# Tier-2 并行 filtered read — 实施与实战记录 2026-09-16

承接：`dsmcFoam_todolist.md` Tier-2 条目（构造期按 cellOwner 过滤初始
cloud，读取时过滤）；串行保守版记录见
`detail_mix/tier2_owned_filtered_initial_read_20260914.md`。

## 驱动问题

超算 ocm-moss3d/230wcell-1bparticle（2,293,832 cells，初始 cloud
116,359,657 parcels，m24o16）：串行版过滤读取单线程构造 **3.5h+ 无输出**
（job 4638184），用户报告"停滞"。实测 24 ranks 每线程满载推进——非死锁，
是构造期单线程磨：per-particle ~105µs（parse + locate，locate 随 mesh
规模劣化），且 93% 构造完即删除（非 owned）。

## 实现（用户指定两项叠加：OMP 分块并行 + 精确预过滤）

### 1. OMP 分块解析（positions 流分块）

- `IOPosition::readHeader(entriesStart)`（新公有方法）：读 header、返回
  条目数 N、给出首条目字节偏移（readStream 后 stdStream().tellg()）。
  readStream 重载选择 `readStream(this->type())`——两棵 OpenFOAM 树的
  无参 readStream 访问级别不同（本地 public/超算 inaccessible），有参
  重载两边都可用且正确校验 positions 的 class="Cloud"。
- 条目区域单遍行扫描（memchr，7GB ≈ 5s）定位每 chunk 首条目字节偏移；
  chunk 按 ceil(N/T) 条目划分。
- 每线程独立 IFstream + seekg 到 chunk 边界；并行段内每条目
  `new ParticleType(mesh, is, false)`（流构造器读 `(x y z) cellI`，
  cell/tet 语义与原路径一致）。
- **mmap**：文件映射共享页缓存（同节点 12 ranks 只读一次盘），避免
  24×17GB 的重复文件 I/O 与逐 IFstream 读放大。

### 2. cellI 精确预过滤（比计划的均匀网格查表更优）

positions 文件每条目自带 `cellI`（流构造器读 `position_ >> cellI_`）——
**无需 findCell 即可精确判定归属**：`cellOwner[cellI] == filterRank`
不满足的条目**连构造都不发生**（计划的均匀网格 O(1) 查表是为"位置无
cell 信息"设计的保守方案；cellI 直接可用使其不必要）。

### 3. 其他

- 进度打印：atomic 计数，每 2000 万条 `Cloud filtered read: parsed
  X / N`（OMP critical 内 Info）。
- keep（文件索引）按 chunk 前缀合并，云布局确定性（chunk 序 merge）。
- 并行路径跳过 initCloudPostRead 的逐粒子 tet 重初始化
  （initCloudPostRead(reLocate=false)，线程内已做）。
- 小文件（<100MB）或单线程自动回退串行 readDataFiltered。

## 调试历程（四个坑，全部实证）

1. **`regIOobject::readData` 默认 `return false`**（GeometricField 未
   重载）——`field.read()` 静默无效（sigmaTcRMax 修复时已踩，本次未再犯）。
2. **裸 `typeName` 名字查找陷阱**：IOPosition 内 `readStream(typeName)`
   解析到继承的 `regIOobject::typeName()`（="regIOobject"），与
   positions header 的 "Cloud" 不匹配 → `unexpected class name` FATAL。
   修复：checkClass=false（与 dsmcCloud 原路径一致）或显式 `this->type()`。
3. **不带 -parallel 时构造期 `Pstream::nProcs() == 1`**：replicated
   mesh 的 MPI 在 `dsmcReplicatedMesh::initialize()` 内才手动 init——
   helper 的 nProcs>=2 条件在构造期恒假（early-return 编号诊断
   `8 × "#4 nProcs<2"` 一锤定音）——删除该条件（单 rank 过滤全保留，
   等价无害）。
4. **`getNewParticleID` 数据竞争**：并行构造使 `particleCount_++` 非原
   子 → origId 竞态。fork 内 particleI.H 改 `#pragma omp atomic capture`。

另有 most vexing parse（`List<char> buf(label(...))`）、`List<off_t>`
两参构造歧义（int 匹配 InputIterator 模板，需 off_t(-1) 显式化）、
regIOobject::readStream 无参重载访问级别两树不同等编译期坑。

## 本地回归（zb mpi4omp2，300 步，阈值临时 20MB 触发并行路径）

| run | 路径 | 结果 |
|---|---|---|
| tier2-par4 | 并行（2 线程） | parsing 1,306,745 / kept 201,681；End main；stuck=0；NoColl=0 |
| tier2-final-off | 串行 | End main；particles 1,958,118；energy 1.96261e-3 |

物理与串行版、6/25 历史基线一致（噪声带内）。

## 超算实战（230wcell-1bparticle，job 4641127，m24o16 双节点）

```
Cloud filtered read: parsing 116359657 entries with 16 threads
Cloud filtered read: parsed 20000000 / 116359657   ← 每 2000 万一条
...
Cloud filtered read: kept ... / 116359657
```

| 指标 | 串行版（job 4638184） | 并行版（job 4641127） |
|---|---|---|
| 116M 条目解析 | **3.5h+ 无任何输出** | **~10 分钟**（每 2000 万可见进度） |
| 每 rank RSS | （未完成，预测 ~87GB） | **实测 15.6 GB**（Threads=17、R 状态、CPU 推进确认） |
| 24 ranks 总 RSS | （未完成） | **实测 205 GB**（节点 1536GB 的 13%） |

对照：若不过滤（原始全量物化路径），116M × 24 ranks ≈ 2.2TB——物理不可
行。**Tier-2（过滤 + 并行）是该规模 case 的使能项，而非可选项。**

## 遗留

1. **bulk 字段读取仍串行 I/O**（17GB/rank × 24 = 408GB Lustre）——当前
   解析后的下一瓶颈（构造+字段合计 ~40-60 分钟）；可并行化（字段间独立）
   或纳入 Tier-2 完全体（分块流式，消全量 IOField 瞬态）。
2. `readHeader` 的 Info 诊断打印已移除；保留 `Cloud filtered read:
   parsing/kept` 两级正式输出。
3. 本地阈值测试补丁（20MB）已恢复 100MB。
4. 时间循环运行状态待观察（日志 stdio 缓冲滞后，监控用 sstat/ps/CPU）。

## 追加：优化 3b（in-solver 并行填充）实施 2026-09-17

设计见 `distributed_initial_cloud_omp_fill_design_20260917.md`（P1-P3）。

### 实现（P1-P2）

- `dsmcConfiguration.H/.C`：虚方法 `setInitialConfigParallel(cellOwner,
  filterRank)`（基类默认 FATAL，提示走串行路径）
- `dsmcMeshFill.H/.C`：override 实现——OMP dynamic 遍历 owned cells，
  **per-cell 确定性 RNG**（seed = cellI×黄金比 ^ rank 混合，与线程/调度
  无关），采样链全部走 Random& 重载（equipartitionLinearVelocity/
  Rotational/Vibrational/Electronic 四个重载：dsmcCloud.H/.C 新增，
  原无参版委托），sigmaTcRMax 初值内算（串行版尾部公式提取），
  addNewParcel 直接复用（内部 critical addParticle 已 OMP-ready）
- `dsmcAllConfigurations.H/.C`：`setInitialConfigParallel` 编排
- `dsmcCloud.H/.C`：`initialFillFromDict`（ctor-2 填充逻辑接入 ctor-1
  路径）+ `dsmcInSolverFill` 触发 helper + buildConstProps 前移
  （填充分支在 constProps 构建前执行的 FATAL 修复）

### 本地验证（P3，zb mpi4omp2 300 步，0/lagrangian 移除触发填充）

```
dsmcMeshFill parallel fill: rank 0 fills 15022 owned cells
dsmcMeshFill parallel fill: rank 0 inserted 201247 particles
Replicated mesh: rank 0 in-solver fill: 201247 particles
Total Iterations = 300 / End main / stuck=0 / NoColl=0
particles 1,957,357 / energy 1.96138e-3（基线噪声带内）
```

Σ 插入 804,988 ≈ dsmcInitialise+ 生成的 806,724（差 0.2% =
particlesRequired 取整的 Poisson 噪声）✓ 守恒按构造精确。

### 超算 3b 实测（230wcell-1bparticle，job 4644800，排队中）

- 重建 0 错误；0/lagrangian（17GB 串行产物）已备份
  `0/lagrangian.serial-bak`；slurm 的 dsmcInitialise+ 步骤注释
  （run-hpc-fill.slurm）；controlDict `dsmcInSolverFill true`
- 启动后预期：填充 1.16 亿粒子（384 核）分钟级，**初始化 6.5h → <5 min**
  （消除 408GB 字段 I/O 与全部文件往返）

### 教训（工程过程）

- 多文件 python 批量 patch 的顺序依赖：前一步插入使后一步的
  锚点失效 → 部分改动静默未落盘。**每次批量 patch 后必须 grep 验证
  每个插入点实际存在**（本次靠编译错误逐一暴露补齐）。
- dsmcCloud.C 的 equipartition 变换曾破坏文件（参数插入位置错误），
  从 git HEAD 干净文本 + 损坏文件完好前缀（L1-3452）重建。

## 追加：3b 超算验收实测（2026-09-17，job 4644800，m24o16 双节点）

230wcell-1bparticle（2,293,832 cells、1.16 亿初始粒子）：

| 指标 | 旧路径（Tier-2 过滤读，job 4641127） | **3b in-solver 填充（4644800）** |
|---|---:|---:|
| Iteration 10 @ ExecutionTime | 23,300 s（6.47 h） | **200.76 s（3.3 min）** |
| ClockTime @ iter 10 | 23,225 s | **36 s** |
| rank 0 填充 | kept 9,111,898（过滤读后） | 95,576 cells / **8,962,116 particles** |
| MaxRSS（24 ranks） | 205 GB | 218 GB |
| 1.16 亿粒子物理状态 | stuck=0、碰撞正常 | stuck=0、碰撞 33,561@iter10 窗口 ✓ |

**初始化 6.47h → 3.3 min（×118）**——远超 <5 min 验收目标。
ExecutionTime/ClockTime = 5.6 ≈ OMP 并行填充并行度 ✓（旧路径 ≈1.0 串行）。
rank0 8.96M 份额为均值 1.85×（METIS 初始失衡，DLB 后续拉平，与
100wcell 观察一致）。作业验证完成后由用户 scancel（验证目的达成）。

**结论：230w/1.16 亿粒子的初始化问题闭环。1B 粒子级（优化 3b 必需项）
已验证。优化 1/2 的原始动机（初始化耗时）由 3b 根除。**

## 追加：log 重复输出去重（2026-09-17）

### 问题

.replication mesh 模式（mpirun 无 -parallel）下，启动期回显（reaction
模型选择 ×192、dsmcVolFields 配置 ×144、porous/timeStep/patchBoundary/
Constructing 等）每 rank 重复打印到单一 .out。

### 修复

1. **10 个回显站点 master 门控**（reaction/field/timeData/timeFluxData/
   timeDataMeas/porous/constantTimeStepModel/variableTimeStepModel/
   patchBoundary/dsmcFoam+ Constructing）
2. 门控谓词 **`dsmcIsPrintingRank()`**（新头 `basic/particle/
   dsmcMasterInfo.H`）：MPI 已初始化 → MPI_Comm_rank；未初始化 →
   **PMI_RANK 环境变量**（mpirun 启动的构造期，Pstream 未 init、
   Pstream::master() 每 rank 恒真——不能用）；真串行 → true
3. 头文件分发：basic/particle（lnInclude 平铺）+ time/timeData +
   time/timeFluxData + time/timeDataMeas 同目录副本（各库 include 路径
   可见）
4. slurm 模板建议：去掉 I_MPI_DEBUG=4（MPI startup ×24 回显源）

### 两个工程坑（实测教训）

1. **wmake 依赖跟踪不穿透 lnInclude symlink**：头文件内容变化（symlink
   目标重写）不触发依赖重编——touch 源文件强制。build 脚本的
   `wmake -j lnInclude 2>/dev/null || true` 是无效语法（从未真正更新
   lnInclude）——已改为 `wmakeLnInclude .`
2. **批量 patch 的锚点顺序依赖**：先插入 initialFillFromDict 会使
   helper patch 的锚点（尾部含 Constructors 标记）失配——**每步 patch
   后 grep 验证**（本次靠编译错误逐一暴露）

### 本地验证（zb mpi4omp2）

| 回显 | 门控前 | 门控后 |
|---|---:|---:|
| Selecting the reaction model | 48 | **12**（12 模型各 1） |
| Selecting field: | 24 | **6**（6 字段各 1） |
| TimeData Statistics | 24 | **6** |
| Constant time-step model | 4 | **1** |
| Constructing dsmcCloud | 4 | **1** |
| Selecting the porous | 4 | **1** |
| Initialising dsmcVolFields | 6 | **6**（本已单份） |

物理不变（1,957,946 / 1.96203e-3 噪声带内）。门控函数返回值经 stderr
诊断实证（MPI 路径 rank=0→1、rank=1/2/3→0；PMI 路径同）。

## 遗留（待 4647240 启动后验证）

- 超算 230wcell-1bparticle 重跑（run-hpc-fill.slurm，3b + 去重）：
  初始化 <5 min + 重复消除双重验证
- dsmcInitialise+ 残留：0/lagrangian.serial-bak 保留可回退

## 追加：log 去重第二轮（2026-09-20）

超算 test case（job 4738920，m16）实测：第一轮门控站点全部生效
（TimeData Statistics 6=实例数 ✓，含 libgeneralMolecule 补编），
但仍有 6 类未门控回显（nSamples/nAverages/measurement option/
total no. of sampling steps ×96、time-step/coordinate-system/
generalBoundary/BinaryCollision 的 Selecting ×16、reactions created ×16）。

### 本轮新增门控（7 个站点，全部仓库代码）

- timeDataMeas.C：setInitialData 内 measurement option/nSamples/nAverages
  块 + total no. of sampling/averaging steps
- timeData.C：total no. of sampling/averaging/control steps
- dsmcTimeStepModel.C / dsmcCoordinateSystem.C / dsmcGeneralBoundary.C /
  BinaryCollisionModel.C 的 RTS New 选择打印
- dsmcReactions.C：Number of reactions created / no reactions defined

### OpenFOAM 基树（argList banner/sigFpe/dynamicFvMesh）——已回滚

曾短暂给基树 argList.C/sigFpe.C/dynamicFvMeshNew.C 加 PMI 门控并验证
有效（banner/nProcs/sigFpe/dynamicFvMesh 全部单份），按决定**不修改
OpenFOAM v1706 源码**，已全部逆向恢复并重建 libOpenFOAM/libOSspecific/
libdynamicFvMesh 至原始状态（基树 PMI 串归零）。这三处 ×N 回显为基树
行为，接受现状。

### 工程坑（新）

- **基树 OSspecific 的产物是可重定位对象 libOSspecific.o**（被
  src/OpenFOAM/Make/options 直接链入 libOpenFOAM），`wmake libso` 产出的
  libOSspecific.so 不被任何程序加载；且链接器合并 .rodata 同名字符串，
  `strings | grep -c PMI_RANK` 无法区分门控是否编入——验证必须以
  运行行为为准。
- wmake 不因依赖库更新而重链 exe（动态链接运行时自然生效，无需重链）。

### 验证（zb mpi4omp2 300 步）

| 回显 | 前 | 后 |
|---|---:|---:|
| Selecting the reaction model | 48 | 12（12 模型） |
| Number of reactions created | 4 | **1** |
| Selecting field:/TimeData/nSamples/measurement option/total steps | 24 | **6**（6 实例） |
| time-step model/coordinate system/generalBoundary/BinaryCollision | 4 | **1** |
| Build banner/sigFpe/Selecting dynamicFvMesh（基树） | 4 | 4（接受） |

物理不变（1,957,627 / 1.96166e-3 噪声带内）。7 文件已同步 bjm8 并重编
（generalMolecule/dsmc/exe，0 错误，.o 时间戳 16:35 验证）。
