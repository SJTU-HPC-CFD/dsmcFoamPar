# DLB 进一步优化（1300 万 cell 级）— 2026-09-20 ～ 2026-09-21

> **状态总览（2026-09-21 收口）**：本档案记录 1200w-cell DSMC（dsmcFoam+，
> replicated mesh + ParMETIS DLB）的 DLB 与热路径优化全周期。**优化阶段已
> 收口**：全部可行动方向已落地或证伪，当前生产配置 = staged 基线 + P0–P5
> + back-off/wall 触发 + fastRng（生产 40000 步 19h12m 正常完赛）。

## 最终结论索引

**已落地（生产配置）**

| 项 | 内容 | 效果 |
|---|---|---|
| F1–F5（§4） | DLB OMP 化、findIndex 位图化、迁移记账去除、两次挂死根因修复 | DLB #1 从 2617s → 30s 量级 |
| §5 | back-off + 每步完整 walltime 触发指标 | 误触发削减（§7 A/B） |
| §6.3/§6.4 → §8 P0–P2 + ① | post cache 持久化+稀疏复位、build 串行段 OMP 化、boundary accumulation OMP、build-fill atomic-cursor | 同窗 −57% build、−25% post（§8.9） |
| §7 | 迁移同步等待相关（delayed receive 等） | 同段 A/B 确认 |

**已证伪关闭（均有 A/B 数据）**

| 项 | 章节 | 判定 |
|---|---|---|
| 双约束权重（ncon=2） | §9.3 | 更差（wall +18s，分布畸变） |
| adaptive alpha | §9.3 | 中性，未激活 |
| cost model 权重（file-static 版） | §10.8/§10.9 | 分区质量回归 +51%→+275%；根因 act 固定项占权重 74% |
| post-move 时序重叠 | §9.2 | 无可重叠窗口 |
| NoAlltoall | §6.2 方案A | +45% |
| collision-aware DLB | §9.3 | 碰撞热点非 wall 瓶颈（0.03s/步） |
| 边界 critical 免锁 | §11.6.4 | 墙钟代价不显著，可动份额 <8% |
| per-segment cell 排序 | §11.7 | move +37%，排序成本 > 局部性收益 |

**机理闭环（§11）**：move 相位 per-parcel 成本空间差异 1.63× = interior
trackToFace 内部成本差 1.57×（主因，与粒子密度负相关，cache/访存属性，
corr=0.997）+ 复杂几何切割 cell 放大 3–4×（占比 2–8%）。均为物理/网格/
访存属性而非 cell 归属属性 ⇒ DLB 权重原理上无法修复——§9/§10 失败的
根本原因。tracking 拓扑优化被"不改 OpenFOAM 源码"约束封死。

**后续工作建议**：转向生产任务推进；本文档不再新增优化方向。

## 目录

- §0–§1 背景与根因（4738842 DLB 开销事件）
- §2–§4 F1–F5 修复与超算验收（4746677 完赛）
- §5 back-off + wall 触发指标设计（已实施，§7）
- §6 100wcell 分析与三项优化设计（优化一无效；二/三 → §8）
- §7 §5/§6 实施与 A/B（4748993/4749270）
- §8 post/build P0–P2 + 追加三项（完成，含累计成果）
- §9 权重重测 + 时序分析（双关闭）
- §10 cost model（实施→回归→归一化修复→证伪关闭）
- §11 move per-parcel 成本研究（分账实验→两方向处理→收口）

---

## 0. 背景与驱动问题

bjm8 超算 `ocm-moss3d/1200w-rcj`（12,958,378 cells，初始 1.5 亿 parcels），
作业 4738842（ocmmoss-m12o48，12 rank × 48 线程，576 核 3 节点，3 万+ 步长跑），
用户报告"实时 CPU 利用率 100%，卡在 DLB"。

日志取证结论：**不是死锁，是 DLB 真实开销**。

| 阶段 | 耗时（ClockTime） |
|---|---:|
| step 10–40（失衡 3.98 未修） | ~1.25 s/步 |
| step 50：DLB #1（threshold 触发，loadImbalance=3.984） | **2617 s** |
| step 60–240（重平衡后） | ~1.08–1.10 s/步，CPU 利用率 36%→66% |
| step 250：DLB #2（threshold 触发，loadImbalance=1.595） | 进行中（remap 后仅 217,938 cells 变更） |

关键运行数据：

- DLB #1 remap（首次生产验证）：changedBefore=10,273,985 → changedAfter=5,775,458，
  **saved=4,498,527 cells 的迁移量**
- DLB #2 remap：changedBefore=10,977,995 → changedAfter=217,938，
  **saved=10,760,057**（划分已接近最优，remap 把伪迁移几乎全部消除）

## 1. 根因分析（代码证据）

### 1.1 DLB 划分侧没有 OMP

`dsmcReplicatedMesh.C` 全部 3 处 `#pragma omp` 都在
`migrateParticlesByCellOwner()`（行 3439/3526/3872，即 M1 的
owner 分类 / pack / unpack 并行化）。DLB 其余环节每 rank 单线程执行：

- 全局每 cell 粒子数构建（13M 循环）+ `MPI_Allreduce`（52 MB）
- ParMETIS 权重构造 vwgt/adjwgt/vsize（13M/1.08M 循环）
- `ParMETIS_V3_AdaptiveRepart` 本体：**纯 MPI 库，12 rank 并行度，
  每 rank 108 万顶点，48 个 OMP 线程闲置**
- greedyOverlap remap：overlap 统计（13M 循环）+ changed 统计（13M×2）
- `cellOwner_` 更新 + nChanged（13M）
- `rebuildMyCells()`（13M）
- `validateCellOwnerMap()`（13M×2 + 68 万 HashSet）

### 1.2 `dsmcLocalMesh::build()` 的 O(N²) 热点（最大单项串行开销）

`dsmcLocalMesh.C:86-93`：对 nTotal（rank0 ≈ 73.4 万）个局部 cell 逐个调用
`findIndex(myCells, globalI)`，在 68 万元素的 `myCells` 数组里线性查找：

```text
73.4 万 × 34 万（平均扫描半长） ≈ 2.5×10^11 次标量比较
单线程、内存带宽受限 → 估算 250–800 s
```

这正是 `todo.md` 已登记的待优化项（"全局 List<bool> 和线性 findIndex()"）。
每次 DLB repartition 后必然执行一次（`dsmcReplicatedMesh.C:1756`）。

### 1.3 迁移接收端每 parcel 的冗余记账

flat POD 接收循环（`dsmcReplicatedMesh.C:3894-3898`）对每个接收 parcel 调
`cloud_.addParticle(received[j])`。`dsmcCloud::addParticle`（dsmcCloud.C:1370）
每次做分支 + `moveAppendedParcels_.append` + 4 个 cache 标志写。而该路径
收尾逻辑（`dsmcReplicatedMesh.C:3952-3959`）必然以
`setMoveOrderedParcels(kept)`（内部 clear moveAppendedParcels_ 并重置全部
cache 标志，dsmcCloud.C:2936-2966）或 `clearMoveOrderedParcels()` 收尾——
**每 parcel 的记账在最终状态上完全冗余**，只浪费串行时间。

### 1.4 附带发现（当时不修，记录在案）——后三项均已处理

> 处置状态：①CPU 时间触发语义偏移 → 已由 §5.6 每步完整
> walltime 指标取代（§7.1 实施）；②阈值抖动/MinGapSteps →
> 已由 §5 back-off 三层节奏实施；③remap 循环 FatalError 与
> OMP 并行区的冲突 → 本次实现已处理。

- 触发指标 `productiveTime_` 来自 CPU 时间（`dsmcCloud.C:2658-2661` 传
  `elapsedCpuTime`），48 线程 barrier 自旋等待计入 → loadImbalance 语义
  偏移；`replicatedMeshDLBParticleGate false` 时日志打印的
  `particleMaxMin=1` 是未参与计算的占位值（`particleGateImbalance`
  初始化值，dsmcReplicatedMesh.C:1935），不反映真实粒子失衡。
- DLB #2 在划分已最优时仍触发（1.595 vs 阈值 1.5，MinGapSteps=50），
  长步数 case 存在阈值抖动风险；6 月报告建议长跑 MinGapSteps 300–1000。
- remap 后 post-remap 循环（dsmcReplicatedMesh.C:1699-1716）内有
  FatalError，OMP 化需把错误检查移出并行区（本次实现已处理）。

## 2. 修复方案

约束：**不修改 OpenFOAM v1706 源码**
（`~/code/OpenFoam/OF-1706/OpenFOAM-v1706` 只读；核心容器
`DLListBase::first_/last_/nItems_` 私有、`transfer()` 为替换语义，
不修改核心即无法对非空侵入式链表做 O(线程数) splice）。全部改动限于
本仓库 `src/lagrangian/dsmc/replicatedMesh/`。

### F1 `dsmcLocalMesh::build()` findIndex → O(N) 成员位图

- 新增 `List<bool> isOwnedCell(nGlobalCells)`，标记 `myCells`；
  halo 判定由 `findIndex(myCells, globalI) < 0` 改为
  `!isOwnedCell[localToGlobalCell_[localI]]`。
- O(nOwned × nTotal) → O(nGlobalCells + nTotal)，单线程约 250–800 s →
  ~30 ms。行为逐位等价（同一集合的成员测试）。

### F2 13M-cell 级串行循环 OMP 化（M1 模式）

线程数口径照抄 M1/现有接收路径：`cloud_.ompNumThreads()`（controlDict
`openmpThreads`），并用 `if(nT>1)` 子句在单线程/无 OMP 时自动退化为
串行路径，避免重复循环体。逐项：

| 循环 | 位置 | 并行方式 |
|---|---|---|
| localCellParticles 构建 | reassignByParMetisAdaptiveRepart 开头 | OMP for，写索引不相交 |
| rank0 统计（particles/activeCells/max） | 同上 | reduction(+,+:)(max:) |
| vwgt 构建（ncon 两分支） | 权重段 | OMP for，按 i 不相交 |
| part 初始化 / vsize 构建 | ParMETIS 参数段 | OMP for |
| changedBeforeRemap 统计 | remap 前 | reduction(+) |
| remap overlap 统计（13M） | greedyOverlap | 每线程扁平 nProcs×nProcs 矩阵 + 串行合并（M1 两遍模式） |
| post-remap fullPart 重写 + changedAfter | remap 后 | OMP for + reduction(+)；FatalError 改为记录坏索引、并行区外报错 |
| cellOwner_ 更新 + nChanged | 收尾 | OMP for，写不相交 + reduction(+) |
| rebuildMyCells | 独立函数 | 两遍法：并行计数 → 前缀 → 并行填充（保持升序） |

明确不并行（记录原因）：

- `localDeg`/`adjncy`/`adjwgt` 构建：按 face 遍历、每 face 对
  `off[]`/计数器自增，多 face 写同一 cell 槽位 → 需要原子或逐线程重排，
  改动风险大于收益（每遍 26M face 串行 ≈ 0.3–0.5 s）；
- `validateCellOwnerMap`：FatalError 在循环内 + HashSet 去重，串行
  13M ≈ 0.1 s，保持原样；
- ParMETIS 本体：外部库，纯 MPI，不可 OMP。

### F3 迁移接收端：去除每 parcel 冗余记账

- flat POD 接收循环中 `cloud_.addParticle(received[j])` 改为基类裸追加
  `cloud_.Cloud<dsmcParcel>::addParticle(received[j])`（跳过
  dsmcCloud 层每 parcel 的分支/标志/append 记账）；
- 最终 cache 状态由既有收尾 `setMoveOrderedParcels(kept)` /
  `clearMoveOrderedParcels()` 完整重建（语义不变，见 §1.3）；
- `kept` 逻辑与 `kept.size() == cloud_.size()` 一致性检查不变。
- **真并行 splice 的取舍**：DLListBase 头/尾/计数私有、`transfer()`
  是替换而非合并，不改 OpenFOAM 核心无法 O(1) 合并非空链表；且实测
  该串行段本身仅 ~0.1–0.3 s/DLB（67M parcels 分摊到 12 rank 后每
  parcel ~20–50 ns）。故本项定位为"去除冗余开销"，不是全并行。

## 3. 验证计划

1. 编译：`source doc/scripts/build-dsmcFoam.sh`，产物
   `platforms/linux64IccDPInt32Opt/bin/dsmcFoam+`（`which` 校验）。
2. 本地回归：`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react-validate/mpi4omp2`
   （mpirun -np 4 × OMP 2，300 步），门检查：
   `Total Iterations = 300`、`End main`、stuck=0、碰撞非零
   （NoColl=0，防 sigmaTcRMax 类假阳性）、末态 particles/energy 落在
   6/25 基线噪声带（1,957,715 ± / 1.9619e-3 ±）。
3. 对照运行日志确认 `Replicated mesh: initialized with 4 MPI ranks`、
   `OpenMP max threads = 2`（AGENTS §4.5）。
4. 超算同步（dlb-new）：重跑 1200w-rcj 观察 DLB 步耗时（预期 F1 生效后
   单次 DLB 显著缩短）与 remap/触发行为不变。

## 4. 实时记录

### 4.1 修复实施（2026-09-20 完成）

- [x] 根因定位与方案设计（本文档 §1–§2）
- [x] F1 `dsmcLocalMesh.C`：halo 判定改为 `isOwnedCell` 位图，删除
      `findIndex(myCells, globalI)` 线性查找
- [x] F2 `dsmcReplicatedMesh.C`：
  - `reassignByParMetisAdaptiveRepart()` 内 `nT`（controlDict
    `openmpThreads`，`openmpEnabled()` 门控）+ 8 处循环 OMP 化：
    localCellParticles 构建（`cellOccupancy()` 引用先串行物化再并行遍历，
    规避惰性物化并发）、rank0 统计（reduction +/max）、vwgt、part、
    vsize（两分支）、changedBeforeRemap、remap overlap（每线程扁平
    nProcs×nProcs 矩阵 + 串行合并）、post-remap fullPart 重写 +
    changedAfter（FatalError 移出并行区，坏索引 critical 记录后区外报错，
    防御条件 `remap[newPart]<0` 保留）、cellOwner_ 更新（reduction）
  - `rebuildMyCells()` 两遍法（并行计数 → 串行前缀 → 并行填充），
    升序保持，`nCells >= 2*nT` 以下走原串行路径
  - 未并行并记录原因：localDeg/adjncy/adjwgt（face 遍历存在跨 face
    写冲突）、validateCellOwnerMap（FatalError 在循环内，串行仅 ~0.1 s）
- [x] F3 `dsmcReplicatedMesh.C` flat POD 接收循环：
      `cloud_.addParticle(...)` → `cloud_.Cloud<dsmcParcel>::addParticle(...)`
      （基类裸追加；dsmcCloud 层每 parcel 记账由收尾
      `setMoveOrderedParcels(kept)` / `clearMoveOrderedParcels()` 完整重建，
      语义不变）
- [x] 编译：`source doc/scripts/build-dsmcFoam.sh` exit=0；
      dsmcLocalMesh.o / dsmcReplicatedMesh.o 重编、libdsmcFoam+.so 重链
      （2026-09-20 17:14）；无新增告警（仅有既有 OpenFOAM 头告警与既有
      ctor reorder 告警）
- [x] 二进制校验：`which dsmcFoam+ dsmcInitialise+` → 本仓库
      `platforms/linux64IccDPInt32Opt/bin/`；`mpirun` → oneAPI 2021.16

### 4.2 本地回归（zb-cylinder-react-validate/mpi4omp2，300 步）

运行命令：`source doc/scripts/env.sh; OMP_NUM_THREADS=2 OMP_DYNAMIC=false
OMP_PROC_BIND=close OMP_PLACES=cores /usr/bin/time -p mpirun -np 4 dsmcFoam+`
（real 116.9 s，exit 0）

Run 1 `log.dlb-opt-20260920`（case 原样，remap off）门检查：

| 门 | 结果 |
|---|---|
| `Total Iterations = 300` / `End main` | ✓ |
| `Replicated mesh: initialized with 4 MPI ranks` / `OpenMP max threads = 2` | ✓ |
| stuck particles | 0 |
| "No collisions" 计数 | 0（碰撞非零，末窗口 182,121） |
| final particles / Total energy | 1,958,335 / 1.96255e-3 |
| 6/25 基线带（1,957,715–1,958,118 / 1.9619e-3–1.9626e-3） | 带内 ✓ |
| DLB | 2 次 rebalance 完整执行（SAR 触发 step 100/200，step 300 被
  MinRemainingSteps 拦截）；owner map validated ×2 |

Run 2 `log.dlb-opt-remap-20260920`（临时 controlDict 追加
`replicatedMeshDLBRemap true;`，跑完已恢复备份）：exit 0，同一门全过
（1,958,080 / 1.96248e-3）；remap 路径两次执行：
`changedBefore=30000→changedAfter=3223`、`59998→6453`，OMP remap
（每线程 overlap 矩阵 + post-remap 重写）行为正确。

### 4.3 计时对比（单次运行，仅作回归参考，不构成性能结论）

| run | total_loop | move | coll | build |
|---|---:|---:|---:|---:|
| tier2-final-off（基线） | 108.9 | 78.2 | 22.5 | 5.17 |
| sigma-fix3（基线） | 126.7 | 87.5 | 23.3 | 5.85 |
| dedup4（基线） | 96.5 | 59.8 | 15.1 | 4.39 |
| **本次 dlb-opt** | **89.0** | **62.7** | **22.5** | **4.25** |

- 各 phase 均落在近期噪声带内，无回归迹象。
- 注意：`log.threadoffset-fix-20260914` 的 total_loop/move 不可比——该
  run 是 sigmaTcRMax 零碰撞 bug 时期的假阳性验证（coll=0.25 s 为零碰撞
  特征），不能作为相位基线。
- 本地 case 仅 6 万 cells，F1/F2 的收益在该规模不可测（F1 的 O(N²) 项在
  6 万 cells 下仅毫秒级）；性能收益验证必须依赖超算 1200w-rcj 重跑。

### 4.4 结果分析与超算验证清单

分析：

1. 正确性：两轮 300 步回归（含 remap on/off）全部门检查通过，末态
   particles/energy 落在 6/25 基线噪声带；DLB 触发行为（SAR 触发步、
   MinRemaining 拦截）与修复前一致。
2. 预期收益归因（1200w-rcj 量级）：
   - F1：单次 DLB 的 localMesh 重建从 ~250–800 s 降到毫秒级——DLB #1
     的 2617 s 中预计削减数百秒；
   - F2：13M 级循环（localCellParticles/vwgt/vsize/overlap/owner 更新/
     rebuildMyCells）在 48 线程下从串行 ~1–3 s 量级降到 ~0.1 s 级；
     ParMETIS 本体（纯 MPI）与 face 三循环（串行保留）不受影响；
   - F3：接收端每 parcel 记账消除，~0.1–0.3 s/DLB。
   - 均为 DLB 步内收益；每步 move/coll/build 路径未改动，正常迭代耗时
     预期不变。
3. 遗留（本次不修，见 §1.4）：触发指标 CPU 时间语义、阈值抖动
   （MinGapSteps=50 偏小）、face 三循环串行、validateCellOwnerMap 串行。

超算同步与验证（dlb-new 对应位置）：

```text
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C   # F2 + F3
src/lagrangian/dsmc/replicatedMesh/dsmcLocalMesh.C        # F1
```

### 4.5 超算同步记录（2026-09-20 17:31，已完成）

- 目标：`bjm8:/publicfs01/fs1-m8/home/m8s000774/users/xiao_chen_xiang/
  611/dsmcFoamPar/v1706/dlb-new/src/lagrangian/dsmc/replicatedMesh/`
- 同步前远端 md5（修复前版本，已备份为同目录
  `*.C.bak-20260920`，非编译后缀，不参与 wmake）：
  - dsmcReplicatedMesh.C `c6ab7d77e3ba5a94f3811afaaee19501`
  - dsmcLocalMesh.C `3c6a97bd877e8d91c6fc4758ef7889a0`
- 同步后远端 md5 与本地一致：
  - dsmcReplicatedMesh.C `510c7208e39e1c147d9259eddb246e81`
  - dsmcLocalMesh.C `db8313f6e25de0cbedb136a811968af8`
- 待办：远端 `dlb-new/build.sh` 重建 → 重跑 1200w-rcj 验证（作业
  4738842 当前仍在跑旧二进制，重跑前确认作业状态/杀旧作业由用户决策）。

重建后重跑 1200w-rcj（建议同时记录）：

- DLB #1（step 50）迭代耗时：修复前 2617 s → 预期显著下降（F1 主导）；
- `Phase C ParMETIS profile` 开关（`replicatedMeshDLBProfile true`）可
  拆分 weights/AdaptiveRepart/remap 分项，确认剩余瓶颈是否已转移到
  ParMETIS 本体与迁移；
- remap/触发行为与本次日志逐项对照（changedBefore/After、触发步、
  rebalance 次数）；
- 正确性门：`Total Iterations`、`End main`、stuck=0、碰撞非零、
  粒子数/能量与 4738842 对照。

### 4.6 重新诊断（4744935 DLB #1 挂死）与 F4 修复（2026-09-20 18:00）

#### 现象

新二进制（F1–F3 已同步重建，libdsmcFoam+.so 17:33）作业 4744935
（同 case，m12o48）在 step 50 DLB #1 触发后"卡住"，top 显示每 rank
进程 ~95–100% CPU（单线程）。

#### 取证（gdb 计算节点逐 rank 抓栈，12/12 全采样）

- **12 个 rank 全部处于 `validateCellOwnerMap()` 内部的
  `MPI_Allreduce` 自旋**（栈：PMPI_Allreduce → hcoll
  `hmca_bcol_ucx_p2p_allreduce_knomial_progress` → UCX
  `uct_rc_mlx5_iface_progress_cyclic`），无掉队 rank。
- 集体操作锁步 + 全员已进入 + 20 分钟不完成 ⇒ **MPI/hcoll 集体通信
  本身挂死**，不是应用层串行代码，也**不是 OMP 未开启**：
  move 阶段每迭代 CPU 47 s / wall 1.2 s（48 线程 ~65% 占用）证明
  OMP 正常；DLB 期间 top 单线程 100% 是全员等集体时的自旋表现。
- 对照旧作业 4738842 日志：**DLB #2（step 250）同样卡在此处**（日志
  终止于 remap/"rank 0 owns"行，无 "owner map validated"；.out 16:18
  后停更，16:36 作业终止）——用户最初报告的"卡在 DLB"即此挂死，
  当时被误判为"纯慢"。DLB #1（旧）能通过 ⇒ 挂死为间歇性。

#### 根因定位

`validateCellOwnerMap()` 的全局 map 一致性校验用分块 MIN/MAX 全量
Allreduce：13,958,378 cells / 256K 块 × 2（MIN+MAX）= **~100 次
1 MB 集体通信**（todo.md 登记的"逐元素分块归约 DLB 检查成本"）。
该序列在 bjm8 的 hcoll/UCX（knomial bcol）路径上间歇性挂死；
挂死暴露面随集体次数线性放大。

#### F4 修复

`dsmcReplicatedMesh.C` `validateCellOwnerMap()`：分块 MIN/MAX 全量
校验（~100 次集体 + 2×1MB/块缓冲）替换为 **64 位包装和校验和 +
2 次 MIN/MAX Allreduce**：

```text
checksum = Σ_i  (cellOwner_[i]+1) * (0x9E3779B97F4A7C15 ^ i)   （uint64 环绕加法）
MIN(checksum) == MAX(checksum)  ⇒  各 rank map 一致（碰撞概率 ~2^-64 量级）
```

- 环绕加法满足结合/交换律，与 rank 求和顺序无关；
- 本地校验（owner 范围、myCells 一致性、coverage SUM）不变；
- 集体通信次数 ~101 → 3，挂死暴露面降 ~30 倍，同时消除大 case 的
  校验耗时（正常路径下 100×1MB allreduce 也是可测成本）。
- 编译坑（同 mixed 报告 §3.7）：OFv1706 Ostream 无
  `unsigned long long` 重载，FatalError 打印 cast 到 `long`。

#### F4 验证

- 重建 exit=0（libdsmcFoam+.so 2026-09-20 18:00）。
- 本地回归 `zb-cylinder-react-validate/mpi4omp2` 300 步
  （`log.dlb-f4-20260920`）：exit 0；2 次 DLB 的 checksum 校验路径
  各执行一次（"owner map validated"×2）；全部门检查通过
  （300 步/End main/stuck=0/NoColl=0；末态 1,957,813 / 1.96192e-3，
  6/25 基线带内）。
- 同步：`dsmcReplicatedMesh.C`（md5 `1f286b1184db124c580bcc74823a9fac`，
  2026-09-20 18:0x 已同步至 bjm8 dlb-new 并校验一致；`dsmcLocalMesh.C`
  无 F4 改动，沿用 17:31 同步版本）。作业 4744935 已不在队列
  （squeue 无记录，应由用户终止）。

### 4.7 根因更正与 F5 修复（2026-09-20 19:00，最终版）

#### 4746365 / 4746420 仍"挂死"推翻 §4.6 的 hcoll 结论

- F4（集体次数 101→3）+ `I_MPI_COLL_EXTERNAL=0`（4746420 environ 已确认
  生效）双管齐下后两作业仍卡在同一位置（Local mesh 之后、validated 之前）
  ⇒ 挂死与 hcoll 无关。
- 逐 rank 抓栈发现**发散**：多数 rank 在 validate 的第一个 Allreduce 等
  待，**至少一个 rank（60720，CPU ~190%）仍在 validate 的本地代码里
  "磨"**——PC 在相距 76 字节的两条指令间振荡，15 分钟未到达集合点。
- 反汇编该热点（validateCellOwnerMap+1339..+1436）：

```text
mov 0x8(%rax),%rax     ; rax = rax->next_   ← HashTable 桶链 next 指针
cmp (%rax),%r15d       ; 比较节点 key
je  → 命中退出
（编译器展开的链遍历循环）
```

  即 **HashSet<label> 的桶链走查**。

#### 真正根因（F5）

`validateCellOwnerMap()` 的去重校验用 `HashSet<label>`：

1. `Foam::Hash<label>` 单参版是**恒等哈希**（`return p`，Hash.H:97-100）；
2. `HashTableCore::canonicalSize` 把桶数强制为**2 的幂**；
3. 桶索引 = cellIndex 的低位；DLB 后各 rank 的 owned cell 集合
   （ParMETIS 空间连续分区 + greedyOverlap remap 的并集，且按粒子加权，
   远场 rank 可拥有数百万 cell）与 2^22 取模相互作用，至少一个 rank 的
   key 分布使桶链病态拉长 → `insert/found` 退化为 O(myCells²) 链走查，
   单 rank 磨数十分钟到小时级，其余 rank 在集合点"看起来挂死"。

这与 F1 修掉的 `dsmcLocalMesh::findIndex` O(N²) 属同一类问题——
**validate 的 HashSet 是漏网的第二个 O(N²) 点**。历史行为全部自洽：
旧作业 DLB#1（rank0 68 万 cells，磨 ~5-10 min）藏在 2617 s 内完成；
DLB#2/F1-F3/F4 各次则等穿了用户耐心。

**注**：§4.6 的"hcoll/UCX 间歇挂死"结论据此**更正为误判**——所有
"全员卡在 Allreduce"的表象都是等待磨链的 straggler；F4 本身（100 次
集体 → 3 次）仍是有效的成本削减，予以保留，但不是挂死的解药。

#### F5 修复

`validateCellOwnerMap()`：`HashSet<label> listedOwned(...)` →
`List<bool> listedOwned(nCells, false)` 位图（insert = 置位、found =
读位）。13 MB 瞬态对 ~90 GB/rank 占用可忽略；O(1) 成员测试与 key 分布
无关，彻底消除该类退化。三项校验逻辑（范围、归属、去重）不变。

#### F5 验证与状态

- 重建 exit=0（libdsmcFoam+.so 2026-09-20 18:42）。
- 本地回归 `zb-cylinder-react-validate/mpi4omp2` 300 步
  （`log.dlb-f5-20260920`）：exit 0；2 次 DLB 位图校验路径各执行一次；
  全部门通过（300 步/End main/stuck=0/NoColl=0；末态
  1,958,263 / 1.96296e-3，基线带内）。
- 作业 4746365 / 4746420 已按用户指示 scancel（分析完成后）。
- F5 已同步 bjm8 dlb-new（md5 双端一致，见同步清单）。

#### 重跑清单（下一次）——已执行，验收见 §4.8

1. 远端 `build.sh` 重建（源已含 F1–F5）→ 重跑 1200w-rcj。
2. 预期：DLB #1 整步耗时从 2617 s 大幅下降（F1 findIndex、F5 hashset
   两个 O(N²) 项均消除，剩余大头应为 ParMETIS 本体 + 迁移 +
   occupancy）；`replicatedMeshDLBProfile true` 可拆分验证。
3. 判读要点：日志应迅速通过 `owner map validated`（不再有分钟级
   停滞）；若个别 rank 仍慢，用同款 gdb 法直接定位剩余热点。

### 4.8 超算重跑验收（作业 4746677，2026-09-20 19:00–）

- 远端 `build.sh` 重建 exit=0（libdsmcFoam+.so 19:00；slurm 脚本已移除
  临时 `I_MPI_COLL_EXTERNAL=0`，恢复 hcoll 默认）。
- 作业 4746677（m12o48，m4cm0602-0604，与 2617 s 基线同节点）。

| 事件 | 修复前（4738842） | **修复后（4746677）** |
|---|---:|---:|
| DLB #1（step 50，触发失衡 5.41） | **2617 s wall**，其后每步仍待验证 | **~30 s wall**（ClockTime 172→202，含 ParMETIS + remap 省 661 万 cells + 635 万 cells 迁移 + occupancy 重建）≈ **87×** |
| DLB #2（step 100，触发失衡 1.503） | 挂死 >17 min 被杀 | **~15 s wall**，remap saved 12,952,343/12,958,175 cells（仅 5832 变更） |
| `owner map validated` | 分钟级停滞/挂死 | 立即通过（两次 DLB 均秒级） |
| 每迭代 wall（重平衡后） | ~1.08–1.1 s | ~1.0–1.2 s（正常） |

- 作业状态：iteration 170 / ClockTime 345 s 正常推进，目标 40,000 步
  （预计 ~12–13 h）。
- 观察项（不阻塞生产）：DLB 触发基于 `loadImbalance > 1.5`，step 100
  触发时失衡 1.503 贴线；若后续每 50 步都触发，单次仅 15–30 s 可接受，
  但未来长跑可考虑 `replicatedMeshDLBImbalanceThreshold 1.8`（todo
  遗留项，见 §1.4）。

**结论：F1–F5 全部生效，DLB 在 1300 万 cell / 1.5 亿粒子 case 上从
"挂死/2617 s"修复为秒级完成；作业 4746677 作为生产运行继续。**

## 5. 设计：DLB 触发/执行优化（back-off + wall 触发指标）— 2026-09-20【已实施并验证，实施记录见 §7.1/§7.4/§7.5】

### 5.1 动机（4746677 实测）

- 作业 4746677 前 4870 步触发 **96 次 rebalance（每 50 步检查全中）**，
  而各次 remap 显示划分早已稳定（step 100：1295.8 万 cell 仅变 5832），
  loadImbalance 稳定在 1.82–1.88 —— **rebalance 无收益的直接信号**，
  但 ParMETIS 每次 ~15–20 s 照付，累计 ~30 分钟（约 6% 已运行时间）。
- 该 case 的失衡属"粒子数已平衡、CPU 负载失衡"类（collision ∝ 密度²、
  per-particle move 成本差异、owner-gated inflow），ParMETIS 的粒子数
  权重压不住它（6 月报告 §5 已记录的客观局限），继续重试没有意义。

### 5.2 现有机制盘点（为什么还需要本设计）

| 已有机制 | 覆盖 | 缺口 |
|---|---|---|
| `nChanged > 0` 守卫（:1918/2392/2429/2451） | 归属无变化时跳过重建/迁移/后处理 | **ParMETIS 本体 + 集体通信照跑**（当前每次成本主体） |
| `replicatedMeshDLBParticleGate`（:2092） | 粒子失衡低于阈值时跳过全局检查 | 省检查不省 rebalance；对"粒子已平衡的 CPU 失衡"无效；当前 case 已关闭 |
| `replicatedMeshDLBMinRemainingSteps`（:2321） | 尾部护栏（剩余步数不回本则跳过） | 只覆盖 run 结尾 |
| adaptive alpha worsen 检测（:2563） | 检测上次 rebalance 是否更差 | 只反馈 alpha 自适应，不跳过后续 rebalance |

缺的正是用户设想的反馈环：**"上次 rebalance 无收益（失衡未改善）或
归属几乎无变化 → 退避跳过后续 ParMETIS"**。

### 5.3 设计

**状态记录**（每次实际执行的 rebalance 后更新，各 rank 天然一致，
无需新增通信）：

```text
lastExecChangedRatio_ = nChanged / nCells      // remap 后归属变化比例
lastExecImbalance_    = 触发时的 loadImbalance  // 全局 max/min（已有）
futileCount_          // 连续无收益计数
backOffRemaining_     // 剩余退避检查次数
```

**退避判定**（在每次 check 得到 loadImbalance、threshold 触发成立之后）：

```text
futile := (lastExecChangedRatio_ < eps)            // 归属几乎无变化
       && (loadImbalance >= lastExecImbalance_ - tol)  // 失衡未见改善
futileCount_ = futile ? futileCount_+1 : 0

if (backOffRemaining_ > 0 && loadImbalance < escalate)
    打印 "Phase C auto DLB skipped by back-off (remaining=M,
          imbalance=..., lastChangedRatio=...)"
    backOffRemaining_--; triggered = false; return;   // 跳过 ParMETIS
else if (futileCount_ >= K)
    backOffRemaining_ = M;                            // 进入退避
```

**退出/安全阀**：`loadImbalance >= escalate`（默认 threshold×1.5，
即显著恶化）时无视退避、立即执行——保证物理剧变（如激波形成、
大规模入流）不会被退避压制。

**一致性与正确性**：`nChanged` 来自 allgatherv 后的 fullPart，所有
rank 相同；`loadImbalance` 来自 gatherLoadExtrema 全局归约；计数与
判定在各 rank 上确定性重复，不引入新的集体点，无 mismatch 风险。
本设计只改变"何时调用 ParMETIS"，不改变 ParMETIS/remap/迁移语义。

**controlDict 开关**（默认关闭，A/B 用）：

```text
replicatedMeshDLBBackOff            false;  // 总开关
replicatedMeshDLBBackOffChangedRatio 1e-3;  // eps：归属变化比例下限
replicatedMeshDLBBackOffTol          0.05;  // 失衡改善容忍带
replicatedMeshDLBBackOffFutileRuns   2;     // K：连续无收益次数
replicatedMeshDLBBackOffSkipChecks   6;     // M：退避跳过的检查次数
                                            //   （check 每 50 步 → 300 步）
replicatedMeshDLBBackOffEscalate     1.5;   // 安全阀倍数（×threshold）
```

**改动点**：`dsmcReplicatedMesh.H` 新增 4 个状态成员 +
`initialize()` 读取开关；`autoRebalance()` 触发段加判定与计数；
`reassignByParMetisAdaptiveRepart()` 返回值处记录
`lastExecChangedRatio_`（nChanged 已有返回值，无额外计算）。

### 5.4 与现有旋钮的关系

- 与 `replicatedMeshDLBImbalanceThreshold 2.0` 互补：抬高阈值是静态
  手调，可能错过真实恶化；back-off 是自适应版本（无收益就退避、
  恶化即恢复），两者可叠加。
- 与 remap 正交：remap 已把"无收益 rebalance"的**迁移**成本降到零，
  back-off 进一步把 **ParMETIS 本体**成本降到零。

### 5.5 预期收益与验收标准

预期（以 4746677 前 4870 步为参照）：96 次 → 预计 2–4 次
（初始失衡 + 物理剧变点），节省 ~25–30 分钟/12 h 量级 run（~4–6%），
长步数 run 收益线性放大。

验收：

1. **功能**：退避生效时日志打印跳过原因与剩余计数；安全阀触发
   （失衡 ≥ escalate）时立即执行；`Total Iterations`/`End main` 正常。
2. **性能**：同 case 同配置 A/B（back-off off/on），DLB 次数显著
   下降（96 → 个位数）、总 wall 不升、每迭代 wall 不变。
3. **正确性**：末态 particles/energy/碰撞数落在基线噪声带；
   无 collective mismatch / rank 分歧。
4. **物理响应**：构造失衡剧变场景（或用 4746677 的 step 50 初始失衡
   5.41 段）验证安全阀能突破退避及时 rebalance。

### 5.6 设计二：触发指标改为每步完整 walltime（用户指定口径）

#### 现状与问题

- `dsmcCloud.C:2658-2661` 以 **CPU 时间**喂触发指标：
  `addEvolveTime(elapsedCpuTime() - moveAndCollideCpuStart)` →
  `evolveStepTime_` → `productiveTime_`（dsmcReplicatedMesh.C:1897-1904）
  → `loadImbalance`（max/min）与 SAR 趋势。
- 缺陷：① `elapsedCpuTime` 计入全部 48 线程忙等/自旋，读数与真实
  负载偏离；② 口径只覆盖 move+collide 段，migration 本地段、
  occupancy、fields 等**每步真实工作不进指标**；③ CPU−wall 混搭
  相减，不自洽。

#### 主口径定义（按用户指定）

```text
stepFullWallTime(rank) = 一次完整迭代的 walltime
                         （步起点 → 步终点，含 move/collision/
                          occupancy/migration/fields/info 全部段）
loadImbalance = max_rank(stepFullWallTime) / min_rank(stepFullWallTime)
```

实现载体：主循环 `dsmcFoam+.C` 已有 per-rank 迭代计时
（`loopWallStart`/`loopWall`，:274-275/514）——改为**无条件**记录每步
δ 并传入 cloud（`setLastStepWallTime()`），`autoRebalance` 读上一步
值（滞后一步，对 50 步窗口无影响）。每次 check 打印 per-rank
`stepFullWallTime` 全列——直接可视化"为什么有的进程 CPU 只有 2000%"。

#### 必须写明的口径风险（设计一级项）：集合点拉平效应

`migrateInterval 1` 下每步存在迁移集合点（size exchange Alltoall、
updateParticleCounts Allgather；delayedReceive 时同步点移到主循环
migrateFinish，效果相同）。轻载 rank 的等待发生在**步内**：

```text
rank_i 步 walltime = 自身集合点前工作 + 等待慢 rank + 自身集合点后工作
```

稳态下所有 rank 的步 walltime 被集合点收敛到同一周期，
**理论 max/min ≈ 1.0**，失衡信号只在负载演化/暂态（DLB 后、激波
发展期）短暂出现。推论：以完整 walltime 做触发指标，阈值触发可能
长期不激活（等效抑制失衡触发）——这是口径的物理属性，不是实现
缺陷；是否如此必须由 A/B 实测判定。

#### 并列备选口径（同时实现、同时打印，A/B 定夺）

```text
stepWorkWall(rank) = stepFullWallTime(rank)
                   − migrationSizeExchangeWall Δ
                   − migrationWaitWall Δ
                   − updateParticleCountsWall Δ
```

即每步完整 walltime 扣除纯等待段 = per-rank **工作墙钟**（本地
pack/localPrep/deserialize/collision/fields 全保留）。三个等待
累加器已存在，snapshot-δ 照抄 `productiveTime_` 模式。理论预期：
workWall 的 max/min 保留失衡信号（预计与 4746677 的 rank 忙闲分布
一致），而 stepFullWallTime ≈ 1.0——若实测如此，触发改用 workWall，
stepFullWallTime 保留为周期一致性诊断。

#### 实现要点与影响面

1. 主循环 per-rank 迭代 walltime 无条件记录 + `setLastStepWallTime()`
   传递（dsmcFoam+.C 与 dsmcCloud 各 ~5 行）；
2. 等待 δ 三处（复用已有累加器）；
3. `autoRebalance` 步级记账（:1897-1904）替换 stepEvolve 来源；
   `productiveTime_`/`gatherLoadExtrema` 结构不变；
4. 每次 check 打印两列 per-rank 值（`stepFullWallTime` 全列 +
   `stepWorkWall` max/min），`report()` 的 per-rank evolve 打印
   （:5163）同步改口径；
5. 开关：`replicatedMeshDLBTriggerWallTime`（默认 false，A/B 后定
   默认）；SAR 趋势（同源）随之变化。

#### 阈值语义与配套

- 历史 `threshold 1.5` 是 CPU 口径调的；新口径读数体系完全不同
  （stepFullWallTime 可能 ≈1.0，workWall 预计 1.8+）——必须与
  §5.3 back-off 配套上线并重校 threshold；
- 上线顺序：① off/on A/B → ② 观测两口径的 per-rank 分布与
  4746677 的 4700%/2000% 现象对照 → ③ 定触发口径与阈值 →
  ④ 叠加 back-off。

#### A/B 验证方案

1. **口径实证门（关键）**：实测 stepFullWallTime 的 max/min 是否被
   集合点拉平（验证本设计的风险预判）；workWall 与 rank 忙闲
   （4700%/2000% 现象、report per-rank 值）的相关性；
2. **正确性门**：指标为纯观测量，末态物理量带内；
3. **性能门**：DLB 触发次数与总 wall 不劣化（配合 back-off）；
4. **留档**：新旧口径完整日志各一份（触发行 + per-rank 打印 +
   report 段），供论文口径说明。

**工作量**：~40 行（主循环传递 5 行 + 等待 δ 3 处 + 触发段替换 +
两列打印 + 开关），复用已有累加器与 snapshot-δ 模式，风险低。

## 6. 100wcell-5000wparticle 性能分析与三项优化设计 — 2026-09-20【优化一判定无效（6.2 方案A）；优化二/三已实施并验证 → §7/§8】

### 6.1 基线（作业 4633878，m24o16/384 核 2 节点，F1–F5 前二进制）

- 规模：1,079,807 cells；粒子 50M 初始化 → 入流充填 → 稳态 79.29M；
  20,000 步，作业 2:57:16，求解循环 wall 9867 s（**0.49 s/步**，
  稳态 ~0.41 s/步）。正确性门全过；2 次写出 + reconstructPar +
  foamToEnsight（8.2 GB/121 s）+ 气动力系数（CD=1.709）全流程闭环。
- 相位 max（占 accounted 7534 s）：move 2881 / comm 2923 / post 2834 /
  build_occupancy 967 / coll 776 / io_info 472 / dlb_noncomm 65；
  other_residual 1433 s（17.5%，未归账）。
- **migration 分解（关键发现）**：wall 3200 s 中 **size exchange
  2619 s（82%）**，pack 340 / wait 208 / deserialize 121，
  updateParticleCounts 224 另计。size exchange 折合 0.13 s/步
  （每步 32%），而其数据量仅 24×24 int ⇒ **成本≈集合延迟 + 到达
  时间差**（per-step 负载失衡被吸收进该计时）——§5.6 拉平效应的
  直接实测。
- DLB 健康：27 触发 / 400 check（6.75%，无抖动，对比 1200w 的 96/96）；
  remap 末期 changedAfter=105/183 cells（0.01%）；dlb_noncomm 65 s
  （0.8%）——也解释了 100w 从不"挂死"而 1200w 会（O(N²) 项随规模
  放大）。

### 6.2 优化一：migration size exchange（2619 s，最大单项）——判定无效

> 实施与判定：方案 A（NoAlltoall）A/B +45% 更差，见本节中部
> 「方案 A 实测结果」小节；设计文字保留作档案。


**根因定性**：每步 Alltoall 交换 24×24 个 int，数据量可忽略；
2619 s = 集体延迟 + 各 rank 到达集合点的时间差（失衡等待被记账到
size exchange 名下）。

**方案 A（零代码，先行 A/B）——`replicatedMeshNoAlltoall true`**：

- 该路径已完整实现：按 peer 预测接收容量
  （`perPeerRecvCapacity_ = prevRecvSizes_×2 + minBuf`，dsmcReplicatedMesh.C:4346-4353），
  64 位 P2P header 携带实际字节数，`prevRecvSizes_` 每步自适应更新
  （:4379），完全消除每步 Alltoall 集合点。
- 历史判定：8-rank zb 上劣于 Alltoall（6 月报告 §7.1）——但当时
  无 2619 s 量级的 skew 可消除；24 rank + 失衡 1.5× 的前提不同，
  **必须重测**。
- A/B：同 case 各 3000 步，对比 size-exchange max、每步 wall、
  迁移正确性（粒子数守恒）。
- 预期：若失配率低（`prevRecvSizes×2` 对相邻步粒子数变化足够），
  2619 s 中的等待部分消失，size-exchange 项趋近纯通信毫秒级。

**方案 B（若 A 失配率高）**：预测改进——滑动平均（最近 N 步）替代
单步记忆 + per-peer 步数记忆；失配逃逸路径（二次接收）单列计时。

**方案 C（结构性，记录不动手）**：size exchange 与 collide 重叠——
不可行：reactions 在 collide 中增删 parcel，迁移分类必须 post-collide。
与 buildCellOccupancy 重叠需重排 migration/occupancy 顺序，风险大，
仅在 A/B 后仍有大头时再议。

**风险**：预测容量 ×2 的接收缓冲内存翻倍（当前每 peer 数 MB 量级，
可忽略）；NoAlltoall 关闭 shared 集合点后，失衡等待转移到
updateParticleCounts Allgather（224 s）——预期净收益仍为正。

#### 方案 A 实测结果（2026-09-21，dlbopt 1200w m12o48，200 步）——**判定无效**

作业 4749119（NoAlltoall true + wall/backoff）vs 4748993（Alltoall +
wall/backoff），其余配置相同；正确性门均过（End main/NoColl=0，
2 次 rebalance）：

| migration 子段（200 步 max） | Alltoall | NoAlltoall |
|---|---:|---:|
| migration wall | 99.0 s | 127.8 s（+29%） |
| **size exchange** | **73.1 s** | **106.2 s（+45%）** |
| pack / wait / deserialize | 7.5 / 13.4 / 11.5 | 8.8 / 14.0 / 11.8 |

每 10 迭代步 wall 相当（10–13 s vs 10–14 s）。**NoAlltoall 更慢**，
与 6 月 8-rank 结论一致。机制解释（本批数据支持）：size-exchange
成本的主体是**同步等待（到达时间差），不是集合算法**——NoAlltoall
只把 Alltoall 换成 2×(n−1) 条 P2P header 消息，同步语义不变
（接收方在发送方完成 move+pack 前无法得知大小），延迟反而更高。
**结论：§6.2 的"换交换算法"路线（A/B/C）整体受制于 skew-bound
本质，全部关闭**；每步迁移成本的进一步压缩只能来自更均衡的负载
（work 口径 1.07–1.36×，已由 wall 指标诚实度量）或减小每步迁移
数据量——size-exchange 等待应视为同步成本基线，不再是优化目标。

### 6.3 优化二：post fields 计算（2834 s / 38%）

**现状核实**：该 case 配置 6 个 dsmcVolFields entry（N2/O2/NO/N/O
+mixture，组分+mixture 全算）；post 子相位计时成员**已存在但从未
打印**（`profilePostReactions/FieldCalc/FieldWrite/Controllers/
Boundaries/BoundaryMeas/Clean*`，printProfileSummary 无输出点）——
先诊断后动手。

**第一步（诊断，必做，~20 行）**：report 增加 post 子相位打印 +
per-rank post 值。要回答两个问题：

1. 子相位分布：reactions / fieldCalc / fieldWrite / controllers
   (forceMoment 每步积分) / boundaries / boundaryMeas / clean 三项
   各占多少；
2. per-rank 分布：DLB 按粒子数加权 → 各 rank owned **cells** 数差异
   大（远场 rank 拥有大量空 cell），而 fields 采样按 owned cells
   展开——post max 是否来自 cell 数失衡（若是，post 失衡与 §5.6
   的 workWall 指标将互相印证）。

**候选（按诊断数据决策，不预设）**：

- a. fieldCalc 主导：确认共享 sample cache 对"5 组分 + mixture"配置
  的命中率（T10/T12 的共享 cache 本为此设计；若 6 个 entry 仍在重复
  遍历 parcel，则是一次遍历派生全部场的改造点）；
- b. controllers/forceMoment 主导：force/ moment 每步积分 → 增加
  采样间隔开关（物理上等价于力系数时间平均）；
- c. boundaries/boundaryMeas 主导：sparse clean 已有（T13），查剩余
  项；clean 三项的 OMP 并行化曾被证伪（D3/D4，OMP 报告 §18.6），
  不重复；
- d. 物理层面：若论文只需 mixture 宏观场，组分 field entry 可裁剪
  （用户决策，非代码问题）。

### 6.4 优化三：buildCellOccupancy（967 s / 12.8%）

**现状核实**：单一 wall/cpu 计时（dsmcCloud.C:855/964），**无子相位**
——OMP 报告 §14.1 "先做 subphase profile 再改" 的遗留事项至今未做。

**第一步（诊断，必做）**：子相位计时——ordered path 命中率、fallback
gather 次数、appended parcel 数、per-thread active cells、
reset/count/reduce/fill 分项、DynamicList resize 行为。

**候选（诊断后决策）**：

- a. **增量更新**（结构性候选）：`migrateInterval 1` 下 migration
  精确知道哪些 cell 的粒子集合变化（每 peer 发送/接收的 cell 明单）
  → 只 patch 这些 cell。**主要工程风险**：flat ordered 视图的
  `occupancyCellOffsets_` 是全局前缀和，局部 patch 无法保持——
  可行形态是 per-cell `cellOccupancy_` 增量维护 + flat 视图仅按需
  重建（或容忍 flat 视图滞后一步，需验证消费方容忍度）。设计论证
  通过后才立项；
- b. 若子相位显示 OMP 并行度/内存带宽瓶颈：当前 gathered thread
  offsets 路径已 OMP（48/16 线程），重点是 NUMA/带宽而非加线程；
- c. v2506 gap（zb omp8 上 7.9 vs 3.5 s）对照：v2506 无 flat 视图
  双重维护成本，若诊断证实，考虑合并两套视图。

### 6.5 实施顺序与统一验收

顺序（诊断先行，零代码项优先）：

```text
① 6.2 方案 A（NoAlltoall A/B，零代码，立即）
② 6.3 诊断打印（~20 行）→ 6.4 子相位计时（~30 行）——可同批
③ 依据 ②数据立项实施 6.3/6.4 的候选
④ 全部生效后重测 20000 步基线，更新本节表格
```

统一验收：correctness 门不变（Total Iterations/End main/stuck=0/
NoColl=0/末态带内）；A/B 需 3000+ 步独立结果目录；每项收益以
相位 max 变化归因（严禁只看总 wall）；与 §5（back-off + wall 触发
指标）联动——size exchange 的等待吸收实测是该设计的直接论据。

## 7. §5/§6 实施与冒烟验证 — 2026-09-21

### 7.1 实施内容（5 文件，均已同步 bjm8 dlb-new）

| 项 | 文件 | 内容 |
|---|---|---|
| §5.6 wall 指标 | dsmcReplicatedMesh.H/C + dsmcCloud.H/C + dsmcFoam+.C | 主循环每迭代完整 wall 无条件记录 → `setLastStepWallTime()` 传入；autoRebalance 记账双口径（开关 `replicatedMeshDLBTriggerWallTime`，默认 false）；等待 δ（sizeExchange+wait+updateParticleCounts）滞后一步配对；**work = 完整 wall − 集体等待**；check 时 per-rank step wall/work 全列打印；report 增加累计口径 |
| §5.3 back-off | dsmcReplicatedMesh.H/C | `replicatedMeshDLBBackOff`（默认 false）+ eps/tol/futileRuns/skipChecks/escalate 五参数；退避门置于 minRemainingSteps 之后、全局一致性 Allreduce 之前（全输入跨 rank 一致，无新集体点）；执行后记录 lastExecChangedRatio/Imbalance 并复位；安全阀 escalate=1.5×threshold |
| §6.3 post 子相位 | 无需代码 | printProfileSummary 已含全部 post 子相位 + per-rank 明细（本次核对确认） |
| §6.4 build occupancy 子相位 | dsmcCloud.C | reset/count/reduce/fill 四段计时 + ordered/fallback 计数（文件级累加器，OMP 路径；printProfileSummary MAX 归约打印） |
| 调试标记 | dsmcFoam+.C/dsmcCloud.C | 已加已清（定位挂死用） |

### 7.2 实施中发生并修复的挂死（重要教训）

本地 300 步 ON 回归出现确定性退出挂死（End stage 后、Total
Iterations 前，rank 自旋）。定位过程：开关二分（wall ON 即挂，OFF
干净）→ **根因 = report() 中新加的 per-rank step wall Allreduce 位于
`if (myRank_ != 0) return;` 之后——只有 rank0 进集合点，其余 rank
提前返回并退出进程，rank0 永久等待**。修复：归约移到 early return
之前（与 allTimes 同段，全 rank 执行）。修复后 300+ON 回归 exit=0
全门通过。教训：**report()/printProfileSummary 类 rank0 早退函数中
新增任何集体通信都必须放在早退之前**（与 :5432 的既有模式一致）。
另：本地 10 步迷你运行（endTime=10 步）在 F5 二进制上同样存在退出
挂死（预存在、与本批改动无关，未深究——短窗测试请用 ≥300 步）。

### 7.3 本地回归（zb mpi4omp2，300 步）

- 开关 OFF：exit=0，全门过（1,957,847/1.96191e-3），build occupancy
  子相位首次有数据（count 1.71 / fill 1.73 s）。
- 开关 ON（wall+backoff）：exit=0，全门过（1,957,944/1.96184e-3），
  2 次 DLB；per-rank step wall 打印 6 次；**拉平效应实证**：
  cumulative max/min 仅 1.0026（理论预判 ✓）。
- F5 基线对照：300+OFF 干净、300+ON（修复前）挂、修复后干净——
  修复有效性直接验证。

### 7.4 超算冒烟（4748993，1200w-rcj-dlbopt，m12o48，200 步）

- 重建（.so 01:01）→ 提交 → **End main、NoColl=0、stuck=0 全门过**。
- back-off 配置打印生效；rebalance #1（step 50，wall 口径失衡
  2.806 vs CPU 口径 3.98–5.41）；lastChangedRatio=0.469（真实重平衡，
  futileCount=0）。
- **per-rank step wall/work 数据（本次最有价值产出）**：

| check | step wall max/min | step work min–max | 解读 |
|---|---:|---|---|
| step 50（DLB#1 前） | 1.14 | 0.28–1.34（**4.9×**） | 初始失衡（METIS 均分 cells，粒子按密度集中） |
| step 100（DLB#1 后瞬态） | 1.10 | **0–1.08** | 重分布瞬态：快 rank 的步内等待被 size exchange 吸收，work≈0 |
| step 150 | 1.09 | 0.80–1.09（1.36×） | 收敛中 |
| step 200 | 1.06 | 1.02–1.09（**1.07×**） | 均衡 ✓ |

  - **step wall 全程拉平（1.06–1.14）——§5.6 拉平预判实证**；
  - **work 口径完整捕获重平衡动态**，且与触发决策自洽
    （step 150/200 的 work max/min 1.36/1.07 < 1.5 → 未触发 ✓）；
  - 这组数据正是 4746677 上"有的进程 CPU 只有 2000%"现象的
    wall 版解释：轻载 rank 的步内等待被记账进 size exchange
    （100wcell 分析中 2619 s 的同源现象）。
- 作业在 End main 后被 slurm 尾部 reconstructPar 拖住（200 步窗
  无写出，latestTime 无新目录）——求解数据完整，可 scancel。

### 7.5 A/B 同比对比（dlbopt 4748993 vs 基线 4746677 前 200 步，同 m12o48 同物理）

两者前 40 步逐 stamp 一致（45.5/45.4/46.4-49.0/47.4-48.0 s，噪声带内）；
DLB#1 均为 30 s wall。分歧从 step 100 的 check 开始——**wall/work 口径
正确避开了基线的两次无效 rebalance**：

| 事件 | 基线 4746677 | dlbopt（wall+backoff） |
|---|---:|---|
| iter 50（DLB#1，真实） | 30 s | 30 s（相同） |
| iter 100（基线 futile #2，5832 cells） | **57 s** | **14 s**（work 1.36×<1.5 → 不触发） |
| iter 200（基线 futile #3） | **137 s** | **12 s**（work 1.07× → 不触发） |
| ClockTime @ iter 200 | 396 s | **374 s（领先 22 s，~5.5%）** |
| 每 10 迭代 wall（100–200 段） | ~13.7 s | ~11.7 s（**~15%**） |

- 提升 = 被避开的无效果平衡机制（每 50 步 ~60–125 秒的无效 ParMETIS +
  迁移重建），**结构性收益、随步数累积放大**；
- 基线 4746677 运行 6h16m 时已累计 **260 次 rebalance**（每 50 步全中
  的抖动模式持续）——按每次 ~15–30 s 估算，全程浪费 1–2 h 量级；
- 注意事项：单次对比含 RNG 轨迹噪声；确定性成分是被避开的 DLB 步。

### 7.6 profileDetail 诊断（作业 4749270，dlbopt 200 步 m12o48，全门过）

§6.3/§6.4 的"诊断先行"步骤完成，结论：

1. **post 定位闭环**：field calculate = 90.5/101.6 s（**89%**），clean
   11%（tracking 8.5），write/controllers/boundaries/reactions ≈0，
   子相位和闭合（残差 0.004）。候选 (a) 确认：fieldCalc 是唯一值得
   动的 post 项（6-entry 配置 0.45 s/步）。
2. **build_occupancy 子相位闭合**（44.2 ≈ 44.0 s）：reduce 16.6 + fill
   17.3 = **77%**（线程合并 + flat 散射，内存流量型），reset 5.5、
   count 4.8；**ordered 命中率 99.5%**（402/2）——stage25 机制健康，
   v2506 gap 定位到 reduce+fill。
3. **move+collide 窗口含常规迁移**：move 62.4 + coll 7.4 + build 44.0
   = 113.8 vs moveAndCollide max 210.8 → ~97 s 常规迁移计入该窗口
   （0.49 s/步），与 size exchange 113 s（0.57 s/步）互证。
4. **碰撞候选热点（新量化）**：rank0 finalCandidates **2.54M** vs 其它
   rank 57–124K（~44×）；rank0 coll CPU 7.4 s vs 1.2–2.3 s（6.15×）。
   rank0 粒子最少（10.4M）却候选最多 → **碰撞负载 ∝ 局部密度结构，
   与粒子数弱相关**——粒子数权重 ParMETIS 原理上看不见它。
   collision-aware DLB（6 月报告 §10.2）的立项依据至此齐备。
5. **back-off 首次真实运作**：step 150 rebalance #2 changedRatio=
   0.00052 < eps → 记 futile #1（机制正确）；step 200 触发被
   near-end guard 拦截（窗口边界效应），arm→skip 完整周期需 ≥500 步。
6. **wall/work 指标自洽**：触发读数 1.559@150 与 per-step work 分布
   （~1.47×）一致；step wall 拉平 1.05–1.15 ✓。

### 7.7 待办（下一步）——已由 §8 承接并完成

1. dlbopt 长窗 A/B：endTime 0.01（2000 步，~40 次 check）量化
   back-off 对 96-次触发模式的抑制收益（对照 4746677 前 2000 步）。
2. §6.2 方案 A：`replicatedMeshNoAlltoall true` A/B（零代码）。
3. §6.3/6.4 诊断数据采集（20000 步 100wcell 复跑或 dlbopt 长窗）
   → 按数据立项 post/build_occupancy 候选。
4. wall 指标口径的阈值重校（当前 1.5 在 wall 口径下的触发行为
   需 2000 步数据确认）+ back-off 参数标定。

## 8. post/build 瓶颈优化计划 — 2026-09-21（基于 4749270 profileDetail）

用户判定：post 与 build 是当前真正瓶颈（合计占关键路径 ~60%：
post max 0.51 s/步 + build max 0.22 s/步，迭代 ~1.1–1.3 s/步）。
profileDetail 的实测分解给出了两个低风险高收益的修复点。

### 8.1 实测分解（4749270，200 步 m12o48）

**post fieldCalc 90.5 s 分解**（既有 BuildProfile 设施输出）：

| 子项 | 时间 | 说明 |
|---|---:|---|
| **cache allocate** | **25.7 s（128 ms/步）** | 共享 cache 数组每步重新分配+清零 |
| parcel accumulate | 12.1 s | 真实采样遍历 |
| field combine + cell reduction（×6 entry） | ~14 s | 派生 |
| boundary accumulation（×6 entry） | ~20 s | 组分边界量 |
| 后续 5 entry 的 cache 复用 | 0.0001 s | 共享机制健康 ✓ |

**build 44.0 s 分解**：reduce 16.6 + fill 17.3 = 77%（含每步
`totalCounts` 52MB 清零 + `nextCellOffsets` 52MB 拷贝 +
`occupancyCellOffsets_` 52MB——**~156 MB/步 的分配churn**）；
ordered 命中率 99.5%（reset 已稀疏化 ✓）。

### 8.2 OMP 覆盖核查（2026-09-21，代码行号实证）——post/build 远未充分利用 OMP

**build 44 s：约 70% 时间在串行段**

| 子相位 | 时间 [s] | OMP 状态 |
|---|---:|---|
| count | 4.8 | ✓ OMP（`#pragma omp parallel`，:649 段） |
| fill 散射段 | ~7 | ✓ OMP（:805 段） |
| **reset** | **5.5** | **✗ 串行**（:608-647 对 48 个线程列表逐个清零） |
| **reduce**（含 totalCounts 52MB 分配清零） | **16.6** | **✗ 串行**（:708-719 合并循环） |
| **fill 前缀段**（13M 前缀和 :781-785 + 52MB 拷贝 :788 + active-cell 更新 :790-803） | **~10** | **✗ 串行** |

**post fieldCalc 90.5 s：约 50% 时间在串行段**

| 子项 | 时间 [s] | OMP 状态 |
|---|---:|---|
| **cache allocate** | **25.7** | **✗ 串行**（`initScalarFields` 的 setSize+双重清零——纯 memset 型，且持久化可整体消除） |
| parcel accumulate | 12.1 | ✓ OMP（:1495/:1667，activeCells 门控） |
| combine（:3561）/ reduction（:3807） | ~14 | ✓ OMP（`if (useOpenMPSampling)` 门，本 run 门开 ✓） |
| **boundary accumulation** | **~20** | **✗ 串行**（dsmcVolFields.C 全部 11 处 omp 站点无一覆盖） |

**合并：两相 ~134 s 中 ~76 s（57%）为串行**，且均为 memset/合并/前缀型
易并行或可消除操作。这也与 m8-new 报告的历史观察吻合（post 扩展性
最差、OMP64 下反而变慢）。

### 8.3 修订收益表（持久化 + 串行段 OMP 化叠加）

| 项 | 现状 [s] | 措施 | 预期 [s] |
|---|---:|---|---:|
| cache allocate | 25.7 | 持久化 + touchedCells 稀疏复位（直接消除） | ~0 |
| build reduce | 16.6 | totalCounts 持久化 + merge 按 cell 并行（每 cell 求 48 线程和，无冲突写） | ~5 |
| build fill 前缀段 | ~10 | nextCellOffsets 持久化 + 分块并行前缀 | ~4 |
| build reset | 5.5 | 按线程分块并行（各线程清自己的列表，无共享写） | ~1 |
| boundary accumulation | ~20 | 按 patch/cell OMP | ~5 |
| **合计** | **~78** | | **~15（−80%）** |

**叠加预期：迭代关键路径 1.1–1.3 s/步 → ~0.75–0.9 s/步（−30% 上下）**，
约为原 P0/P1（仅持久化，−15%）的两倍。全部为无冲突写或前缀和型
低风险 OMP 模式。

### 8.4 P0：共享 cache 持久化 + 稀疏复位

现状（dsmcVolFields.C:714-760）：`allocateAllFields` 每步对每个数组
执行 `setSize(nCells, 0.0)` **加** `fields[typei] = 0.0`（双重清零），
fresh 分配又触发内核页清零。规模：~10 个数组 × 5 组分 × 13M cells
≈ 2.6–4 GB/步 的分配+清零（实测 128 ms/步）。

修复（模式 = T13 sparse clean + 结构体已有 `touchedCells`）：

1. 数组容量持久化：`size()==nCells` 时跳过 setSize；
2. 删除冗余的 `fields[typei] = 0.0` 第二遍；
3. 累计前的清零改为**按 touchedCells 稀疏复位**（上一步被写的 cell
   才需要清）——touchedCells 基建已在结构体中；
4. 验证：gate = 4749270 同窗复跑，cache allocate 25.7 → ≈0，
   sample accumulation 不变，末态场逐 cell 一致（mixture 修复的
   验证模式可复用）。

### 8.5 P1：build 串行段持久化 + 并行化

1. **reduce**：`totalCounts` 提升为持久成员（容量保留）；合并循环改为
   按 cell 并行（cellI 的 48 线程值求和，各线程写不同 cell 槽无冲突，
   或按线程外层、cell 内层分块）；52MB 分配消失；
2. **fill 前缀段**：`nextCellOffsets`/`occupancyCellOffsets_` 持久化；
   13M 前缀和改分块并行（各块局部和 + 串行块间前缀 + 平移）；
   52MB 拷贝因持久化消失；
3. **reset**：按线程分块并行（threadI 列表只被 threadI 清，无冲突）；
4. 验证：子相位计时逐段前后对比（reset→~1、reduce→~5、fill→~10），
   occupancy 内容逐 cell 一致（旧实现 vs 新实现同步骤运行比对）。

### 8.6 P2：boundary accumulation OMP 化

~20 s 串行（6 entry × 每 entry 的 boundary 循环）：按 patch/face 维度
并行（各 patch 写各自的测量槽，无冲突）。位于 post 第二大子项，
模式与 field sampling 的 OMP 化一致。

### 8.7 残余问题（暂不动）

- post/build 的 per-rank 失衡（post 1.28×、build 3.7×——rank0 低密度
  大 footprint 区域的 locality 成本）：粒子数权重 DLB 看不见，属
  footprint/collision-aware DLB 范畴（与 §7.6 候选热点同根）；
- per-entry 派生 ~14 s（combine/reduction）：已 OMP，剩余为访存本质；
- back-off 长窗验证与阈值标定：200 步 A/B 见 §7.5；40000 步
  生产 run（19h12m 完赛）已包含 back-off 长窗行为，无异常
  报告——视为完成。

### 8.8 P0/P1/P2 实施与验证记录（2026-09-21，全部完成）

**P0（cache 持久化 + 稀疏复位）**：initScalar/VectorFields 容量保留、
删除每步全量 operator= 清零；build() 的 !validFor 路径改为
"allocateFields（容量保留）+ touchedCells=当前 activeCells +
resetFields（稀疏）"。正确性验证：本地强制写出（writeInterval 5e-6，
3 个写出时刻）+ mixture == Σ组分 **逐 cell 精确相等**（0.000e+00）
+ 瞬时场平稳性（N2 t3/t1 均值 1.044，无累积增长）。

**P1（build 串行段）**：reset 按线程分块并行；reduce 用持久
`occupancyTotalCounts_`（dsmcCloud.H 新成员）+ 并行清零 + 按线程
atomic merge（active lists 线程分离、cell 可跨线程 → atomic）；
fill 前缀和分块并行（块内和 → 块间扫描 → 平移加，3 段式）。
nextCellOffsets 52MB 拷贝与 cursor 更新按计划保留（收益/风险比
不足）。

**P2（boundary accumulation OMP）**：species×patch×face 累积循环按
sampled patch 并行（dynamic schedule）；写 per-(patch,face) 分离 ✓、
boundaryFlux 访问器为 const 纯读 ✓（已核验 boundaryMeasurementsI.H）。

**bjm8 200 步逐段归档（4749270 修复前 → 4749668/4749760/4749840
P0/P1/P2 后，同 dlbopt case 同 m12o48）**：

| 指标（200 步 max） | 修复前 | 修复后 | 变化 |
|---|---:|---:|---|
| cache allocate | 25.7 s | **0.45 s** | −98% |
| cache reset | 0 s | 6.2–11.5 s | 稀疏复位代价（新增，符合设计） |
| sample accumulation | 38.8 s | 13.4 s | −66%（持久页局部性） |
| **boundary accumulation** | **~20 s** | **2.31 s** | **−88%** |
| **post max** | **101.6 s** | **78.7 s** | **−23%** |
| **build max**（reset 0.6/count 4.8/reduce 7.5/fill 16.9） | **44.0 s** | **29.6 s** | **−33%** |
| move max（未触碰，对照） | 62.4 s | 62.8 s | 不变 ✓ |
| rank step wall（full iter, cumulative） | 281–297 s | 241–256 s | **−40 s（−14%）** |

每步关键路径 1.40 → 1.21 s（−14%），与 §8.3 的保守预期一致
（move 62 s 与迁移同步 ~0.5 s/步不随本批改动变化，稀释了全局占比；
build/post 内部的串行段本身 −70~80%）。

正确性门（每阶段 bjm8 200 步）：End main、NoColl=0、stuck=0 全过；
本地另有 mixture==Σ组分逐 cell 精确验证。作业号：P0=4749668、
P1=4749760、P2=4749840（同 dlbopt case）。

### 8.9 三项追加尝试（2026-09-21）：①完成 ②审计否决 ③评估后缓行

#### ① build-fill atomic-cursor（实施并验证 ✓，作业 4750177）

删除串行 cursor 更新循环（对全部 active cell 写 per-thread 游标，
14M 触点串行），scatter 改为对共享游标数组 `nextCellOffsets` 的
`#pragma omp atomic capture` fetch-add（游标初始 = 排他前缀）。
cell 内 slot 分配从"线程序"变为"到达序"——**统计等价、非逐位**
（occupancy 消费方为无序求和/独立轨迹，T13/Tier-2 同类先例）。

| 指标（200 步 max） | P1 后 | ①后 | 变化 |
|---|---:|---:|---|
| build fill | 15.76 s | **8.23 s** | **−48%** |
| **build max** | **26.02 s** | **19.10 s** | **−27%** |
| rank step wall | 221–236 s | 227–242 s | 噪声带内持平（build 非 step-wall 最大段） |

正确性门全过（End main/NoColl=0/stuck=0）。P0–① 累计：build
44.0 → 19.1 s（**−57%**）。

#### ② post-cache reset 移除（审计否决，不实施）

审计 build() 后的 per-cell 赋值语义（dsmcVolFields.C:1405-1440）：

1. 赋值循环遍历 **`activeTypes`**（cell 上实际存在的组分，
   :1408-1410）——cell 上缺席的组分**不被赋值**，其 cache 值依赖
   resetFields 清零；移除 reset 会让缺席组分的场出现幻影密度
   （派生读取全部 active cell 的各组分场）；
2. **totals 是 `+=` 累计**（totalDsmcN[celli] += localN，
   :1425-1430），不是赋值——依赖清零后的 cache（T12 错误类同源）。

**结论：reset load-bearing，保留。** 若要消除其成本（6.2–11.5 s），
需将赋值循环改为全 type 覆盖（缺席 type 显式赋 0）——把 reset 的
清零量换成交付量的加宽，净收益不明确且触及 T12 类正确性边界，
不做。

#### ③ post-entry fusion（评估后缓行，需独立工作包）

6 个 field entry 为独立 dsmcField 对象（虚 calculateField），entry
fusion = 一遍 activeCells 遍历派生全部 entry（共享同一份 cache 读），
可省 ~2/3 的 combine/reduction 时间（~10 s 量级）。但需重构
dsmcFieldProperties 的 per-entry 对象设计（批处理 entry 组、统一
派生上下文）——架构级改动，风险在 per-entry 场语义（typeIds 子集、
resetAtOutput 时序、writeFields 白名单的独立门控）。**缓行**：登记
为独立设计任务（优先级低于 collision/footprint-aware DLB）。

#### P0–① 累计成果（同 200 步窗对比 4749270 修复前）

| 指标 | 修复前 | 现状 | 变化 |
|---|---:|---:|---|
| **build max** | **44.0 s** | **19.1 s** | **−57%** |
| **post max** | **101.6 s** | **75.9 s** | **−25%** |
| rank step wall | 281–297 s | 227–242 s | **−55 s（−19%）** |

**实施顺序：P0 → P1（reduce → fill 前缀 → reset）→ P2 → dlbopt 同窗
复跑归因（子相位逐段对比）→ 长窗 A/B。**

## 9. 迁移同步等待削减：① 权重方案重测 + ② 时序分析 — 2026-09-21

背景：负载分布分析（4750177）确认 work 失衡 ~1.4–1.5×（candidate
热点 22×、move 2×、build footprint 1.7×），迁移 size-exchange 等待
~0.5 s/步 为其直接影子。两项工作并行推进。

### 9.1 ① 双约束/adaptive-alpha 权重方案重测（实施：A/B）

**重测理由**：6 月判定"无效"的上下文与现在不同——当时（8 rank、
无热点量化）热点未被发现；现在 candidate 热点已量化（1.81M vs
80K，22×），而**双约束的第二权重 N(N-1) 恰好就是碰撞负载代理**，
adaptive-alpha 则动态调整 N^α 的平衡刚度——两者都直接瞄准已确诊
的失衡维度。

A/B 矩阵（dlbopt 200 步 m12o48，其余配置同 4750177 基线）：

| run | 配置 | 验证点 |
|---|---|---|
| 基线 | dual=false, adaptive=false | 4748993/4749668/4749760/4750177 已有 |
| B | `DLBDualConstraint true`（ubvec1=1.5） | 候选热点是否收敛（rank6 1.81M → ?）、coll max/min |
| C | `DLBAdaptiveAlpha true` | work 失衡与 candidates 分布 |
| D | B+C 组合（视 B/C 结果） | — |

判定指标：candidates max/min、coll max/min、per-rank step work
max/min、rank full max/min、正确性门。

### 9.2 ② post-move 时序分析（结论：无可重叠窗口，方案关闭）

取证（dsmcCloud.C:2594-2629）：每步迁移序列为
`migrateBegin(); migrateFinish();` **背靠背**（delayedReceive 的
异步机制在每步路径上没有可延迟的间隙）；时序为
move → 迁移 → collide。

关键约束：collide 只遍历 owned collision cells（occupancy 由迁移后
的 cloud 重建）——**跨 owner parcel 必须在 collide 前到达其 owner**，
否则该步碰撞不可见（migrateInterval=10 的错误机制）。因此：

1. 迁移点（move 后、collide 前）已是最小 skew 位置——若后移到
   fields 之后，fields 的 0.39 s/步（1.28× 失衡）会**叠加**进 skew；
2. begin 与 finish 之间无可重叠的本地工作（collide 需要迁移后的
   occupancy）；
3. 等待量 = move 相位 skew（≈0.15–0.3 s/步）+ 周期内 fields/post
   贡献——**唯一削减手段 = ①的负载均衡**；
4. 理论终态（ghost 副本供碰撞 + 真实迁移异步）需重构 DSMC 所有权
   模型（双计风险），超出范围，登记为远期方向。

**结论：② 无独立实施方案，其收益由 ① 承载。**

### 9.3 ①的实施与结果记录（2026-09-21，A/B 完成）

**run B（dual=true，4750453）——更差，确认关闭**：

| 指标（200 步） | 基线（4750177） | dual=true |
|---|---:|---:|
| rank step wall | 227–242 s | **245–260 s（+18 s）** |
| rank collision max/min | 5.33 | **6.34（更差）** |
| candidates max | 1.81M | 1.25M（略降但新增 r0 1.25M/r2 692K/r6 535K 三处次热点） |
| owned cells 分布 | 0.60–0.89M | **0.012–1.08M（畸变）** |
| 正确性 | — | End main/NoColl=0 ✓ |

**run C（adaptive=true，4750481）——短窗中性，确认关闭**：

| 指标（200 步） | 基线 | adaptive=true |
|---|---:|---:|
| rank step wall | 227–242 s | 229–244 s（噪声带内） |
| rank collision max/min | 5.33 | 5.48（同带） |
| adaptive 状态 | — | lastImbalance=1e15（仅 1 次 rebalance，从未激活） |

**综合判定（重要认知修正）**：

1. 两个权重方案均未收敛 candidate 热点（dual 后热点分裂为三处，
   adaptive 短窗未激活）——**热点是物理性的**（激波密集区的
   cell 密度结构），ParMETIS 只能移动 cell，移动密集 cell 必然
   破坏其它维度的分布。
2. **碰撞热点不是 wall-time 瓶颈**：coll max 5.8 s / 200 步 =
   0.029 s/步，其失衡的墙钟代价仅 ~0.02 s/步（且热点 rank6 的
   full 并非关键路径最大值）——collision-aware DLB 的预期收益
   被实测数据**降级为低优先级**。
3. **真正的关键路径失衡是 move 相位的 2× 分散**（move 31–61 s，
   粒子数仅差 1.34×→ per-parcel tracking 成本随区域变化 ~1.5×）；
   step-wall 最大的 rank 是 move-轻 rank（等待 move-重 rank 的
   迁移同步）。**下一步的正确诊断工具 = move 子相位剖分**
   （tracking/访存按区域分解），而非 DLB 权重改造。
4. 本轮 A/B 同时验证了 back-off 与 wall 指标在 dual/adaptive
   配置下行为正常（dual 的 3 次 rebalance、lastChangedRatio 7e-5
   均被正确记录）。

配置已恢复基线（dual=false, adaptive=false）。§6.5 的
collision-aware DLB 方向降级；move 子相位剖分升级为下一步首选。

## 10. Per-cell 成本模型 DLB 权重（计数器 × 窗口拟合系数）— 2026-09-21【已实施（10.8）→ bjm8 回归（10.8.6）→ 归一化修复（10.9）→ 证伪关闭（10.9.1）】

### 10.1 定位（与 §9.3 结论的对接）

§9.3 实测降级了 collision-aware DLB（热点墙钟代价仅 ~0.02 s/步），
并把真正的关键路径失衡定位为 **move 相位的 2× 分散**（move 31–61 s，
粒子数仅差 1.34× → per-parcel tracking 成本随区域变化 ~1.5×）。

本模型的定位据此修正：它**不是** collision-aware 的翻版——其最有价值
的项是 **move 项（moveItersPerCell）**，恰好直接捕获 §9.3 要求度量的
"per-parcel tracking 成本随区域变化"。模型精神：**不预设哪个驱动重要，
拟合系数由数据决定**——碰撞项若真的不重要，拟合会自动给它小系数。

### 10.2 最终形态（已与用户确认）

```
窗口内每步（×50 步累计）：
  build 时：N(c)        ← occupancy 计数，每 cell 每步（已有，零开销）
  collide 时：cand(c)   ← NTC 候选计数，每 cell 每步（已有 nCandidatesPerCell_）
  move 时：iters(c)     ← 跟踪迭代数，每 cell 每步（已有
                           moveItersPerCellCumulative_ 成员）
  → 累进 *_Cumulative_ 版本（窗口内累加，平均单步涨落）

窗口末（DLB check，每 50 步）：
  每 rank 实测分相 wall（profile 计时器现成，建议用 §5.6 work 口径）
  ↓ 最小二乘（12 rank = 12 方程）
  解出分相系数
  ↓
  权重：vwgt(c) = Σ_相 系数_相 × 计数(c)   （c = 本 rank owned cells，
  feed 给 ParMETIS vwgt）
```

**关键口径修正（与用户确认的表述的差异）**：per-cell 层面测量的是
**计数**（N/cand/iters），不是耗时——cell 级计时不可行也不必要；
耗时只在 **rank 层**作为拟合目标出现，cell 级成本是方程的推断值。

### 10.3 拟合设计

- **分相拟合**（move/coll/build/post 各自小回归）而非复合单方程：
  每相的物理驱动不同（move ~ N+iters、coll ~ cand、post/build ~
  N+active），12 方程对 2–3 系数，条件数好；vwgt = 分相权重之和。
- **共线性风险**：N 与 iters 可能强相关（iters ≈ N × crossing-rate）。
  处理：逐窗检查各列条件数；若共线则合并 N 与 iters 项（crossing-rate
  均匀时合并无损，不均匀时 iters 项优先保留——它才是 2× 失衡的来源）。
- **非负约束**：系数物理上非负（NNLS 或残差投影）。
- **拟合目标口径**：用 §5.6 work wall（扣除集体等待）——与触发指标
  同源，避免 CPU 自旋伪影污染系数。

### 10.4 前提条件（第一道门）

**moveItersPerCell 的生产可用性审计**：该计数器当前是否只在
moveDetailProfile 开启时累加？若是，需轻量化常开（每 parcel 一次
计数，预计开销 <1%，需实测）。这一步同时产出 move 子相位剖分的
基础数据——**与 §9.3 的"move 子相位剖分首选"合并为一个工作包**。

### 10.5 权重接线与冷启动

- vwgt[c] = scale × cost(c)，正值 clamp 复用现有 positiveWeight；
- 冷启动：第一窗沿用现有 N^α（无计数数据），第二窗起切 fitted；
- 与 §5.3 back-off 协同：新权重生效初期失衡读数波动由 back-off
  吸收；与 §5.6 wall 指标协同：拟合目标与触发指标同口径。

### 10.6 验证方案

1. **拟合优度门**：逐窗 R²、rank 级残差 ≤10%（目标 5%）；
2. **A/B**（dlbopt ≥2000 步）：weight=particle vs fitted——
   per-rank step work max/min（1.4–1.5 → 目标 ≤1.2）、迁移同步等待
   （rank step wall 与 move+collide 之差）、step wall、正确性门；
3. **机理确认**：fitted 系数的 move 项应显著非零（验证 §9.3 的
   per-parcel 成本随区域变化假设）；
4. **热点跟随**：候选/cell 密度热点（§7.6 rank6 1.81M）在新权重下
   的 cell 应被 ParMETIS 主动减载。

### 10.7 工作量与风险

| 阶段 | 内容 | 量级 |
|---|---|---|
| 前提审计 | moveItersPerCell 生产开销 + 候选累计版 | 0.5 天 |
| 拟合例程 | 分相 NNLS/简化回归 + 接线 | 1 天 |
| A/B 与标定 | dlbopt 长窗 + 系数审视 | 1 天 |

风险：①共线性（有预案）；②计数器常开的热路径开销（需实测，
超预期则退回采样版计数器）；③新权重的 DLB 行为漂移（back-off
安全网兜底）。

### 10.8 实施记录（2026-09-21，file-static 重新实施版）

#### 10.8.1 为什么重写：第一版实施与堆腐蚀事件

第一版将拟合状态（`fitCostModel` 相关成员）加为 `dsmcReplicatedMesh`
的类成员（改 .H），13M-cell 规模 bjm8 作业（4750806-4751264）在
MPI_Init 内部 malloc 处报 `corrupted double-linked list` 挂死。
排查结论：

- 腐蚀**检测点**在 MPI_Init，但写入发生在更早的构造阶段，新代码
  彼时未执行；
- pre-cost-model 二进制同阶段干净 ⇒ 腐蚀是"潜在越界写"被对象布局
  偏移（.H 改动 → 成员偏移变化）暴露，而非新代码直接越界；
- 该潜在越界属于存量代码（P0/P1/P2 之前就存在），在原布局下落在
  填充区，未触发检测。

处置（用户指示）：不用备份、不改暂存区，重新思考实现。方案：
**零 .H 改动**的 file-static 状态，使对象布局与"200 步正常"的
staged 版本逐字节一致，存量潜在越界回到其原（良性）位置。
注：潜在越界本身仍待后续定位（独立任务），本次不扩大改动面。

#### 10.8.2 实现形态（全部在 dsmcReplicatedMesh.C 内，无 .H 改动）

1. **状态**：`namespace Foam` 内结构体 `dsmcReplicatedCostModelState`
   （enabled/fitted/7 个系数/4 个相位 wall 快照/move 残差）+
   全局唯一实例 `costModelState_`；拟合函数 `dsmcFitCostModel`
   为文件内 static。
2. **开关**：`initialize()` 读 `replicatedMeshDLBCostModel`
   （默认 false），开启时打印确认行。
3. **拟合触发**：`autoRebalance()` 两个 auto 分支
   （legacyWindow 与默认 sarEval 窗口）的 `maxT > SMALL` 内、
   `printPerRankStepWall()` 之后调用 `dsmcFitCostModel`。
   此刻相位 wall（move/coll/build/post）与 occupancy 均为最新。
4. **拟合内容**（按 §10.3 设计，公式计数器版，零热路径插桩）：
   - 每 rank 相位 wall 增量 ΔW（对上次快照，负值截 0）；
   - 每 rank 计数器和：ΣN、ΣN/ℓ（ℓ=cbrt(V) 跟踪代理）、
     Σcand（闭式 0.5·nC(nC−1)·nParticles·σTcR·dt/V，仅 nC>1）、
     Σactive；
   - `MPI_Allgather` 8 标量/rank，rank 内最小二乘（非负截断）：
     move = a·N + b·N/ℓ（2×2）、coll = c·cand（1×1）、
     build = d·N + e·act、post = f·N + g·act（2×2）；
   - **共线性降级**：det 判据失败 → 该相合并为 N-only 系数；
   - **拟合门**：coll 相 Σcand=0（无有效窗）→ 本窗放弃（不 fitted）；
   - 质量：move 相相对残差（最大相位）记录并打印。
5. **权重接线**（`reassignByParMetisAdaptiveRepart`）：
   - `useCostModel = enabled && fitted` 时强制 `ncon=1`；
   - 每个活跃 owned cell 权重
     `w = (moveN+buildN+postN)·N + moveIt·N/ℓ + collCand·cand
     + (buildAct+postAct)·[N>0]`，
     N 取 occupancy（0 时回退 globalCellParticles）；
   - rawW 为秒量级 ⇒ 按 **全局最大值归一**（Allreduce MAX，
     目标 2^24，正权重截断照旧），避免全 1 退化；
   - 双约束分支整体退化为 else 分支（useDualConstraint 在
     cost model 下自动失效，ubvec[1] 块由 ncon==2 守卫不受影响）；
   - Info 打印权重模式（costModel weights / N^alpha weights）。
6. **可观测性**：`report()` 在 DLB 块内输出 fitted/残差/7 系数。

#### 10.8.3 与设计的偏差

| §10 设计 | 实施版 | 原因 |
|---|---|---|
| NNLS 逐相拟合 | 非负截断最小二乘（clamp 后重解） | 避免 NNLS 依赖，2×2/1×1 闭式解 |
| moveItersPerCell 计数器 | ΣN/ℓ 跟踪代理（cbrt(V)） | 零插桩；ℓ 即"每 cell 有效线性尺寸"，N/ℓ ∝ tracking 工作量 |
| cand 计数器 | 闭式公式（同 NTC :527-533） | 零插桩，公式与实际选择完全一致 |
| R² 逐窗门 | move 相相对残差 | 简化；coll 窗空门保留 |

#### 10.8.4 验证状态（后三项结果见 10.8.5 / 10.8.6 / 10.9.1）

- [x] 编译通过（linux64IccDPInt32Opt）
- [x] 本地回归：OFF 门通过（10.8.5）
- [x] bjm8 500 步：堆检测干净 + 拟合触发 + 权重生效（10.8.5 前段、
      10.8.6 前段），但同段 wall 回归 +51%→+275%，用户取消
- [x] A/B：max-anchored（4768139）与 mean-anchored+P99 cap
      （4770718）均回归，归一化假设证伪 → 方向关闭（10.9.1）；
      ≥2000 步长窗 A/B 不再需要

#### 10.8.5 本地回归结果（2026-09-21，zb 300 步 mpirun -np 8，独立目录 zb-dlbcost-off/on）

| 项 | OFF（开关缺省） | ON（cost model） |
|---|---|---|
| 300 步正常结束/无 fatal/stuck | ✓ | ✓ |
| 开关确认行 | 无（0 次） | "per-cell cost model DLB weights enabled" |
| 拟合触发 | — | 7 次成功；move it=0 共线性合并按设计降级 |
| ParMETIS 权重模式 | N^alpha | costModel weights（4 次重分区） |
| 全局粒子 step270/280 | 1901473 / 1920575 | 1901877 / 1920755（0.02% 内） |
| Total energy（末窗） | 0.0019619 | 0.0019623 |
| acceptance rate | 0.5251-0.5256 | 0.5209-0.521（物理等价） |
| rank0 local particles | 279123 | 182127 |
| rank0 local Collisions | 71430 | 3463 |
| ExecTime step270 | 101.9 s | 96.9 s |

注：Collisions/local particles 为 rank0 局部诊断（global reduction
skipped）。两 run 差异 = cost model 权重重分区后 rank0 负载下降，
即权重生效的直接证据；全局量（粒子/能量/acceptance）一致 ⇒ 物理
等价。本地小 case 上 ON 略快 5%（单机噪声量级，不作结论）。

结论：file-static 实现的功能门全部通过；MPI_Init 堆腐蚀回归需
bjm8 13M 规模 500 步确认（本地 case 无法复现该规模）。

#### 10.8.6 bjm8 500 步验证：回归发现，用户取消（作业 4768139，2026-09-21）

file-static 版通过堆腐蚀门 + 拟合门（残差 3.9%）+ DLB 触发
（step50/100，imbalance 4.21→2.72），但**同段 wall 对比显示显著
回归且逐步恶化**（m12o48，对照 staged 基线 4764853）：

| step | cost model (4768139) | 基线 (4764853) | 慢 |
|---|---|---|---|
| 50 | 298 s | 197 s | +51% |
| 100 | 617 s | 240 s | +157% |
| 150 | 872 s | 287 s | +204% |
| 200 | 1267 s | 338 s | +275% |

基线 step250/400 的 imbalance 已在 1.50-1.52 附近，而 cost model
版 step100 后步时持续变差 ⇒ 重分区#1/#2 产出的 cellOwner 分区
劣于静态 Scotch 初分区（或逐步碎片化）。用户指示停止，作业于
27:13（step ~240）取消。

**初步归因（首要嫌疑）**：权重归一化方案。本实施用"全局最大值
→ 2^24"，若个别 cell 的 cand 项极端（σTcR·N²/V 热点），全局 max
被其主导，典型 cell 权重被压到 1-几的量级 → ParMETIS 眼中典型
cell 近乎等权 → 分区质量劣化 + 界面/迁移负担上升。次嫌疑：
cand 项把碰撞热点权重放得过大，导致热点 cell 集中而拉长尾部。

**教训**：拟合优度门（残差 3.9%）只保证"模型能解释 wall 分解"，
不保证"权重分布适合 ParMETIS"——两者之间还隔着权重标度/分布
形态。后续若重启本方向，先把归一化改为 mean-anchored（典型 cell
→ 4096 级）+ cand 项截断/压缩（如 log1p 或 cap 分位数），单独
A/B 归一化，不与模型本身捆绑验证。

### 10.9 归一化修复：mean-anchored + cand 全局 P99 cap（2026-09-21）

按"改动收敛到一处，单独 A/B 验证归一化"的原则重写
`reassignByParMetisAdaptiveRepart` 的 cost model 权重分支。
**模型公式与拟合系数完全不动**，仅替换权重后处理：

1. **cand 全局 P99 cap**：cand 原始值先存 `candW`，128 个 log2 桶
   （偏移 +64，覆盖 2^-64..2^63）直方图 → `MPI_Allreduce` 全局直方
   图 → 从高桶向低扫描累计，超过总数 1% 的首个桶的下界即 cap
   （最多 ~1% cell 被钉住，且钉到同一值，保留"热点/典型"比值）。
2. **mean-anchored 归一化**：`rawW` 用 capped cand 计算，
   Allreduce(SUM) 得全局 ΣrawW 与活跃 cell 数 → 全局均值；
   `wScale = 4096/mean` —— 典型活跃 cell 权重 ~4096，热点 cell
   ≤ ~4e6（仍远离 ParMETIS 上限 1.3e8），空 cell 权重 1。
3. rank0 打印 scale/candCap/全局均值/活跃数/cand 数供诊断。
4. 删除 max→2^24 路径（已证伪，不留开关）。

A/B 语义（三组同段 wall，50/100/150/200 步）：
- 基线 N^alpha：4764853（已有，500 步完赛）
- max-anchored：4768139（已有，+51%..+275% 回归，step240 取消）
- mean-anchored + P99 cap：本次作业

本地回归（zb-dlbcost-mean，300 步）结果待补。

#### 10.9.1 A/B 结果：归一化假设证伪，真正根因 = act 固定项主导权重

作业 4770718（mean-anchored + P99 cap，500 步计划）step150 前的
同段 wall（ClockTime）：

| step | mean-anchored | max-anchored(4768139) | 基线(4764853) |
|---|---|---|---|
| 50 | 262 s | 298 s | 197 s |
| 100 | 537 s | 617 s | 240 s |
| 150 | 840 s | 872 s | 287 s |

mean-anchored 只带来 ~10% 改善 ⇒ **归一化不是主因**。真正的
失败信号是 imbalance 历史：cost model 版 4.08(step50) → 3.27
(step100) → **8.62(step150)**，基线 step250 即收敛 1.51 ——
重分区在把分区越调越差。作业已于 19:51 取消。

**根因（数值实证，step50 窗拟合系数）**：权重四项中
- act 固定项 = buildAct+postAct = 5.4e-5/活跃 cell —— 占全局
  mean rawW（7.27e-5）的 **74%**；
- N 项 = moveN·N ≈ 8.4e-7×17.5 ≈ 1.5e-5，仅 ~20%；
- cand 项被 P99 cap（cap=1-2，cand 分布极度长尾，99% cell <1）
  压至 ~5% 以下。

⇒ 权重语义退化为"均衡活跃 cell 数"，对粒子数差异（N 项）几乎
不敏感（20%）。粒子分布高度不均（起步 imbalance 4.08）时，
ParMETIS 按 cell 数均分 → 密集区高 N cell 无法被减载 → 失衡
持续并漂移恶化，且每窗都触发无效重分区（ParMETIS+迁移白付）。
基线 N^1 权重直接均衡粒子数 ⇒ 收敛 1.5。

act 项来自拟合：build/post 的 wall 大头确实是"每活跃 cell 固定
遍历成本"，**模型能解释 wall 分解 ≠ 固定项该进分区权重**——
分区权重的语义应当只保留随负载可迁移的部分（per-particle /
per-candidate 增量），固定项对所有 rank 是公共底数，进权重只会
稀释信号。这是继 §10.8.6 之后的第二层教训。

**方向处置**：
- cost model 保持默认 off（开关控制），生产路径（N^alpha）不受
  影响；两份 bjm8 日志（4768139/4770718）与三组同段数据已留档。
- 若未来重启：权重只保留 a·N + b·N/ℓ + c·cand（剔除 act 项），
  等效于"N 增量 + 碰撞增量"的纯可迁移成本；但鉴于基线已收敛
  1.5、收益上限（move 相失衡贡献）有限，本方向标记为**暂停**，
  优先级让位给其他工作。

## 11. move 相位 per-parcel 成本空间差异研究（2026-09-21）

### 11.1 问题定义（承接 §9.3 结论 3）

§9.3 实测 move 相位 rank 级 2× 分散而粒子数仅差 1.34×，判明
per-parcel tracking 成本随区域变化 ~1.5×，且 step-wall 最大的 rank
是 move-轻 rank（等迁移同步）。本节目标：定位该成本差异的来源。

### 11.2 诊断工具：per-rank move 子相位表

发现既有 `moveDetailProfile` 设施（trackingData 的 moveTrackCalls/
moveTrackWallTime/faceHits 分类/boundaryWallTime，chrono 插桩），
但只输出本 rank 累计、无 per-rank 汇总。本次在
`printProfileSummary()` 增加 per-rank Gather 表（10 值/rank：
parcels、trackCalls、faceHits、cyclicHits、patchHits、procHits、
trackWall、trackerWall、bndWall、moveWall），rank = replicated mesh
下的空间区域。开关：`profileDetail true` + `moveDetailProfile true`
（仅诊断 run 使用，符合 AGENTS §3）。

实施中踩了三个工程坑（均已修复并记录）：
1. **Gather 放进 rank0 门控块** → 其余 rank 不可达，集体通信挂死
   （与 §5.6 report() Allreduce 同款错误，第 2 次发生——规则升级：
   在该函数内新增任何 MPI 集合通信，必须先核对函数的 rank 门控结构）；
2. **label 累计 int32 溢出**：200 步 × 1.5e8 粒子，累计计数远超
   2^31。改 file-static double 累计数组 `dsmcMoveDetailCum[12]`
   （零 .H 改动，精确到 2^53；单步内 trackingData label 计数器
   不溢出，逐步步加到 double 安全）；
3. **`Info<< long(x)` 静默降级**：OF-1706 Ostream 无 long 重载，
   long 隐式转 label(int32) 打印 → >2^31 的计数回绕显示为负
   （与此前 unsigned long 无重载同族）。计数列全部改 scalar 打印。
   （注：tc/parcel 等 scalar 列始终正确，第 2 次诊断 run 的
   max/min=1.053 结论已成立，第 3 次 run 仅修正显示。）

### 11.3 数据（作业 4775054，1200w-rcj-dlbopt，200 步，m12o48，cost model off）

关键 per-rank 量（完整 12 行见 4775054.out）：

| 量 | max/min | 备注 |
|---|---|---|
| parcels（累计 move 调用） | **2.29×** | r6 最少 1.34e9，r10 最多 3.06e9 |
| trackCalls/parcel（跨 cell 率） | **1.045** | 1.6495-1.7242，基本均匀 |
| faceHits/trackCall | 1.070 | 命中模式均匀 |
| trackWall/parcel | **1.63×** | 6.40e-7（r8）→ 1.04e-6（r6）s |
| trackWall/trackCall | **1.56×** | 3.86e-7（r8）→ 6.04e-7（r6）s |
| moveWall | **1.95×** | 35.2（r8）→ 68.7（r9）s |
| bndWall | r0/r6/r7/r9 = 59-74 s；内部 r4/r5/r8/r10/r11 ≈ 0.1-0.2 s | 强空间结构 |
| patchHits/parcel | 125× | 近壁 rank 7e-3 vs 内部 1e-4 |

相关性：corr(moveWall, trackWall)=0.97、corr(moveWall,
parcels)=0.86、corr(per-parcel 成本, bndWall/parcel)=0.87、
corr(per-parcel 成本, tc/parcel)=0.96（r0/r6 双高驱动）、
corr(per-parcel 成本, patchHits/parcel)=0.68。

### 11.4 归因

1. **跨 cell 率不是来源**：tc/parcel 均匀（1.045）。ocm-moss3d
   网格均匀，§10 cost model 里的 N/ℓ 项在该 case 无空间信号
   （与其分区失败互相印证——该模型没有可用的空间成本信号）。
2. **per-parcel 成本差异 1.63× 的主体 = 每次 trackToFace 调用
   成本差 1.56×**（不是调用次数差），集中在近壁/对称面 rank
   （r0/r6/r7/r9：patchHits/parcel 高 7e-3-1.1e-2，内部 rank
   ~1e-4）。候选机理：近壁网格 tet 形状/尺寸（拉伸层）加深
   tet 遍历、boundary 面 test、RWF。
3. **边界处理 critical 串行化是独立确认的热点**：
   `#pragma omp critical(dsmcMoveBoundary)` 包住 patch/cyclic
   控制（含 cloud RNG 与 wall measurements 共享写）。近壁 rank
   bndWall（thread-sum，含等锁）59-74 s/200 步，内部 rank≈0；
   bndWall/parcel 差 1118×、单次边界处理成本差 13×（r0 vs r5）
   → 锁争用放大。这是明确的、可削减的工程热点。
4. **moveWall 分散 = 粒子数失衡（2.29×，DLB 可修）× per-parcel
   成本空间差（1.63×，DLB 不可修）**。per-parcel 成本是 rank 的
   空间属性（近壁就是贵），cell 所有权再怎么调，近壁 cell 的
   单位成本不降——**这正是 §9/§10 所有 DLB 权重方案失败的根本
   原因**：权重只能均衡"N×f"里的 N，f 是位置属性。

### 11.5 下一步（按收益排序）——后已全部处理：方向1/2 结论见 §11.6/§11.7

1. **削减边界 critical 串行化**（近壁 rank 直接受益）：per-thread
   RNG（fastRng 路径已部分存在）、wall measurements 分片累加 +
   串行合并、cyclic 控制无锁化。预期收益：r0/r6/r7/r9 的
   bndWall/48 ≈ 1.2-1.5 s/thread/200 步的串行段削减。
2. **近壁 per-call 成本机理**：离线从 constant/polyMesh 计算
   per-rank 平均 tet 体积/面数/拉伸比，验证"近壁 tet 遍历更深"
   假设；若成立，评估 trackToFace 的近壁路径优化。
3. DLB 权重方向维持关闭（§9/§10/§11.4-4 的完整证据链）。

### 11.6 nearWall/interior 分账实验 + 方向 1/2 结论（2026-09-21）

#### 11.6.1 实验：按 cell 贴壁性分账 trackToFace

新增 per-track 分账：parcel 所在 cell 贴任意 boundary 面（nearWall）
与否（interior），分别累计 trackCalls/chrono。实现：
trackingData 加 4 计数器 + nearWall boolList 指针（dsmcCloud.C
file-static 惰性构建，boundary 面 owner 集合，零 .H 数据成员）；
Cloud.C 的 td 合并适配器同步扩展。修复时序 bug：类别必须在
trackToFace **前**取（调用后 cell() 已变化）。

#### 11.6.2 数据（作业 4775404，1200w-rcj-dlbopt，200 步，m12o48）

| rank | wall per-call μs | interior per-call μs | 比值 | wall 占 trackWall |
|---|---|---|---|---|
| r0 | 1.89 | 0.509 | 3.7× | 6.1% |
| r2 | 1.93 | 0.585 | 3.3× | 7.4% |
| r7 | 1.22 | 0.401 | 3.0× | 7.7% |
| r9 | 1.03 | 0.420 | 2.5× | 6.6% |
| r4/r5/r8/r10/r11 | 0.38-0.49 | 0.37-0.41 | 1.0-1.2× | 2-3% |
| r1/r3/r6 | 0.57-0.74 | 0.40-0.41 | 1.4-1.8× | 3-5% |

#### 11.6.3 发现

1. **per-parcel 成本差异的主因是 interior per-call 的 rank 差
   1.57×**（0.373-0.585 μs）：interior 占 trackWall 的 92-98%，
   corr(trackWall/parcel, interior per-call) = **0.997**。与跨面率
   （均匀 1.07）、tc/parcel（均匀 1.05）均无关——是 trackToFace
   调用**内部**工作（tet 遍历深度/访存）的区域差异。
2. **corr(per-parcel 成本, parcels) = -0.66（负）**：粒子稀疏的
   rank（r0/r2/r6）per-call 反而贵——支持 **cache 效应假设**：
   粒子密 → cell 邻域数据热 → track 快；稀疏区冷数据访存多。
3. **nearWall cell 在复杂几何 rank（持 Part58 切割层的 r0/r2/r7/r9）
   per-call 贵 3-4×**（snappyHexMesh 贴体切割 cell 的 tet 更碎），
   在只有平面端面的 rank 与 interior 相同——网格属性；但 wall
   cell track 只占 2-8%，非主因。
4. 跨面率/调用次数均匀 ⇒ 跨 cell 机制整体排除。

#### 11.6.4 方向 1 初步决策：评估后不做（后按用户指示实施验证 → §11.7 证伪关闭）

- ratio 证据（moveWall/(trackWall/48)）：近壁 rank 1.30-1.32 vs
  内部 1.30-1.41 —— 边界 critical 串行化的墙钟代价不显著；
- FreeStreamInflow 的工作在串行 controlParcelsBeforeMove（不在锁内）；
  Deletion 天然线程安全（porousMeas no-op + keepParticle=false）；
- interior 占 92-98% ⇒ 即使 bndWall 全免，可动份额 <8%。
  thread-RNG/分级免锁方案保留在记录中，不实施。

#### 11.6.5 总结论

move 相位 per-parcel 成本的空间差异（§9.3 发现的 ~1.5×）分解为：
interior trackToFace 内部成本差 1.57×（主因，与粒子密度负相关，
cache/访存属性 + 局部位移场差异）+ 复杂几何切割 cell 放大（局部，
占比小）。两者都是**物理/网格/访存属性，非 cell 归属属性**——
DLB 权重原理上无法修复，§9/§10/§11 的失败证据链闭环。

可行动方向的残余空间（均低优先）：
- 粒子排序/布局优化以改善稀疏区 cache 命中（moveOrderedParcels
  基础已存在，需验证稀疏-密集 rank 的收益差）；
- trackToFace 热路径的访存审计（cell 定位查表结构）。
当前生产配置（staged 基线 + 已落地的 P0-P5）维持不变。

### 11.7 方向 1 实验：per-segment cell 排序——证伪关闭（2026-09-21）

#### 11.7.1 实施

- `openmpMoveSortByCell`（controlDict，默认 false）：move 遍历数组
  在 Cloud.C 汇合点（threadOffsets 定型后）按线程段 std::sort by
  cell()。fallback（链表序）与 reuse（持久数组）两路径都覆盖；
  reuse 路径的排序序持久化。dsmcCloud 开关用 file-static 缓存
  （零数据成员，布局不变）。
- 正确性：统计等价（遍历/RNG 消耗序变化，物理不变）。

#### 11.7.2 本地 A/B（zb omp8，300 步，OMP 8）

| | OFF | ON |
|---|---|---|
| move max [s] | 40.70 | **55.78（+37%）** |
| Total energy | 0.0019628 | 0.0019630（统计等价 ✓） |

#### 11.7.3 判定：证伪，方向关闭

排序成本（指针比较排序本身 + 每 cell 拓扑数据仍须逐个读取）远超
局部性收益。与 §11.6 的 cache 假设对照：假设方向正确（稀疏 rank
per-call 贵），但**排序不是经济的干预**——冷数据的量不变，改变的
只是访问序；而每步排序的成本是确定性的。周期化（每 K 步排）会
同时衰减收益（粒子每步跨 ~0.65 cell，排序序 1-2 步内半失效），
净期望仍为负。1200w 规模的排序成本更高（1.25e7 指针/rank），
不投入 bjm8 机时。

#### 11.7.4 优化阶段收口

至此 §11 盘点的全部可行动方向均已处理：
- 已落地：back-off/wall 触发、P0-P5（build/post OMP 化）、
  fastRng、延迟接收等（生产基线 40000 步 19h12m 正常完赛）；
- 已证伪关闭：DLB 权重×3、collision-aware、NoAlltoall、时序重叠、
  边界 critical 免锁、wall-cell track、**cell 排序（本轮）**；
- 被约束封死：tracking 拓扑数据结构优化（OF 核心类，不改源码）。

结论：当前生产配置（staged 基线 + 已落地优化）接近该代码架构在
本硬件上的实际极限。move 相位（wall 70%+）由 OF trackToFace 拓扑
结构与密度分布的物理属性锁定；post/build/通信已完成可用优化。
后续工作建议转向生产任务推进。

## 12. 开关清理（2026-09-21，证伪开关代码级删除）

### 12.1 分析原则

- 证伪 + **零对象布局风险**（file-static / 函数局部变量）→ 直接删除；
- 证伪但删除需改 **.H 数据成员** → 保留（默认 false），原因：
  §10.8.1 定位的存量潜在越界写尚未修复，dsmcReplicatedMesh.H /
  dsmcCloud.H 的任何数据成员增删都会再次偏移对象布局，可能重新
  触发 13M 规模堆腐蚀（4750806 事件）；布局越界定位并修复前，
  .H 数据成员只减不加不删。

### 12.2 处置清单

| 开关 | 证伪依据 | 形态 | 处置 |
|---|---|---|---|
| `replicatedMeshDLBCostModel` | §10.9.1 | file-static（.C） | **已删除**（state/fit/vwgt 分支/report，~300 行） |
| `openmpMoveSortByCell` | §11.7 +37% | file-static（.C） | **已删除**（访问器/排序块/开关，3 文件） |
| `replicatedMeshDLBDualConstraint`（+Ubvec1） | §9.3 更差 | 函数局部变量 | **已删除**（vwgt 双分支并入单 N^alpha 分支） |
| `replicatedMeshDLBAdaptiveAlpha`（+5 子参数） | §9.3 中性 | .H 数据成员 ×8 | **保留**（布局风险） |
| `replicatedMeshNoAlltoall` | +45% | .H 成员+迁移路径分支 | **保留**（布局风险+大范围 diff 风险） |

保留理由补充：adaptive alpha 证伪结论为"中性"（无害），NoAlltoall
虽有明确负收益但删除需贯穿迁移路径大改；两者默认 false 零运行成本，
controlDict 不写即不激活。若未来完成存量越界定位修复，可一并清理。

### 12.3 验证

- 编译通过（BUILD=0，dsmcReplicatedMesh.C / dsmcCloud.C / Cloud.C）；
- 本地回归（zb mpi8，300 步，独立目录 zb-clean-check）：正常
  End、已删开关输出 0 命中、Total energy = 0.0019617（基线
  0.0019619 波动内）、全局粒子 1958111 一致。
- 生产路径（N^alpha 权重 + 全部保留开关）行为不变。

## 13. 500wcell-1bparticle DLB 卡死诊断 + back-off 激活缺陷修复（2026-09-22）

### 13.1 事件

作业 4776404（500wcell-1bparticle：5,822,470 cells、1.98 亿粒子、
36 ranks m4cl0505-0507）报"卡在 DLB"。实际：step 16100 的
rebalance #160 进入 ParMETIS AdaptiveRepart 后 12 小时无输出
（.out 最后写入 02:42:15，用户 14:4x 取消）。

取证：
- gdb 抓栈（0505 ×2 rank、0506 ×1 rank）：全部停在
  `ParMETIS_V3_AdaptiveRepart → CommSetup → MPI_Alltoall(count=1,
  hcoll ucx p2p) → uct progress 轮询`；
- 进程 RSS 43GB/rank（ps），作业步聚合 MaxRSS **521 GB**；
- 0507 未抓到栈（用户随后取消，pam_slurm_adopt 拒绝再连）。

### 13.2 代码级根因：back-off 激活缺陷（已修复）

该 case 的 DLB 稳态：每 100 步触发、ParMETIS 输出与现状 100% 不同
（changedBefore=5822470）、greedyOverlap remap 全盘否决（net
changedAfter=50，比率 8.6e-6）、失衡 1.66 物理性存在 → `futile`
判定每次都成立——但 **159 次重分区全部执行、back-off 0 次跳过**。

缺陷（`dsmcReplicatedMesh.C` autoRebalance 执行段）：`futileCount_`
在每次执行后被**无条件清零**，而触发检查时的累积在下次执行前必然
又被清零——`futileCount_` 恒为 1 < `futileRuns_=2`，back-off 在
"remap 持续否决 ParMETIS"的稳态下**数学上不可能激活**。

修复：清零改为条件式——仅当本次执行有效
（`lastExecChangedRatio_ >= backOffChangedRatio_`）时重置；
`backOffRemaining_ = 0` 保留（执行后重置跳过状态）。
修复后稳态行为：第 2 次连续无效 → 启动 6 窗跳过 → 每逢 remaining
耗尽重置 6 窗 → ParMETIS 调用趋零；imbalance 突涨 ≥
threshold×1.5=2.25 时 escalate 安全阀仍强制执行。

**教训**：back-off 当初的 A/B（§7.5）在 1200w 上做——该 case
remap 有效、futile 永不成立，激活路径从未被测过。"设计已实现"
≠"实现可达设计语义"——对条件触发的防御机制，必须构造触发条件
本身做一次验证。

### 13.3 相关联审计（同日）

- `nChanged > 0` 守卫（rebuild/localMesh :1947、迁移 :2508、
  后处理 :2544/:2566）：在位且正确；本 case 每次 nChanged=50>0，
  全量链照付——但 rebuild/migrate 是正确性必需（owner 变了就必须
  跟），**不可加比例阈值**；省 ParMETIS 调用只能靠 back-off。
- `printPerRankStepWall`：单一口径实现，§5.6 的"备选口径"无双实现
  残留，无需清理。
- 层级关系：back-off（外层，修好）→ nChanged 守卫（内层，正确）。

### 13.4 重启配置建议（500wcell-1bparticle）

1. `export I_MPI_COLL_EXTERNAL=0`（run-hpc.slurm）：本次死锁的直接
   机制在 hcoll/ucx 集合点（36 rank 全员等 Alltoall 12h），禁用
   hcoll 回退 Intel MPI 原生集合通信；
2. back-off 修复版二进制（已同步编译）；
3. 可选：`replicatedMeshDLBMinGapSteps 100→300`，进一步降低
   检查频率（该 case 失衡 1.66 物理性存在，检查本身便宜，非必需）。

### 13.5 back-off 修复的生产验证（作业 4783381，2026-09-22）

用户在 500wcell-1bparticle-syncproblem（粒子数缩小至 ~1000 万级）
以修复版二进制 + `I_MPI_COLL_EXTERNAL=0` +
`replicatedMeshDLBImbalanceThreshold 1.8` 完整跑完 40000 步。

**正确性门（AGENTS §6 全过）**：COMPLETED、Total Iterations=40000 /
End main、无 fatal/NaN/abort/corrupt、stuck=0、粒子 972 万 → 1326 万
（物理发展正常）、Total energy 8233.3。

**back-off 修复的直接验证（首次真实激活）**：
- `skipped by back-off` ×6，remaining 序列 6→5→4… 严格符合设计；
- lastChangedRatio=2.97e-4 < 1e-3，futile 判定正确成立；
- **9 次 rebalance / 40000 步**（每 ~4400 步一次）vs 缺陷版
  500wcell 的每 100 步一次——重分区频率降 ~44 倍；
- DLB wall max = 19.9 s（全程累计）vs 缺陷版单次卡死 12 h；
- report 尾行 futile count=0 / lastChangedRatio=1.61e-4（探测窗
  行为正常：每 6 窗跳过后 1 次验证执行）。

**性能**：ClockTime 7026 s（2h01m）≈ 0.176 s/步（13M 粒子 36
ranks）；40000 步累计 max：post 2996 s（30%，最大项）、comm 1889 s、
move 1152 s、build 998 s、coll 232 s。I_MPI_COLL_EXTERNAL=0 下
40000 步集合通信零故障——§13.1 的 hcoll 卡死机制未再现。

**残余观察**：粒子数在 +36% 增长中（注入>删除）；若目标形态是
10 亿粒子，内存/DLB 压力将回升至 4776404 的量级——缩小粒子数的
本次验证覆盖了 DLB/通信路径，未覆盖极限内存形态。
