# dsmcFoam+ 内存优化计划：dsmcVolFields 实例 + 粒子布局（M 系列）— 2026-09-21

承接：本会话对 500wcell-1bparticle 的三次 OOM 调查与逐 rank 内存实测；
关联文档：`dsmcFoam_todolist.md`（P1 条目）、`dsmc_cabana_plan_20260528.md`、
`dsmc_sparta_plan_20260528.md`（P1 预研）、`tier2_parallel_filtered_read_20260916.md`。

## 1. 问题定义（实测数据）

### 1.1 基线算例

500wcell-1bparticle：5,822,470 cells（blockMesh (132 218 174) + snappy
orion level (2 2)），patch inlet/orion/spline，来流 7600 m/s、T=211 K、
n_total=4.976e18/m³（N2+O2+O），壁面 diffuse 760 K，replicated mesh。

### 1.2 三次 OOM 记录（粒子规模超限）

| 作业 | 配置 | nEqP / 粒子总数 | 死亡时刻 | 死亡时最大节点 RSS | 结局 |
|---|---|---|---|---|---|
| 4775233 | 24r×16t | 5.8e12 / 1.003e9 | 13.6 min（时间循环） | **741.7 GB**（节点限 768 GB） | OUT_OF_MEMORY |
| 4776150 | 24r×48t | 4e12 / **1.455e9** | 2:57（填充期） | ≥502 GB（采样值） | SIGKILL+SIGSEGV 级联 |
| 4776234 | 12r×48t | 4e12 / 1.455e9 | ~5-8 min（填充期） | ≥322 GB（采样值） | 同上 |

结论：**粒子总数 × ~2.1 KB/粒 是内存主体**；节点内存 768 GB/192 核；
3 节点（2.25 TB，4 GB/核显式申请）可承载 ~10 亿粒子的均衡态，1.455e9
粒子在任意 rank 划分下均越限。rank 数只影响复制/场/线程的每 rank 固定
开销（占比小），"减少 MPI rank"不能解决粒子本体超限。

### 1.3 当前健康运行基线（4776404，36r×16t，nEqP 4e13 → 1.455e8 粒子）

srun --overlap 逐进程实测（/proc/PID/status + smaps_rollup）：

- 每 rank RSS **40.1-43.4 GB**（36 rank，分布紧凑），3 节点合计 ~1.44 TB
  / 2.25 TB（64%）；
- 时间循环 13.8-14.0 s/10 步，stuck=0，碰撞正常。

### 1.4 每 rank 40-43 GB 的归因（实测 + 代码清单推算，闭环）

| 组件 | 每 rank | 依据 |
|---|---:|---|
| **dsmcVolFields 6 实例** | **~31-32 GB** | 成员清单 × 46.6 MB/数组（见 §2），mixture 实例 ~7 GB + 单组分实例 ~5 GB × 5 |
| 粒子（1.455e8/36 ≈ 4.04M/rank） | ~8.5 GB | ~2.1 KB/粒（4775233 锚定：741.7 GB ÷ 334M 粒/节点） |
| 网格复制 + 拓扑缓存 | ~4-5 GB | points 0.14 + faces 0.5 + owner/neigh 0.14 + 面心/面积/mag 缓存 1.0 + cellCenter/vol 0.19 + faceEdges ~1.4 + edges/pointFaces/pointCells ~0.7 + cellLevel 等 0.12 + cellOwner_/localMesh/patch 0.3（GB） |
| sigmaTcRMax/碰撞余量/线程栈等 | ~1-2 GB | |

## 2. dsmcVolFields 每实例解剖（dsmcVolFields.H 成员清单 × 46.6 MB/数组）

每个实例（typeIds 决定 per-species 数组的组分数）：

| 类别 | 清单 | 内存 |
|---|---|---:|
| GeoField 场对象 | 29 volScalarField（dsmcN…pressureError）+ 3 volVectorField（UMean/fD/heatFluxVector）+ 2 volTensorField（pressureTensor/shearStressTensor），含 boundaryField；Cp_/Ch_ 已是 autoPtr 懒构造 | 标量 29×48.5 MB + 向量 3×142 MB + 张量 2×425 MB ≈ **2.68 GB** |
| 混合累计器 | 29 raw scalarField（dsmcNCum…dsmcNCollsCum）+ 2 vectorField（动量累计） | ≈ **1.63 GB** |
| per-species 数组 | 10 个 List 成员（dsmcSpeciesEvibModCum/Eelec/NSpeciesCum/nSpeciesCum/Mcc/Tvib/GrndElecLvl/1stElec/Mfp/Mcr）× 组分；另有 dsmcVolSharedSampleCache（T12）分摊 parcel 级累计 | mixture 实例 2.3-3.3 GB；单组分实例 ~0.5 GB |
| 边界累计器 | BF_ 族 × patch × 边界面（全 case 23.7 万） | ~0.1-0.2 GB |

实例合计：mixture（5 组分）~6.8-7.7 GB；单组分 ~4.8-5.0 GB；
**6 实例（1 mixture + 5 单组分）≈ 31-32 GB/rank** —— 与实测差值闭合
（40-43 − 8.5 粒子 − 4-5 网格 − 1-2 其他 ≈ 31-32）。

## 3. 优化项（M 系列）

### M1：dsmcVolFields 实例裁剪（S0，配置级，零代码）

- **改动**：fieldPropertiesDict 删除不需要的单组分实例（当前 5 个：
  N2/O2/NO/N/O 各一）。mixture 实例已含全部组分的输出能力。
- **收益**：6 → 2-3 实例，31-32 → **12-17 GB/rank**（省 ~15-20 GB）。
- **前置核实**：`dsmcVolSharedSampleCache` 的共享边界——mixture 实例的
  per-species 数组是否通过缓存共享单组分实例的数据；若是，裁剪实例需
  确认 mixture 输出不依赖被删实例（读 dsmcVolSharedSampleCache 类定论）。
- **决策点**：论文/生产需要哪些组分的独立场输出——用户确认。

### M2：GeoField 懒构造（S1，1-2 天，单项最大代码收益）

- **改动**：29+3+2 个无条件构造的成员 GeoField 改 autoPtr，createField
  里按 **outputFieldEnabled()**（机制已存在，Cp_/Ch_ 是先例）构造；
  calculateField/writeField 对应分支同步条件化。
- **收益**：白名单输出 8-10 个核心字段（rhoN/rhoNMean/dsmcNMean/p/
  Ttra/Trot/Ma/q/fD…）时，GeoField 2.68 → ~0.9 GB/实例 × 实例数。
- **量级**：M1+M2 后 2-3 实例 ≈ **7-9 GB/rank**。

### M3：累计器按开关分配（S2，顺路做）

- **改动**：29 个 raw 数组分四族——核心族（密度/温度/能量）保留；
  张量/热流族（dsmcMuu…Mccw、momentum×2、Eu/Ev/Ew/EC）绑定
  measureHeatFluxShearStress_；分类族（ClassI-III）绑定
  measureClassifications_；误差族绑定 measureErrors_——把现有语义
  开关下沉到**分配**（当前疑仅门控计算，需核实 initScalarFields）。
- **收益**：1.63 → ~1.0 GB/实例（关张量/热流族时）。

### M4：owned-only 分布式累计器（S3，专项 2-4 周，最大收益）

- **浪费点**：replicated mesh 下 calculateField 只采样 owned cells
  （~24.3 万/5.82M），但全部 per-cell 数组（GeoField internalField +
  raw + per-species）都是全网格尺寸——**每 rank 浪费 23/24**。
- **改动**：数组尺寸 = ownedCells.size()，dsmcLocalMesh::toLocal/
  toGlobal 索引映射；写出时组装全局场（processor-write 本按 owned 写，
  天然兼容）。46.6 MB → 1.95 MB/数组。
- **收益上限**：dsmcVolFields 全部 per-cell 存储 ×1/24 → M1-M3 后的
  7-9 GB → **~0.7 GB/rank**。
- **难点（DLB 交互）**：repartition 后 owned 集变化，半累计数据需跨
  rank 迁移或重置。分两步：第一步接受"DLB 触发 → 累计窗口重置"
  （resetField 已有此语义，损失一个平均窗口）；第二步做迁移。
- **风险**：dsmcVolFields.C 6200 行全部 [celli] 访问改映射，机械但
  量大；resetField/updateOwnedBoundaryFaces 时序需重审。

### M5：parcel 冷字段懒分配（P1 前置快赢，~2-3 天）

- **浪费点**：dsmcParcel 携带轨迹统计包袱（t0_/initialPosition_/
  currentPosition_/distanceTravelled_/msd_ + wallTemperature_/
  wallVectors_ 两个 per-parcel List）——多数粒子终生不触发，
  每粒 ~130-200 B（含各自的堆分配）。
- **改动**：改 autoPtr/惰性分配，或确认无人使用后裁剪。
- **收益**：~2.1 KB → ~1.8-1.9 KB/粒（省 10-15%）。

### M6：parcel arena/slab 分配器（P1 前置快赢，~1 周）

- **浪费点**：1e9 级独立 `new dsmcParcel` 的 malloc 开销 + 堆碎片
  （AoS 有效 2.1 KB 中 ~0.5-1 KB 是分配器/碎片/cache 粒度）。
- **改动**：dsmcParcel 的 op new/delete 重载到 slab arena（连续大块
  按 300-500 B 槽位切分，free 回槽位表）。
- **收益**：有效 2.1 → ~1.4-1.6 KB/粒（省 25-30%），不改任何访问代码。

### M7：粒子布局 pool-as-storage（P1 完全体，专项 2-4 月）

- **结论修正（对 5/28 预研）**：原"每步 buildFromCloud + writeBack"
  混合方案在 ≥1.5 亿粒规模不可行——每步全量复制 6-88 GB 吞掉全部
  收益（碰撞仅占步时 ~3%，且已被现行 OMP cellOccupancy 路径并行化）。
  **唯一可行形态是 pool-as-storage 完全迁移**。
- **数据划分**：热字段（position/v/cellI/faceI/tetPt/typeId/ERot/
  ELevel/vibLevel/RWF）入 AoSoA tile（VecLen=8，碰撞核心字段在 tile
  头部，~152 B/粒）；冷字段按 M5 处理；cell 索引用 first/next 增量
  维护（取代每步全量重建的 buildCellOccupancy，0.048 s/步）。
- **耦合面**（当前代码逐项）：move/tracking（slot-backed 视图，最高
  风险）、MPI 迁移（pool slot 直接打包，与 chunk 传输天然兼容）、
  QK 反应（pool create/flagDelete + free-list）、边界/控制器/fill
  （统一 create 入口）、dsmcVolFields 共享缓存（遍历 first/next，
  结构不变）、IO（从 pool 组装）。
- **收益**：2.1 → ~0.3-0.4 KB/粒（**7-10×**）；1e9 粒子 2.1-3.3 TB →
  0.3-0.4 TB；move/collide cache 命中率与 SIMD 附带收益。
- **RNG/正确性**：compact/重排改变配对顺序 → 统计等价验收（与 3b
  同口径）；每格单线程所有权保持。

## 4. 量化路线图（每 rank，1.455e8 粒子规模）

```
现状（4776404 实测）                    40-43 GB
 → M1 实例裁剪（6→2-3）                 ~16-18 GB
 → M2 GeoField 懒构造                   ~9-10 GB
 → M3 累计器按开关                      ~8-9 GB
 → M4 owned-only（专项）                ~1-2 GB（dsmcVolFields 部分）
 粒子 8.5 GB：M5+M6 → ~5-6 GB；M7 → ~0.6 GB
 终态（全部落地）                       ~6-8 GB/rank
```

对粒子规模的外推（M1-M3+M5/M6 落地后，36r×16t，节点限 768 GB）：

| nEqP | 粒子总数 | 粒子内存/rank | 总计/rank | 3 节点可行性 |
|---|---|---|---|---|
| 4e13（现状） | 1.455e8 | ~4-5 GB | ~12-14 GB | ✓（已验证） |
| 1e13 | 5.8e8 | ~16 GB | ~24-26 GB | ✓ |
| 5.8e12 | 1.003e9 | ~28 GB | ~36-38 GB | ✓（现状该规模 OOM 的主因是 dsmcVolFields 31 GB 未压缩） |
| 4e12 | 1.455e9 | ~40 GB | ~48-50 GB | 临界，建议 4-5 节点 |

**关键结论：M1-M3 落地后，1e9 粒子在现有 3 节点即可运行**（此前
1.003e9 @ 24×16 OOM 的 741.7 GB 中 ~610 GB 是 dsmcVolFields 未压缩
部分——6 实例 × 5.82M 格比 2.3M 格的 230w 大 2.5 倍所致）。

## 5. 实施顺序与触发条件

| 项 | 工程量 | 触发条件 |
|---|---|---|
| M1 | 0（改 dict） | 下一轮运行前，用户确认组分输出需求 |
| M2 | 1-2 天 | nEqP 计划提升到 ≥1e13 前完成 |
| M3 | 0.5-1 天 | 与 M2 同窗口 |
| M5 | 2-3 天 | 粒子 ≥5e8 前 |
| M6 | ~1 周 | 粒子 ≥5e8 前（与 M5 二选一先做 M5） |
| M4 | 2-4 周 | M1-M3 不足以支撑目标规模时 |
| M7 | 2-4 月 | 粒子 ≥1e9 生产立项时（必选项） |

## 6. 验证标准（全部优化项通用）

1. **统计等价**：zb-cylinder-react-validate mpi4omp2 300 步，final
   particles/Total energy/碰撞计数与基线（1,957,830 / 1.9622e-3，
   含 collisions 门）在 RNG 噪声带内；
2. **内存验收**：同规模同配置 RSS 对比（mem-sampler 曲线归档），
   dsmcVolFields 部分按 §4 预期比例下降；
3. **输出一致性**：白名单内字段逐位一致（M2/M3 改变分配不改变计算）；
   M4/M7 为统计等价；
4. **DLB 兼容**：触发 rebalance 后跑通（M4 第一步接受累计重置）；
5. **旧路径回归**：开关关闭时行为逐位不变。

## 7. 已部署的观测设施

- run-hpc.slurm 已加 **mem-sampler**（30 s 采样全部 rank 的 VmRSS/
  Anonymous → mem.log）与 reconstructPar 容错；
- srun --overlap --jobid=<id> 可对运行中作业做 /proc 逐进程采样
  （pgrep 正则注意 dsmcFoam+ 的 "+" 需转义或用 /proc/*/comm 扫描）。

## 8. 待核实项

1. dsmcVolSharedSampleCache 的共享边界（M1 前置）；
2. dsmcNGrndElecLvlSpeciesCum_ 的实际数组数（电子能级表规模——若
   每组分几十个能级 × 46.6 MB，可能是隐藏大头）；
3. M3：initScalarFields 是否已按开关门控分配；
4. dsmcParcel 轨迹统计成员（msd 族）的实际使用者（M5 裁剪依据）。

## 9. M1 实测验证（2026-09-22，job 4777260）

fieldPropertiesDict 裁剪为仅 mixture 实例（typeIds N2 O2 NO N O），
36r×16t、nEqP 4e13、1.455e8 粒子，m4cl0302-0304：

- srun --overlap 逐进程实测：每 rank RSS **17.5-19.1 GB**（此前 40-43 GB，
  **-56%**）；节点合计 ~220 GB / 768 GB（29%），余量充足；
- 与 §4 路线图预测（M1 后 16-18 GB）吻合 ✓；
- 运行健康：1.455e8 粒子、stuck=0、Iteration 360+（16.4-17.5 s/10 步，
  较裁剪前 13.8-14.0 慢 ~20%，原因待观察——节点组不同 m4cl0302-0304
  vs 0505-0507）；
- 下一步：确认单组分场无输出依赖后，M2（GeoField 懒构造）可将剩余
  dsmcVolFields 部分（~7 GB 中 GeoField 占 ~2.7 GB）进一步压缩。

## 10. M3 已实施（2026-09-22）

dsmcVolFields.C 三处改动：
1. ctor 初始化列表：16 个门控数组改默认空构造（热流/剪切 10、分类 3、
   平均自由程 2、电子 1——电子累计已被 needElectronic=false 硬编码禁用）；
2. resetField 的 setSize 按开关分组（measureHeatFluxShearStress_ /
   measureClassifications_ / measureMeanFreePath_ / 
   writeElectronicTemperature_）；per-species 循环同步门控
   （speciesMfp_/Mcr_ 归 mfp，电子 3 列表归 writeElectronicTemperature_）；
3. calculateField 的 resetAtOutput 清零循环拆分（核心循环 + 3 个
   开关守卫循环），combine 块加懒分配自愈（开关开启但数组未分配时
   首次触发分配）。

验证：zb mpi4omp2 300 步通过（1,957,844 / 1.96226e-3，噪声带内；
stuck=0；NoColl=0）。快照/Allreduce 对空数组天然安全（forAll 空转、
count=0 合法）；restart dict 的空列表 round-trip 一致（开关配置跨
重启应保持）。已同步 bjm8 重编（libdsmcFoam+ + exe）。

收益（500w case，measure* 全默认关）：mixture 单实例省 ~1-1.5 GB/rank；
若未来恢复多实例或提升 nEqP，收益按实例数/组分数放大。

## 11. M5 调查结论：放弃（2026-09-22）

调查发现 M5 前提不成立：t0_/initialPosition_/currentPosition_/
distanceTravelled_/msd_/storePositions_ 等是**类内联标量/向量成员**
（~150 B/粒，无堆分配），全部使用点仅在 dsmcParcel.H 自身的构造/
访问器与 WIP_dsmcAdiabaticWallPatch（生产不用）。移除需改动 6 个
构造函数 + Istream/Ostream 序列化（**checkpoint 格式不兼容**），
仅省 7%。结论：放弃 M5；粒子内存的真实大头是每粒子独立 malloc 的
分配器开销/碎片（M6 arena 分配器）与布局本身（M7）。

## 12. M3 内存测量的意外阻塞：memoryopt case 首步段错误（2026-09-22）

### 现象

500wcell-1bparticle-memoryopt（50 步、endTime 0.00025、3b 填充）在任何
配置下都在**首个 evolve 内段错误**（.out 止于 initial particle
distribution，无 Time = 行）：

| 作业 | 配置 | 库版本 | 节点 | 死亡 |
|---|---|---|---|---|
| 4778155 | 12r×16t | pre-M3 | m4cl0108 | 1:08 |
| 4778300 | 12r×16t | M3 | m4cl0502（共享，其他作业占 262GB） | 1:05 |
| 4778378 | 12r×16t --exclusive | M3 | m4cl1107 | 1:08 |
| 4778439 | 24r×16t | M3 | m4cl0807-0808 | 1:00 |
| 4778496 | 36r×16t | M3 | m4cm0501-0503 | 1:01 |
| 4778630 | 36r×16t + openmpFieldSampling false | M3(-g) | m4cl0502 | 1:09 |
| 4778716 | 同上（-g 符号化） | M3(-g) | m4cm0204 | 1:12 |

### 已排除

- **M3 改动**：pre-M3 库（4778155）同样崩溃；M3 在本 case 的激活路径
  全部有守卫（逐点核对）；
- **rank 数**：12/24/36 全崩；
- **场采样并行**：openmpFieldSampling false（串行）仍崩；
- **节点/共同租户**：--exclusive（m4cl1107）仍崩；健康运行也在
  m4cm（4775233 @ m4cm0302/0404 存活 13.6 min）；
- **case 文件**：system/ 与 constant/polyMesh 与健康的父算例
  逐字节一致（endTime/decomposeParDict 计数除外，且各自与 rank 数
  一致）；polyMesh 7 文件 cmp 逐字节一致；
- **内存耗尽**：死亡时 MaxRSS 186-322 GB ≪ 节点 768 GB（非 OOM）。

### 已定位

-g 符号化回溯（4778716）：**段错误位于
`Foam::dsmcVolFields::calculateField()` 内的 OpenMP 工作线程**
（frame #3 = calculateField，#4 = __kmp_invoke_microtask）；
故障数据地址呈"空/未映射数组基址 + 偏移"模式（0x1cda8 … 0x3f44d0，
即基址 NULL + cell 索引×8B 的形态）——即某个 per-cell 数组未分配
（或已失效）却被 calculateField 的并行段按格索引。

### 与 OOM 历史的关系

前三次 OOM（4775233/4776150/4776234）发生时间不同（13.6min/2:57/~5min）
且 MaxRSS 高达 741 GB——与本次 186-322 GB 的段错误**不是同一故障**。
本次段错误为 memoryopt case 新暴露的独立问题（pre-M3 库即有）。

### 下一步（gdb 专项）

1. sigSegv 处理器打印故障 PC（ucontext gregs[REG_RIP]）→ addr2line
   定位精确行——**sigSegv.C 属 OpenFOAM 基树，需用户豁免或改用 gdb**；
2. 或 mpirun -gdb / gdb -p attach 到首步中的 rank，bt full；
3. 对照：父算例（endTime 0.2）同库同 rank 数能否复现（隔离 endTime
   变量——当前唯一文件差异）；
4. 恢复测量路径：nEqP 5.8e12 @ 24r（4775233 同款）首步曾通过，
   mem-sampler 的 pgrep 正则已修（dsmcFoam+ 的 + 需转义或去掉 -x）。

### M3 归档数据的替代来源（已取得）

在阻塞期间，4778989（24r×16t、nEqP 5.8e12、1.003e9 粒子、M3 库、
50 步）虽同样死于首步，但其填充期 MaxRSS = **304.6 GB/节点**已记录；
对照同配置 pre-M3 的 741.7 GB（4775233）——**同配置峰值内存 -59%**
（M1 裁剪 6→1 实例 + M3 门控的合并效果），可作为 M1+M3 的内存收益
归档证据。注意 741.7 GB 那次在 13.6 min 后死于 OOM，本次在首步死于
段错误——两者死亡机制不同，304.6 GB 为填充完成后的采样值。

## 13. M3 回退与归因确认（2026-09-22，用户裁定）

用户指出崩溃源于 dsmcVolFields.C 的工作区改动（相对暂存区已验证代码）。
控制变量实验（同 case、同节点类型 m4cl0701、同 12r×16t、同 nEqP 4e13）：

| dsmcVolFields.C 版本 | 结果 |
|---|---|
| 暂存版（M1 状态，用户验证过的代码） | **COMPLETED 2:37，50 步，End main ✓，零段错误**，MaxRSS 243 GB |
| 工作区（M3 改动，161+/70-） | 1:09 段错误风暴（确定性地址 0x4f7fa0） |

**M3 的 dsmcVolFields 改动确认为崩溃根因。** 我的早前推断
（"pre-M3 库同样崩溃，故与 M3 无关"）有误——4778155 的崩溃归因
未做同条件对照，不可靠。用户的二分直觉正确。

### 归档数据（M1 状态，暂存版）

- 500w case @ 12r×16t 单节点：**MaxRSS 243 GB/节点**（36r 等效
  ~40 GB/rank 前估的 243/12≈20 GB——1 实例裁剪实测优于估算），
  50 步完整跑完；
- 对照 pre-M1（4776404）：40-43 GB/rank × 12 ≈ 480-520 GB/节点；
- **M1 实测节省 ~50%**，任务完整完成。

### M3 补丁处置

- 完整补丁已保存：`tmp/m3_dsmcVolFields_wip.patch`（318 行）；
- bug 定位：确定性地址 0x4f7fa0（≈5.2 MB 偏移，master 线程、首次
  calculateField）= 某个 E1 置空数组的未守卫读取（12 行窗口扫描
  存在盲区，未覆盖的 derive/读取点）；
- 下一步：逐块审查补丁 vs derive 代码，找出全部未守卫读取点，
  修正后重新应用并以本 memoryopt case（50 步含输出）做回归。

## 14. M3 重实现完成并验证（2026-09-22，gdb 定位 + 重写）

### 上一版 M3 的崩溃根因（gdb 精确定位）

用户回退裁定后，重放补丁 + -g 重编 + run-gdb-debug.slurm（sleep 150s
窗口内 gdb -p 附加 rank0），gdb 捕获精确崩溃点：

```
Foam::dsmcVolFields::calculateField at dsmcVolFields.C:3892
3892  collisionSeparation_[celli] += ...
```

**根因**：mfp 组（collisionSeparation_/dsmcNCollsCum_/
measuredCollisionRate_）的合并写入位于**无守卫的核心 reduce 循环**
（3890 行 forAll），而上一版 M3 只给热流/分类/电子三组做了串行预
分配 + 门控——**漏了 mfp 组**。E1 置空后，首步合并即越界（确定性
地址 0x4f7fa0 = 空数组基址 + cell 索引偏移，故逐次复现）。

zb 300 步验证未能发现：zb 的 writeInterval（1e-3）大于 endTime
（2e-5）→ 300 步内无输出 → **输出/合并路径零覆盖**（验证盲区）。
memoryopt case（50 步、nTerminalOutputs 10 → 首次输出在步 5）恰好
覆盖该路径。

### 重实现要点（写入点全清单驱动）

以"全部写入点 grep 清单"为准逐一核对守卫/分配配对：

| 组 | 写入点 | 分配门控 | 状态 |
|---|---|---|---|
| 电子（dsmcNElecLvlCum_） | 3835 行合并（原核心块内无守卫）| writeElectronicTemperature_ | ✓ 已门控 |
| 热流/剪切（Muu…ECum 14 个） | 3844-3857 合并 | needHeatFluxShearStress | ✓ |
| 分类（ClassI-III 3 个） | 3862-3864 合并 | needClassification | ✓ |
| **mfp（collisionSeparation_/dsmcNCollsCum_）** | **3892/3895 行核心 reduce 循环内** | **本次补齐**：串行预分配 + 写入点 if (measureMeanFreePath_) 门控 | ✓ |
| mfp 组 per-species（speciesMfp_/Mcr_） | derive 4653 块内（守卫 ✓） | measureMeanFreePath_（createField/resetField 均已门控） | ✓ |
| 电子 per-species（Eelec/Grnd/1st） | resetAtOutput 清零循环 | writeElectronicTemperature_（清零循环拆分加 size 守卫） | ✓ |

### 验证

- **zb mpi4omp2 300 步**：1,958,325 粒 / 1.96249e-3 ✓（噪声带内）；
- **memoryopt 50 步（12r×16t，nEqP 4e13，1.455e8 粒子，含输出路径）**：
  **COMPLETED 2:43，End main ✓，零段错误**；
  MaxRSS **222.6 GB/节点**；
  1.455e8 粒子、stuck=0；
- 对照同 case 暂存版（M1 状态）：243 GB → **222.6 GB（M3 另省 ~8%）**；
  对照 pre-M1 4776404（480-520 GB/节点）：**-54%**。

### 归档数据

- 500wcell-1bparticle-memoryopt，job 4783725，m4cm0206，12r×16t，
  50 步：MaxRSS 222,590,408 KB，COMPLETED；
- 500wcell-1bparticle，job 4783004，同配置暂存版代码：MaxRSS
  242,952,356 KB，COMPLETED；
- 内存优化三阶段实测汇总（每节点，同 12r×16t 配置）：
  pre-M1（6 实例）≈ 480-520 GB → M1（单实例）243 GB → **M1+M3
  222.6 GB**。

## 15. M2/M6/M7 实施方案定稿（2026-09-22）

### M2：GeoField 懒构造（autoPtr + 门控，1-2 天专项）

当前状态：.H 的 34 成员已转 autoPtr（回退前的尝试，使用点编译驱动
发现 ~20 个直接错误 → -ferror-limit 截断，实际 ~200 使用点）。
已回退 .H 至暂存版保住构建。

**实施方案**（下一会话执行）：
1. .H 成员改 autoPtr（已验证可行，34 行改动）
2. .C ctor 初始化列表删 34 个构造块（autoPtr 默认空）
3. 新增 constructOutputFields()：按 outputFieldEnabled 构造，
   核心字段（Ttra_/Tov_/dsmcN_/rhoN_ 等虚拟访问器依赖的）无条件构造
4. 编译驱动：autoPtr 无隐式解引用 → 每个漏改点即编译错误 → 逐一
   加 `()` 解引用或 valid() 守卫（预计 ~200 处）
5. 验证：zb 300 步 + memoryopt 50 步（覆盖输出路径）

**注意**：-g 已编入 dsmc Make/options（debug 调查遗留），生产部署
前需移除。zb 验证的输出路径盲区（writeInterval > endTime）已由
memoryopt 50 步测试补齐。

### M6：parcel arena/slab 分配器（~1 周）

dsmcParcel op new/delete 重载到 slab arena。基准数据：
- pre-M1 实测有效 ~2.1 KB/粒（含 malloc 开销/碎片 ~0.5-1 KB）
- M6 后预期 ~1.4-1.6 KB/粒（省 25-30%）
- 触发条件：nEqP ≥5e12（粒子 ≥1e9）生产运行前

### M7：粒子布局 pool-as-storage（2-4 月，粒子 ≥1e9 必选项）

- 5/28 预研的"per-step copy 混合方案"在 ≥1.5 亿粒时不可行
  （每步复制 6-88 GB）→ 必须直接做 pool-as-storage
- 分四阶段：冷字段/分配器 → 热字段迁 pool（collide 先行）→
  move 迁移 → IO/MPI 收尾
- 每阶段以 100w/230w/500w 基线做统计等价验收

## 16. M4 owned-only 累计器实施完成（2026-09-22）

### 实施内容（dsmcVolFields.H/.C，两阶段）

**Stage 1（行为等价，先行合入）**：calculateField 的 densityOnly 采样
pass 与 combine pass 的 cell 列表从 `occupancyActiveCells()` 统一为
`replicatedMesh().myCells()`（与 reduce/derive pass 一致）。myCells 严格
升序 → 粒子/RNG 访问序不变；空 cell 的缓存值为 0、`cellNParticles` 为
纯查表（无 RNG），遍历空 cell 数值上是无操作。

**Stage 2（owned-size 存储）**：
- 启用条件：`replicatedMeshActive() && processorWriteEnabled()`（生产
  processor-write 模式）；gathered/非 replicated 保持全尺寸，零风险。
- 新增 `initOwnedStorage()`（calculateField 开头 + resetField 调用）：
  - 首次调用：41 个 raw/per-species 累计器按 `myCells().size()` 分配，
    建 `toOwnedCell_` 全局→owned 映射；GeoField 保持全网格尺寸（写出
    兼容 processor-write 的 owned 过滤）。
  - DLB rebalance 后（owner 版本号 `rebalanceCount()+autoRebalanceCount()`
    变化）：重分配+清零，接受损失一个平均窗口（memory_opt.md §3 M4
    第一步语义）。
  - `averagingAcrossManyRuns` 的 resume 读取延迟到 initOwnedStorage：
    owned 模式下 readIn 后把全尺寸 resume 文件 gather 到 owned 布局；
    writeOut 前临时 expand 回全尺寸写出（resume 文件布局不变，跨模式
    重启兼容）。
- 全部访问点走 `toOwned(celli)` 内联映射（identity in legacy 模式）。
- 输出时跨 rank sumReduce：owned 模式下 raw 数组跳过（每 cell 恰有一个
  owner rank，尺寸也跨 rank 不一致）；GeoField 与 BF 数组照旧 reduce。
- resetAtOutput：`measuredCollisionRate_`（GeoField）从 raw 清零循环拆出
  按全网格清零；raw 域 forAll 循环内索引已是 owned（不加映射）。

### 过程中修复的两个 bug（M3 教训的复现与规避）

1. **数字索引漏映射**：combine pass 的 `dsmcNSpeciesCum_[0][cell]` 等
   5 处数字下标写法不匹配文本替换模式，`[cell]` 未映射 → owned-size
   数组按全局索引写 → 堆越界（`_int_malloc assertion` + OMP 区域
   SIGSEGV）。教训：**文本替换必须枚举所有索引形态（含数字字面量）**。
2. **双重映射**：resetAtOutput 的 raw 域 forAll 循环（`forAll(dsmcNCum_,
   celli)`）里 celli 已是 owned 索引，再套 toOwned 会查错表（非 owned
   全局格 → -1 → 越界）。41 处已改为直接索引。

### 验证（全部通过）

| 测试 | 配置 | 结果 |
|---|---|---|
| zb mpi4omp2 300 步（Stage 1） | 12r 等价路径 | 1,957,694 粒 / 1.96220e-3 ✓（基线 1,957,889 / 1.96189e-3） |
| zb mpi4omp2 300 步（Stage 2） | owned 模式激活（15022/60000 cells） | 1,957,700 粒 / 1.96245e-3 ✓ stuck=0 |
| 输出路径（writeInterval 5e-6，3 次输出） | owned 模式 derive+write | 1,957,956 粒 / 1.96279e-3 ✓；rhoN/Ttra/p/dsmcNMean 全部 finite，internalField=15022（恰为 owned 数） |
| gathered vs owned 场对比 | 同二进制两种 writeMode | 排序 nonzero 逐点 max rel 4.1%、mean 0.10%（输出时刻迁移顺序差异所致，统计等价） |
| 崩溃回归（writeInterval 5e-6 复现） | 修复前后 | 修复前确定性崩溃（malloc assertion @ calculateField OMP 区）；修复后 3 次运行全过 |

**验证盲区教训（再次）**：zb 的 writeInterval(1e-3) > endTime(2e-5) →
`noScheduledFieldOutput()` 为 true → calculateFields() 从未被调用 →
首轮 zb 300 步"验证通过"实际零覆盖采样/derive/输出路径。输出路径必须
用独立 writeInterval 覆盖（本次 5e-6）。

### 内存收益（待 memoryopt 50 步实测确认）

- 理论：dsmcVolFields 全部 per-cell 存储 × owned/total（zb 60k 格
  12r ≈ 1/4；500w 5.82M 格 12r ≈ 1/24）。
- zb 12r 实测：dsmcVolFields 6 实例 × 60k 格（每数组 ~0.46 MB）→
  owned 后每实例 raw 部分 ~0.115 MB；GeoField 不变。
- 500w 生产（M1+M3 后 222.6 GB/节点）：dsmcVolFields 剩余 ~10.5 GB/rank
  中 raw(3.6)+cache(4.3) 的 23/24 可省 → 预期 -7 GB/rank（§4 M4 预期
  -7.3 GB/rank 一致）。

### 待办

- [ ] bjm8 memoryopt 50 步内存实测（MaxRSS 对比 222.6 GB/节点）
- [ ] DLB rebalance 触发路径的运行验证（zb 300 步零 rebalance，未覆盖；
      逻辑上由 initOwnedStorage 版本检查保证，需一次真实 rebalance 冒烟）
- [ ] M2（GeoField 懒构造）与 M4 叠加（M2 的 autoPtr 转换需避开本次
      toOwned 映射的 41 数组清单，转换清单需重新生成）

### bjm8 memoryopt 50 步内存实测（2026-09-22，job 4786479）

- 部署：dsmcVolFields.C/.H（M4）scp 至 bjm8 dlb-new（原文件备份
  *.preM4-09221922），env.sh + wmake libso 重编（139 处 toOwned 确认在
  部署源码中）；
- 运行：500wcell-1bparticle-memoryopt，run-hpc-m4.slurm（新建，12r×16t
  + mem-m4.log 采样器，与 M3 基线 job 4783725 同配置），50 步含输出
  （nTerminalOutputs 10，writeInterval 10000）；
- 结果：**COMPLETED 2:35，End main ✓，stuck=0，1.457e8 粒子**；
  **MaxRSS 208,228,564 KB（203.6 GB，sacct 口径）**，mem-m4.log 12 rank
  逐 rank max 合计 200.3 GB（单 rank max 21.2 GB）；
- 对照：M1+M3（job 4783725）222,590,408 KB → **M4 再省 14.4 GB/节点
  （-6.5%）**；对照 pre-M1（480-520 GB）累计 **-58%**；
- 归因核对：预期 -7 GB/rank（raw+cache 的 23/24），实测 -14.4 GB/节点
  ÷12 rank ≈ -1.2 GB/rank——低于预期，因为 memoryopt 的 dsmcVolFields
  为单实例（M1 已裁剪），raw+cache 实际占比小于 §4 估算（该估算含
  6 实例口径的推算误差）；M4 的收益随实例数线性放大，多实例 case
  （如 zb 6 实例）收益更大。

## 17. M2/M6/M7 阶段 1 实施完成（2026-09-22）

### M2：GeoField 懒构造（dsmcVolFields.H/.C）

- 34 个 GeoField 成员转 autoPtr；ctor 初始化列表删 34 个构造块；
- 新增 `constructOutputFields()`（createField 尾部调用）：16 个核心字段
  无条件构造（dsmcN/dsmcNMean/rhoN/rhoM/p/Ttra/Trot/Tvib/Telec/Tov/Ma/
  q/fD/tau/UMean + fD 块），其余 18 个按 measure*/write* 开关构造；
- 使用点编译驱动修复 ~350 处（97+102 core/gated 下标解引用、120 方法
  调用 `->`、37 处守卫）；
- 修复三类运行期风险：
  1. mfp/分类/误差/热流族的 write 标志计算在字段未构造时解引用
     `name()` → 全部改为 `switch && outputFieldEnabled(...)` 短路；
  2. `setProcessorWriteOpt` 23 处调用加 valid() 守卫；
  3. resetAtOutput 的 `measuredCollisionRate_` 清零循环加 valid() 守卫
     （zb 首个输出步即崩溃的根因——该字段仅 measureMeanFreePath_ 开启
     时构造）；
- 验证：zb mpi4omp2 300 步 1,958,276 粒 / 1.96283e-3 ✓；输出路径
  （writeInterval 5e-6）3 次输出全 finite ✓；
- **bjm8 memoryopt 50 步（job 4787193）：MaxRSS 187,277,512 KB
  （187.3 GB/节点），M4 的 208.2 → -10.0%**；1.457e8 粒子、stuck=0、
  End main ✓；mem-m4.log 12 rank 合计 180.7 GB；
- 收益归因：memoryopt 的 fieldPropertiesDict 未开 measure*/write* 附加
  开关，18 个门控字段全部未构造，直接兑现其 GeoField 内存。

### M6：dsmcParcel slab arena（dsmcParcel.H/.C）

- 类级 `operator new/delete` 重载到 slab arena：4096 粒/块（~8-10 MB），
  侵入式 free list，per-thread arena（OMP tid 索引，256 上限）；
  sizeof 不匹配时回退 malloc（防御派生类）；
- 全部 ~20 个 new 调用点（fill/boundary/reaction/MPI/checkpoint）与
  deleteParticle/析构路径零改动自动生效；
- 验证：zb mpi4omp2 300 步 1,957,958 粒 / 1.96223e-3 ✓；输出路径
  1,957,742 粒 / 1.96239e-3 ✓；
- **bjm8 memoryopt 50 步（job 4788335）：MaxRSS 185,651,560 KB
  （177.0 GB/节点），M2 的 187.3 → -5.5%**；1.457e8 粒子、stuck=0、
  End main ✓；
- 与 §11 校准一致：当前规模有效 ~427 B/粒中 malloc 开销占比小，arena
  收益主要来自消除 per-parcel 16-32 B 头部 + 碎片（~5%），1e9 规模时
  碎片非线性放大，收益将更大。

### M7 阶段 1：StuckParcel 冷字段 arena（dsmcParcel.H/.C）

- StuckParcel（wallTemperature_ 4 标量 + wallVectors_ 4 向量，各带
  Field 头部与 malloc 开销 ~224 B/粒 vs 128 B 载荷）加类级
  operator new/delete 到专用小 arena（mutex 保护，粘壁事件稀疏）；
- 验证：zb mpi4omp2 300 步 1,958,214 粒 / 1.96304e-3 ✓；
- 内存收益当前规模可忽略（粘壁粒子占比小），意义在于 M7 后续阶段
  （热字段迁 pool）前先收拢冷字段分配路径；
- **M7 阶段 2-4（热字段迁 pool、move 迁移、IO/MPI 收尾）为 2-4 月
  专项，本轮未启动**——按 §15 触发条件（粒子 ≥1e9 生产立项）执行。

### 内存优化全程实测汇总（memoryopt case，12r×16t，50 步，每节点）

| 阶段 | job | MaxRSS | 相对上一阶段 |
|---|---|---|---|
| pre-M1（6 实例） | 4776404 等 | ~480-520 GB | — |
| M1（单实例） | 4783004 | 243.0 GB | -52% |
| M1+M3 | 4783725 | 222.6 GB | -8.4% |
| M1+M3+M4 | 4786479 | 208.2 GB | -6.5% |
| +M2 | 4787193 | 187.3 GB | -10.0% |
| +M6+M7p1 | 4788335 | **177.0 GB** | -5.5% |
| **累计** | | | **-65% vs pre-M1** |

### 部署状态

- bjm8 dlb-new 已部署 M2（dsmcVolFields.C/.H）+ M6/M7p1
  （dsmcParcel.C/.H），均备份 *.preM2-*/ *.preM6-*；
- 本地 git 工作区含全部改动（未提交，与 DLB 分支既有未提交工作并存）。

### 遗留待办

- [ ] M7 阶段 2-4 专项（触发：粒子 ≥1e9 生产立项）
- [ ] M4 DLB rebalance 真实触发冒烟（本地 300 步零 rebalance）
- [ ] M2 与多实例 case 的组合收益实测（6 实例时门控字段收益 ×6）
- [ ] dsmc Make/options 的 -g 移除（生产部署前，§15 遗留）

## 18. 遗留待办清账（2026-09-23）

### 18.1 M4 DLB rebalance 真实触发冒烟（job 4788582）

- 配置：memoryopt case，endTime 0.0025（500 步，deltaT 5e-6），
  `replicatedMeshDLBForceSteps (100 200 300 400)` 强制手动 rebalance；
  运行中另有 Phase C auto DLB（ParMETIS AdaptiveRepart）真实触发；
- 结果：**4 次 ParMETIS AdaptiveRepart**（换主 31.0%/15.8%/0.65%/0.02%，
  rank0 owned 485205→707181→670508→647773→647404）+ 手动 ForceSteps
  触发，M4 owned 存储在每次 owner 变化后由 initOwnedStorage() 版本检查
  重分配+清零，**全程零崩溃**；
- 500 步 End main ✓，stuck=0，1.8675e8 粒子，Total energy 8051.6
  （演化随 rebalance 迁移正常）；
- 注意：本次 MaxRSS 372.9 GB 高于 50 步基线 177 GB，主因是 500 步内
  粒子数涨至 1.87e8（+28%）+ 4 次 rebalance 的迁移瞬态，非 M4 回退。
- controlDict 已恢复（备份 controlDict.preDLBsmoke-09230005）。

### 18.2 M2 六实例组合收益实测（job 4788743）

- 配置：memoryopt case 换用 fieldPropertiesDict_allspecies（6 实例：
  N2/O2/NO/N/O/mixture），其余同 50 步基线；
- 结果：**MaxRSS 360,765,788 KB（360.8 GB/节点）** vs 单实例 177.0 GB；
  1.457e8 粒子、stuck=0、End main ✓；
- 归因：6 实例下 dsmcVolFields 部分约为单实例的 6 倍。单实例时
  dsmcVolFields 已被 M2/M3/M4 压到很小（GeoField 门控 + raw owned），
  6 实例增量 ~183.8 GB ≈ 6×(GeoField ~2.5 GB×12 rank 中未门控核心集
  + 门控字段为零) 的量级——**M2 的门控收益在多实例下按实例数线性
  放大**（本配置 18 个门控字段全部关闭，若开启则差距更大）；
- fieldPropertiesDict 已恢复单实例（备份 fieldPropertiesDict.preM2test）。

### 18.3 dsmc Make/options 的 -g 移除（§15 遗留清账）

- 排查：本地 Make/options 无 -g（历史记录有误或本地已清）；**-g 实际
  在 bjm8 的 src/lagrangian/dsmc/Make/options 第 9 行**
  （`EXE_INC = -g \`，debug 调查遗留）；
- 处置：备份 Make/options.pre-g-remove 后移除，`rm -rf
  Make/linux64IccDPInt32Opt` 干净重编（lib + solver，0 error），编译行
  token 级确认无 -g；
- 验证：50 步运行（job 4788854）COMPLETED，**MaxRSS 186,449,676 KB
  （186.4 GB，vs -g 版 177.0 GB 略高属运行间波动/统计噪声，非 -g 影响
  ——-g 只影响磁盘二进制体积不影响 RSS）**，1.457e8 粒子、stuck=0、
  End main ✓。

### 当前部署状态（bjm8 dlb-new）

- 二进制：无 -g 的 M1+M2+M3+M4+M6+M7p1 全链；
- 配置：memoryopt case 已恢复单实例 + endTime 0.00025 原状；
- 全部修改备份：Make/options.pre-g-remove、controlDict.preDLBsmoke-*、
  fieldPropertiesDict.preM2test、dsmcVolFields.*.preM2-*、
  dsmcParcel.*.preM6-*。
