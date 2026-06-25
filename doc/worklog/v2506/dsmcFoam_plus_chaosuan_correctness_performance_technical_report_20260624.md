# dsmcFoam+ 超算正确性与性能回修技术报告 - OFv1706 hyStrath_dlb

日期：2026-06-24

工作目录：

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

参考实现：

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
```

本文是 `dsmcFoam+` 在 OFv1706 `hyStrath_dlb` 分支进入超算运行后的阶段性
技术报告。报告范围覆盖：

- 超算环境下 ParMETIS/METIS 构建、ABI 和运行时绑定问题；
- raw-MPI replicated mesh 与 mixed MPI+OpenMP 在超算上的崩溃、卡死和正确性修复；
- replicated mesh 输出正确性，从 rank0 gather-write 过渡到真正 processor-write；
- `replicatedMeshMigrateInterval=1` 作为当前正确性基线的确认；
- FastRNG scope 修正；
- post/output 热点定位、`dsmcVolFields` 共享 cache、输出字段白名单和 processor cloud
  可选字段懒分配；
- 最新超算 64 核完整 8750 步性能结果。

本报告是以下报告之后的超算阶段 follow-up：

```text
doc/worklog/v2506/dsmcFoam_plus_omp_technical_report_20260606.md
doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_dlb_technical_report_20260608.md
doc/worklog/v2506/dsmcFoam_plus_mpi_omp_mixed_technical_report_20260610.md
doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_further_optimization_technical_report_20260617.md
```

详细工作日志位于：

```text
doc/worklog/v2506/detail_chaosuan/chaosuan_correctness_performance_worklog_20260621.md
```

## 1. 结论摘要

本轮超算阶段的结论可以压缩为八点：

1. 超算 early failure 的第一层问题不是 DSMC 算法本身，而是 ParMETIS/METIS 构建、
   ABI 和运行时动态库绑定。修正后，`mpi64` 和 mixed replicated mesh 能进入长时间运行。
2. replicated mesh 在超算上的写出问题已经从 rank0 gather-write 修到真正
   processor-write：每个 MPI rank 写自己的 `processorN/` mesh/addressing、
   volFields 和 lagrangian cloud，标准 `reconstructPar` 可恢复 root 结果。
3. `replicatedMeshMigrateInterval=10` 与当前 owned-cell collision 架构不自洽。
   跨 owner 粒子会在旧 rank 滞留最多 9 步，导致 bulk flow collision/sampling 和
   wall flux 都可能偏离。因此当前正确性基线必须是：

```text
replicatedMeshMigrateInterval 1;
```

4. 早期正确性修复曾带来性能退化，主要来自两个保守化：
   - migration 后过度清空 `moveOrderedParcels_`，破坏 build occupancy 快路径；
   - OpenMP move lazy mesh data 预热被计入 move timer。
   后续已通过 ordered-cache 完整性 guard 和 move timer 外预热修回大部分问题。
5. FastRNG scope 已修正。`collisionFastRng true` 现在覆盖 OpenMP collision、
   mixed MPI+OpenMP collision、pure MPI replicated single-thread collision，并通过
   `cloud_.setCollisionRngContext(...)` 使 collision/reaction 热路径使用同一 fast RNG
   context。
6. post/output 的主要开销不是文件写入，而是 `dsmcVolFields::calculateField()` 的
   宏观场计算。本轮已加入共享 sample cache、输出字段白名单和 processor cloud
   可选字段懒分配，降低重复遍历和不必要输出。
7. 最新超算 64 核 8750 步完整系列中，当前最优配置是：

```text
MPI32 x OMP2
ClockTime = 1680 s
```

   次优是 `MPI16 x OMP4`，`ClockTime = 1759 s`。`MPI64` 可以正常跑完，但
   `ClockTime = 2291 s`，不是性能最优；`OMP64` 为 `3140 s`，在该超算 case
   也不是最优。
8. 后续优化不应再回到 `migrate=10`，而应在 `migrate=1` 正确性基线下继续优化
   每步 migration、occupancy 构建、collision rank imbalance 和 post 字段输出量。

## 2. 环境、算例和复现实验规则

### 2.1 本地构建环境

激活环境：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
```

编译：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
```

主要产物：

```text
platforms/linux64IccDPInt32Opt/lib/libdsmcFoam+.so
platforms/linux64IccDPInt32Opt/bin/dsmcFoam+
platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+
```

### 2.2 超算环境约束

超算上必须保证三件事一致：

- OpenFOAM 编译器和 MPI wrapper 来自同一套 oneAPI/OpenFOAM 环境；
- ParMETIS/METIS/GKlib 的头文件和动态库来自同一个安装前缀；
- OpenFOAM 编译时和运行时看到的 `metis.h` / `parmetis.h` ABI 一致。

现场确认过的 ParMETIS/METIS 头文件宽度：

```text
IDXTYPEWIDTH  = 32
REALTYPEWIDTH = 64
```

运行时需要通过 `ldd` 确认：

```text
libparmetis.so
libmetis.so
libGKlib.so
```

都来自用户编译的 `PARMETIS_DIR`，而不是系统或第三方残留库。

### 2.3 主要超算 case

超算完整算例：

```text
results/palphd3.3.1react-m8-fixoutput
```

物理口径：

- Palharini PhD 3.3.1 reactive cylinder case；
- `endTime 0.0035`；
- `deltaT 4e-07`；
- 8750 iterations；
- 总核数 64；
- `replicatedMeshMigrateInterval=1`；
- mixed/pure MPI 使用 raw-MPI replicated mesh，不使用 OpenFOAM decomposed
  `-parallel` processor mesh。

### 2.4 本地回归 case

本地快速回归主要使用：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react
```

严格 `MPI4xOMP2` 运行口径：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
cd run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi4omp2
export OMP_NUM_THREADS=2
export OMP_PROC_BIND=close
export OMP_PLACES=cores
/usr/bin/time -p mpirun -np 4 "$FOAM_USER_APPBIN/dsmcFoam+"
```

该口径必须满足：

```text
Replicated mesh: initialized with 4 MPI ranks
OpenMP max threads = 2
Total Iterations = 300
End main
```

不能使用 OpenFOAM `-parallel`。

## 3. 超算现场正确性问题与修复

### 3.1 ParMETIS/METIS 初始化异常

现场问题：

- `mpi64` 初始化时触发 `Floating point exception`；
- 调用栈进入 `METIS_PartGraphKway()` /
  `libmetis__SetupKWayBalMultipliers()` /
  `dsmcReplicatedMesh::computeCellOwnerScotch()`；
- 相同网格、相同代码本地可运行，超算失败。

处理过程：

- 写最小 MPI/METIS 验证程序，只允许用户修改 `PARMETIS_DIR`；
- 自动推导 OpenFOAM `owner/neighbour` 的 cell 数，避免固定 cell 数错误；
- 重写 ParMETIS 构建脚本，从 OpenFOAM 编译环境抓取 `icc/icpc/mpiicc`；
- 避免直接 `source` 构建脚本中断超算连接，改为普通 `bash` 子 shell 执行；
- 用 `ldd` 核对运行时库绑定。

结论：

- 该问题本质是超算环境 ABI 和动态库绑定问题；
- 只有当 ParMETIS 头文件、库文件、OpenFOAM 编译环境和运行时 `LD_LIBRARY_PATH`
  一致后，后续 DLB 与 replicated mesh 问题才有可诊断基础。

### 3.2 ParMETIS `tpwgts` 参数错误

现场问题：

```text
PARMETIS ERROR: The sum of tpwgts for constraint #0 is not 1.0
```

触发路径：

```text
Phase C ParMETIS: AdaptiveRepart (..., 1 constraints, ubvec=1.05)
ParMETIS_V3_AdaptiveRepart(...)
```

判断：

- `tpwgts` 按 `real_t` 传入；
- 若编译时头文件和运行时库对 `real_t` 的宽度理解不同，ParMETIS 会按错误内存布局解释
  `1/nparts`；
- `REALTYPEWIDTH=64` 要同时作用于编译和运行。

结果：

- ParMETIS 重新构建和绑定后，此类初始化期错误消失；
- `mpi64` 能进入完整运行阶段。

### 3.3 OpenMP/mixed 崩溃

现场问题包括：

- `OMP64` 在 `dsmcCloud::buildCellOccupancy()` 附近 segfault；
- `MPI8+OMP8` 在 initial particle distribution 之后或 SAR/DLB 触发附近 segfault/卡死；
- pure MPI 修复后 mixed 反而暴露新的崩溃。

根因分为两类：

- OpenMP move region 内部触发 OpenFOAM mesh lazy data 并发初始化；
- migration 后 `moveOrderedParcels_` 与真实 cloud size 不一致或生命周期不完整。

修复：

- 新增 `Cloud<ParticleType>::prepareOpenMPMoveMeshData()`，提前预热 move 必需 mesh data；
- 将预热移到 `dsmcCloud::evolve_moveAndCollide()` move timer 之前；
- `autoMap()` 后重置预热标志；
- migration 后只在 `kept.size() == cloud_.size()` 时恢复 ordered cache，否则清空；
- `migrateFinish()` append received batch 时保持 occupancy cache 失效。

效果：

- mixed/OMP 崩溃消失；
- move timer 不再被一次性 lazy data 预热污染；
- build occupancy 快路径可恢复使用。

### 3.4 DLB trigger 一致性

为避免某些 rank 进入 ParMETIS、其他 rank 不进入 collective，`autoRebalance()` 增加了
跨 rank trigger 一致性保护：

- local trigger 后通过 `MPI_Allreduce` 汇总；
- 若 rank 间不一致，打印 mismatch；
- 使用全局 OR 结果进入 DLB。

同时 DLB decision timer 从 local wall 改为全局最大 wall，以反映关键路径。

该修改偏正确性和诊断，会带来少量 collective 开销，但避免了最危险的 rank collective
顺序不一致。

## 4. replicated mesh 输出修复

### 4.1 rank0 gather-write 阶段

早期 replicated mesh 输出问题：

- 多 rank 持有完整网格，但只拥有部分 parcels；
- 普通 OpenFOAM write 会让多 rank 写同一套 root cloud；
- cloud 字段和 parcel 数可能不匹配；
- `foamToEnsight` 后流场不对；
- write 后粒子 ownership 状态可能被破坏。

第一阶段修复采用 rank0 gather-write：

1. `gatherParcelsToRank0()`；
2. rank0 暂时关闭 cloud 自动写；
3. `runTime.write()` 写非 cloud 对象；
4. rank0 `writeGatheredCloudOnRank0()` 写 gathered cloud；
5. 写 `cellOwner` 场；
6. 输出后 `migrateParticlesByCellOwner()` 恢复计算态。

该阶段能修复写出正确性，但超算大粒子数 case 中 rank0 write 成本达到约
`20+ s/output`，不适合作为最终大规模输出路径。

### 4.2 真正 processor-write

后续实现的 processor-write 闭环：

- `replicatedMeshWriteMode processor` 下不再由 rank0 串行写 root 全局场/root cloud；
- 每个 MPI rank 写自己的 `processorN/` mesh/addressing、volFields 和 lagrangian cloud；
- processor mesh/addressing 按 owner 分区版本复用，避免每个 output time 重写；
- DLB 后 processor patch list 固定为 `nProcs-1` 个 processor patch，非邻居 rank
  使用 0-face patch，保证 `reconstructPar` 跨多个 output time 稳定读取；
- 标准 `reconstructPar` 可重构 root volFields 和 lagrangian cloud。

验证：

- 2-rank smoke 中，root 粒子数等于 processor 粒子数之和；
- `mpi4omp2` DLB 验证完成 1000 steps、500-step output、250-step forced DLB；
- 4 次 DLB 和 2 次 processor output 均成功；
- `3.320025e-05` 与 `6.64005e-05` 两个 output time 均可 `reconstructPar`。

生产使用注意：

- processor-write 运行阶段不生成 root time 场；
- `foamToEnsight` 或 root case 后处理必须在 `reconstructPar` 之后；
- no-write compute、processor-write、reconstructPar、foamToEnsight 必须分开计时。

## 5. `replicatedMeshMigrateInterval=1` 正确性基线

### 5.1 问题现象

在 `replicatedMeshMigrateInterval=10` 时，超算结果出现：

- 通信边界处流场不连续；
- 圆柱壁面附近 `wallHeatFlux` 偏高；
- 部分 mixed 结果与 OMP64 对照不一致。

典型对比：

| quantity | MPI16xOMP4, migrate=10 | OMP64 | difference |
|---|---:|---:|---:|
| `p_mixture` | 307.50 | 301.84 | +1.9% |
| `rhoN_mixture` | 1.866e22 | 1.902e22 | -1.9% |
| `dsmcNMean_mixture` | 125.54 | 132.32 | -5.1% |
| `wallHeatFlux_mixture` | 92501.56 | 69468.62 | +33.2% |

### 5.2 根因

当前 replicated mesh collision 是 owner-cell 模式：

- move 后粒子可能跨到另一个 rank owner cell；
- migration 只在 `stepCounter % replicatedMeshMigrateInterval == 0` 时执行；
- collision 只遍历本 rank `occupancyOwnedCollisionCells()`。

当 `migrateInterval=10` 时，跨 owner 粒子可在旧 rank 滞留最多 9 步：

```text
particle moves from rank A owned cell to rank B owned cell
particle still physically resides on rank A until next migration

rank A:
  can see the parcel in a B-owned cell
  but collision only loops over A-owned cells
  -> parcel does not collide in B-owned cell

rank B:
  owns B-owned cell
  but does not have parcel yet
  -> parcel also missing from rank B collision
```

因此这不是单纯 wall flux 归属问题，而是 bulk flow collision/reaction/sampling 的
owner 一致性问题。

### 5.3 结论

当前实现必须保持：

```text
replicatedMeshMigrateInterval 1;
```

评估过的替代方案：

| direction | correctness | implementation risk | conclusion |
|---|---|---|---|
| `migrationInterval=1` | 高 | 低 | 当前正确性基线 |
| 只修 wall flux face-owner 归属 | 只修壁面事件 | 中 | 不足以修 bulk flow |
| 非 owner active cell 也参与 collision | 不严格 | 中 | 易造成同一物理 cell 粒子分裂碰撞 |
| ghost/halo parcel collision | 理论可行 | 高 | 接近重设 DSMC 跨分区 collision |
| 每步迁移但优化 migration 成本 | 高 | 中 | 后续推荐方向 |

后续性能优化应在 `migrate=1` 基线下优化 migration 成本，而不是恢复 `migrate=10`。

## 6. FastRNG scope 修正

### 6.1 修正内容

修正后，`collisionFastRng true` 覆盖：

- OpenMP collision 线程循环；
- mixed MPI+OpenMP collision；
- pure MPI replicated mesh single-thread collision；
- collision/reaction 内部通过 `cloud_.setCollisionRngContext(...)` 进入同一个
  fast RNG context。

关键源码：

```text
src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.H
```

种子设计：

```text
OpenMP/mixed: rank + thread + timeIndex
pure MPI single-thread: rank + timeIndex
```

### 6.2 before/after 对比

超算 `MPI32xOMP2` 单次对比：

| case | ClockTime [s] | move+collide [s] | move [s] | buildCellOccupancy [s] | collision [s] | migration [s] | DLB 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| before FastRNG scope fix | 1569 | 1251.37 | 716.93 | 123.53 | 344.43 | 220.61 | 13 |
| after FastRNG scope fix | 1680 | 1370.49 | 756.26 | 126.55 | 399.90 | 196.28 | 11 |

该差异不能解释为“DLB 次数增加”，因为修正后 DLB 次数从 13 降到 11，
migration wall 也从 220.61 s 降到 196.28 s。增加主要体现在 move 和 collision。

判断：

- FastRNG scope 修正改变了随机序列；
- 随机序列变化会改变粒子分布、候选对、反应路径和局部负载；
- DSMC 随机程序不应跨 RNG 修正版本期待 bitwise 复现旧性能日志；
- `collisionFastRng false` 会回到原始 RNG 热路径，通常更慢，不建议作为性能配置。

## 7. post/output 优化

### 7.1 热点定位

profile 表明，`post fields/output` 的主要开销不是 field write，而是：

```text
dsmcVolFields::calculateField()
```

即宏观场计算。文件写出可优化，但不是最主要的 post 成本来源。

### 7.2 `dsmcVolFields` 共享 sample cache

优化内容：

- 在单次 field calculation 中先构建 per-cell/per-species 共享累积量；
- 混合场 `rhoN/rhoM/p/Ttra/Trot/Tvib/U` 等从共享 cache 派生；
- 减少多场重复遍历 parcel；
- 预计算 type mass、rotational zeta、vibrational quantum；
- 对 active cell 局部 reset，避免每次清空全域大数组。

相关源码：

```text
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H
```

### 7.3 输出字段白名单

`dsmcVolFieldsProperties` 支持：

```text
writeFields (p Ttra U wallHeatFlux);
```

语义：

- 不配置时保持旧行为，写全部字段；
- 配置后只写指定宏观场或字段组；
- 只控制 volFields 输出数量，不改变求解器物理推进。

### 7.4 processor cloud 可选字段懒分配

processor cloud 输出中，以下字段改为按需分配和写出：

```text
radialWeight
ERot
ELevel
stuckToWall
wallTemperature
wallVectors
isTracked
inPatchId
tracerInitialTime
tracerInitialPosition
tracerCurrentPosition
tracerDistanceTravelled
```

实现方式：

- 本 rank 先扫描自身 parcels 判断是否需要；
- `MPI_Allreduce(MPI_MAX)` 得到全局是否需要；
- 只有全局需要时才分配对应 `IOField` 并写出；
- 基础字段仍保留：

```text
position
U
vibLevel
typeId
newParcel
classification
```

新增控制项：

```text
replicatedMeshProcessorWriteCloud true;
replicatedMeshProcessorWriteCloud false;
```

默认 `true`。若设为 `false`，只写 processor mesh/volFields，跳过 lagrangian cloud。
这只适合纯宏观场输出或诊断 post/write 成本；若后续需要 reconstruct lagrangian cloud，
不能关闭。

### 7.5 本地验证

10-step output-filter/cloud-skip 验证：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi4omp2/log.codex_mpi4omp2_10step_outputfilter_cloudskip_np4_omp2_foreground_20260624_155259
```

结果：

- `writeFields (p Ttra U wallHeatFlux)`；
- `replicatedMeshProcessorWriteCloud false`；
- 写出 `scalar=19, vector=6`；
- 跳过 `processor0 lagrangian cloud`；
- `processor fields/cloud write = 0.399130284 s`；
- `post field calculate = 0.439340028 s`；
- `post field write = 0.001786365 s`；
- 10 steps 正常 `End main`。

300-step no-write 回归：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi4omp2/log.codex_mpi4omp2_300step_outputfilter_nowrite_np4_omp2_foreground_20260624_154408
```

结果：

| metric | value |
|---|---:|
| mode | 4 MPI x 2 OpenMP |
| iterations | 300 |
| real | 68.23 s |
| move+collide wall | 67.85252556 s |
| move only | 46.46549126 s |
| buildCellOccupancy | 3.738757474 s |
| collision phase | 10.53373378 s |
| migration wall time | 20.01529794 s |
| DLB rebalances | 2 |
| post field calculate | 0 |
| post field write | 0 |

该结果证明 no-write compute 口径没有引入 post/output 开销。

## 8. 最新超算 64 核性能结果

### 8.1 结果表

结果目录：

```text
results/palphd3.3.1react-m8-fixoutput
```

口径：

- 8750 iterations；
- 64 cores；
- `replicatedMeshMigrateInterval=1`；
- mixed/pure MPI 使用 processor-write；
- `ClockTime` 取最终 `Stage 1.0`；
- profile wall 取最终 `DSMC solver profile summary`。

| case | log | ClockTime [s] | move+collide [s] | move [s] | buildCellOccupancy [s] | collision [s] | migration [s] | DLB 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP64 | `omp64_3637765.out` | 3140 | 1786.35 | 892.42 | 421.54 | 430.05 | - | - |
| MPI2xOMP32 | `mpi2omp32_3637766.out` | 3096 | 2124.54 | 810.85 | 190.96 | 414.33 | 828.14 | 43 |
| MPI4xOMP16 | `mpi4omp16_3637767.out` | 2390 | 1731.05 | 814.76 | 133.05 | 382.48 | 502.73 | 30 |
| MPI8xOMP8 | `mpi8omp8_3637768.out` | 1923 | 1420.56 | 755.70 | 111.97 | 336.05 | 272.75 | 14 |
| MPI16xOMP4 | `mpi16omp4_3637769.out` | 1759 | 1417.55 | 760.07 | 114.38 | 381.68 | 249.56 | 18 |
| MPI32xOMP2 | `mpi32omp2_3637770.out` | 1680 | 1370.49 | 756.26 | 126.55 | 399.90 | 196.28 | 11 |
| MPI64 | `mpi64_3637771.out` | 2291 | 1882.19 | 930.64 | 417.37 | 401.50 | 386.18 | 42 |

### 8.2 性能解释

当前排序：

```text
MPI32xOMP2 1680 s
MPI16xOMP4 1759 s
MPI8xOMP8  1923 s
MPI64      2291 s
MPI4xOMP16 2390 s
MPI2xOMP32 3096 s
OMP64      3140 s
```

关键判断：

- `MPI32xOMP2` 是当前完整超算系列最优；
- `MPI16xOMP4` 接近最优，可作为 repeat3 对照；
- `MPI8xOMP8` collision 最低，但总时间不最优，说明总性能受 migration/DLB/move
  共同限制；
- `MPI64` 可以跑完，但 buildCellOccupancy 和 migration 明显偏高；
- `OMP64` 不再是最佳，说明单进程 64 线程在该 case 上受 OpenMP 调度、内存访问、
  build/move 和 post 阶段限制。

## 9. 当前推荐配置

超算 full-run 生产建议：

```text
replicatedMesh true;
replicatedMeshMigrateInterval 1;
replicatedMeshDecompMethod metis;
replicatedMeshAutoDLB true;
replicatedMeshDLBDualConstraint false;
replicatedMeshDLBAlpha 1;
replicatedMeshDLBAdaptiveAlpha false;
collisionFastRng true;
replicatedMeshWriteMode processor;
replicatedMeshProcessorWriteTimeMesh false;
```

64 核优先配置：

```text
MPI32 x OMP2
```

候选备选：

```text
MPI16 x OMP4
```

输出建议：

- 需要完整 lagrangian 重构时，保持 `replicatedMeshProcessorWriteCloud true` 或不写该项；
- 只需要宏观场时，可以临时设为 `false`，并用 `writeFields (...)` 限制 volFields；
- 正式性能对比继续区分 no-write compute、processor-write、reconstructPar、
  foamToEnsight。

## 10. 剩余问题和后续方向

### 10.1 必须继续保持的边界

后续不能把以下口径混在一起：

- no-write compute；
- write-on run；
- processor-write；
- reconstructPar；
- foamToEnsight；
- `migrate=1` 正确性基线；
- `migrate=10` 历史错误/诊断口径。

### 10.2 需要 repeat 的性能结论

当前超算 64 核最优来自单次完整 8750-step 系列。建议下一步对：

```text
MPI32xOMP2
MPI16xOMP4
MPI8xOMP8
```

做 repeat3，以确认 `MPI32xOMP2` 对 `MPI16xOMP4` 的优势是否稳定。

### 10.3 继续优化方向

优先级：

1. 每步 migration 成本：
   - pack/unpack；
   - size exchange；
   - received parcel apply；
   - `updateParticleCounts()`；
   - 只迁移跨 owner 粒子的扫描成本。
2. `buildCellOccupancy()`：
   - pure MPI64 中该项达到 `417.37 s`，是 MPI64 明显慢于 MPI32xOMP2 的主要原因之一；
   - 继续检查 ordered occupancy、active cell reset 和 replicated owner filtering。
3. collision rank imbalance：
   - 继续区分 local loop、reduce wait、sigmaBC、reaction；
   - 不应只用总粒子数判断 collision 负载。
4. post/output：
   - 用 `writeFields` 控制变量数量；
   - 用 `sampleInterval` 降低后处理频率；
   - 只在需要时写 lagrangian cloud。
5. wall heat flux：
   - 保持 `migrate=1`；
   - 统一 processor reconstruct 流程；
   - 单独检查采样面积、面 owner 事件归属和采样步数；
   - 不把热流偏差和性能口径混在同一组结论里。

## 11. 证据路径

主工作日志：

```text
doc/worklog/v2506/detail_chaosuan/chaosuan_correctness_performance_worklog_20260621.md
```

ParMETIS 最小验证和重构脚本：

```text
doc/worklog/v2506/detail_mpi/metis_minimal_validation_20260618
doc/worklog/v2506/detail_mpi/parmetis_rebuild_20260619/makeParMETIS.sh
```

processor-write 验证：

```text
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/processor_write_full_20260621.md
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/log.true_parallel_final_smoke_run
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/log.true_parallel_final_smoke_reconstruct
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/log.zb_mpi4omp2_dlb250_write500_1000step_fixedpatches_run
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/log.zb_mpi4omp2_dlb250_write500_1000step_fixedpatches_reconstruct
```

超算 64 核系列：

```text
results/palphd3.3.1react-m8-fixoutput/omp64_3637765.out
results/palphd3.3.1react-m8-fixoutput/mpi2omp32_3637766.out
results/palphd3.3.1react-m8-fixoutput/mpi4omp16_3637767.out
results/palphd3.3.1react-m8-fixoutput/mpi8omp8_3637768.out
results/palphd3.3.1react-m8-fixoutput/mpi16omp4_3637769.out
results/palphd3.3.1react-m8-fixoutput/mpi32omp2_3637770.out
results/palphd3.3.1react-m8-fixoutput/mpi64_3637771.out
results/palphd3.3.1react-m8-fixoutput/mpi32omp2beforefastrng_3634454.out
```

FastRNG 与 post/output 构建：

```text
doc/worklog/v2506/detail_mix/build_fastrng_scope_20260622.log
doc/worklog/v2506/detail_mix/build_fastrng_scope_final_20260622.log
doc/worklog/v2506/detail_mix/build_output_filter_cloudskip_20260624_154051.log
doc/worklog/v2506/detail_mix/build_output_filter_cloudskip_solver_20260624_154110.log
```

本地 post/output 验证：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi4omp2/log.codex_mpi4omp2_300step_outputfilter_nowrite_np4_omp2_foreground_20260624_154408
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi4omp2/log.codex_mpi4omp2_10step_outputfilter_cloudskip_np4_omp2_foreground_20260624_155259
```

核心源码：

```text
applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C
src/lagrangian/basic/Cloud/Cloud.C
src/lagrangian/basic/Cloud/Cloud.H
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.H
src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H
```
