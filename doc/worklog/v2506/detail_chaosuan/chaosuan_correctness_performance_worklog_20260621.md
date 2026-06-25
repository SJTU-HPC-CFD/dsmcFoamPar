# 超算正确性与性能回修工作日志

日期：2026-06-21

工作目录：

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

相关提交与当前状态：

```text
HEAD = 6e1e671 V-1.7: finish output error, run correctly in Beijingchaosuan. But make performance worse. Will optimize later.
```

当前 `HEAD` 已包含超算正确性修复和部分性能诊断日志；其后还有一组未提交的
性能回修补丁，主要修改：

```text
src/lagrangian/basic/Cloud/Cloud.C
src/lagrangian/basic/Cloud/Cloud.H
src/lagrangian/basic/Cloud/CloudIO.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.H
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C
```

本文记录“开始超算运行之后”围绕正确性和性能展开的工作，不替代此前三份
阶段报告：

```text
doc/worklog/v2506/dsmcFoam_plus_omp_technical_report_20260606.md
doc/worklog/v2506/dsmcFoam_plus_mpi_omp_mixed_technical_report_20260610.md
doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_further_optimization_technical_report_20260617.md
```

## 1. 本轮目标

超算工作开始后，目标从“本地 8 核性能筛选”转为三件事：

1. 让 replicated-mesh raw-MPI / mixed MPI+OpenMP 在超算环境下稳定运行；
2. 修复有 replicatedMesh 参与时的写出正确性，使超算完整算例能得到可后处理结果；
3. 修复正确性补丁引入的性能退化，恢复接近 2026-06-17 本地最佳性能的计算口径。

本轮特别需要分清两条线：

- **正确性线**：不崩溃、不死锁、写出的 cloud / field 能用于后处理；
- **性能线**：no-write compute 口径下，move/build/collision/DLB 不因保守修复长期退化。

## 2. 超算现场问题与定位

### 2.1 ParMETIS/METIS 环境问题

现场表现：

- `mpi64` 在初始化阶段触发 `Floating point exception`，调用栈进入
  `METIS_PartGraphKway()` / `libmetis__SetupKWayBalMultipliers()` /
  `dsmcReplicatedMesh::computeCellOwnerScotch()`；
- 最小验证程序读取 `constant/polyMesh/owner` / `neighbour` 时先暴露过
  cell 数推导错误，自动推导后又出现 rank 被 `SIGKILL`；
- 超算上同一网格、同一代码，本地正常而超算异常，问题集中到 ParMETIS/METIS
  编译、ABI 和运行时动态库绑定。

处理：

- 写最小 METIS/MPI 验证程序，固定从 OpenFOAM 编译环境抓取编译器和 MPI wrapper，
  只允许用户修改 `PARMETIS_DIR`；
- 重写 ParMETIS 构建脚本，避免直接 `source` 脚本导致超算连接中断，改为子 shell
  中解析 OpenFOAM 环境；
- 确认安装头文件中：

```text
IDXTYPEWIDTH  = 32
REALTYPEWIDTH = 64
```

- 运行时 `ldd` 确认 `libparmetis.so` / `libmetis.so` / `libGKlib.so` 都来自用户
  编译安装目录，并链接到 oneAPI 2022.1 的 MPI/编译器运行时。

结论：

- 原始问题并不只是 `controlDict`，超算上的 ParMETIS/METIS 构建和运行时绑定是
  第一层不确定因素；
- 后续 `tpwgts` 报错说明 ABI 修正后进入了真正的 ParMETIS 参数检查阶段，
  问题从“库不可用/不稳定”转为“传入参数需要严格满足 ParMETIS 约束”。

### 2.2 ParMETIS `tpwgts` sum 错误

现场表现：

```text
PARMETIS ERROR: The sum of tpwgts for constraint #0 is not 1.0
```

触发点：

```text
Phase C ParMETIS: AdaptiveRepart (..., 1 constraints, ubvec=1.05)
ParMETIS_V3_AdaptiveRepart(...)
```

分析：

- 当前源码中 `tpwgts` 使用 `real_t`，按 `ncon * nparts` 填充 `1/nparts`；
- 如果头文件、库、调用侧的 `real_t` 宽度不一致，或者编译时使用的 ParMETIS 头
  与运行时动态库不一致，ParMETIS 会把 `tpwgts` 内存按错误宽度解释；
- 用户检查到 `REALTYPEWIDTH=64` 后，必须保证 OpenFOAM 编译时和运行时都使用同一套
  `metis.h` / `parmetis.h` / 动态库。

结果：

- 重新构建 ParMETIS 后，`tpwgts` 这一类初始化期错误消失；
- 超算 `mpi64` 以及 pure MPI replicated mesh 可以进入长时间运行。

### 2.3 MPI64 / OMP64 / mixed 运行崩溃

现场出现过三类崩溃：

1. `OMP64` 在 `dsmcCloud::buildCellOccupancy()` 处 segfault；
2. `MPI64` 在 ParMETIS/METIS 初始化阶段 FPE；
3. `MPI8+OMP8` 在初始 particle distribution 之后或 SAR/DLB 触发附近 segfault/卡死。

定位后分成两类处理：

- ParMETIS/METIS 问题按 2.1 和 2.2 处理；
- mixed/OpenMP 崩溃与 OpenMP move 中 lazy mesh data 的并发初始化、以及 migration 后
  `moveOrderedParcels_` 缓存生命周期有关。

关键判断：

- `MPI4+OMP16`、`MPI8+OMP8`、`MPI16+OMP4` 曾经能跑，说明 mixed 总体路径可用；
- pure MPI 修复后 mixed 反而出问题，说明后续稳定性修复改变了缓存/预热/迁移时序，
  不是简单的 ParMETIS 库问题。

## 3. 正确性修复内容

### 3.1 DLB trigger 决策一致性

`HEAD` 中 `dsmcReplicatedMesh::autoRebalance()` 增加了跨 rank trigger 决策一致性保护：

- 当某一步已经完成全局 trigger 判断后，用 `MPI_Allreduce` 汇总所有 rank 的
  `triggered`；
- 如果部分 rank 触发、部分 rank 未触发，打印 mismatch，并把结果提升为全局 OR；
- 这样避免部分 rank 进入 ParMETIS/rebalance、另一部分 rank 不进入，导致后续 collective
  顺序不一致而卡死。

同时，`tdecps_` 不再只累积本 rank local decision wall，而是通过 `MPI_Allreduce`
取全局最大值：

```text
localDecWall -> MPI_Allreduce(MPI_MAX) -> globalDecWall -> tdecps_
```

这会让计时更接近关键路径，但也增加了一次 collective。该修改偏正确性/诊断，
对性能不利。

### 3.2 ParMETIS profile 点

`HEAD` 中为 `reassignByParMetisAdaptiveRepart()` 增加了可选 profile：

```text
replicatedMeshDLBProfile true/false
```

覆盖：

- particle-count allreduce；
- ParMETIS AdaptiveRepart；
- partition allgatherv。

用途：

- 分清 DLB 时间到底在 `MPI_Allreduce`、ParMETIS 本身，还是 repartition 结果分发；
- 现场排查时确认卡住位置，而不是只看到 `Phase C auto DLB triggered`。

### 3.3 OpenMP move lazy mesh data 预热

mixed/OpenMP 崩溃暴露出 `Cloud::move()` 进入 OpenMP region 后，部分 mesh 访问器仍会
lazy 初始化。多个 OpenMP worker 同时触发这些 demand-driven 数据时，可能出现
segfault 或不可复现卡死。

当前未提交性能回修补丁中引入：

```cpp
void Cloud<ParticleType>::prepareOpenMPMoveMeshData() const;
```

预热数据缩减为 OpenMP tracking 必要集合：

```text
polyMesh_.tetBasePtIs()
polyMesh_.cells()
polyMesh_.cellVolumes()
polyMesh_.cellCentres()
boundaryMesh.patchID()
boundaryMesh[patchI].faceCells()
cellHasWallFaces()
```

同时：

- `Cloud` 新增 `openmpMoveMeshDataReady_`；
- constructor 初始化为 `false`；
- `autoMap()` 后重置为 `false`；
- `dsmcCloud::evolve_moveAndCollide()` 在 move timer 开始前预热，避免把一次性预热成本
  计入 `move only`。

这一修改兼顾正确性和性能：

- 正确性：避免 OpenMP worker 并发构造 lazy mesh data；
- 性能：预热只做一次，且不污染 move profile。

### 3.4 replicated-mesh write 修复

问题：

有 replicatedMesh 参与时，所有 MPI rank 持有完整网格但只拥有部分 parcel。
如果直接走普通 OpenFOAM write，容易出现：

- 多 rank 写同一套 root cloud 文件；
- cloud 字段和 parcel 数不匹配；
- `foamToEnsight` / 后处理流场不对；
- write 后 particle ownership 状态被破坏。

修复思路：

1. 输出时先 `gatherParcelsToRank0()`；
2. rank0 暂时把 `dsmc.writeOpt()` 设为 `NO_WRITE`，先让 `runTime.write()` 写非 cloud
   对象；
3. rank0 调用 `writeGatheredCloudOnRank0()` 写 gathered cloud 字段；
4. rank0 写 `cellOwner` 场；
5. 输出后调用 `migrateParticlesByCellOwner()` 和 `updateParticleCounts()`，恢复按
   `cellOwner_` 分布的计算状态。

代码路径：

```text
applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C
```

已见运行标记：

```text
Replicated mesh: wrote gathered cloud on rank0 at output time ...
Written cellOwner field at output time ...
Replicated mesh output timing:
    gather parcels [s]
    rank0 write [s]
    migrate-back [s]
```

超算 `pal-phd-3.3.1-react` 多组输出日志已经出现这些标记，例如：

```text
results/palphd3.3.1react-m8/mpi8omp8_3629855.out
results/palphd3.3.1react-m8/mpi4omp16_3629839.out
results/palphd3.3.1react-m8/mpi16omp4_3629840.out
results/palphd3.3.1react-m8/mpi32omp2_3629843.out
results/palphd3.3.1react-m8/mpi64_3629848.out
```

其中最终输出处的 rank0 write 成本约为 `22-24 s`，migrate-back 约为 `1.1-1.6 s`。
这说明 write 正确性已经可跑通，但 write 成本不可忽略，性能对比仍应区分
write-on 与 no-write compute。

## 4. 正确性修复引入的性能退化

`HEAD` 提交信息已经明确：

```text
run correctly in Beijingchaosuan. But make performance worse.
```

性能退化主要来自两处保守化：

### 4.1 migration 后禁用 `moveOrderedParcels_` 复用

为避免 dangling pointer / 不完整缓存导致 mixed 崩溃，`HEAD` 在多个 migration 路径
中把原先的：

```cpp
cloud_.setMoveOrderedParcels(kept);
```

改成：

```cpp
cloud_.clearMoveOrderedParcels();
```

并强制 migration partition 阶段不走 ordered traversal：

```cpp
const bool useOrderedTraversal = false;
```

效果：

- 正确性更稳；
- 但每次 migration 后，后续 `buildCellOccupancy()` 失去 move-ordered 快路径，
  容易退回链表扫描或 fallback gather；
- 对 OMP/mixed 的 move/build profile 影响明显。

### 4.2 OpenMP move 预热最初计入 move profile

最初的 OpenMP lazy mesh data 预热在 `Cloud::move()` 内部、move timer 之后触发。
这会把一次性 mesh data 构造成本计入 `move only`，造成 OMP8 / mixed 的 move
时间异常偏大。

后续补丁把预热提前到 `dsmcCloud::evolve_moveAndCollide()` 的 move timer 前，并把预热
集合缩小到 tracking 必需数据。

## 5. 性能回修补丁

当前未提交性能回修补丁做了两类修复。

### 5.1 恢复 migration 后的 ordered cache 复用，但加完整性 guard

迁移路径改为：

```cpp
const bool useOrderedTraversal =
    cloud_.hasMoveOrderedParcels()
 && cloud_.moveOrderedParcels().size() == cloud_.size();
```

迁移结束时不再无条件清空 cache，而是：

```cpp
if (kept.size() == cloud_.size())
{
    cloud_.setMoveOrderedParcels(kept);
}
else
{
    cloud_.clearMoveOrderedParcels();
}
```

异步 migration 中：

- `migrateBegin()` 先保存本 rank kept parcels；
- `migrateFinish()` 反序列化 received parcels 后，如果已有完整 ordered cache，
  则 `appendBatchToMoveOrdered(received)`；
- 否则清空 cache，避免使用不完整指针列表。

`dsmcCloud` 中所有写入/替换 moveOrdered cache 的接口同步失效 occupancy cache：

```text
occupancyOrderedParcelsValid_ = false
cellOccupancyMaterialized_   = false
```

意义：

- 不再为正确性长期牺牲 `moveOrderedParcels_ -> buildCellOccupancy` 快路径；
- guard 保证 cache size 与 cloud size 一致才复用，避免旧问题复发。

### 5.2 OpenMP move 预热移出 move profile

新增：

```cpp
Cloud<ParticleType>::prepareOpenMPMoveMeshData()
```

并在 `dsmcCloud::evolve_moveAndCollide()` 中：

```cpp
if (openmpEnabled_ && openmpMoveEnabled_ && ompNumThreads_ > 1)
{
    prepareOpenMPMoveMeshData();
}
```

这段位于 move timer 前。`Cloud::move()` 内仍保留一次 fallback 调用，保证如果有其他
入口直接调用 move，也不会在 OpenMP region 内首次构造 lazy mesh data。

## 6. 验证与性能数据

### 6.1 编译验证

最后一次回修编译通过：

```text
doc/worklog/v2506/detail_mix/build_perf_cache_repair_20260621_reduced_prewarm.log
```

构建命令口径：

```bash
source doc/scripts/env.sh
wmake src/lagrangian/basic
wmake src/lagrangian/dsmc
wmake applications/solvers/discreteMethods/dsmc/dsmcFoam+
```

### 6.2 2026-06-17 本地最佳基线

来源：

```text
doc/worklog/v2506/detail_mix/zb_full_series_current_mpi8cfg_repeat3_20260617
doc/worklog/v2506/detail_mix/ourmesh_full_series_current_mpi8cfg_repeat3_20260618
```

`zb-cylinder-react` 300-step repeat3 均值：

| mode | real mean | full evolve mean | move mean | buildOcc mean | collision mean |
|---|---:|---:|---:|---:|---:|
| OMP8 | 55.41 | 50.85 | 37.06 | 3.92 | 8.68 |
| MPI2xOMP4 | 59.26 | 53.35 | 39.78 | 4.79 | 8.54 |
| MPI4xOMP2 | 60.48 | 54.02 | 40.46 | 3.83 | 9.61 |
| MPI8 | 69.32 | 60.67 | 46.39 | 4.98 | 10.75 |
| MPI8origin | 107.68 | 107.65 | 79.64 | 9.87 | 23.71 |

`ourmesh` 500-step repeat3 均值：

| mode | real mean | full evolve mean | move mean | buildOcc mean | collision mean |
|---|---:|---:|---:|---:|---:|
| OMP8 | 63.23 | 57.55 | 47.41 | 6.21 | 2.80 |
| MPI2xOMP4 | 67.43 | 61.79 | 50.32 | 7.45 | 3.23 |
| MPI4xOMP2 | 66.09 | 59.50 | 48.91 | 6.56 | 3.11 |
| MPI8 | 84.25 | 75.28 | 60.46 | 10.43 | 4.85 |
| MPI8origin | 126.63 | 126.50 | 109.11 | 15.21 | 11.49 |

该组数据仍是目前本地 no-write compute 的主要性能基线。

### 6.3 `HEAD` 正确性修复后的回归状态

2026-06-20 no-write 性能复测中：

```text
doc/worklog/v2506/detail_mix/zb_full_series_nowrite_perf_recheck_20260620
```

结果：

| mode | exit | real | 状态 |
|---|---:|---:|---|
| OMP8 | 139 | 10.27 | 未跑完，日志无完整 profile |
| MPI2xOMP4 | 255 | 18.56 | segfault / killed |
| MPI4xOMP2 | 255 | 20.54 | segfault / killed |
| MPI8 | 0 | 68.86 | 正常 |
| MPI8origin | 0 | 112.45 | 正常 |

此时 pure MPI8 已能正常完成，mixed/OMP 出现新的稳定性问题。

MPI8 分项：

| metric | value |
|---|---:|
| real | 68.86 |
| full evolve | 60.07 |
| move | 47.80 |
| buildCellOccupancy | 4.48 |
| collision | 10.94 |
| migration wall max | 3.62 |
| DLB wall max | 5.97 |
| rank wall max/min | 1.30 |

结论：

- pure MPI8 相比 2026-06-17 均值没有明显坏化，甚至接近当时均值；
- 主要问题转为 OMP/mixed 稳定性与后续 cache/prewarm 修复。

### 6.4 性能回修后的单轮五组测试

来源：

```text
doc/worklog/v2506/detail_mix/zb_perf_cache_repair_single_20260621
```

五组均 `exit 0`：

| mode | real | full evolve | move | buildOcc | collision | migration | DLB wall | rank wall max/min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 62.38 | 60.11 | 44.31 | 4.15 | 10.23 | - | - | - |
| MPI2xOMP4 | 65.45 | 61.66 | 45.44 | 3.48 | 8.64 | 2.54 | 2.09 | 1.07 |
| MPI4xOMP2 | 66.12 | 60.32 | 46.29 | 3.03 | 11.44 | 3.43 | 4.24 | 1.10 |
| MPI8 | 83.64 | 74.21 | 55.89 | 5.80 | 12.66 | 4.24 | 6.82 | 1.35 |
| MPI8origin | 116.58 | 124.04 | 92.27 | 10.97 | 26.96 | - | - | - |

这轮证明五组能跑完，但 OMP8 的 move 明显偏大，说明 OpenMP 预热仍在污染 move
profile。

### 6.5 预热移出 move profile 后的直接复测

来源：

```text
doc/worklog/v2506/detail_mix/zb_reduced_prewarm_omp8_20260621
doc/worklog/v2506/detail_mix/zb_reduced_prewarm_mpi8_20260621
```

| mode | real | full evolve | move | buildOcc | collision | migration | DLB wall | rank wall max/min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 56.30 | 54.90 | 40.24 | 3.72 | 9.65 | - | - | - |
| MPI8 | 78.23 | 72.05 | 48.57 | 5.14 | 16.40 | 3.91 | 6.79 | 1.52 |

对比 2026-06-17 最佳单次：

| mode | historical best real | current real | 主要差异 |
|---|---:|---:|---|
| OMP8 | 55.07 | 56.30 | 基本修回；move 仍高约 3.37 s，collision 高约 0.95 s |
| MPI8 | 66.82 | 78.23 | 未修回；move 接近，但 collision 与 DLB/rank imbalance 明显变差 |

判断：

- OMP8 的大头回归已修复；
- MPI8 退化不来自 OpenMP 预热，因为 MPI8 `useOpenMP false`；
- 当前 MPI8 剩余问题集中在 collision phase、DLB check/trigger overhead 和 rank wall
  imbalance。当前 `collision phase = 16.40 s`，而 2026-06-17 最佳单次为
  `10.62 s`。

### 6.6 `replicatedMeshMigrateInterval=10` 的正确性复核

背景：

- write 修复后，mixed 超算结果仍出现局部流场异常；
- `MPI16+OMP4` 在 1/4 圆柱驻点附近与 `OMP64` 对照相比，`p`、`rhoN`、
  `dsmcNMean` 大体接近，但 `wallHeatFlux` 明显偏高；
- 用户随后将 `replicatedMeshMigrateInterval` 从 `10` 改为 `1`，同一类流场和
  壁面结果恢复正常。

典型点值对比：

| quantity | MPI16xOMP4, migrate=10 | OMP64 | difference |
|---|---:|---:|---:|
| `p_mixture` | 307.50 | 301.84 | +1.9% |
| `rhoN_mixture` | 1.866e22 | 1.902e22 | -1.9% |
| `dsmcNMean_mixture` | 125.54 | 132.32 | -5.1% |
| `wallHeatFlux_mixture` | 92501.56 | 69468.62 | +33.2% |

最初判断是 boundary flux 归属问题：粒子跨到非本 rank owner cell 后，在下一次
migration 前仍由旧 rank 执行 move；如果这期间撞到壁面，`dsmcPatchBoundary` 会把
`qBF/fDBF/rhoNBF` 记在旧 rank 的 `boundaryFluxMeasurements()` 中，而
`dsmcVolFields` 在 replicated mesh 输出时会对 boundary 累积量做全局 `sumReduce`。
因此非 owner rank 的壁面事件也会进入最终 `wallHeatFlux`。

进一步分析后，问题被提升为更底层的 owner 一致性问题：

1. `dsmcCloud::evolve_moveAndCollide()` 中 `Cloud<dsmcParcel>::move()` 先执行；
2. replicated mesh 的 `migrateParticlesByCellOwner()` 只在
   `stepCounter % replicatedMeshMigrateInterval == 0` 时执行；
3. move 和可选 migration 之后才 `buildCellOccupancy()`；
4. `noTimeCounter` 默认只遍历 `occupancyOwnedCollisionCells()`，即只对本 rank
   owner cell 做 collision。

这意味着 `replicatedMeshMigrateInterval=10` 时，跨 owner 粒子最多会在旧 rank 上
滞留 9 步：

```text
particle moves from rank A owned cell to rank B owned cell
particle still physically resides on rank A until the next migration

rank A:
  buildCellOccupancy can see the parcel in a B-owned cell
  but collision loops only over A-owned cells
  -> the parcel does not collide in that B-owned cell

rank B:
  owns the B-owned cell
  but does not have the parcel yet
  -> the parcel is also missing from rank B collision
```

因此 `migrationInterval=10` 不是单纯的 write 或 wall flux 归并问题，而是会破坏
当前 raw-MPI replicated mesh 的核心不变量：

```text
before collision/reaction/sampling, every parcel must reside on the owner rank
of its current cell
```

这也解释了用户观察到的通信边界不连续：在 processor-owner 边界附近，跨 owner 粒子
在若干步内没有参与正确 cell 的 collision/reaction/sampling，局部分布函数自然会出现
断裂或偏差。

评估过的修复方向：

| direction | correctness | implementation risk | conclusion |
|---|---|---|---|
| `migrationInterval=1` | 高 | 低 | 当前正确性基线 |
| 只修 wall flux face-owner 归属 | 只修壁面事件 | 中 | 不足以修 bulk flow |
| 非 owner active cell 也参与 collision | 不严格 | 中 | 会把同一物理 cell 的粒子拆成多份局部碰撞 |
| ghost/halo parcel collision | 理论可行 | 高 | 接近重新设计跨分区 DSMC collision |
| 每步迁移但优化 migration 成本 | 高 | 中 | 后续推荐方向 |

结论：

- `replicatedMeshMigrateInterval=10` 与当前 owned-cell collision 架构不自洽；
- `replicatedMeshMigrateInterval=1` 不是单纯 workaround，而是当前实现保证物理一致性的
  必要条件；
- 后续性能工作不应继续尝试让 `migrate=10` 通过局部补丁“看起来正常”，而应在保持
  `migrate=1` 的前提下优化每步 migration 的扫描、打包、通信和计数更新成本。

## 7. 当前结论

### 7.1 正确性

已解决或基本解决：

- 超算 ParMETIS/METIS 构建与运行时绑定问题；
- ParMETIS `tpwgts` / datatype 类错误；
- pure MPI replicated mesh 在超算上初始化和长时间运行问题；
- DLB trigger 跨 rank 决策不一致风险；
- OpenMP move lazy mesh data 并发初始化风险；
- replicatedMesh write 只写 rank0 gathered cloud，并输出 `cellOwner` 场；
- write 后 migrate-back 恢复计算态；
- 明确 `replicatedMeshMigrateInterval=1` 是当前 raw-MPI replicated mesh 的正确性基线。

仍需继续验证：

- full write-on 流场结果是否在 `replicatedMeshMigrateInterval=1` 口径下与 OMP8
  对照完全一致，尤其是 `foamToEnsight` 后的物理场分布；
- 大规模 MPI64 / OMP64 / mixed64 在 write-on 口径下的稳定性和输出后继续计算能力；
- processor-write 路径是否值得作为替代 rank0 gather-write 的大规模输出方案。

### 7.2 性能

已修回：

- OMP8 因 OpenMP mesh data 预热计入 move timer 导致的主要性能退化；
- migration 后完全禁用 ordered cache 造成的部分 build/move 快路径退化。

未完全修回：

- MPI8 当前仍比 2026-06-17 最佳慢，主要不是 move，而是 collision/DLB/rank imbalance；
- write-on 口径下，rank0 gather/write 本身在大粒子数超算 case 中约 `20+ s/次`，
  会显著改变端到端总时间；
- 若从 `replicatedMeshMigrateInterval=10` 改为 `1`，migration 通信和 pack/unpack
  开销会增加，必须重新建立性能基线。

当前建议：

1. 正式性能比较继续使用 no-write compute 口径；
2. write 正确性单独用短步数/指定 output time 验证；
3. 大规模超算生产输出优先考虑减少输出频率，或继续完善 processor-write/reconstructPar
   路径，避免所有 parcel 长期集中到 rank0 写出；
4. replicated mesh 正确性测试统一使用 `replicatedMeshMigrateInterval=1`；
5. MPI8 下一步不应再优先做 OpenMP prewarm，而应回到每步 migration 成本、
   collision rank imbalance、DLB check cadence、以及 collision 权重 proxy 的针对性分析。

## 8. 后续待办

1. `zb-cylinder-react`：
   - 对当前未提交性能回修补丁做 repeat3；
   - 分别记录 OMP8、MPI2xOMP4、MPI4xOMP2、MPI8、MPI8origin；
   - 与 2026-06-17 repeat3 均值比较。
2. MPI8 collision：
   - 开启 detail profile，只看 collision localLoop / reduce / sigmaBC；
   - 记录 rank wall max/min 与 particle max/min；
   - 判断是否需要把 DLB weight 从纯粒子数 proxy 扩展到 active collision-cell 或候选数
     的低开销近似。
3. 每步 migration 成本：
   - 在 `replicatedMeshMigrateInterval=1` 下重建 OMP8、MPI2xOMP4、MPI4xOMP2、MPI8
     的 no-write 性能基线；
   - profile `migrateBegin/migrateFinish`、pack、size exchange、receive/apply、
     `updateParticleCounts()`；
   - 优先优化“每步只迁移跨 owner 粒子”的扫描和通信开销，而不是恢复
     `replicatedMeshMigrateInterval=10`。
4. DLB cadence：
   - 重新验证 `replicatedMeshSARSteps=50`、`MinGapSteps=50`、`allgather` 的当前最优性；
   - 若当前 `checks=300` 仍导致 check wall 偏大，继续把 expensive collective 与
     per-step entry 分开优化。
5. write：
   - 对比 rank0 gather-write 与 processor-write；
   - 检查 `foamToEnsight` 后的场和 OMP8 对照；
   - 如果 rank0 write 成本成为超算总时间主要部分，优先推进 processor-write。

### 8.1 processor-write 后续实现状态

已在 2026-06-21 补齐 raw-MPI replicated-mesh 真正 processor 写出闭环：

```text
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/processor_write_full_20260621.md
```

当前验证结论：

- `replicatedMeshWriteMode processor` 下，运行时不再由 rank0 串行写 root 全局场、
  root cloud 或 root `cellOwner`；
- 每个 MPI rank 写自己的 `processorN/` mesh/addressing、volFields 和
  lagrangian cloud；
- processor mesh/addressing 按 owner 分区版本复用，避免每个 output time 重写；
- 标准 `reconstructPar` 可重构 root volFields 和 lagrangian cloud；
- 2-rank clean smoke 中，`1e-07` 和 `2e-07` 的 root 粒子数均等于 processor
  粒子数之和；
- DLB owner 变化验证已按 `mpi4omp2` 完成：1000 steps，500-step output，
  250-step forced DLB，4 次 DLB 与 2 次 processor output 均成功；
- 为适配 DLB 后 processor patch 拓扑变化，processor-write 路径现在固定写出
  `nProcs-1` 个 processor patch；非邻居 rank 使用 0-face patch，保证
  `reconstructPar` 连续读取多个 output time 时 patch list 稳定；
- `mpi4omp2` DLB 验证中，`3.320025e-05` 与 `6.64005e-05` 两个 output time
  均 `reconstructPar` 成功，root 粒子数等于各 processor 粒子数之和；
- 编译、运行、重构验证统一使用 `doc/scripts/env.sh` 和
  `doc/scripts/build-dsmcFoam.sh`，避免旧用户 OpenFOAM `reconstructPar` /
  `libreconstruct.so` 与当前库混用。

最新验证路径：

```text
doc/worklog/v2506/detail_mix/build_processor_write_true_parallel_20260621.log
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/log.true_parallel_final_smoke_run
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/log.true_parallel_final_smoke_reconstruct
doc/worklog/v2506/detail_mix/build_processor_write_fixed_proc_patches_20260621.log
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/log.zb_mpi4omp2_dlb250_write500_1000step_fixedpatches_run
doc/worklog/v2506/detail_mpi/processor_write_full_20260621/log.zb_mpi4omp2_dlb250_write500_1000step_fixedpatches_reconstruct
```

生产使用注意：

- processor-write 的目标是降低大网格、大粒子数下 root 串行 write 成本；
- 因为运行阶段不再生成 root time 场，`foamToEnsight` 或 root case 后处理必须放在
  `reconstructPar` 之后；
- 该路径和 no-write compute 性能基线应分开评估。

## 9. 证据路径

构建：

```text
doc/worklog/v2506/detail_mix/build_perf_cache_repair_20260621_reduced_prewarm.log
```

历史本地最佳：

```text
doc/worklog/v2506/detail_mix/zb_full_series_current_mpi8cfg_repeat3_20260617/aggregate.csv
doc/worklog/v2506/detail_mix/ourmesh_full_series_current_mpi8cfg_repeat3_20260618/aggregate.csv
```

正确性修复后回归状态：

```text
doc/worklog/v2506/detail_mix/zb_full_series_nowrite_perf_recheck_20260620/results.csv
doc/worklog/v2506/detail_mix/zb_full_series_nowrite_perf_recheck_20260620/aggregate.csv
```

性能回修测试：

```text
doc/worklog/v2506/detail_mix/zb_perf_cache_repair_single_20260621
doc/worklog/v2506/detail_mix/zb_reduced_prewarm_omp8_20260621
doc/worklog/v2506/detail_mix/zb_reduced_prewarm_mpi8_20260621
```

`replicatedMeshMigrateInterval` 正确性复核相关源码：

```text
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C
src/lagrangian/dsmc/boundaries/basic/dsmcPatchBoundary/dsmcPatchBoundary.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
```

超算 write-on 输出日志：

```text
results/palphd3.3.1react-m8/mpi4omp16_3629839.out
results/palphd3.3.1react-m8/mpi8omp8_3629855.out
results/palphd3.3.1react-m8/mpi16omp4_3629840.out
results/palphd3.3.1react-m8/mpi32omp2_3629843.out
results/palphd3.3.1react-m8/mpi64_3629848.out
```

ParMETIS 最小验证和重构脚本：

```text
doc/worklog/v2506/detail_mpi/metis_minimal_validation_20260618
doc/worklog/v2506/detail_mpi/parmetis_rebuild_20260619/makeParMETIS.sh
```

## 10. 2026-06-24 续补：超算 64 核系列、FastRNG 与 post/output 优化

### 10.1 本轮新增范围

本轮在前述正确性修复基础上继续补齐三类工作：

- 超算 `palphd3.3.1react-m8-fixoutput` 完整 8750 步 64 核系列性能复核；
- `collisionFastRng true` 在 OpenMP、mixed MPI+OpenMP、pure MPI replicated mesh collision 路径中的实际生效修正；
- `dsmcVolFields` 后处理计算和 replicated processor-write 输出字段数量优化。

该节只记录这轮最新状态，不覆盖前文关于 ParMETIS、processor-write 正确性和
`replicatedMeshMigrateInterval=1` 的结论。

### 10.2 最新超算性能结果

结果目录：

```text
results/palphd3.3.1react-m8-fixoutput
```

算例口径：

- Palharini PhD 3.3.1 reactive cylinder case；
- 总核数 64；
- `replicatedMeshMigrateInterval=1`；
- 8750 iterations，对应 `endTime 0.0035`、`deltaT 4e-07`；
- mixed/pure MPI replicated mesh 走 processor-write 输出；
- 表中 `ClockTime` 取最终 `Stage 1.0` 行，profile wall 取最终
  `DSMC solver profile summary`。

| case | log | ClockTime [s] | move+collide [s] | move [s] | buildCellOccupancy [s] | collision [s] | migration [s] | DLB 次数 | 结论 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| OMP64 | `omp64_3637765.out` | 3140 | 1786.35 | 892.42 | 421.54 | 430.05 | - | - | 纯 OMP 在该超算 case 不再是最优，collision/build/move 总体偏大。 |
| MPI2xOMP32 | `mpi2omp32_3637766.out` | 3096 | 2124.54 | 810.85 | 190.96 | 414.33 | 828.14 | 43 | migration 和 DLB 过重，性能最差之一。 |
| MPI4xOMP16 | `mpi4omp16_3637767.out` | 2390 | 1731.05 | 814.76 | 133.05 | 382.48 | 502.73 | 30 | 比 MPI2xOMP32 明显好，但 migration 仍重。 |
| MPI8xOMP8 | `mpi8omp8_3637768.out` | 1923 | 1420.56 | 755.70 | 111.97 | 336.05 | 272.75 | 14 | mixed 中间档，collision 最低，但不是总时间最优。 |
| MPI16xOMP4 | `mpi16omp4_3637769.out` | 1759 | 1417.55 | 760.07 | 114.38 | 381.68 | 249.56 | 18 | 接近最优，migration 低于 MPI8xOMP8，但 collision 偏高。 |
| MPI32xOMP2 | `mpi32omp2_3637770.out` | 1680 | 1370.49 | 756.26 | 126.55 | 399.90 | 196.28 | 11 | 当前这组完整 8750 步超算系列最优。 |
| MPI64 | `mpi64_3637771.out` | 2291 | 1882.19 | 930.64 | 417.37 | 401.50 | 386.18 | 42 | pure MPI64 可正常跑完，但性能明显不如 MPI32xOMP2。 |

直接结论：

- 当前 64 核生产配置优先级应是 `MPI32xOMP2`，其次 `MPI16xOMP4`；
- `MPI64` 不是最佳，慢点主要不是单一 collision，而是 `move/buildCellOccupancy/migration/DLB`
  的组合开销；
- `OMP64` 在该超算 case 也不是最佳，说明大核数 OpenMP collision 的线程收益已经被
  build/move、线程调度、内存访问和后处理阶段限制；
- `MPI8xOMP8` 的 collision profile 最低，但总时间输给 `MPI32xOMP2`，说明不能只看
  collision，需要同时看 migration 和 DLB 触发代价。

### 10.3 FastRNG 修正

这轮修正后，`collisionFastRng true` 不再只对一部分 OpenMP collision 生效，而是在：

- OpenMP collision 线程循环中为每个线程建立 `FastRng`；
- mixed MPI+OpenMP 中用 `rank + thread + timeIndex` 生成线程本地种子；
- pure MPI replicated mesh 单线程 collision 路径中也建立 `FastRng`；
- collision/reaction 调用链通过 `cloud_.setCollisionRngContext(...)` 使用同一个
  fast RNG context，避免 reaction 内部仍退回原始 `rndGen_` 热路径。

关键源码位置：

```text
src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C
  OpenMP path: threadFastRng + setCollisionRngContext
  replicated single-thread path: FastRng + setCollisionRngContext
```

`mpi32omp2beforefastrng_3634454.out` 与 `mpi32omp2_3637770.out` 的单次对比：

| case | ClockTime [s] | move+collide [s] | move [s] | buildCellOccupancy [s] | collision [s] | migration [s] | DLB 次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| before FastRNG scope fix | 1569 | 1251.37 | 716.93 | 123.53 | 344.43 | 220.61 | 13 |
| after FastRNG scope fix | 1680 | 1370.49 | 756.26 | 126.55 | 399.90 | 196.28 | 11 |

这个差异不能简单解释为“DLB 次数增加”：修正后 DLB 次数反而从 13 降到 11，
migration wall 也从 220.61 s 降到 196.28 s。性能增加主要体现在 move 和 collision
本身，说明 FastRNG scope 修正改变了随机序列和微观演化轨迹后，粒子分布、候选对、
反应路径和局部负载状态都可能发生真实变化。由于这是随机 DSMC 程序，修正 RNG
热路径后不应期待 bitwise 复现旧日志。

当前判断：

- `collisionFastRng false` 会退回原始 RNG 热路径，通常更慢，不建议作为性能配置；
- FastRNG 修正后的性能需要用新基线评估，不能再拿修正前日志作为“同一物理轨迹”的严格对照；
- 若要区分纯波动和 RNG 轨迹差异，需要固定同一版本做 repeat3，而不是跨修正版本单次比较。

### 10.4 post/output 开销定位与优化

此前确认 `post fields/output` 的主要开销不是文件写入本身，而是
`dsmcVolFields::calculateField()` 中宏观场计算。为降低大算例 output 成本，本轮做了三类优化：

1. `dsmcVolFields` 共享采样 cache：
   - 在单次 field calculation 中先构建 per-cell/per-species 的共享累积量；
   - 混合场 `rhoN/rhoM/p/Ttra/Trot/Tvib/U` 等从共享 cache 派生，减少多场重复遍历 parcel；
   - 预计算 type mass、rotational zeta、vibrational quantum 等常量；
   - 对 active cell 做局部 reset，避免每次清空全域大数组。

2. processor cloud 可选字段懒分配：
   - `radialWeight`、`ERot`、`ELevel`、`stuckToWall/wallTemperature/wallVectors`、
     tracked 相关字段先在本 rank 检测，再用 `MPI_Allreduce` 得到全局是否需要；
   - 只有全局确实存在对应信息时才分配和写出这些字段；
   - 基础字段 `position/U/vibLevel/typeId/newParcel/classification` 仍保持写出。

3. 输出字段白名单：
   - `dsmcVolFieldsProperties` 支持 `writeFields (...)`；
   - 不配置时保持旧行为，仍写全部字段；
   - 可配置只写关心的宏观场，例如 `writeFields (p Ttra U wallHeatFlux);`；
   - 该白名单只控制 volFields 写出数量，不改变求解器物理推进。

新增控制项：

```text
replicatedMeshProcessorWriteCloud true;   // 默认 true，写 processor lagrangian cloud
replicatedMeshProcessorWriteCloud false;  // 只写 processor mesh/volFields，跳过 cloud
```

`replicatedMeshProcessorWriteCloud false` 的用途是做纯宏观场输出或调试后处理开销。
如果后续需要 `reconstructPar` 恢复 lagrangian cloud，则不能关闭。

本地验证结果：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi4omp2/log.codex_mpi4omp2_10step_outputfilter_cloudskip_np4_omp2_foreground_20260624_155259
```

- 临时配置 `writeFields (p Ttra U wallHeatFlux)`；
- 临时配置 `replicatedMeshProcessorWriteCloud false`；
- 输出日志显示 rank0 写出 `scalar=19, vector=6`，并跳过 `processor0 lagrangian cloud`；
- `processor fields/cloud write = 0.399130284 s`；
- `post field calculate = 0.439340028 s`；
- `post field write = 0.001786365 s`；
- 10 steps 正常 `End main`。

no-write compute 口径回归：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi4omp2/log.codex_mpi4omp2_300step_outputfilter_nowrite_np4_omp2_foreground_20260624_154408
```

- 严格运行模式：`4 MPI x 2 OpenMP`，`mpirun -np 4`，`OMP_NUM_THREADS=2`；
- 300 steps 正常 `End main`；
- `real 68.23 s`；
- `post field calculate = 0`、`post field write = 0`，证明 no-write 口径没有引入后处理开销；
- `move+collide wall = 67.85252556 s`；
- `move only = 46.46549126 s`；
- `buildCellOccupancy = 3.738757474 s`；
- `collision phase = 10.53373378 s`；
- `migration wall time = 20.01529794 s`；
- DLB rebalances = 2。

### 10.5 当前配置建议

超算 full-run 生产：

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

64 核建议优先尝试：

```text
MPI32 x OMP2
```

若关注壁面热流和通信边界正确性，继续保持：

```text
replicatedMeshMigrateInterval 1;
```

不建议为了性能回到 `migrate=10`。之前已经确认 `migrate=10` 可能导致通信边界处粒子
重定位滞后，进而影响局部流场连续性和壁面事件统计。生产上应先接受 `migrate=1`
作为正确性基线，再优化每步 migration 的实现成本。

输出配置：

- 需要完整 lagrangian 重构时：保持 `replicatedMeshProcessorWriteCloud true` 或不写该项；
- 只需要宏观场时：可以临时设为 `false`，同时用 `writeFields (...)` 限制 volFields；
- 任何性能对比都应区分 no-write compute、processor-write、reconstructPar 和
  foamToEnsight 四个阶段。

### 10.6 仍需跟进的问题

1. 超算 `MPI32xOMP2` 需要 repeat3，确认 1680 s 是否稳定优于 `MPI16xOMP4`。
2. FastRNG 修正后应重新建立 mixed 系列基线，旧的 before-fastrng 单次日志只能作历史参考。
3. post 计算已经减少重复遍历，但宏观场计算本身不可避免；后续收益主要来自字段白名单、
   sampleInterval、减少输出频率和只在需要时写 cloud。
4. `MPI64` 的 buildCellOccupancy 和 migration 仍明显偏高，后续若继续优化 pure MPI，
   应优先看 occupancy 数据结构、迁移扫描和每步通信，而不是继续调 DLB 权重。
5. 壁面热流偏高/不平滑问题和并行模式相关时，必须保持 `migrate=1`、统一 processor
   reconstruct 流程，并单独检查采样面积、面 owner 事件归属和采样步数，不应混入性能口径。

### 10.7 新增证据路径

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

相关源码：

```text
src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.H
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H
```
