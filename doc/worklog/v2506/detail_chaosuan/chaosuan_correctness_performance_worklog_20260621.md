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

## 7. 当前结论

### 7.1 正确性

已解决或基本解决：

- 超算 ParMETIS/METIS 构建与运行时绑定问题；
- ParMETIS `tpwgts` / datatype 类错误；
- pure MPI replicated mesh 在超算上初始化和长时间运行问题；
- DLB trigger 跨 rank 决策不一致风险；
- OpenMP move lazy mesh data 并发初始化风险；
- replicatedMesh write 只写 rank0 gathered cloud，并输出 `cellOwner` 场；
- write 后 migrate-back 恢复计算态。

仍需继续验证：

- full write-on 流场结果是否与 OMP8 对照完全一致，尤其是 `foamToEnsight`
  后的物理场分布；
- 大规模 MPI64 / OMP64 / mixed64 在 write-on 口径下的稳定性和输出后继续计算能力；
- processor-write 路径是否值得作为替代 rank0 gather-write 的大规模输出方案。

### 7.2 性能

已修回：

- OMP8 因 OpenMP mesh data 预热计入 move timer 导致的主要性能退化；
- migration 后完全禁用 ordered cache 造成的部分 build/move 快路径退化。

未完全修回：

- MPI8 当前仍比 2026-06-17 最佳慢，主要不是 move，而是 collision/DLB/rank imbalance；
- write-on 口径下，rank0 gather/write 本身在大粒子数超算 case 中约 `20+ s/次`，
  会显著改变端到端总时间。

当前建议：

1. 正式性能比较继续使用 no-write compute 口径；
2. write 正确性单独用短步数/指定 output time 验证；
3. 大规模超算生产输出优先考虑减少输出频率，或继续完善 processor-write/reconstructPar
   路径，避免所有 parcel 长期集中到 rank0 写出；
4. MPI8 下一步不应再优先做 OpenMP prewarm，而应回到 collision rank imbalance、
   DLB check cadence、以及 collision 权重 proxy 的针对性分析。

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
3. DLB cadence：
   - 重新验证 `replicatedMeshSARSteps=50`、`MinGapSteps=50`、`allgather` 的当前最优性；
   - 若当前 `checks=300` 仍导致 check wall 偏大，继续把 expensive collective 与
     per-step entry 分开优化。
4. write：
   - 对比 rank0 gather-write 与 processor-write；
   - 检查 `foamToEnsight` 后的场和 OMP8 对照；
   - 如果 rank0 write 成本成为超算总时间主要部分，优先推进 processor-write。

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

