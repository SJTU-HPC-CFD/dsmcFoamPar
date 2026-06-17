# dsmcFoam+ MPI replicated mesh 进一步优化技术报告 - OFv1706 hyStrath_dlb

日期：2026-06-17

工作目录：

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

参考实现：

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
```

本文是 `doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_dlb_technical_report_20260608.md`
之后的 follow-up 收尾报告，覆盖 2026-06-11 到 2026-06-17 这轮 MPI replicated-mesh
进一步优化工作。主线包括：

- collision owned-cell 迭代裁剪和单线程 collision kernel 热路径优化；
- MPI replicated collision 的串行基线、OMP 对照、subphase profile 和负载分析；
- `ourmesh` 上 alpha / Dual / adaptive-alpha / forced-vs-auto 的 DLB 行为复核；
- `zb-cylinder-react` 上当前 MPI8 replicated-mesh DLB 最优配置的矩阵筛选、repeat3
  确认和回归复测；
- 已保留、默认关闭、以及明确证伪的后续优化方向。

本文与以下报告互为同阶段技术文档：

```text
doc/worklog/v2506/dsmcFoam_plus_omp_technical_report_20260606.md
doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_dlb_technical_report_20260608.md
doc/worklog/v2506/dsmcFoam_plus_mpi_omp_mixed_technical_report_20260610.md
```

本轮详细分析和单次/多次实验日志主要位于：

```text
doc/worklog/v2506/detail_mix
```

## 1. 结论摘要

当前这轮 follow-up 的结论可以压缩为五点：

1. MPI replicated-mesh collision 的第一层问题确实存在于“遍历空间和串行热路径”，
   但这不是最终瓶颈。
   - 通过 owned collision-cell 迭代和串行 `noTimeCounter` kernel 优化，
     `ourmesh` MPI8 500-step `collision phase` 从 `32.80 s` 降到 `25.96 s`，
     `real` 从 `101.45 s` 降到 `92.69 s`。
2. MPI replicated collision 的主导问题随后转移为“同步等待和 phase skew”，
   而不是继续缺少某个局部 kernel 微优化。
   - 在 `ourmesh` 详细 subphase profile 中，MPI8 的最大 local collision loop
     只有 `5.29 s`，但 collision phase 可被记为 `34.28 s`，其中最大
     `reduce` 等待达到 `32.40 s`。
3. OMP8 对 collision 的“超级加速”不能被 pure MPI8 完整复现，根因不是单一的
   `particles max/min` 或总墙钟不均衡，而是：
   - OMP 具有单进程共享内存下对 active collision cells 的步内动态调度；
   - MPI replicated mesh 以 rank-owned cells 为粗粒度静态切分，关键路径由最慢 rank
     决定；
   - DLB 当前更容易平衡 move+collision 总和，而不是 collision 候选负载本身。
4. `zb-cylinder-react` 300-step 上，当前确认的 MPI8 replicated-mesh 最优配置为：

```text
replicatedMeshNoAlltoall false;
replicatedMeshDLBDualConstraint false;
replicatedMeshDLBAlpha 1;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshSARSteps 50;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBParticleGate false;
replicatedMeshDLBAdaptiveAlpha false;
replicatedMeshDLBMinRemainingSteps 50;
replicatedMeshGatherCandidates false;
replicatedMeshOverlapSizeExchange false;
```

   - 2026-06-17 首次 repeat3 均值：`real = 65.88 s`
   - 同日回归复测 repeat3 均值：`real = 66.54 s`
   - 相比此前确认的最佳单次 `67.65 s`，平均改进约 `2.6%`
5. 本轮有两个方向被明确证伪为“不应默认开启”：
   - `DualConstraint true` 在 `ourmesh` 和 `zb` 都没有带来端到端收益；
   - `replicatedMeshOverlapSizeExchange true` 的 `MPI_Ialltoall` overlap
     实验在 `zb` 上反而慢于当前最优。

## 2. 环境、算例和复现实验规则

### 2.1 环境

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
```

### 2.2 主要算例

`ourmesh` 500-step：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh
```

`zb-cylinder-react` 300-step：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8
```

对照路径还用到了：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/omp8
```

### 2.3 运行口径

本轮正式性能分析继续采用 no-write compute 口径：

- `ourmesh`：500 steps；
- `zb-cylinder-react`：300 steps；
- `profileSummary true`；
- 按需要开启 `profileDetail true` 做 subphase / rank 级分析；
- replicated mesh 使用 raw-MPI 路径，不使用 `-parallel` processor mesh。

报告中所有 DLB cadence 解释继续严格区分三件事：

1. `autoRebalance()` 是否每步进入；
2. expensive collective check 是否受 `replicatedMeshSARSteps` 窗口控制；
3. 真正 ParMETIS repartition 是否发生。

## 3. 当前源码保留状态

相对 2026-06-08 报告，本轮新增并保留的关键实现分为四组。

### 3.1 collision owned-cell 与单线程 kernel

保留状态：

- `dsmcCloud` 新增 `occupancyOwnedCollisionCells_`；
- `buildCellOccupancy()` 在 replicated mesh 下显式构建 owned + active collision cells；
- `noTimeCounter` 的 OMP 与串行 collision 路径都改为遍历
  `occupancyOwnedCollisionCells()`；
- 串行 collision kernel 复用了 scratch list，缓存 `constProps` /
  `pairModelAddressing` / `binaryCollision` 等热路径引用，并使用 `randomIndex(n)`。

对应明细：

```text
doc/worklog/v2506/detail_mix/collision_owned_cells_iterator_20260611.md
doc/worklog/v2506/detail_mix/mpi8_serial_collision_kernel_rngindex_20260611
doc/worklog/v2506/detail_mix/mpi8_collision_cached_refs_20260611
```

### 3.2 collision / migration profiling 与默认行为调整

保留状态：

- `noTimeCounter` 保留 collision subphase profile；
- `dsmcCloud::printProfileSummary()` 保留 raw-MPI replicated mesh 的 rank profile 汇总；
- `collisionReduceOnlyOnOutput` / `collisionOutputGlobalReduce` 已接入；
- replicated mesh 默认路径已避免把每步 collision global reduce 继续放在关键路径上；
- 最终 summary 仍保留一次全局 collisions/candidates 汇总。

对应明细：

```text
doc/worklog/v2506/detail_mix/ourmesh_mpi8_collision_subphase_20260612
doc/worklog/v2506/detail_mix/ourmesh_mpi8_profile_opt_20260612
```

### 3.3 DLB 控制和自动 alpha 基础设施

保留状态：

- `replicatedMeshDLBAdaptiveAlpha` 及其上下界/步长/增益控制；
- `replicatedMeshDLBParticleGate` 仍存在，默认值仍是 `true`，但本轮正式最优配置
  在需要稳定对比时会显式关掉；
- `replicatedMeshDLBCheckCollective`、`replicatedMeshSARSteps`、`replicatedMeshDLBMinGapSteps`
  仍是当前 check cadence 的主要控制项。

对应明细：

```text
doc/worklog/v2506/detail_mix/ourmesh_mpi8_alpha08_dualfalse_auto_nogate_adaptive_alpha_20260613
```

### 3.4 replicated-mesh 迁移/DLB follow-up 开关

当前源码新增并保留以下控制项：

```text
replicatedMeshGatherCandidates false
replicatedMeshOverlapSizeExchange false
replicatedMeshDLBMinRemainingSteps 0/自定义
```

语义：

- `replicatedMeshGatherCandidates false`
  - 跳过当前路径实际上未消费的 per-rank candidate `MPI_Allgather`
- `replicatedMeshOverlapSizeExchange false`
  - 默认关闭 `MPI_Ialltoall` size-exchange overlap 实验路径
- `replicatedMeshDLBMinRemainingSteps`
  - near-end guard，避免在剩余步数太少时再触发一次几乎无法回本的 DLB

## 4. collision 方向：已完成的优化和证据链

### 4.1 owned collision-cell 迭代裁剪

问题起点是 replicated mesh collision 仍然按全 `nCells` 或过大的 active-cell
集合遍历，导致总迭代量接近 `nRanks x nCells`。

引入 `occupancyOwnedCollisionCells_` 后，`ourmesh` 500-step 单次快速对比为：

| 模式 | real [s] | move | buildOcc | collision | 结论 |
|---|---:|---:|---:|---:|---|
| OMP8 | 62.31 | 47.75 | 6.58 | 3.49 | 基本持平 |
| MPI2xOMP4 | 69.89 | 52.17 | 7.58 | 16.04 | 温和改善 |
| MPI8 | 97.59 | 68.05 | 10.22 | 24.82 | 显著改善 |

这一步证明：owned-cell 裁剪是必要条件，但还不是充分条件。

### 4.2 单线程 collision kernel 热路径优化

在 `noTimeCounter` 串行 fallback 中，继续做了：

- scratch list 复用；
- 遍历 `occupancyCollisionCells()` / `occupancyOwnedCollisionCells()`；
- parcel/type/charge/model lookup 缓存；
- `randomIndex(n)`；
- 精确 candidate clear 和热点引用缓存。

`ourmesh` MPI8 500-step 三次均值从 clean mean 到 cached-ref mean 的改善为：

| metric | clean MPI8 mean [s] | cached-ref mean [s] | delta |
|---|---:|---:|---:|
| real | 101.446667 | 92.690000 | -8.63% |
| full evolve | 98.629420 | 92.260071 | -6.46% |
| collision | 32.798558 | 25.957368 | -20.86% |

这一步说明 kernel 优化对 MPI8 仍然有效，但改善量级已经不可能解释 OMP8 对 MPI8
近 10 倍的 collision 差距。

### 4.3 单核串行基线：OMP8 和 MPI8 的 collision speedup 并不对称

`ourmesh` 500-step 串行基线：

| mode | collision [s] | speedup vs serial1 |
|---|---:|---:|
| serial1 | 14.636285 | 1.00x |
| OMP8 | 3.003647 | 4.87x |
| MPI8 replicated current mean | 25.957368 | 0.56x |

即：

- OMP8 collision 真正比单核串行快 `4.87x`；
- MPI8 replicated collision 不仅没有接近 `8x`，反而比单核串行慢 `1.77x`。

`zb-cylinder-react` 300-step 串行基线同样支持这一点：

| mode | collision [s] | speedup vs serial |
|---|---:|---:|
| serial1 | 49.2992 | 1.00x |
| OMP8 current | 11.39 | 4.33x |
| MPI8 replicated current | 42.58 | 1.16x |

结论很明确：MPI replicated collision 的剩余差距，已经不能再归因于“缺某个串行 kernel
优化”或者“collision 负载太小”。

## 5. 为什么 OMP8 的 collision 超级加速无法在 MPI8 上完全复现

### 5.1 `ourmesh`：collision phase 很大，但 local collision loop 很小

`ourmesh` MPI8 subphase detail run 的关键数据：

| 指标 | 数值 |
|---|---:|
| max local collision loop | 5.285824 s |
| max collision reduce | 32.396364 s |
| max collision phase | 34.279284 s |
| rank full max/min | 1.00007 |

这说明 collision phase 变大主要是因为：

- 快 rank 更早到达 `noTimeCounter::collide()` 里的 collective/reduce；
- 等待慢 rank 的 move/build/local collision 完成；
- 这些等待时间被记到了 collision bucket。

也就是说，MPI8 的“大 collision”在很大程度上并不是“有用碰撞计算”。

### 5.2 去掉 per-step collision reduce 后，等待只会迁移位置

实验控制 `collisionReduceOnlyOnOutput true` 后：

- `collision phase` 可以从 `25.96 s` 级别塌缩到约 `5.30 s`；
- 但 `real` 不会同步改善，等待转移到下一处 DLB check / migration collective。

因此这不是一个可独立成立的端到端优化。它证明了等待来源，但没有消除等待。

### 5.3 `zb`：DLB 能平衡总时间，却平衡不了 collision 候选负载

`zb` 300-step detail run 中：

| 指标 | 数值 |
|---|---:|
| full wall max/min | 1.000463931 |
| collision wall max/min | 1.935047402 |
| cumulative candidates max/min | 2.944458304 |

同时 inter-DLB 实际负载块显示 move 和 collision 在 rank 间呈反相关：

- move 重的 rank，collision 往往轻；
- collision 重的 rank，move 往往轻；
- 总 `move+collision` 看起来很平，但 collision 单相并不平。

这正是 pure MPI8 无法复制 OMP8 collision 加速的根因：

1. OMP8 是单进程共享 active collision cells，步内可动态调度；
2. MPI8 replicated mesh 是 rank-owned static cells，步内不能把 hot rank 的碰撞
   工作再切给别的 rank；
3. 当前 DLB 主要按粒子数 proxy、move/collision 混合代价近似或时间窗口做决策，
   不足以直接平衡 collision candidates 的极端偏斜。

## 6. `ourmesh`：DLB alpha / Dual / adaptive-alpha 的复核结论

### 6.1 forced8 严格对比下，`DualConstraint true` 没有回本

`ourmesh` 500-step forced8 严格比较：

| Variant | real [s] | full evolve [s] | move [s] | build | collision | 结论 |
|---|---:|---:|---:|---:|---:|---|
| `alpha=1, Dual=false` | 77.85 | 71.70 | 57.15 | 6.06 | 4.09 | 最优 |
| `alpha=0.8, Dual=false` | 77.79 | 75.00 | 59.23 | 8.84 | 4.94 | external 接近，但内部更差 |
| `alpha=0.8, Dual=true` | 85.27 | 82.70 | 65.13 | 10.37 | 4.41 | 明显更慢 |

结论：

- `DualConstraint true` 虽然有时改善 rank-wall 比值，但无法回收额外的 repartition /
  migration / build / move 成本；
- `alpha=1, Dual=false` 仍是更稳健的单约束基线。

### 6.2 adaptive alpha 可以改善 auto no-gate，但没有打赢 forced8 alpha1

`ourmesh` 500-step auto no-gate 比较：

| Variant | real [s] | full evolve [s] | move [s] | build | collision | DLB rebalances |
|---|---:|---:|---:|---:|---:|---:|
| `alpha=0.8, adaptive=true` | 77.78 | 73.76 | 57.59 | 9.10 | 3.87 | 7 |
| `alpha=0.8, adaptive=false` | 79.49 | 74.16 | 58.25 | 8.30 | 3.98 | 7 |
| `alpha=1.0, adaptive=false` | 82.52 | 79.21 | 62.17 | 9.40 | 4.31 | 6 |
| forced8 `alpha=1.0, adaptive=false` | 77.85 | 71.70 | 57.15 | 6.06 | 4.09 | 8 |

解释：

- adaptive alpha 的确可以改善 auto no-gate `alpha=0.8`；
- 但它更像是“让非 forced 路径更接近 forced8 alpha1”，而不是已经超越了后者；
- 从 `full evolve` 看，forced8 `alpha=1, Dual=false` 仍然更干净。

因此 adaptive alpha 目前属于“已实现、可继续研究、但不是当前默认推荐路径”。

## 7. `zb-cylinder-react`：本轮 DLB 调优和当前最优配置

### 7.1 筛选结论

`zb` 上的矩阵和 repeat 结果最终把最优配置收敛到：

```text
replicatedMeshNoAlltoall false;
replicatedMeshDLBDualConstraint false;
replicatedMeshDLBAlpha 1;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshSARSteps 50;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBParticleGate false;
replicatedMeshDLBAdaptiveAlpha false;
replicatedMeshDLBMinRemainingSteps 50;
replicatedMeshGatherCandidates false;
replicatedMeshOverlapSizeExchange false;
```

关键判断：

- `SARSteps=50` 是真正有效的 cadence lever；
- `allgather` 比 `allreduce` 更合适；
- `Dual=false` 比 `Dual=true` 更稳；
- `alpha=1` 比 `alpha=0.8` 更稳；
- `NoAlltoall false` 在这个 8-rank `zb` case 上优于 no-Alltoall 路径；
- step 300 near-end DLB 应被 guard 掉；
- `allProcCandidates_` 当前不参与决策，因此 candidate gather 应默认关闭。

### 7.2 repeat3 确认

2026-06-17 当前最优配置 repeat3：

| run | real [s] | full evolve [s] | sizeX max [s] | wait max [s] | DLB rebalances | 说明 |
|---|---:|---:|---:|---:|---:|---|
| rep1 | 66.18 | 56.94 | 2.595 | 0.077 | 5 | step 300 skipped |
| rep2 | 65.80 | 56.53 | 2.234 | 0.071 | 5 | step 300 skipped |
| rep3 | 65.66 | 54.87 | 2.745 | 0.072 | 5 | step 300 skipped |
| mean | 65.88 | 56.11 | 2.525 | 0.073 | 5 | 当前确认最优 |

相对历史最佳单次 `67.65 s`，平均提升约 `2.62%`。

### 7.3 回归复测

在 `replicatedMeshOverlapSizeExchange` 实验之后，对当前最优再做了 repeat3 回归：

| run | real [s] |
|---|---:|
| rep1 | 67.94 |
| rep2 | 67.25 |
| rep3 | 64.42 |
| mean | 66.54 |

解释：

- 这组均值比 `65.88 s` 略慢，但仍维持同一最优配置；
- migration `sizeX max`、`wait max` 与前一组基本同一量级；
- 差异主要来自 DLB/check 波动以及 move/collision 随机负载波动。

因此当前最优配置仍成立，只是推荐口径应同时记录：

- 稳定确认均值：`65.88 s`
- 同日复测均值：`66.54 s`

## 8. 本轮明确放弃或默认关闭的方向

### 8.1 `DualConstraint true`

本轮在 `ourmesh` 和 `zb` 都没有看到端到端收益：

- collision 有时会下降；
- 但 total / full evolve / move / build / migration 往往变差；
- 当前不适合作为默认推荐。

### 8.2 `replicatedMeshOverlapSizeExchange true`

该实验在普通 Alltoall 路径中：

- pack 后先发起 `MPI_Ialltoall`；
- 用本地 delete/local prep 阶段尝试覆盖 size exchange；
- 再 `MPI_Wait`。

`zb` 300-step repeat3 结果：

| 配置 | real mean [s] |
|---|---:|
| current optimum | 65.88 |
| overlapSizeExchange=true | 68.49 |

原因：

- 可重叠的本地工作窗口太小；
- MPI 非阻塞 collective 不保证后台推进；
- 实测 `wait` 反而上升。

因此这个方向可以作为默认关闭的实验代码保留，但不应纳入推荐配置。

### 8.3 跳过非输出 write / 继续压 check 之外的激进节流

本轮也验证过若干更激进的 squeeze 方向，例如：

- 改更稀的 `SARSteps=100`；
- 更大的 `MinGapSteps=100`；
- 更高 threshold；
- 跳过部分 post diagnostic。

这些都没有稳定超过当前最优，说明 `zb` 这一路径已经比较接近“可由控制项榨出的上限”。

## 9. 当前推荐状态和后续方向

### 9.1 当前推荐状态

对当前本地源码树，推荐分两层理解。

源码默认保留：

- owned collision-cell 遍历；
- collision/migration subphase profiling；
- adaptive alpha 基础设施；
- `replicatedMeshGatherCandidates`；
- `replicatedMeshOverlapSizeExchange`；
- `replicatedMeshDLBMinRemainingSteps`。

但 `zb` 当前推荐运行配置仍是：

```text
replicatedMeshNoAlltoall false;
replicatedMeshDLBDualConstraint false;
replicatedMeshDLBAlpha 1;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshSARSteps 50;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBParticleGate false;
replicatedMeshDLBAdaptiveAlpha false;
replicatedMeshDLBMinRemainingSteps 50;
replicatedMeshGatherCandidates false;
replicatedMeshOverlapSizeExchange false;
```

### 9.2 后续真正值得做的方向

这轮工作已经把“继续简单切配置”的收益空间基本榨干。若要再向前推进，应优先考虑：

1. collision-aware DLB
   - 直接把 cumulative local candidates / collision wall 纳入约束或目标，
     而不是继续只靠粒子数 proxy。
2. 分离 move-balance 与 collision-balance 判据
   - 避免 move/collision 反相关时，总时间平衡掩盖 collision 失衡。
3. MPI 内再引入小粒度 shared-memory 并行
   - 纯 MPI owned-cell 静态切分天然缺少 OMP 那种步内动态调度能力；
   - 若目标是尽量复现 OMP8 collision speedup，`MPI x small-OMP` 比继续 pure MPI
     更现实。
4. 若继续 pure MPI
   - 重点不再是补更多 collision kernel 微优化，而是减少 phase barrier、
     缩小 move/build arrival skew，或者让 DLB 真正看到 collision 负载。

## 10. 对本轮结果的最终判断

这轮 follow-up 已经把 replicated-mesh MPI8 的问题边界收紧得比较清楚：

- collision kernel 本身已经不再是主要未知项；
- OMP8 对 collision 的大幅加速，本质上来自共享内存步内动态调度；
- pure MPI8 的剩余问题主要是静态 owner 切分、collision 候选负载偏斜，以及
  move/build/collision 之间的同步等待；
- `zb` 当前最优配置已经被 repeat3 和回归复测确认，但进一步收益不会再主要来自
  简单配置切换。

因此，后续若继续追求“真正的 mpi8 replicated-mesh coll 和总时间加速”，
应从 collision-aware DLB 或 MPI+OMP 混合分工入手，而不是继续期待单纯 cadence、
Dual 或 alpha 微调能复制 OMP8 的 collision 行为。
