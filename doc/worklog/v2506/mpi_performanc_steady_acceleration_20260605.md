# 当前有效优化与源码核对技术总结

日期：2026-06-05

本文根据 `hyStrath_xcx/docs/detail` 下已有工作日志、`ourmesh` 三方法最新确认日志，以及当前工作树源码/配置核对，梳理已经证明“实际有效并仍在使用”的工作。核对对象是当前工作树文件内容，不以 git index/staged 中早期残留为准。

## 1. 当前稳定性能结论

最终以 `ourmesh_three_methods_perf_20260605` 的 Run 3 confirmation 为准。三套正式 case 均为 `profileSummary true` / `profileDetail false`。

| 方法 | 运行形态 | external wall [s] | main loop [s] | full evolve [s] | move only [s] | move kernel [s] | collision [s] | DLB |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `omp8` | OpenMP 8 | 103.66 | 104.8437967 | 102.1225033 | 48.82365912 | 40.14260107 | 9.876225069 | 0 |
| `mpi8replicatedmesh` | MPI 8 replicated mesh, no `-parallel` | 117.31 | 118.3501321 | 114.9207628 | 47.1652636 | 45.76155029 | 10.89824478 | 8 rebalances |
| `mpi8origin` | MPI 8 decomposed mesh, `-parallel` | 160.95 | 166.6245026 | 163.0536101 | 120.3658396 | 56.6983455 | 18.73261305 | off |

最新排序明确：

```text
omp8               103.66 s
mpi8replicatedmesh 117.31 s
mpi8origin         160.95 s
```

三条 Run 3 日志均达到 `End`，未检出 `FOAM FATAL`、`MPI_ABORT`、`Abort(`、`nan` 或段错误标记。最终粒子数、碰撞数和平均总能量也在日志中记录。完整日志和对应 `controlDict` 已统一保存到 `docs/detail/ourmesh_three_methods_perf_20260605/confirm_pdFalse_20260605_024610_bundle/`：

```text
mpi8origin.log
mpi8origin.controlDict
mpi8replicatedmesh.log
mpi8replicatedmesh.controlDict
omp8.log
omp8.controlDict
```

## 2. 实际有效且仍在使用的工作

### 2.1 生产运行统一为 `profileDetail false`

**日志证据**

`ourmesh_three_methods_perf_20260605/README.md` 记录 Run 2 开始三方法统一为 `profileDetail false`；Run 3 复测前再次确认三套 production case 已经是 `profileSummary true` / `profileDetail false`。

统一关闭详细 profiling 后，Run 2 中两个 MPI 方法相对 Run 1 明显下降：

| 方法 | Run 1 wall [s] | Run 2 wall [s] | 变化 |
|---|---:|---:|---:|
| `mpi8origin` | 176.45 | 159.69 | -9.50% |
| `mpi8replicatedmesh` | 141.97 | 127.37 | -10.28% |

`omp8` 在 Run 1 已经是 `profileDetail false`，所以 Run 1 到 Run 2 的小差异不能解释为配置收益。

**当前配置证据**

当前三套正式 case 均仍保持该设置：

```text
case/.../ourmesh/mpi8origin/system/controlDict: profileSummary true; profileDetail false;
case/.../ourmesh/mpi8replicatedmesh/system/controlDict: profileSummary true; profileDetail false;
case/.../ourmesh/omp8/system/controlDict: profileSummary true; profileDetail false;
```

**当前代码证据**

`src/lagrangian/dsmc/parcels/dsmcParcel.C` 的 `dsmcParcel::move()` 先读取 `cloud.profilingDetailEnabled()`。当 `profileDetail false` 时，代码直接执行普通 tracking/control 流程，不进入详细计时、move-iteration 累加、face/patch hit 计数等路径。`moveItersPerCell_` 和 `moveItersPerCellCumulative_` 更新只在详细 profiling 分支中发生。

**结论**

这是当前生产性能口径的前提配置。性能比较必须使用 `profileDetail false`，详细 profiling 日志只能用于归因，不能作为 clean baseline。

### 2.2 `replicatedMeshActive()` 语义收紧

**日志证据**

`omp8_activefix_followup_20260605.md` 记录：将若干单 rank 敏感的 `replicatedMesh_.valid()` 检查改为 `replicatedMeshActive()` 后，`omp8` 恢复了一部分慢速回归：

```text
external wall: 125.46 -> 118.35 s
main loop:     127.2241534 -> 120.6044536 s
move kernel:   61.01169953 -> 58.98261391 s
move only:     72.52985026 -> 69.50859184 s
```

这不是最终最大收益，但它是正确性/运行形态修复：避免单进程 OpenMP case 因为存在 inactive replicated mesh 对象而误入 replicated migration/report/rank-time 路径。

**当前代码证据**

`src/lagrangian/dsmc/clouds/dsmcCloudI.H`：

```cpp
return replicatedMesh_.valid() && replicatedMesh_->active();
```

`src/lagrangian/dsmc/clouds/dsmcCloud.C` 的 evolve/reporting 路径继续使用 `replicatedMeshActive()` 保护初始迁移、延迟接收、迁移、rebalance、rank-time/report 等 replicated mesh 分支。

**当前配置关联**

`omp8/system/controlDict` 里仍有 `replicatedMesh true` 和 `replicatedMeshAutoDLB true`，但最新 Run 3 中 OMP DLB 统计为 0。当前代码通过 `replicatedMeshActive()` 区分“配置存在”和“并行 replicated mesh 实际激活”，这正是该修复的价值。

**结论**

保留。它是 OpenMP/单 rank 不误走 replicated 路径的基础修复，也避免后续 profiling 或 report 误判。

### 2.3 移除 `moveTrackCallCount` profiling 链

**日志证据**

`omp8_activefix_followup_20260605.md` 明确记录：与 `src/gitbkp` 对照后，发现额外 `moveTrackCallCount` profiling 链是 OMP move 回归主因。移除后：

```text
external wall: 118.35 -> 103.77 s
main loop:     120.6044536 -> 103.4617452 s
full evolve:   117.680763 -> 100.914042 s
move only:      69.50859184 -> 48.20069174 s
move kernel:    58.98261391 -> 39.593751 s
```

随后正式 `ourmesh` Run 3 也确认 `omp8` 回到最快：

```text
Run 2 omp8 wall: 125.46 s
Run 3 omp8 wall: 103.66 s
```

**当前源码核对**

当前源码检索未发现以下符号：

```text
moveTrackCallCount
trackNsPerCall
moveThreadTrackCallCounts
```

核对范围包括：

```text
src/lagrangian/basic/Cloud/Cloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.H
src/lagrangian/dsmc/clouds/dsmcCloudI.H
src/lagrangian/dsmc/parcels/dsmcParcel.C
src/lagrangian/dsmc/parcels/dsmcParcel.H
```

**结论**

这是本轮把 OMP 性能拉回来的最大单项有效工作。后续不应把 per-track/per-call 计数重新放回生产 move 热路径。

### 2.4 `profileDetail false` 下的 parcel false-path 整理

**日志证据**

`omp8_activefix_followup_20260605.md` 记录：`dsmcParcel::move()` 拆分后，`profileDetail false` 绕过详细 move timing/counting work。这个改动只带来小幅改善：

```text
move kernel: 58.98261391 -> 58.39925487 s
move only:   69.50859184 -> 68.68627768 s
```

它不是最终主因，但方向正确，并与生产 `profileDetail false` 口径一致。

**当前代码证据**

`src/lagrangian/dsmc/parcels/dsmcParcel.C` 当前保留：

- `const bool profileMoveDetail = cloud.profilingDetailEnabled();`
- `if (!profileMoveDetail)` 直接走普通 free/stuck parcel tracking/control；
- `moveItersPerCell_` 只在 detail 分支里更新；
- `hitPatch()` / `hitWallPatch()` 在 `profileDetail false` 下只执行必要 boundary control，不做计时累加。

**结论**

保留。它是小幅生产路径清理，不应夸大为 OMP 回归主因；主因仍是 `moveTrackCallCount` 链移除。

### 2.5 replicated mesh DLB 当前有效组合

**日志证据**

replicated mesh/DLB 的有效线索不是单个参数，而是一组“可运行且避免慢路径”的组合。

`replicated_mesh_retest_20260602.md` 中的关键结论：

- `legacyWindow` 自动触发模式是 post-restart 后最有效的恢复方向之一，把 case 拉回约 130 s 级别；
- `lowcheck_reduce` 的 `Allreduce` 检查实验是负结果，虽然减少了显式 post diagnostic，但改变 DLB 轨迹并导致 `141.86 s`；
- 显式恢复 `allgather` 和 `replicatedMeshDLBSkipFixedKPostDiag false` 后，case 回到 130 s 级别；
- `replicatedMeshDLBMinGapSteps 50` 有效消除了 `90 -> 120` 这类短间隔 rebound，三次 minGap 复测从旧慢路径的 `140+ s` 收敛到约 `134 s` 均值，并将 DLB count 限制在 7-8 左右；
- 方向 2 的 skip-minGap-closed SAR collective 降低了 `check max`，但改变轨迹并未带来端到端收益，已回退。

后续正式 Run 3 在当前源码和 `profileDetail false` 下给出更强最终结果：

```text
mpi8replicatedmesh external wall = 117.31 s
mpi8origin external wall         = 160.95 s
```

**当前配置证据**

`mpi8replicatedmesh/system/controlDict` 当前保留：

```text
replicatedMesh true;
replicatedMeshDelayedReceive true;
replicatedMeshNoAlltoall true;
replicatedMeshFlatTransfer true;
replicatedMeshMigrateInterval 10;
replicatedMeshDLBDualConstraint true;
replicatedMeshAutoDLB true;
replicatedMeshDLBSteps 50;
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshDLBSkipFixedKPostDiag false;
replicatedMeshDLBImbalanceThreshold 1.5;
replicatedMeshDLBFixedK 1;
replicatedMeshDLBUbvec 1.05;
replicatedMeshDLBUbvec1 1.5;
replicatedMeshDLBProfile false;
replicatedMeshDLBInitialAlpha 0.8;
replicatedMeshDLBVsizeExp 0;
```

`mpi8origin/system/controlDict` 当前为 `replicatedMesh false` / `replicatedMeshAutoDLB false`，作为非 replicated MPI 对照。

**当前代码证据**

`src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C` 当前仍读取并实现：

- `replicatedMeshDLBTriggerMode`，默认 `cumulative`，当前 case 用 `legacyWindow`；
- `replicatedMeshDLBMinGapSteps`；
- `replicatedMeshDLBCheckCollective`，可选 `allgather` / `allreduce`，当前 case 用 `allgather`；
- `replicatedMeshDLBSkipFixedKPostDiag`，当前 case 为 `false`；
- `replicatedMeshDLBForceSteps` 只作为可选 forced cadence，不在当前正式 case 中启用。

`legacyWindow` 分支中仍按 `ndecps_ % sarSteps_ == 0` 做 SAR 采样，并用：

```text
triggered = minGapSatisfied && (sarTriggered || thresholdTriggered)
```

完成当前 minGap 保护下的 `SAR || threshold` 触发。

**结论**

这是 MPI replicated mesh 当前仍在使用的有效 DLB 组合。它的价值是避免非 replicated MPI 的重 move 代价，并把 DLB 轨迹从不稳定慢路径拉回可用范围。最新 117.31 s 结果还叠加了后续 profiling 热路径清理和 `profileDetail false` 生产口径。

### 2.6 `dsmcVolFields` per-cell local accumulator

**日志证据**

`dsmcVolFields_percell_20260604/README.md` 的决策：

- 保留 per-cell local accumulator；
- 回退 per-type constant cache；
- 判断口径以 `parcel accumulate` 和 `fields.calculate` 目标热点为主，不把单次总 wall 当作局部优化判据。

三次稳定复测记录在 `dsmcVolFields_percell_stability_20260604/README.md`：

| 指标 | 基线摘要 | per-cell 三次均值 | 均值相对基线 |
|---|---:|---:|---:|
| `parcel accumulate` | 23.141 | 22.307 | -3.61% |
| `fields.calculate` | 30.602 | 29.765 | -2.73% |

回退第1项后的三次均值更差：

| 指标 | 回退后三次均值 | 相对 per-cell |
|---|---:|---:|
| `parcel accumulate` | 24.499 | +9.83% |
| `fields.calculate` | 31.868 | +7.07% |

**当前源码核对**

`src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C` 当前仍有：

- per-cell local scalar/vector accumulators；
- `activeTypes`；
- parcel loop 中直接读取 `cloud.constProps(typeId)`；
- 对 `activeTypes` 写回每个 type/cell 的 field；
- 未发现已回退的 per-type constant cache 符号，例如 `massByType`、`thetaVByType`、`electronicEnergyByType`。

**结论**

保留，但应表述为 fields/post 侧小幅稳定热点优化，不能写成端到端 10% 以上稳定加速，也不是 `omp8` 从 125 s 回到 103 s 的主因。

### 2.7 `particle` tracking 中非 debug `position()` 计算移除

**日志证据**

`move_tracking_optimization_20260604/README.md` 记录：移除 `trackToStationaryTri()` / `trackToMovingTri()` 中非 debug 的多余 `position()` 计算是可接受优化。100-step clean comparison：

| run | move only [s] | move kernel [s] | full evolve [s] | external wall [s] |
|---|---:|---:|---:|---:|
| old baseline, fast path off | 8.187566989 | 7.618559921 | 18.92929133 | 24.72 |
| optimized, fast path off | 7.532625702 | 6.954936849 | 18.1020032 | 23.90 |

对应改善：

```text
move only              8.00% faster
move parallel kernel   8.62% faster
full evolve            4.37% faster
external wall          3.32% faster
```

**当前源码核对**

`src/lagrangian/basic/particle/particle.C` 中：

- `trackToStationaryTri()` 只有在 `debug` 时才计算 `x0 = position()` 并打印；
- `trackToMovingTri()` 同样只在 `debug` 时计算 `x0 = position()`；
- 非 debug 常规 tracking 不再为日志输出提前构造全局位置。

**结论**

保留。它是 move tracking 常规路径的低风险小优化，但需要与 `moveTrackCallCount` 大收益区分开。

## 3. 保留但主要用于诊断的工作

### 3.1 move detail / hotspot profiling

`move_detail_profile_full_20260604/README.md` 和 `move_cell_hotspot_profile_20260604/README.md` 的主要价值是归因：

- 500-step detail run 显示 `trackToAndHitFace` 占 observed move time 约三分之二；
- boundary control 仅约 `0.228 s`，不是主要瓶颈；
- hotspot 100-step 报告显示 move/track iterations 不集中在少数 top cells，也不是前三层 DSMC boundary layer 主导。

当前源码仍保留 `moveItersPerCell_` / `moveItersPerCellCumulative_` 和 hotspot report，但这些只在 `profileDetail true` 下工作。生产 case 当前为 `profileDetail false`，因此不计入 clean 性能收益。

**结论**

保留为诊断工具。不要把 detail profiling run 当作性能 baseline。

## 4. 已验证无效、已撤回或未作为当前默认的工作

| 工作项 | 日志结论 | 当前核对结论 |
|---|---|---|
| same-current-tet fast path (`moveSameTetFastPath`) | 5-step 能减少 track calls，但 100-step/500-step clean run 端到端变慢；500-step fast-on external wall `129.46 -> 138.92 s`，DLB 从 5 次增至 8 次 | 当前源码/正式 case 检索无 `moveSameTetFastPath` / `sameTetFastPath` / `SameTetFastPath` |
| `hitPatch()` / `hitWallPatch()` 额外 boundary-hit fastpath | `log.activefix_fastpath_unsandbox_20260605` 比 kept fix 慢，已撤回 | 当前只保留 `profileDetail false` 下不做计时/计数的普通 boundary control，不把额外 fastpath 作为有效项 |
| `SAR && threshold` 触发 | DLB 次数减少但 wall 变慢；已回退到 `SAR || threshold` 带 minGap/eligibility | 当前 `legacyWindow` 分支仍是 `sarTriggered || thresholdTriggered` |
| single-constraint / 过紧或过松 ubvec / adaptive alpha 等配置扫参 | 多数不优于 best-known 或不稳定 | 当前正式 case 使用 dual constraint、`ubvec 1.05`、`ubvec1 1.5`、`fixedK 1` |
| QualityGate / ROI gate / reject cooldown | 机制上可运行，但 proxy 过乐观或全拒绝，端到端不稳定；reject cooldown 也不足以修复无有效候选的问题 | 当前源码检索无 `replicatedMeshDLBQualityGate` / `RejectCooldown` |
| `replicatedMeshDLBUseSAR false` no-SAR threshold-only | `threshold=1.5` 太晚，`threshold=1.2` 太频繁，`DLBSteps=25` 仍不如 stronger comparison | 当前源码/正式 case 不使用该开关 |
| low-check `Allreduce` 默认化 | `lowcheck_reduce_20260604_1` 变慢到 `141.86 s`，并改变触发轨迹 | 源码保留 opt-in `allreduce`，但当前正式 case 显式 `replicatedMeshDLBCheckCollective allgather` |
| skip minGap-closed SAR collectives | 降低 check max，但改变 DLB trajectory，端到端没有稳定收益；已回退 | 当前源码检索未见 skipgap/prewarm 运行路径；正式 case只保留 minGap50 |
| async SAR overlap | 日志中显示主要是计时归属迁移，未带来相同轨迹下端到端收益 | 当前源码/正式 case 未检出 `replicatedMeshDLBAsyncSAR`，不作为当前有效默认 |
| per-type constant cache | `parcel accumulate` 不优于 per-cell local accumulator，已回退 | 当前 `dsmcVolFields.C` 无 `massByType`、`thetaVByType`、`electronicEnergyByType` 等符号 |

## 5. 当前推荐运行形态

### 5.1 环境和日志

继续遵守 `docs/detail/dsmcFoam++tips.md`：

```bash
source /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/docs/env.sh
```

保留完整 `dsmcFoam+` 日志。每次性能结论至少记录：

- external wall；
- `main loop wall time`；
- `full evolve wall`；
- `move only` / `move parallel kernel wall`；
- `buildCellOccupancy`；
- `collision phase`；
- `post fields/output`；
- DLB trigger steps、rebalances、check/migration/post diagnostic；
- 最终 particles、last-step collisions/candidates、average total energy；
- `End` 和异常标记检查。

### 5.2 三方法命令

`mpi8origin`：

```bash
cd case/pal-phd-ch3.3.1-cylinder-react/timecompare/ourmesh/mpi8origin
mpirun -np 8 dsmcFoam+ -parallel
```

`mpi8replicatedmesh`：

```bash
cd case/pal-phd-ch3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh
mpirun -np 8 dsmcFoam+
```

replicated mesh 模式不要运行 `decomposePar`，不要加 `-parallel`。

`omp8`：

```bash
cd case/pal-phd-ch3.3.1-cylinder-react/timecompare/ourmesh/omp8
OMP_NUM_THREADS=8 dsmcFoam+
```

### 5.3 配置底线

性能比较默认使用：

```text
profileSummary true;
profileDetail false;
```

`mpi8replicatedmesh` 当前有效 DLB 默认：

```text
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshDLBSkipFixedKPostDiag false;
replicatedMeshDLBProfile false;
```

不要把 `profileDetail true`、same-tet fast path、QualityGate/ROI、Allreduce low-check、skipgap 或 async SAR overlap 当成当前生产默认。

## 6. 总结

当前性能稳定不是来自单个大改动，而是几类工作叠加：

1. `mpi8replicatedmesh` 的基础收益来自 replicated mesh 运行形态和 DLB 轨迹控制：no `-parallel`、no `decomposePar`、dual-constraint fixedK DLB、`legacyWindow`、`minGap50`、`allgather` 检查路径。
2. `omp8` 当前最快，主要因为移除了 `moveTrackCallCount` 这条侵入 OpenMP move 热路径的 profiling 链；`replicatedMeshActive()` 和 `profileDetail false` false-path 清理是必要但较小的修复。
3. `profileDetail false` 是所有正式性能比较的前提。详细 profiling 能解释热点，但本身会改变性能。
4. `dsmcVolFields` per-cell local accumulator 和 `particle` 非 debug `position()` 移除是仍在源码里的小幅有效优化，但不应被写成三方法排序改变的主因。
5. 已回退/负结果方向已经明确：same-tet fast path、ROI/QualityGate、Allreduce low-check 默认化、skipgap、async SAR overlap、额外 boundary fastpath 和 per-type constant cache 都不能作为当前有效优化宣传。

因此，当前最可靠的技术表述是：

> 当前源码和 `profileDetail false` 配置下，`omp8` 通过移除 move profiling 热路径开销恢复为最快；`mpi8replicatedmesh` 通过 replicated mesh 运行形态、DLB 轨迹控制和同一批 profiling 清理稳定优于传统 `mpi8origin`；`mpi8origin` 仍受非 replicated MPI move/通信路径拖累，是三者中最慢。
