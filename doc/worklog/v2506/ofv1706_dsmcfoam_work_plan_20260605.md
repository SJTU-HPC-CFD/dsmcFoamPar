# OFv1706 dsmcFoam+ 性能优化与 DLB 复现工作计划

日期：2026-06-05

## 1. 目标和准则

本轮工作的目标是在

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

中复现并继续推进 OF-v2506 版本 `dsmcFoam+` 的性能优化工作，最终覆盖三种运行形态：

1. `omp8`：OpenMP 8 线程优化和 OpenMP 侧负载平衡能力。
2. `mpi8replicatedmesh`：MPI 8 rank replicated mesh DLB。
3. `mpi+omp`：MPI+OpenMP 混合并行。

执行准则：

- `doc/worklog/v2506` 是历史背景和决策记录。
- `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx` 是参考代码源，后续实现以它的实际源码为准。
- 不修改 `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`。
- 每个阶段必须保留完整日志和对应 `controlDict`。
- 每个测试必须检查粒子数、碰撞数、能量一致性。
- 正式性能比较默认使用 `profileSummary true`、`profileDetail false`。

## 2. 当前已知状态

### 2.1 当前 OFv1706 工作树

当前 OFv1706 树中已存在：

```text
src/lagrangian/basic/Cloud/Cloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C
src/lagrangian/dsmc/collisions/derived/VariableHardSphere/VariableHardSphere.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
src/lagrangian/dsmc/dynamicLoadBalancing/dsmcDynamicLoadBalancing.C
```

当前 OFv1706 树中未发现：

```text
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H
```

这说明 OFv1706 当前不能直接按 v2506 replicated mesh DLB 日志运行，必须先做源码级移植和编译集成。

### 2.2 参考代码中确认存在的关键模块

参考代码 `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx` 中确认存在：

```text
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H
src/lagrangian/basic/Cloud/Cloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C
src/lagrangian/dsmc/collisions/derived/VariableHardSphere/VariableHardSphere.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
src/lagrangian/dsmc/dynamicLoadBalancing/dsmcDynamicLoadBalancing.C
```

参考代码中还确认存在这些能力：

- `replicatedMeshActive()` / `dsmcReplicatedMesh`。
- `openmpMoveSchedule`、`openmpMoveChunk`、`openmpCollisionSchedule`、`openmpCollisionChunk`。
- `profileDetail` 控制。
- `moveOrderedParcels_` / flat occupancy 相关路径。
- `collisionFastRng` 和 `FastRng`。
- replicated mesh 输出侧的 `isMyCell` 和 MPI 汇总逻辑。

## 3. 基准和测试目录

主测试 case：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh
```

基准日志和恢复来源：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp
```

阶段记录目录：

```text
doc/worklog/v2506/detail
```

正式比较至少保留三类数据：

- external wall / main loop / full evolve。
- move only / move kernel / buildCellOccupancy / collision / post fields/output。
- 粒子数、碰撞数、平均总能量、异常标记检查。

## 4. 阶段 0：环境、源码、基准和完整性能剖析

目标：在做任何功能移植前，把当前 OFv1706、参考代码、已有基准日志三者对齐，并先完成当前代码的完整性能剖析。后续所有优化必须能回到阶段 0 的剖析数据上做 before/after 量化，否则不进入默认保留路径。

步骤：

1. 激活环境：

   ```bash
   source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
   ```

2. 编译当前基线：

   ```bash
   source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
   ```

3. 建立源码对照清单：

   - `Cloud.C`
   - `dsmcCloud.C/H/I.H`
   - `noTimeCounter.C/H`
   - `VariableHardSphere.C/H`
   - `dsmcVolFields.C/H`
   - `dsmcDynamicLoadBalancing.C/H`
   - `dsmcFoam+.C`
   - `Make/files`、`Make/options`

4. 解析 `ourmeshbkp` 三组基准日志：

   - `omp8`
   - `mpi8origin`
   - `mpi8replicatedmesh`

5. 建立当前 OFv1706 的运行时调用链性能地图：

   ```text
   dsmcFoam+.C
     -> dsmcCloud::evolve()
       -> controllers/boundaries controlBeforeMove
       -> Cloud<dsmcParcel>::move()
       -> buildCellOccupancy()
       -> coordSystem().evolve()
       -> collisions()
       -> reactions / post-collision controls
       -> fields / properties / output
       -> loadBalanceCheck()
   ```

6. 做完整性能剖析，至少拆出以下层级：

   - solver 总体：external wall、main loop、full evolve、write、non-DSMC wall。
   - move：pre-control/partition、reset/setup、extract parcels、parallel kernel、transfer/delete finalize、commit/survivor rebuild。
   - buildCellOccupancy：extract/gather、count/reduce、allocate/fill、fallback 次数。
   - collision：candidate 数、active cell 数、selection/collide、reaction、per-rank max/avg/sum。
   - post fields/output：field calculate、parcel accumulate、field combine、boundary accumulation、writeField。
   - MPI：每 rank compute/comm/wait、move transfer 总量、DLB check/trigger/rebalance/migration/post 诊断。
   - OpenMP：每线程粒子数、wall time、threadWorkMax/threadWorkAvg、fallback reason。
   - 正确性：最终粒子数、last-step collisions/candidates、平均总能量、异常标记。

7. 区分两种剖析口径：

   - clean baseline：`profileSummary true`、`profileDetail false`，用于性能结论。
   - diagnostic run：必要时打开 detail 或临时计时，只用于定位热点，不作为性能 baseline。

8. 如果 OFv1706 当前缺少必要计时点，阶段 0 可以先加入低扰动、可开关的 profiling instrumentation，但这些计时必须满足：

   - 默认关闭或由 `profileSummary/profileDetail` 控制。
   - 不改变物理路径和并行语义。
   - 先在小步数 case 上验证不引入明显 wall time 污染。
   - 后续正式性能比较仍以 clean baseline 为准。

9. 生成可复用的日志解析脚本或命令片段，输出统一表格：

   - 三方法总览表。
   - 单 case 分项时间表。
   - MPI per-rank 表。
   - OpenMP per-thread 表。
   - 正确性检查表。

10. 将阶段 0 结果写入：

   ```text
   doc/worklog/v2506/detail/stage0_baseline_source_audit_20260605.md
   doc/worklog/v2506/detail/stage0_full_performance_profile_20260605.md
   ```

完成条件：

- 当前 OFv1706 能编译。
- 明确当前已有功能、缺失功能和参考代码对应位置。
- 得到可复核的三方法基准表。
- 得到当前代码的完整性能剖析表，覆盖 move、buildCellOccupancy、collision、post、MPI、OpenMP 和 correctness。
- 明确后续每一阶段的主要优化目标、Amdahl 上限和预期收益来源。
- 后续优化结果能直接和阶段 0 指标做量化对比。

## 5. 阶段 1：最小 OpenMP 性能基础移植

目标：先让 OFv1706 具备 v2506 中已经证明有效的 OpenMP 基础能力，避免一上来就引入 replicated mesh 的复杂度。

优先移植内容：

1. `profileSummary` / `profileDetail` 生产口径。

   - 正式性能比较默认 `profileDetail false`。
   - detail profiling 只用于归因，不作为 clean baseline。

2. `Cloud::move()` OpenMP 基础路径。

   - 并行 `particles[]` extract。
   - `p.reset()` 并入并行路径。
   - `move commit/survivor rebuild` 两阶段化。
   - 避免把高扰动 per-track 计数链放入生产热路径。

3. `buildCellOccupancy()` flat/ordered 数据流。

   - 先确认 OFv1706 当前 `cellOccupancy_` 的使用面。
   - 再移植 reference 中可闭合的 `moveOrderedParcels_` / flat occupancy 路径。
   - 保留旧 `cellOccupancy_` 作为兼容层，避免一次性改坏 collisions/reactions/boundaries/fields。

4. `post fields/output` 小幅稳定优化。

   - `dsmcVolFields` 的 per-cell/local accumulator。
   - 按 field 配置裁剪无用累积项。
   - 复用 occupancy 连续视图。

5. `collision` 基础优化。

   - `collisionFastRng`。
   - `VariableHardSphere::sigmaTcR()` cr2 查表优化是否可移植。
   - OpenMP 反应路径线程安全：`addNewParcel()` deferred buffer 和反应计数 atomic。

不默认移植的内容：

- 已证明负收益的 same-tet fast path。
- 高扰动 `moveTrackCallCount` / per-track profiling 生产路径。
- getter 级微缓存。
- 无明确收益的 collision partition 复杂化。

测试：

- 先跑小步数 smoke test。
- 再跑 `ourmesh/omp8` 正式对比。
- 每次保留日志和 `controlDict`。

完成条件：

- `omp8` 正确完成。
- 粒子数、碰撞数、能量与基准同量级。
- 记录 OpenMP 移植后分项时间变化。

## 6. 阶段 2：OpenMP DLB / 负载诊断能力

目标：建立 OpenMP 侧负载诊断和必要的线程级负载均衡基础，而不是盲目加入复杂 DLB。

步骤：

1. 复核 reference 中实际保留的 OpenMP 负载相关代码。
2. 提供低扰动的线程级统计：

   - thread move wall。
   - collision per-thread wall。
   - candidates / active cell 统计。

3. 只在 `profileDetail true` 下打开重型诊断。
4. 对 `openmpMoveSchedule` 和 `openmpCollisionSchedule` 做最小必要支持。

判断标准：

- 如果调度改变不能稳定降低 main loop，不保留为默认。
- 如果诊断本身改变性能，不纳入 clean baseline。

输出：

```text
doc/worklog/v2506/detail/stage1_omp_dlb_diagnostics_*.md
```

## 7. 阶段 3：MPI replicated mesh DLB 移植

目标：在 OFv1706 中建立 replicated mesh 运行形态，优先实现能正确运行和输出，再追求 DLB 性能。

核心移植内容：

1. 新增 `dsmcReplicatedMesh` 模块。

   - `dsmcReplicatedMesh.C/H`
   - `Make/files`
   - `Make/options`
   - METIS/ParMETIS 依赖检查。

2. 接入 `dsmcCloud`。

   - `replicatedMesh_` 成员。
   - `replicatedMeshActive()`。
   - `replicatedMeshRef()` / `replicatedMesh()` accessors。
   - `evolve()` 中 migration、build occupancy、collision、rebalance 的顺序。

3. 接入 `dsmcFoam+.C`。

   - replicated mesh 模式下不使用 `-parallel`。
   - 按参考代码处理 gather/write/migrate/cellOwner 输出。

4. 接入 `Cloud::move()`。

   - replicated mesh 模式下避免普通 OpenFOAM processor patch transfer 慢路径。
   - 保证 migration 后 `moveOrderedParcels_` 状态一致。

5. 输出正确性修复。

   - `dsmcVolFields` 对 replicated mesh 使用 owned-cell 采样。
   - 输出前 MPI_Allreduce。
   - reduce 后派生场 recomputation。
   - Tvib / species Evib 累积必须使用 `isMyCell` 过滤。

6. DLB 参数先采用当前稳定口径。

   - `replicatedMesh true`
   - `replicatedMeshDelayedReceive true`
   - `replicatedMeshNoAlltoall true`
   - `replicatedMeshFlatTransfer true`
   - `replicatedMeshMigrateInterval 10`
   - `replicatedMeshDLBDualConstraint true`
   - `replicatedMeshAutoDLB true`
   - `replicatedMeshDLBSteps 50`
   - `replicatedMeshDLBTriggerMode legacyWindow`
   - `replicatedMeshDLBMinGapSteps 50`
   - `replicatedMeshDLBCheckCollective allgather`
   - `replicatedMeshDLBSkipFixedKPostDiag false`
   - `replicatedMeshDLBProfile false`

不作为默认方向：

- external reconstruct/decompose/restart mesh DLB。
- Allreduce low-check 默认化。
- QualityGate / ROI gate。
- skip-minGap-closed SAR collective。
- async SAR overlap。
- weighted decomposition 作为默认初始方案。

测试顺序：

1. 2 rank smoke test。
2. `mpi8replicatedmesh` 小步数正确性测试。
3. `mpi8replicatedmesh` 正式基准对比。
4. 输出字段检查：非零 cell 数、温度场、粒子数、能量。

输出：

```text
doc/worklog/v2506/detail/stage2_replicated_mesh_port_*.md
```

完成条件：

- replicated mesh 模式能不加 `-parallel` 运行。
- DLB rebalances 数量和触发日志可解释。
- 输出字段不出现 owned-cell 缺失或 Tvib 异常。
- 相比 `mpi8origin` 有明确性能收益，或能明确定位未达成原因。

## 8. 阶段 4：MPI+OpenMP 混合并行

目标：在 replicated mesh 或普通 MPI 形态上启用 OpenMP，使 mixed 模式具备可测性能。

移植和验证重点：

1. mixed move 数据流。

   - `stagebuf1`：`moveAppendedParcels_` 只表示 pre-move inflow。
   - `stagebuf2`：`moveAppendToPending_` 在整个 move 多 pass 中持续有效。
   - 不引入已验证破坏 `moveOrderedParcels_ -> buildCellOccupancy` 的双缓冲或过细列表拆分。

2. OpenMP guard cell / processor boundary 交互。

   - 先保持参考代码默认策略。
   - 不盲目启用 guided/dynamic move，除非能证明 guard cell 语义正确。

3. mixed collision。

   - 默认使用 reference 中当前有效的 `openmpCollisionSchedule`。
   - FastRng 和 sigmaTcR 优化必须保持线程安全。
   - 反应 case 必须使用 deferred `addNewParcel()`。

4. mixed 输出和 DLB。

   - 复用阶段 3 的 replicated mesh 输出规则。
   - 若普通 MPI 模式下测试 mixed，必须区分 `Pstream::parRun()` 和 replicated raw MPI 模式。

测试：

- 小步数 `mpi2 x omp4` 或 `mpi4 x omp2` smoke。
- `ourmesh` mixed 正式测试。
- 与 `omp8` 和 `mpi8replicatedmesh` 比较。

输出：

```text
doc/worklog/v2506/detail/stage3_mpi_omp_mixed_*.md
```

完成条件：

- mixed 模式稳定完成。
- 无链表损坏、无 PstreamBuffers 未消费错误、无反应 OpenMP crash。
- 性能分项能解释与 `omp8` / `mpi8replicatedmesh` 的差距。

## 9. 阶段 5：统一性能确认

目标：以同一套源码、同一批 case、同一性能口径做最终确认。

正式对比对象：

```text
ourmesh/omp8
ourmesh/mpi8origin
ourmesh/mpi8replicatedmesh
ourmesh/<mixed-case>
```

每个 case 保存：

- 完整 solver log。
- `/usr/bin/time` 输出。
- `system/controlDict`。
- `system/loadBalanceDict`。
- 当前 git diff 或源码改动摘要。

必须统计：

```text
external wall
main loop wall
full evolve wall
move only
move parallel kernel
move transfer/delete finalize
buildCellOccupancy
collision phase
post fields/output
DLB trigger/rebalance/check/migration/post diagnostics
final Number of molecules
last-step Collisions
last-step Collision candidates
Average total energy
End / FOAM FATAL / MPI_ABORT / nan / segfault 标记
```

最终总结写入：

```text
doc/worklog/v2506/detail/final_three_modes_performance_*.md
```

并在 `doc/worklog/v2506/` 下补一份总览总结。

## 10. 风险和回退策略

### 10.1 主要风险

1. OFv1706 和 OF-v2506 API 差异。

   - `Cloud.C` 模板接口可能不同。
   - `PstreamBuffers` / MPI 包装接口可能不同。
   - `fvMesh` / `polyMesh` 构造和分发接口可能不同。

2. replicated mesh 依赖较重。

   - 需要 METIS/ParMETIS。
   - Make/options 需要确认库路径。
   - raw MPI 与 OpenFOAM Pstream 必须分清。

3. correctness 风险。

   - 反应 OpenMP 路径可能产生链表 race。
   - replicated mesh 输出容易漏 reduce 或漏 `isMyCell`。
   - DLB 后粒子数/碰撞数/能量必须逐步验证。

4. profiling 污染性能。

   - `profileDetail true` 只能用于定位。
   - clean baseline 必须 `profileDetail false`。

### 10.2 回退策略

- 每个阶段完成后保留可编译 checkpoint。
- 每个负优化必须记录并回退。
- 如果某个移植点阻塞编译，先拆成更小补丁，不跨模块堆叠。
- 如果 correctness 失败，优先回到最近通过的 correctness checkpoint，而不是继续优化性能。

## 11. 建议的执行顺序

建议按以下顺序执行：

1. 阶段 0：环境、源码、基准和完整性能剖析。
2. 阶段 1：OpenMP 基础优化移植。
3. 阶段 2：OpenMP 负载诊断和低扰动 DLB 支撑。
4. 阶段 3：MPI replicated mesh DLB。
5. 阶段 4：MPI+OpenMP mixed。
6. 阶段 5：统一性能确认和总结。

原因：

- 当前 OFv1706 缺 replicated mesh 模块，不能直接从 MPI DLB 开始。
- OpenMP/flat occupancy 是 replicated mesh 和 mixed 的基础数据流。
- 输出正确性必须在 replicated mesh 性能调参之前解决。
- mixed 模式最容易同时触发 OpenMP、MPI、replicated mesh 三类问题，应该放在前两类能力稳定后做。

## 12. 近期下一步

下一步建议执行阶段 0：

1. 编译当前 OFv1706 基线。
2. 提取 `ourmeshbkp` 三组基准日志的关键指标。
3. 对当前 OFv1706 代码做完整性能剖析，先形成 clean baseline，再按需形成 diagnostic profile。
4. 生成 OFv1706 与 `hyStrath_xcx` 的关键源码差异表。
5. 写入：

   ```text
   doc/worklog/v2506/detail/stage0_baseline_source_audit_20260605.md
   doc/worklog/v2506/detail/stage0_full_performance_profile_20260605.md
   ```

阶段 0 完成后，再决定第一批最小源码移植点。
