# OFv1706 dsmcFoam+ 全量执行计划

## Summary

目标是在 `/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb` 中，按现有 `doc/worklog/v2506/ofv1706_dsmcfoam_work_plan_20260605.md` 推进完整工作：先完成阶段 0 的源码审计、基准整理和完整性能剖析，再依次实施 OpenMP 基础优化、OpenMP 负载诊断、MPI replicated mesh DLB、MPI+OpenMP mixed，最后做统一性能确认。

本文件是执行模式下的落盘计划。后续工作以参考代码 `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx` 的实际源码为准，参考代码只读。

## Key Changes And Implementation Plan

### 阶段 0：基准、源码审计、完整性能剖析

- 激活和编译：`source doc/scripts/env.sh`，`source doc/scripts/build-dsmcFoam.sh`。
- 若编译失败，记录失败日志到 `doc/worklog/v2506/detail/stage0_build_*.log`，只修复编译阻塞，不做性能优化。
- 生成源码差异表，覆盖 `Cloud.C`、`dsmcCloud.C/H/I.H`、`noTimeCounter.C/H`、`dsmcVolFields.C/H`、`VariableHardSphere.C/H`、`dsmcDynamicLoadBalancing.C/H`、`dsmcFoam+.C`，并核对参考代码中的 `dsmcReplicatedMesh.C/H`、OpenMP、FastRng、profile、replicated output 相关实现。
- 解析 `ourmeshbkp/omp8`、`ourmeshbkp/mpi8origin`、`ourmeshbkp/mpi8replicatedmesh` 现有基准日志，提取 external wall、main loop、full evolve、move、buildCellOccupancy、collision、post、DLB 次数、粒子数、碰撞数、能量和异常标记。
- 建立 clean baseline 与 diagnostic profile 的口径：clean 使用 `profileSummary true`、`profileDetail false`；detail 只用于热点定位。
- 输出 `doc/worklog/v2506/detail/stage0_baseline_source_audit_20260605.md` 与 `doc/worklog/v2506/detail/stage0_full_performance_profile_20260605.md`。

### 阶段 1：OpenMP 基础优化移植

- 从参考代码移植最小 OpenMP 基础能力：`profileSummary/profileDetail`、`Cloud::move()` 两阶段化、flat/moveOrdered occupancy、`dsmcVolFields` per-cell accumulator、`collisionFastRng`、`VariableHardSphere::sigmaTcR()` cr2 查表、反应路径 deferred `addNewParcel()` 和 reaction counter atomic。
- 不移植默认负优化：same-tet fast path、生产路径 `moveTrackCallCount`、getter 微缓存、过细 profiling。
- 阶段输出 `doc/worklog/v2506/detail/stage1_omp_base_port_*.md`，并保留每次编译日志、case log、controlDict。

### 阶段 2：OpenMP 负载诊断与 OMP-DLB 支撑

- 实现低扰动线程级诊断：thread move wall、thread collision wall、threadWorkMax/threadWorkAvg、candidate/active-cell 分布、fallback reason。
- 重型统计只在 `profileDetail true` 或明确诊断开关下启用。
- 支持 `openmpMoveSchedule`、`openmpMoveChunk`、`openmpCollisionSchedule`、`openmpCollisionChunk`。
- 调度变化必须降低 clean baseline main loop 才能保留默认。
- 输出 `doc/worklog/v2506/detail/stage1_omp_dlb_diagnostics_*.md`。

### 阶段 3：MPI replicated mesh DLB

- 移植 `dsmcReplicatedMesh.C/H`，接入 `Make/files`、`Make/options`，确认 METIS/ParMETIS 链接。
- 接入 `dsmcCloud`：`replicatedMesh_`、`replicatedMeshActive()`、migration/build occupancy/collision/rebalance 顺序、`cellOwner` 输出。
- 接入 `dsmcFoam+.C`：replicated mesh 模式不加 `-parallel`，支持 gather/write/migrate/cellOwner 流程。
- 接入 `Cloud::move()`：replicated mesh 模式避免普通 processor patch transfer 慢路径，migration 后保持 `moveOrderedParcels_` 状态一致。
- 修复 replicated mesh 输出：owned-cell 采样、`MPI_Allreduce`、reduce 后派生场 recomputation、Tvib/species Evib 使用 `isMyCell`。
- 默认参数采用当前稳定组合：`replicatedMesh true`、`replicatedMeshDelayedReceive true`、`replicatedMeshNoAlltoall true`、`replicatedMeshFlatTransfer true`、`replicatedMeshMigrateInterval 10`、`replicatedMeshDLBDualConstraint true`、`replicatedMeshDLBTriggerMode legacyWindow`、`replicatedMeshDLBMinGapSteps 50`、`replicatedMeshDLBCheckCollective allgather`、`replicatedMeshDLBSkipFixedKPostDiag false`、`replicatedMeshDLBProfile false`。
- 输出 `doc/worklog/v2506/detail/stage2_replicated_mesh_port_*.md`。

### 阶段 4：MPI+OpenMP mixed

- 在 replicated mesh 或普通 MPI 形态稳定后启用 mixed。
- 移植并验证 mixed move 数据流：`stagebuf1` 中 `moveAppendedParcels_` 只表示 pre-move inflow；`stagebuf2` 中 `moveAppendToPending_` 在整个 move 多 pass 中持续有效。
- 不引入已验证破坏 `buildCellOccupancy` 快路径的双缓冲或过细列表拆分。
- 保持 OpenMP guard cell 策略与参考代码一致，不默认启用 guided/dynamic move。
- mixed collision 默认使用参考代码当前有效调度；FastRng、sigmaTcR 和反应线程安全必须同时验证。
- 输出 `doc/worklog/v2506/detail/stage3_mpi_omp_mixed_*.md`。

### 阶段 5：统一性能确认和总结

- 统一源码、统一 case、统一 `profileDetail false` 生产口径测试：`ourmesh/omp8`、`ourmesh/mpi8origin`、`ourmesh/mpi8replicatedmesh`、mixed case。
- 每个 case 保存 solver log、`/usr/bin/time` 输出、`controlDict`、`loadBalanceDict`、当前源码改动摘要。
- 最终报告包含三方法/四方法性能表、correctness 表、per-rank/per-thread 分析、DLB trigger/rebalance/migration 统计、负优化和回退列表。
- 输出 `doc/worklog/v2506/detail/final_three_modes_performance_*.md` 与 `doc/worklog/v2506/ofv1706_dsmcfoam_final_summary_*.md`。

## Test Plan

- 编译测试：每个阶段至少运行一次 `source doc/scripts/build-dsmcFoam.sh`。
- smoke 测试：每个新功能先短步数运行，确认无 crash、无 `FOAM FATAL`、无 `MPI_ABORT`、无 `nan`。
- 正式性能测试：OMP 用 `OMP_NUM_THREADS=8 dsmcFoam+`；MPI origin 用 `mpirun -np 8 dsmcFoam+ -parallel`；MPI replicated 用 `mpirun -np 8 dsmcFoam+`，不运行 `decomposePar`，不加 `-parallel`；mixed 按阶段 4 最终 case 配置运行。
- correctness 检查：final particles、last-step collisions/candidates、average total energy、replicated output 非零 cell、温度场和 Tvib 合理性。
- 性能验收：所有优化必须和阶段 0 clean baseline 做 before/after；“局部更平衡但 main loop 更慢”的改动记录为负优化并回退；detail profiling 只用于解释，不作为正式速度结论。

## Assumptions And Defaults

- 执行范围只限 `/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb`。
- `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx` 只读参考，不修改。
- `doc/worklog/v2506` 是历史依据，具体实现以参考代码实际源码为准。
- 若 OFv1706 API 与参考代码不兼容，按“最小可编译补丁”拆分移植，不跨阶段堆叠。
- 阶段 0 的完整性能剖析是后续所有优化排序和收益判断的强制前置条件。
