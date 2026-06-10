# replicated-mesh processor 写出与 reconstructPar 验证记录（2026-06-10）

## 目标

在 raw-MPI replicated-mesh 模式下，不再只由 rank0 写 gathered cloud，而是输出 OpenFOAM 原生 `processorN/` case，使标准 `reconstructPar` 可以直接重建 lagrangian DSMC 数据。

## 代码状态

- 默认旧路径保持为 `replicatedMeshWriteMode gathered`。
- 新路径由 `replicatedMeshWriteMode processor` 显式启用。
- `dsmcFoam+` 在 processor 模式下跳过 rank0 gather/write/migrate-back，改为调用 `dsmcReplicatedMesh::writeProcessorOutput()`。
- `domainDecomposition` 支持外部 `cellToProc_`，并支持限制当前调用只写指定 `processorN`。
- `writeProcessorOutput()` 中所有 raw-MPI rank 都执行 mesh decomposition/write，但每个 rank 只写自己的 `processorN`，避免多个 rank 写同一目录。
- processor mesh 默认写到 `processorN/constant/polyMesh`。此前写到 `processorN/<time>/polyMesh` 会导致 `reconstructPar -noFields` 在时间 `0` 查找 `constant/polyMesh/points` 时失败。

## 本轮关键修正

失败现象：

```text
reconstructPar -noFields
FOAM FATAL ERROR:
Cannot find file "points" in directory "polyMesh" in times 0 down to constant
```

原因：

- smoke case 设置 `replicatedMeshProcessorWriteTimeMesh true`。
- processor mesh 被写入 `processorN/1e-07/polyMesh`、`processorN/2e-07/polyMesh`。
- 标准静态 processor case 需要 `processorN/constant/polyMesh`。

修正：

- `replicatedMeshProcessorWriteTimeMesh` 默认值改为 `false`。
- `writeProcessorOutput()` 显式把目标 processor mesh instance 传给 `domainDecomposition`。
- smoke case 改为 `replicatedMeshProcessorWriteTimeMesh false`。

## 验证 case

case：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi2replicatedmesh_processor_write_smoke_20260610
```

关键设置：

```text
numberOfSubdomains 2
endTime 2.e-07
deltaT 1.e-07
writeInterval 1.e-07
replicatedMesh true
replicatedMeshWriteMode processor
replicatedMeshProcessorWriteTimeMesh false
replicatedMeshAutoDLB false
```

## 命令与日志

构建：

```bash
source doc/scripts/build-dsmcFoam.sh
```

运行：

```bash
source doc/scripts/env.sh
timeout 300 mpirun -np 2 dsmcFoam+ -case "$CASE"
```

日志：

```text
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/processor_write_smoke_constant_mesh_rerun7_20260610_070906.log
```

结果：

```text
Exit status: 0
Elapsed wall clock: 0:13.58
processor write [s] at 1e-07 = 2.95
processor write [s] at 2e-07 = 2.74
```

## 产物检查

已生成：

```text
processor0/constant/polyMesh/{points,faces,owner,neighbour,boundary}
processor0/constant/polyMesh/{pointProcAddressing,faceProcAddressing,cellProcAddressing,boundaryProcAddressing}
processor1/constant/polyMesh/{points,faces,owner,neighbour,boundary}
processor1/constant/polyMesh/{pointProcAddressing,faceProcAddressing,cellProcAddressing,boundaryProcAddressing}
processor0/{1e-07,2e-07}/lagrangian/dsmc/{positions,U,ERot,vibLevel,typeId,newParcel,classification}
processor1/{1e-07,2e-07}/lagrangian/dsmc/{positions,U,ERot,vibLevel,typeId,newParcel,classification}
```

## reconstructPar 验证

`reconstructPar -noFields`：

```bash
timeout 120 reconstructPar -case "$CASE" -noFields
```

日志：

```text
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/reconstructPar_noFields_constant_mesh_rerun7_20260610_070941.log
```

结果：

```text
Exit status: 0
Reconstructing lagrangian fields for cloud dsmc
End
```

普通 `reconstructPar`：

```bash
timeout 120 reconstructPar -case "$CASE"
```

日志：

```text
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/reconstructPar_full_constant_mesh_rerun7_20260610_071001.log
```

结果：

```text
Exit status: 0
No FV fields
No point fields
Reconstructing lagrangian fields for cloud dsmc
End
```

## 粒子数一致性

`1e-07`：

```text
root positions       = 2058583
processor0 positions = 1195805
processor1 positions = 862778
sum                  = 2058583
```

`2e-07`：

```text
root positions       = 2059436
processor0 positions = 1196060
processor1 positions = 863376
sum                  = 2059436
```

## 当前结论

2-rank replicated-mesh smoke case 已经可以输出标准静态 `processorN/constant/polyMesh`、processor addressing 和 rank-local lagrangian cloud；OpenFOAM 原生 `reconstructPar -noFields` 与普通 `reconstructPar` 均可直接完成 lagrangian DSMC 重建。

## 4-rank 与 8-rank 扩展验证

按 2-rank smoke case 复制了两个独立输入 case，只复制 `0/`、`constant/`、`system/`、`boundaries/`、`fieldMeasurements/`，未复制已有 processor/time/log 输出：

```text
run/.../ourmesh/mpi4replicatedmesh_processor_write_smoke_20260610
run/.../ourmesh/mpi8replicatedmesh_processor_write_smoke_20260610
```

共同设置：

```text
endTime 2.e-07
deltaT 1.e-07
writeInterval 1.e-07
replicatedMeshWriteMode processor
replicatedMeshProcessorWriteTimeMesh false
replicatedMeshAutoDLB false
```

### 4-rank

命令：

```bash
timeout 600 mpirun -np 4 dsmcFoam+ -case "$CASE4"
timeout 180 reconstructPar -case "$CASE4" -noFields
timeout 180 reconstructPar -case "$CASE4"
```

日志：

```text
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/processor_write_smoke_constant_mesh_mpi4_20260610_083853.log
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/reconstructPar_noFields_constant_mesh_mpi4_20260610_083920.log
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/reconstructPar_full_constant_mesh_mpi4_20260610_083930.log
```

结果：

```text
dsmcFoam+ exit status = 0, elapsed = 0:12.29, max RSS = 1725276 KB
reconstructPar -noFields exit status = 0, elapsed = 0:10.13
reconstructPar exit status = 0, elapsed = 0:11.68
```

产物完整性：

```text
processor0..processor3 均有 constant/polyMesh/{points,faces,owner,neighbour,boundary}
processor0..processor3 均有 constant/polyMesh/{pointProcAddressing,faceProcAddressing,cellProcAddressing,boundaryProcAddressing}
processor0..processor3 的 1e-07/2e-07 均有 lagrangian/dsmc/{positions,U,ERot,vibLevel,typeId,newParcel,classification}
```

粒子数一致性：

```text
1e-07 root=2058615 processors=637593+559410+448164+413448 sum=2058615
2e-07 root=2059408 processors=637660+559572+448472+413704 sum=2059408
```

### 8-rank

命令：

```bash
timeout 900 mpirun -np 8 dsmcFoam+ -case "$CASE8"
timeout 240 reconstructPar -case "$CASE8" -noFields
timeout 240 reconstructPar -case "$CASE8"
```

日志：

```text
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/processor_write_smoke_constant_mesh_mpi8_20260610_083959.log
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/reconstructPar_noFields_constant_mesh_mpi8_20260610_084031.log
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/reconstructPar_full_constant_mesh_mpi8_20260610_084044.log
```

结果：

```text
dsmcFoam+ exit status = 0, elapsed = 0:19.36, max RSS = 1709936 KB
reconstructPar -noFields exit status = 0, elapsed = 0:12.85
reconstructPar exit status = 0, elapsed = 0:12.11
```

产物完整性：

```text
processor0..processor7 均有 constant/polyMesh/{points,faces,owner,neighbour,boundary}
processor0..processor7 均有 constant/polyMesh/{pointProcAddressing,faceProcAddressing,cellProcAddressing,boundaryProcAddressing}
processor0..processor7 的 1e-07/2e-07 均有 lagrangian/dsmc/{positions,U,ERot,vibLevel,typeId,newParcel,classification}
```

粒子数一致性：

```text
1e-07 root=2058580 processors=251174+200034+206879+204602+265451+293023+327772+309645 sum=2058580
2e-07 root=2059433 processors=251261+200241+207027+204742+265594+293071+327767+309730 sum=2059433
```

## 更新后结论

2-rank、4-rank、8-rank replicated-mesh smoke case 均已验证通过。当前 processor 写出路径可以生成标准静态 `processorN/constant/polyMesh`、processor addressing 与 rank-local lagrangian DSMC 数据；OpenFOAM 原生 `reconstructPar -noFields` 与普通 `reconstructPar` 在 2/4/8 rank 下均可直接完成重建，且重建后的 root `positions` 数量等于各 processor cloud 数量之和。

下一步可扩展到更长步数或真实 500-step benchmark，确认多输出时刻、大 I/O 量和 DLB 触发后的 `cellOwner_` 更新状态仍保持同样的 reconstructPar 兼容性。

## MPI4 OMP2 DLB 后写出验证

为确认 DLB 触发并更新 `cellOwner_` 后仍能写出可重建的 processor case，新增专用 smoke case：

```text
run/.../ourmesh/mpi4omp2_processor_write_dlb_smoke_20260610
```

该 case 从 4-rank processor-write smoke 输入目录复制而来，只复制 `0/`、`constant/`、`system/`、`boundaries/`、`fieldMeasurements/`，未复用旧输出。原 `mix-mpi4omp2` 性能 case 未修改。

关键设置：

```text
numberOfSubdomains 4
useOpenMP true
openmpThreads 2
endTime 3.e-06
writeInterval 3.e-06
replicatedMeshWriteMode processor
replicatedMeshProcessorWriteTimeMesh false
replicatedMeshAutoDLB true
replicatedMeshDLBForceSteps (30)
replicatedMeshDLBMinGapSteps 0
```

代码执行顺序确认：

```text
dsmc.evolve()
  -> replicatedMesh_->autoRebalance()
solver output block
  -> writeProcessorOutput()
```

`autoRebalance()` 内部有 `currentStep < 30` 保护，因此本验证强制第 30 步 DLB，并只在 `3e-06` 写出一次，保证该输出发生在 DLB 完成之后。

运行命令：

```bash
OMP_NUM_THREADS=2 timeout 900 mpirun -np 4 dsmcFoam+ -case "$CASE"
timeout 240 reconstructPar -case "$CASE" -noFields
timeout 240 reconstructPar -case "$CASE"
```

注意：第一次在 workspace sandbox 内启动 MPI 失败，日志为：

```text
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/processor_write_dlb_after_mpi4omp2_20260610_211937.log
```

失败原因是 sandbox 禁止 Hydra 打开本地通信 socket：

```text
cannot open socket (Operation not permitted)
```

随后用同一 case、同一命令在允许 MPI socket 的环境下重跑成功。

有效日志：

```text
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/processor_write_dlb_after_mpi4omp2_20260610_212011.log
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/reconstructPar_noFields_dlb_after_mpi4omp2_20260610_212107.log
doc/worklog/v2506/detail_mpi/processor_write_20260610/logs/reconstructPar_full_dlb_after_mpi4omp2_20260610_212114.log
```

DLB 证据：

```text
Phase C auto DLB triggered (forced) at step 30
Phase C ParMETIS AdaptiveRepart: 97767 / 104151 cells changed (0.9387043811)
Phase C auto DLB complete: rebalance #1
Replicated mesh: processor output begin at time 3e-06
```

运行结果：

```text
dsmcFoam+ exit status = 0, elapsed = 0:15.15, max RSS = 1725556 KB
OpenMP enabled = 1
OpenMP max threads = 2
Phase C auto DLB rebalances = 1
Total Iterations = 30
processor write [s] at 3e-06 = 1.85
```

reconstructPar 结果：

```text
reconstructPar -noFields exit status = 0, elapsed = 0:07.21
reconstructPar exit status = 0, elapsed = 0:07.10
```

产物完整性：

```text
processor0..processor3 均有 constant/polyMesh/{points,faces,owner,neighbour,boundary}
processor0..processor3 均有 constant/polyMesh/{pointProcAddressing,faceProcAddressing,cellProcAddressing,boundaryProcAddressing}
processor0..processor3 的 3e-06 均有 lagrangian/dsmc/{positions,U,ERot,vibLevel,typeId,newParcel,classification}
```

粒子数一致性：

```text
3e-06 root=2082048 processors=603430+525200+534124+419294 sum=2082048
```

### DLB 后写出结论

MPI4 OMP2 混合并行下，强制 DLB 后的 processor 写出已验证通过。该验证中 ParMETIS 改变了 97767 个 cell 的归属，随后 `writeProcessorOutput()` 使用更新后的 `cellOwner_` 生成 `processor0..processor3/constant/polyMesh` 与 rank-local lagrangian cloud；OpenFOAM 原生 `reconstructPar -noFields` 和普通 `reconstructPar` 均可重建，且重建后 root 粒子数与各 processor 粒子数之和一致。
