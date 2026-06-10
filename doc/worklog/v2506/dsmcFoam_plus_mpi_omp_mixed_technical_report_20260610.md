# dsmcFoam+ MPI+OpenMP mixed-parallel 技术报告 - OFv1706 hyStrath_dlb

日期：2026-06-10

工作目录：

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

参考实现：

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
```

本文是本轮 `dsmcFoam+` 在 OFv1706 上的 MPI+OpenMP mixed-parallel 方向收尾
报告。报告范围包括 replicated-mesh raw-MPI mixed path、OpenMP move/collision
与 replicated mesh lifecycle 的组合、`zb-cylinder-react` raw-MPI 修复、raw-MPI
输出路径修复、replicated-mesh ParMETIS DLB 权重污染修复，以及最终 no-write
8-core 对比测试。

本报告与以下报告互为同阶段技术文档：

```text
doc/worklog/v2506/dsmcFoam_plus_omp_technical_report_20260606.md
doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_dlb_technical_report_20260608.md
```

完整日志、controlDict 快照、CSV 解析结果和自动生成明细报告位于：

```text
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean
```

## 1. 结论摘要

当前保留的 mixed 技术状态：

- mixed MPI+OpenMP 使用 replicated-mesh raw-MPI 路径，不使用 OpenFOAM
  decomposed `-parallel` processor mesh；
- patch boundary 的 `controlParticle()` 在 OpenMP move 中进入保守 critical，
  避免 diffuse/specular wall 更新共享测量状态和共享 RNG 时并发崩溃；
- `Cloud::move()` 中 OpenFOAM processor-patch transfer 已对 replicated-mesh
  raw-MPI 模式加保护，避免 raw-MPI replicated mesh 错走标准 `-parallel` transfer；
- replicated-mesh initial distribution / migration 后会重建 cell occupancy，
  避免第一步 move/collision 使用 stale occupancy；
- raw-MPI replicated-mesh 输出路径使用 rank0 gather/write/re-migrate 方式，
  避免多个 MPI rank 并发写同一套 lagrangian cloud 文件；
- replicated-mesh ParMETIS DLB 使用全局每 cell 粒子数构造权重，修复 DLB 后
  本 rank `cellOccupancy()[globalCell]` 与 ParMETIS 连续 `vtxdist` 分片不一致
  导致的权重污染和 MPI8 卡死；
- profile/timing 能在 OMP、raw-MPI replicated mesh、mixed MPI+OMP 和标准
  `MPI8origin` 中统一提取外部 `real`、full evolve、move、build occupancy、
  collision、post fields/output、DLB 统计和正确性指标。

最终正式对比采用 no-write compute 口径：

- `ourmesh`：500 steps；
- `zb-cylinder-react`：300 steps；
- 每个 case 五组 8-core 模式；
- 每组重复三次；
- 共 30 次正式运行；
- 所有运行 `exit 0`，全部通过 correctness gates。

最终性能结论：

1. `OMP8` 仍是两个 case 的最佳 8-core 默认配置。
2. mixed 模式已经可用，但没有超过 `OMP8`：
   - `ourmesh` 最好 mixed 是 `MPI4xOMP2`，平均 `70.10 s`，比 `OMP8`
     慢 `12.66%`；`MPI2xOMP4` 为 `70.30 s`，差距只有 `0.20 s`；
   - `zb-cylinder-react` 最好 mixed 是 `MPI2xOMP4`，平均 `63.66 s`，
     比 `OMP8` 慢 `8.38%`。
3. pure replicated-mesh `MPI8` 已稳定跑完，不再出现 ParMETIS 卡死，但不适合
   作为默认 8-core 性能路径：
   - `ourmesh` 比 `OMP8` 慢 `63.04%`；
   - `zb-cylinder-react` 比 `OMP8` 慢 `35.86%`。
4. 标准 OpenFOAM decomposed `MPI8origin` 在两个 case 上都明显更慢：
   - `ourmesh` 平均 `135.07 s`，比 `OMP8` 慢 `117.07%`；
   - `zb-cylinder-react` 平均 `114.80 s`，比 `OMP8` 慢 `95.44%`。
5. mixed 模式的剩余 gap 主要来自 collision、build occupancy 和 MPI rank 间
   duplicated/owner-filtered work；no-write 口径下 post 已不是主导项。

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

主要可执行文件：

```text
platforms/linux64IccDPInt32Opt/bin/dsmcFoam+
platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+
```

### 2.2 正式 case 和五组模式

正式比较包括两个 case：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react
```

每个 case 的五组模式：

| mode | execution form | mesh/parallel meaning |
|---|---|---|
| `OMP8` | `OMP_NUM_THREADS=8 dsmcFoam+` | single process, 8 OpenMP threads |
| `MPI2xOMP4` | `OMP_NUM_THREADS=4 mpirun -np 2 dsmcFoam+` | replicated mesh raw-MPI, 2 ranks x 4 threads |
| `MPI4xOMP2` | `OMP_NUM_THREADS=2 mpirun -np 4 dsmcFoam+` | replicated mesh raw-MPI, 4 ranks x 2 threads |
| `MPI8` | `OMP_NUM_THREADS=1 mpirun -np 8 dsmcFoam+` | replicated mesh raw-MPI, 8 ranks |
| `MPI8origin` | `OMP_NUM_THREADS=1 mpirun -np 8 dsmcFoam+ -parallel` | standard OpenFOAM decomposed processor mesh |

`zb-cylinder-react/mpi8origin` 本轮没有现成 case，因此从 `zb-cylinder-react/omp8`
复制创建，修改 `decomposeParDict` 为 8 分区，并执行 `decomposePar` 生成：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8origin
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean/logs/zb_MPI8origin_decomposePar.log
```

### 2.3 no-write 控制

正式性能测试关闭场文件写出，但保留完整 stdout log：

```text
writeControl runTime;
writeInterval 1.e-3;
profileSummary true;
profileDetail false;
runTimeModifiable no;
```

`ourmesh`：

```text
endTime 5.e-05;
deltaT 1.e-07;
```

`zb-cylinder-react`：

```text
endTime 1.9920146682e-05;
deltaT 6.640048894e-08;
```

`writeInterval 1.e-3` 大于两个 case 的 `endTime`，因此正式步数范围内不触发
OpenFOAM write。已对正式 log 目录搜索写出标记，没有发现 `Writing` /
`runTime.write` / `writeTime` 等输出标记。

早期失败目录不参与统计：

```text
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610
```

该目录的失败原因是把 `writeControl` 改为 `timeStep` 后触发 `timeDataMeas`
除零/浮点异常。

另一个旧目录：

```text
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_postwrite
```

曾被中断运行污染，已经清除，不参与任何正式统计。

本报告只使用 ParMETIS 修复后的 clean 重测目录：

```text
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean
```

### 2.4 运行脚本和解析脚本

正式运行脚本：

```text
doc/worklog/v2506/detail_mix/run_mix_perf_nowrite_repeats_20260610.sh
doc/worklog/v2506/detail_mix/run_mpi8origin_nowrite_repeats_20260610.sh
doc/worklog/v2506/detail_mix/run_zb_mpi8origin_nowrite_repeats_20260610.sh
```

解析脚本：

```text
doc/worklog/v2506/detail_mix/parse_mix_perf_nowrite_repeats_20260610.py
```

输出结果：

```text
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean/run_manifest.tsv
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean/results.csv
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean/aggregate.csv
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean/performance_correctness_summary_20260610.md
```

## 3. 源码修改和技术实现

### 3.1 源码文件

本 mixed 阶段涉及的核心源码文件：

| file | purpose |
|---|---|
| `src/lagrangian/dsmc/parcels/dsmcParcel.C` | mixed OpenMP move 中保守保护 patch-boundary `controlParticle()` |
| `src/lagrangian/basic/Cloud/Cloud.C` | OpenMP move detail aggregation；replicated mesh raw-MPI 下跳过 OpenFOAM processor-patch transfer |
| `src/lagrangian/dsmc/clouds/dsmcCloud.C` | replicated-mesh initial distribution / migration 后重建 occupancy |
| `applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C` | replicated raw-MPI 输出时 gather cloud、rank0 write、re-migrate |
| `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C/H` | rank0 gathered cloud 写出 API、flat POD gather 路径、ParMETIS DLB 全局 cell 权重和权重饱和保护 |

当前工作树同时包含 mixed 运行脚本、ParMETIS 修复和若干后续优化线；本报告不再用
单一 diff 规模作为交付判据，正式结论以 clean 重测目录和定向 MPI8 验证日志为准。

```text
doc/worklog/v2506/detail_mix/mpi8_parmetis_fix_20260610
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean
```

### 3.2 构建证据

关键 build logs：

| stage | log | result |
|---|---|---|
| boundary critical fix | `detail_mix/stage_mix_boundary_critical_build_20260608.log` | build passed |
| move-detail aggregation | `detail_mix/stage_mix_move_detail_omp_aggregation_build_20260608.log` | build passed |
| raw-MPI move guard | `detail_mix/stage_zb_final_rawmpi_guard_build_20260608.log` | build passed |
| raw-MPI gathered cloud write | `detail_mix/stage_write_flat_gather_cloudprops_build_20260608.log` | build passed |

### 3.3 mixed OpenMP boundary crash 修复

最早的 `MPI4xOMP2` smoke 在第一步 OpenMP move 中崩溃，stack 指向
`dsmcPatchBoundary::measurePropertiesBeforeControl()`，调用链来自
`dsmcDiffuseWallPatch::controlParticle()`。

原因是 patch boundary control 会更新共享 wall measurement 状态，diffuse wall
路径还使用共享 cloud RNG。mixed OpenMP move 允许 diffuse/specular wall 并发后，
这些共享状态不再安全。

最终采用保守修复：所有 patch-boundary `controlParticle()` 调用进入既有
`dsmcMoveBoundary` critical section。修复后 10-step smoke 通过。

### 3.4 raw-MPI replicated mesh move 修复

`zb-cylinder-react` 最初的 raw-MPI replicated mesh 运行在第一步 move 中失败。
原因是 `Cloud::move()` 仍会进入 OpenFOAM standard `-parallel` processor-patch
transfer 逻辑，而 replicated mesh raw-MPI 模式没有标准 decomposed processor mesh。

修复方式是在 replicated mesh active 时跳过 OpenFOAM processor-patch transfer，
粒子跨 owner 的迁移由 replicated mesh owner map/migration 路径处理。

### 3.5 replicated initial distribution occupancy 修复

`zb-cylinder-react` 后续暴露 initial distribution 后 occupancy 不一致的问题。
修复方式是在 replicated-mesh initial distribution / migration 点之后重建
`cellOccupancy`，保证第一步 move/collision 使用当前 owner-rank 粒子状态。

### 3.6 raw-MPI 输出修复

raw-MPI replicated mesh 不使用 OpenFOAM `-parallel` processor directories。早期
输出路径在共享 cloud 文件写出时可能 hang 或失败。

最终输出顺序：

1. 所有 rank 把 parcels gather 到 output rank；
2. output rank 临时抑制 automatic cloud write；
3. output rank 调用 `runTime.write()` 写非 cloud 对象；
4. output rank 显式写 gathered lagrangian cloud 和 `cellOwner`；
5. 所有 rank 按 `cellOwner_` 把 parcels migrate 回 owner ranks。

这不是并行文件 I/O，而是 parallel communication + rank0 serial cloud write。
正式性能比较采用 no-write 口径，因此该路径主要作为功能正确性修复保留。

### 3.7 ParMETIS DLB 权重污染修复

MPI8 定向测试暴露 ParMETIS DLB 第二次以后可能不返回。根因不在 MPI8 本身，
而在 replicated-mesh DLB 的权重构造：DLB 后粒子实际归属由 `cellOwner_`
决定，本 rank 持有的 `cellOccupancy()` 不再等于 ParMETIS 连续 `vtxdist`
分片。如果继续用本 rank 的 `cellOccupancy()[globalCell]` 直接构造
`vwgt`、`adjwgt` 和 `vsize`，后续 AdaptiveRepart 会收到污染的权重输入。

修复方式：

- 在 `reassignByParMetisAdaptiveRepart()` 开始处对每个 cell 粒子数做
  `MPI_Allreduce`，得到全局每 cell 粒子数；
- `vwgt`、`adjwgt` 和 `vsize` 全部改用全局 cell 粒子数构造；
- 对权重转换加入正值保护和 `idx_t` 饱和保护，避免极端粒子数下整数溢出；
- 修正 `Ostream << long long` 的重载歧义，保证 OFv1706 编译通过。

定向验证：

| case | run | exit | DLB result | final marker |
|---|---|---:|---|---|
| `ourmesh` replicated `MPI8` | 320-step no-write | 0 | 2 次 DLB 返回，`Phase C ParMETIS max = 0.132687133 s` | `Total Iterations = 320`, `End main` |
| `ourmesh` replicated `MPI8` | 500-step no-write | 0 | 3 次 DLB 返回，`Phase C ParMETIS max = 0.18440913 s` | `Total Iterations = 500`, `End main` |

对应日志：

```text
doc/worklog/v2506/detail_mix/mpi8_parmetis_fix_20260610/ourmesh_MPI8_320step_parmetis_fix_20260610.log
doc/worklog/v2506/detail_mix/mpi8_parmetis_fix_20260610/ourmesh_MPI8_500step_parmetis_fix_20260610.log
```

该修复是 replicated-mesh + ParMETIS DLB 的通用修复。OMP-only 不走这条路径；
纯 MPI 和 mixed 只要启用 replicated-mesh auto DLB，都应使用修复后的全局权重。

## 4. Bring-up 和诊断结果

### 4.1 `ourmesh/mix-mpi4omp2` formal baseline

早期 mixed baseline：

```text
run/.../ourmesh/mix-mpi4omp2
OMP_NUM_THREADS=2 mpirun -np 4 dsmcFoam+
```

500-step 结果：

| metric | value |
|---|---:|
| real | 104.65 s |
| total iterations | 500 |
| particles | 2463810 |
| stuck particles | 0 |
| collisions | 35271 |
| candidates | 65578 |
| total energy | 1.241652653 |
| full evolve wall | 100.2682179 s |
| move only | 68.54304325 s |
| buildCellOccupancy | 7.287947673 s |
| collision phase | 4.94542299 s |
| post fields/output | 25.13833162 s |
| DLB rebalances | 3 |

该阶段证明 mixed path 可完成 500 steps，但后续 no-write repeat 显示当前最终
性能口径下 mixed move 已显著改善，剩余 gap 转向 collision/build/post。

### 4.2 move-detail diagnostic

200-step move-detail 诊断显示 boundary critical 不是主要瓶颈：

| metric | value |
|---|---:|
| real | 45.11 s |
| total iterations | 200 |
| move only | 27.47583734 s |
| buildCellOccupancy | 2.280479583 s |
| collision phase | 1.211962116 s |
| post fields/output | 8.968296275 s |
| move detail track max | 29.28668551 s |
| move detail boundary max | 0.100746499 s |

结论：boundary critical 修复正确但不是 mixed 主瓶颈；active hot path 仍是
same-tet/tet-walk tracking 和 replicated-mesh move data flow。

### 4.3 `zb-cylinder-react` 300-step bring-up

`zb-cylinder-react` 用于验证反应流和 raw-MPI replicated path。关键问题链：

1. first-move segfault：raw-MPI replicated mesh 错走 OpenFOAM processor-patch
   transfer；
2. initial distribution 后 occupancy stale；
3. output gather/write cleanup 失败或 hang；
4. no-write 正式性能需要保持 `writeControl runTime`，只把 `writeInterval`
   设成大于 `endTime`。

最终四个 replicated/raw-MPI 模式和新增 `MPI8origin` 标准 decomposed baseline
均完成 300 steps。

## 5. 最终 no-write 对比测试

### 5.1 完整性和正确性

最终 manifest 包含 30 次运行：

```text
2 cases x 5 modes x 3 repeats = 30 runs
```

所有运行满足：

- exit code `0`；
- expected iteration count；
- `End main` present；
- no fatal / segmentation fault / MPI abort / killed / NaN markers；
- final particle count nonzero；
- stuck particles `0`；
- finite total energy。

### 5.2 `ourmesh` 500-step 性能

| mode | ok/runs | real mean | stdev | min-max | full evolve | move | buildOcc | collision | post | delta vs OMP8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 3/3 | 62.22 | 1.60 | 61.10-64.06 | 60.46 | 49.11 | 6.62 | 3.52 | 0.65 | 0.00% |
| MPI2xOMP4 | 3/3 | 70.30 | 1.48 | 68.60-71.18 | 68.12 | 54.41 | 7.54 | 16.87 | 1.07 | +12.99% |
| MPI4xOMP2 | 3/3 | 70.10 | 1.49 | 68.42-71.24 | 67.91 | 51.95 | 6.84 | 12.67 | 2.20 | +12.66% |
| MPI8 | 3/3 | 101.45 | 1.43 | 100.01-102.87 | 98.63 | 69.07 | 10.61 | 32.80 | 5.00 | +63.04% |
| MPI8origin | 3/3 | 135.07 | 2.31 | 132.44-136.78 | 144.27 | 117.47 | 16.31 | 20.56 | 0.48 | +117.07% |

`ourmesh` 结论：

- `OMP8` 最快；
- `MPI4xOMP2` 是平均最快 mixed split，但只比 `MPI2xOMP4` 快 `0.20 s`；
- mixed move 已低于 pure replicated `MPI8`，但仍慢于 `OMP8`；
- replicated `MPI8` 已稳定，但 move/build/collision 都高于 mixed split；
- standard decomposed `MPI8origin` 的 move 极高，是最慢组。

### 5.3 `zb-cylinder-react` 300-step 性能

| mode | ok/runs | real mean | stdev | min-max | full evolve | move | buildOcc | collision | post | delta vs OMP8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 3/3 | 58.74 | 1.25 | 57.31-59.61 | 57.76 | 41.69 | 4.26 | 10.46 | 0.18 | 0.00% |
| MPI2xOMP4 | 3/3 | 63.66 | 1.22 | 62.49-64.92 | 63.59 | 42.89 | 5.11 | 13.43 | 0.27 | +8.38% |
| MPI4xOMP2 | 3/3 | 67.72 | 0.17 | 67.58-67.91 | 67.42 | 44.66 | 3.73 | 24.39 | 0.55 | +15.28% |
| MPI8 | 3/3 | 79.81 | 1.59 | 78.20-81.38 | 81.33 | 54.01 | 6.51 | 33.73 | 1.29 | +35.86% |
| MPI8origin | 3/3 | 114.80 | 1.20 | 113.42-115.60 | 124.21 | 90.73 | 11.23 | 28.69 | 0.13 | +95.44% |

`zb-cylinder-react` 结论：

- `OMP8` 最快；
- `MPI2xOMP4` 是平均最快 mixed split，比 `OMP8` 慢 `8.38%`；
- `MPI4xOMP2` 的 build occupancy 更低，但 collision 明显高于 `MPI2xOMP4`；
- replicated `MPI8` 已稳定，但 move/collision 高于 mixed split；
- standard decomposed `MPI8origin` 比 replicated `MPI8` 还慢，主要因为 move 和
  build occupancy 高。

### 5.4 Correctness metric ranges

| case | mode | particles min-max | collisions mean | total energy min-max | stuck max |
|---|---|---:|---:|---:|---:|
| ourmesh | OMP8 | 2463545-2463681 | 35591 | 1.2428346060-1.2430942860 | 0 |
| ourmesh | MPI2xOMP4 | 2463659-2463786 | 35550 | 1.2430114520-1.2431917970 | 0 |
| ourmesh | MPI4xOMP2 | 2463788-2463938 | 35614 | 1.2423908420-1.2428044570 | 0 |
| ourmesh | MPI8 | 2463776-2463829 | 35101 | 1.2417871210-1.2424040990 | 0 |
| ourmesh | MPI8origin | 2463640-2463765 | 35517 | 1.2432013240-1.2433580940 | 0 |
| zb | OMP8 | 1957484-1958390 | 241823 | 0.0019609871-0.0019620457 | 0 |
| zb | MPI2xOMP4 | 1958061-1958320 | 226886 | 0.0019510445-0.0019613185 | 0 |
| zb | MPI4xOMP2 | 1959318-1960011 | 222185 | 0.0019434562-0.0019474198 | 0 |
| zb | MPI8 | 1957998-1958771 | 206645 | 0.0019326233-0.0019393615 | 0 |
| zb | MPI8origin | 1957683-1958417 | 240788 | 0.0019622125-0.0019627143 | 0 |

粒子数和总能量在各模式间保持可接受一致性。碰撞数随并行分解和 RNG 消耗顺序
变化，这是 DSMC 并行调度下预期差异；本轮没有任何 correctness gate 失败。

### 5.5 DLB diagnostics

Replicated-mesh raw-MPI 模式的 DLB 统计：

| case | mode | migration wall mean | particles max/min mean | rank wall max/min mean | rebalances mean |
|---|---|---:|---:|---:|---:|
| ourmesh | MPI2xOMP4 | 1.05 | 1.1623 | 1.0005 | 2.67 |
| ourmesh | MPI4xOMP2 | 1.53 | 1.2190 | 1.0018 | 2.33 |
| ourmesh | MPI8 | 2.96 | 1.3181 | 1.0030 | 2.67 |
| zb | MPI2xOMP4 | 2.06 | 1.2048 | 1.0008 | 2.33 |
| zb | MPI4xOMP2 | 2.47 | 1.2532 | 1.0013 | 2.00 |
| zb | MPI8 | 2.39 | 1.3424 | 1.0020 | 2.33 |

解释：

- `particles max/min` 不是单独性能判据；
- 本轮 rank-wall ratio 通常只有 `1.0005-1.0030`，说明 DLB 后 critical-path balance
  并不差；
- pure `MPI8` 仍然慢，主要不是 DLB cadence 本身，而是 collision/move/build
  代价在更多 MPI ranks 下上升。

## 6. 与 OMP/MPI 单项报告的关系

### 6.1 相对 OMP 技术报告

OMP 报告的主结论是 stage25 OMP8 已经超过历史 `ourmeshbkp/omp8`，并且 OMP
路径的核心收益来自 collision 并行化、move 热路径、post fields/output 和
`buildCellOccupancy()`。

本 mixed 报告的最终对比延续这个结论：在 8-core 范围内，`OMP8` 仍然是最强
baseline。mixed 要超过 OMP8，不能只依赖 MPI 分解；必须进一步降低 mixed
collision 和 build occupancy 开销。

### 6.2 相对 MPI replicated mesh DLB 技术报告

MPI replicated mesh DLB 报告证明 pure-MPI replicated mesh 已从 correctness/port
阶段进入可运行性能状态，核心收益来自 post/output cleanup 和 same-tet area reuse。

本 mixed 报告显示：mixed split 可以显著优于 pure replicated `MPI8`，并且
ParMETIS 全局 cell 权重修复后 pure replicated `MPI8` 已稳定通过 500-step 和
30 组 clean repeat。replicated mesh + OMP 的组合是有效的，但 pure `MPI8`
在 8-core 下仍不是合理默认性能模式。

### 6.3 相对 standard OpenFOAM decomposed MPI8

本轮补齐了两个 case 的 `MPI8origin` baseline。结论一致：

| case | MPI8origin vs OMP8 | MPI8origin vs replicated MPI8 |
|---|---:|---:|
| ourmesh | +117.07% real | +33.14% real |
| zb | +95.44% real | +43.85% real |

标准 decomposed MPI8 的主要问题是 move/build 开销显著高于 OMP8 和 replicated
MPI8。因此除非必须保持 standard OpenFOAM `-parallel` 工作流，否则不应把
`MPI8origin` 作为当前 8-core 性能默认路径。

## 7. 剩余优化方向

优先级按当前 no-write repeat 数据排序：

1. mixed collision 路径。`ourmesh` 中 collision 从 `OMP8 3.52 s` 增至
   `MPI2xOMP4 16.87 s`、`MPI4xOMP2 12.67 s`、`MPI8 32.80 s`；`zb`
   中 collision 从 `OMP8 10.46 s` 增至 `MPI2xOMP4 13.43 s`、
   `MPI4xOMP2 24.39 s`、`MPI8 33.73 s`。应检查 collision
   是否仍有 per-rank 全网格/空 cell 遍历、owner filtering、RNG 或反应调用开销。
2. `buildCellOccupancy()`。`ourmesh MPI8 10.61 s`、`ourmesh MPI8origin 16.31 s`、
   `zb MPI8 6.51 s`、`zb MPI8origin 11.23 s`，说明更多 MPI ranks 下该路径
   仍有可优化空间。可评估 owned/active
   cell list 或迁移后增量更新。
3. OMP8 move。若目标是继续提升当前最佳配置，`OMP8` 的 move 仍是最大项：
   `ourmesh 49.11 s`、`zb 41.69 s`。mixed 是否能超过 OMP8，最终仍取决于
   move/collision shared baseline 是否继续下降。
4. `zb-cylinder-react` 反应路径。反应 case 的 collision/reaction 成本对 split
   很敏感，应继续检查反应调用、owner filtering 和 rank-local 统计开销。
5. 运行时绑核和波动控制。clean repeat 中 `zb MPI8` 的 `real` 范围已收敛到
   `78.20-81.38 s`，但 mixed split 的相对排序仍值得用 `OMP_PROC_BIND`、
   `OMP_PLACES`、Intel MPI pinning 做低成本 A/B。
6. DLB threshold/cadence 不是下一步第一优先级。本轮 rank-wall ratio 已很低，
   继续调 ParMETIS/DLB 不太可能单独改变端到端排序。

## 8. 最终建议

当前 8-core 正式默认配置仍应使用 `OMP8`。

如果需要 mixed MPI+OpenMP：

- `ourmesh`：`MPI4xOMP2` 平均最快，`MPI2xOMP4` 只慢 `0.20 s`，两者都可作为
  mixed 候选；
- `zb-cylinder-react`：优先用 `MPI2xOMP4`。

不要把 pure replicated `MPI8` 或 standard decomposed `MPI8origin` 作为当前
8-core 默认性能路径。它们可以作为 MPI path 功能验证和对照基准保留，但不是
性能最佳选择。ParMETIS 修复后 replicated `MPI8` 可以作为稳定的 DLB 功能验证
路径保留。
