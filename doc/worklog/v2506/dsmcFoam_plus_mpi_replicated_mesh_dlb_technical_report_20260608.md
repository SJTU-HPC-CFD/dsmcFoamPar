# dsmcFoam+ MPI replicated mesh DLB 技术报告 - OFv1706 hyStrath_dlb

日期：2026-06-08

工作目录：

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

参考实现：

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
```

本文是本轮 `dsmcFoam+` 在 OFv1706 上的 MPI replicated mesh DLB 方向收尾
报告。报告范围是 `mpi8replicatedmesh` pure-MPI 路径，包括 replicated mesh
port、owner-filter、particle migration、ParMETIS AdaptiveRepart DLB、replicated
field/output、profile、DLB 触发策略，以及围绕 move/post 热路径的优化和回退。
OpenMP2 曾产生过较快结果，但用户已明确本轮正式口径必须是 pure MPI：
`useOpenMP false`、`openmpThreads 1`、`mpirun -np 8`，因此 OpenMP2 结果只作为
历史旁证，不作为本文正式结论。

## 1. 结论摘要

当前保留的技术状态：

- `dsmcReplicatedMesh` 已完成可运行 port，包含 owner map、local mesh、
  synchronous/async migration、ParMETIS AdaptiveRepart DLB、rank-reduced profile
  report、rank0 gather/write/re-migrate 输出路径。
- `dsmcCloud` 已接入 replicated mesh lifecycle：初始化、初始迁移、周期迁移、
  manual/auto DLB、output gather、profile report。
- `dsmcVolFields` 已补齐 replicated owned-cell / owned-boundary-face sampling
  和 output-rank reduce/write 语义，避免早期 replicated rank 重复采样导致的
  post/output 过高。
- accepted source 优化包括：
  - `boundaryMeasurements` sparse clean；
  - `particleTemplates.C` same-tet area reuse；
  - backup suffix hardening，避免 `lnInclude` 误链接 `particle/bkp` 中的备份
    `.C/.H` 文件。
- barycentric full tracking port 已证明可构建、可 smoke、可 500-step correctness，
  但最终没有留在 active source；当前主路径没有 active barycentric tracking 代码，
  只有 `.Cbkp/.Hbkp` 备份。
- 本轮最后一次用户要求后，DLB trigger-control runtime 改动已回退，当前 case
  为 area-reuse source + 非 forced legacy auto/SAR DLB 配置：

```text
useOpenMP false;
openmpThreads 1;
profileSummary true;
profileDetail false;
replicatedMeshDLBForceSteps ();
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBImbalanceThreshold 1.5;
```

注意：当前非 forced 配置是在最后回退后恢复的状态，尚未重新做 500-step retest。
因此本文的正式性能结论仍以已完成的 strict 500-step 日志为准，不能把当前
`ForceSteps ()` 状态直接等同于已有 forced8 性能。

关键性能节点：

| stage | source/control meaning | real | full evolve | move only | post fields/output | DLB |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| historical `ourmeshbkp/mpi8replicatedmesh` | preserved pure MPI reference | 117.31 s | 114.9207628 s | 47.1652636 s | 30.35563996 s | 8 |
| owner-filter early run | correctness fix, post still high | 184.13 s | n/a | 83.71 s | 75.75 s | 6 |
| first replicated DLB implementation | field/output and async path fixed | 152.07 s | n/a | 82.58246361 s | 43.16658299 s | 4 |
| pure-MPI scope-corrected forced8 | OpenMP2 removed from formal scope | 131.16 s | 129.1197774 s | 70.58583509 s | 41.44612826 s | 8 |
| strict same-tree baseline | retained forced8 repeat before post cleanup | 133.83 s | 132.7056254 s | 72.9567669 s | 40.37082178 s | 8 |
| sparse boundary clean repeat | accepted post cleanup | 113.24 s | 110.5389487 s | 70.66798424 s | 18.95231379 s | 8 |
| same-tet area reuse formal | accepted move/local geometry reuse | 104.32 s | 104.4499775 s | 67.43927759 s | 17.64348673 s | 8 |
| same-tet area reuse repeat | accepted repeat | 107.45 s | 105.2572367 s | 67.84144046 s | 17.2290821 s | 8 |
| area-reuse rebuild retest | after barycentric rollback/lnInclude repair | 108.69 s | 100.2458092 s | 62.51511822 s | 17.10465238 s | 8 |
| area-reuse rebuild repeat | same rebuilt state, slower repeat | 116.04 s | 108.3720354 s | 70.26858009 s | 18.14279997 s | 8 |
| current forced8 3-run mean | suffix-hardened area reuse, forced8 | 109.89 s | 101.86325 s | 64.5141319 s | 18.2834816 s | 8 |
| current auto/SAR 3-run mean | same source, `ForceSteps ()` | 115.816667 s | 108.660507 s | 70.9917698 s | 20.3781463 s | 4 |

结论：

1. 本轮 pure-MPI replicated mesh DLB 已从早期 correctness/port 阶段的
   `real 184.13 s` 降到 accepted area-reuse forced8 的 `real 104-110 s` 级别。
2. 最大稳定收益不是来自 ParMETIS 本身，而是来自 post/output 和 tracking hot
   path：
   - sparse boundary clean 将 post 从约 `40.37 s` 降到约 `19 s`；
   - same-tet area reuse 将 move 从约 `70.67 s` 降到 `67.44-67.84 s`，后续
     rebuild retest 还出现过 `62.52 s` 的 move-only 单次结果。
3. `particles per rank max/min` 不是单独的性能判据。DLB 使用当前
   `cellOccupancy()`、dual constraints 和 ParMETIS AdaptiveRepart；最终粒子
   max/min 变大不等价于总粒子数错误，也不必然等价于 rank-wall imbalance 变差。
4. 500-step 下 forced8 比旧 auto/SAR 更稳：3-run mean 中 forced8 比 auto 快
   `5.93 s real`、`6.80 s full evolve`。旧 auto 的 4 次 SAR 触发太稀疏。
5. `threshold=1.08` 的 threshold-only trigger 虽然在 500-step benchmark 中
   明显快于旧 auto/SAR，但这个阈值对上万步 DSMC 过小，已按用户要求回退到
   旧 auto/SAR 配置和旧源码触发逻辑。

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
platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+
```

最后一次 trigger rollback rebuild：

```text
doc/worklog/v2506/detail_mpi/build_dlb_trigger_revert_20260608.log
```

结果：构建通过，`libdsmcFoam+.so` 在 `2026-06-08 10:08` 重建。

### 2.2 正式算例

正式 MPI replicated mesh DLB case：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh
```

历史参考 case/log：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610
```

正式 pure-MPI 执行形式：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
cd run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.<name> 2>&1
```

正式 500-step controls：

```text
endTime 5.e-05;
deltaT 1.e-07;
profileSummary true;
profileDetail false;
useOpenMP false;
openmpThreads 1;
```

strict full-fields 配置要求保留一个 `dsmcVolFields` mixture entry。清空
`fieldPropertiesDict` 的 `dsmcFields ();` 是 no-fields/no-output 性能模式，
不能和 full-fields strict 结果混为一组。

### 2.3 控制文件备份规则

本轮所有 smoke/signal/formal 配置修改都先备份 `controlDict`，结束后恢复。
最后一次 trigger rollback 前的 forced8 备份为：

```text
run/.../mpi8replicatedmesh/system/controlDict.codex_forced_backup_before_trigger_revert_20260608_rollback
```

常用 strict forced8 hash：

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

常用 strict full-fields `fieldPropertiesDict` hash：

```text
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
```

## 3. Replicated Mesh DLB 实现结构

### 3.1 模块与 build 集成

核心源文件：

```text
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C
src/lagrangian/dsmc/replicatedMesh/dsmcLocalMesh.H
src/lagrangian/dsmc/replicatedMesh/dsmcLocalMesh.C
```

主要状态：

- `cellOwner_`：全局 cell owner map；
- `myCells_`：当前 rank 负责的 owned cells；
- `localMesh_`：owned + halo 的本地 mesh 视图；
- `allParticleCounts_`：迁移和 report 中使用的 rank 粒子数；
- `autoRebalanceChecks_` / `autoRebalanceCount_`：Phase C check/repartition 统计。

### 3.2 Phase A/B/C

Phase A：静态 owner decomposition 和 local mesh 构造。

Phase B：根据 `cellOwner_` 进行 particle migration。当前实现支持
synchronous flat transfer，也保留 async begin/finish API；`dsmcCloud` 已接入
`replicatedMeshDelayedReceive`，但正式 accepted path 仍需按实际控制项判断。

Phase C：ParMETIS AdaptiveRepart 动态负载均衡。

关键控制项：

```text
replicatedMeshAutoDLB true;
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshDLBDualConstraint true;
replicatedMeshDLBInitialAlpha 0.8;
replicatedMeshDLBUbvec 1.05;
replicatedMeshDLBUbvec1 1.5;
```

DLB 的 dual constraints：

- `N^alpha`：近似 move balance；
- `N*(N-1)`：近似 collision/candidate balance。

因此 DLB 目标不是最终 particle count 完全均匀，而是更接近 move/collision
rank-wall balance。

### 3.3 DLB cadence 的三个层次

报告中必须分开三件事：

1. `dsmcCloud::evolve()` 每步调用 `replicatedMesh_->autoRebalance()`；
2. `legacyWindow` 模式下 full collective check 受 `sarSteps` 等窗口控制；
3. 真正 ParMETIS repartition 还受 forced steps、`minGapSteps`、SAR/threshold
   条件控制。

一个日志里 `checks=500` 不代表 repartition 500 次；本轮 strict forced8
典型结果是 `checks=500`、`rebalances=8`。

### 3.4 Replicated output/profile

当前 solver 输出分支：

- `dsmc.evolve()` 后，replicated mesh 模式下把 parcels gather 到 output rank；
- 只有 output rank 调用 `runTime.write()`；
- 输出后再把 particles migrate 回 owner ranks 并更新 counts。

`dsmcVolFields` 当前已使用 owner cells / owned boundary faces 进行 replicated
sampling，再在 output-time 做 reduce/write。这个路径修复了早期 replicated rank
重复采样和 post/output 过高的问题。

## 4. Port 与正确性修复主线

### 4.1 初始 gap

早期 gap analysis 结论：

- replicated mesh module 已 port，不再是缺模块状态；
- 最大缺口不在 `dsmcReplicatedMesh` 类本身，而在 `dsmcVolFields`、
  output-rank gating、delayed receive wiring、rank-reduced profile 和 collision
  diagnostics。

早期 owner-filter 500-step run：

| metric | value |
| --- | ---: |
| real | 184.13 s |
| move only | 83.71 s |
| buildCellOccupancy | 6.78 s |
| collision phase | 14.25 s |
| post fields/output | 75.75 s |
| final particles | 2463897 |
| DLB rebalances | 6 |

owner-filter 修复了 replicated run 的粒子数爆炸问题：

```text
unfiltered replicated final particles = 8089972
owner-filtered final particles        = 2463897
mpi8origin current final particles    = 2463565
```

### 4.2 DLB implementation formal run

在补齐 `dsmcVolFields` owned sampling、output-rank reduce/write、
`replicatedMeshDelayedReceive` wiring、raw MPI profile reduction 后：

| metric | value |
| --- | ---: |
| real | 152.07 s |
| move only | 82.58246361 s |
| buildCellOccupancy | 9.800232796 s |
| collision phase | 15.99929084 s |
| post fields/output | 43.16658299 s |
| final particles | 2463839 |
| DLB rebalances | 4 |

相对 owner-filter 早期 run，post 从 `75.75 s` 降到 `43.17 s`，粒子数保持正确。
但它仍慢于历史 `117.31 s` 参考。

### 4.3 OpenMP2 结果的范围纠正

OpenMP2 replicated-mesh DLB 曾得到：

| metric | value |
| --- | ---: |
| real | 113.48 s |
| move only | 69.34624257 s |
| collision phase | 5.523175253 s |
| post fields/output | 33.6372609 s |

但用户明确本轮必须是 pure MPI replicated mesh DLB，OpenMP2 被排除在正式
比较之外。后续所有正式结论均以：

```text
useOpenMP false;
openmpThreads 1;
mpirun -np 8;
```

为边界。

## 5. Pure-MPI DLB 调参与触发策略

### 5.1 Pure-MPI scope correction

scope corrected pure-MPI forced8：

| metric | current pure MPI forced8 | historical pure MPI reference | gap |
| --- | ---: | ---: | ---: |
| real | 131.16 s | 117.31 s | +16.85 s |
| move only | 70.58583509 s | 47.1652636 s | +24.42057149 s |
| buildCellOccupancy | 7.037232476 s | 10.7084157 s | -3.671183224 s |
| collision phase | 14.06214838 s | 10.89824478 s | +3.16490360 s |
| post fields/output | 41.44612826 s | 30.35563996 s | +11.09048830 s |
| full evolve wall | 129.1197774 s | 114.9207628 s | +14.1990146 s |
| DLB rebalances | 8 | 8 | 0 |

### 5.2 Forced DLB steps

forced8 schedule：

```text
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

这个 schedule 在 500-step case 中更稳定地贴近历史 reference 的 8 次 DLB。
早期 DLB objective 调参结论：

| candidate | result |
| --- | --- |
| dense forced 13 steps | rejected；更多 rebalances 让 wall/move 变慢 |
| move-only ParMETIS constraint | rejected；move 和 collision 均变慢 |
| dual constraint alpha=1.0 | rejected；move/collision/post/end-to-end 均变慢 |
| no470 | 单次有改善，但后续正式口径仍回到 forced8 用于一致比较 |

### 5.3 Auto/SAR 与 forced8 对比

同一 source、strict 500-step、3 次 retest：

| metric | forced8 mean | auto/SAR mean | auto - forced |
| --- | ---: | ---: | ---: |
| real | 109.89 s | 115.816667 s | +5.926667 s |
| full evolve | 101.86325 s | 108.660507 s | +6.797257 s |
| move | 64.5141319 s | 70.9917698 s | +6.477638 s |
| post | 18.2834816 s | 20.3781463 s | +2.094665 s |
| migration wall | 2.24147209 s | 1.78783609 s | -0.453636 s |
| DLB rebalances | 8 | 4 | -4 |
| particles max/min | 2.54273483 | 1.88848103 | -0.654254 |
| rank wall max/min | 1.04319302 | 1.06797824 | +0.024785 |

解释：

- auto/SAR 用更少 DLB，migration wall 少约 `0.45 s`，但 move/build/collision/post
  增量远大于这点节省；
- auto/SAR 粒子 max/min 往往更好，但 rank-wall max/min 和 full evolve 更差；
- 对性能而言，rank wall 和 full evolve 比最终 particle max/min 更关键。

### 5.4 Threshold-only trigger rollback

临时测试过：

```text
replicatedMeshDLBForceSteps ();
replicatedMeshDLBMinStartStep 120;
replicatedMeshDLBRemainingGuardSteps 20;
replicatedMeshDLBUseSAR false;
replicatedMeshDLBImbalanceThreshold 1.08;
```

3 次均值：

| metric | optimized trigger mean | vs forced8 | vs old auto |
| --- | ---: | ---: | ---: |
| real | 109.42 s | -0.47 s | -6.396667 s |
| full evolve | 102.116092 s | +0.252842 s | -6.544415 s |
| move | 65.1392083 s | +0.625076 s | -5.852562 s |
| DLB rebalances | 6.333333 | -1.666667 | +2.333333 |

结论：

- 500-step 下它明显强于旧 auto/SAR；
- 但它未在 `full evolve` 和 `move` 上稳定超过 forced8；
- `threshold=1.08` 对上万步 DSMC 过小，容易被随机粒子占据波动频繁触发；
- 已按用户要求回退源码和 case 配置，回到旧 auto/SAR。

## 6. Post Fields/Output 优化

### 6.1 结构审计

`post fields/output [s]` 不是单纯磁盘输出；在当前 solver 中它主要覆盖
`dsmcCloud::evolve_fields()`：

- reactions output；
- fields calculate/write；
- controllers；
- boundaries；
- boundary measurements；
- `trackingInfo_.clean()`；
- `boundaryMeas_.clean()`；
- `cellMeas_.clean()`。

正式 500-step case：

```text
endTime 5.e-05;
writeInterval 1.e-3;
```

因此没有 `5e-05` 输出目录，formal `post fields/output` 中的大头不是 disk write。

### 6.2 No-fields mode

清空 `dsmcFields ();` 的 no-fields mode：

| metric | strict full fields | no fields | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 119.69 s | -14.14 s |
| move only | 72.9567669 s | 64.72452806 s | -8.23223884 s |
| post fields/output | 40.37082178 s | 30.60678571 s | -9.76403607 s |
| full evolve wall | 132.7056254 s | 118.6090831 s | -14.09654230 s |

它是合法的 no-field-measurement/no-field-output 性能模式，但不是 strict
full-fields 默认配置，不能作为 full-fields baseline。

### 6.3 Sparse boundary clean

诊断显示 `boundaryMeas_.clean()` 是可移除 post 子成本。accepted change：

- 为 boundary flux entries 维护 sparse touched tracking；
- wall boundary 写 boundary flux 时标记 `(species, patch, face)`；
- `boundaryMeasurements::clean()` 只清 touched entries；
- sticking boundary particle count full clear 保持不变。

500-step formal：

| metric | strict baseline | sparse clean | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 121.40 s | -12.43 s |
| move only | 72.9567669 s | 80.52274982 s | +7.56598292 s |
| post fields/output | 40.37082178 s | 19.56710696 s | -20.80371482 s |
| full evolve wall | 132.7056254 s | 121.1536038 s | -11.5520216 s |

500-step repeat：

| metric | strict baseline | sparse clean repeat | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 113.24 s | -20.59 s |
| move only | 72.9567669 s | 70.66798424 s | -2.28878266 s |
| post fields/output | 40.37082178 s | 18.95231379 s | -21.41850799 s |
| full evolve wall | 132.7056254 s | 110.5389487 s | -22.1666767 s |

决策：保留。它不禁用 field sampling，不改变 strict full-fields contract。

## 7. Move/Tracking 优化

### 7.1 same-tet non-normalised helper

早期 retained same-tet fast path 使用 unnormalised area vectors 快速判断
`endPosition` 是否仍在当前 tet。历史 formal pair：

```text
real 127.91-129.93 s
move only 67.56621687-69.11276280 s
full evolve 118.7495988-120.5815986 s
DLB rebalances 8
```

后续 same-day rebuild/retest 显示环境和 DLB partition 有明显波动：

```text
retained same-tet rerun: real 139.34-139.66 s
clean rebuild without same-tet diff: real 144.89-145.79 s
```

解释：same-tet source 仍比 clean source 快，但历史 `127-130 s` 与 same-day
`139 s` 不应混作同一基准。

### 7.2 same-tet area reuse

accepted change：

- same-tet fast path 失败时，第一轮 legacy tet-walk 仍在同一个 tet；
- 原代码会重算同一个 tet 的 4 个面积向量；
- 新实现复用 failed same-tet check 中已经计算的 `initialTetAreas[0..3]`；
- 后续 topology 变化后的 loop 仍重新计算 fresh areas；
- 不改变 `findTris()` / outer `tetLambda()` selection。

500-step strict formal：

| metric | strict baseline | sparse-clean repeat | area reuse formal | area reuse repeat |
| --- | ---: | ---: | ---: | ---: |
| real | 133.83 s | 113.24 s | 104.32 s | 107.45 s |
| move only | 72.9567669 s | 70.66798424 s | 67.43927759 s | 67.84144046 s |
| post fields/output | 40.37082178 s | 18.95231379 s | 17.64348673 s | 17.2290821 s |
| full evolve wall | 132.7056254 s | 110.5389487 s | 104.4499775 s | 105.2572367 s |
| DLB rebalances | 8 | 8 | 8 | 8 |

决策：保留。

### 7.3 Barycentric full port

本轮尝试过 reference-style barycentric tracking：

- v1706-compatible barycentric primitives；
- persistent barycentric state；
- static mesh wrapper；
- zero-progress fallback；
- explicit particle binary IO field ordering；
- processor patch transfer coverage；
- runtime switch `barycentricTracking true|false`。

最终 500-step replicated comparison中，barycentric hybrid correctness 通过，但未成为默认。
后来 active source 回退到 area-reuse Cartesian path；当前主路径没有 active
`barycentricTracking` / `trackToAndHitFaceBarycentric` 代码，只有 `.Cbkp/.Hbkp`
备份文件。

参考 formal 结果：

| replicated metric | barycentric on | final default Cartesian in that experiment | old strict baseline |
| --- | ---: | ---: | ---: |
| real | 118.32 s | 114.42 s | 114.9207628 s |
| move only | 69.68940261 s | 65.4460646 s | 47.1652636 s |
| full evolve wall | 108.2546436 s | 103.9147289 s | 114.9207628 s |
| DLB rebalances | 8 | 8 | 8 |

结论：barycentric 方向对 correctness 和未来 reference-style port 有价值，但它没有
关闭 move gap，且最终被回退，不作为当前 production default。

## 8. Rejected/Not-Retained 方案

### 8.1 DLB/objective 类

| candidate | reason |
| --- | --- |
| dense forced DLB steps | 13 次 rebalances 增加 overhead，real/move 变慢 |
| move-only ParMETIS constraint | particle balance 和 move/collision 都变差 |
| alpha=1.0 dual constraint | move/collision/post/end-to-end 均变慢 |
| old auto/SAR as formal default | 3-run mean 比 forced8 慢约 `5.93 s real` |
| threshold-only `1.08` trigger | 500-step 可用，但长步数 DSMC 阈值过低，已回退 |

### 8.2 Post/output 类

| candidate | reason |
| --- | --- |
| direct sampling | 局部 dsmcVolFields timers 下降，但总 post 没改善 |
| post flat serial path | targeted post sub-timers 变差 |
| dsmcN owned reset | 10-step smoke 未改善，回退 |
| final reset skip | 200-step signal 不够，未进 formal |
| sparse clean last-touch fast path | 200-step signal real/move/post 均未改善 |
| no-fields mode | 可作为 no-output 性能模式，但不是 strict full-fields default |

### 8.3 Move/tracking 类

| candidate | reason |
| --- | --- |
| pointer snapshot / cache-only move optimizations | formal pure-MPI move 未改善或 full run 变慢 |
| Pstream no-change fast path | added reduction/branch 不划算，smoke 变慢 |
| tracking API compatibility wrapper | 10-step 通过但 500-step move 大幅回退 |
| quad edge fast path / fixed tris array | 10-step tracking timer 变差 |
| same-face tet index increment | 单次 real 可能下降，但 move/full evolve 变差 |
| static lambda helper | 10-step 信号未在 500-step 保持 |
| direct same-tet normal calculation | 10-step 全面变慢 |
| outside-plane tri prefilter | 10-step move 好看，但 500-step formal 变差 |
| `tetNeighbour(0)` early return | formal real/move/full evolve 均回退 |
| same-tet sqrt-free tolerance | 10-step 全面回退 |
| direct Utracking | 500-step formal required metrics 全回退 |
| transition cache / trackingData cache / tetNeighbour cached refs | smoke/signal 不稳定，formal 不过 gate |

### 8.4 Migration/output 类

| candidate | reason |
| --- | --- |
| `replicatedMeshFlatTransfer false` | smoke 卡住到 `real 30882.88`，不可用 |
| `dsmcCloud` no-wall override | 500-step formal external wall/move/full evolve 回退 |

## 9. 正确性与验收口径

正式保留/对比的 500-step 日志均检查：

- `Total Iterations = 500`；
- `Number of stuck particles = 0`；
- final DSMC particles 在 `2.463e6` 量级；
- collision count、candidate count、total energy 与同类 run 的 DSMC 随机波动一致；
- 无 `Fatal`、`NaN`、`BAD TERMINATION`、`Segmentation`、`abort` 等标记；
- `OpenMP enabled = 0`、`OpenMP max threads = 1`；
- case `controlDict` 临时修改后恢复或明确记录当前状态。

性能 gate：

1. smoke 只用于排除明显错误或明显信号；
2. 200-step signal 可用于决定是否进入 formal；
3. retained source candidate 必须通过 strict full-fields 500-step；
4. 接受优化不能只看 `real`，还要看 `full evolve`、`move only`、post、DLB count、
   rank-wall balance；
5. particle max/min 只作为辅助解释，不单独决定保留或回退。

## 10. 当前状态

当前 active case：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/system/controlDict
```

当前关键配置：

```text
useOpenMP false;
openmpThreads 1;
endTime 5.e-05;
profileSummary true;
profileDetail false;
replicatedMeshDLBForceSteps ();
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBImbalanceThreshold 1.5;
```

当前 source 状态：

- active `particleTemplates.C` 含 same-tet area reuse；
- active boundary measurement path 含 sparse clean；
- active replicated mesh DLB trigger 逻辑已回到旧 `currentStep < 30` 和
  `sar_ > 0.0`；
- active source 不含 `replicatedMeshDLBMinStartStep` /
  `replicatedMeshDLBRemainingGuardSteps` / `replicatedMeshDLBUseSAR`；
- active source 不含 active barycentric tracking switch/path；
- `src/lagrangian/basic/particle/bkp` 备份文件已改为 `.Cbkp/.Hbkp`，避免
  `lnInclude` 误链接。

最后构建：

```text
doc/worklog/v2506/detail_mpi/build_dlb_trigger_revert_20260608.log
```

当前状态没有新的 500-step 性能结论。若要把当前非 forced auto/SAR 作为正式
baseline，必须重新做 strict 500-step retest。

## 11. 建议后续工作

1. 先重测当前状态：

```text
area-reuse source + non-forced legacy auto/SAR
useOpenMP false
openmpThreads 1
mpirun -np 8
strict full-fields 500 steps
```

2. 长步数 DSMC 的 DLB trigger 不应沿用 `threshold=1.08`。建议重新设计为
成本收益型：

```text
trigger if expected saved imbalance wall time > DLB cost * safety factor
```

工程近似可以从更保守参数开始：

```text
replicatedMeshDLBImbalanceThreshold 1.15-1.25
replicatedMeshDLBMinGapSteps 300-1000
minStart / remainingGuard 使用总步数比例，而不是 500-step 固定数
```

3. 若继续优化 move，优先围绕 active Cartesian same-tet/tet-walk path，而不是
重启已回退的 barycentric wrapper。

4. 若目标是生产长步数，不应只看 500-step benchmark。至少增加 2000/5000/10000
step 的 DLB cadence 和 partition stability 观测：

- actual DLB count；
- migration wall；
- rank-wall max/min；
- particle max/min；
- move/full-evolve trend by window；
- final particle/collision/energy consistency。

## 12. 主要工作日志索引

Core port / DLB：

```text
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_port_gap_analysis_20260606.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_dlb_implementation_20260606.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_scope_correction_pure_mpi_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_forced_dlb_steps_20260606.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_current_auto_vs_forced8_20260607.md
```

Post/output：

```text
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_post_output_structure_audit_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_no_fields_no_output_mode_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_post_boundary_sparse_clean_20260608.md
```

Move/tracking：

```text
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_tracking_same_tet_nonorm_20260606.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_tracking_tet_walk_followup_20260606.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_move_window_profile_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_barycentric_full_port_20260608.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_same_tet_area_reuse_20260608.md
```

Rejected candidate groups：

```text
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_move_candidate_rejections_20260606.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_move_trackingdata_cache_rejection_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_move_transition_cache_rejection_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_move_tetneigh_cached_refs_rejection_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_post_direct_sampling_rejection_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_post_flat_serial_rejection_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_post_final_reset_skip_rejection_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_post_dsmcN_owned_reset_rejection_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_flattransfer_false_rejection_20260607.md
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_dsmccloud_nowall_override_rejection_20260607.md
```

