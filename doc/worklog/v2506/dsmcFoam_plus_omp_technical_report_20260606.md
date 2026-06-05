# dsmcFoam+-OMP 技术报告 - OFv1706 hyStrath_dlb

日期：2026-06-06

工作目录：

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

参考实现：

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
```

本文是本轮 `dsmcFoam+` 在 OFv1706 上的 OMP 方向收尾报告。报告范围是
单进程 8 线程 OpenMP 路径，包括 profile 基础设施、OpenMP move、
OpenMP collision、post fields/output、`buildCellOccupancy()` 以及对应的
正确性/性能验证。MPI replicated mesh、MPI DLB、MPI+OMP 混合并行只在必要
的接口兼容处提到，不作为本文主线结论。

## 1. 结论摘要

最终保留状态是 stage25：

- `Cloud::move()` 保留 stage21 的无删除 fast path；
- `dsmcVolFields` 保留 stage20 的 shared sample cache 和 boundary patch
  mask；
- `buildCellOccupancy()` 保留 stage25 的 sparse reset、`DynamicList`
  复用和 appended-aware ordered path；
- OpenMP collision 保留 stage5 引入的 per-thread NTC buffers 和 fast RNG
  路径；
- profile/timer/控制项保留 stage1 后续逐步完善的实现。

最终 `ourmesh/omp8` 正式 500-step OMP 结果：

| metric | historical `ourmeshbkp/omp8` | final stage25 | delta | pct |
|---|---:|---:|---:|---:|
| real [s] | 103.66 | 90.34 | -13.32 | -12.849701% |
| full evolve wall [s] | 102.1225033 | 90.00671836 | -12.11578494 | -11.863972% |
| move only [s] | 48.82365912 | 49.24385297 | +0.42019385 | +0.860636% |
| buildCellOccupancy [s] | 6.490896398 | 6.447003770 | -0.043892628 | -0.676218% |
| collision phase [s] | 9.876225069 | 3.376808391 | -6.499416678 | -65.808714% |
| post fields/output [s] | 36.92670834 | 30.39917036 | -6.527537980 | -17.677010% |

最终 stage25 相对 stage21 稳定性均值：

| metric | stage21 stability mean | final stage25 | delta | pct |
|---|---:|---:|---:|---:|
| real [s] | 98.38333333 | 90.34 | -8.04333333 | -8.175504% |
| full evolve wall [s] | 98.16965099 | 90.00671836 | -8.16293263 | -8.315129% |
| move only [s] | 53.73855140 | 49.24385297 | -4.49469843 | -8.364011% |
| buildCellOccupancy [s] | 8.408399703 | 6.447003770 | -1.961395933 | -23.326626% |
| collision phase [s] | 3.608682147 | 3.376808391 | -0.231873756 | -6.425441% |
| post fields/output [s] | 31.83707414 | 30.39917036 | -1.437903780 | -4.516444% |

最终结论：

1. 相对历史 `ourmeshbkp/omp8` 基准，最终 OMP 结果已经超过历史基准，
   端到端 `real` 降低约 12.85%。
2. 相对本轮 OFv1706 最早可靠 profile 阶段 stage4，最终 `full evolve
   wall` 从 272.2832057 s 降到 90.00671836 s，降低约 66.94%。
3. 主要收益来源分阶段不同：
   - stage5 主要来自 collision 路径并行化；
   - stage7/stage8/stage9/stage15 主要来自 move 热路径；
   - stage17/stage20 主要来自 post fields/output；
   - stage25 主要来自 `buildCellOccupancy()`；
   - stage21 的 no-delete fast path 是 stage20 后端到端跨过历史基准的
     关键一步。
4. 正确性方面，最终 `ourmesh/omp8` stage25 500-step run 完成 500/500
   iterations，无 stuck particles，无 fatal/NaN/segfault/MPI abort 标记；
   粒子数、碰撞数、候选数、能量接近 stage21 稳定性带，偏差处于很小的
   DSMC 随机波动量级。
5. `zb-cylinder-react/omp8` 300-step 验证也通过，但它暴露了一个剩余问题：
   v1706 stage25 的 `buildCellOccupancy = 7.902615635 s` 仍显著慢于
   v2506 reference 的 `3.491106571 s`。这个 gap 不应视为已关闭。

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

### 2.2 正式算例

正式 OMP 500-step 算例：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8
```

历史基准与配置恢复来源：

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/omp8
```

正式 OMP run 使用 `ourmeshbkp/omp8/system/controlDict` 中的 OMP/profile 配置：

```text
useOpenMP true;
openmpThreads 8;
openmpMoveSchedule static;
openmpMoveChunk 64;
openmpCollisionSchedule dynamic;
openmpCollisionChunk 8;
profileSummary true;
profileDetail false;
collisionFastRng true;
```

每次正式测试的执行规则：

1. 不重新初始化，直接运行；
2. 临时把 `ourmeshbkp/omp8/system/controlDict` 复制到
   `ourmesh/omp8/system/controlDict`；
3. 运行完成后恢复原 `ourmesh/omp8/system/controlDict`；
4. 校验恢复 hash：

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d
```

正式运行命令形式：

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
cd run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8
OMP_NUM_THREADS=8 /usr/bin/time -p dsmcFoam+ > log.<stage-name> 2>&1
```

### 2.3 辅助验证算例

辅助反应流 OMP 300-step 算例：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/omp8
```

该 case 原始 `controlDict` 有 v2506 风格：

```text
writeControl none;
```

OFv1706 不接受该枚举。本轮只做兼容性单行修改：

```diff
-writeControl none;
+writeControl runTime;
```

备份：

```text
system/controlDict.codex_bak_20260606_before_v1706_writecontrol
```

## 3. 起点分析和历史基准

### 3.1 源码起点

stage0 审计确认：当前 OFv1706 源码不是历史 `ourmeshbkp` profiling log 的
等价源码。它缺少大量 v2506/reference 能力：

- `Cloud.C` 缺少 reference OpenMP move path、ordered parcel extraction、
  move profiling helpers；
- `dsmcCloud.C` 缺少 profile controls、OpenMP move/collision controls、
  flat/moveOrdered occupancy、replicated mesh lifecycle、DLB reporting；
- `dsmcVolFields.C/H` 与 reference 有大量 post-field/replicated-output 差异；
- replicated mesh 源码目录最初不存在；
- METIS/ParMETIS 链接环境最初也不存在。

因此，历史 `ourmeshbkp` 日志只能作为性能参考，不能直接假设当前源码已经
实现同一代码路径。本轮 OMP 工作必须分阶段移植、编译、smoke、formal。

### 3.2 历史 OMP 基准

`ourmeshbkp/omp8/log.omp8.confirm_pdFalse_20260605_024610`：

| metric | value |
|---|---:|
| external wall [s] | 103.66 |
| full evolve wall [s] | 102.1225033 |
| move only [s] | 48.82365912 |
| buildCellOccupancy [s] | 6.490896398 |
| collision phase [s] | 9.876225069 |
| post fields/output [s] | 36.92670834 |
| final particles | 2463391 |
| final collisions | 35847 |
| collision candidates | 75838 |

时间占比：

| phase | wall [s] | share |
|---|---:|---:|
| move only | 48.82365912 | 47.81% |
| buildCellOccupancy | 6.490896398 | 6.36% |
| collision phase | 9.876225069 | 9.67% |
| post fields/output | 36.92670834 | 36.16% |

基准解释：

- OMP8 是 move-heavy + post-heavy；
- collision 只占约 10%，所以只做 collision 并行不能解决端到端；
- `buildCellOccupancy` 占比不大，但在后续 OFv1706 移植阶段多次成为
  剩余 gap；
- post fields/output 是第二大耗时，必须单独优化。

## 4. Profile 和 OpenMP 控制基础设施

stage1 首先移植低扰动 profile/timer 设施，而不是直接上大规模 OMP 改动。

新增/保留能力：

- `profileSummary`；
- `profileDetail`；
- `full evolve wall`；
- `move+collide wall`；
- `move only`；
- `buildCellOccupancy`；
- `collision phase`；
- `post fields/output`；
- CPU timing 汇总；
- stage end 的 `dsmc.printProfileSummary()`；
- 并行时对阶段时间做 rank max reduction。

关键实现细节：

- wall timer 后续改为 `std::chrono::steady_clock`，避免 OFv1706
  `clockTime/gettimeofday()` 粒度和非单调问题；
- `buildCellOccupancy()` 只统计 evolve 内部调用，构造期 initial occupancy
  不纳入 solver-step summary；
- `profileDetail false` 是正式比较默认配置，避免诊断 profile 本身污染
  生产路径。

## 5. 阶段演进总表

下表只列正式 OMP 500-step 或与正式 500-step 可比的关键阶段。stage3 的
phase timer 当时尚不可靠，因此只在后文说明，不放入同口径阶段表。

| stage | real [s] | full evolve wall [s] | move [s] | build [s] | collision [s] | post [s] | 说明 |
|---|---:|---:|---:|---:|---:|---:|---|
| historical bkp | 103.66 | 102.1225033 | 48.82365912 | 6.490896398 | 9.876225069 | 36.92670834 | 历史参考 |
| stage4 | 265.13 | 272.2832057 | 144.8259732 | 7.895838180 | 68.18518398 | 50.78234158 | steady timer + flat post，首个可靠 phase profile |
| stage5 | 163.13 | 167.3012128 | 104.1127577 | 7.730683037 | 3.632583082 | 51.26087370 | retained move/collision，collision 大幅下降 |
| stage6 | 161.79 | 167.0823348 | 103.6495219 | 7.483010573 | 3.613524884 | 51.76373504 | full move outer-control port，小幅收益 |
| stage7 | 151.06 | 156.8234055 | 95.58241050 | 7.194255854 | 3.516419771 | 49.98115219 | tracking hot-path cache，move 明显下降 |
| stage8 | 137.46 | 141.1479884 | 78.62462603 | 7.449130794 | 3.631726566 | 50.86713931 | OMP move complete，move 主收益 |
| stage9 | 127.02 | 128.4992481 | 64.61371202 | 7.589475879 | 3.711549217 | 51.98935525 | inline reset + stationary-tet fast path |
| stage15 | 114.53 | 116.6882995 | 52.99979159 | 8.663262117 | 3.628983484 | 50.80901580 | uniform-dt fast path |
| stage17 | 106.84 | 110.2731630 | 53.31126355 | 8.160171479 | 3.692066583 | 44.52651830 | post local accumulation |
| stage20 | 103.34 | 102.6520492 | 56.13548084 | 9.115554564 | 3.760305061 | 33.03582209 | shared sample cache + boundary mask |
| stage21 best | 100.65 | 93.28696662 | 51.17943521 | 8.087915571 | 3.314541395 | 30.15687469 | no-delete move fast path，首次明显超历史 |
| stage21 mean | 98.38333333 | 98.16965099 | 53.73855140 | 8.408399703 | 3.608682147 | 31.83707414 | 三次稳定性均值 |
| stage25 final | 90.34 | 90.00671836 | 49.24385297 | 6.447003770 | 3.376808391 | 30.39917036 | final kept OMP result |

stage25 相对 stage4：

| metric | stage4 | stage25 | delta | pct |
|---|---:|---:|---:|---:|
| real [s] | 265.13 | 90.34 | -174.79 | -65.926149% |
| full evolve wall [s] | 272.2832057 | 90.00671836 | -182.2764873 | -66.943713% |
| move only [s] | 144.8259732 | 49.24385297 | -95.58212023 | -65.997913% |
| buildCellOccupancy [s] | 7.895838180 | 6.447003770 | -1.448834410 | -18.349343% |
| collision phase [s] | 68.18518398 | 3.376808391 | -64.80837559 | -95.047592% |
| post fields/output [s] | 50.78234158 | 30.39917036 | -20.38317122 | -40.138305% |

## 6. OpenMP move 技术路线

### 6.1 stage5：reference-retained move 基础

stage5 首先移植 reference 中仍保留且能适配 OFv1706 的 move hot-path 改动：

- `dsmcCloud::trackerActive_`；
- `refreshTrackerUsage()`；
- 无 tracker 时 gate 掉 face-transition tracker 调用；
- wall/boundary path 使用 patch-to-model 直接映射，减少扫描；
- 简单 wall patch model 走直接 OMP path；
- 非线程安全 patch model 仍保留 critical 保护；
- OpenMP move kernel 使用 per-thread move RNG。

stage5 formal：

| metric | stage4 | stage5 | delta |
|---|---:|---:|---:|
| real [s] | 265.13 | 163.13 | -102.00 |
| full evolve wall [s] | 272.2832057 | 167.3012128 | -104.9819929 |
| move only [s] | 144.8259732 | 104.1127577 | -40.7132155 |
| collision phase [s] | 68.18518398 | 3.632583082 | -64.552600898 |

stage5 的端到端收益很大，但真正最大的是 collision；move 仍远慢于历史。

### 6.2 stage6：完整 move outer-control port

stage6 将 `Cloud<ParticleType>::move()` 改为 reference-style
extract/kernel/commit 路径，并适配 OFv1706 的 `PstreamBuffers`、
`IDLList`、`finishedSends` 机制：

- pre-move append capture；
- `moveOrderedParcels()` 和 appended parcels 复用；
- per-thread survivor/delete 分类；
- 非 MPI 和 MPI+OMP 路径均能启动；
- MPI transfer 分支仍保持旧 OFv1706 通信机制。

stage6 相对 stage5：

| metric | stage5 | stage6 | delta |
|---|---:|---:|---:|
| full evolve wall [s] | 167.3012128 | 167.0823348 | -0.2188780 |
| move only [s] | 104.1127577 | 103.6495219 | -0.4632358 |
| buildCellOccupancy [s] | 7.730683037 | 7.483010573 | -0.247672464 |

结论：outer-control port 正确但收益很小，剩余 move gap 在更底层 tracking
核心、tet 路径或版本差异中。

### 6.3 stage7：tracking hot path cache

stage7 在不复制 v2506 新 API 的前提下优化 OFv1706 tracking 路径：

- 缓存 cloud、mesh、boundary mapping、coordSystem type、tracker gate；
- boundary face path 只计算一次 patch index；
- `trackToFace(..., true)` 缓存 owner/neighbour/tet-base/cell-volume 等引用；
- 当 `hasWallImpactDistance()` 为 false 时跳过 `hitWallFaces()`；
- 保持 DSMC 专用路径，不改普通 `trackToFace()` 模板。

stage7 相对 stage6：

| metric | stage6 | stage7 | delta |
|---|---:|---:|---:|
| full evolve wall [s] | 167.0823348 | 156.8234055 | -10.2589293 |
| move only [s] | 103.6495219 | 95.58241050 | -8.0671114 |
| post fields/output [s] | 51.76373504 | 49.98115219 | -1.78258285 |

结论：这是 stage5 后第一个对 formal move 有明显收益的 move 级改动。

### 6.4 stage8：OMP move complete

stage8 完成 retained move-OMP port：

- append capture 提前到 pre-move controls/boundaries 前；
- `Cloud::move()` 保持幂等 begin；
- OMP move path 不再限制于 `!Pstream::parRun()`；
- MPI+OMP 下 OMP kernel 分类 keep/switch/delete，然后交给旧 v1706 transfer。

stage8 相对 stage7：

| metric | stage7 | stage8 | delta |
|---|---:|---:|---:|
| full evolve wall [s] | 156.8234055 | 141.1479884 | -15.6754171 |
| move only [s] | 95.58241050 | 78.62462603 | -16.95778447 |
| buildCellOccupancy [s] | 7.194255854 | 7.449130794 | +0.25487494 |

结论：stage8 的收益几乎完全来自 move，副作用是 build/collision/post 有
小幅波动但远小于 move 收益。

### 6.5 stage9：inline reset + stationary-tet fast path

stage9 两个关键点：

- `stepFraction = 0` 从串行前扫移动到 OMP kernel 内；
- DSMC stationary-mesh `trackToFace(endPosition, td, true)` fast path：
  如果目标点仍在当前 tet 且无 wall-impact-distance 处理，直接更新
  `position_` 并返回。

stage9 相对 stage8：

| metric | stage8 | stage9 | delta |
|---|---:|---:|---:|
| full evolve wall [s] | 141.1479884 | 128.4992481 | -12.6487403 |
| move only [s] | 78.62462603 | 64.61371202 | -14.01091401 |
| buildCellOccupancy [s] | 7.449130794 | 7.589475879 | +0.140345085 |
| post fields/output [s] | 50.86713931 | 51.98935525 | +1.12221594 |

结论：收益仍然 move-dominated，说明 particle tracking 仍是当时最大瓶颈。

### 6.6 stage15：uniform-dt fast path

stage15 引入时间步模型能力：

- `dsmcTimeStepModel::uniformDeltaT()`；
- `dsmcConstantTimeStepModel` 返回 true；
- `dsmcVariableTimeStepModel` 返回 false；
- `dsmcCloud` 缓存 `uniformDeltaT_`；
- constant time-step case 中 parcel move 直接使用已有 `trackTime`，
  避免内层反复调用 `cloud.deltaTValue(cell)`。

stage15 相对 stage12：

| metric | stage12 | stage15 | delta |
|---|---:|---:|---:|
| real [s] | 116.55 | 114.53 | -2.02 |
| full evolve wall [s] | 118.4147219 | 116.6882995 | -1.7264224 |
| move only [s] | 54.62342851 | 52.99979159 | -1.62363692 |

stage16 后续缓存 branch/reduced-D 尝试在 smoke 变慢，已回退。

### 6.7 stage21：no-delete fast path

stage21 针对 pure OMP 非 MPI move 常见无删除情况：

- OMP move kernel 内统计 `keepParticle == false`；
- 如果删除数为 0，直接复用 ordered parcel path；
- 跳过后续 survivor/delete 计数扫描和临时 survivor/delete list；
- MPI transfer branch 不变；
- 如果有删除，仍回到原完整 commit path。

stage21 相对 stage20：

| metric | stage20 | stage21 best | delta |
|---|---:|---:|---:|
| real [s] | 103.34 | 100.65 | -2.69 |
| full evolve wall [s] | 102.6520492 | 93.28696662 | -9.36508258 |
| move only [s] | 56.13548084 | 51.17943521 | -4.95604563 |
| buildCellOccupancy [s] | 9.115554564 | 8.087915571 | -1.027638993 |
| post fields/output [s] | 33.03582209 | 30.15687469 | -2.87894740 |

stage21 是当时保留的最优配置，并通过三次 500-step 稳定性复测确认：

| metric | stage21 stability mean | range | CV |
|---|---:|---:|---:|
| real [s] | 98.38333333 | 96.01-99.84 | 2.107090% |
| full evolve wall [s] | 98.16965099 | 96.50163862-99.04519771 | 1.472082% |
| move only [s] | 53.73855140 | 52.42366711-54.42676369 | 2.119780% |
| buildCellOccupancy [s] | 8.408399703 | 8.283163669-8.529617478 | 1.466106% |
| post fields/output [s] | 31.83707414 | 31.55022843-32.00615842 | 0.784424% |

## 7. OpenMP collision 技术路线

collision 的核心收益主要来自 stage5：

- `noTimeCounter` 使用 per-thread subcell buffers；
- per-thread parcel pointer/type/velocity/charge buffers；
- OpenMP cell-level collision loop；
- per-thread fast RNG；
- flat occupancy 可用时复用；
- serial collision fallback 保留。

stage4 到 stage5：

| metric | stage4 | stage5 | delta |
|---|---:|---:|---:|
| collision phase [s] | 68.18518398 | 3.632583082 | -64.552600898 |

之后 collision phase 基本保持在 3.3-3.8 s：

| stage | collision phase [s] |
|---|---:|
| stage5 | 3.632583082 |
| stage8 | 3.631726566 |
| stage15 | 3.628983484 |
| stage20 | 3.760305061 |
| stage21 best | 3.314541395 |
| stage25 final | 3.376808391 |

结论：

- collision 不是最终 OMP 瓶颈；
- stage5 后继续深入 collision 的收益空间有限；
- stage25 相对历史基准的总优势中，collision 仍贡献很大，但这是 stage5
  及后续保留路径的结果，不是 stage25 occupancy patch 的直接效果。

## 8. post fields/output 技术路线

### 8.1 stage4：flat occupancy sampling

stage4 先引入轻量 OMP sampling：

- 如果 flat occupancy ordered parcels 可用，则按 cell 并行采样；
- 每个 cell 由单线程写入；
- 不引入 atomics；
- `openmpFieldSampling` 默认 true。

同时修复 profile wall timer，保证正式 profile 内部一致。

stage4 证明：

- flat field sampling 相对关闭版本可降低 post；
- 但整体仍很慢，正式 500-step `real = 265.13 s`，远慢于历史。

### 8.2 stage17：local accumulation

stage17 对 `dsmcVolFields::calculateField()` 做局部累加优化：

- 缓存 `openmpFieldSampling`；
- 构建 `typeId -> species index` map；
- flat occupancy view 上按 cell 并行；
- 跳过 empty cells；
- 使用局部 scalar/vector 累加后一次写回；
- `measureHeatFluxShearStress false` 时跳过高阶热流/剪切累加。

stage17 相对 stage15：

| metric | stage15 | stage17 | delta |
|---|---:|---:|---:|
| real [s] | 114.53 | 106.84 | -7.69 |
| full evolve wall [s] | 116.6882995 | 110.2731630 | -6.4151365 |
| post fields/output [s] | 50.80901580 | 44.52651830 | -6.2824975 |

结论：stage17 是 post 主线第一个正式有效阶段。

### 8.3 stage20：shared sample cache + boundary patch mask

stage20 保留 reference-style shared sample cache，进一步解决 boundary
accumulation 过多的问题：

- 保留 shared sample cache 和 active-cell reset/combine；
- 将 `cloud.nParticles(celli)` 移出 per-parcel OMP cache build loop；
- 新增 `sampledBoundaryPatches_`；
- 跳过 `processor`、`empty`、`symmetry`、`wedge` patches；
- 在 `ourmesh/omp8` 中避免 `BaseAndTop` symmetry patch 208302 faces 的
  无效 boundary measurement accumulation；
- 缓存 `cloud_.boundaryFluxMeasurements()` 引用。

stage20 相对 stage17：

| metric | stage17 | stage20 | delta |
|---|---:|---:|---:|
| real [s] | 106.84 | 103.34 | -3.50 |
| full evolve wall [s] | 110.2731630 | 102.6520492 | -7.6211138 |
| post fields/output [s] | 44.52651830 | 33.03582209 | -11.49069621 |

stage20 相对历史基准：

| metric | historical bkp | stage20 | delta |
|---|---:|---:|---:|
| post fields/output [s] | 36.92670834 | 33.03582209 | -3.89088625 |

结论：post fields/output 在 stage20 后已经快于历史 OMP 基准，不再是
正式主要 gap。

## 9. buildCellOccupancy 技术路线

### 9.1 早期失败尝试

stage18 尝试：

- member `occupancyTotalCounts_` 复用；
- per-thread cell counts active-cell sparse clear；
- 避免 `nextCellOffsets(occupancyCellOffsets_)` 全量复制。

10-step smoke：

| metric | stage17 kept smoke | stage18 smoke | delta |
|---|---:|---:|---:|
| buildCellOccupancy [s] | 0.106092180 | 0.113830229 | +0.007738049 |
| real [s] | 6.72 | 6.89 | +0.17 |

结论：stage18 拒绝并回退。

stage22 在 stage21 后再次尝试 active-cell reset：

| metric | stage21 kept | stage22 candidate | delta |
|---|---:|---:|---:|
| full evolve wall [s] | 93.28696662 | 98.63267366 | +5.34570704 |
| buildCellOccupancy [s] | 8.087915571 | 8.651435857 | +0.563520286 |

结论：stage22 拒绝并回退。

### 9.2 stage25 最终 occupancy 改动

stage25 不是简单重复 stage18/stage22，而是结合了三个条件：

1. `occupancyOrderedParcels_` 从 `List<dsmcParcel*>` 改为
   `DynamicList<dsmcParcel*>`，用 `resize()` 复用容量；
2. ordered path 支持 `moveOrderedParcels_ + moveAppendedParcels_`，适配
   inflow/reaction 追加粒子；
3. 只有在 ordered path 命中时，对 per-thread `localCounts` 做上轮 active
   cells sparse reset；fallback gather 仍保持全 clear。

stage25 对 `ourmesh/omp8`：

| metric | stage21 stability mean | stage25 final | delta | pct |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 98.16965099 | 90.00671836 | -8.16293263 | -8.315129% |
| buildCellOccupancy [s] | 8.408399703 | 6.447003770 | -1.961395933 | -23.326626% |

stage25 对 `zb-cylinder-react/omp8`：

| metric | stage21 | stage25 | delta | pct |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 75.24167335 | 71.98270041 | -3.25897294 | -4.331340% |
| buildCellOccupancy [s] | 8.394238035 | 7.902615635 | -0.491622400 | -5.856665% |

结论：

- stage25 是保留的 occupancy 优化；
- 对 `ourmesh/omp8`，它将 occupancy 拉回历史参考水平附近；
- 对 `zb-cylinder-react/omp8`，它有收益但仍远慢于 v2506 reference；
- 下一个 occupancy 工作应先做 subphase profile，而不是继续盲改 clear/fill。

## 10. 拒绝项和回退原因

| stage | candidate | 结果 | 回退原因 |
|---|---|---|---|
| stage16 | cached reduced-D tracking flag / dt branch | rejected | smoke 中 `move only` 和 `full evolve wall` 均变慢 |
| stage18 | early buildCellOccupancy sparse clear | rejected | 10-step `buildCellOccupancy` 从 0.1061 s 变 0.1138 s |
| stage22 | buildCellOccupancy active-cell reset | rejected | formal `buildCellOccupancy` 和 `full evolve wall` 都变慢 |
| stage23 | move delete-index recording | rejected | external real 变快但 internal profiled path 变慢，formal 不可信 |
| stage24 | uninitialized keep/switch flags | rejected | move 小幅改善但 real 和 build 回归，收益太小 |

技术经验：

1. smoke positive 不能替代 500-step formal；
2. external `real` 不能单独作为保留依据，必须同时看 internal
   `full evolve wall` 和目标 phase；
3. active-cell sparse reset 只有和正确的 ordered/appended 判定一起做才有价值；
4. move flag 写入类微优化很容易被内存布局、cache、后续 phase 波动抵消；
5. `profileDetail false` 必须保持生产路径干净，不能把诊断计数放回 move
   热路径。

## 11. 正确性验证总结

正式保留阶段均检查：

- `Total Iterations = 500`；
- `End main`；
- stuck particles 为 0；
- 无 `FOAM FATAL`；
- 无 `NaN` / `nan`；
- 无 segmentation marker；
- 无 `MPI_ABORT`；
- 粒子数、碰撞数、候选数、总能量处于同量级随机波动范围内，并结合
  stuck particles 和错误扫描判断正确性。

最终 `ourmesh/omp8` stage25：

| metric | value |
|---|---:|
| iterations | 500 |
| real [s] | 90.34 |
| full evolve wall [s] | 90.00671836 |
| particles | 2463726 |
| stuck particles | 0 |
| collisions | 35977 |
| collision candidates | 69118 |
| total energy | 1.243107069 |

stage21 三次稳定性复测给出的随机波动带：

| metric | mean | min | max |
|---|---:|---:|---:|
| particles | 2463658.67 | 2463620 | 2463727 |
| collisions | 35717.33 | 35678 | 35748 |
| collision candidates | 68942.00 | 68636 | 69204 |
| total energy | 1.243132700 | 1.243124401 | 1.243145619 |

stage25 与这个稳定性带接近；碰撞数略高，但仍在 DSMC 随机统计可接受范围，
没有数值异常或 stuck-particle 问题。

## 12. 性能收益来源归因

### 12.1 相对历史 `ourmeshbkp/omp8`

stage25 总收益：

| component | historical bkp | stage25 | delta | share of full-wall gain |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 102.1225033 | 90.00671836 | -12.11578494 | 100.000% |
| move only [s] | 48.82365912 | 49.24385297 | +0.420193850 | -3.468% |
| buildCellOccupancy [s] | 6.490896398 | 6.447003770 | -0.043892628 | 0.362% |
| collision phase [s] | 9.876225069 | 3.376808391 | -6.499416678 | 53.644% |
| post fields/output [s] | 36.92670834 | 30.39917036 | -6.527537980 | 53.876% |

解释：

- 相对历史基准，主要收益来自 stage5 collision 和 stage20/21 后保留的
  post-field 路径；
- stage25 occupancy 的作用是把 build 与历史参考持平，而不是贡献大部分
  总收益；
- move 与历史参考基本持平，仍略慢 0.42 s。

### 12.2 相对 stage21 最优

stage25 对 stage21 best：

| component | stage21 best | stage25 | delta | share of full-wall gain |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 93.28696662 | 90.00671836 | -3.28024826 | 100.000% |
| move only [s] | 51.17943521 | 49.24385297 | -1.935582240 | 59.007% |
| buildCellOccupancy [s] | 8.087915571 | 6.447003770 | -1.640911801 | 50.024% |
| collision phase [s] | 3.314541395 | 3.376808391 | +0.062266996 | -1.898% |
| post fields/output [s] | 30.15687469 | 30.39917036 | +0.242295670 | -7.387% |

解释：

- 直接可归因到 stage25 代码的是 `buildCellOccupancy` 的下降；
- `move only` 同时下降，但 stage25 没有直接改 move kernel，因此应视为
  run-to-run 波动或内存/cache 副作用，除非后续复跑确认；
- collision/post 略慢，抵消一部分收益。

### 12.3 `zb-cylinder-react/omp8` 对照

stage25 对当前 v1706 stage21：

| component | stage21 | stage25 | delta |
|---|---:|---:|---:|
| full evolve wall [s] | 75.24167335 | 71.98270041 | -3.25897294 |
| move only [s] | 41.62369781 | 40.34492359 | -1.27877422 |
| buildCellOccupancy [s] | 8.394238035 | 7.902615635 | -0.491622400 |
| collision phase [s] | 10.26836592 | 9.837595081 | -0.430770839 |
| post fields/output [s] | 13.75492210 | 12.73351032 | -1.021411780 |

stage25 对 v2506 reference：

| component | v2506 reference | stage25 v1706 | delta |
|---|---:|---:|---:|
| full evolve wall [s] | 74.71342017 | 71.98270041 | -2.73071976 |
| move only [s] | 38.10127114 | 40.34492359 | +2.24365245 |
| buildCellOccupancy [s] | 3.491106571 | 7.902615635 | +4.411509064 |
| collision phase [s] | 20.61886309 | 9.837595081 | -10.781268009 |
| post fields/output [s] | 12.49893663 | 12.73351032 | +0.23457369 |

解释：

- `zb` 的 v1706 stage25 总 wall 与 v2506 reference 接近甚至略快，但来源
  不是 occupancy；
- v1706 的 `buildCellOccupancy` 仍比 v2506 reference 慢 126%；
- 总 wall 能接近，是 collision phase 大幅更低抵消了 build/move/post 的
  劣势；
- 因此 `zb` 后续如果继续优化，应优先做 `buildCellOccupancy` subphase
  profile。

## 13. 当前保留代码结构

主要保留模块：

### `src/lagrangian/basic/Cloud/Cloud.C`

保留内容：

- OMP extract/kernel/commit move path；
- ordered parcel reuse；
- append capture 配合；
- no-delete fast path；
- MPI transfer 分支兼容 OFv1706 `PstreamBuffers`；
- 非 OpenMP path 保持原兼容逻辑。

### `src/lagrangian/basic/particle/particleTemplates.C`

保留内容：

- DSMC stationary tet fast path；
- OpenMP 安全的局部 `DynamicList<label> tris(4)`；
- mesh/tet/boundary 引用缓存；
- `hasWallImpactDistance()` gate。

### `src/lagrangian/dsmc/clouds/dsmcCloud.H/C`

保留内容：

- profile/timer counters；
- OpenMP controlDict controls；
- OpenMP move/collision schedule/chunk 控制；
- move ordered/appended parcel data；
- flat occupancy view；
- stage25 `buildCellOccupancy` sparse/dynamic/appended path；
- profile summary 输出。

### `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter`

保留内容：

- per-thread buffers；
- OpenMP cell-level collision loop；
- fast RNG；
- serial fallback。

### `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields`

保留内容：

- flat occupancy sampling；
- local accumulation；
- shared sample cache；
- active-cell reset/combine；
- sampled boundary patch mask；
- species map/cache；
- `openmpFieldSampling` 控制。

### time-step model headers

保留内容：

- `uniformDeltaT()` API；
- constant/variable model 区分；
- move loop constant-dt fast path。

## 14. 当前遗留问题

### 14.1 `zb-cylinder-react/omp8` occupancy 仍慢于 v2506

当前数据：

```text
v2506 reference buildCellOccupancy = 3.491106571 s
v1706 stage25 buildCellOccupancy   = 7.902615635 s
```

需要下一步补充：

- `buildCellOccupancy` subphase profile；
- ordered path 命中率；
- fallback gather 次数；
- appended parcel 数量；
- per-thread active cell 数；
- reset/count/reduce/fill 分项；
- `occupancyOrderedParcels_` resize/capacity 行为。

### 14.2 move 与历史参考基本持平但未稳定领先

最终 stage25：

```text
historical move = 48.82365912 s
stage25 move    = 49.24385297 s
```

move 已从 stage4 的 144.826 s 大幅下降，但最终仍略慢于历史参考。进一步
优化需要 profileDetail 级别的 move subphase，而不是继续盲改 tracking
分支。

### 14.3 stage25 后还缺少多次稳定性复跑

stage21 做过 3 次稳定性复跑。stage25 目前有一次有效 `ourmesh/omp8`
formal run。虽然结果正确且明显更快，但若要作为长期发布基线，建议再做
2-3 次 500-step formal 复跑，确认：

- `full evolve wall` 是否稳定在 90-92 s；
- `buildCellOccupancy` 是否稳定在 6.4-6.6 s；
- move 下降是否可重复；
- 粒子数/碰撞数/能量是否仍在 stage21 稳定性带附近。

### 14.4 controlDict 配置风险

`ourmesh/omp8/system/controlDict` 当前恢复 hash `ff112...`，该文件本身不含
OMP/profile entries。正式 OMP 回归必须临时使用 `ourmeshbkp/omp8` 的
`9e7d...` OMP/profile controlDict。否则会产生无效 no-OMP run：

```text
OpenMP enabled = 0
real = 365.19 s
full evolve wall = 385.4477822 s
```

该日志已在 stage25 记录中标为 invalid no-OMP，不可用于 OMP 性能结论。

### 14.5 OpenMP schedule 策略测试状态

截至 stage25 正式 `ourmesh/omp8` 500-step，生产配置统一使用：

```text
openmpMoveSchedule static;
openmpMoveChunk 64;
openmpCollisionSchedule dynamic;
openmpCollisionChunk 8;
```

该组合来自历史 v2506/reference 工作中的生产推荐，并在 OFv1706 stage5
到 stage25 的正式 500-step 过程中持续作为比较基准使用。

stage26 已在 `zb-cylinder-react/omp8` 上补做当前 OFv1706 代码的
schedule/chunk sweep。该 sweep 是两个独立 300-step 矩阵，而不是四维
全组合：

- move sweep：固定 collision 为 `dynamic / 8`，扫描 move
  `static/dynamic x chunk 8,16,32,64,128,256`；
- collision sweep：固定 move 为 `static / 64`，扫描 collision
  `static/dynamic x chunk 8,16,32,64,128,256`。

stage26 关键结论：

- move：`static / 256` 在 `zb` 上最好；相对同组 `static / 64`，
  full wall `78.06102765 -> 77.11135257 s`，`-1.22%`，
  move-only `43.85303096 -> 42.88355957 s`，`-2.21%`；
- move：`dynamic` 对单 rank OMP8 没有收益，所有 chunk 的 move-only
  都慢于 `static / 64`；
- collision：`dynamic / 8` 仍是最优，collision phase
  `10.35458841 s`；
- collision：`static / 8` 的 collision phase 已比 `dynamic / 8`
  慢 `11.75%`，`static / 256` 慢 `35.12%`；
- collision：`dynamic / 64` 的 full wall 与 `dynamic / 8` 接近
  (`+0.10%`)，但 collision phase 已慢 `3.07%`，不建议放大为默认
  chunk。

stage26 边界：

- 该结论来自 `zb-cylinder-react/omp8` 300-step single-run sweep；
- `ourmesh/omp8` 正式 500-step 尚未复测 `move static / 256`；
- 因此正式 `ourmesh/omp8` 结果仍以 `move static / 64 +
  collision dynamic / 8` 为已验证配置；若要改 move chunk，需要在
  `ourmesh/omp8` 500-step 上单独验证。

stage26 前，本轮 OFv1706 中出现过但不构成完整 sweep 的 schedule/chunk
变化：

- stage8/stage9 的 MPI+OMP smoke 用过 `openmpMoveSchedule static`、
  `openmpMoveChunk 1`、`openmpCollisionChunk 1`，只用于启动/通信正确性；
- stage6 有一个非正式 500-step run 使用 `openmpCollisionChunk 1`，
  记录为不可与 stage5/formal bkp controlDict 直接比较；
- stage25 有一次 no-OMP 误配置 run，`OpenMP enabled = 0`，不可用于
  schedule 结论。

旧 v2506 工作日志中的 schedule 结论：

- OMP8 单 rank：`move` 改 `dynamic/guided` 没有本质收益；主要问题不是
  粒子数调度不均，而是单位粒子 tracking 成本不同。因此继续在
  `move schedule/chunk` 层反复调参性价比低；
- OMP8 collision：`dynamic` 是更稳妥的选择；candidate-only partition
  或 weighted partition 没有稳定超过 dynamic；
- collision chunk 扫描做过 `chunk=4/8/16/64`，结论是 chunk 调优只带来
  轻微波动，不能改变主趋势；
- 2MPI x 4OMP / 4MPI x 2OMP 等混合配置中，`move dynamic chunk=256`
  曾有小幅收益，但这是混合并行配置，不应直接外推到单 rank OMP8；
- guided schedule 在某些混合配置下与 guard-cell/固定线程范围假设冲突，
  曾出现 segfault 风险。

因此，当前 OFv1706 报告中的 OMP8 schedule 结论应表述为：

```text
ourmesh/omp8 正式 500-step 已沿用并验证 v2506/reference 推荐配置：
move static chunk 64 + collision dynamic chunk 8。
zb-cylinder-react/omp8 stage26 300-step sweep 进一步支持 collision dynamic
chunk 8；move chunk 在 zb 上显示 static chunk 256 更好，但尚未完成
ourmesh/omp8 正式 500-step 复测。
```

如果后续继续补 schedule 证据，优先级应是：

| test | move schedule/chunk | collision schedule/chunk | 用途 |
|---|---|---|---|
| ourmesh move retest | static / 256 | dynamic / 8 | 验证 zb 上的 move chunk 收益能否迁移到正式 500-step |
| ourmesh baseline rerun | static / 64 | dynamic / 8 | 与 retest 同批比较，控制运行时漂移 |
| collision guard | static / 64 | dynamic / 4,16 | 只有怀疑 chunk=8 不稳时再补 |
| guided smoke only | guided / 256 | dynamic / 8 | 只做 smoke，先查安全性 |

## 15. 报告证据索引

主要阶段文档：

| file | 用途 |
|---|---|
| `doc/worklog/v2506/detail/stage0_baseline_source_audit_20260605.md` | 起点源码差距和 reference 审计 |
| `doc/worklog/v2506/detail/stage0_full_performance_profile_20260605.md` | historical benchmark profile |
| `doc/worklog/v2506/detail/stage1_omp_base_port_profile_timers_20260605.md` | profile/timer 基础设施 |
| `doc/worklog/v2506/detail/stage3_formal_ourmesh_500step_results_20260605.md` | 早期 formal 结果和初始 OMP 回归 |
| `doc/worklog/v2506/detail/stage4_omp_flat_post_steadywall_results_20260605.md` | steady wall timer 和 flat post |
| `doc/worklog/v2506/detail/stage5_omp_ref_retained_move_collision_results_20260605.md` | retained OMP move/collision |
| `doc/worklog/v2506/detail/stage6_omp_full_move_port_results_20260605.md` | full move outer-control port |
| `doc/worklog/v2506/detail/stage7_move_tracking_results_20260605.md` | tracking hot path |
| `doc/worklog/v2506/detail/stage8_move_omp_complete_results_20260605.md` | complete OMP move |
| `doc/worklog/v2506/detail/stage9_move_inline_reset_fastpath_results_20260605.md` | inline reset / stationary-tet fast path |
| `doc/worklog/v2506/detail/stage15_move_uniform_dt_fastpath_results_20260605.md` | uniform-dt fast path 和 stage16 回退 |
| `doc/worklog/v2506/detail/stage17_post_fields_omp_local_accum_results_20260605.md` | post local accumulation |
| `doc/worklog/v2506/detail/stage18_buildCellOccupancy_sparse_clear_failed_20260605.md` | early occupancy failed attempt |
| `doc/worklog/v2506/detail/stage20_post_shared_cache_boundarymask_results_20260606.md` | shared sample cache + boundary mask |
| `doc/worklog/v2506/detail/stage21_move_nodelete_fastpath_results_20260606.md` | no-delete move fast path |
| `doc/worklog/v2506/detail/stage21_stability_rerun_results_20260606.md` | stage21 三次稳定性复跑 |
| `doc/worklog/v2506/detail/stage22_24_rejected_candidates_20260606.md` | rejected candidates |
| `doc/worklog/v2506/detail/stage25_buildCellOccupancy_sparse_dynamic_appended_results_20260606.md` | final occupancy stage 和性能归因 |
| `doc/worklog/v2506/detail/stage26_zb_move_coll_schedule_sweep_results_20260606.md` | `zb-cylinder-react/omp8` move/collision schedule sweep |
| `doc/worklog/v2506/detail/zb_cylinder_react_omp8_stage21_current_test_20260606.md` | `zb-cylinder-react/omp8` 当前版本验证 |

主要正式日志：

| log | 用途 |
|---|---|
| `ourmeshbkp/omp8/log.omp8.confirm_pdFalse_20260605_024610` | historical OMP baseline |
| `ourmesh/omp8/log.codex_omp_flat_post_steadywall_500step_20260605` | stage4 |
| `ourmesh/omp8/log.codex_omp_refretained_steadywall_500step_20260605` | stage5 |
| `ourmesh/omp8/log.codex_omp_move_fullport_bkpctrl_500step_20260605` | stage6 |
| `ourmesh/omp8/log.codex_omp_stage7_move_tracking_bkpctrl_500step_20260605` | stage7 |
| `ourmesh/omp8/log.codex_omp_stage8_move_omp_complete_bkpctrl_500step_20260605` | stage8 |
| `ourmesh/omp8/log.codex_omp_stage9_inline_reset_fastpath_bkpctrl_500step_20260605` | stage9 |
| `ourmesh/omp8/log.codex_omp_stage15_uniform_dt_fastpath_bkpctrl_500step_20260605` | stage15 |
| `ourmesh/omp8/log.codex_omp_stage17_post_fields_omp_local_accum_bkpctrl_500step_20260605` | stage17 |
| `ourmesh/omp8/log.codex_omp_stage20_post_shared_cache_boundarymask_bkpctrl_500step_20260605` | stage20 |
| `ourmesh/omp8/log.codex_omp_stage21_move_nodelete_fastpath_bkpctrl_500step_20260606` | stage21 best |
| `ourmesh/omp8/log.codex_omp_stage21_stability_run*_bkpctrl_500step_20260606` | stage21 stability |
| `ourmesh/omp8/log.codex_omp_stage25_buildocc_sparse_dynamic_appended_validomp_500step_20260606` | final stage25 |
| `zb-cylinder-react/omp8/log.codex_omp8_zb_buildocc_sparse_dynamic_appended_300step_20260606` | final `zb` verification |
| `zb-cylinder-react/omp8/schedule_sweep_stage26_20260606/log.codex_zb_*_300step_20260606` | stage26 schedule sweep 24-run matrix |

## 16. 最终建议

1. 将 stage25 作为当前 OMP 最终保留状态。
2. 如果要发布为稳定 OMP 基线，补做 stage25 的 2-3 次 `ourmesh/omp8`
   500-step formal 复跑。
3. 后续不要继续做无 profile 支撑的 move branch 微优化；move 已接近历史
   参考，必须先拆 subphase。
4. 后续 occupancy 工作只应围绕 `zb-cylinder-react/omp8` 的 v2506 gap 做
   subphase profile 后再改。
5. 保持 `profileSummary true`、`profileDetail false` 作为正式性能比较
   配置；诊断 profile 只用于定位，不混入正式基准。
6. 所有正式结论继续使用 500-step `ourmesh/omp8`，smoke 只作为编译/启动/
   correctness bring-up。
