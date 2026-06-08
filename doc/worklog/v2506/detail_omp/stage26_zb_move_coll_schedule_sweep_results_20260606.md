# Stage26: zb-cylinder-react/omp8 move/collision schedule sweep

日期：2026-06-06

目标：在 `zb-cylinder-react/omp8` 上分别测试 `move` 和 `collision`
的 `static` / `dynamic` OpenMP schedule，并扫描 chunk size
`8,16,32,64,128,256`。该 case 中 collision 占比相对更容易观察
schedule/chunk 影响。

## 测试对象

case：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/omp8
```

代码状态：stage25 当前代码，即 `buildCellOccupancy` sparse reset +
`DynamicList` reuse + appended-aware ordered path 后的版本。

运行方式：

```text
source doc/scripts/env.sh
OMP_NUM_THREADS=8 /usr/bin/time -p dsmcFoam+
```

本次没有重新初始化 case，直接从 `startTime 0` 运行。`controlDict`
为 300-step 配置：

```text
endTime 1.9920146682e-05;
deltaT 6.640048894e-08;
useOpenMP true;
openmpThreads 8;
profileSummary true;
profileDetail false;
```

## 测试矩阵

本次不是完整 `move schedule x move chunk x collision schedule x
collision chunk` 四维矩阵，而是两个独立 sweep：

1. move sweep：固定 collision 为 `dynamic / 8`，只扫描 move
   `static/dynamic x chunk 8..256`。
2. collision sweep：固定 move 为 `static / 64`，只扫描 collision
   `static/dynamic x chunk 8..256`。

共 24 个 300-step run。

批量脚本：

```text
doc/worklog/v2506/detail/run_zb_move_coll_schedule_sweep_20260606.sh
```

日志和表格目录：

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/omp8/schedule_sweep_stage26_20260606
```

关键文件：

```text
manifest.tsv
summary.tsv
controlDict.before_sweep
controlDict.after_restore
log.codex_zb_*_300step_20260606
```

文件 hash：

```text
summary.tsv sha256 = 31017a93c2e229b62abce268d98b4447c1601f4f9a9e0e36842a35f89b41cb4e
run script sha256 = 86985c35aabedd6f65d3963bb9ed1014b6980a56acc38a329f882397993f4fd9
```

`controlDict` sweep 前后已恢复一致：

```text
system/controlDict sha256 = 1f1f20be1740bae0c26637c8ce1e233951963114f62d486ef1b77da439f52b75
openmpMoveSchedule static;
openmpMoveChunk 64;
openmpCollisionSchedule dynamic;
openmpCollisionChunk 8;
```

## 正确性检查

结果：

- `manifest.tsv` 共有 25 行，即 24 个组合 + 表头。
- 24 个 run 的 `exitCode` 全部为 `0`。
- 24 个 run 均到达 `Iteration 300`。
- 错误关键字扫描无命中：
  `FOAM FATAL` / `nan` / `NaN` / `MPI_ABORT` / `Segmentation` /
  `SIGSEGV` / `core dumped` / `Floating point exception`。
- final `Collisions` 范围：`241056 -> 242836`，span `1780`，
  mean `242109.4`，span/mean `0.735%`。

碰撞数用于 sanity check，不作为 bitwise 一致性标准；不同 schedule 会改变
并行执行顺序和随机碰撞路径。

## Move sweep 结果

固定：

```text
openmpCollisionSchedule dynamic;
openmpCollisionChunk 8;
```

同组基准：

```text
move static / 64
full evolve wall = 78.06102765 s
move only        = 43.85303096 s
```

按 `full evolve wall` 排序：

| move schedule/chunk | full evolve wall [s] | move only [s] | buildCellOccupancy [s] | collision phase [s] | real [s] |
|---|---:|---:|---:|---:|---:|
| static / 256 | 77.11135257 | 42.88355957 | 8.49646293 | 10.50747748 | 80.89 |
| static / 128 | 77.50246610 | 43.22449083 | 8.755019848 | 10.43021164 | 81.89 |
| static / 64 | 78.06102765 | 43.85303096 | 8.335356222 | 10.51135249 | 84.64 |
| static / 32 | 78.43745217 | 44.35011103 | 8.320832266 | 10.51031655 | 82.61 |
| dynamic / 128 | 78.71659021 | 45.22702633 | 8.125458627 | 10.36533293 | 82.03 |
| dynamic / 256 | 79.01221285 | 45.26499047 | 8.286185276 | 10.42439840 | 82.23 |
| static / 16 | 79.28067769 | 45.41806557 | 8.232832785 | 10.52402565 | 83.40 |
| static / 8 | 80.29099797 | 46.73535309 | 7.939946204 | 10.43405242 | 84.78 |
| dynamic / 16 | 80.75199882 | 46.72304560 | 8.442430365 | 10.47823073 | 87.59 |
| dynamic / 64 | 81.06350078 | 47.11674827 | 8.208773984 | 10.49659521 | 84.87 |
| dynamic / 32 | 82.83674084 | 48.17373349 | 8.582965528 | 10.61897487 | 87.01 |
| dynamic / 8 | 83.06845838 | 48.40619232 | 8.836215193 | 10.60091996 | 87.33 |

相对 `static / 64`：

| move schedule/chunk | full delta | move-only delta |
|---|---:|---:|
| static / 8 | +2.86% | +6.57% |
| static / 16 | +1.56% | +3.57% |
| static / 32 | +0.48% | +1.13% |
| static / 64 | +0.00% | +0.00% |
| static / 128 | -0.72% | -1.43% |
| static / 256 | -1.22% | -2.21% |
| dynamic / 8 | +6.41% | +10.38% |
| dynamic / 16 | +3.45% | +6.54% |
| dynamic / 32 | +6.12% | +9.85% |
| dynamic / 64 | +3.85% | +7.44% |
| dynamic / 128 | +0.84% | +3.13% |
| dynamic / 256 | +1.22% | +3.22% |

Move 结论：

- `dynamic` move 在这个单 rank OMP8 case 上没有收益，所有 chunk 都慢于
  `static / 64` 的 move-only 时间。
- `static` chunk 从 8 增大到 256 时，move-only 基本改善：
  `46.73535309 -> 42.88355957 s`。
- 本次 `zb` 上的 move 最优是 `static / 256`：
  full wall 比 `static / 64` 低 `1.22%`，move-only 低 `2.21%`。
- 这只说明 `zb-cylinder-react/omp8` 300-step 上 `static / 256` 更好；
  若要改正式生产配置，还需要在 `ourmesh/omp8` 500-step 上复测。

## Collision sweep 结果

固定：

```text
openmpMoveSchedule static;
openmpMoveChunk 64;
```

同组基准：

```text
collision dynamic / 8
full evolve wall = 75.45810589 s
collision phase  = 10.35458841 s
```

按 `collision phase` 排序：

| collision schedule/chunk | full evolve wall [s] | move only [s] | buildCellOccupancy [s] | collision phase [s] | real [s] |
|---|---:|---:|---:|---:|---:|
| dynamic / 8 | 75.45810589 | 42.06382167 | 7.970480739 | 10.35458841 | 81.17 |
| dynamic / 16 | 76.72571503 | 42.60142367 | 8.545208794 | 10.43977883 | 80.13 |
| dynamic / 32 | 76.17704143 | 42.35313915 | 8.199493522 | 10.54607559 | 79.58 |
| dynamic / 64 | 75.53186369 | 41.88465382 | 7.765691990 | 10.67298286 | 81.64 |
| dynamic / 128 | 76.65006902 | 42.31881036 | 8.295323981 | 10.89192718 | 80.18 |
| dynamic / 256 | 77.21817424 | 42.45125659 | 8.316535153 | 11.44430031 | 80.76 |
| static / 8 | 76.86783591 | 42.18998927 | 8.138643928 | 11.57134530 | 82.97 |
| static / 16 | 78.31582568 | 42.86002794 | 8.595644551 | 11.76955654 | 80.55 |
| static / 32 | 77.40826851 | 42.38330788 | 8.064784944 | 11.82185272 | 83.26 |
| static / 64 | 78.83915353 | 42.87104388 | 8.635456359 | 12.25675103 | 82.23 |
| static / 128 | 78.97288225 | 42.61511125 | 8.512498083 | 12.87052694 | 82.44 |
| static / 256 | 79.33333096 | 42.32312361 | 7.946598657 | 13.99113719 | 82.92 |

相对 `dynamic / 8`：

| collision schedule/chunk | full delta | collision-phase delta |
|---|---:|---:|
| static / 8 | +1.87% | +11.75% |
| static / 16 | +3.79% | +13.67% |
| static / 32 | +2.58% | +14.17% |
| static / 64 | +4.48% | +18.37% |
| static / 128 | +4.66% | +24.30% |
| static / 256 | +5.14% | +35.12% |
| dynamic / 8 | +0.00% | +0.00% |
| dynamic / 16 | +1.68% | +0.82% |
| dynamic / 32 | +0.95% | +1.85% |
| dynamic / 64 | +0.10% | +3.07% |
| dynamic / 128 | +1.58% | +5.19% |
| dynamic / 256 | +2.33% | +10.52% |

Collision 结论：

- `dynamic / 8` 是本组 collision phase 最优：`10.35458841 s`。
- `dynamic` chunk 增大后 collision phase 逐步变差：
  `dynamic / 256` 比 `dynamic / 8` 慢 `10.52%`。
- `static` collision 明显不适合这个 case：
  最好的 `static / 8` 的 collision phase 仍比 `dynamic / 8` 慢
  `11.75%`；`static / 256` 慢 `35.12%`。
- `dynamic / 64` 的 full wall 和 `dynamic / 8` 很接近
  (`+0.10%`)，但 collision phase 已慢 `3.07%`，不建议仅凭
  full wall 把 collision chunk 放大。

## 总结

本次 `zb-cylinder-react/omp8` 300-step schedule sweep 的结论：

```text
move:      static / 256 在本 case 上最好，但需要 ourmesh/omp8 500-step 复测后才能改正式配置。
collision: dynamic / 8 明确最好，应继续作为当前 OMP8 collision 推荐配置。
```

边界：

- 这是 single-run sweep，未做每个组合多次重复。
- `moveSweep` 和 `collisionSweep` 中同名基准的绝对 wall time 存在运行顺序/
  系统状态漂移，因此本文只在各自 sweep 内部比较百分比。
- 当前 `ourmesh/omp8` 正式 500-step 配置仍保持
  `move static / 64 + collision dynamic / 8`，除非后续专门复测
  `move static / 256`。
