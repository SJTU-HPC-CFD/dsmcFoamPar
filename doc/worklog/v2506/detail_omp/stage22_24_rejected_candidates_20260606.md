# Stage22-24 rejected OMP candidates - 2026-06-06

## Context

Stage21 remains the kept OMP baseline before these experiments:

- source change: pure OMP `Cloud::move()` no-delete fast path;
- formal log: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage21_move_nodelete_fastpath_bkpctrl_500step_20260606`;
- key formal profile: `real = 100.65 s`, `full evolve wall = 93.28696662 s`, `move only = 51.17943521 s`, `buildCellOccupancy = 8.087915571 s`.

All candidates below were tested with the same `ourmesh/omp8` case, temporarily
using `ourmeshbkp/omp8/system/controlDict` for 500 steps, then restoring
`ourmesh/omp8/system/controlDict`.

Restored controlDict hash after the runs:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d  ourmesh/omp8/system/controlDict
```

## Stage22: buildCellOccupancy active-cell reset

Candidate:

- in `dsmcCloud::buildCellOccupancy()`, reset per-thread `localCounts` only for the previous active cells when `moveOrderedReady`;
- intended to reduce `ompNumThreads * nCells` zeroing.

Build log:

- `doc/worklog/v2506/detail/stage22_buildCellOccupancy_active_reset_build_20260606.log`

Formal log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage22_buildocc_active_reset_bkpctrl_500step_20260606`

Formal result:

| metric | stage21 kept | stage22 candidate | delta |
|---|---:|---:|---:|
| real | 100.65 s | 99.41 s | 1.24 s faster |
| full evolve wall | 93.28696662 s | 98.63267366 s | 5.34570704 s slower |
| move only | 51.17943521 s | 54.12781802 s | 2.94838281 s slower |
| buildCellOccupancy | 8.087915571 s | 8.651435857 s | 0.563520286 s slower |
| collision phase | 3.314541395 s | 3.566299609 s | 0.251758214 s slower |
| post fields/output | 30.15687469 s | 31.71996632 s | 1.56309163 s slower |

Decision: rejected and reverted.  The external `real` value improved, but the
solver's internal profiled path and the targeted `buildCellOccupancy` metric
both regressed.

Revert build log:

- `doc/worklog/v2506/detail/stage22_buildCellOccupancy_active_reset_revert_build_20260606.log`

## Stage23: move delete-index recording

Candidate:

- replace full `keepParticleFlags` writes in the pure OMP non-MPI move branch
  with per-thread deleted-particle index recording;
- intended to avoid writing one flag per particle on no-delete steps.

Build log:

- `doc/worklog/v2506/detail/stage23_move_delete_indices_build_20260606.log`

Formal log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage23_move_delete_indices_bkpctrl_500step_20260606`

Formal result:

| metric | stage21 kept | stage23 candidate | delta |
|---|---:|---:|---:|
| real | 100.65 s | 95.81 s | 4.84 s faster |
| full evolve wall | 93.28696662 s | 98.47757234 s | 5.19060572 s slower |
| move only | 51.17943521 s | 53.54094492 s | 2.36150971 s slower |
| buildCellOccupancy | 8.087915571 s | 8.709473857 s | 0.621558286 s slower |
| collision phase | 3.314541395 s | 3.608371929 s | 0.293830534 s slower |
| post fields/output | 30.15687469 s | 32.04237918 s | 1.88550449 s slower |

Decision: rejected and reverted.  The smoke signal was positive, but the formal
profile showed the delete path is not free enough for this approach.

Revert build log:

- `doc/worklog/v2506/detail/stage23_move_delete_indices_revert_build_20260606.log`

## Stage24: uninitialized move flags

Candidate:

- construct `keepParticleFlags` and `switchProcessorFlags` without an initial
  fill value because every element is written before read;
- intended to remove redundant flag-array initialization.

Build log:

- `doc/worklog/v2506/detail/stage24_move_uninitialized_flags_build_20260606.log`

Formal log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage24_move_uninit_flags_bkpctrl_500step_20260606`

Because recent 500-step runs were slower than the older stage21 formal run, a
fresh current-condition stage21 rerun was made after reverting stage24:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage21_rerun_after_stage24_revert_bkpctrl_500step_20260606`

Current-condition comparison:

| metric | stage21 rerun | stage24 candidate | delta |
|---|---:|---:|---:|
| real | 97.31 s | 99.53 s | 2.22 s slower |
| full evolve wall | 99.18206278 s | 99.01781108 s | 0.16425170 s faster |
| move only | 54.12212295 s | 53.97329809 s | 0.14882486 s faster |
| buildCellOccupancy | 8.666368709 s | 8.872791784 s | 0.206423075 s slower |
| collision phase | 3.584365182 s | 3.554273378 s | 0.030091804 s faster |
| post fields/output | 32.23474339 s | 32.04972285 s | 0.18502054 s faster |

Decision: rejected and reverted.  The internal move-only metric improved only
slightly, while external `real` and `buildCellOccupancy` regressed.

Revert build log:

- `doc/worklog/v2506/detail/stage24_move_uninit_flags_revert_build_20260606.log`

## Correctness notes

All formal candidate runs reached `End main` and `Total Iterations = 500`.
No `FOAM FATAL`, `Fatal`, `NaN`, `Floating point`, `Segmentation`, or
`MPI_ABORT` marker was found in the candidate formal logs.

Last-step values were in the same stochastic band as stage21:

| run | collisions | candidates | particles | stuck | total energy |
|---|---:|---:|---:|---:|---:|
| stage21 rerun | 35674 | 68747 | 2463809 | 0 | 1.243072026 |
| stage22 candidate | 35471 | 68577 | 2463621 | 0 | 1.242988757 |
| stage23 candidate | 35620 | 68797 | 2463656 | 0 | 1.24306114 |
| stage24 candidate | 35564 | 68634 | 2463704 | 0 | 1.24306551 |

## Current state

Stage22, stage23, and stage24 source changes are not kept.  The source and
rebuilt local binaries are back on the stage21 kept implementation.
