# zb 300-step serial collision baseline, 2026-06-12

## Scope

- Case source: `run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/omp8`
- Run mode: single process, no MPI, `useOpenMP false`, `openmpThreads 1`, `OMP_NUM_THREADS=1`
- Steps: 300 (`endTime 1.9920146682e-05`, `deltaT 6.640048894e-08`)
- Output control: `writeControl runTime`, `writeInterval 1.e-3`, no scheduled field write in this 300-step window
- Source state: local dirty tree with owned collision-cell iterator related edits; `src/lagrangian/dsmc` rebuilt before the run
- Temporary copied case was removed after the run; retained artifacts are the log and controlDict snapshots.

## Artifacts

- Serial log: `doc/worklog/v2506/detail_mix/zb_serial_baseline_20260612/logs/zb_SERIAL_300step_20260612.log`
- Current OMP8 log: `doc/worklog/v2506/detail_mix/zb_serial_baseline_20260612/logs/zb_OMP8_current_300step_20260612.log`
- Current MPI8 replicated log: `doc/worklog/v2506/detail_mix/zb_serial_baseline_20260612/logs/zb_MPI8_current_300step_20260612.log`
- Exit codes: matching `.exit` files under `doc/worklog/v2506/detail_mix/zb_serial_baseline_20260612/logs`
- Run controlDicts: `doc/worklog/v2506/detail_mix/zb_serial_baseline_20260612/controlDicts`
- ControlDict diff vs copied OMP8 case: only `useOpenMP true -> false` and `openmpThreads 8 -> 1`.

## Serial result

| metric | serial |
| --- | ---: |
| exit | 0 |
| external real | 342.78 s |
| iterations | 300 |
| full evolve wall | 370.9618039 s |
| move+collide wall | 370.7683617 s |
| move only | 285.3681810 s |
| buildCellOccupancy | 35.11854775 s |
| collision phase | 49.29919808 s |
| post fields/output | 0.191958561 s |
| OpenMP enabled | 0 |
| final particles | 1958171 |
| stuck particles | 0 |
| final collisions / candidates | 241613 / 475278 |
| total energy | 0.001962662871 |

Note: for total wall-time comparison use external `real`. The serial profile `full evolve wall` is larger than external `real`, so phase totals are treated as hotspot indicators rather than a replacement for `/usr/bin/time`.

## Current-Source Single-Run Comparison

These three runs were all made after rebuilding `src/lagrangian/dsmc` in the current local tree.

| mode | real | full evolve | move | buildOcc | collision | post | final particles | collisions / candidates | energy | stuck |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Serial, 1 core | 342.78 | 370.96 | 285.37 | 35.12 | 49.30 | 0.19 | 1958171 | 241613 / 475278 | 0.001962662871 | 0 |
| OMP8 current | 74.05 | 72.01 | 52.78 | 6.19 | 11.39 | 0.26 | 1957853 | 242094 / 477721 | 0.001961627640 | 0 |
| MPI8 replicated current | 100.53 | 103.37 | 69.87 | 7.85 | 42.58 | 1.60 | 1959041 | 198452 / 385483 | 0.001932343654 | 0 |

Current-source speedup against serial:

| mode | real speedup | full-evolve speedup | move speedup | buildOcc speedup | collision speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| OMP8 current | 4.63x | 5.15x | 5.41x | 5.67x | 4.33x |
| MPI8 replicated current | 3.41x | 3.59x | 4.08x | 4.47x | 1.16x |

Current-source OMP8 vs MPI8 replicated:

| metric | MPI8 / OMP8 |
| --- | ---: |
| real | 1.36x slower |
| full evolve | 1.44x slower |
| move | 1.32x slower |
| buildOcc | 1.27x slower |
| collision | 3.74x slower |

MPI8 replicated balance diagnostics in this current run:

- rank evolve wall max/min: `1.004031819`
- final particles per rank max/min: `1.422579507`
- Phase C checks / rebalances / triggered checks: `300 / 2 / 2`
- migration wall: `2.866175552 s`

Interpretation for current source: owned collision-cell traversal is helping relative to the old MPI8 mean, but not enough to reproduce OMP's collision speedup. The current pure MPI8 replicated collision phase is `42.58 s`, only `1.16x` faster than serial collision, while OMP8 is `4.33x` faster than serial collision. Since rank evolve wall max/min is almost flat (`1.004`), this run does not support "gross rank imbalance" as the main explanation.

## Historical 2026-06-10 Three-Run Mean

The OMP8/MPI8 rows below are the existing three-run means from
`doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_rerun/logs`.
They are not a strict same-source rerun of the 2026-06-12 local tree. They are still useful for the load-effect question because the serial `zb` collision workload is now measured directly.

| mode | real | full evolve | move | buildOcc | collision | post |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Serial, 1 core | 342.78 | 370.96 | 285.37 | 35.12 | 49.30 | 0.19 |
| OMP8 mean | 72.78 | 66.87 | 37.38 | 7.58 | 9.40 | 11.37 |
| MPI8 replicated mean | 107.31 | 101.31 | 51.03 | 14.19 | 49.16 | 14.97 |

Speedup against this serial baseline:

| mode | real speedup | full-evolve speedup | move speedup | buildOcc speedup | collision speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| OMP8 mean | 4.71x | 5.55x | 7.63x | 4.63x | 5.25x |
| MPI8 replicated mean | 3.19x | 3.66x | 5.59x | 2.47x | 1.00x |

## Interpretation

1. The `zb` case really has a long enough serial collision section: 49.30 s over 300 steps, about 14.4% of external real time and much larger than the smaller `ourmesh` OMP collision numbers.
2. OMP8 gives a real collision-kernel speedup: current-source `49.30 / 11.39 = 4.33x`, historical mean `49.30 / 9.40 = 5.25x`. It is not 8x because collision has scheduling, random-number, memory, per-cell serial work, and load-balance losses, but the critical path is clearly shortened.
3. Current pure MPI8 replicated collision is only `49.30 / 42.58 = 1.16x` faster than serial. Historical pure MPI8 replicated mean was essentially no collision speedup: `49.30 / 49.16 = 1.00x`.
4. Therefore the OMP8 vs MPI8 collision gap cannot be explained by "the collision load is too small" or by missing serial baseline. The case has enough collision work, and MPI8 still does not shorten the collision critical path much.
5. The current MPI8 run has rank evolve wall max/min `1.004`, so gross rank wall imbalance is not the main explanation. The remaining gap is in the replicated collision path itself: ownership partitioning does not equal intra-rank kernel parallelization, and per-rank traversal/bookkeeping plus DLB/candidate-statistics overhead still sit on the critical path.
6. MPI8 total still improves over serial (`3.41x` current) because move and build occupancy are reduced, but collision is not carrying the same speedup as OMP.

## Next Check

The next MPI optimization should target the pure replicated collision critical path, not scheduling configuration alone. The current evidence points at:

- clearing/storing `nCandidatesPerCell` without stale full-array side effects for DLB accounting;
- rechecking active owned collision-cell traversal in all MPI collision branches;
- reducing replicated collision bookkeeping that remains serial per rank;
- rerunning at least three current-source `zb` MPI8 reps, because historical MPI8 collision ranged from `34.14 s` to `62.80 s`.
