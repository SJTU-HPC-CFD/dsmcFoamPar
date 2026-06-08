# Stage 5 OMP Reference-Retained Move and Collision Results

Date: 2026-06-05

Scope:
- OMP-only line.
- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8`
- Reference config/log source: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/omp8`
- Reference code stayed read-only: `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`

## Retained OMP Changes Ported

This stage stopped treating the tracker gate as a standalone proof target and directly ported OMP hot-path optimizations that are still retained in the reference code and fit the current OFv1706 interfaces.

Move path:
- added `dsmcCloud::trackerActive_` and `refreshTrackerUsage()`;
- gated face-transition tracker calls in `dsmcParcel::move()` and new-parcel insertion when no flux surface tracker is active;
- replaced boundary-model scans in the move/wall path with direct patch-to-model addressing;
- kept simple wall patch models (`dsmcDiffuseWallPatch`, `dsmcSpecularWallPatch`) on the direct OMP path and used a critical section only for patch models that are not known thread-safe;
- used per-thread move RNG when the OpenMP move kernel is active.

Collision path:
- added per-thread `noTimeCounter` buffers for subcells, parcel pointers, velocities, type ids, and charges;
- added an OpenMP cell-level collision loop with per-thread fast RNG for candidate selection;
- reused flat occupancy when available;
- retained the original serial collision logic as the fallback when OpenMP is disabled.

## Build

Command used before the smoke and formal runs:

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
```

Status:
- incremental build passed after the move-path tracker/boundary changes;
- incremental build passed after the parallel `noTimeCounter` collision changes;
- the only compile fix needed in this stage was using the OFv1706-compatible `autoPtr` access form instead of direct dereference.

## Logs And Case Cleanup

Logs:
- `ourmesh/omp8/log.codex_omp_refretained_smoke_10step_20260605`
- `ourmesh/omp8/log.codex_omp_refretained_steadywall_500step_20260605`

The formal run temporarily copied the `ourmeshbkp/omp8/system/controlDict` settings into `ourmesh/omp8/system/controlDict`.

Cleanup:
- original `ourmesh/omp8/system/controlDict` backup: `/tmp/hystrath_dlb_ourmesh_omp8_controlDict.pre_stage5_20260605`
- restored after the 500-step run;
- restored hash matched the backup:
  `ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d`

Error scan:
- no `Fatal`, `NaN`, `Floating point`, or `Segmentation` strings were found in the 10-step or 500-step retained-OMP logs.

## Smoke Result

10-step smoke:

| metric | value |
| --- | ---: |
| real | 7.36 s |
| user | 19.37 s |
| sys | 1.62 s |
| full evolve wall | 1.95439844 s |
| move+collide wall | 1.230075391 s |
| move only | 1.096597142 s |
| buildCellOccupancy | 0.104245321 s |
| collision phase | 0.019392573 s |
| post fields/output | 0.72418286 s |
| collisions | 2897 |
| collision candidates | 3001 |
| collision acceptance rate | 0.965344885 |
| particles | 2065859 |
| total energy | 1.150631839 |

OpenMP settings reported by the log:
- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

## Formal 500-Step Result

Current stage5 retained-OMP OMP8:

| metric | value |
| --- | ---: |
| real | 163.13 s |
| user | 1254.28 s |
| sys | 24.50 s |
| ClockTime | 162 s |
| full evolve wall | 167.3012128 s |
| move+collide wall | 116.0347194 s |
| move only | 104.1127577 s |
| buildCellOccupancy | 7.730683037 s |
| collision phase | 3.632583082 s |
| post fields/output | 51.2608737 s |
| move+collide cpu | 869.71 s |
| move only cpu | 780.19 s |
| buildCellOccupancy cpu | 57.33 s |
| collision phase cpu | 27.26 s |
| post fields/output cpu | 385.22 s |
| full evolve cpu | 1254.93 s |
| collisions | 35466 |
| collision candidates | 68727 |
| collision acceptance rate | 0.5160417303 |
| particles | 2463658 |
| stuck particles | 0 |
| total energy | 1.242958948 |

OpenMP settings reported by the log:
- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

## Performance Comparison

Against previous current-tree formal runs:

| metric | stage3 current | stage4 current | stage5 retained OMP |
| --- | ---: | ---: | ---: |
| real | 248.80 s | 265.13 s | 163.13 s |
| full evolve wall | unreliable older timer | 272.2832057 s | 167.3012128 s |
| move only | unreliable older timer | 144.8259732 s | 104.1127577 s |
| buildCellOccupancy | unreliable older timer | 7.89583818 s | 7.730683037 s |
| collision phase | unreliable older timer | 68.18518398 s | 3.632583082 s |
| post fields/output | unreliable older timer | 50.78234158 s | 51.2608737 s |

Stage5 vs stage4:
- `real`: `265.13 s -> 163.13 s`, 102.00 s faster, 38.47% lower, 1.625x speedup;
- `full evolve wall`: `272.2832057 s -> 167.3012128 s`, 38.56% lower;
- `move only`: `144.8259732 s -> 104.1127577 s`, 28.11% lower;
- `collision phase`: `68.18518398 s -> 3.632583082 s`, 94.67% lower;
- `post fields/output`: `50.78234158 s -> 51.2608737 s`, essentially unchanged and slightly slower by 0.48 s.

Stage5 vs stage3:
- `real`: `248.80 s -> 163.13 s`, 85.67 s faster, 34.43% lower, 1.525x speedup.

Against historical `ourmeshbkp/omp8`:

| metric | historical bkp | stage5 retained OMP | gap |
| --- | ---: | ---: | ---: |
| external/real wall | 103.66 s | 163.13 s | 59.47 s slower |
| full evolve wall | 102.1225033 s | 167.3012128 s | 65.1787095 s slower |
| move only | 48.82365912 s | 104.1127577 s | 55.28909858 s slower |
| buildCellOccupancy | 6.490896398 s | 7.730683037 s | 1.239786639 s slower |
| collision phase | 9.876225069 s | 3.632583082 s | 6.243641987 s faster |
| post fields/output | 36.92670834 s | 51.2608737 s | 14.33416536 s slower |
| collisions | 35847 | 35466 | -381 |
| particles | 2463391 | 2463658 | +267 |

The main remaining historical gap is now move, not collision. The current collision phase is faster than the historical bkp profile, while move is still 2.13x historical and post is 1.39x historical.

## Correctness Check

Completion and error checks:
- the retained-OMP formal log reached `Total Iterations = 500` and `End main`;
- final simulated time reached `Time = 5e-05`;
- final stuck particles remained `0`;
- no `Fatal`, `NaN`, `Floating point`, or `Segmentation` failure string was found in the retained-OMP 500-step log.

Current-tree 500-step comparison:

| metric | stage3 current | stage4 current | stage5 retained OMP | stage5 vs stage4 |
| --- | ---: | ---: | ---: | ---: |
| particles | 2463725 | 2463724 | 2463658 | -66 (-0.002679%) |
| stuck particles | 0 | 0 | 0 | 0 |
| collisions | 35283 | 35753 | 35466 | -287 (-0.802730%) |
| total energy | 1.243229348 | 1.24318294 | 1.242958948 | -0.000223992 (-0.018018%) |
| avg component energy | 1.009227367520e-18 | 1.009190104030e-18 | 1.009035303071e-18 | -1.548009590000e-22 (-0.015339%) |

Historical `ourmeshbkp/omp8` comparison:

| metric | historical bkp | stage5 retained OMP | difference |
| --- | ---: | ---: | ---: |
| particles | 2463391 | 2463658 | +267 (+0.010839%) |
| collisions | 35847 | 35466 | -381 (-1.062850%) |
| collision candidates | 75838 | 68727 | -7111 (-9.376566%) |
| collision acceptance rate | 0.4726786044 | 0.5160417303 | +0.0433631259 |
| avg total/component energy | 1.032208342000e-18 | 1.009035303071e-18 | -2.317303892900e-20 (-2.244996%) |

Assessment:
- stage5 is correct enough against the current-tree stage3/stage4 OMP baselines: particle drift is below 0.003%, collision-count drift is below 1%, and total-energy drift is below 0.022%;
- against historical bkp, particle and collision counts are still close, but collision candidates and average energy are not identical;
- the historical bkp log reports `Average total energy`, while the current-tree log reports individual average components plus `Total energy`, so the bkp energy comparison above uses the sum of current average components and should be treated as an approximate cross-log check, not strict identity;
- there is no NaN, stuck-particle growth, or early-exit signal in the retained-OMP run.

## Assessment

This stage is a real end-to-end OMP improvement:
- stage5 is much faster than both stage3 and stage4 current-tree OMP8;
- the retained parallel collision path is effective and should stay;
- the tracker/boundary move-path changes reduced move time, but move still dominates the remaining gap to historical bkp;
- flat post remains locally useful, but post is now the second-largest historical gap after move.

Next OMP-only direction:
1. Continue porting retained move-path optimizations from the reference tree.
2. Recheck whether the current move kernel still has shared writes, avoidable boundary criticals, or ordering/rebuild costs that are absent in the reference implementation.
3. After move is closer to historical bkp, revisit post fields/output; do not spend the next iteration on collision unless new evidence shows the 3.63 s result is unstable.
