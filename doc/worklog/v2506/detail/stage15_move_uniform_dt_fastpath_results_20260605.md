# Stage15 Move Uniform-dt Fast Path Results - 2026-06-05

## Scope

This stage continues the pure OMP move path from the current stage12+stage13+stage14 code state.

Implemented change:

- added `dsmcTimeStepModel::uniformDeltaT()` with `dsmcConstantTimeStepModel` returning true and `dsmcVariableTimeStepModel` returning false;
- cached `uniformDeltaT_` in `dsmcCloud`;
- in `dsmcParcel::move()`, used the already supplied `trackTime` for constant time-step cases instead of calling `cloud.deltaTValue(cell)` inside the parcel move loop;
- preserved the existing variable time-step path.

The target case uses `dsmcConstantTimeStepModel`, confirmed in both smoke and formal logs.

## Build

- `doc/worklog/v2506/detail/stage15_move_uniform_dt_fastpath_build_20260605.log`
- `doc/worklog/v2506/detail/stage15_move_uniform_dt_fastpath_rebuild_after_stage16_revert_20260605.log`

Both builds completed. The logs contain the existing OpenFOAM template-instantiation warnings, but no compile/link error markers.

## 10-Step OMP Smoke

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage15_uniform_dt_fastpath_smoke_10step_20260605`

Control setup:

- started from `ourmeshbkp/omp8/system/controlDict`;
- temporarily changed `endTime` to `1.e-06`;
- restored `ourmesh/omp8/system/controlDict` after the run.

Results:

| metric | stage14 smoke | stage15 smoke |
|---|---:|---:|
| move only | 0.860485809 s | 0.798181885 s |
| move+collide wall | 1.003752111 s | 0.93172987 s |
| full evolve wall | 1.739159432 s | 1.605409188 s |
| real | 6.81 s | 6.94 s |

Correctness smoke values:

| metric | value |
|---|---:|
| total iterations | 10 |
| collisions | 2856 |
| collision candidates | 2971 |
| DSMC particles | 2065867 |
| total energy | 1.150630665 |

No `FOAM FATAL`, `Fatal`, `NaN`, `Floating point`, `Segmentation`, or `MPI_ABORT` marker was found.

## 500-Step Formal OMP Run

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage15_uniform_dt_fastpath_bkpctrl_500step_20260605`

Control setup:

- formally used `ourmeshbkp/omp8/system/controlDict`;
- restored `ourmesh/omp8/system/controlDict` after the run.

Results:

| metric | stage10 | stage12 | stage15 |
|---|---:|---:|---:|
| real | 119.99 s | 116.55 s | 114.53 s |
| full evolve wall | 119.4492059 s | 118.4147219 s | 116.6882995 s |
| move only | 55.72693395 s | 54.62342851 s | 52.99979159 s |
| buildCellOccupancy | 7.66116264 s | 7.704320872 s | 8.663262117 s |
| collision phase | 3.641443767 s | 3.618495509 s | 3.628983484 s |
| post fields/output | 51.83930591 s | 51.88702306 s | 50.8090158 s |

Stage15 compared with stage12:

- real: `116.55 s -> 114.53 s`, 2.02 s faster;
- move only: `54.62342851 s -> 52.99979159 s`, 1.62363692 s faster;
- full evolve wall: `118.4147219 s -> 116.6882995 s`, 1.7264224 s faster.

Against the historical `ourmeshbkp/omp8` baseline:

| metric | historical bkp | stage15 | gap |
|---|---:|---:|---:|
| external wall / real | 103.66 s | 114.53 s | 10.87 s slower |
| full evolve wall | 102.1225033 s | 116.6882995 s | 14.5657962 s slower |
| move only | 48.82365912 s | 52.99979159 s | 4.17613247 s slower |
| buildCellOccupancy | 6.490896398 s | 8.663262117 s | 2.172365719 s slower |
| collision phase | 9.876225069 s | 3.628983484 s | 6.247241585 s faster |
| post fields/output | 36.92670834 s | 50.8090158 s | 13.88230746 s slower |

Correctness values at `Time = 5e-05`:

| metric | stage12 | stage15 |
|---|---:|---:|
| collisions | 35773 | 35458 |
| collision candidates | 68765 | 68600 |
| DSMC particles | 2463594 | 2463707 |
| total energy | 1.242990897 | 1.243139771 |

The stage15 formal log reached `Total Iterations = 500` and `End main`. No failure marker was found.

## Stage16 Failed Follow-up

An additional micro-optimization was tested after stage15:

- cached a reduced-D tracking flag in `dsmcCloud`;
- avoided constant-dt path `cell()` reads and replaced one division with a reciprocal multiply.

Build:

- `doc/worklog/v2506/detail/stage16_move_cached_dt_branch_build_20260605.log`

Smoke log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage16_cached_dt_branch_smoke_10step_20260605`

Stage16 smoke result:

| metric | stage15 smoke | stage16 smoke |
|---|---:|---:|
| move only | 0.798181885 s | 0.803218327 s |
| move+collide wall | 0.93172987 s | 0.939831001 s |
| full evolve wall | 1.605409188 s | 1.636234363 s |

Conclusion: stage16 was slower on smoke and was reverted. The current code and rebuilt binary are back to the stage15 state.

## ControlDict Restore

After smoke/formal runs and the stage16 revert, both controlDict files were restored to:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d  ourmesh/omp8/system/controlDict
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d  ourmesh/mpi8origin/system/controlDict
```

## Conclusion

Stage15 is a real move improvement but does not fully close the historical gap.

Current best formal pure OMP result:

- `real = 114.53 s`;
- `move only = 52.99979159 s`.

Remaining formal gap to historical `ourmeshbkp/omp8`:

- `real +10.87 s`;
- `move only +4.17613247 s`.

The remaining end-to-end gap is now dominated more by post fields/output and buildCellOccupancy than by collision. The remaining move gap is still present, but further blind micro-branch edits are not paying off; next work should use detailed move/tracking breakdown before trying another tracking change.
