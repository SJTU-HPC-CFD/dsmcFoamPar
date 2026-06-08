# Stage21 Move No-Delete Fast Path Results - 2026-06-06

## Scope

This stage continues from stage20.  It keeps the shared post-field sample cache
and only changes the pure OMP `Cloud::move()` path.

Implemented change:

- in the non-MPI OpenMP move kernel, count `keepParticle == false` while parcels
  are moved;
- when the count is zero, return through the existing ordered-parcel reuse path
  without the later survivor/delete counting scan and temporary survivor/delete
  lists;
- leave the MPI transfer branch unchanged.

The target `ourmesh/omp8` case has no stuck particles in the tested runs, so
this is a low-risk fast path for the common no-delete OMP move step.  If any
parcel is deleted, the original survivor/delete commit path is still used.

## Build

Build log:

- `doc/worklog/v2506/detail/stage21_move_nodelete_fastpath_build_20260606.log`

Result:

- `dsmcFoam+` rebuilt at
  `platforms/linux64IccDPInt32Opt/bin/dsmcFoam+`;
- `dsmcInitialise+` rebuilt at
  `platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+`;
- no compile/link error markers were found.  The build log contains the same
  OpenFOAM v1706 template-instantiation warnings seen in previous stages.

## 10-Step OMP Smoke

Control setup:

- temporarily used `ourmeshbkp/omp8/system/controlDict`;
- changed `endTime` to `1.e-06`;
- restored `ourmesh/omp8/system/controlDict` after the runs.

Logs:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage21_move_nodelete_fastpath_smoke_10step_20260606`
- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage21_move_nodelete_fastpath_smoke2_10step_20260606`
- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage21_move_nodelete_fastpath_smoke3_10step_20260606`

Smoke comparison:

| metric | stage20 smoke | stage21 smoke 1 | stage21 smoke 2 | stage21 smoke 3 | stage21 mean |
|---|---:|---:|---:|---:|---:|
| real | 6.58 s | 7.07 s | 6.62 s | 6.56 s | 6.75 s |
| full evolve wall | 1.519459887 s | 1.523759716 s | 1.440944734 s | 1.449909909 s | 1.471538120 s |
| move only | 0.868335554 s | 0.851717874 s | 0.825542506 s | 0.850265246 s | 0.842508542 s |
| buildCellOccupancy | 0.123102951 s | 0.113586003 s | 0.107953341 s | 0.107228537 s | 0.109589294 s |
| collision phase | 0.020993557 s | 0.022261412 s | 0.016799953 s | 0.016667453 s | 0.018576273 s |
| post fields/output | 0.496770673 s | 0.525342541 s | 0.479969612 s | 0.465018729 s | 0.490110294 s |

The smoke signal was noisy in external `real`, but the profiled move path was
consistently lower on average, so a 500-step formal run was warranted.

No `FOAM FATAL`, `Fatal`, `NaN`, `Floating point`, `Segmentation`, or
`MPI_ABORT` marker was found in the smoke logs.

## 500-Step Formal OMP Run

Control setup:

- temporarily used the full `ourmeshbkp/omp8/system/controlDict`;
- restored `ourmesh/omp8/system/controlDict` after the run.

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage21_move_nodelete_fastpath_bkpctrl_500step_20260606`

Formal profile versus stage20:

| metric | stage20 | stage21 | delta |
|---|---:|---:|---:|
| real | 103.34 s | 100.65 s | 2.69 s faster |
| full evolve wall | 102.6520492 s | 93.28696662 s | 9.36508258 s faster |
| move only | 56.13548084 s | 51.17943521 s | 4.95604563 s faster |
| buildCellOccupancy | 9.115554564 s | 8.087915571 s | 1.027638993 s faster |
| collision phase | 3.760305061 s | 3.314541395 s | 0.445763666 s faster |
| post fields/output | 33.03582209 s | 30.15687469 s | 2.87894740 s faster |

Formal profile versus historical `ourmeshbkp/omp8`:

| metric | historical bkp | stage21 | delta |
|---|---:|---:|---:|
| real | 103.66 s | 100.65 s | 3.01 s faster |
| full evolve wall | 102.1225033 s | 93.28696662 s | 8.83553668 s faster |
| move only | 48.82365912 s | 51.17943521 s | 2.35577609 s slower |
| buildCellOccupancy | 6.490896398 s | 8.087915571 s | 1.597019173 s slower |
| collision phase | 9.876225069 s | 3.314541395 s | 6.561683674 s faster |
| post fields/output | 36.92670834 s | 30.15687469 s | 6.76983365 s faster |

Correctness values at `Time = 5e-05`:

| metric | stage20 | stage21 |
|---|---:|---:|
| collisions | 35624 | 36005 |
| collision candidates | 68756 | 69029 |
| collision acceptance rate | 0.5181220548 | 0.5215923742 |
| DSMC particles | 2463612 | 2463537 |
| stuck particles | 0 | 0 |
| total energy | 1.24300724 | 1.243069478 |
| total iterations | 500 | 500 |

The formal run reached `End main`.  No `FOAM FATAL`, `Fatal`, `NaN`,
`Floating point`, `Segmentation`, or `MPI_ABORT` marker was found.

## ControlDict Restore

After smoke and formal runs, `ourmesh/omp8/system/controlDict` was restored to:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d  ourmesh/omp8/system/controlDict
```

## Conclusion

Stage21 is kept.  It is a real OMP move-path improvement and also improves the
current formal end-to-end result:

- `real = 100.65 s`;
- `full evolve wall = 93.28696662 s`;
- `move only = 51.17943521 s`;
- `post fields/output = 30.15687469 s`.

The overall result is now faster than the historical `ourmeshbkp/omp8` profile.
The remaining formal gaps are still in `move only` and `buildCellOccupancy`;
post fields/output is now clearly ahead of both stage20 and the historical
baseline.
