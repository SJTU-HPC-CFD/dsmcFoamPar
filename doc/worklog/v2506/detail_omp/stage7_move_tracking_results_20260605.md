# Stage7 OMP Move Tracking Hot-Path Results - 2026-06-05

## Scope

This stage continued only on the OMP line after stage6.  It did not change the
reference tree under `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`.

Changes applied in this stage:

- `dsmcParcel::move()` now caches the `dsmcCloud`, mesh, boundary mapping,
  coordinate-system type, and tracker-active gate outside the per-parcel
  tracking loop;
- the boundary face path now computes the patch index once and uses direct
  `cyclicBoundaryToModelIds()` / `patchToModelIds()` lookup;
- the DSMC-specific `particle::trackToFace(..., true)` path now caches mesh
  owner/neighbour/tet-base/cell-volume/boundary references and the moving-mesh
  flag;
- the DSMC `trackToFace(..., true)` path now skips `hitWallFaces()` entirely
  when `hasWallImpactDistance()` is false.  For this current `dsmcCloud`, that
  is the inherited default, so this removes a repeated no-op call from the move
  hot path.

The v2506 `trackToAndHitFace()` API was still not copied into OFv1706; this
stage stayed inside the old OFv1706 tracking API.

## Build

Command:

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
```

Log:

- `doc/worklog/v2506/detail/stage7_move_tracking_build_20260605.log`

Result:

- build completed;
- `dsmcFoam+` rebuilt at `platforms/linux64IccDPInt32Opt/bin/dsmcFoam+`;
- `dsmcInitialise+` rebuilt at `platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+`;
- only the known OpenFOAM v1706 template-instantiation warnings were present.

## 10-Step OMP Smoke

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage7_move_tracking_smoke_10step_20260605`

Temporary control settings:

- copied `ourmeshbkp/omp8/system/controlDict` into `ourmesh/omp8/system/controlDict`;
- changed `endTime` to `1.e-06`;
- restored the original `ourmesh/omp8/system/controlDict` after the run.

Result:

| metric | value |
| --- | ---: |
| real | 7.28 s |
| user | 18.73 s |
| sys | 1.48 s |
| full evolve wall | 1.896343598 s |
| move+collide wall | 1.210692305 s |
| move only | 1.081026197 s |
| buildCellOccupancy | 0.102229993 s |
| collision phase | 0.017659035 s |
| post fields/output | 0.685546801 s |
| collisions | 2804 |
| collision candidates | 2939 |
| particles | 2065825 |
| stuck particles | 0 |
| total energy | 1.150609378 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

No `Fatal`, `NaN`, `Floating point`, or `Segmentation` entries were reported.

## 500-Step Formal Run

Formal comparison used the same control settings as stage5 and stage6.  For
this run, `ourmeshbkp/omp8/system/controlDict` was temporarily copied to
`ourmesh/omp8/system/controlDict`, then the original
`ourmesh/omp8/system/controlDict` was restored.

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage7_move_tracking_bkpctrl_500step_20260605`

Result:

| metric | value |
| --- | ---: |
| real | 151.06 s |
| user | 1172.88 s |
| sys | 19.71 s |
| ClockTime | 151 s |
| full evolve wall | 156.8234055 s |
| move+collide wall | 106.8374108 s |
| move only | 95.5824105 s |
| buildCellOccupancy | 7.194255854 s |
| collision phase | 3.516419771 s |
| post fields/output | 49.98115219 s |
| move+collide cpu | 796.60 s |
| move only cpu | 712.71 s |
| buildCellOccupancy cpu | 52.76 s |
| collision phase cpu | 26.16 s |
| post fields/output cpu | 373.29 s |
| full evolve cpu | 1169.91 s |
| collisions | 35601 |
| collision candidates | 68574 |
| collision acceptance rate | 0.5191617814 |
| particles | 2463537 |
| stuck particles | 0 |
| total energy | 1.242846706 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

No `Fatal`, `NaN`, `Floating point`, or `Segmentation` entries were reported.

## Stage6 Comparison

| metric | stage6 full move port | stage7 move tracking | delta |
| --- | ---: | ---: | ---: |
| real | 161.79 s | 151.06 s | 10.73 s faster, 6.63% lower |
| full evolve wall | 167.0823348 s | 156.8234055 s | 10.2589293 s faster, 6.14% lower |
| move only | 103.6495219 s | 95.5824105 s | 8.0671114 s faster, 7.78% lower |
| buildCellOccupancy | 7.483010573 s | 7.194255854 s | 0.288754719 s faster, 3.86% lower |
| collision phase | 3.613524884 s | 3.516419771 s | 0.097105113 s faster, 2.69% lower |
| post fields/output | 51.76373504 s | 49.98115219 s | 1.78258285 s faster, 3.44% lower |

This is the first post-stage5 move-only change that materially reduced the
formal 500-step move time.  The main gain is still move, not collision.

## Correctness Check

Against stage6 full move port:

| metric | stage6 | stage7 | difference |
| --- | ---: | ---: | ---: |
| particles | 2463609 | 2463537 | -72 (-0.002923%) |
| collisions | 35898 | 35601 | -297 (-0.827344%) |
| collision candidates | 68794 | 68574 | -220 (-0.319795%) |
| total energy | 1.242916724 | 1.242846706 | -0.000070018 (-0.005633%) |
| stuck particles | 0 | 0 | unchanged |

The run completed 500 iterations and reached `End main`.  The correctness drift
relative to stage6 is small for the current-tree OMP regression level.

## Historical Gap

Against historical `ourmeshbkp/omp8`:

| metric | historical bkp | stage7 move tracking | gap |
| --- | ---: | ---: | ---: |
| external/real wall | 103.66 s | 151.06 s | 47.40 s slower |
| full evolve wall | 102.1225033 s | 156.8234055 s | 54.7009022 s slower |
| move only | 48.82365912 s | 95.5824105 s | 46.75875138 s slower |
| buildCellOccupancy | 6.490896398 s | 7.194255854 s | 0.703359456 s slower |
| collision phase | 9.876225069 s | 3.516419771 s | 6.359805298 s faster |
| post fields/output | 36.92670834 s | 49.98115219 s | 13.05444385 s slower |
| particles | 2463391 | 2463537 | +146 (+0.005927%) |
| collisions | 35847 | 35601 | -246 (-0.686250%) |
| collision candidates | 75838 | 68574 | -7264 (-9.578312%) |

The remaining historical gap is still dominated by move.  Collision remains
faster than the historical bkp run.

## Cleanup

`ourmesh/omp8/system/controlDict` was restored after the smoke and formal runs.

Restored hash:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d
```

No `dsmcFoam+` process was left running by the stage7 commands.
