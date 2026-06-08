# Stage 20 - Shared-cache post fields with boundary patch mask

## Scope

This stage continues the stage19 `dsmcVolFields` shared sample cache port from
the reference tree:

`/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`

The kept requirement is the reference-style shared sample cache data flow.  This
stage does not revert to the earlier stage17 direct local-accumulation path.

Files changed:

- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`
- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H`

## Diagnostic basis

Stage19 active-cell shared-cache profile log:

`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_stage19_post_shared_cache_activecells_profile_smoke_10step_20260605`

Key 10-step profile-detail numbers:

| item | time |
|---|---:|
| shared cache build | 0.294076519 s |
| parcel accumulate | 0.29040268 s |
| field combine | 0.058142003 s |
| cell reduction | 0.004870337 s |
| boundary accumulation | 0.196422763 s |
| post fields/output | 0.694367741 s |
| real | 7.40 s |

The active-cell path was already working, so the remaining cost was not an
all-cell empty scan.  The biggest concrete avoidable cost was boundary
accumulation over non-DSMC boundary patches.  In `ourmesh/omp8`,
`BaseAndTop` is a `symmetry` patch with 208302 faces.

## Implementation

- Kept the shared sample cache and active-cell reset/combine path.
- Moved `cloud.nParticles(celli)` out of the per-parcel OMP cache build loop.
- Added `sampledBoundaryPatches_`, populated once in `createField()`.
- Skipped `processor`, `empty`, `symmetry`, and `wedge` patches during boundary
  measurement accumulation.
- Cached `cloud_.boundaryFluxMeasurements()` as a local const reference in the
  boundary accumulation loop.

The vibrational energy convention remains the current OFv1706 convention:
`constantProperties::eVib_m()`.  This intentionally does not copy the reference
tree's `vibLevel + 0.5` hot-path expression.

## Build

Build log:

`doc/worklog/v2506/detail/stage20_post_shared_cache_boundarymask_build_20260605.log`

Result: build passed.

## Smoke test

Smoke log:

`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage20_post_shared_cache_boundarymask_smoke_10step_20260605`

The smoke run temporarily used `ourmeshbkp/omp8/system/controlDict` with
`endTime 1.e-06`, then restored the original `ourmesh/omp8/system/controlDict`.

| metric | stage17 kept smoke | stage20 smoke | delta |
|---|---:|---:|---:|
| real | 6.72 s | 6.58 s | 0.14 s faster |
| full evolve wall | n/a | 1.519459887 s | n/a |
| move only | 0.826445053 s | 0.868335554 s | 0.041890501 s slower |
| buildCellOccupancy | 0.106092180 s | 0.123102951 s | 0.017010771 s slower |
| collision phase | 0.019019326 s | 0.020993557 s | 0.001974231 s slower |
| post fields/output | 0.632757959 s | 0.496770673 s | 0.135987286 s faster |

Smoke correctness tail:

```text
Collisions                = 2836
Collision candidates      = 2942
Number of DSMC particles  = 2065870
Number of stuck particles = 0
Number of free particles  = 2065870
Total energy              = 1.150626777
OpenMP enabled            = 1
OpenMP max threads        = 8
End main
```

## Formal OMP8 500-step validation

Formal log:

`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage20_post_shared_cache_boundarymask_bkpctrl_500step_20260605`

The formal run temporarily used `ourmeshbkp/omp8/system/controlDict`, then
restored the original `ourmesh/omp8/system/controlDict`.

Restored hash:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d
```

Formal profile versus stage17:

| metric | stage17 | stage20 | delta |
|---|---:|---:|---:|
| real | 106.84 s | 103.34 s | 3.50 s faster |
| full evolve wall | 110.2731630 s | 102.6520492 s | 7.6211138 s faster |
| move only | 53.31126355 s | 56.13548084 s | 2.82421729 s slower |
| buildCellOccupancy | 8.160171479 s | 9.115554564 s | 0.955383085 s slower |
| collision phase | 3.692066583 s | 3.760305061 s | 0.068238478 s slower |
| post fields/output | 44.5265183 s | 33.03582209 s | 11.49069621 s faster |

Formal profile versus historical `ourmeshbkp/omp8`:

| metric | ourmeshbkp/omp8 | stage20 | delta |
|---|---:|---:|---:|
| real | 103.66 s | 103.34 s | 0.32 s faster |
| full evolve wall | 102.1225033 s | 102.6520492 s | 0.5295459 s slower |
| move only | 48.82365912 s | 56.13548084 s | 7.31182172 s slower |
| buildCellOccupancy | 6.490896398 s | 9.115554564 s | 2.624658166 s slower |
| collision phase | 9.876225069 s | 3.760305061 s | 6.115920008 s faster |
| post fields/output | 36.92670834 s | 33.03582209 s | 3.89088625 s faster |

Formal correctness tail:

```text
Collisions                = 35624
Collision candidates      = 68756
Number of DSMC particles  = 2463612
Number of stuck particles = 0
Number of free particles  = 2463612
Average linear KE         = 9.805389794e-19
Average rotational energy = 2.423940447e-20
Average vibrational energy= 4.314963522e-21
Average electronic energy = 0
Total energy              = 1.24300724
Total Iterations          = 500
End main
```

## Conclusion

Stage20 is a successful post-fields/output stage.  The reference-style shared
sample cache is retained, and `post fields/output` is now faster than both the
stage17 kept implementation and the historical `ourmeshbkp/omp8` profile.

End-to-end `real` time is now slightly better than the historical bkp run, while
`full evolve wall` remains 0.53 s slower.  The remaining formal gap is no longer
post fields/output; it is mainly move and `buildCellOccupancy`.
