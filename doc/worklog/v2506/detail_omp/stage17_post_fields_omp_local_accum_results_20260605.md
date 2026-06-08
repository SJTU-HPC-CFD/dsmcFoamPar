# Stage 17: post fields OMP local accumulation

Date: 2026-06-05

## Scope

Current kept baseline before this stage: stage15, uniform-dt move fast path.

This stage starts the `post fields/output` optimization line.  It does not
change move, collision, or the `buildCellOccupancy()` algorithm.

Changed files:

- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H`
- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`

Main changes:

- Cache `openmpFieldSampling` instead of reading `controlDict` in every
  `calculateField()` call.
- Build a `typeId -> species index` map once in `createField()` and use it in
  parcel sampling instead of inner-loop `findIndex(speciesIds_, p.typeId())`.
- Keep the OpenMP field-sampling path over the flat occupancy view, but skip
  empty cells before allocating local accumulators.
- Accumulate core per-cell quantities into local scalars/vectors and write them
  back once per cell.
- Skip heat-flux/shear-stress high-order parcel accumulators when
  `measureHeatFluxShearStress` is false.  The formal case only enables
  `measureMeanFreePath true`.

## Build

Build logs:

- `doc/worklog/v2506/detail/stage17_post_fields_omp_local_accum_build_20260605.log`
- `doc/worklog/v2506/detail/stage17_post_fields_omp_local_accum_build2_20260605.log`

Both builds completed successfully.  `build2` is the kept revision with the
empty-cell fast skip.

`git diff --check` passed for the touched `dsmcVolFields` files.

## 10-step smoke

Case:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8
```

Smoke used `ourmeshbkp/omp8/system/controlDict` with `endTime 1.e-06`.
The original `ourmesh/omp8/system/controlDict` was restored after the run.

Logs:

- first attempt: `log.codex_omp_stage17_post_fields_omp_local_accum_smoke_10step_20260605`
- kept revision: `log.codex_omp_stage17_post_fields_omp_local_accum_v2_smoke_10step_20260605`

The first attempt initialized local accumulators for every cell, including
empty cells, and was rejected:

| run | real | move only | buildCellOccupancy | collision | post fields/output |
|---|---:|---:|---:|---:|---:|
| stage15 smoke | 6.94 s | 0.798181885 s | 0.104578723 s | 0.019113147 s | 0.673572436 s |
| stage17 attempt 1 | 7.25 s | 0.918349019 s | 0.119848878 s | 0.020275160 s | 0.692359891 s |
| stage17 kept | 6.72 s | 0.826445053 s | 0.106092180 s | 0.019019326 s | 0.632757959 s |

Kept smoke vs stage15:

- real: `6.94 -> 6.72 s`, `0.22 s` faster.
- post fields/output: `0.673572436 -> 0.632757959 s`, `0.040814477 s`
  faster, about `6.1%`.

Smoke correctness tail for kept run:

```text
Collisions                = 2882
Collision candidates      = 2998
Number of DSMC particles  = 2065874
Total energy              = 1.150627055
End main
```

## Formal 500-step result

Formal run:

```text
OMP_NUM_THREADS=8 dsmcFoam+
```

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage17_post_fields_omp_local_accum_bkpctrl_500step_20260605
```

The run temporarily used `ourmeshbkp/omp8/system/controlDict`, then restored the
original `ourmesh/omp8/system/controlDict`.

Restored hash:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d
```

Formal profile:

| metric | stage15 | stage17 | delta |
|---|---:|---:|---:|
| real | 114.53 s | 106.84 s | 7.69 s faster |
| full evolve wall | 116.6882995 s | 110.2731630 s | 6.4151365 s faster |
| move only | 52.99979159 s | 53.31126355 s | 0.31147196 s slower |
| buildCellOccupancy | 8.663262117 s | 8.160171479 s | 0.503090638 s faster |
| collision phase | 3.628983484 s | 3.692066583 s | 0.063083099 s slower |
| post fields/output | 50.8090158 s | 44.5265183 s | 6.2824975 s faster |

Formal correctness tail:

```text
Collisions                = 35552
Collision candidates      = 68374
Number of DSMC particles  = 2463676
Average linear KE         = 9.806150241e-19
Average rotational energy = 2.420874363e-20
Average vibrational energy= 4.333286131e-21
Total energy              = 1.243118007
End main
```

## Baseline comparison

Historical `ourmeshbkp/omp8` profile:

| metric | ourmeshbkp/omp8 | stage17 | remaining delta |
|---|---:|---:|---:|
| real | 103.66 s | 106.84 s | 3.18 s slower |
| move only | 48.82365912 s | 53.31126355 s | 4.48760443 s slower |
| buildCellOccupancy | 6.490896398 s | 8.160171479 s | 1.669275081 s slower |
| collision phase | 9.876225069 s | 3.692066583 s | 6.184158486 s faster |
| post fields/output | 36.92670834 s | 44.5265183 s | 7.59980996 s slower |

Stage17 recovers most of the stage15 end-to-end gap, but post fields/output is
still the largest remaining gap against `ourmeshbkp/omp8`.  The next post step
should target the remaining per-species/vibrational sampling writes or port a
minimal shared sample cache only after profiling confirms the cost.

`buildCellOccupancy()` remains about `1.67 s` slower than historical baseline.
This stage did not structurally change it; a separate stage should inspect
whether the remaining cost is full-cell scans, count-list reset, prefix build,
or fill.
