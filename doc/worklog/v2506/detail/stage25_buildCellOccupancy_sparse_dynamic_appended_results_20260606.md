# Stage25 BuildCellOccupancy Sparse/Dynamic/Appended Results - 2026-06-06

## Scope

This note records the direct implementation and validation of the
`buildCellOccupancy()` changes requested after the `zb-cylinder-react/omp8`
hotspot check:

- sparse reset of per-thread cell counts on the ordered OpenMP path;
- `occupancyOrderedParcels_` changed from `List` to `DynamicList` and resized
  without discarding capacity;
- ordered occupancy path now accepts `moveOrderedParcels_ + moveAppendedParcels_`
  when the move list does not already cover the full cloud.

Source files changed:

- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`

Build command:

```text
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
```

Build result: success.

`git diff --check` result: success.

## Implementation Notes

The OpenMP occupancy path now distinguishes three cases:

- move-ordered list already covers the full cloud;
- move-ordered list plus appended parcels covers the full cloud;
- fallback gather path over the full linked-list cloud.

The appended path is only used when `moveOrderedParcels_.size()` does not
already equal `this->size()`.  This avoids double-counting cases where the move
path has already folded appended parcels into the ordered list.

On the ordered path, per-thread `localCounts` are reset only for cells recorded
in the previous `occupancyThreadActiveCells_`.  The fallback gather path keeps
the full clear because its previous active-cell coverage is not guaranteed to
match the next full-cloud rebuild.

## `zb-cylinder-react/omp8` 300-Step Validation

Case:

- `run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/omp8`

Logs:

- before: `log.codex_omp8_zb_currentctrl_stage21_300step_20260606`
- after: `log.codex_omp8_zb_buildocc_sparse_dynamic_appended_300step_20260606`
- v2506 reference: `log.dsmcFoam+.20260605`

The case was run directly without initialization, using the existing v1706
compatible `writeControl runTime` controlDict.

### Current v1706 Before/After

| metric | stage21 before | stage25 after | delta | pct |
|---|---:|---:|---:|---:|
| real [s] | 80.79 | 76.97 | -3.82 | -4.728308% |
| full evolve wall [s] | 75.24167335 | 71.98270041 | -3.25897294 | -4.331340% |
| move only [s] | 41.62369781 | 40.34492359 | -1.27877422 | -3.072226% |
| buildCellOccupancy [s] | 8.394238035 | 7.902615635 | -0.49162240 | -5.856665% |
| collision phase [s] | 10.26836592 | 9.837595081 | -0.430770839 | -4.195126% |
| post fields/output [s] | 13.75492210 | 12.73351032 | -1.02141178 | -7.425791% |
| DSMC particles | 1958022 | 1958225 | +203 | +0.010368% |
| collisions | 242587 | 242699 | +112 | +0.046169% |
| collision candidates | 479075 | 479261 | +186 | +0.038825% |
| total energy | 0.001961480323 | 0.001962071283 | +5.9096e-07 | +0.030128% |

Correctness scan for the after log: no `FOAM FATAL`, `NaN`, segmentation
marker, `MPI_ABORT`, floating-point exception marker, or kill marker.

### Cross-Version Reference Check

This is not a strict regression comparison because the reference log is from
the v2506 code path.

| metric | v2506 reference | stage25 v1706 | delta | pct |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 74.71342017 | 71.98270041 | -2.73071976 | -3.654925% |
| move only [s] | 38.10127114 | 40.34492359 | +2.24365245 | +5.888655% |
| buildCellOccupancy [s] | 3.491106571 | 7.902615635 | +4.411509064 | +126.364205% |
| collision phase [s] | 20.61886309 | 9.837595081 | -10.781268009 | -52.288373% |
| post fields/output [s] | 12.49893663 | 12.73351032 | +0.23457369 | +1.876749% |
| DSMC particles | 1957421 | 1958225 | +804 | +0.041074% |
| collisions | 240977 | 242699 | +1722 | +0.714591% |
| collision candidates | 518328 | 479261 | -39067 | -7.537119% |

The v1706 stage25 total wall is close to or better than the v2506 reference in
this single log, but `buildCellOccupancy` remains much slower than the v2506
implementation.  The remaining gap should not be considered closed.

## `ourmesh/omp8` 500-Step Regression

Case:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8`

Valid OMP log:

- `log.codex_omp_stage25_buildocc_sparse_dynamic_appended_validomp_500step_20260606`

Comparison baseline:

- stage21 stability rerun mean from
  `doc/worklog/v2506/detail/stage21_stability_rerun_results_20260606`

The valid run used the OMP/profile controlDict from:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/omp8/system/controlDict`

ControlDict hashes:

```text
9e7d21d8c4bac968a736bc48b26c9fc377a7ac74bcaf5d155d1b07f06c62a2da  OMP/profile controlDict used for valid run
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d  restored ourmesh/omp8/system/controlDict after run
```

The pre-run `ff112...` file was preserved as:

- `system/controlDict.codex_bak_20260606_before_stage25_valid_omp`

### Valid OMP Result

The valid run reached `End main`, `Total Iterations = 500`, and reported:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

| metric | stage21 stability mean | stage25 valid OMP | delta | pct |
|---|---:|---:|---:|---:|
| real [s] | 98.38333333 | 90.34 | -8.04333333 | -8.175504% |
| full evolve wall [s] | 98.16965099 | 90.00671836 | -8.16293263 | -8.315129% |
| move only [s] | 53.73855140 | 49.24385297 | -4.49469843 | -8.364011% |
| buildCellOccupancy [s] | 8.408399703 | 6.447003770 | -1.961395933 | -23.326626% |
| collision phase [s] | 3.608682147 | 3.376808391 | -0.231873756 | -6.425441% |
| post fields/output [s] | 31.83707414 | 30.39917036 | -1.43790378 | -4.516444% |
| DSMC particles | 2463658.67 | 2463726 | +67.33 | +0.002733% |
| collisions | 35717.33 | 35977 | +259.67 | +0.727005% |
| collision candidates | 68942.00 | 69118 | +176 | +0.255287% |
| total energy | 1.243132700 | 1.243107069 | -0.000025631 | -0.002062% |

Correctness scan for the valid OMP log: no `FOAM FATAL`, `NaN`, segmentation
marker, `MPI_ABORT`, floating-point exception marker, or kill marker.

The stage25 run improves the main `ourmesh/omp8` stability baseline by about
8.3% in `full evolve wall`.  The largest direct improvement is
`buildCellOccupancy`, down 23.3% versus the stage21 stability mean.

### Invalid No-OMP Run

One preliminary run used the currently restored `ff112...` controlDict, which
does not contain `useOpenMP true` or the profile/OpenMP tuning entries.  That
log is retained but excluded from the performance conclusion:

- `log.codex_omp_stage25_buildocc_sparse_dynamic_appended_bkpctrl_500step_20260606`

It reached 500 steps but reported:

- `OpenMP enabled = 0`
- `real = 365.19 s`
- `full evolve wall = 385.4477822 s`
- `move only = 261.5343933 s`
- `buildCellOccupancy = 27.28651133 s`

This log is useful only as a configuration caution.  It is not an OMP
regression result.

## Conclusion

Stage25 is a positive change on both tested surfaces:

- `zb-cylinder-react/omp8` current-control 300-step run: `full evolve wall`
  improved by 4.33%, and `buildCellOccupancy` improved by 5.86%;
- `ourmesh/omp8` valid OMP 500-step run: `full evolve wall` improved by 8.32%
  versus the stage21 stability mean, and `buildCellOccupancy` improved by
  23.33%.

The code compiles cleanly, passes whitespace checks, and both valid test logs
finish with stable particle, collision, candidate, and energy metrics.  The
stage25 source changes should be kept as the current best OMP occupancy
implementation.

Remaining issue: on `zb-cylinder-react/omp8`, v1706 stage25
`buildCellOccupancy` is still much slower than the v2506 reference
(`7.902615635 s` versus `3.491106571 s`).  That remaining gap needs a separate
subphase-level profile before more occupancy changes are made.

## Performance Source Analysis

This section separates three comparison targets:

- reference baseline logs;
- previous best stage21 logs;
- the direct source-level effect of stage25.

### `zb-cylinder-react/omp8`

Compared with the v2506 reference log, stage25 v1706 is faster in total wall
only because the collision phase is much lower:

| component | v2506 reference | stage25 v1706 | delta | share of full-wall gain |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 74.71342017 | 71.98270041 | -2.73071976 | 100.000% |
| move only [s] | 38.10127114 | 40.34492359 | +2.24365245 | -82.163% |
| buildCellOccupancy [s] | 3.491106571 | 7.902615635 | +4.411509064 | -161.551% |
| collision phase [s] | 20.61886309 | 9.837595081 | -10.781268009 | +394.814% |
| post fields/output [s] | 12.49893663 | 12.73351032 | +0.23457369 | -8.590% |

Conclusion for the reference comparison: the v1706 stage25 full-wall advantage
is not caused by `buildCellOccupancy`.  `buildCellOccupancy` remains the main
deficit against the v2506 reference.  The apparent total advantage is dominated
by collision phase differences, with different candidate counts
(`518328 -> 479261`) also contributing to the collision-time gap.

Compared with the current v1706 stage21 log, the stage25 gain is distributed:

| component | stage21 | stage25 | delta | share of full-wall gain |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 75.24167335 | 71.98270041 | -3.25897294 | 100.000% |
| move only [s] | 41.62369781 | 40.34492359 | -1.27877422 | 39.239% |
| buildCellOccupancy [s] | 8.394238035 | 7.902615635 | -0.491622400 | 15.085% |
| collision phase [s] | 10.26836592 | 9.837595081 | -0.430770839 | 13.218% |
| post fields/output [s] | 13.75492210 | 12.73351032 | -1.021411780 | 31.342% |

Conclusion for the previous-best comparison: the directly attributable stage25
source change explains the `buildCellOccupancy` drop.  The additional drops in
move, collision, and post fields are measured in the same run but should be
treated as run-to-run variation and possible secondary cache/allocation effects,
not as direct evidence that stage25 changed those kernels.

### `ourmesh/omp8`

Compared with the preserved `ourmeshbkp/omp8` reference log, stage25 is faster
mostly because collision and post-field/output are already much lower:

| component | ourmeshbkp reference | stage25 valid OMP | delta | share of full-wall gain |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 102.1225033 | 90.00671836 | -12.11578494 | 100.000% |
| move only [s] | 48.82365912 | 49.24385297 | +0.420193850 | -3.468% |
| buildCellOccupancy [s] | 6.490896398 | 6.447003770 | -0.043892628 | 0.362% |
| collision phase [s] | 9.876225069 | 3.376808391 | -6.499416678 | 53.644% |
| post fields/output [s] | 36.92670834 | 30.39917036 | -6.527537980 | 53.876% |

Conclusion for the reference comparison: relative to the preserved baseline,
stage25's total speedup is not mainly from this occupancy patch.  The large
speedup against the reference is mostly inherited from earlier collision and
post-field/output improvements.  Stage25 brings `buildCellOccupancy` roughly to
reference parity (`6.447003770 s` versus `6.490896398 s`).

Compared with the previous best single stage21 formal run:

| component | stage21 best | stage25 valid OMP | delta | share of full-wall gain |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 93.28696662 | 90.00671836 | -3.28024826 | 100.000% |
| move only [s] | 51.17943521 | 49.24385297 | -1.935582240 | 59.007% |
| buildCellOccupancy [s] | 8.087915571 | 6.447003770 | -1.640911801 | 50.024% |
| collision phase [s] | 3.314541395 | 3.376808391 | +0.062266996 | -1.898% |
| post fields/output [s] | 30.15687469 | 30.39917036 | +0.242295670 | -7.387% |

Compared with the stage21 stability mean:

| component | stage21 stability mean | stage25 valid OMP | delta | share of full-wall gain |
|---|---:|---:|---:|---:|
| full evolve wall [s] | 98.16965099 | 90.00671836 | -8.16293263 | 100.000% |
| move only [s] | 53.73855140 | 49.24385297 | -4.49469843 | 55.062% |
| buildCellOccupancy [s] | 8.408399703 | 6.447003770 | -1.961395933 | 24.028% |
| collision phase [s] | 3.608682147 | 3.376808391 | -0.231873756 | 2.841% |
| post fields/output [s] | 31.83707414 | 30.39917036 | -1.437903780 | 17.615% |

Conclusion for the previous-best comparison: the robust, directly explained
stage25 gain is the `buildCellOccupancy` drop (`8.09/8.41 s -> 6.45 s`).  The
move-time drop is visible in the logs but is not a direct target of this patch,
so it should be treated as run-to-run variation or secondary memory/cache
behavior unless repeated runs confirm it.

### Final Attribution

The source of the performance gain depends on the baseline:

- against preserved/reference baselines, most total speedup comes from earlier
  collision and post-field/output improvements, not from stage25;
- against the previous stage21 best, stage25's directly attributable source is
  `buildCellOccupancy`, especially on `ourmesh/omp8`;
- on `zb-cylinder-react/omp8`, stage25 improves `buildCellOccupancy`, but the
  total-wall gain is spread across several phases and the v2506 reference still
  has a much faster occupancy path;
- on `ourmesh/omp8`, stage25 closes the occupancy gap relative to the preserved
  reference while keeping the earlier collision/post-field gains.
