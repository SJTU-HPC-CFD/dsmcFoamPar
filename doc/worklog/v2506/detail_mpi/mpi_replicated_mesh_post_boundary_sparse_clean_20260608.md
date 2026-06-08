# MPI replicated mesh post boundary sparse clean - 2026-06-08

## Scope

- Mode: pure MPI8 replicated mesh DLB only, no OpenMP.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Strict full-fields formal controls:
  `endTime 5.e-05`, `profileSummary true`, `profileDetail false`,
  one `dsmcVolFields` mixture entry in `fieldPropertiesDict`.
- Current restored case hashes:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

## Starting Point

Strict same-tree full-fields baseline:

```text
log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607
```

| metric | strict baseline |
| --- | ---: |
| real | 133.83 s |
| move only | 72.9567669 s |
| buildCellOccupancy | 6.729516136 s |
| collision phase | 14.8398869 s |
| post fields/output | 40.37082178 s |
| full evolve wall | 132.7056254 s |
| DLB rebalances | 8 |

This is separate from the historical `ourmeshbkp` reference
(`External wall seconds 117.31`, `move only 47.1652636`) and from the
no-fields optional mode (`real 119.69`, `post 30.60678571`).

## Diagnostic Split

Temporary source instrumentation was added only to split
`dsmcCloud::evolve_fields()` into sub-timers.  It was later removed before the
formal candidate runs.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_evolve_fields_subprofile_build_20260608.log
```

10-step diagnostic:

```text
log.codex_mpi8_replicatedmesh_post_evolve_fields_subprofile_smoke_10step_20260608
```

| metric | value |
| --- | ---: |
| real | 12.54 s |
| post fields/output | 0.9103465 s |
| fields calculate | 0.336047609 s |
| tracking clean | 0.098220448 s |
| boundary measurement clean | 0.514156583 s |
| cell measurement clean | 0.005162174 s |

200-step diagnostic:

```text
log.codex_mpi8_replicatedmesh_post_evolve_fields_subprofile_signal_200step_20260608
```

| metric | value |
| --- | ---: |
| real | 56.06 s |
| move only | 28.060843 s |
| post fields/output | 17.06698454 s |
| fields calculate | 8.718456902 s |
| tracking clean | 1.46317059 s |
| boundary measurement clean | 8.409648626 s |
| cell measurement clean | 0.071565471 s |

Interpretation: `boundaryMeas_.clean()` was the largest removable post sub-cost.
It was full-clearing all species/patch/face boundary-flux arrays every step.

## Candidate

Source change:

- files:
  - `src/lagrangian/dsmc/boundaryMeasurements/boundaryMeasurements.H`
  - `src/lagrangian/dsmc/boundaryMeasurements/boundaryMeasurementsI.H`
  - `src/lagrangian/dsmc/boundaryMeasurements/boundaryMeasurements.C`
- add per-step sparse touched tracking for boundary flux entries:
  `(species, patch, face)`;
- mark touched entries through the non-const boundary-flux accessors;
- change `boundaryMeasurements::clean()` to reset only touched entries, while
  still clearing all flux components for each touched face;
- keep `nParticlesOnStickingBoundaries_` full clear unchanged.

The candidate preserves the existing order:

1. wall boundary model writes `boundaryFluxMeasurements()`;
2. `dsmcVolFields::calculateField()` consumes those values;
3. `boundaryMeas_.clean()` resets values for the next step.

Build logs:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_boundary_sparse_clean_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_post_boundary_sparse_clean_formal_build_20260608.log
```

Both builds passed, with only the existing OpenFOAM template warnings.

## 10-step Smoke

Temporary controls:

```text
endTime 1.e-06;
profileSummary true;
profileDetail true;
```

Log:

```text
log.codex_mpi8_replicatedmesh_post_boundary_sparse_clean_smoke_10step_20260608
```

| metric | diagnostic baseline | sparse clean |
| --- | ---: | ---: |
| real | 12.54 s | 8.51 s |
| post fields/output | 0.9103465 s | 0.385807841 s |
| fields calculate | 0.336047609 s | 0.339796802 s |
| boundary measurement clean | 0.514156583 s | 0.003575961 s |
| tracking clean | 0.098220448 s | 0.041921618 s |

Correctness:

```text
Total Iterations = 10
Number of DSMC particles = 2065840
Number of stuck particles = 0
Collisions = 2907
Total energy = 1.150617828
OpenMP enabled = 0
```

No `Fatal`, `NaN`, or `BAD TERMINATION` string was found.

## 200-step Signal

Temporary controls:

```text
endTime 2.e-05;
profileSummary true;
profileDetail true;
```

Log:

```text
log.codex_mpi8_replicatedmesh_post_boundary_sparse_clean_signal_200step_20260608
```

| metric | diagnostic baseline | sparse clean |
| --- | ---: | ---: |
| real | 56.06 s | 44.20 s |
| move only | 28.060843 s | 25.15106083 s |
| post fields/output | 17.06698454 s | 6.963116119 s |
| fields calculate | 8.718456902 s | 6.346710544 s |
| boundary measurement clean | 8.409648626 s | 0.037444006 s |
| full evolve wall | 51.42166433 s | 39.20235374 s |
| DLB rebalances | 2 | 2 |

Correctness:

```text
Total Iterations = 200
Number of DSMC particles = 2220224
Number of stuck particles = 0
Collisions = 15114
Total energy = 1.169532893
OpenMP enabled = 0
```

No `Fatal`, `NaN`, or `BAD TERMINATION` string was found.

## 500-step Formal

The diagnostic `dsmcCloud::evolve_fields()` sub-timers were removed before
these formal runs.  Only the sparse boundary clean source change was retained.
The strict full-fields case hashes were restored.

Formal log:

```text
log.codex_mpi8_replicatedmesh_post_boundary_sparse_clean_formal_500step_20260608
```

| metric | strict baseline | sparse clean | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 121.40 s | -12.43 s |
| move only | 72.9567669 s | 80.52274982 s | +7.56598292 s |
| buildCellOccupancy | 6.729516136 s | 7.583478988 s | +0.853962852 s |
| collision phase | 14.8398869 s | 15.63210682 s | +0.79221992 s |
| post fields/output | 40.37082178 s | 19.56710696 s | -20.80371482 s |
| full evolve wall | 132.7056254 s | 121.1536038 s | -11.5520216 s |
| DLB rebalances | 8 | 8 | 0 |

Correctness:

```text
Total Iterations = 500
Number of DSMC particles = 2463879
Number of stuck particles = 0
Collisions = 34493
Total energy = 1.238780249
OpenMP enabled = 0
```

No `Fatal`, `NaN`, or `BAD TERMINATION` string was found.

Formal repeat log:

```text
log.codex_mpi8_replicatedmesh_post_boundary_sparse_clean_formal_repeat_500step_20260608
```

| metric | strict baseline | sparse clean repeat | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 113.24 s | -20.59 s |
| move only | 72.9567669 s | 70.66798424 s | -2.28878266 s |
| buildCellOccupancy | 6.729516136 s | 7.264959245 s | +0.535443109 s |
| collision phase | 14.8398869 s | 14.50986947 s | -0.33001743 s |
| post fields/output | 40.37082178 s | 18.95231379 s | -21.41850799 s |
| full evolve wall | 132.7056254 s | 110.5389487 s | -22.1666767 s |
| DLB rebalances | 8 | 8 | 0 |

Correctness:

```text
Total Iterations = 500
Number of DSMC particles = 2463782
Number of stuck particles = 0
Collisions = 34583
Total energy = 1.238881763
OpenMP enabled = 0
```

No `Fatal`, `NaN`, or `BAD TERMINATION` string was found.

## Decision

Retain the sparse boundary clean candidate.

Reason:

- Both 500-step strict full-fields formal runs beat the same-tree baseline in
  end-to-end wall time.
- The post aggregate drops from `40.37082178 s` to about `19 s`.
- The first formal run showed extra move cost, but the repeat did not; this
  points to run-to-run/DLB balance variability rather than a stable move
  regression large enough to reject the candidate.
- The candidate does not disable field sampling and does not change the
  no-fields/full-fields contract.

Current accepted source delta is intentionally limited to
`boundaryMeasurements.{H,I.H,C}`.  The temporary `dsmcCloud` post subprofile
instrumentation was removed before formal validation.

## Remaining Space

The repeat formal result is already near the historical reference wall time,
but the source state and contract are different:

- current strict full-fields sparse-clean repeat: `real 113.24 s`,
  `move only 70.66798424 s`, `post 18.95231379 s`;
- historical `ourmeshbkp` reference: `External wall seconds 117.31`,
  `move only 47.1652636 s`, `post 30.35563996 s`.

The sparse clean candidate fixed the current-tree post floor; it did not solve
the remaining move gap to the historical reference.  Further work should not
return to the rejected tracking micro-candidates unless new evidence changes
the cost structure.

## Rejected Follow-up: last touched fast path

After accepting sparse clean, one narrow move-side follow-up was tested:

- add a one-entry cache to `boundaryMeasurements::markTouched()`;
- return early when consecutive boundary-flux accessor writes hit the same
  `(species, patch, face)` and `cleanStamp`;
- goal: reduce repeated stamp-array lookups caused by many boundary-flux
  component writes during one wall hit.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_boundary_sparse_clean_lasttouch_build_20260608.log
```

10-step smoke:

```text
log.codex_mpi8_replicatedmesh_post_boundary_sparse_clean_lasttouch_smoke_10step_20260608
```

| metric | value |
| --- | ---: |
| real | 12.30 s |
| move only | 1.114251679 s |
| post fields/output | 0.375942484 s |
| full evolve wall | 2.151170739 s |
| final particles | 2065850 |
| final collisions | 2890 |
| total energy | 1.150622528 |

200-step signal:

```text
log.codex_mpi8_replicatedmesh_post_boundary_sparse_clean_lasttouch_signal_200step_20260608
```

| metric | sparse clean signal | last-touch signal |
| --- | ---: | ---: |
| real | 44.20 s | 46.51 s |
| move only | 25.15106083 s | 25.46408405 s |
| post fields/output | 6.963116119 s | 7.051579863 s |
| full evolve wall | 39.20235374 s | 39.410004 s |
| DLB rebalances | 2 | 2 |

Decision: reject and revert.  The signal did not improve move or end-to-end
time, and the last-touch run used `profileDetail false` while the previous
sparse-clean signal still included detail profiling overhead.  This is not a
credible positive candidate.

Revert build:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_boundary_sparse_clean_lasttouch_revert_build_20260608.log
```

Post-revert state:

```text
system/controlDict hash restored to 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
system/fieldPropertiesDict hash remains 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
no lastTouched / last-touch / postDetail / post detail symbols remain in touched source files
```
