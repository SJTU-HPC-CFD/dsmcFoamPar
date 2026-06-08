# MPI replicated mesh strict full-fields post/output structure audit - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI8 replicated mesh DLB only, `useOpenMP false`,
  `openmpThreads 1`, `mpirun -np 8`.
- Source/case action in this audit: structural read-only analysis plus this
  worklog note.  No source, binary, or case configuration was changed.

Current restored formal hashes:

```text
controlDict         8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
fieldPropertiesDict 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
```

The current formal case is a strict full-fields run:

```text
endTime       5.e-05;
deltaT        1.e-07;
writeInterval 1.e-3;
profileSummary true;
profileDetail  false;
dsmcFields: one dsmcVolFields mixture field, typeIds (N2 O2 NO N O)
```

## Baseline Numbers

Strict same-tree full-fields baseline:

```text
log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607
```

| metric | value |
| --- | ---: |
| real | 133.83 s |
| move only | 72.9567669 s |
| buildCellOccupancy | 6.729516136 s |
| collision phase | 14.8398869 s |
| post fields/output | 40.37082178 s |
| full evolve wall | 132.7056254 s |
| OpenMP enabled | 0 |
| DLB rebalances | 8 |
| Total Iterations | 500 |

Preserved historical reference, kept separate:

```text
ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610
```

| metric | value |
| --- | ---: |
| external wall | 117.31 s |
| move only | 47.1652636 s |
| post fields/output | 30.35563996 s |
| full evolve wall | 114.9207628 s |
| DLB rebalances | 8 |

No-fields diagnostic mode, not a full-fields baseline:

```text
log.codex_mpi8_replicatedmesh_no_fields_500step_20260607
```

| metric | strict full fields | no fields | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 119.69 s | -14.14 s |
| move only | 72.9567669 s | 64.72452806 s | -8.23223884 s |
| post fields/output | 40.37082178 s | 30.60678571 s | -9.76403607 s |
| full evolve wall | 132.7056254 s | 118.6090831 s | -14.09654230 s |

Interpretation: no-fields proves that the `dsmcVolFields` entry accounts for
roughly a 10 s same-tree post delta over 500 steps, but it also changes runtime
state and must not be mixed into strict full-fields comparisons.

## What The Timer Actually Covers

The solver loop calls `dsmc.evolve()` first, then handles replicated-mesh output
and `runTime.write()` outside the cloud post timer
(`applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C:146-205`).

The printed `post fields/output [s]` is accumulated in
`dsmcCloud::evolve_fields()`:

```text
src/lagrangian/dsmc/clouds/dsmcCloud.C:1808-1837
```

It includes all of these calls:

- `reactions_.outputData()`
- `fields_.calculateFields()`
- `fields_.writeFields()`
- `controllers_.calculateProps()`
- `controllers_.outputResults()`
- `boundaries_.calculateProps()`
- `boundaries_.outputResults()`
- `boundaryMeas_.outputResults()`
- `trackingInfo_.clean()`
- `boundaryMeas_.clean()`
- `cellMeas_.clean()`

Therefore the timer name is broader than `dsmcVolFields` and broader than disk
output.  It is the whole end-of-step field/control/boundary/measurement cleanup
section.

## Current Case Objects

Current dictionaries reduce the active path:

- `system/chemReactDict` has `reactions ();`, so `reactions_.outputData()`
  only increments its counter and loops over an empty list.
- `system/controllersDict` has no state or flux controllers, so controller
  loops are empty.
- `system/fieldPropertiesDict` has one `dsmcVolFields` mixture entry with
  `measureMeanFreePath true`.
- `system/boundariesDict` has one `dsmcFreeStreamInflowPatch`, two
  `dsmcDeletionPatch` entries, and one `dsmcDiffuseWallPatch`.

The active non-field post floor is therefore mostly boundary model property
calculation plus measurement cleanup, not reactions or controllers.

## Per-Step Versus Output-Time Work

Per-step work in the current formal run:

- `dsmcFieldProperties::calculateFields()` calls every field's
  `calculateField()` every step
  (`src/lagrangian/dsmc/macroscopicProperties/basic/dsmcFieldProperties/dsmcFieldProperties.C:150-157`).
- `dsmcVolFields::calculateField()` samples owned cells in replicated mesh mode
  and accumulates boundary measurements
  (`src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:2417-2999`).
- `dsmcBoundaries::calculateProps()` loops patch, cyclic, and general boundary
  models every step
  (`src/lagrangian/dsmc/boundaries/basic/dsmcBoundaries/dsmcBoundaries.C:570-585`).
- `trackingInfo_.clean()` clears per-type face flux arrays
  (`src/lagrangian/dsmc/faceTracker/dsmcFaceTracker.C:86-94`).
- `boundaryMeas_.clean()` clears species and patch boundary measurement arrays
  (`src/lagrangian/dsmc/boundaryMeasurements/boundaryMeasurements.C:364-385`).
- `cellMeas_.clean()` clears cell collision measurement arrays
  (`src/lagrangian/dsmc/cellMeasurements/cellMeasurements.C:82-88`).

Output-time-only work:

- `dsmcFieldProperties::writeFields()` only creates the `uniform` directory and
  re-reads field dictionaries under `runTime.outputTime()`, although it still
  calls each field's `writeField()` every step
  (`dsmcFieldProperties.C:172-195`).
- `dsmcVolFields::writeField()` is empty in the current tree
  (`dsmcVolFields.C:4330-4331`).
- `dsmcVolFields::calculateField()` has a large `outputTime()` block that
  performs replicated reductions, rank0 field computation, field writes, and
  optional reset
  (`dsmcVolFields.C:3001-4083`).
- Controller, boundary, and boundary-measurement `outputResults()` methods are
  guarded by `runTime.outputTime()`.

For the formal 500-step controls, `writeInterval 1.e-3` is larger than
`endTime 5.e-05`.  The case directory has no `5e-05` output directory and the
formal log has no replicated output timing block.  Thus the formal 40 s post
timer is not disk field output.

## Profile Detail Trap

Current `dsmcVolFields` declares and prints these fields:

```text
profileOutputComputeWallTime_
profileFieldWriteWallTime_
profileOutputResetWallTime_
profileOutputTimeWallTime_
```

But in the current source they are only declared, initialized, and printed.
There is no `+=` accumulation for them.  The populated detail timers are the
sampling-body timers:

- `profileSampleAccumWallTime_`
- `profileSharedCacheBuildWallTime_`
- `profileFieldCombineWallTime_`
- `profileCellReduceWallTime_`
- `profileBoundaryAccumWallTime_`

So a detail line of `output compute = 0`, `field writes = 0`, or
`output-time block = 0` should not be interpreted as a measured breakdown of the
output branch.  In this formal case the output branch is also not reached, but
that conclusion comes from `writeInterval/endTime`, case directories, and log
absence, not from those currently-unwired detail counters.

## Existing Detail Evidence

10-step post detail from the strict full-field path showed:

```text
log.codex_mpi8_replicatedmesh_post_dsmcN_owned_reset_smoke_10step_20260607
post fields/output                 0.865637994 s
dsmcVolFields sample accumulation  0.239334675 s
shared cache build                 0.167777043 s
field combine                      0.071534244 s
cell reduction                     0.005973127 s
boundary accumulation              0.003791123 s
```

Another 10-step direct-sampling attempt removed cache/combine internally but did
not improve the total post timer:

| metric | strict detail smoke | direct sampling smoke | delta |
| --- | ---: | ---: | ---: |
| dsmcVolFields sample accumulation | 0.21323199 s | 0.199297762 s | -0.013934228 s |
| dsmcVolFields field combine | 0.061793793 s | 0 s | -0.061793793 s |
| dsmcVolFields shared cache build | 0.151396085 s | 0 s | -0.151396085 s |
| dsmcVolFields direct parcel accum | 0 s | 0.199266246 s | +0.199266246 s |
| post fields/output | 0.834544903 s | 0.836535452 s | +0.001990549 s |

This explains why the post path should not be optimized by only rearranging the
existing `dsmcVolFields` cache/combine structure.

## Conclusions

1. `post fields/output` is an `evolve_fields()` aggregate, not just
   `dsmcVolFields` and not just disk output.
2. Formal strict full-fields post cost is a per-step cost.  With
   `writeInterval 1.e-3 > endTime 5.e-05`, the expensive disk-output branch is
   not the current 500-step bottleneck.
3. `dsmcVolFields` contributes materially: no-fields mode removes about
   `9.76 s` from the 500-step post timer.  But no-fields is a different
   operational mode, not the strict baseline.
4. A large non-field floor remains: no-fields still reports
   `post fields/output 30.60678571 s`, close to the preserved reference
   `30.35563996 s`.
5. Current profile detail is insufficient to split that non-field floor because
   `evolve_fields()` itself has only one aggregate timer.

## Next Valid Step

Before another optimization candidate, add a diagnostic-only sub-timer around
the individual `evolve_fields()` calls:

- reactions output;
- field calculate;
- field write wrapper;
- controller calculate/output;
- boundary calculate/output;
- boundary measurement output;
- face-tracker clean;
- boundary-measurement clean;
- cell-measurement clean.

Run it as pure MPI8 with `profileDetail true` on a 10-step smoke first.  If the
split is stable, promote to a 200-step signal run.  Keep the results separated
from the formal strict full-fields baseline until the instrumentation is removed
and the binary is rebuilt.
