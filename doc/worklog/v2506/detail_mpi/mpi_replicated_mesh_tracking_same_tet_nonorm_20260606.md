# MPI replicated mesh same-tet non-normalised inside check - 2026-06-06

## Scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI, `useOpenMP false`, `openmpThreads 1`
- Retained DLB configuration:
  `replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);`
- Baseline for this candidate: forced8 500-step run,
  `log.codex_mpi8_replicatedmesh_mpionly_forced_dlb_steps_500step_20260606`

Acceptance gate: improve both forced8 `real 132.29 s` and
`move only 69.49945105 s` in 500-step pure-MPI validation.

## Source candidate

File:

```text
src/lagrangian/basic/particle/particleTemplates.C
```

Change:

- Added `dsmcTrackInsideTetNoNormalise()`.
- Replaced the DSMC static-mesh same-tet fast-path call from
  `tet.inside(endPosition)` to `dsmcTrackInsideTetNoNormalise(tet, endPosition)`.

Reasoning:

- The standard `tet.inside()` implementation normalises all four tet face
  normals before the half-space checks.
- For this MPI replicated-mesh move path, most `trackToFace(..., true)` calls do
  not cross a face.
- The replacement keeps the same tolerance relation:
  `(pt-base) & (n/(mag(n)+VSMALL)) > SMALL`
  is evaluated as
  `(pt-base) & n > SMALL*(mag(n)+VSMALL)`.
- The common inside path short-circuits the `mag(n)` work when the signed
  distance is non-positive.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_same_tet_nonorm_build_20260606.log
```

Build result: passed.  Only the existing OpenFOAM template-instantiation
warnings were emitted.

## 10-step diagnostic smoke

Temporary controls:

```text
endTime 1.e-06;
profileDetail true;
moveDetailProfile true;
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_mpionly_same_tet_nonorm_smoke_10step_20260606
```

Comparison to the previous move-detail diagnostic:

| metric | previous diagnostic | same-tet non-normalised |
| --- | ---: | ---: |
| real | 10.27 s | 10.26 s |
| move only | 1.17689206 s | 1.145794634 s |
| move detail track calls | 25,758,269 | 25,758,353 |
| move detail face hits | 5,135,781 | 5,135,912 |
| move detail track max | 0.719672337 s | 0.699361356 s |
| move detail boundary max | 0.001153038 s | 0.001129234 s |
| OpenMP enabled | 0 | 0 |
| Total Iterations | 10 | 10 |

Smoke decision: proceed to 500-step validation.  The move and track timer deltas
match the intended same-tet fast-path mechanism.

## 500-step validation

First formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_same_tet_nonorm_500step_20260606
```

Repeat formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_same_tet_nonorm_confirm_500step_20260606
```

Results:

| metric | forced8 baseline | same-tet run 1 | same-tet confirm |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 127.91 s | 129.93 s |
| move only | 69.49945105 s | 67.56621687 s | 69.11276280 s |
| full evolve wall | 123.0267073 s | 118.7495988 s | 120.5815986 s |
| buildCellOccupancy | 6.216832977 s | 6.460351318 s | 6.357288369 s |
| collision phase | 12.98068576 s | 13.32302832 s | 13.64615726 s |
| post fields/output | 36.15433245 s | 36.72935704 s | 36.30228583 s |
| DLB checks | 500 | 500 | 500 |
| DLB rebalances | 8 | 8 | 8 |

Correctness checks from both formal logs:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- No Fatal/NaN match in the logs.
- First run final state: particles `2463800`, stuck `0`, final total energy
  `1.239005282`.
- Confirm run final state: particles `2463823`, stuck `0`, final total energy
  `1.238912214`.

## Decision

Retain the source candidate.

Both 500-step runs beat the forced8 baseline on the required metrics:

- First run: `real -4.38 s`, `move -1.93323418 s`.
- Confirm run: `real -2.36 s`, `move -0.38668825 s`.

The repeat margin on move is small, but it is still positive and the smoke
detail timer also moved in the expected direction.  Treat the validated range
for the retained candidate as:

```text
real 127.91-129.93 s
move only 67.56621687-69.11276280 s
```

Remaining gap to the preserved `ourmeshbkp` reference is still large:

- Versus `ourmeshbkp real 117.31 s`: at least `+10.60 s`.
- Versus `ourmeshbkp move 47.1652636 s`: at least `+20.40095327 s`.

Next work should continue in the old tracking core, but this candidate shows
that shaving same-tet geometry overhead is measurable.  Further changes should
be validated with the same repeat discipline when the 500-step margin is small.

## Follow-up exclusion: reduced-D constraint helpers

A quick `checkMesh -constant` dimension check was run after retaining this
candidate:

```text
log.codex_checkMesh_dims_20260606
```

Result:

```text
Mesh has 3 geometric (non-empty/wedge) directions (1 1 1)
Mesh has 3 solution (non-empty) directions (1 1 1)
Mesh OK.
```

Therefore the `dsmcParcel::move()` branch guarded by
`mesh.nGeometricD() < 3 || mesh.nSolutionD() < 3` is not active for this case.
Do not spend the next candidate on `meshTools::constrainToMeshCentre()` or
`meshTools::constrainDirection()` unless the target case changes.
