# MPI replicated mesh full barycentric tracking port - 2026-06-08

## Scope

- Current checkout:
  `/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb`
- Reference source:
  `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
- Purpose: start the full reference-style barycentric particle/tracking port for
  the pure-MPI replicated-mesh move gap.

This is not the previous transient local barycentric hit-selection candidate.
That candidate avoided the storage/IO port and was rejected by the strict
500-step formal test. This round starts from the reference source's coherent
particle state model: persistent barycentric coordinates, tet topology state,
tracking API, IO, transfer, and patch handling.

## Initial Source State

The current worktree has existing staged changes. This port must not use
`git add` while working; new changes are kept as unstaged overlays until the
user decides what to do with the index.

Particle-related staged files at start:

- `src/lagrangian/basic/particle/particleTemplates.C`
- `src/lagrangian/dsmc/parcels/dsmcParcel.C`
- `src/lagrangian/dsmc/parcels/dsmcParcel.H`
- several `src/lagrangian/dsmc/*` replicated-mesh and timing files

Untracked backup files also exist:

- `src/lagrangian/basic/particle/particleTemplates.Cbkp`
- `src/lagrangian/dsmc/clouds/dsmcCloud.Cbkp`
- `src/lagrangian/dsmc/clouds/dsmcCloud.Hbkp`
- `src/lagrangian/dsmc/parcels/dsmcParcel.Cbkp`
- `src/lagrangian/dsmc/parcels/dsmcParcel.Hbkp`

## Reference Tracking Core

Reference particle state is not Cartesian `position_` plus legacy tet indices.
It stores:

- `barycentric coordinates_`
- `celli_`
- `tetFacei_`
- `tetPti_`
- `facei_`
- `stepFraction_`
- stuck-track guards `behind_` and `nBehind_`

Reference DSMC move calls `trackToAndHitFace(displacement, fraction, cloud, td)`
with remaining step fraction. The lower-level stationary tracking path projects
the displacement into the current tet's barycentric basis, picks the first
barycentric component that reaches zero, updates `coordinates_`, and then
updates tet/cell topology through `changeTet()`, `changeFace()`, or
`changeCell()`.

## Current Dependency Gaps

The full port cannot begin by swapping only `dsmcParcel.C` or
`particleTemplates.C`. The first missing layers are:

1. v1706-compatible `barycentric` and `barycentricTensor` primitives.
2. A v1706-compatible path for the newer lightweight `tetIndices` semantics, or
   explicit adaptation of the reference particle code to the old v1706
   `tetIndices` API.
3. `particle.H/C/I.H/Templates.C/IO.C` state and IO migration.
4. DSMC parcel move call-site migration from absolute-position tracking to
   displacement/fraction tracking.
5. Replicated flat-transfer packing audit for the changed particle storage.
6. Patch/cyclic/processor transfer semantics audit.

## Important Port Constraint

OpenFOAM-v2506 `Barycentric*.H` cannot be copied verbatim into this v1706
checkout. It uses newer type-trait forms such as `is_contiguous`,
`is_contiguous_label`, and `is_contiguous_scalar`; OF-v1706 uses the older
`contiguous<T>()` function style. The first implementation slice therefore
needs a v1706-compatible barycentric shim or an adapted primitive import.

## Execution Plan

1. Add or adapt the barycentric primitive layer in a way that compiles under
   OF-v1706/C++11.
2. Decide whether to shadow/adapt `tetIndices` locally or keep v1706
   `tetIndices` and rewrite the reference particle code against the old API.
3. Port `particle` state and stationary tracking in a compile-focused slice.
4. Port DSMC `move()` call shape.
5. Audit particle IO and replicated flat-transfer packing.
6. Build `liblagrangian+`, then `libdsmcFoam+`, before running any benchmark.

## Current Status

Implementation slice built and smoke-tested.

### Implemented In This Slice

- Added v1706-compatible barycentric primitives and compiled
  `particle/barycentric.C` into `liblagrangian+`.
- Added persistent barycentric state to the existing v1706 `particle` without
  replacing the legacy Cartesian `position_` storage model.
- Added a static-mesh `trackToAndHitFaceBarycentric(...)` wrapper with
  reference-style stationary tet tracking, patch dispatch reuse, and legacy
  fallback when the barycentric attempt cannot advance.
- Connected `dsmcParcel::move()` to the wrapper only for conservative cases:
  static mesh, uniform DSMC time step, not reduced-D Cartesian constrained, and
  not in wall-impact-distance cells.
- Added move-detail counters:
  `move detail barycentric calls` and `move detail barycentric fb`.

### Build Logs

```text
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_hybrid_lagrangian_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_hybrid_dsmc_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_hybrid_solver_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_hybrid_zero_progress_fallback_lagrangian_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_hybrid_zero_progress_fallback_dsmc_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_hybrid_zero_progress_fallback_solver_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_hybrid_counters_dsmc_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_hybrid_counters_solver_build_20260608.log
```

Latest build products:

```text
liblagrangian+.so  2026-06-08 04:57:32
libdsmcFoam+.so    2026-06-08 05:02:40
dsmcFoam+          2026-06-08 05:02:42
```

Build result: passed. The grep hits on `FatalIOErrorInFunction` in the logs are
OpenFOAM warning context lines, not compile/link failures.

### Smoke Validation

Case:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh
```

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_barycentric_hybrid_smoke_20260608
```

Temporary smoke controls:

```text
endTime 1.e-06;
profileSummary true;
profileDetail true;
moveDetailProfile true;
```

Final smoke log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_hybrid_counters_smoke_10step_20260608
```

Result:

| metric | value |
| --- | ---: |
| real | 11.94 s |
| move only | 1.433015986 s |
| full evolve wall | 2.214138123 s |
| move detail track calls | 30,410,444 |
| move detail barycentric calls | 30,410,444 |
| move detail barycentric fallback | 2,332,756 |
| move detail same-tet no-face | 15,791,800 |
| move detail internal tet only | 9,482,938 |
| move detail face hits | 5,135,706 |

Correctness smoke signals:

```text
Total Iterations = 10
Number of DSMC particles = 2065853
Number of stuck particles = 0
Collisions = 2794
Total energy = 1.150623417
OpenMP enabled = 0
Phase C auto DLB checks = 10
Phase C auto DLB rebalances = 0
```

No `Fatal`, `NaN`, `BAD TERMINATION`, or timeout string was found in the final
smoke log. The case `system/controlDict` was restored to the 500-step formal
settings after the smoke.

### Zero-Progress Fallback Fix

The first post-build smoke attempt was terminated manually after more than two
minutes without reaching the first time output:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_hybrid_smoke_10step_20260608
```

Cause found in the wrapper path: a barycentric attempt could return a completed
fraction of zero without hitting a face. `dsmcParcel::move()` multiplies `dt`
by that value, so `tEnd` did not decrease and the parcel loop retried forever.

Fix applied:

- match the reference stuck-guard behavior by advancing `stepFraction_` on
  stuck barycentric tracks;
- preserve and restore the original particle state before barycentric tracking;
- if the barycentric wrapper makes zero progress and does not hit a face,
  restore the state and call the legacy `trackToFace()` path.

### Current Limitations

This is only a bring-up smoke, not formal performance evidence. The current
slice still syncs barycentric coordinates from `position_` before each wrapper
call and has not yet ported full binary IO/transfer semantics for persistent
barycentric state. The fallback count is non-trivial, so the next stage should
audit the remaining topology/state mismatch before any 500-step performance
claim.

### Review Fixes And Transfer Coverage

Follow-up review found two correctness risks in the hybrid wrapper:

- `particleIO.C` still used old contiguous `&position_` binary block IO. After
  adding `coordinates_`, `behind_`, and `nBehind_` between `position_` and
  `cellI_`, that could corrupt the binary particle stream used by processor
  transfer. The fix keeps the old stream field set and order, but reads/writes
  each field explicitly.
- The barycentric wrapper synchronised coordinates after every boundary patch.
  For processor patches, transfer is completed later by
  `prepareForParallelTransfer()` and `correctAfterParallelTransfer()`. The
  wrapper now skips the immediate barycentric resync on processor patches and
  lets the next tracking entry resync from the received `position_/cell/tet`
  state.

Final build log after the review fix:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_review_formatfix_build_20260608.log
```

Build result: passed. No compile/link error keyword was found in the final
build log.

Final smoke logs after the review fix:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_review_final_smoke_10step_20260608
run/.../ourmesh/mpi8origin/log.codex_mpi8origin_barycentric_transfer_final_smoke_10step_20260608
```

Final smoke result:

| case | iterations | particles | stuck | energy | bary calls | fallbacks | processor hits | real |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| mpi8replicatedmesh | 10 | 2,065,832 | 0 | 1.150600843 | 30,410,908 | 2,333,046 | 0 | 9.33 s |
| mpi8origin -parallel | 10 | 2,065,849 | 0 | 1.150617584 | 30,404,199 | 2,329,476 | 35,088 | 4.83 s |

No `Fatal`, `NaN`, `BAD TERMINATION`, `timeout`, `MPI_ABORT`, or segmentation
string was found in either final smoke log. The `mpi8origin -parallel` smoke is
the targeted processor-transfer coverage for this review because it hit
processor patches 35,088 times.

Control backups used for the temporary 10-step runs:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_barycentric_review_smoke_20260608
doc/worklog/v2506/detail_mpi/mpi8origin_controlDict_before_barycentric_transfer_smoke_20260608
```

Both case `system/controlDict` files were restored to the original 500-step
settings after the smoke tests. These are still smoke tests only; formal
performance/correctness claims still require the 500-step `ourmesh` comparison
against `ourmeshbkp`.

### Formal 500-Step Comparison

Source state for this comparison:

- Current dirty-tree barycentric hybrid port, including the zero-progress
  fallback fix, explicit particle binary IO field ordering, and processor-patch
  resync guard.
- No git index operation was performed for this run.

Build products used:

```text
liblagrangian+.so  2026-06-08 05:15:30 +0800
libdsmcFoam+.so    2026-06-08 05:11:59 +0800
dsmcFoam+          2026-06-08 05:12:01 +0800
```

Control state:

- `mpi8replicatedmesh`: `endTime 5.e-05`, `deltaT 1.e-07`,
  `profileSummary true`, `profileDetail false`, `replicatedMesh true`,
  forced DLB steps `(120 170 220 270 320 370 420 470)`.
- `mpi8origin`: `endTime 5.e-05`, `deltaT 1.e-07`.
- `git diff` for both case `system/controlDict` files was zero after the runs.

Current formal logs:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_review_formal_500step_20260608
run/.../ourmesh/mpi8origin/log.codex_mpi8origin_barycentric_review_formal_500step_20260608
```

Baseline logs:

```text
run/.../ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610
run/.../ourmeshbkp/mpi8origin/log.mpi8
run/.../ourmeshbkp/mpi8origin/log.mpi8origin.confirm_pdFalse_20260605_024610
```

Correctness signals from current formal runs:

| case | iterations | particles | stuck | collisions | candidates | total energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| mpi8replicatedmesh current | 500 | 2,463,647 | 0 | 34,177 | 60,387 | 1.238944357 |
| mpi8origin current | 500 | 2,463,593 | 0 | 35,441 | 68,149 | 1.243171874 |

No `Fatal`, `NaN`, `BAD TERMINATION`, `timeout`, `MPI_ABORT`, or segmentation
string was found in either current formal log.

End-to-end/profile comparison:

| case | metric | current | baseline | change |
| --- | --- | ---: | ---: | ---: |
| mpi8replicatedmesh | full evolve wall | 121.1328433 s | 114.9207628 s | +5.41% slower |
| mpi8replicatedmesh | real | 122.14 s | 114.9207628 s | +6.28% slower |
| mpi8replicatedmesh | move only | 80.62990327 s | 47.1652636 s | +70.95% slower |
| mpi8replicatedmesh | rank wall max | 94.79 s | 99.13380862 s | -4.38% better |
| mpi8origin | final Stage 1 ExecutionTime | 177.96 s | 171.4 s | +3.83% slower |
| mpi8origin | real | 178.68 s | 171.4 s | +4.25% slower |
| mpi8origin | move only | 147.4417899 s | 120.3658396 s | +22.49% slower |
| mpi8origin | full evolve wall | 189.6118895 s | 163.0536101 s | +16.29% slower |

Replicated-mesh DLB comparison:

| metric | current | baseline | change |
| --- | ---: | ---: | ---: |
| migration calls | 59 | 59 | same |
| migration wall time | 2.490682605 s | 4.406303137 s | -43.48% better |
| particles per rank max/min | 1.315786739 | 1.775539064 | -25.89% better |
| rank wall max/min | 1.036749426 | 1.080575306 | -4.05% better |
| auto DLB checks | 500 | 500 | same |
| auto DLB rebalances | 8 | 8 | same |
| auto DLB wall max | 2.087319069 s | 25.67342357 s | -91.87% better |

Interpretation:

- The current replicated-mesh run improves load balance and DLB overhead, but
  does not improve 500-step end-to-end wall time against the preserved
  `ourmeshbkp` baseline.
- The main regression is still in move tracking: replicated-mesh `move only`
  increases from 47.17 s to 80.63 s, and origin `move only` increases from
  120.37 s to 147.44 s.
- Therefore the barycentric hybrid bring-up is correctness-positive under
  500-step testing, but not performance-positive yet. The next optimization
  should focus on reducing fallback/sync overhead or porting the fuller
  persistent barycentric path before claiming a move-speedup.

### Valid-Coordinate Cache and Same-Tet Fast Path

Root cause confirmed after the first formal run:

- The hybrid wrapper was still calling `syncCoordinatesFromPosition()` on every
  barycentric track attempt.  This preserved correctness but did not act like
  the reference implementation, where barycentric coordinates are persistent
  particle state.
- A large share of tracks remained same-tet/no-face.  These tracks do not need
  the full wrapper state-save/fallback path if the target barycentric point is
  strictly inside the current tet.

Implemented changes:

- Added `coordinatesValid_` to `particle`.
- `syncCoordinatesFromPosition()` now becomes an on-demand cache refresh:
  initial Cartesian constructors, binary reads, old tracker fallback, reduced-D
  position constraints, and external shock-reset translations invalidate it;
  successful barycentric updates keep it valid.
- Boundary/patch hits invalidate coordinates instead of eagerly resyncing.
  Processor-patch transfer remains safe because the send-side binary/flat
  transfer path still sends the old Cartesian fields and the receive side starts
  with invalid barycentric coordinates.
- Added a conservative same-tet fast path before the wrapper saves state:
  if the proposed endpoint has all four barycentric components greater than
  `SMALL`, update coordinates/position/stepFraction directly and return a
  completed track.

Source files touched in this optimization layer:

```text
src/lagrangian/basic/particle/particle.H
src/lagrangian/basic/particle/particle.C
src/lagrangian/basic/particle/particleI.H
src/lagrangian/basic/particle/particleIO.C
src/lagrangian/basic/particle/particleTemplates.C
src/lagrangian/dsmc/parcels/dsmcParcel.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
```

Build logs:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_valid_cache_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_valid_cache_sametet_fast_build_20260608.log
```

Both builds passed with only the known OpenFOAM-v1706 template warnings.

#### 10-Step Diagnostics

Valid-cache smoke logs:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_valid_cache_smoke_10step_20260608
run/.../ourmesh/mpi8origin/log.codex_mpi8origin_barycentric_valid_cache_smoke_10step_20260608
```

Same-tet fast smoke logs:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_valid_cache_sametet_fast_smoke_10step_20260608
run/.../ourmesh/mpi8origin/log.codex_mpi8origin_barycentric_valid_cache_sametet_fast_smoke_10step_20260608
```

Smoke comparison:

| case | variant | move only | full evolve | bary calls | fallbacks | processor hits | stuck |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| replicated | review fix | 1.472413352 s | 2.286895246 s | 30,410,908 | 2,333,046 | 0 | 0 |
| replicated | valid cache | 1.344919782 s | 2.150191794 s | 30,412,145 | 2,332,564 | 0 | 0 |
| replicated | valid cache + same-tet fast | 1.150808333 s | 1.890849379 s | 27,467,150 | 0 | 0 | 0 |
| origin | review fix | 2.420838506 s | 2.913935642 s | 30,404,199 | 2,329,476 | 35,088 | 0 |
| origin | valid cache | 1.941524455 s | 2.412688337 s | 30,406,898 | 2,330,433 | 35,101 | 0 |
| origin | valid cache + same-tet fast | 1.80668225 s | 2.236762088 s | 27,462,436 | 0 | 35,079 | 0 |

The same-tet fast path removed the short-run fallback counter because the
previous zero-progress retry cases are now completed by the strict same-tet
inside test.  The origin smoke still hit processor patches, so the transfer
path remains covered.

#### 500-Step Formal Re-Test

Formal logs for valid-cache only:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_valid_cache_formal_500step_20260608
run/.../ourmesh/mpi8origin/log.codex_mpi8origin_barycentric_valid_cache_formal_500step_20260608
```

Formal logs for final same-tet fast version:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_valid_cache_sametet_fast_formal_500step_20260608
run/.../ourmesh/mpi8origin/log.codex_mpi8origin_barycentric_valid_cache_sametet_fast_formal_500step_20260608
```

Control state:

- Formal temporary controls were `endTime 5.e-05`, `deltaT 1.e-07`,
  `profileSummary true`, `profileDetail false`, `moveDetailProfile false`.
- Both `mpi8replicatedmesh/system/controlDict` and
  `mpi8origin/system/controlDict` were restored after the runs.
- `git diff` for the two case `controlDict` files was empty after validation.

Correctness signals:

| case | final variant | iterations | particles | stuck | collisions | total energy |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| replicated | valid cache + same-tet fast | 500 | 2,463,656 | 0 | 34,332 | 1.238835075 |
| origin | valid cache + same-tet fast | 500 | 2,463,670 | 0 | 35,398 | 1.243232007 |

No `Fatal`, `NaN`, `BAD TERMINATION`, `timeout`, `MPI_ABORT`, segmentation, or
floating-point failure string was found in the final formal logs.  No
`dsmcFoam+`, `mpirun`, or `mpiexec` process was left running.

500-step performance comparison:

| case | metric | baseline | review fix | valid cache | final same-tet fast |
| --- | --- | ---: | ---: | ---: | ---: |
| replicated | full evolve wall | 114.9207628 s | 121.1328433 s | 114.8409521 s | 108.2546436 s |
| replicated | external real | 114.9207628 s | 122.14 s | 116.33 s | 118.32 s |
| replicated | move only | 47.1652636 s | 80.62990327 s | 74.13468215 s | 69.68940261 s |
| origin | Stage 500 / external real baseline | 171.4 s | 178.68 s | 172.96 s | 170.86 s |
| origin | full evolve wall profile | 163.0536101 s | 189.6118895 s | 169.1924754 s | 167.1265645 s |
| origin | move only | 120.3658396 s | 147.4417899 s | 130.1203701 s | 127.9978233 s |

Final same-tet fast deltas:

| case | metric | vs review fix | vs valid cache | vs baseline |
| --- | --- | ---: | ---: | ---: |
| replicated | full evolve wall | -10.63% | -5.74% | -5.80% |
| replicated | move only | -13.57% | -6.00% | +47.76% |
| replicated | external real | -3.13% | +1.71% | +2.96% |
| origin | full evolve wall profile | -11.86% | -1.22% | +2.50% |
| origin | move only | -13.19% | -1.63% | +6.34% |
| origin | external real | -4.38% | -1.21% | -0.32% |

Final replicated-mesh DLB state:

| metric | final same-tet fast | baseline |
| --- | ---: | ---: |
| migration calls | 59 | 59 |
| migration wall time | 2.344602125 s | 4.406303137 s |
| particles per rank max/min | 1.528281083 | 1.775539064 |
| rank wall max/min | 1.039600551 | 1.080575306 |
| auto DLB checks | 500 | 500 |
| auto DLB rebalances | 8 | 8 |
| auto DLB wall max | 2.072637642 s | 25.67342357 s |

Interpretation after both optimizations:

- The barycentric hybrid is no longer a clear performance regression.  In the
  replicated-mesh profile summary, full evolve improves from 114.92 s baseline
  to 108.25 s, while DLB overhead remains much lower than baseline.
- The replicated external `real` timer did not improve versus baseline
  (118.32 s vs 114.92 s).  This should be reported separately from the solver
  profile summary; do not claim universal end-to-end speedup for replicated
  mesh from this single run.
- The origin `real` timer is slightly better than the old `log.mpi8` baseline
  (170.86 s vs 171.4 s), but its profile `move only` is still 6.34% slower than
  the profile baseline.
- The remaining bottleneck is still the tracking core itself.  Even after
  same-tet fast completion, replicated `move only` remains 69.69 s versus the
  old 47.17 s baseline.  There is no evidence here that another wrapper-level
  micro-cache will close that gap.
- The next meaningful choice is architectural: either keep this correctness-
  compatible barycentric hybrid because it improves replicated profile full
  evolve and keeps origin real competitive, or benchmark a selectable fallback
  to the earlier optimized Cartesian same-tet path if absolute replicated
  external `real` time is the primary objective.

## Runtime Switch And Final Default

The architectural choice above was tested directly by adding a runtime
`controlDict` switch:

```text
barycentricTracking true|false;
```

If the entry is absent, the production default is now `false`, which keeps the
optimized Cartesian same-tet/tet-walk path as the default move tracker.  The
barycentric hybrid remains available for explicit A/B and future reference
work, but it is no longer the default for the replicated-mesh production path.

Code/control changes:

- `dsmcCloud` stores `barycentricTracking_`, reads it from `controlDict`, and
  prints the active value in the profile summary.
- `dsmcParcel::move()` only enters the barycentric wrapper when the switch is
  true and the existing conservative guards pass.
- When the switch is false, the Cartesian path no longer maintains barycentric
  cache state after each legacy `trackToFace()` call.

Build logs:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_runtime_switch_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_default_cartesian_build_20260608.log
```

10-step smoke logs:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_switch_off_smoke_10step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_switch_on_smoke_10step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_default_cartesian_smoke_10step_20260608
```

Smoke switch verification:

| variant | iterations | stuck | barycentric tracking | bary calls | fallbacks | move only | full evolve | real |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| explicit false | 10 | 0 | 0 | 0 | 0 | n/a | 2.111812417 s | 11.04 s |
| explicit true | 10 | 0 | 1 | 27,467,340 | 0 | n/a | 1.942417104 s | 10.66 s |
| default absent | 10 | 0 | 0 | 0 | 0 | 1.200054566 s | 2.021418635 s | 11.17 s |

500-step formal logs:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_switch_off_formal_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_default_cartesian_formal_500step_20260608
```

Formal controls were the strict full-fields replicated-mesh controls:

```text
endTime 5.e-05;
deltaT 1.e-07;
profileSummary true;
profileDetail false;
```

No `barycentricTracking` entry is left in the restored production
`controlDict`; the absent-entry default is the tested Cartesian default.

Final 500-step comparison:

| replicated metric | barycentric on | explicit false before default cleanup | final default Cartesian | old strict baseline |
| --- | ---: | ---: | ---: | ---: |
| real | 118.32 s | 118.59 s | 114.42 s | 114.9207628 s |
| move only | 69.68940261 s | 68.91372364 s | 65.4460646 s | 47.1652636 s |
| post fields/output | 18.72222119 s | 18.46509146 s | 18.24213585 s | n/a |
| full evolve wall | 108.2546436 s | 108.2008378 s | 103.9147289 s | 114.9207628 s |
| particles per rank max/min | 1.528281083 | 3.635513688 | 1.508825188 | 1.775539064 |
| rank wall max/min | 1.039600551 | 1.030386109 | 1.025803369 | 1.080575306 |
| DLB rebalances | 8 | 8 | 8 | 8 |

Final correctness:

```text
Total Iterations = 500
Number of DSMC particles = 2463663
Number of stuck particles = 0
Collisions = 33856
Total energy = 1.238815995
Barycentric tracking = 0
```

No `Fatal`, `NaN`, `BAD TERMINATION`, `timeout`, `MPI_ABORT`, segmentation, or
floating-point failure string was found in the final default-Cartesian formal
log.  No `dsmcFoam+`, `mpirun`, or `mpiexec` process was left running.

Final decision:

- Keep the runtime switch.
- Use Cartesian tracking as the production default for replicated mesh.
- Keep barycentric tracking as explicit opt-in because it is correct and useful
  for future reference-style work, but it does not beat the final Cartesian
  default in the 500-step replicated-mesh formal run.

Final status versus the strict old baseline:

- `full evolve wall`: `114.9207628 s -> 103.9147289 s`, 9.58% faster.
- external `real`: `114.9207628 s -> 114.42 s`, slightly faster but close
  enough that this should be treated as parity/slight win rather than a large
  end-to-end claim.
- `move only`: still slower than the historical strict baseline
  (`65.4460646 s` vs `47.1652636 s`), so the remaining raw tracking gap is not
  fully closed.  The end-to-end replicated profile is nevertheless now ahead
  because DLB/migration/post-field overheads are much lower.
