# MPI replicated-mesh move tracker atomic candidate - 2026-06-07

## Scope correction

This note is retained only as a record of an out-of-scope attempt.

The user clarified after this run that the current round must consider only
pure MPI8 replicated-mesh DLB with OpenMP completely disabled.  Therefore the
OpenMP2 result below must not be used as the retained result, acceptance
evidence, or comparison basis for this round.

The source candidate was reverted, the case control was restored to
`useOpenMP false; openmpThreads 1;`, and the solver was rebuilt.  The valid
reference baseline for this round is:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610
```

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Reference source:
  `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
- Candidate purpose: port the reference move-path optimization that replaces a
  global OpenMP critical around face-transition tracking with atomic flux
  accumulation inside `dsmcFaceTracker`.
- Performance target requested in this round: reach the historical
  MPI8+OpenMP2 reference class around `113 s`.

Important comparison boundary:

- `113 s` is a MPI8+OpenMP2 target, not a pure-MPI target.
- The latest pure-MPI retained baseline before this candidate was forced8:
  `real 133.83 s`, `move only 72.9567669 s`.

## Source candidate

Changed files:

- `src/lagrangian/dsmc/faceTracker/dsmcFaceTracker.C`
- `src/lagrangian/dsmc/parcels/dsmcParcel.C`

Reference-derived change:

- Added `_OPENMP` handling in `dsmcFaceTracker.C`.
- In `trackFaceTransition()`, route all `parcelIdFlux_` and `massIdFlux_`
  updates through a local `addFlux()` helper.
- When OpenMP move is active inside a parallel region, `addFlux()` uses
  `#pragma omp atomic update` for each scalar flux update.
- With those atomics in place, `trackParcelFaceTransitionThreadSafe()` in
  `dsmcParcel.C` directly calls
  `cloud.tracker().trackParcelFaceTransition(p)` instead of serialising every
  face transition through one `dsmcMoveTracker` critical section.

Reasoning:

- The reference implementation does not wrap tracker transition calls in a
  global critical; it makes the tracker field updates thread-safe at the update
  site.
- The current 10-step move-detail diagnostics showed millions of face-hit
  events per run.  A global critical around those calls is therefore a plausible
  OpenMP2 bottleneck.
- Directly removing the critical without atomic updates would be a data race;
  the tracker atomic update is the required paired change.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_tracker_atomic_build_20260607.log
```

Build result: passed.

## 10-step OpenMP2 smoke

Control backup before changing the forced8 pure-MPI case:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_tracker_atomic_openmp2_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Temporary controls:

```text
useOpenMP true;
openmpThreads 2;
endTime 1.e-06;
profileSummary true;
profileDetail false;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_openmp2_tracker_atomic_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 8.99 s |
| move+collide wall | 1.41770603 s |
| move only | 1.018788372 s |
| buildCellOccupancy | 0.09662842 s |
| collision phase | 0.070577728 s |
| post fields/output | 0.86746907 s |
| full evolve wall | 2.112235493 s |
| final particles | 2065847 |
| stuck particles | 0 |
| final collisions | 2858 |
| final candidates | 2975 |
| final total energy | 1.150620479 |
| DLB checks | 10 |
| DLB rebalances | 0 |

Correctness: no `Fatal`, `Segmentation`, `Floating`, `NaN`, `nan`, or
`BAD TERMINATION` match was found.

Smoke decision: promote to 500-step formal.  The previous OpenMP2 10-step
record in this worklog series was `real 10.70 s`, `move only 1.288977197 s`,
so the smoke signal was strongly positive.

## 500-step OpenMP2 formal

Control backup before restoring 500 steps:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_tracker_atomic_openmp2_formal_20260607
sha256 1c01a2753f8e20d8fb49733020a4905f1d6fd1542c10ae1ebc2596cc71b25bc8
```

Formal controls:

```text
useOpenMP true;
openmpThreads 2;
endTime 5.e-05;
profileSummary true;
profileDetail false;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Formal log:

```text
log.codex_mpi8_replicatedmesh_openmp2_tracker_atomic_500step_20260607
```

Result:

| metric | historical OpenMP2 reference | tracker atomic candidate |
| --- | ---: | ---: |
| real | 113.48 s | 104.84 s |
| move only | 69.34624257 s | 56.24644404 s |
| buildCellOccupancy | 5.975890432 s | 5.90283204 s |
| collision phase | 5.523175253 s | 4.633412965 s |
| post fields/output | 33.63726090 s | 37.47362874 s |
| full evolve wall | 104.7277518 s | 100.7882073 s |
| DLB rebalances | 4 | 8 |

Additional formal metrics:

| metric | value |
| --- | ---: |
| move+collide wall | 68.68705888 s |
| total profiled | 104.2563178 s |
| migration calls | 59 |
| migration wall | 1.955936549 s |
| particles max/min | 1.645531292 |
| rank wall max/min | 1.088521999 |
| DLB checks | 500 |
| final particles | 2463685 |
| stuck particles | 0 |
| final collisions | 34235 |
| final candidates | 61587 |
| final total energy | 1.238876671 |

Correctness:

- `Total Iterations = 500`
- `OpenMP enabled = 1`
- `OpenMP max threads = 2`
- No `Fatal`, `Segmentation`, `Floating`, `NaN`, `nan`, or `BAD TERMINATION`
  match was found.

Current formal OpenMP2 controlDict hash:

```text
bd07cb5ab0c716a246e273fd8f2e957cdab5a1797c5930b3ff6bd1bb44c40bc3
```

## Decision before scope correction

This was initially treated as retained, but that decision is superseded by the
scope correction above.

This candidate beats the requested `113 s` reference class by `8.64 s` in a
500-step run and reduces `move only` by `13.09979853 s` relative to the
historical OpenMP2 reference.

The source change is specifically an OpenMP move optimization.  It does not
explain or close the pure-MPI gap by itself, because pure MPI does not enter the
OpenMP tracker critical/atomic path.  It is therefore not retained for this
pure-MPI round.

## Further optimization directions after the target is reached

The requested target is already reached (`real 104.84 s < 113 s`).  If this
round continues beyond the target, the remaining high-value directions are:

1. Reduce `post fields/output` cost, now `37.47362874 s`, because it is larger
   than in the historical OpenMP2 reference.
2. Continue old tracking-core work only with better diagnostics; prior tet-walk
   micro-candidates had poor 500-step stability.
3. Treat a full barycentric `trackToAndHitFace()` particle API migration as a
   larger port, not a narrow move candidate, because it crosses base particle
   storage, IO, transfer, and cyclic/patch semantics.
