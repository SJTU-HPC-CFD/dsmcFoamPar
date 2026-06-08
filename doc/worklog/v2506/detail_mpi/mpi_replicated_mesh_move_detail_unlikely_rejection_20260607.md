# MPI replicated mesh move-detail branch hint attempt - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI only, `useOpenMP false`, `openmpThreads 1`,
  `mpirun -np 8`
- Retained DLB configuration:
  `replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);`

## Candidate

File:

```text
src/lagrangian/dsmc/parcels/dsmcParcel.C
```

Temporary source change:

- Added an inline `dsmcMoveDetailLikely()` wrapper around
  `__builtin_expect(enabled, false)`.
- Replaced the default-off `moveDetailProfile` branches in
  `dsmcParcel::move()`, `hitWallPatch()`, and `hitPatch()` with that branch
  hint.
- Did not remove the diagnostic counters or timers.

Reasoning:

- `moveDetailProfile` is normally false in formal runs.
- The default-off diagnostics leave several runtime branches inside the DSMC
  move hot path.
- The candidate tried to reduce branch cost without changing diagnostics or
  runtime controls.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_detail_unlikely_build_20260607.log
```

Build result: passed.

## 10-step no-detail smoke

Temporary control:

```text
endTime 1.e-06;
profileDetail false;
useOpenMP false;
openmpThreads 1;
```

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_detail_unlikely_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_mpionly_move_detail_unlikely_smoke_10step_20260607
```

Comparison against the previous same-tet no-detail smoke:

| metric | no-detail baseline | branch-hint smoke |
| --- | ---: | ---: |
| real | 10.86 s | 9.03 s |
| move+collide wall | 1.775343523 s | 1.539431429 s |
| move only | 1.223100367 s | 1.055965168 s |
| buildCellOccupancy | 0.107306958 s | 0.107809421 s |
| collision phase | 0.192508068 s | 0.183138816 s |
| post fields/output | 0.859970621 s | 0.863455345 s |
| full evolve wall | 2.433862566 s | 2.284059652 s |

Correctness:

- `Total Iterations = 10`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- final particles `2065864`
- stuck particles `0`
- final total energy `1.150621535`
- no Fatal/NaN/BAD TERMINATION match was found in the log.

Smoke decision: positive enough for one 500-step validation.

## 500-step validation

Formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_move_detail_unlikely_500step_20260607
```

Comparison against the current pure-MPI forced8 baseline from
`log.codex_mpi8_replicatedmesh_purempi_scope_corrected_rebuild_500step_20260607`:

| metric | current baseline | branch-hint candidate |
| --- | ---: | ---: |
| real | 131.16 s | 134.25 s |
| move+collide wall | 93.59239769 s | 98.63919566 s |
| move only | 70.58583509 s | 74.36248110 s |
| buildCellOccupancy | 7.037232476 s | 7.191969421 s |
| collision phase | 14.06214838 s | 15.23411885 s |
| post fields/output | 41.44612826 s | 40.62043329 s |
| full evolve wall | 129.1197774 s | 134.2197901 s |
| migration calls | 59 | 59 |
| DLB checks | 500 | 500 |
| DLB rebalances | 8 | 8 |

Correctness:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- final particles `2463609`
- stuck particles `0`
- final total energy `1.23883972`
- no Fatal/NaN/BAD TERMINATION match was found in the log.

## Decision

Rejected and reverted.

The 10-step smoke improvement did not survive the formal run.  The candidate
worsened the required metrics:

- `real +3.09 s`
- `move only +3.77664601 s`
- `move+collide wall +5.04679797 s`
- `full evolve wall +5.1000127 s`

The small post-field improvement is not enough to offset the move and evolve
regressions.

## Revert state

The source change was reverted, the case control remained restored to the
500-step forced8 pure-MPI hash, and the binary was rebuilt.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_detail_unlikely_revert_build_20260607.log
```

Restored control hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

There is no remaining unstaged diff in `dsmcParcel.C` or the case
`system/controlDict` from this candidate.
