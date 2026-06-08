# MPI replicated mesh direct Utracking rejection - 2026-06-07

## Scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI8 replicated mesh, `useOpenMP false`, `openmpThreads 1`
- Baseline for this decision: `log.codex_mpi8_replicatedmesh_purempi_scope_corrected_rebuild_500step_20260607`

## Candidate

- File: `src/lagrangian/dsmc/parcels/dsmcParcel.C`
- Function: `dsmcParcel::move()`
- Idea: for 3D Cartesian tracking where `constrainCartesianTracking` is false,
  avoid copying `U_` into `Utracking` on every move loop and use a direct
  displacement for the `trackToFace(..., true)` target.

This was intentionally separate from the previously rejected
Cartesian/reduced-D condition cache.  It did not change the reduced-D branch or
boundary interaction semantics.

## Build

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_dsmcparcel_direct_Utracking_build_20260607.log
```

Build result: passed.

## 10-step smoke

Temporary control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_direct_Utracking_smoke_20260607
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_mpionly_direct_Utracking_smoke_10step_20260607
```

| metric | same-tet no-detail baseline | direct-Utracking |
| --- | ---: | ---: |
| real | 10.86 s | 9.12 s |
| move+collide wall | 1.775343523 s | 1.705514051 s |
| move only | 1.223100367 s | 1.147873963 s |
| post fields/output | 0.859970621 s | 0.94129796 s |
| full evolve wall | 2.433862566 s | 2.552192378 s |

The smoke result improved `real` and `move only` but worsened full evolve and
post time, so it was promoted to a 500-step formal check rather than retained.

## 500-step formal

Formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_direct_Utracking_500step_20260607
```

| metric | current forced8 baseline | direct-Utracking |
| --- | ---: | ---: |
| real | 131.16 s | 133.20 s |
| move only | 70.58583509 s | 75.99879936 s |
| buildCellOccupancy | 7.037232476 s | 7.443867945 s |
| collision phase | 14.06214838 s | 15.31489044 s |
| post fields/output | 41.44612826 s | 42.15341799 s |
| full evolve wall | 129.1197774 s | 132.6370253 s |
| DLB rebalances | 8 | 8 |

Correctness:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- stuck particles remained `0`
- no `Fatal`, `Segmentation`, `Floating`, `NaN`, `nan`, or `BAD TERMINATION`
  match was found.

## Decision

Rejected and reverted.  The 500-step formal result regressed every required
timer relative to the current pure-MPI forced8 baseline.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_dsmcparcel_direct_Utracking_revert_build_20260607.log
```

Revert build result: passed.  After the revert, `dsmcParcel.C` has no residual
diff from this candidate, and the formal `controlDict` hash is again:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```
