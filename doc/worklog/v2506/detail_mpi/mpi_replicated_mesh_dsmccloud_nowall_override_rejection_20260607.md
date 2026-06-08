# MPI replicated mesh dsmcCloud no-wall override rejection - 2026-06-07

## Scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI8 replicated mesh, `useOpenMP false`, `openmpThreads 1`
- Baseline for this decision: `log.codex_mpi8_replicatedmesh_purempi_scope_corrected_rebuild_500step_20260607`

## Candidate

- File: `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- Change: add an inline `hasWallImpactDistance() const final { return false; }`
  override for `dsmcCloud`.
- Purpose: make the DSMC no-wall-impact-distance semantics explicit for the
  templated `particle::trackToFace(..., true)` hot path.

## Build

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_dsmccloud_nowall_override_build_20260607.log
```

Build result: passed, with only the existing OpenFOAM template warnings.

## 10-step smoke

Temporary control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_dsmccloud_nowall_smoke_20260607
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_mpionly_dsmccloud_nowall_smoke_10step_20260607
```

| metric | same-tet no-detail baseline | no-wall override |
| --- | ---: | ---: |
| real | 10.86 s | 11.43 s |
| move+collide wall | 1.775343523 s | 1.58839893 s |
| move only | 1.223100367 s | 1.089370857 s |
| post fields/output | 0.859970621 s | 0.875652929 s |
| full evolve wall | 2.433862566 s | 2.325763264 s |

The smoke result improved the move timer but not external wall time, so the
candidate was promoted to a 500-step formal check rather than retained.

## 500-step formal

Formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_dsmccloud_nowall_500step_20260607
```

| metric | current forced8 baseline | no-wall override |
| --- | ---: | ---: |
| real | 131.16 s | 132.54 s |
| move+collide wall | not recorded in correction note table | 97.48059437 s |
| move only | 70.58583509 s | 74.20316995 s |
| buildCellOccupancy | 7.037232476 s | 6.992181256 s |
| collision phase | 14.06214838 s | 14.62113283 s |
| post fields/output | 41.44612826 s | 40.31207845 s |
| full evolve wall | 129.1197774 s | 130.9484461 s |
| DLB rebalances | 8 | 8 |

Correctness:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- stuck particles remained `0`
- no `Fatal`, `Segmentation`, `Floating`, `NaN`, `nan`, or `BAD TERMINATION`
  match was found.

## Decision

Rejected and reverted.  The 500-step formal result regressed external wall time,
move-only time, collision time, and full evolve wall time relative to the current
pure-MPI forced8 baseline.  The source change was removed.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_dsmccloud_nowall_override_revert_build_20260607.log
```

Revert build result: passed, with only the existing OpenFOAM template warnings.
After the revert, `dsmcCloud.H` has no residual diff from this candidate, and
the formal `controlDict` hash is again:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```
