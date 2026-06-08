# MPI replicated mesh post flat-serial sample-cache attempt - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI only, `useOpenMP false`, `openmpThreads 1`,
  `mpirun -np 8`
- Formal retained control after this attempt:
  `replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);`

## Candidate

File:

```text
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
```

Temporary source change:

- In the serial branch of `dsmcVolSharedSampleCache::build()`, prefer the
  already-built flat occupancy data:
  `occupancyCellOffsets_ + occupancyOrderedParcels_`.
- Keep the old cloud-list scan as fallback.
- Motivation: the current pure-MPI post timer is high, and
  `dsmcVolFields` already receives cell-ordered occupancy from
  `dsmcCloud::buildCellOccupancy()`.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_flat_serial_build_20260607.log
```

Build result: passed.

## 10-step smoke

Temporary controls:

```text
endTime 1.e-06;
profileDetail true;
profilePostDetail true;
useOpenMP false;
openmpThreads 1;
```

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_post_flat_serial_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_mpionly_post_flat_serial_smoke_10step_20260607
```

Comparison against the recent same-control post-detail smoke
`log.codex_mpi8_replicatedmesh_retained_forced8_move_detail_smoke_10step_20260607`:

| metric | retained smoke | flat-serial candidate |
| --- | ---: | ---: |
| real | 11.69 s | 11.32 s |
| move only | 1.322484047 s | 1.143610321 s |
| post fields/output | 0.945544449 s | 0.896456407 s |
| full evolve wall | 2.763975522 s | 2.298740528 s |
| sample accumulation | 0.255075814 s | 0.369975288 s |
| shared cache build | 0.175329391 s | 0.264224101 s |
| parcel accumulate | 0.155086605 s | 0.246493887 s |
| field combine | 0.079714644 s | 0.105713237 s |
| cell reduction | 0.009494015 s | 0.010890903 s |
| boundary accumulation | 0.010482669 s | 0.005722558 s |

Correctness:

- `Total Iterations = 10`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- final particles `2065846`
- stuck particles `0`
- final total energy `1.150623817`
- no Fatal/NaN/BAD TERMINATION match was found in the log.

## Decision

Rejected and reverted.

Although the 10-step external wall and total post timer moved slightly lower,
the targeted post sub-timers worsened:

- `sample accumulation` increased by `0.114899474 s`;
- `shared cache build` increased by `0.088894710 s`;
- `parcel accumulate` increased by `0.091407282 s`;
- `field combine` increased by `0.025998593 s`.

The small `post fields/output` improvement is not a reliable signal because it
does not come from the targeted sample-cache work.  Do not promote this
candidate to 500-step validation.

## Revert state

The source change was reverted, the case control was restored, and the binary
was rebuilt.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_flat_serial_revert_build_20260607.log
```

Restored control hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

There is no remaining unstaged diff in `dsmcVolFields.C` or the case
`system/controlDict` from this candidate.
