# MPI replicated mesh trackingData per-step cache rejection - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI only, `useOpenMP false`, `openmpThreads 1`,
  `mpirun -np 8`
- Field semantics: strict full-fields, unchanged `fieldPropertiesDict`
- Baseline control hash:
  `8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2`
- Baseline fieldProperties hash:
  `73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9`

The candidate only tried to cache per-time-step invariants in
`dsmcParcel::trackingData`:

- `constrainCartesianTracking`
- `trackerActive`
- `uniformDeltaT`

The values were filled once in `dsmcCloud::evolve_moveAndCollide()` and read in
`dsmcParcel::move()`.  It did not alter tet/tri tracking geometry, migration,
post fields/output, no-fields mode, or OpenMP runtime settings.

Build log for the candidate:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_trackingdata_cache_build_20260607.log
```

## 10-step smoke

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_trackingdata_cache_smoke_10step_20260607
```

Result:

| metric | current strict smoke reference | candidate | delta |
| --- | ---: | ---: | ---: |
| real | 12.10 s | 12.21 s | +0.11 s |
| move only | 1.138835794 s | 1.031074239 s | -0.107761555 s |
| full evolve wall | 2.633190181 s | 2.516925162 s | -0.116265019 s |
| OpenMP enabled | 0 | 0 | same |
| DLB rebalances | 0 | 0 | same |

Correctness markers:

| item | value |
| --- | ---: |
| Total Iterations | 10 |
| final particles | 2065839 |
| stuck particles | 0 |
| final collisions | 2786 |
| final total energy | 1.150604113 |

The smoke showed a move-path signal, so the candidate was promoted to a
200-step signal test.

## 200-step signal test

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_trackingdata_cache_smoke_200step_20260607
```

Result:

| metric | vsize025 screen | itr100 screen | candidate |
| --- | ---: | ---: | ---: |
| real | 55.00 s | 54.93 s | 53.52 s |
| move only | 28.3369108 s | 25.19764393 s | 25.35173156 s |
| buildCellOccupancy | 2.607412446 s | 2.43696367 s | 2.279840184 s |
| collision phase | 5.402109378 s | 4.780814507 s | 4.929787418 s |
| post fields/output | 17.03748328 s | 17.1258495 s | 15.67096643 s |
| full evolve wall | 49.89051853 s | 47.81540611 s | 46.3701355 s |
| OpenMP enabled | 0 | 0 | 0 |
| DLB rebalances | 2 | 2 | 2 |

Correctness markers:

| item | value |
| --- | ---: |
| Total Iterations | 200 |
| final particles | 2220220 |
| stuck particles | 0 |
| final collisions | 15415 |
| final total energy | 1.169670371 |

The 200-step result had end-to-end signal, so the candidate was promoted to the
formal 500-step gate.

## 500-step formal gate

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_trackingdata_cache_500step_20260607
```

Result:

| metric | strict full-fields baseline | candidate | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 131.33 s | -2.50 s |
| move only | 72.9567669 s | 74.64296316 s | +1.68619626 s |
| move+collide wall | n/a | 99.46871365 s | n/a |
| buildCellOccupancy | n/a | 6.594258188 s | n/a |
| collision phase | n/a | 15.70982533 s | n/a |
| post fields/output | n/a | 38.50132308 s | n/a |
| full evolve wall | n/a | 133.3647574 s | n/a |
| OpenMP enabled | 0 | 0 | same |
| DLB checks | 500 | 500 | same |
| DLB rebalances | 8 | 8 | same |

Correctness markers:

| item | value |
| --- | ---: |
| Total Iterations | 500 |
| final particles | 2463655 |
| stuck particles | 0 |
| final collisions | 34246 |
| final total energy | 1.23918131 |
| migration calls | 59 |
| migration wall time | 2.608120396 s |
| particles max/min | 1.716807519 |
| rank wall max/min | 1.051932793 |

## Decision

Rejected.  Do not retain this source candidate.

The formal result improved end-to-end `real` time by `2.50 s`, but it missed the
strict acceptance gate because `move only` regressed by `1.68619626 s` versus
the current strict full-fields baseline.  This is not a valid retained
move-path optimization under the current gate, even though the 10-step and
200-step runs had partial signal.

## Revert state

The source candidate was reverted in:

```text
src/lagrangian/dsmc/parcels/dsmcParcel.H
src/lagrangian/dsmc/parcels/dsmcParcel.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
```

The revert only removed the three trackingData cache fields and restored the
existing direct reads in `dsmcParcel::move()`.  It did not touch
`src/lagrangian/basic/particle/particleTemplates.C`.

Rebuild after revert:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_trackingdata_cache_revert_build_20260607.log
```

The rebuild completed with `=== Done ===`; the remaining messages are the
known OpenFOAM template warnings.

Restored case hashes:

```text
system/controlDict
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2

system/fieldPropertiesDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
```
