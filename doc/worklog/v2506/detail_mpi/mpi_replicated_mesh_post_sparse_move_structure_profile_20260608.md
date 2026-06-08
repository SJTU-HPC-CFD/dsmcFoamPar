# MPI replicated mesh post-sparse move structure profile - 2026-06-08

## Scope

- Mode: pure MPI8 replicated mesh DLB only, no OpenMP.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Source state: accepted `boundaryMeasurements` sparse clean retained.
- Purpose: diagnostic-only split after sparse clean, to decide whether the
  remaining move cost is dominated by wall measurement writes, replicated-mesh
  migration/DLB, or the original tet-walk tracking path.

This is not a formal 500-step performance comparison.  The temporary patch
adds detail timers and changes timing overhead.

Strict full-fields hashes after cleanup:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

## Temporary Instrumentation

Existing `moveDetailProfile` already times:

- `trackToFace()` in the parcel move loop;
- face tracker transitions;
- boundary model calls.

This run temporarily added cloud-level timers for replicated-mesh structure:

- pre-move migration/distribution;
- post-move periodic migration;
- manual DLB reassignment;
- auto DLB checks/rebalances.

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_post_sparse_move_structure_signal_20260608
```

Temporary controls:

```text
endTime 2.e-05;
profileSummary true;
profileDetail true;
moveDetailProfile true;
moveStructureProfile true;
```

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_sparse_move_structure_profile_build_20260608.log
```

Build result: passed, with only the known OpenFOAM-v1706 template warnings.

## 200-step Signal

Run command shape:

```text
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_post_sparse_move_structure_signal_200step_20260608
```

Aggregate profile:

| metric | value |
| --- | ---: |
| external real | 49.77 s |
| solver profile steps | 200 |
| move+collide wall | 36.10882215 s |
| move only | 28.04566287 s |
| buildCellOccupancy | 2.348702429 s |
| collision phase | 4.791171992 s |
| post fields/output | 6.769591856 s |
| full evolve wall | 41.58454539 s |

Move detail:

| metric | value |
| --- | ---: |
| move detail parcels | 428,041,919 |
| move detail track calls | 533,643,758 |
| same-tet no-face | 232,595,387 |
| internal tet only | 195,284,306 |
| face hits | 105,764,065 |
| processor hits | 0 |
| cyclic hits | 0 |
| patch hits | 439,207 |
| track max | 14.84928018 s |
| tracker max | 0 s |
| boundary max | 0.04610543 s |

Move structure timers:

| metric | value |
| --- | ---: |
| pre migration calls | 1 |
| post migration calls | 20 |
| manual DLB calls | 0 |
| auto DLB checks | 200 |
| auto DLB rebalances | 2 |
| pre migration | 0.408016217 s |
| post migration | 0.323874592 s |
| manual DLB | 0 s |
| auto DLB | 0.596339039 s |

Replicated-mesh summary:

| metric | value |
| --- | ---: |
| migration calls | 23 |
| migration wall time | 0.744608988 s |
| particles per rank max/min | 2.689175003 |
| rank wall time max/min | 1.034134008 |
| Phase C checks | 200 |
| Phase C rebalances | 2 |
| Phase C triggered checks | 2 |
| Phase C auto DLB wall max | 0.59499556 s |
| Phase C migration max | 0.462937643 s |

Correctness/check lines:

```text
Total Iterations = 200
Number of DSMC particles = 2220258
Number of stuck particles = 0
Collisions = 15433
Total energy = 1.169687546
OpenMP enabled = 0
```

No `Fatal`, `NaN`, or `BAD TERMINATION` string was found.

## Interpretation

The post-sparse-clean remaining move cost is not explained by boundary
measurement writes:

- boundary model timing is only `0.04610543 s` in the 200-step diagnostic;
- patch hits are `439,207` against `533,643,758` track calls;
- the rejected last-touch fast path is consistent with this result.

It is also not dominated by replicated-mesh migration or DLB:

- total migration wall time is `0.744608988 s`;
- post-move migration measured in the cloud path is `0.323874592 s`;
- auto DLB checks/rebalances are `0.596339039 s`;
- together these are small relative to `move only = 28.04566287 s` and
  `full evolve wall = 41.58454539 s`.

The remaining move-side work is still the original tet-walk tracking path:

- `trackToFace()` max timer is `14.84928018 s` under detail profiling;
- previous sampled tet-walk profile still applies: `findTris()` is about
  `39.62%`, topology transition about `24.88%`, area-normal/plane-base setup
  about `16.29%`, metadata reconstruction about `9.47%`, and outer lambda
  selection about `9.74%` of the sampled tet-walk split.

This points away from:

- wall-measurement accessor micro-optimisation;
- migration/auto-DLB tuning as the primary move solution;
- full fixed tet-geometry cache;
- narrow barycentric partial tracking patches.

The plausible remaining directions are larger:

1. a current-tree algorithmic change that reduces `findTris()` /
   `tetLambda()` / `tetNeighbour()` work without changing particle storage; or
2. a full barycentric particle/tracking-stack port with IO, transfer, cyclic,
   and replicated flat-transfer semantics audited together.

The second direction is too broad for this narrow pure-MPI move round.  The
first direction needs a new invariant or algorithmic simplification; the old
one-off tracking micro-candidates should remain rejected.

## Cleanup

Temporary source instrumentation was removed from:

```text
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/clouds/dsmcCloud.H
```

Cleanup build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_sparse_move_structure_profile_cleanup_build_20260608.log
```

Cleanup build result: passed, with only the known OpenFOAM-v1706 template
warnings.

Post-cleanup checks:

- no `moveStructureProfile`, `move structure`, or `moveStructure` symbols
  remain in `src/lagrangian` or the case `system/controlDict`;
- strict full-fields `controlDict` and `fieldPropertiesDict` hashes are
  restored to the values shown above;
- `git diff --check` over the files touched in this diagnostic passed.

