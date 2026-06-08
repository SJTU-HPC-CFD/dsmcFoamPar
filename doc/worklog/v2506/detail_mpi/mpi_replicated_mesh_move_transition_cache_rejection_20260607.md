# MPI replicated mesh move transition-cache rejection - 2026-06-07

## Scope

- Mode: pure MPI8 replicated mesh, no OpenMP evidence.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Candidate feasibility note:
  `doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_move_transition_cache_feasibility_20260607.md`
- Formal comparison target: current strict full-fields forced8 baseline:
  `log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607`

## Candidate

Implemented a fixed-mesh topology table for non-face tet transitions only:

```text
(cellI, tetFaceI, tetPtI, triI=1/2/3) -> (next tetFaceI, next tetPtI)
```

Files touched during the candidate:

- `src/lagrangian/basic/particle/particle.H`
- `src/lagrangian/basic/particle/particleTemplates.C`
- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
- case `system/controlDict`

The candidate did not alter tri0 face-hit handling, `findTris()`,
`tetLambda()`, wall/boundary/tracker handling, replicated-mesh migration, or
the retained same-tet non-normalised helper.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_transition_cache_build_20260607.log
```

Build result: passed after fixing a helper-scope `Foam::FatalError`
qualification issue.  The retained log ends with `=== Done ===`.

Runtime cache size reported on each rank:

```text
states=1249812
transitions=3749436
table=28.60592651 MiB
offsets=3.184906006 MiB
```

## 10-step smoke

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_transition_cache_smoke_20260607
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_transition_cache_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| real | 10.30 s |
| move+collide wall | 1.903644187 s |
| move only | 1.155141187 s |
| post fields/output | 1.006511989 s |
| full evolve wall | 2.809113787 s |
| particles per rank max/min | 1.638461346 |
| DLB rebalances | 0 |

Correctness/check lines:

```text
Total Iterations = 10
Number of DSMC particles = 2065850
Collisions = 2775
Total energy = 1.150620014
```

No addressing failure or fatal error occurred.  The signal was not strong, but
it was clean enough to run the 200-step gate.

## 200-step signal

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_transition_cache_signal_20260607
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_transition_cache_signal_200step_20260607
```

Result:

| metric | value |
| --- | ---: |
| real | 55.45 s |
| move+collide wall | 36.80751543 s |
| move only | 28.07271777 s |
| buildCellOccupancy | 2.585665653 s |
| collision phase | 5.405701958 s |
| post fields/output | 17.67378473 s |
| full evolve wall | 50.71722401 s |
| particles per rank max/min | 3.320858831 |
| DLB rebalances | 2 |

Correctness/check lines:

```text
Total Iterations = 200
Number of DSMC particles = 2220201
Collisions = 15513
Total energy = 1.169917056
```

The run covered forced DLB steps 120 and 170.  No correctness anomaly was
observed, so the candidate was promoted to one 500-step formal run.

## 500-step formal

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_transition_cache_formal_20260607
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_transition_cache_500step_20260607
```

Formal comparison:

| metric | current forced8 baseline | transition cache | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 137.97 s | +4.14 s |
| move+collide wall | 97.10060854 s | 99.82211041 s | +2.72150187 s |
| move only | 72.9567669 s | 75.98126057 s | +3.02449367 s |
| buildCellOccupancy | 6.729516136 s | 6.750326008 s | +0.020809872 s |
| collision phase | 14.8398869 s | 14.07148121 s | -0.76840569 s |
| post fields/output | 40.37082178 s | 43.02448373 s | +2.65366195 s |
| full evolve wall | 132.7056254 s | 138.6072552 s | +5.9016298 s |
| particles per rank max/min | 1.412911358 | 1.797055588 | +0.38414423 |
| DLB rebalances | 8 | 8 | 0 |

Final check lines:

| metric | current forced8 baseline | transition cache |
| --- | ---: | ---: |
| Total Iterations | 500 | 500 |
| Number of DSMC particles | 2463567 | 2463716 |
| Collisions | 34522 | 34092 |
| Total energy | 1.239014073 | 1.238625899 |

## Decision

Reject and revert.

Reason:

- The candidate missed the strict full-fields formal gate.
- End-to-end `real` regressed by `+4.14 s`.
- `move only` regressed by `+3.02449367 s`.
- `post fields/output` and `full evolve wall` also regressed.
- The 200-step signal did not predict the 500-step behavior well enough to
  justify retaining the code.

This is another short/signal-stage candidate that does not survive formal
validation.  The transition table also added about `31.79 MiB` per rank of
static table/offset storage without producing a formal benefit.

## Revert

The transition-cache source and control switch were removed from:

- `particle.H`
- `particleTemplates.C`
- `dsmcCloud.H`
- `dsmcCloud.C`
- case `system/controlDict`

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_transition_cache_revert_build_20260607.log
```

Revert build result: passed, with only existing OpenFOAM template warnings.

Post-revert checks:

```text
system/controlDict hash restored to 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
system/fieldPropertiesDict hash remains 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
git diff for the candidate source files and controlDict is empty
no moveTetTransition / tetTransition / cachedTetNeighbour / cachedCrossEdge symbols remain
```

The retained same-tet non-normalised helper remains in `particleTemplates.C`.
