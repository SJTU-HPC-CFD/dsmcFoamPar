# MPI replicated mesh move tetNeighbour cached refs rejection - 2026-06-07

## Scope

- Mode: pure MPI8 replicated mesh, no OpenMP.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Formal comparison target: current strict full-fields forced8 baseline:
  `log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607`

## Candidate

Files:

- `src/lagrangian/basic/particle/particle.H`
- `src/lagrangian/basic/particle/particleI.H`
- `src/lagrangian/basic/particle/particleTemplates.C`

Idea:

- Add a `tetNeighbour(triI, faceOwner, faces, tetBasePtIs)` overload.
- Keep the original `tetNeighbour(triI)` as a wrapper for ordinary callers.
- In the DSMC `trackToFace(..., true)` overload only, call the new overload
  with the already cached `faceOwner`, `pFaces`, and `tetBasePtIs` references.

This was a data-layout candidate only.  It did not intentionally change the tet
geometry, lambda tests, same-tet fast path, or boundary handling.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_tetneigh_cached_refs_build_20260607.log
```

Build result: passed, with only existing OpenFOAM template warnings.

## 10-step smoke

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_tetneigh_cached_refs_smoke_20260607
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_tetneigh_cached_refs_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| real | 7.85 s |
| move only | 1.04021043 s |
| move+collide wall | 1.554254451 s |
| full evolve wall | 2.324689906 s |
| post fields/output | 0.887244899 s |
| particles per rank max/min | 1.638504107 |
| DLB rebalances | 0 |

Correctness/check lines:

```text
Total Iterations = 10
Number of DSMC particles = 2065852
Collisions = 2928
Total energy = 1.150631355
```

The smoke result was positive enough to run a 200-step signal.

## 200-step signal

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_tetneigh_cached_refs_signal_20260607
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_tetneigh_cached_refs_signal_200step_20260607
```

Result:

| metric | value |
| --- | ---: |
| real | 58.45 s |
| move only | 27.67124267 s |
| move+collide wall | 36.78781236 s |
| full evolve wall | 51.15396182 s |
| post fields/output | 17.16831739 s |
| particles per rank max/min | 1.692796917 |
| DLB rebalances | 2 |

Correctness/check lines:

```text
Total Iterations = 200
Number of DSMC particles = 2220208
Collisions = 15125
Total energy = 1.169917014
```

The 200-step run did not show a correctness anomaly.  Because the 10-step move
signal was strong and the 200-step move value did not obviously regress, the
candidate was promoted to one 500-step formal validation.

## 500-step formal

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_tetneigh_cached_refs_formal_20260607
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_tetneigh_cached_refs_500step_20260607
```

Formal comparison:

| metric | current forced8 baseline | cached refs | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 135.50 s | +1.67 s |
| move+collide wall | 97.10060854 s | 100.3172837 s | +3.21667516 s |
| move only | 72.9567669 s | 76.02605933 s | +3.06929243 s |
| buildCellOccupancy | 6.729516136 s | 7.250746091 s | +0.521229955 s |
| collision phase | 14.8398869 s | 15.07103127 s | +0.23114437 s |
| post fields/output | 40.37082178 s | 43.33847728 s | +2.96765550 s |
| full evolve wall | 132.7056254 s | 136.6747296 s | +3.96910420 s |
| particles per rank max/min | 1.412911358 | 2.285750483 | +0.872839125 |
| DLB rebalances | 8 | 8 | 0 |

Final check lines:

| metric | baseline | cached refs |
| --- | ---: | ---: |
| Total Iterations | 500 | 500 |
| Number of DSMC particles | 2,463,567 | 2,463,765 |
| Collisions | 34,522 | 34,662 |
| Total energy | 1.239014073 | 1.239225394 |

## Decision

Reject and revert.

Reason:

- The formal run missed the strict acceptance gate.
- Both end-to-end `real` and `move only` regressed.
- The final particles-per-rank ratio was materially worse, even though the
  number of forced DLB rebalances was unchanged.

The short-run improvement did not survive formal validation.  This is another
case where 10-step smoke and 200-step signal were insufficient to justify
retaining a move micro-candidate.

## Revert

The cached-ref candidate was removed from:

- `particle.H`
- `particleI.H`
- `particleTemplates.C`

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_tetneigh_cached_refs_revert_build_20260607.log
```

Revert build result: passed, with only existing OpenFOAM template warnings.

Post-revert checks:

```text
system/controlDict hash restored to 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
system/fieldPropertiesDict hash remains 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
git diff for the three candidate source files and controlDict is empty
```

The retained same-tet non-normalised helper remains in `particleTemplates.C`.
