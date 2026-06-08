# MPI replicated mesh move subprofile - 2026-06-07

## Scope

- Mode: pure MPI8 replicated mesh, no OpenMP.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Purpose: diagnostic-only split of the current `move` path below the existing
  aggregate `moveDetailProfile`.
- This is not a formal 500-step performance comparison.  The diagnostic patch
  adds counters in the tracking hot path and changes timing overhead.

Formal case control hash before and after the diagnostic runs:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

## Temporary source instrumentation

Added temporary counters under `moveDetailProfile true`:

- `same-tet fast`: retained non-normalised same-tet fast path returns.
- `tet-walk loops`: iterations of the internal tet-walk loop after the
  same-tet fast path fails.
- `findTris cand`: candidate triangle count produced by `findTris()`.
- `lambda tests`: second-stage `tetLambda(position_, endPosition, ...)`
  tests over the candidate triangles.
- `internal tet x`: internal tet transitions selected by `triI > 0`.
- `cell-face x`: cell-face crossings selected by `triI == 0`.
- `wall checks`: `hitWallFaces()` checks gated by wall-impact distance.
- `no-tri end`: non-fast calls that found no candidate triangle and returned
  directly to `endPosition`.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_subprofile_build_20260607.log
```

Build result: passed, with only the existing OpenFOAM template warnings.

## Case controls

Original formal `controlDict` backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_subprofile_smoke_20260607
```

10-step smoke temporary changes:

```text
endTime 1.e-06;
profileDetail true;
moveDetailProfile true;
```

200-step signal temporary changes:

```text
endTime 2.e-05;
profileDetail true;
moveDetailProfile true;
```

The formal `controlDict` was restored after the diagnostic runs.

## 10-step smoke

Log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_subprofile_smoke_10step_20260607
```

Summary:

| metric | value |
| --- | ---: |
| external real | 13.44 s |
| solver profile steps | 10 |
| move only | 1.327636155 s |
| move detail track max | 0.858761831 s |
| move detail boundary max | 0.001324257 s |
| move detail tracker max | 0 s |
| track calls | 25,758,317 |
| same-tet fast | 11,139,367 |
| internal tet only | 9,483,026 |
| face hits | 5,135,924 |
| tet-walk loops | 27,917,095 |
| findTris candidates | 20,316,952 |
| lambda tests | 20,316,952 |
| internal tet transitions | 13,298,145 |
| cell-face crossings | 5,135,924 |
| wall checks | 0 |
| no-tri end | 9,483,026 |

Derived:

| derived metric | value |
| --- | ---: |
| same-tet fast / track calls | 43.245710% |
| internal tet only / track calls | 36.815394% |
| face hits / track calls | 19.938896% |
| non-fast calls | 14,618,950 |
| tet-walk loops / non-fast call | 1.909651 |
| findTris candidates / tet-walk loop | 0.727760 |
| internal tet transitions / non-fast call | 0.909651 |

Correctness/check lines:

```text
Total Iterations = 10
Number of DSMC particles = 2065827
Collisions = 2864
Total energy = 1.1506003
```

## 200-step signal

Log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_subprofile_signal_200step_20260607
```

Summary:

| metric | value |
| --- | ---: |
| external real | 58.38 s |
| solver profile steps | 200 |
| move only | 30.19626312 s |
| move detail track max | 16.2060381 s |
| move detail boundary max | 0.048134089 s |
| move detail tracker max | 0 s |
| track calls | 533,640,621 |
| same-tet fast | 232,603,167 |
| internal tet only | 195,273,674 |
| face hits | 105,763,780 |
| tet-walk loops | 574,273,748 |
| findTris candidates | 417,477,672 |
| lambda tests | 417,477,672 |
| internal tet transitions | 273,236,294 |
| cell-face crossings | 105,763,780 |
| wall checks | 0 |
| no-tri end | 195,273,674 |

Derived:

| derived metric | value |
| --- | ---: |
| same-tet fast / track calls | 43.587980% |
| internal tet only / track calls | 36.592730% |
| face hits / track calls | 19.819290% |
| non-fast calls | 301,037,454 |
| tet-walk loops / non-fast call | 1.907649 |
| findTris candidates / tet-walk loop | 0.726966 |
| internal tet transitions / non-fast call | 0.907649 |
| track timer / move only | 53.669019% |

DLB/check lines:

```text
Phase C auto DLB rebalances = 2
Phase C auto DLB triggered checks = 2
particles per rank max/min = 2.191442775
rank wall time (evolve, 200 steps) max/min = 1.079306072
```

Correctness/check lines:

```text
Total Iterations = 200
Number of DSMC particles = 2220202
Collisions = 15056
Total energy = 1.169487359
```

## Interpretation

The 10-step and 200-step distributions match closely.  The retained same-tet
non-normalised helper accounts for all observed same-tet no-face returns:

```text
same-tet no-face == same-tet fast
```

That means the current same-tet path is already covered by the retained helper.
The next meaningful move space is not another same-tet micro-check.

The non-fast path is the main remaining local tracking workload:

- about 56% of track calls enter the tet-walk path;
- each non-fast call performs about 1.91 tet-walk loop iterations;
- wall checks are zero in this case;
- tracker time is zero and boundary time is tiny.

For the 200-step signal, `findTris()` internally evaluates four
`tetLambda(Ct, endPosition, ...)` calls per tet-walk loop.  The outer selection
then evaluates `tetLambda(position_, endPosition, ...)` over the candidates.
That is approximately:

```text
4*574,273,748 + 417,477,672 = 2,714,572,664 tetLambda-style evaluations
```

These are not direct duplicate calculations: `findTris()` uses the line from
the tet centre to `endPosition`, while the outer pass uses the actual parcel
segment from `position_` to `endPosition`.  Therefore the candidate lambdas
cannot simply be reused.

## Optimization direction

Do not spend more time on:

- wall-impact-distance bypasses;
- boundary/tracker paths;
- another same-tet normal/sqrt/outside-plane variant;
- same-face increment or `tetNeighbour(0)` early return;
- direct `Utracking` or `trackingData` cache changes.

The remaining plausible current-tree direction is a larger static-mesh data
layout change:

1. Static tet-neighbour transition table:
   `(cellI, tetFaceI, tetPtI, triI) -> next tetFace/tetPt or cell-face`.
2. Static tet-geometry cache:
   per-tet base point ids, oriented face point ids, plane base ids, and possibly
   normalised tet face areas.

Both must be generated from the existing mesh/tet logic, not from a new
heuristic.  The acceptance gate should be:

- 10-step detail counts exactly match baseline categories;
- 200-step signal shows no particle/collision/energy anomaly;
- only then consider 500-step formal.

The current diagnostic does not justify a low-risk one-line move optimization.
It narrows the next work to static tet transition/geometry layout or a broader
tracking API migration.

## Cleanup

The temporary internal counter patch was removed after the 10-step and 200-step
diagnostic runs.  The final rebuilt binary therefore does not retain the
additional `same-tet fast`, `tet-walk loops`, `findTris cand`, `lambda tests`,
`internal tet x`, `cell-face x`, `wall checks`, or `no-tri end` counters.

Cleanup build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_subprofile_cleanup_build_20260607.log
```

Cleanup build result: passed, with only the existing OpenFOAM template warnings.

Post-cleanup checks:

```text
system/controlDict hash restored to 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
system/fieldPropertiesDict hash remains 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
No dsmcMoveDetail* or temporary internal-counter symbols remain in the source.
```
