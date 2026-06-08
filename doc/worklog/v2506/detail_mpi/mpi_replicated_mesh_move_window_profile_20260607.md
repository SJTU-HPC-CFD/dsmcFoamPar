# MPI replicated-mesh move window profile - 2026-06-07

## Scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI8 replicated mesh, OpenMP disabled.
- Purpose: diagnostic-only window profile of `evolve_moveAndCollide()`.
- Formal controls were restored after the diagnostic runs.

Formal control hashes after restore:

```text
controlDict         8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
fieldPropertiesDict 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
```

## Diagnostic source

A temporary `moveWindowProfile` block was added to
`dsmcCloud::evolve_moveAndCollide()` to accumulate per-window wall time for:

- move+collide;
- move only;
- buildCellOccupancy;
- collision;
- local parcel count.

This block is diagnostic-only and is not retained as production source.

## First 500-step diagnostic

Backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_window_profile_20260607
```

Temporary controls:

```text
moveWindowProfile true;
moveWindowProfileInterval 500;
```

Log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_window_profile_500step_20260607
```

Result:

| metric | value |
| --- | ---: |
| real | 127.89 s |
| Total Iterations | 500 |
| OpenMP enabled | 0 |
| move+collide wall | 93.84286564 s |
| move only | 71.45264026 s |
| buildCellOccupancy | 6.419074039 s |
| collision phase | 13.56931447 s |
| post fields/output | 37.52243611 s |
| full evolve wall | 125.2529162 s |
| DLB checks / rebalances | 500 / 8 |
| final particles | 2463745 |
| stuck particles | 0 |
| final total energy | 1.238762473 |

Window output:

| steps | move+collide min | move+collide max | move+collide mean | move min | move max | move mean | build max | collision max | local particles min/max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1-500 | 87.7193194 | 93.84286564 | 91.40195048 | 65.34166537 | 71.45264026 | 69.05893555 | 6.419074039 | 13.56931447 | 232816 / 519394 |

Interpretation: this run is valid as a full-run rank aggregation, but it is not
a segmented profile.  The temporary source checked `replicatedMesh::rebalanceSteps()`,
whereas this case uses Phase C forced DLB steps from
`replicatedMeshDLBForceSteps`.  With `moveWindowProfileInterval 500`, the only
guaranteed split point was the end of the run.

## 50-step segmented diagnostic

Backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_window50_profile_20260607
```

Temporary controls:

```text
moveWindowProfile true;
moveWindowProfileInterval 50;
```

Log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_window50_profile_500step_20260607
```

Result:

| metric | value |
| --- | ---: |
| real | 140.40 s |
| Total Iterations | 500 |
| OpenMP enabled | 0 |
| OpenMP max threads | 1 |
| move+collide wall | 103.6679114 s |
| move only | 78.45228037 s |
| buildCellOccupancy | 7.909960787 s |
| collision phase | 15.52829072 s |
| post fields/output | 44.4835671 s |
| total profiled | 146.374099 s |
| full evolve wall | 138.4330278 s |
| migration calls / wall | 59 / 2.434911923 s |
| particles per rank | min 218420 max 433677 max/min 1.985518725 |
| rank wall time | min 84.42 max 93.71 max/min 1.110045013 |
| DLB checks / rebalances / triggered checks | 500 / 8 / 8 |
| final particles | 2463700 |
| stuck particles | 0 |
| final collisions / candidates | 35046 / 62552 |
| final total energy | 1.239344272 |

No `Fatal`, `Segmentation`, `Floating`, word-boundary `NaN`, or
`BAD TERMINATION` entries were found.

### Window table

All times are wall seconds.  Values are reduced across ranks as
`min/max/mean`.

| steps | move+collide max | move+collide mean | move max | move mean | build max | build mean | collision max | collision mean | local particles min/max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1-50 | 7.793032751 | 7.372890022 | 5.713631821 | 5.302939529 | 0.581164681 | 0.4921873732 | 1.062149551 | 0.9639238727 | 209445 / 321258 |
| 51-100 | 8.364809214 | 7.536651436 | 6.556180014 | 5.744892427 | 0.661837895 | 0.5552430188 | 1.24128694 | 1.101982373 | 218000 / 322006 |
| 101-150 | 11.19590263 | 10.1244859 | 8.646236638 | 7.57504912 | 0.82514423 | 0.601376227 | 1.701418763 | 1.509640491 | 123846 / 517051 |
| 151-200 | 10.87368173 | 10.15704402 | 8.314509534 | 7.604909176 | 0.944758897 | 0.6017394296 | 1.674952868 | 1.529395048 | 156561 / 410105 |
| 201-250 | 10.48077418 | 9.811495967 | 7.993546703 | 7.327554569 | 0.922790867 | 0.6313094349 | 1.651303029 | 1.471065632 | 206836 / 405472 |
| 251-300 | 11.89242099 | 11.20357358 | 9.067166302 | 8.371224749 | 1.045816202 | 0.723531134 | 1.89784883 | 1.718594808 | 219437 / 525580 |
| 301-350 | 11.88873386 | 10.73704434 | 9.160224189 | 7.998961262 | 0.950184807 | 0.75944577 | 1.798849468 | 1.571845527 | 242444 / 383428 |
| 351-400 | 10.99519041 | 9.877031835 | 8.414966302 | 7.294988527 | 1.00272472 | 0.7632289234 | 1.625477597 | 1.42131078 | 255957 / 364633 |
| 401-450 | 11.15935712 | 10.75597966 | 8.439029083 | 8.028162787 | 0.844672999 | 0.7362065206 | 1.667130313 | 1.554069651 | 209368 / 467460 |
| 451-500 | 11.693686 | 10.95212894 | 8.859273953 | 8.101872308 | 0.898102238 | 0.7709131763 | 1.795639152 | 1.62182772 | 214758 / 439741 |

Forced Phase C DLB triggers:

```text
120 170 220 270 320 370 420 470
```

## Interpretation

- The first 100 steps are materially cheaper than the later windows.  Move mean
  is `5.30-5.74 s` per 50 steps before the first forced DLB, then rises to
  `7.29-8.37 s` per 50 steps after the particle population grows.
- The worst move window by rank max is `301-350` at `9.160224189 s`; the worst
  move+collide windows are `251-300` and `301-350`.
- DLB reduces some rank-particle extremes, but large local-particle spread
  persists.  The `101-150` and `251-300` windows reach `517051` and `525580`
  local particles on the heaviest rank.
- buildCellOccupancy and collision also grow with particle count, but the
  dominant segment growth remains move.
- The diagnostic run itself is not an acceptance baseline.  It only localizes
  where the 500-step move cost accumulates.  Strict candidate acceptance remains
  the current formal full-fields baseline: improve both `real 133.83 s` and
  `move only 72.9567669 s`.

## Follow-up

The temporary `moveWindowProfile` source should be reverted after this note and
the DSMC library rebuilt so the retained tree is clean.
