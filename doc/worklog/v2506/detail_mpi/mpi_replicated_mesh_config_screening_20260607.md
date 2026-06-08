# MPI replicated mesh config screening - 2026-06-07

Scope:

- pure MPI8 replicated mesh DLB only;
- OpenMP remains disabled;
- formal comparison target remains forced8 with controlDict hash
  `8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2`.

## Delayed receive / flat transfer path check

Current formal controls include:

```text
replicatedMeshDelayedReceive true;
replicatedMeshNoAlltoall true;
replicatedMeshFlatTransfer true;
```

Source check:

- `dsmcCloud.C` calls `migrateBegin(); migrateFinish();` when
  `replicatedMeshDelayedReceive true`.
- `dsmcReplicatedMesh::migrateBegin()` immediately falls back to
  `migrateParticlesByCellOwner()` when `useFlatTransfer_` is true.
- Therefore, with the current accepted `replicatedMeshFlatTransfer true`,
  `replicatedMeshDelayedReceive true` is not an overlap path and does not make
  `replicatedMeshNoAlltoall true` active.
- Turning `replicatedMeshFlatTransfer false` would route to the stream transfer
  fallback.  Earlier port validation recorded that the stream fallback failed on
  OFv1706 binary parcel reads, so this is not a low-risk config candidate.

Decision: do not spend a formal run on delayedReceive/noAlltoall alone while
flat transfer is enabled.

## Candidate: `replicatedMeshDLBVsizeExp 0.25`

Rationale:

- current forced8 uses `replicatedMeshDLBVsizeExp 0`, so ParMETIS sees unit
  migration sizes for 4+ ranks;
- logs show very large cell-owner changes at forced DLB steps, often around
  70-98% of cells;
- `vsizeExp=0.25` should penalize particle-heavy cell movement mildly without
  changing DSMC math.

Temporary config:

```text
replicatedMeshDLBVsizeExp 0.25;
```

Control backup before the candidate:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_vsize025_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

### 10-step smoke

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_vsize025_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 12.10 s |
| move only | 1.138835794 s |
| full evolve wall | 2.633190181 s |
| OpenMP enabled | 0 |
| DLB rebalances | 0 |
| final particles | 2065859 |
| stuck particles | 0 |
| final collisions | 2792 |
| final candidates | 2930 |
| final total energy | 1.150623258 |

Smoke passed; no Fatal/NaN/BAD TERMINATION string was found.

### 200-step DLB signal test

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_vsize025_smoke_200step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 200 |
| real | 55.00 s |
| move only | 28.3369108 s |
| buildCellOccupancy | 2.607412446 s |
| collision phase | 5.402109378 s |
| post fields/output | 17.03748328 s |
| full evolve wall | 49.89051853 s |
| OpenMP enabled | 0 |
| DLB rebalances | 2 |
| changed cells at step 120 | 72410 / 104151 = 0.6952405642 |
| changed cells at step 170 | 71549 / 104151 = 0.6869737208 |
| Phase C ParMETIS max | 0.20446931 s |
| Phase C migration max | 0.551284359 s |
| particles max/min | 2.753409602 |
| rank wall max/min | 1.127638191 |
| final particles | 2220342 |
| stuck particles | 0 |
| final collisions | 15366 |
| final candidates | 25696 |
| final total energy | 1.170357806 |

Interpretation:

- `vsizeExp=0.25` reduced DLB movement and DLB accounting cost.
- The same run also left worse particle/rank imbalance and did not provide a
  strong move-path signal.
- Since current accepted forced8 already spends only about 2 s in total Phase C
  DLB accounting over 500 steps, this candidate has too small a likely upside
  and too much imbalance risk for a 500-step formal run.

Decision: reject `replicatedMeshDLBVsizeExp 0.25` for now; retain forced8
`replicatedMeshDLBVsizeExp 0`.

The production `controlDict` was restored after the test:

```text
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

## Baseline contract and `replicatedMeshMigrateInterval`

Two comparison layers must stay separate in later notes:

- preserved historical reference:
  `ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610`
  reports `External wall seconds 117.31`, older detailed-profile
  `move only 47.1652636 s`, `migration calls 59`, `Phase C auto DLB checks 500`,
  and `Phase C auto DLB rebalances 8`;
- current same-tree forced8 repeats use the current profile-summary path and
  current retained source state.  The best near-time repeat for strict config
  screening remains
  `log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607`:
  `real 133.83 s`, `move only 72.9567669 s`, `DLB rebalances 8`.

The historical reference remains useful for the larger performance gap, but its
older detailed profile is not a like-for-like timer for every current profile
field.  Candidate retention in this config-screening round should therefore
beat the current same-tree forced8 run first, then report the remaining gap to
the preserved `ourmeshbkp` reference separately.

Source check for `replicatedMeshMigrateInterval`:

- `dsmcReplicatedMesh.C` still has a stale comment claiming every-step migration
  is required, but the implementation reads
  `replicatedMeshMigrateInterval` from `controlDict`;
- `dsmcCloud.C` uses it only to gate regular particle migration after move:
  step 0, then `stepCounter % migrateInterval == 0`;
- this is distinct from the DLB cadence: current runs still call the
  auto-rebalance entry each step and report `Phase C auto DLB checks 500`, while
  actual forced repartitions are the eight configured DLB force steps.

History check:

- older replicated-mesh tuning already scanned interval values on this case
  family and retained `replicatedMeshMigrateInterval 10`;
- the 2026-06-01 matrix recorded dual-constraint `interval=1` around `138 s`,
  `interval=2` around `130 s`, `interval=5` around `125 s`, and
  `interval=10` around `121 s`;
- current `detail_mpi` evidence has only interval-1 10-step smoke logs, not a
  current pure-MPI forced8 formal run.

Decision: do not spend another formal run on `replicatedMeshMigrateInterval`
right now.

Rationale:

- decreasing the interval would increase regular migration frequency, and old
  data already ranked it below interval 10;
- increasing the interval can only save a small regular migration budget in the
  current retained case, while it risks stale ownership/occupancy effects;
- current 500-step forced8 logs spend only about `1.7-2.2 s` in migration-wall
  summary, whereas the current move path is `72-75 s`, so the maximum plausible
  upside is too small for the risk.

Retain the current production setting:

```text
replicatedMeshMigrateInterval 10;
```

The production `controlDict` remains at forced8 hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

## Candidate: `replicatedMeshDLBItr 100`

Rationale:

- current forced8 writes `replicatedMeshDLBItr 1000`;
- source default is `100`;
- `itr` is passed only to `ParMETIS_V3_AdaptiveRepart`, so this is a pure DLB
  partitioning tradeoff test and does not change DSMC tracking/collision/post
  math.

Temporary config:

```text
replicatedMeshDLBItr 100;
```

Control backup before the candidate:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_itr100_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

### 10-step smoke

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_itr100_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 10.50 s |
| move only | 1.300827125 s |
| full evolve wall | 2.80487811 s |
| OpenMP enabled | 0 |
| DLB rebalances | 0 |
| final particles | 2065843 |
| stuck particles | 0 |
| final collisions | 2830 |
| final candidates | 2932 |
| final total energy | 1.150616813 |

Smoke passed; no Fatal/NaN/BAD TERMINATION string was found.

### 200-step DLB signal test

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_itr100_smoke_200step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 200 |
| real | 54.93 s |
| move only | 25.19764393 s |
| buildCellOccupancy | 2.43696367 s |
| collision phase | 4.780814507 s |
| post fields/output | 17.1258495 s |
| full evolve wall | 47.81540611 s |
| OpenMP enabled | 0 |
| DLB rebalances | 2 |
| changed cells at step 120 | 101363 / 104151 = 0.973231174 |
| changed cells at step 170 | 80831 / 104151 = 0.7760943246 |
| Phase C ParMETIS max | 0.174268232 s |
| Phase C migration max | 0.34043422 s |
| particles max/min | 2.74916965 |
| rank wall max/min | 1.083991385 |
| final particles | 2220161 |
| stuck particles | 0 |
| final collisions | 15117 |
| final candidates | 25384 |
| final total energy | 1.169867889 |

The 200-step signal was not strong, but it had enough DLB-accounting reduction
and rank-wall balance to justify one 500-step formal run.

### 500-step formal

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_itr100_500step_20260607
```

Result:

| metric | current forced8 better repeat | `DLBItr=100` | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 133.85 s | +0.02 s |
| move only | 72.9567669 s | 74.02791549 s | +1.07114859 s |
| move+collide wall | 97.10060854 s | 98.23853773 s | +1.13792919 s |
| full evolve wall | 132.7056254 s | 135.3508952 s | +2.6452698 s |
| post fields/output | 40.37082178 s | 44.19082667 s | +3.82000489 s |
| DLB rebalances | 8 | 8 | 0 |
| Phase C ParMETIS max | 0.68273936 s | 0.72637036 s | +0.043631 s |
| Phase C migration max | 1.687160342 s | 1.747275407 s | +0.060115065 s |
| particles max/min | 1.412911358 | 2.325948261 | +0.913036903 |
| rank wall max/min | 1.053059781 | 1.082145034 | +0.029085253 |

Correctness:

- `Total Iterations = 500`;
- `OpenMP enabled = 0`, `OpenMP max threads = 1`;
- final particles `2463712`, stuck `0`;
- final collisions `34380`, candidates `62113`, acceptance `0.5535073173`;
- final total energy `1.238712934`;
- no Fatal/NaN/BAD TERMINATION string was found.

Decision: reject `replicatedMeshDLBItr 100` and restore forced8
`replicatedMeshDLBItr 1000`.  The 500-step formal missed the strict target and
regressed the move path, full evolve wall, post/output time, and final particle
balance.

The production `controlDict` was restored after the test:

```text
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```
