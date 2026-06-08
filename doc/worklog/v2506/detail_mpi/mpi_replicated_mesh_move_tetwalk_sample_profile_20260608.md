# MPI replicated mesh tet-walk sampled profile - 2026-06-08

## Scope

- Mode: pure MPI8 replicated mesh DLB only, no OpenMP evidence.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Purpose: diagnostic-only split of the non-fast tet-walk loop after the
  fixed tet-geometry cache feasibility audit.
- This is not a formal 500-step performance comparison.  The temporary patch
  adds sampled timers in the tracking hot path and changes overhead.

Formal case hashes after restoring controls:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

## Temporary Instrumentation

Added diagnostic-only fields to `particle::TrackingData` and DSMC profile
plumbing:

- `tetWalkLoops`
- `tetWalkSamples`
- sampled `tetNeighbour(triI)` time
- sampled metadata reconstruction time
- sampled area-normal plus plane-base setup time
- sampled `findTris()` time
- sampled outer candidate `tetLambda()` selection time

Control switch:

```text
profileDetail true;
moveDetailProfile true;
moveTetWalkDetailProfile true;
moveTetWalkDetailSampleStride 1024;
```

The sampling is every roughly 1024 tet-walk loop iterations.  The printed
`scaled` values are aggregate rank-time estimates produced from summed samples
and summed loop counts.  They are useful for relative section ranking, but they
are not the same quantity as `move detail track max [s]`, which is a max-rank
wall timer.

Build logs:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_tetwalk_sample_profile_build_20260608.log
doc/worklog/v2506/detail_mpi/stage_mpi_move_tetwalk_sample_profile_build_retry_20260608.log
```

The first build failed because the local `tetWalkDetailProfile` variable was
inserted into the non-DSMC `trackToFace()` overload while the use site was in
the DSMC overload.  The retry build passed and ended with `=== Done ===`.

## 10-step Smoke

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_tetwalk_sample_profile_smoke_20260608
```

Temporary changes:

```text
endTime 1.e-06;
profileDetail true;
moveDetailProfile true;
moveTetWalkDetailProfile true;
moveTetWalkDetailSampleStride 1024;
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_tetwalk_sample_profile_smoke_10step_20260608
```

Summary:

| metric | value |
| --- | ---: |
| external real | 9.20 s |
| solver profile steps | 10 |
| move only | 1.369392983 s |
| move detail track max | 0.904237983 s |
| move detail boundary max | 0.00138664 s |
| move detail tracker max | 0 s |
| track calls | 25,758,232 |
| same-tet no-face | 11,139,301 |
| internal tet only | 9,483,161 |
| face hits | 5,135,770 |
| tet-walk loops | 27,917,330 |
| tet-walk samples | 27,224 |
| sample scale | 1025.467602 |

Sampled scaled split:

| section | scaled aggregate rank-time | share of sampled split |
| --- | ---: | ---: |
| tetNeighbour | 2.295453852 s | 23.310002% |
| metadata reconstruction | 0.9022658735 s | 9.162380% |
| area normals + plane bases | 1.675107481 s | 17.010474% |
| findTris | 3.936139462 s | 39.970928% |
| outer lambda selection | 1.038539238 s | 10.546216% |

Correctness/check lines:

```text
Total Iterations = 10
Number of DSMC particles = 2065856
Collisions = 2803
Total energy = 1.150631787
```

## 200-step Signal

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_move_tetwalk_sample_profile_signal_20260608
```

Temporary changes:

```text
endTime 2.e-05;
profileDetail true;
moveDetailProfile true;
moveTetWalkDetailProfile true;
moveTetWalkDetailSampleStride 1024;
```

Run log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_move_tetwalk_sample_profile_signal_200step_20260608
```

Summary:

| metric | value |
| --- | ---: |
| external real | 58.53 s |
| solver profile steps | 200 |
| move+collide wall | 39.17950336 s |
| move only | 30.93683107 s |
| buildCellOccupancy | 2.434135007 s |
| collision phase | 5.170508854 s |
| post fields/output | 17.36747305 s |
| full evolve wall | 53.75868322 s |
| move detail track max | 15.53090221 s |
| move detail boundary max | 0.045771432 s |
| move detail tracker max | 0 s |
| track calls | 533,643,543 |
| same-tet no-face | 232,602,591 |
| internal tet only | 195,276,408 |
| face hits | 105,764,544 |
| processor hits | 0 |
| cyclic hits | 0 |
| patch hits | 439,790 |
| tet-walk loops | 574,273,170 |
| tet-walk samples | 560,016 |
| sample scale | 1025.458505 |

Sampled scaled split:

| section | scaled aggregate rank-time | share of sampled split |
| --- | ---: | ---: |
| tetNeighbour | 47.46400934 s | 24.881215% |
| metadata reconstruction | 18.06732779 s | 9.471115% |
| area normals + plane bases | 31.07408760 s | 16.289417% |
| findTris | 75.57627232 s | 39.618008% |
| outer lambda selection | 18.58072770 s | 9.740245% |

DLB/check lines:

```text
Phase C auto DLB checks = 200
Phase C auto DLB rebalances = 2
Phase C auto DLB triggered checks = 2
particles per rank max/min = 1.90216623
rank wall time (evolve, 200 steps) max/min = 1.077318031
```

Correctness/check lines:

```text
Total Iterations = 200
Number of DSMC particles = 2220162
Collisions = 15307
Total energy = 1.169617757
```

## Interpretation

The 10-step and 200-step sampled splits are consistent:

- `findTris()` is the largest sampled section at about `40%` of the sampled
  tet-walk split.
- `tetNeighbour()` and edge-connected topology transition work is about
  `23-25%`.
- area-normal construction plus normalisation is only about `16-17%`.
- metadata reconstruction is under `10%`.
- outer candidate lambda selection is about `10%`.

This does not support a full per-state static geometry cache as the next
candidate.  Such a cache would mainly target the area-normal plus plane-base
setup section, which is smaller than `findTris()` and topology transition work.
The earlier feasibility note also estimated about `133.49 MiB/rank` for
normals plus plane bases, before optional vertex-label storage.

The result points away from a narrow fixed geometry cache and toward one of two
larger directions:

1. Reduce or restructure `findTris()` / lambda screening itself.
2. Audit the broader barycentric/tracking API migration path rather than
   continuing current-tree micro-caches.

## Cleanup

The temporary instrumentation was removed after this note and the solver was
rebuilt.

Cleanup build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_tetwalk_sample_profile_cleanup_build_20260608.log
```

Cleanup build result: passed, with only the existing OpenFOAM template
warnings.  The retained log ends with `=== Done ===`.

Post-cleanup checks:

```text
system/controlDict hash restored to 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
system/fieldPropertiesDict hash remains 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
no tetWalkDetail / moveTetWalkDetail / tetWalkSample / dsmcTrackElapsedSeconds symbols remain
```

The retained same-tet non-normalised helper remains in
`particleTemplates.C`.
