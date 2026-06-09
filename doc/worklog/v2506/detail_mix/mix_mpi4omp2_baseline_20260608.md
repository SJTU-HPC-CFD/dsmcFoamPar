# mix-mpi4omp2 baseline worklog - 2026-06-08

## Scope

User requested starting the MPI+OpenMP mixed-parallel phase.

Work directory:

```text
/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
```

Test case:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mix-mpi4omp2
```

Log directory:

```text
doc/worklog/v2506/detail_mix
```

## Initial Case Audit

The case name is `mix-mpi4omp2`, but the initial `controlDict` still used:

```text
useOpenMP false;
openmpThreads 1;
endTime 5.e-05;
profileSummary true;
profileDetail false;
replicatedMesh true;
replicatedMeshDLBForceSteps ();
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBImbalanceThreshold 1.5;
```

The case uses the replicated-mesh raw-MPI path, so the mixed baseline should be
run as `mpirun -np 4 dsmcFoam+` with `OMP_NUM_THREADS=2`, not as an OpenFOAM
decomposed `-parallel` run.

There were no existing `log*` files in this case and no matching `mix-*`
reference case/log under `ourmeshbkp`.

## Source Baseline

This mixed run inherits the current replicated-mesh source baseline after the
pure-MPI report:

```text
same-tet area reuse: active in particleTemplates.C
barycentric tracker port: inactive; backup files only under src/lagrangian/basic/particle/bkp
temporary trigger controls: reverted; no replicatedMeshDLBMinStartStep,
  replicatedMeshDLBRemainingGuardSteps, or replicatedMeshDLBUseSAR lookup remains
```

The active tracker path computes the initial tet face areas once, tests the
same-tet no-normalise fast path, and reuses the same `initialTetAreas` when the
particle has to enter the normal tet-walk path.

## Control Backups

Before changing `controlDict`, the original file was copied to:

```text
run/.../mix-mpi4omp2/system/controlDict.codex_mix_start_backup_20260608_103121
doc/worklog/v2506/detail_mix/mix_mpi4omp2_controlDict_before_smoke_20260608_103121
```

## Smoke Configuration

Temporary 10-step smoke settings:

```text
useOpenMP true;
openmpThreads 2;
endTime 1.e-06;
deltaT 1.e-07;
profileSummary true;
profileDetail false;
```

Run form:

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/env.sh
cd run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mix-mpi4omp2
OMP_NUM_THREADS=2 /usr/bin/time -p mpirun -np 4 dsmcFoam+
```

## First Smoke Result

Log:

```text
doc/worklog/v2506/detail_mix/stage_mix_mpi4omp2_smoke_np4_omp2_20260608.log
```

Result: failed with exit code `255`.

The run reached replicated-mesh initialisation:

```text
Replicated mesh: METIS decomposition complete, 104151 cells -> 4 ranks
Phase C auto DLB: enabled
Replicated mesh: initialized with 4 MPI ranks, migrate interval 10
```

The crash happened during the first OpenMP move step.  The stack pointed to
`dsmcPatchBoundary::measurePropertiesBeforeControl()` called from
`dsmcDiffuseWallPatch::controlParticle()` inside `Cloud<dsmcParcel>::move()`.

## Mixed OMP Boundary Fix

`src/lagrangian/dsmc/parcels/dsmcParcel.C` was changed conservatively so that
all patch-boundary `controlParticle()` calls enter the existing
`dsmcMoveBoundary` OpenMP critical section.  The previous special case allowed
`dsmcDiffuseWallPatch` and `dsmcSpecularWallPatch` to run concurrently, but wall
patches update shared wall measurements and the diffuse wall path uses the
shared cloud RNG.

Build log:

```text
doc/worklog/v2506/detail_mix/stage_mix_boundary_critical_build_20260608.log
```

Result: build passed.

## Passing Smoke Result

Log:

```text
doc/worklog/v2506/detail_mix/stage_mix_mpi4omp2_smoke_np4_omp2_after_boundary_critical_20260608.log
```

Result: passed with exit code `0`.

Key 10-step results:

| metric | value |
|---|---:|
| real | 8.24 s |
| total iterations | 10 |
| OpenMP max threads | 2 |
| collisions | 2845 |
| collision candidates | 2935 |
| DSMC particles | 2065852 |
| stuck particles | 0 |
| total energy | 1.150612748 |
| Phase C auto DLB checks | 10 |
| Phase C auto DLB rebalances | 0 |

## Formal 500-Step Baseline

Before switching from smoke to formal, the smoke `controlDict` was backed up to:

```text
run/.../mix-mpi4omp2/system/controlDict.codex_mix_smoke_before_formal_20260608_105248
doc/worklog/v2506/detail_mix/mix_mpi4omp2_controlDict_before_formal_20260608_105248
```

Formal `controlDict` copy:

```text
doc/worklog/v2506/detail_mix/mix_mpi4omp2_controlDict_formal_500step_20260608_105248
```

Formal setting:

```text
useOpenMP true;
openmpThreads 2;
endTime 5.e-05;
deltaT 1.e-07;
replicatedMesh true;
replicatedMeshAutoDLB true;
```

Formal log:

```text
doc/worklog/v2506/detail_mix/stage_mix_mpi4omp2_formal_np4_omp2_500step_20260608.log
```

Result: passed with exit code `0`.

Key 500-step results:

| metric | value |
|---|---:|
| real | 104.65 s |
| total iterations | 500 |
| OpenMP max threads | 2 |
| collisions | 35271 |
| collision candidates | 65578 |
| DSMC particles | 2463810 |
| stuck particles | 0 |
| total energy | 1.241652653 |
| move+collide wall | 81.20954127 s |
| move only | 68.54304325 s |
| buildCellOccupancy | 7.287947673 s |
| collision phase | 4.94542299 s |
| post fields/output | 25.13833162 s |
| full evolve wall | 100.2682179 s |
| migration calls | 54 |
| migration wall time | 1.360225229 s |
| particles per rank max/min | 1.191738949 |
| rank wall time max/min | 1.083234328 |
| Phase C auto DLB checks | 500 |
| Phase C auto DLB rebalances | 3 |
| Phase C triggered checks | 3 |
| Phase C auto DLB wall max | 0.98730365 s |

No fatal error, segmentation fault, MPI abort, or NaN marker was found in the
passing smoke or formal baseline logs.  `sigFpe` lines are the normal OpenFOAM
floating-point exception trapping banner.

The case is currently left in the formal mixed baseline state
(`useOpenMP true`, `openmpThreads 2`, `endTime 5.e-05`).

## Position Versus Existing Baselines

This is an orientation comparison only; it is not a same-run A/B sweep.

| mode | source/control label | real | full evolve | move only | collision | post |
|---|---|---:|---:|---:|---:|---:|
| OMP8 | stage25 valid OMP | 90.34 s | 90.00671836 s | 49.24385297 s | 3.376808391 s | 30.39917036 s |
| pure MPI8 replicated | same-tet area reuse formal | 104.32 s | 104.4499775 s | 67.43927759 s | 14.09394481 s | 17.64348673 s |
| pure MPI8 replicated | same-tet area reuse repeat | 107.45 s | 105.2572367 s | 67.84144046 s | 14.86660294 s | 17.2290821 s |
| mixed 4MPI x 2OMP | this baseline, auto DLB | 104.65 s | 100.2682179 s | 68.54304325 s | 4.94542299 s | 25.13833162 s |

Interpretation:

- mixed 4x2 starts correctly and improves the collision phase versus pure MPI8;
- mixed move time remains essentially pure-MPI-like, so the OpenMP move path is
  not yet converting 2 threads/rank into a useful move speedup;
- mixed post/output is lower than OMP8 but higher than pure-MPI replicated
  same-tet, which is consistent with combining OpenMP field/post paths with the
  replicated mesh ownership model;
- mixed auto DLB performed 3 rebalances in this run, while the accepted pure-MPI
  area-reuse formal pair used 8 rebalances in the forced schedule.  DLB policy
  is therefore another comparison axis and should not be collapsed into the
  OpenMP-vs-MPI effect.

Next useful gates:

1. Run a 200-step mixed diagnostic with `profileDetail true` and
   `moveDetailProfile true` to split mixed move into tracking, boundary, and
   append/ordered-list costs after the new boundary critical.
2. If boundary dominates, replace the coarse boundary critical with a real
   thread-local path: thread-local wall RNG plus per-thread boundary measurement
   accumulation/reduction.
3. If tracking dominates, keep using the current same-tet/tet-walk path and
   avoid the barycentric branch unless a new correctness/performance gate
   justifies re-opening it.
4. Keep schedule defaults at `openmpMoveSchedule static`,
   `openmpMoveChunk 64`, `openmpCollisionSchedule dynamic`,
   `openmpCollisionChunk 8` until a mixed-specific sweep beats them on this
   same `mix-mpi4omp2` case.

## 200-Step Move-Detail Diagnostic

The first 200-step diagnostic completed, but did not print move-detail counters:

```text
doc/worklog/v2506/detail_mix/stage_mix_mpi4omp2_detail_np4_omp2_200step_20260608.log
```

Reason: the OpenMP `Cloud::move()` path creates per-thread `localTd` objects.
The diagnostic flag `moveDetailProfile` was set on the outer `td`, but was not
copied into `localTd`, and local counters were not accumulated back into the
outer `td`.

Instrumentation fix:

- `src/lagrangian/basic/Cloud/Cloud.C`
- added no-op-safe `cloudOpenMP` helpers to copy/accumulate move-detail fields
  when the `TrackData` type supports them;
- copied `moveDetailProfile` into per-thread `localTd`;
- accumulated per-thread move-detail counters at the end of both OpenMP move
  branches.

Build log:

```text
doc/worklog/v2506/detail_mix/stage_mix_move_detail_omp_aggregation_build_20260608.log
```

Result: build passed.

Diagnostic control backups:

```text
doc/worklog/v2506/detail_mix/mix_mpi4omp2_controlDict_before_detail_200step_20260608_105744
doc/worklog/v2506/detail_mix/mix_mpi4omp2_controlDict_detail_200step_20260608_105744
doc/worklog/v2506/detail_mix/mix_mpi4omp2_controlDict_after_detail_restore_20260608_105744
```

The case was restored to the formal 500-step mixed baseline control after the
diagnostic.

Rerun log with valid move-detail counters:

```text
doc/worklog/v2506/detail_mix/stage_mix_mpi4omp2_detail_np4_omp2_200step_after_aggregation_20260608.log
```

Result: passed with exit code `0`.

Key 200-step diagnostic results:

| metric | value |
|---|---:|
| real | 45.11 s |
| total iterations | 200 |
| OpenMP max threads | 2 |
| collisions | 15693 |
| DSMC particles | 2220195 |
| stuck particles | 0 |
| total energy | 1.169994588 |
| move+collide wall | 31.8647929 s |
| move only | 27.47583734 s |
| buildCellOccupancy | 2.280479583 s |
| collision phase | 1.211962116 s |
| post fields/output | 8.968296275 s |
| move detail parcels | 428039784 |
| move detail track calls | 533628552 |
| same-tet no-face | 232615890 |
| internal tet only | 195261685 |
| face hits | 105750977 |
| patch hits | 438401 |
| move detail track max | 29.28668551 s |
| move detail boundary max | 0.100746499 s |
| Phase C auto DLB rebalances | 2 |

Diagnostic interpretation:

- `boundary max` is tiny compared with `move only`, so the conservative
  boundary critical is not the dominant mixed bottleneck in this case;
- `track calls`, `same-tet no-face`, and `internal tet only` dominate the move
  profile, keeping mixed move time close to the pure-MPI replicated path;
- the next source-level work should stay on the active Cartesian
  same-tet/tet-walk path.  Re-opening the barycentric branch is not justified by
  this diagnostic.
