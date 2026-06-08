# MPI replicated mesh forced DLB steps - 2026-06-06

## Scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI, `useOpenMP false`, `openmpThreads 1`
- Formal run length: 500 steps (`endTime 5.e-05`, `deltaT 1.e-07`)
- Current acceptance gate: improve both current MPI-only baseline `real 135.30 s` and `move only 72.59819434 s`.
- Final target remains `ourmeshbkp/mpi8replicatedmesh`, not the older OpenMP2 113 s run.

## References

Current MPI-only baseline after post113 rollback:

- Log: `log.codex_mpi8_replicatedmesh_mpionly_post113_revert_500step_20260606`
- `real 135.30 s`
- `move only 72.59819434 s`
- `full evolve wall 125.5070652 s`
- `buildCellOccupancy 7.432408543 s`
- `collision phase 12.97955973 s`
- `post fields/output 38.10104302 s`
- `Phase C auto DLB checks 500`
- `Phase C auto DLB rebalances 5`

Preserved `ourmeshbkp` reference:

- Log: `run/.../ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610`
- `External wall seconds 117.31`
- `move only 47.1652636 s`
- `buildCellOccupancy 10.7084157 s`
- `collision phase 10.89824478 s`
- `post fields/output 30.35563996 s`
- `Phase C auto DLB checks 500`
- `Phase C auto DLB rebalances 8`

## Candidate A: force reference-like DLB steps

Configuration delta:

```text
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Run:

```text
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_forced_dlb_steps_500step_20260606 2>&1
```

Result:

| metric | current MPI baseline | forced 8-step candidate | delta |
| --- | ---: | ---: | ---: |
| real | 135.30 s | 132.29 s | -3.01 s |
| move only | 72.59819434 s | 69.49945105 s | -3.09874329 s |
| full evolve wall | 125.5070652 s | 123.0267073 s | -2.4803579 s |
| buildCellOccupancy | 7.432408543 s | 6.216832977 s | -1.215575566 s |
| collision phase | 12.97955973 s | 12.98068576 s | +0.00112603 s |
| post fields/output | 38.10104302 s | 36.15433245 s | -1.94671057 s |
| DLB rebalances | 5 | 8 | +3 |

Correctness checks from the log:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- Final particles: `2463524`, stuck `0`
- Final collisions: `34199`, candidates `60209`, acceptance `0.5680047833`
- Final total energy: `1.238930476`
- No Fatal/NaN match in the log.

Decision: keep this configuration as the current best MPI-only candidate. It improves both move and end-to-end wall time versus the current source baseline, but it is still far slower than `ourmeshbkp` (`117.31 s`, `move 47.1652636 s`).

## Candidate B: dense forced DLB steps

Configuration delta:

```text
replicatedMeshDLBForceSteps (90 120 150 180 210 240 270 300 330 360 390 420 450);
```

Run:

```text
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_dense_forced_dlb_steps_500step_20260606 2>&1
```

Result:

| metric | forced 8-step retained | dense forced candidate | delta |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 134.81 s | +2.52 s |
| move only | 69.49945105 s | 71.78133755 s | +2.28188650 s |
| full evolve wall | 123.0267073 s | 125.5469243 s | +2.5202170 s |
| buildCellOccupancy | 6.216832977 s | 6.206141795 s | -0.010691182 s |
| collision phase | 12.98068576 s | 13.19992068 s | +0.21923492 s |
| post fields/output | 36.15433245 s | 36.42859915 s | +0.27426670 s |
| DLB rebalances | 8 | 13 | +5 |

Correctness checks from the log:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- Final particles: `2463865`, stuck `0`
- Final collisions: `34469`, candidates `61856`, acceptance `0.5572458614`
- Final total energy: `1.238143942`
- No Fatal/NaN match in the log.

Decision: reject and revert. More frequent forced rebalancing improved final particle-count balance (`max/min 1.942` vs `2.953`) but worsened move and total wall time. The extra rebalances increased migration/DLB overhead and did not reduce the critical move path.

## Current retained case state

`system/controlDict` is restored to Candidate A:

```text
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Next candidate should stay in the DLB/owner-distribution area, but simple force-step densification is not useful. A more targeted test is to change the ParMETIS objective, for example testing move-only `replicatedMeshDLBDualConstraint false` with the retained 8 forced steps.

## Candidate C: move-only ParMETIS constraint with retained 8 forced steps

Configuration delta:

```text
replicatedMeshDLBDualConstraint false;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Run:

```text
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_moveonly_constraint_forced8_500step_20260606 2>&1
```

Result:

| metric | forced 8-step retained | move-only candidate | delta |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 136.75 s | +4.46 s |
| move only | 69.49945105 s | 75.65805056 s | +6.15859951 s |
| full evolve wall | 123.0267073 s | 127.3023569 s | +4.2756496 s |
| buildCellOccupancy | 6.216832977 s | 6.192013971 s | -0.024819006 s |
| collision phase | 12.98068576 s | 14.32712324 s | +1.34643748 s |
| post fields/output | 36.15433245 s | 36.0820299 s | -0.07230255 s |
| DLB rebalances | 8 | 8 | 0 |

Correctness checks from the log:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- Final particles: `2463818`, stuck `0`
- Final collisions: `35301`, candidates `63004`, acceptance `0.5602977589`
- Final total energy: `1.239467981`
- No Fatal/NaN match in the log.

Decision: reject and revert. Single-constraint partitioning made both move and collision slower. The final particle distribution was also worse (`max/min 3.770`) than the retained dual-constraint forced8 run (`max/min 2.953`). The next lower-risk objective test is to keep dual constraints but increase the move-weight exponent from `alpha=0.8` toward linear particle weighting.

## Candidate D: dual-constraint alpha=1.0 with retained 8 forced steps

Configuration delta:

```text
replicatedMeshDLBDualConstraint true;
replicatedMeshDLBInitialAlpha 1.0;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Run:

```text
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_dual_alpha1_forced8_500step_20260606 2>&1
```

Result:

| metric | forced 8-step retained | alpha=1.0 candidate | delta |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 135.52 s | +3.23 s |
| move only | 69.49945105 s | 76.01587475 s | +6.51642370 s |
| full evolve wall | 123.0267073 s | 126.4185434 s | +3.3918361 s |
| buildCellOccupancy | 6.216832977 s | 6.889014191 s | +0.672181214 s |
| collision phase | 12.98068576 s | 14.16367361 s | +1.18298785 s |
| post fields/output | 36.15433245 s | 37.89520374 s | +1.74087129 s |
| DLB rebalances | 8 | 8 | 0 |

Correctness checks from the log:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- Final particles: `2463588`, stuck `0`
- Final collisions: `34317`, candidates `59877`, acceptance `0.5731249061`
- Final total energy: `1.238951262`
- No Fatal/NaN match in the log.

Decision: reject and revert. Raising `alpha` toward linear particle weighting worsened move, collision, post/output, and end-to-end time. The best retained configuration remains dual constraint with `alpha=0.8` and forced steps `(120 170 220 270 320 370 420 470)`.

## Candidate E: remove the final 470-step forced DLB

Configuration delta:

```text
replicatedMeshDLBDualConstraint true;
replicatedMeshDLBInitialAlpha 0.8;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420);
```

Run:

```text
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_forced8_no470_500step_20260606 2>&1
```

Result:

| metric | previous best forced8 | no470 candidate | delta |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 130.89 s | -1.40 s |
| move only | 69.49945105 s | 68.57542261 s | -0.92402844 s |
| full evolve wall | 123.0267073 s | 121.6518645 s | -1.3748428 s |
| buildCellOccupancy | 6.216832977 s | 6.971536085 s | +0.754703108 s |
| collision phase | 12.98068576 s | 12.78095883 s | -0.19972693 s |
| post fields/output | 36.15433245 s | 37.82747705 s | +1.67314460 s |
| DLB rebalances | 8 | 7 | -1 |

Correctness checks from the log:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- Final particles: `2463743`, stuck `0`
- Final collisions: `35054`, candidates `64255`, acceptance `0.5455450938`
- Final total energy: `1.239748536`
- No Fatal/NaN match in the log.

Decision: keep. Removing the step-470 forced DLB reduced both move and end-to-end time. The final 470-step rebalance was too late in the run to pay back its migration/repartition cost. Current best configuration is:

```text
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420);
```

Next local test: remove the 420-step forced DLB as well and verify whether the 420 rebalance is still useful for the final 80 steps.

## Candidate F: remove both 420-step and 470-step forced DLB

Configuration delta:

```text
replicatedMeshDLBForceSteps (120 170 220 270 320 370);
```

Run:

```text
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_forced8_no420_no470_500step_20260606 2>&1
```

Result:

| metric | current best no470 | no420/no470 candidate | delta |
| --- | ---: | ---: | ---: |
| real | 130.89 s | 137.01 s | +6.12 s |
| move only | 68.57542261 s | 77.04926169 s | +8.47383908 s |
| full evolve wall | 121.6518645 s | 127.8606428 s | +6.2087783 s |
| buildCellOccupancy | 6.971536085 s | 7.415596836 s | +0.444060751 s |
| collision phase | 12.78095883 s | 14.60510432 s | +1.82414549 s |
| post fields/output | 37.82747705 s | 38.17155245 s | +0.34407540 s |
| DLB rebalances | 7 | 6 | -1 |

Correctness checks from the log:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- Final particles: `2463924`, stuck `0`
- Final collisions: `34127`, candidates `63199`, acceptance `0.5399927214`
- Final total energy: `1.239477642`
- No Fatal/NaN match in the log.

Decision: reject and revert. Removing the 420-step DLB caused the post-370 segment to become the critical path. The 420 rebalance is useful; only the final 470 rebalance should be removed.

## No470 repeat verification

Because Candidate E only improved the previous best by `1.40 s`, an immediate repeat was run with the same no470 configuration:

```text
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_forced8_no470_confirm_500step_20260606 2>&1
```

Repeat result:

| metric | no470 first run | no470 repeat | previous best forced8 |
| --- | ---: | ---: | ---: |
| real | 130.89 s | 137.65 s | 132.29 s |
| move only | 68.57542261 s | 74.5510391 s | 69.49945105 s |
| full evolve wall | 121.6518645 s | 128.369011 s | 123.0267073 s |
| DLB rebalances | 7 | 7 | 8 |

Correctness checks from the repeat:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- Final particles: `2463645`, stuck `0`
- Final collisions: `34384`, candidates `62741`, acceptance `0.5480307933`
- Final total energy: `1.238969612`
- No Fatal/NaN match in the log.

Decision update: no470 is not retained. The repeat was worse than the 8-step forced configuration and worse than the current MPI-only baseline in real wall time. Treat the first no470 result as non-repeatable run-to-run variance for now.

## Current retained configuration after forced-step sweep

```text
useOpenMP false;
openmpThreads 1;
replicatedMeshDLBDualConstraint true;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
replicatedMeshDLBInitialAlpha 0.8;
```

Retained 500-step result:

- `real 132.29 s`
- `move only 69.49945105 s`
- `full evolve wall 123.0267073 s`
- `DLB checks 500`
- `DLB rebalances 8`

Remaining gap to `ourmeshbkp`: `+14.98 s` real wall and `+22.33418745 s` move-only. The force-step sweep improved the current source baseline but did not close the inherited tracking/move gap.
