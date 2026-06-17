# zb mpi8 current-best retest, 2026-06-16

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Configuration under test:

- `mpirun -np 8 dsmcFoam+`
- 300 steps: `endTime 1.99201e-05`
- `useOpenMP false`, `openmpThreads 1`
- `replicatedMeshAutoDLB true`
- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBAlpha 1`
- `replicatedMeshDLBForceSteps ()`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBCheckCollective allgather`
- `replicatedMeshSARSteps 50`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha false`
- `replicatedMeshDLBVsizeExp 0`

The first launch attempt used a too-strict shell wrapper (`set -u`) and failed
before solver startup while sourcing the environment.  It produced an empty
solver log and is not a performance sample.  The formal sample is
`current_best_retest_rep1`.

The production case `controlDict` was restored after the run.

## Result

| run | exit | real [s] | execution@300 [s] | full evolve [s] | move+collide [s] | move [s] | buildOcc [s] | collision [s] | post [s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| current_best_retest_rep1 | 0 | 77.18 | 75.25 | 57.96273238 | 57.75175882 | 44.46065634 | 4.366002621 | 10.98005539 | 0.228299266 |

## DLB and Balance

| run | checks | rebalances | trigger steps | DLB wall max [s] | check max [s] | ParMETIS max [s] | DLB migration max [s] | migration wall max [s] | particles max/min | rank wall max/min |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| current_best_retest_rep1 | 300 | 5 | 50,100,150,200,250 | 4.941800583 | 4.187642619 | 0.192254026 | 0.60617592 | 3.635295055 | 1.136633076 | 1.21532185 |

## Correctness Gate

| run | iterations | final particles | final energy | cumulative global collisions | cumulative global candidates | stuck |
|---|---:|---:|---:|---:|---:|---:|
| current_best_retest_rep1 | 300 | 1958143 | 0.001954403543 | 40666925 | 74472952 | 0 |

## Comparison With Previous Best-Candidate Single Run

Previous candidate:

`doc/worklog/v2506/detail_mix/zb_mpi8_single_alpha1_check_cadence_20260613/logs/single_alpha1_allgather_sar50.log`

| metric | 2026-06-13 allgather_sar50 | 2026-06-16 retest | delta |
|---|---:|---:|---:|
| real [s] | 67.65 | 77.18 | +9.53 |
| execution@300 [s] | 65.90 | 75.25 | +9.35 |
| full evolve [s] | 61.78422282 | 57.96273238 | -3.82149044 |
| move+collide [s] | 61.59120205 | 57.75175882 | -3.83944323 |
| move [s] | 49.65921506 | 44.46065634 | -5.19855872 |
| buildOcc [s] | 4.617690398 | 4.366002621 | -0.251687777 |
| collision [s] | 14.19315533 | 10.98005539 | -3.21309994 |
| DLB wall max [s] | 5.34415917 | 4.941800583 | -0.402358587 |
| DLB check max [s] | 4.428155121 | 4.187642619 | -0.240512502 |
| rebalances | 6 | 5 | -1 |
| rank wall max/min | 1.316795741 | 1.21532185 | -0.101473891 |

## Interpretation

1. The current-best DLB configuration is still valid.  The retest keeps the low
   check-cost behavior: `Phase C auto DLB check max = 4.19 s`, slightly better
   than the previous `4.43 s`.

2. This run avoided the previous step-300 terminal DLB.  Rebalances are
   `5` instead of `6`, with triggers only at `50,100,150,200,250`.  This is the
   desired direction for the no-final-DLB concern, although it happened
   naturally here because step 300 did not satisfy the trigger.

3. The solver-core profile improved versus the previous `allgather_sar50`
   sample: full evolve `61.78 -> 57.96 s`, collision `14.19 -> 10.98 s`,
   move `49.66 -> 44.46 s`.

4. The external end-to-end time did not reproduce the earlier best.  `real`
   worsened from `67.65` to `77.18 s`, and `ExecutionTime@300` from `65.90` to
   `75.25 s`.  This extra time is outside the reported full-evolve profile.
   It is already visible by iteration 10 (`ExecutionTime = 12.17 s` here versus
   `5.18 s` in the previous sample), so it should not be attributed to DLB or
   collision kernel performance.

5. For optimization decisions, this retest strengthens the internal DLB
   conclusion but does not replace the external-time best.  The current best
   candidate remains:

   - `DualConstraint false`
   - `DLBAlpha 1`
   - `DLBCheckCollective allgather`
   - `SARSteps 50`
   - `DLBParticleGate false`
   - `DLBAdaptiveAlpha false`

   A clean repeated run should be used before claiming a new end-to-end best,
   because this sample has large non-profile startup/setup overhead.

## Artifacts

- Formal log: `logs/current_best_retest_rep1.log`
- Formal time: `logs/current_best_retest_rep1.time`
- Formal exit: `logs/current_best_retest_rep1.exit`
- Formal control snapshot: `controlDicts/current_best_retest_rep1_controlDict`
- Restored control snapshot: `controlDicts/zb_mpi8_controlDict_restored_after_rep1`
