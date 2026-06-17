# zb mpi8 single_alpha1 check/cadence follow-up, 2026-06-13

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Baseline for comparison:

`doc/worklog/v2506/detail_mix/zb_mpi8_single_alpha1_repeat3_20260613/analysis.md`

Common controls:

- `mpirun -np 8 dsmcFoam+`
- 300 steps
- `replicatedMeshAutoDLB true`
- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBAlpha 1`
- `replicatedMeshDLBForceSteps ()`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha false`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBVsizeExp 0`

Tested variants:

| Variant | check collective | SAR steps |
|---|---|---:|
| repeat3 baseline mean | allgather | 10 |
| single_alpha1_allreduce_sar10 | allreduce | 10 |
| single_alpha1_allgather_sar50 | allgather | 50 |
| single_alpha1_allreduce_sar50 | allreduce | 50 |

Note: `Phase C auto DLB checks` increments every step that enters
`autoRebalance()`.  It is not the number of expensive global load collectives.
In `legacyWindow` mode, `replicatedMeshSARSteps` controls how often the global
load check is evaluated.

Source check:

- `dsmcCloud.C` calls `replicatedMesh_->autoRebalance()` once per evolve step
  when replicated mesh is active.
- `dsmcReplicatedMesh.C::autoRebalance()` increments
  `autoRebalanceChecks_` before any cadence gate, so the summary still reports
  `checks=300` for all 300-step runs.
- In `legacyWindow`, the expensive load-extrema collective is under
  `if (ndecps_ % sarSteps_ == 0)`.  Therefore `replicatedMeshSARSteps 50`
  reduces the global timing checks by about 5x relative to the default
  `SARSteps=10`.
- `replicatedMeshDLBCheckCollective allreduce` currently performs two scalar
  collectives (`MPI_MAX` and `MPI_MIN`).  `allgather` performs one scalar
  allgather and computes min/max locally.  At 8 ranks this is not automatically
  cheaper than allgather.

## Performance

| Variant | real [s] | execution@300 [s] | full evolve [s] | move+collide [s] | move [s] | buildOcc [s] | collision [s] | post [s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| repeat3 baseline mean | 68.31333333 | 67.62333333 | 67.56856274 | 67.29541664 | 45.34930456 | 4.824498188 | 10.87385769 | 0.277891099 |
| allreduce_sar10 | 69.56 | 68.03 | 68.06842613 | 67.79557865 | 45.84640637 | 4.612909571 | 14.84666563 | 0.277374337 |
| allgather_sar50 | 67.65 | 65.90 | 61.78422282 | 61.59120205 | 49.65921506 | 4.617690398 | 14.19315533 | 0.194665781 |
| allreduce_sar50 | 69.04 | 67.42 | 63.32075528 | 63.13013652 | 47.87971708 | 6.100416443 | 9.322817921 | 0.197944849 |

## DLB and Balance

| Variant | checks | rebalances | trigger steps | DLB wall max [s] | check max [s] | ParMETIS max [s] | DLB migration max [s] | migration wall max [s] | particles max/min | rank wall max/min |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| repeat3 baseline mean | 300 | 5.333333333 | mixed | 18.94386891 | 18.07361404 | 0.186567323 | 0.73262435 | 3.66107507 | 1.099823693 | 1.312828772 |
| allreduce_sar10 | 300 | 6 | 50,100,150,200,250,300 | 19.648361 | 18.7339347 | 0.2113502 | 0.765731706 | 3.505031427 | 1.094527363 | 1.239125535 |
| allgather_sar50 | 300 | 6 | 50,100,150,200,250,300 | 5.34415917 | 4.428155121 | 0.23355827 | 0.756277535 | 3.807108542 | 1.120322481 | 1.316795741 |
| allreduce_sar50 | 300 | 6 | 50,100,150,200,250,300 | 5.276535654 | 4.420695816 | 0.208027368 | 0.709637563 | 3.949020055 | 1.138100078 | 1.26273096 |

## Correctness Gate

| Variant | exit | iterations | final particles | final energy | cumulative global collisions | cumulative global candidates | stuck |
|---|---:|---:|---:|---:|---:|---:|---:|
| allreduce_sar10 | 0 | 300 | 1960729 | 0.001954837635 | 39749954 | 72518939 | 0 |
| allgather_sar50 | 0 | 300 | 1958911 | 0.001955262483 | 40960793 | 74884607 | 0 |
| allreduce_sar50 | 0 | 300 | 1958919 | 0.001954723025 | 39862128 | 72210192 | 0 |

## Interpretation

1. Changing `allgather -> allreduce` alone does not solve the problem.
   `allreduce_sar10` has similar or slightly worse check cost than the
   repeat3 baseline: check max `18.73 s` vs baseline mean `18.07 s`.

2. Increasing `replicatedMeshSARSteps` from 10 to 50 is the effective change.
   Both `sar50` variants reduce check max to about `4.42 s`, and DLB wall max
   to about `5.3 s`.

3. `allgather_sar50` is the best external real time in this single-run
   follow-up: `67.65 s`, slightly faster than repeat3 mean `68.31 s`.  It also
   has the lowest full-evolve wall in this group: `61.78 s`.

4. `allreduce_sar50` gives a similar check-cost reduction but worse external
   real time (`69.04 s`) than `allgather_sar50`.  There is no evidence here
   that allreduce is better at 8 ranks.

5. All `sar50` runs still rebalance at step 300.  The next DLB test should keep
   `SARSteps=50` and remove the terminal step-300 rebalance, either via a
   forced no-final schedule for measurement or a real min-remaining-steps gate
   in code.

6. The wall-time gain is smaller than the DLB-check reduction because this
   300-step `zb` run still has large run-to-run movement/collision variation.
   The clean signal is the DLB component itself: check max drops from
   `18.07 s` mean to `4.42 s`, while ParMETIS remains only about `0.2 s`.
   This confirms the target is cadence/check overhead, not partitioner compute.

7. The best variant in this group is the single-item cadence change
   `allgather_sar50`.  The combination `allreduce_sar50` keeps the same low
   check cost, but its external `real` time is worse and it has higher
   `buildCellOccupancy`/migration spread in this run.  So there is no reason to
   switch the check collective away from `allgather` for MPI8.

## Current Recommendation

For `zb` MPI8 replicated mesh, use this as the next candidate control state:

- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBAlpha 1`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha false`
- `replicatedMeshDLBCheckCollective allgather`
- `replicatedMeshSARSteps 50`

Then test no-final DLB on top of it.

## Artifacts

- Logs: `logs/`
- Control snapshots: `controlDicts/`
- Restored case controlDict: `controlDicts/zb_mpi8_controlDict_restored`
