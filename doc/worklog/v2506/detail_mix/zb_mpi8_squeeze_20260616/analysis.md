# zb mpi8 squeeze attempts, 2026-06-16

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Baseline candidate:

`single_alpha1_allgather_sar50`

Common candidate controls unless noted:

- `mpirun -np 8 dsmcFoam+`
- 300 steps
- `replicatedMeshAutoDLB true`
- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBAlpha 1`
- `replicatedMeshDLBCheckCollective allgather`
- `replicatedMeshSARSteps 50`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha false`
- `replicatedMeshDLBVsizeExp 0`

The production case `controlDict` was restored after each run.

## Source Changes Tested

### Kept: near-end DLB guard

Added opt-in control:

```text
replicatedMeshDLBMinRemainingSteps 50;
```

Default is `0`, so existing control files keep previous behavior.

When a DLB trigger is detected, `autoRebalance()` estimates the planned total
step count from `startTime`, `endTime`, and `deltaT`.  If the remaining steps
are less than `replicatedMeshDLBMinRemainingSteps`, it logs a skip and avoids
ParMETIS/migration for that trigger.

### Rejected: skip non-output `runTime.write()`

An opt-in `replicatedMeshSkipNonOutputWrite` was tested and then removed.
The result was effectively unchanged in `real` time (`73.03 s` vs `73.09 s`)
and the internal partition/load trajectory was worse, so the code was not kept.

## Results

| Variant | changed controls | exit | real [s] | execution@300 [s] | full evolve [s] | move [s] | buildOcc [s] | collision [s] |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| previous best-candidate single run | historical `allgather_sar50` | 0 | 67.65 | 65.90 | 61.78422282 | 49.65921506 | 4.617690398 | 14.19315533 |
| current-best retest | same candidate, before squeeze | 0 | 77.18 | 75.25 | 57.96273238 | 44.46065634 | 4.366002621 | 10.98005539 |
| best_terminal300_skipdiag | `nTerminalOutputs 300`, `DLBSkipPostDiag true` | 0 | 74.24 | 73.29 | 69.79375438 | 45.02454614 | 4.475151299 | 10.60595723 |
| best_terminal300_skipdiag_profileoff | above + `profileSummary false` | 0 | 79.61 | 79.04 | - | - | - | - |
| forced5_terminal300_skipdiag | forced DLB `(50 100 150 200 250)`, terminal300, skipPostDiag | 0 | 74.48 | 73.40 | 70.48311165 | 46.97016698 | 5.616683353 | 13.22883834 |
| current_best_minremaining50 | `DLBMinRemainingSteps 50` | 0 | 73.09 | 72.01 | 64.92950539 | 52.34311918 | 5.27835767 | 11.52732143 |
| current_best_minremaining50_skipwrite | minRemaining + rejected skip-write switch | 0 | 73.03 | 71.94 | 66.27682377 | 54.31273268 | 5.456147976 | 10.68482634 |

## DLB / Balance

| Variant | rebalances | trigger/skip behavior | DLB wall max [s] | check max [s] | ParMETIS max [s] | particles max/min | rank wall max/min |
|---|---:|---|---:|---:|---:|---:|---:|
| previous best-candidate single run | 6 | 50,100,150,200,250,300 | 5.34415917 | 4.428155121 | 0.23355827 | 1.120322481 | 1.316795741 |
| current-best retest | 5 | 50,100,150,200,250 | 4.941800583 | 4.187642619 | 0.192254026 | 1.136633076 | 1.21532185 |
| best_terminal300_skipdiag | 6 | 50,100,150,200,250,300 | 6.479168797 | 5.514981517 | 0.224438033 | 1.090984681 | 1.056910569 |
| forced5_terminal300_skipdiag | 5 | forced 50,100,150,200,250 | 6.226060479 | 5.452481337 | 0.187308249 | 1.28007781 | 1.061028352 |
| current_best_minremaining50 | 5 | 50,100,150,200,250; step 300 skipped | 5.984115494 | 5.253919627 | 0.182551748 | 1.25815154 | 1.394597091 |
| current_best_minremaining50_skipwrite | 4 | 50,100,150,200; step 300 skipped | 6.014259126 | 5.452252599 | 0.151419402 | 1.351496778 | 1.48709369 |

## Correctness Gate

| Variant | iterations | final particles | final energy | cumulative collisions | stuck |
|---|---:|---:|---:|---:|---:|
| best_terminal300_skipdiag | 300 | 1959719 | 0.001951808166 | 39781991 | 0 |
| best_terminal300_skipdiag_profileoff | 300 | 1958231 | 0.001954441389 | not profiled | 0 |
| forced5_terminal300_skipdiag | 300 | 1959517 | 0.001956177384 | 40662847 | 0 |
| current_best_minremaining50 | 300 | 1958789 | 0.001949936782 | 40393919 | 0 |
| current_best_minremaining50_skipwrite | 300 | 1960063 | 0.00195392772 | 40409736 | 0 |

## Interpretation

1. No tested squeeze variant beats the historical external best
   `allgather_sar50` sample (`67.65 s real`).  The best new external result is
   `73.03 s`, but that used the rejected skip-write switch and had worse
   internal balance.

2. Against the noisy 2026-06-16 current-best retest (`77.18 s real`), the
   near-end guard improved external time to `73.09 s` and correctly skipped
   the step-300 DLB:

   ```text
   Phase C auto DLB skipped near end at step 300:
   remainingSteps=0 minRemainingSteps=50
   ```

   However, its full-evolve profile was worse (`64.93 s` vs `57.96 s`), so the
   gain is not a clean solver-core improvement.

3. `nTerminalOutputs 300` is not a free win.  It reduced terminal output volume,
   but in this run it changed the DLB trajectory enough to produce a step-300
   DLB and a large migration wall max (`22.03 s`).  Do not use it as a default
   performance setting without more repeats.

4. `profileSummary false` should not be used for performance claims.  It removed
   profile visibility and produced a slower external time (`79.61 s`).

5. Forced DLB at `50,100,150,200,250` did not improve end-to-end time.  It
   avoided the final DLB but worsened move/build/collision and particle balance,
   so auto DLB remains better.

6. `replicatedMeshDLBMinRemainingSteps 50` is still worth keeping as a guardrail:
   it is opt-in, default-off, and solves the specific terminal-DLB pathology
   without changing earlier triggers.  It should be enabled when benchmarking
   cases where the last DLB can land at the final step.

## Current Recommendation

Keep the previously selected candidate as the main configuration:

- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBAlpha 1`
- `replicatedMeshDLBCheckCollective allgather`
- `replicatedMeshSARSteps 50`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha false`

Optionally add for production/benchmark safety:

- `replicatedMeshDLBMinRemainingSteps 50`

Do not adopt these attempted squeeze controls as defaults:

- `nTerminalOutputs 300`
- `profileSummary false`
- forced DLB `(50 100 150 200 250)`
- `replicatedMeshSkipNonOutputWrite`

## Artifacts

- Logs: `logs/`
- Control snapshots: `controlDicts/`
- Restored case controlDict: `controlDicts/zb_mpi8_controlDict_restored_after_skipwrite`
