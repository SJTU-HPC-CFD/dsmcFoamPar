# zb mpi8 current optimum retest, 2026-06-17

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Purpose: retest the current best configuration after the default-off
`replicatedMeshOverlapSizeExchange` experiment, to make sure the recommended
configuration still performs correctly.

The live case `controlDict` before this retest was not the best configuration
(`replicatedMeshNoAlltoall true`, `replicatedMeshDLBAlpha 0.8`).  Therefore the
retest used the saved best control from
`zb_mpi8_current_opt_repeat3_20260617/controlDicts/current_opt_rep1_controlDict`
and explicitly added:

```text
replicatedMeshOverlapSizeExchange false;
```

The live case `controlDict` was restored after each run and after the full
repeat set.

## Configuration

```text
replicatedMeshNoAlltoall false;
replicatedMeshDLBDualConstraint false;
replicatedMeshDLBAlpha 1;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshSARSteps 50;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBParticleGate false;
replicatedMeshDLBAdaptiveAlpha false;
replicatedMeshDLBMinRemainingSteps 50;
replicatedMeshGatherCandidates false;
replicatedMeshOverlapSizeExchange false;
```

Runtime confirmed:

```text
Replicated mesh: gatherCandidates=0
Replicated mesh: overlapSizeExchange=0
Phase C auto DLB ... sarSteps=50 ... checkCollective=allgather ... adaptiveAlpha=0
```

## Results

| run | exit | real [s] | execution@300 [s] | full evolve [s] | buildOcc [s] | migration max [s] | sizeX max [s] | wait max [s] | candidate gather max [s] | updateParticleCounts max [s] | DLB rebalances | DLB wall max [s] | DLB check max [s] | particles max/min | final energy | stuck | trigger/skip |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| current_opt_retest_rep1 | 0 | 67.94 | 67.20 | 56.79943256 | 4.778274201 | 3.499779839 | 2.713016930 | 0.067694473 | 0.000003214 | 0.298138501 | 5 | 6.009258826 | 5.360506197 | 1.236263362 | 0.001953006938 | 0 | 50,100,150,200,250; step 300 skipped |
| current_opt_retest_rep2 | 0 | 67.25 | 66.50 | 53.90763881 | 3.894422764 | 3.401282620 | 2.575646076 | 0.060293785 | 0.000003225 | 0.234906634 | 5 | 5.315277191 | 4.660296880 | 1.153013592 | 0.001955106044 | 0 | 50,100,150,200,250; step 300 skipped |
| current_opt_retest_rep3 | 0 | 64.42 | 63.69 | 55.20683303 | 3.883461674 | 3.294992279 | 2.358409811 | 0.075156681 | 0.000004990 | 0.202147937 | 5 | 4.909811394 | 4.235731028 | 1.199406975 | 0.001954317644 | 0 | 50,100,150,200,250; step 300 skipped |

## Statistics

| metric | retest mean | retest pstdev | previous current-opt mean | delta |
|---|---:|---:|---:|---:|
| real [s] | 66.5366666667 | 1.5229868315 | 65.88 | +0.6566666667 |
| execution@300 [s] | 65.7966666667 | 1.5168021917 | 64.99 | +0.8066666667 |
| full evolve [s] | 55.3046348000 | 1.1825936645 | 56.1130466333 | -0.8084118333 |
| buildCellOccupancy [s] | 4.1853862130 | 0.4192589980 | 4.4483185007 | -0.2629322877 |
| migration max [s] | 3.3986849127 | 0.0836243475 | 3.3353183650 | +0.0633665477 |
| migration sizeX max [s] | 2.5490242723 | 0.1459865120 | 2.5247513317 | +0.0242729406 |
| migration wait max [s] | 0.0677149797 | 0.0060677692 | 0.0734097753 | -0.0056947956 |
| candidate gather max [s] | 0.0000038097 | 0.0000008346 | 0.0000032643 | +0.0000005454 |
| updateParticleCounts max [s] | 0.2450643573 | 0.0398407798 | 0.2167314677 | +0.0283328896 |
| DLB wall max [s] | 5.4114491370 | 0.4539698568 | 5.0518810413 | +0.3595680957 |
| DLB check max [s] | 4.7521780350 | 0.4637610054 | 4.4037939930 | +0.3483840420 |

## Interpretation

1. The current best configuration is still valid.  All three runs exit `0`, all
   have `stuck=0`, runtime confirms `gatherCandidates=0` and
   `overlapSizeExchange=0`, and all runs perform five real DLB rebalances at
   steps `50,100,150,200,250` with the step-300 trigger skipped by the near-end
   guard.

2. End-to-end performance is not as tight as the previous repeat3.  The retest
   `real` mean is `66.54 s`, compared with previous `65.88 s`.  The spread is
   larger (`64.42-67.94 s`), so this retest shows more machine/runtime noise
   than the earlier stable set.

3. There is no evidence that the default-off overlap experiment damaged the
   recommended path.  The migration submetrics are essentially the same:
   `sizeX max` mean `2.549 s` versus previous `2.525 s`, and `wait max` is
   slightly lower.  Candidate gather remains near zero.

4. The main difference is run-to-run DLB/check variability and stochastic
   collision/move balance.  Rep1 has `DLB check max=5.36 s`, while rep3 has
   `4.24 s` and reaches the best external time, `real=64.42 s`.

## Current Recommendation

Keep the same current best configuration:

```text
replicatedMeshNoAlltoall false;
replicatedMeshDLBDualConstraint false;
replicatedMeshDLBAlpha 1;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshSARSteps 50;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBParticleGate false;
replicatedMeshDLBAdaptiveAlpha false;
replicatedMeshDLBMinRemainingSteps 50;
replicatedMeshGatherCandidates false;
replicatedMeshOverlapSizeExchange false;
```

Use `real ~= 66.5 s` as the latest retest mean under the current machine state,
with the caveat that this set is noisier than the previous `65.88 s` repeat.

## Artifacts

- Logs: `logs/current_opt_retest_rep{1,2,3}.log`
- Timings: `logs/current_opt_retest_rep{1,2,3}.time`
- Exits: `logs/current_opt_retest_rep{1,2,3}.exit`
- Test control: `controlDicts/current_opt_retest_controlDict`
- Pre-test live control: `controlDicts/zb_mpi8_controlDict_before_retest`
- Restored live control: `controlDicts/zb_mpi8_controlDict_restored_after_retest`
