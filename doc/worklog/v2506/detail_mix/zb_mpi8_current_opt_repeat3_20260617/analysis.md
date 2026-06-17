# zb mpi8 current optimum repeat3, 2026-06-17

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Current optimum configuration tested:

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
```

The case `controlDict` was restored after the repeat set.

## Results

| run | exit | real [s] | execution@300 [s] | full evolve [s] | buildOcc [s] | migration max [s] | sizeX max [s] | wait max [s] | candidate gather max [s] | updateParticleCounts max [s] | DLB rebalances | DLB wall max [s] | DLB check max [s] | particles max/min | final energy | stuck | trigger/skip |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| current_opt_rep1 | 0 | 66.18 | 65.43 | 56.93766103 | 4.078883636 | 3.39763305 | 2.595145288 | 0.077398671 | 0.000003379 | 0.231650022 | 5 | 5.124735585 | 4.451945135 | 1.242859515 | 0.001953692983 | 0 | 50,100,150,200,250; step 300 skipped |
| current_opt_rep2 | 0 | 65.80 | 64.91 | 56.53195911 | 4.665883792 | 3.058505831 | 2.234161471 | 0.070885106 | 0.000003488 | 0.232906396 | 5 | 4.673208182 | 4.069275326 | 1.204576176 | 0.001954946958 | 0 | 50,100,150,200,250; step 300 skipped |
| current_opt_rep3 | 0 | 65.66 | 64.63 | 54.86951976 | 4.600188074 | 3.549816214 | 2.744947236 | 0.071945549 | 0.000002926 | 0.185637985 | 5 | 5.357699357 | 4.690161518 | 1.238312392 | 0.001957968272 | 0 | 50,100,150,200,250; step 300 skipped |

## Statistics

| metric | mean | pstdev | min | max |
|---|---:|---:|---:|---:|
| real [s] | 65.88 | 0.2196967607 | 65.66 | 66.18 |
| execution@300 [s] | 64.99 | 0.3314614105 | 64.63 | 65.43 |
| full evolve [s] | 56.1130466333 | 0.8947691790 | 54.86951976 | 56.93766103 |
| buildCellOccupancy [s] | 4.4483185007 | 0.2626030861 | 4.078883636 | 4.665883792 |
| migration max [s] | 3.3353183650 | 0.2053595434 | 3.058505831 | 3.549816214 |
| migration sizeX max [s] | 2.5247513317 | 0.2143859545 | 2.234161471 | 2.744947236 |
| migration wait max [s] | 0.0734097753 | 0.0028536061 | 0.070885106 | 0.077398671 |
| candidate gather max [s] | 0.0000032643 | 0.0000002433 | 0.000002926 | 0.000003488 |
| updateParticleCounts max [s] | 0.2167314677 | 0.0219923944 | 0.185637985 | 0.232906396 |
| DLB wall max [s] | 5.0518810413 | 0.2841512280 | 4.673208182 | 5.357699357 |
| DLB check max [s] | 4.4037939930 | 0.2557522437 | 4.069275326 | 4.690161518 |

## Comparison

Previous confirmed references:

- historical best single run, `single_alpha1_allgather_sar50`: `real=67.65 s`, `full evolve=61.78422282 s`
- old `alltoall_guard50` repeat pair before candidate-gather removal: `real=70.86/70.93 s`
- same config after candidate-gather removal, noisy first checks: `76.51 s` then `68.72 s`

This repeat set is the first stable confirmation that the current configuration beats the previous historical best:

```text
real mean = 65.88 s
best repeat = 65.66 s
vs historical best single = 67.65 s
```

Relative to the historical best single run, mean external improvement is about:

```text
(67.65 - 65.88) / 67.65 = 2.62%
```

## Interpretation

1. The current optimum is now confirmed by three stable repeats.  The `real`
   spread is only `0.52 s` from min to max, and all runs exit `0`.

2. `replicatedMeshGatherCandidates false` is mechanically effective.  Candidate
   gather max is reduced to about `3e-6 s`, compared with old alltoall_guard50
   values around `0.19-0.23 s`.

3. The end-of-run DLB guard is active.  All three runs detect a step-300 trigger
   and skip it via `replicatedMeshDLBMinRemainingSteps 50`, leaving five actual
   rebalances at steps `50,100,150,200,250`.

4. Regular Alltoall-size-exchange remains better than the previous no-Alltoall
   path on this 8-rank case.  Wait time is very small (`~0.07 s` max), while the
   remaining migration cost is mostly explicit size exchange.

5. Correctness gates are clean for this comparison: all runs have `stuck=0`,
   similar final energy, and normal completion.

## Current Recommendation

Use this as the current `zb` MPI8 replicated-mesh DLB best configuration:

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
```

## Artifacts

- Logs: `logs/current_opt_rep{1,2,3}.log`
- Timings: `logs/current_opt_rep{1,2,3}.time`
- Control snapshots: `controlDicts/current_opt_rep{1,2,3}_controlDict`
