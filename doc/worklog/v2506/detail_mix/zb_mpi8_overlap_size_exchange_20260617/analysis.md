# zb mpi8 overlap size-exchange attempt, 2026-06-17

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Starting point was the current confirmed best `single_alpha1` configuration from
`zb_mpi8_current_opt_repeat3_20260617`, with one additional option:

```text
replicatedMeshOverlapSizeExchange true;
```

The code path is opt-in and defaults to `false`.  It starts `MPI_Ialltoall` for
the flat-transfer size exchange immediately after packing, then overlaps it with
local outgoing-parcel deletion.  The no-Alltoall path is unchanged.

The live case `controlDict` was restored to its pre-test state after the repeat
set.

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
replicatedMeshOverlapSizeExchange true;
```

## Results

| run | exit | real [s] | execution@300 [s] | full evolve [s] | move only [s] | collision [s] | migration max [s] | sizeX max [s] | wait max [s] | candidate gather max [s] | DLB rebalances | DLB wall max [s] | DLB check max [s] | particles max/min | final energy | stuck | trigger/skip |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| overlap_size_exchange_rep1 | 0 | 68.50 | 67.58 | 56.03496006 | 39.08004813 | 10.3044029 | 3.065039225 | 2.208704504 | 0.133768223 | 0.000003689 | 5 | 5.083444019 | 4.440773875 | 1.191506352 | 0.001952365074 | 0 | 50,100,150,200,250; step 300 skipped |
| overlap_size_exchange_rep2 | 0 | 70.37 | 69.56 | 60.09292996 | 47.2291355 | 13.23147504 | 3.569241599 | 2.731578073 | 0.112720102 | 0.000003958 | 5 | 5.896441072 | 5.288353173 | 1.176259794 | 0.001957156926 | 0 | 50,100,150,200,250; step 300 skipped |
| overlap_size_exchange_rep3 | 0 | 66.61 | 65.86 | 57.43833693 | 44.79432422 | 10.52786289 | 3.79059604 | 2.905803859 | 0.177340955 | 0.000002833 | 5 | 5.396473711 | 4.748388507 | 1.198430812 | 0.001955707176 | 0 | 50,100,150,200,250; step 300 skipped |

## Statistics

| metric | overlap mean | overlap pstdev | current-opt mean | delta |
|---|---:|---:|---:|---:|
| real [s] | 68.4933333333 | 1.5350208106 | 65.88 | +2.6133333333 |
| execution@300 [s] | 67.6666666667 | 1.5117612980 | 64.99 | +2.6766666667 |
| full evolve [s] | 57.8554089833 | 1.6827045196 | 56.1130466333 | +1.7423623500 |
| migration max [s] | 3.4749589547 | 0.3036171775 | 3.3353183650 | +0.1396405897 |
| migration sizeX max [s] | 2.6153621453 | 0.2962166820 | 2.5247513317 | +0.0906108136 |
| migration wait max [s] | 0.1412764267 | 0.0269102645 | 0.0734097753 | +0.0678666514 |
| DLB wall max [s] | 5.4587862673 | 0.3348165598 | 5.0518810413 | +0.4069052260 |
| DLB check max [s] | 4.8258385183 | 0.3503298872 | 4.4037939930 | +0.4220445253 |

## Interpretation

This attempt does not beat the current best configuration.

1. End-to-end time regresses.  The overlap repeat mean is `68.49 s`, compared
   with the current best repeat mean `65.88 s`.  Even the best overlap repeat
   (`66.61 s`) is slower than the best current-opt repeat (`65.66 s`).

2. The intended overlap is not reliable.  Rep1 reduced size-exchange max, but
   rep2/rep3 did not.  Across repeat3, `sizeX max` is slightly worse than the
   current optimum, and `wait max` is almost doubled.

3. The theoretical overlap window is small.  The overlapped local-prep phase is
   only about `0.22 s` max, while size exchange is about `2.5 s`.  Even perfect
   progress could only hide a small fraction of the collective, and typical MPI
   nonblocking collectives do not guarantee useful progress while the process is
   busy deleting parcels outside MPI.

4. The current critical path remains move/collision imbalance and DLB check
   variability, not only migration size exchange.  Rep2 is a clear example:
   migration is not catastrophic, but move-only, collision, rank wall imbalance,
   and DLB check time all worsen.

## Recommendation

Do not enable `replicatedMeshOverlapSizeExchange` in the current best
configuration.

Keep the code as a default-off experimental switch if useful for future MPI
progress experiments, but the current recommended `zb` MPI8 configuration
remains:

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

## Artifacts

- Logs: `logs/overlap_size_exchange_rep{1,2,3}.log`
- Timings: `logs/overlap_size_exchange_rep{1,2,3}.time`
- Exits: `logs/overlap_size_exchange_rep{1,2,3}.exit`
- Test control: `controlDicts/overlap_size_exchange_controlDict`
- Pre-test live control: `controlDicts/zb_mpi8_controlDict_before_overlap`
- Restored live control: `controlDicts/zb_mpi8_controlDict_restored_after_overlap`
