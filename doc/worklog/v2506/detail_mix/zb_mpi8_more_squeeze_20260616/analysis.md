# zb mpi8 additional squeeze attempts, 2026-06-16

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Baseline control snapshot:

`doc/worklog/v2506/detail_mix/zb_mpi8_single_alpha1_check_cadence_20260613/controlDicts/single_alpha1_allgather_sar50_controlDict`

Common baseline controls:

- `mpirun -np 8 dsmcFoam+`
- 300 steps
- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBAlpha 1`
- `replicatedMeshDLBCheckCollective allgather`
- `replicatedMeshSARSteps 50`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha false`

The production case `controlDict` was restored after the runs.

## Important Correction

The first four variants with `minremaining50` in the filename did not actually
contain `replicatedMeshDLBMinRemainingSteps 50`; the insertion command missed the
end-of-file marker.  Treat them as unguarded one-factor tests.

The guarded variants are:

- `alltoall_guard50`
- `alltoall_guard50_rep2`
- `alltoall_guard50_skipdiag`

## Results

| Variant | Main change | guard50 active | exit | real [s] | execution@300 [s] | full evolve [s] | buildOcc [s] | migration max [s] | sizeX max [s] | wait max [s] | DLB rebalances |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| historical best single run | `NoAlltoall true`, SAR50 | no | 0 | 67.65 | 65.90 | 61.78422282 | 4.617690398 | 3.807108542 | 0 | 2.595128201 | 6 |
| current-best retest | `NoAlltoall true`, SAR50 | no | 0 | 77.18 | 75.25 | 57.96273238 | 4.366002621 | 3.635295055 | 0 | 2.712906627 | 5 |
| previous guarded no-alltoall | `NoAlltoall true`, guard50 | yes | 0 | 73.09 | 72.01 | 64.92950539 | 5.27835767 | 4.104941978 | 0 | 3.037809792 | 5 |
| sar100_minremaining50 | SAR100, no actual guard | no | 0 | 75.14 | 74.94 | 69.80916433 | 7.588610289 | 3.660696632 | 0 | 2.603656365 | 2 |
| mingap100_minremaining50 | minGap100, no actual guard | no | 0 | 78.74 | 76.84 | 69.03662309 | 6.868383783 | 4.33565509 | 0 | 3.295877961 | 3 |
| threshold2_minremaining50 | threshold2.0, no actual guard | no | 0 | 76.95 | 76.85 | 69.29664993 | 5.877969695 | 3.9442013 | 0 | 2.866534521 | 3 |
| alltoall_minremaining50 | `NoAlltoall false`, no actual guard | no | 0 | 72.42 | 71.16 | 67.3848837 | 4.791773081 | 3.720807226 | 2.551709879 | 0.071689674 | 6 |
| alltoall_guard50 | `NoAlltoall false`, guard50 | yes | 0 | 70.86 | 69.54 | 62.49474379 | 5.231935936 | 3.720826666 | 2.79109105 | 0.070676603 | 5 |
| alltoall_guard50_rep2 | repeat of alltoall_guard50 | yes | 0 | 70.93 | 69.70 | 62.42188458 | 4.534247594 | 3.360279869 | 2.250693878 | 0.078226887 | 5 |
| alltoall_guard50_skipdiag | alltoall_guard50 + skipPostDiag | yes | 0 | 74.59 | 73.19 | 65.42705046 | 4.743621858 | 4.252296919 | 3.23402168 | 0.092672523 | 5 |

`alltoall_guard50` two-run mean:

- `real = 70.895 s`
- population stdev over the two runs: `0.035 s`

## DLB and Correctness

| Variant | trigger/skip behavior | DLB wall max [s] | DLB check max [s] | ParMETIS max [s] | DLB migration max [s] | particles max/min | final energy | stuck |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| sar100_minremaining50 | 200,300 | 3.395839158 | 3.023023314 | 0.097238926 | 0.324337869 | 1.119855608 | 0.001955938595 | 0 |
| mingap100_minremaining50 | 50,150,250 | 6.614978391 | 6.149613668 | 0.129703478 | 0.386046213 | 1.217863068 | 0.001955424815 | 0 |
| threshold2_minremaining50 | 100,200,300 | 5.957546529 | 5.41646368 | 0.109522403 | 0.463660813 | 1.108631902 | 0.001955247389 | 0 |
| alltoall_minremaining50 | 50,100,150,200,250,300 | 5.76216278 | 4.857318737 | 0.215896652 | 0.754396207 | 1.173866671 | 0.00195653673 | 0 |
| alltoall_guard50 | 50,100,150,200,250; step 300 skipped | 4.506035192 | 3.778185699 | 0.180620973 | 0.596189798 | 1.200627186 | 0.001959286903 | 0 |
| alltoall_guard50_rep2 | 50,100,150,200,250; step 300 skipped | 5.43824471 | 4.659641272 | 0.169417976 | 0.645470535 | 1.222443149 | 0.001953648664 | 0 |
| alltoall_guard50_skipdiag | 50,100,150,200,250; step 300 skipped | 6.642212407 | 5.849849832 | 0.172482465 | 0.64729265 | 1.195674865 | 0.001954974374 | 0 |

## Interpretation

1. Reducing DLB frequency further did not help.  `SARSteps=100`,
   `minGapSteps=100`, and `imbalanceThreshold=2.0` all reduced actual
   repartitions, but full-evolve time grew to roughly `69 s`.  The saved DLB
   overhead was smaller than the cost of running longer under poorer balance.

2. The useful new lever is `replicatedMeshNoAlltoall false`.  On this 8-rank
   case, the regular Alltoall-size-exchange path moves cost from `MPI_Waitall`
   into the explicit size exchange, but the total critical path is better and
   more stable after adding the near-end guard.

3. `alltoall_guard50` is a stable improvement over the same-day no-alltoall
   guarded run (`73.09 s -> 70.90 s mean real`), but it still does not beat the
   historical fastest single run (`67.65 s`).  The best guarded Alltoall
   full-evolve time (`62.42-62.49 s`) is close to the historical best
   full-evolve time (`61.78 s`).

4. `replicatedMeshDLBSkipPostDiag true` should not be adopted.  In the Alltoall
   guarded setting it worsened `real` to `74.59 s` and increased DLB check max.

5. `updateParticleCounts` is not a worthwhile target here.  Its reported max is
   about `0.02-0.04 s`, much smaller than DLB/migration and move/build/collision
   variation.

6. Do not increase `replicatedMeshMigrateInterval` as a time shortcut.  The
   source comment in `dsmcReplicatedMesh::initialize()` states that larger
   intervals can leave particles in non-owned cells and lose collisions; the
   current case already uses interval 10, so pushing it higher would trade
   correctness for speed.

## Current Best Candidate From This Round

Use this only as the new candidate to retest, not as a confirmed replacement for
the historical best:

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
```

Do not adopt from this round:

- `replicatedMeshSARSteps 100`
- `replicatedMeshDLBMinGapSteps 100`
- `replicatedMeshDLBImbalanceThreshold 2.0`
- `replicatedMeshDLBSkipPostDiag true`

## Artifacts

- Logs: `logs/`
- Control snapshots: `controlDicts/`
- Recommended snapshot: `controlDicts/alltoall_guard50_controlDict`
