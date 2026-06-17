# zb mpi8 single_alpha1 repeat3, 2026-06-13

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Configuration:

- `mpirun -np 8 dsmcFoam+`
- 300 steps: `endTime 1.9920146682e-05`, `deltaT 6.640048894e-08`
- `useOpenMP false`, `openmpThreads 1`
- replicated mesh raw MPI
- `replicatedMeshAutoDLB true`
- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBAlpha 1`
- `replicatedMeshDLBForceSteps ()`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha false`
- `replicatedMeshDLBCheckCollective allgather`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBVsizeExp 0`

The production case `controlDict` was restored after the three runs.

## Per-Run Results

| run | exit | real [s] | execution@300 [s] | full evolve [s] | move+collide [s] | move [s] | buildOcc [s] | collision [s] | post [s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rep1 | 0 | 69.12 | 69.77 | 69.1902746 | 68.91241458 | 45.99229342 | 4.622833415 | 10.60012566 | 0.281043256 |
| rep2 | 0 | 66.87 | 65.60 | 65.62081459 | 65.35283173 | 43.64349234 | 4.758788429 | 12.6954848 | 0.270564739 |
| rep3 | 0 | 68.95 | 67.50 | 67.89459903 | 67.62100362 | 46.41212793 | 5.091872721 | 9.325962611 | 0.282065301 |

## Statistics

| metric | mean | stdev | min | max |
|---|---:|---:|---:|---:|
| real [s] | 68.31333333 | 1.252850084 | 66.87 | 69.12 |
| execution@300 [s] | 67.62333333 | 2.087734019 | 65.60 | 69.77 |
| full evolve [s] | 67.56856274 | 1.806927209 | 65.62081459 | 69.1902746 |
| move+collide [s] | 67.29541664 | 1.801988534 | 65.35283173 | 68.91241458 |
| move [s] | 45.34930456 | 1.492116537 | 43.64349234 | 46.41212793 |
| buildOcc [s] | 4.824498188 | 0.241325086 | 4.622833415 | 5.091872721 |
| collision [s] | 10.87385769 | 1.701357359 | 9.325962611 | 12.6954848 |
| migration wall max [s] | 3.66107507 | 0.340171037 | 3.26832111 | 3.862409834 |

## DLB and Balance

| run | checks | rebalances | trigger steps | DLB wall max [s] | check max [s] | ParMETIS max [s] | DLB migration max [s] | particles max/min | rank wall max/min |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| rep1 | 300 | 5 | 100,150,200,250,300 | 20.25654738 | 19.39734315 | 0.174367132 | 0.739651917 | 1.050202923 | 1.31013236 |
| rep2 | 300 | 6 | 50,100,150,200,250,300 | 17.99137675 | 17.06098416 | 0.20107497 | 0.762658504 | 1.109332284 | 1.310643275 |
| rep3 | 300 | 5 | 100,150,200,250,300 | 18.5836826 | 17.76251481 | 0.184259868 | 0.695562628 | 1.139935873 | 1.31771068 |
| mean | 300 | 5.333333333 | - | 18.94386891 | 18.07361404 | 0.186567323 | 0.73262435 | 1.099823693 | 1.312828772 |

## Correctness Gate

| run | iterations | final particles | final energy | cumulative global collisions | cumulative global candidates | stuck |
|---|---:|---:|---:|---:|---:|---:|
| rep1 | 300 | 1958921 | 0.001952100143 | 40196870 | 73110906 | 0 |
| rep2 | 300 | 1958449 | 0.001954851589 | 40939438 | 74545204 | 0 |
| rep3 | 300 | 1958437 | 0.001955759067 | 40998611 | 74493988 | 0 |

## Comparison

| Source | real [s] | full evolve [s] | collision [s] | note |
|---|---:|---:|---:|---|
| 2026-06-10 clean MPI8 best single run | 78.20 | 78.87242265 | 30.2890735 | previous formal repeated best |
| 2026-06-10 clean MPI8 mean | 79.81 | 81.33 | 33.73 | 3-run mean |
| 2026-06-12 current-source rank-detail | 77.61 | 79.5623 | 31.7615 | no candidate-load gather |
| 2026-06-13 matrix `single_alpha1` | 71.87 | 70.1149057 | 13.18486084 | first single run |
| 2026-06-13 repeat3 `single_alpha1` mean | 68.31333333 | 67.56856274 | 10.87385769 | this retest |

Repeat3 improves over the previous clean repeated best single run:

- real: `78.20 -> 68.31 s`, `-9.89 s` / `-12.6%`
- full evolve: `78.87 -> 67.57 s`, `-11.30 s` / `-14.3%`
- collision: `30.29 -> 10.87 s`, `-19.42 s` / `-64.1%`

Repeat3 also improves over the 2026-06-13 matrix single run:

- real: `71.87 -> 68.31 s`, `-3.56 s` / `-5.0%`
- full evolve: `70.11 -> 67.57 s`, `-2.55 s` / `-3.6%`
- collision: `13.18 -> 10.87 s`, `-2.31 s` / `-17.5%`

## Interpretation

1. `single_alpha1` is now confirmed as a stable fast `zb` MPI8 setting in the
   current source/control state.  The three-repeat mean is faster than the
   previous matrix single run and much faster than the previous clean repeated
   MPI8 best.

2. The remaining DLB concern is not ParMETIS compute time.  ParMETIS max is
   about `0.19 s`, while DLB check max is still about `18.07 s` mean.  The next
   DLB optimization should target check cadence/collective overhead and the
   zero-payoff final-step rebalance at step 300.

3. Rebalance count is stable at 5-6 actions, and all runs still rebalance at
   step 300.  A no-final-DLB guard or forced schedule without the terminal
   rebalance should be the next direct check.

## Artifacts

- Logs: `logs/`
- Control snapshots: `controlDicts/`
- Restored case controlDict: `controlDicts/zb_mpi8_controlDict_restored`
