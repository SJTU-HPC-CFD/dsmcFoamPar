# zb mpi8 DLB strategy matrix, 2026-06-13

## Scope

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Run mode:

- `mpirun -np 8 dsmcFoam+`
- 300 steps: `endTime 1.9920146682e-05`, `deltaT 6.640048894e-08`
- replicated mesh raw MPI
- `useOpenMP false`, `openmpThreads 1`
- `profileSummary true`, `profileDetail false`
- `replicatedMeshAutoDLB true`
- `replicatedMeshDLBForceSteps ()`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha false`
- `replicatedMeshDLBVsizeExp 0`

Matrix:

| Variant | DualConstraint | alpha | ncon |
|---|---|---:|---:|
| single_alpha08 | false | 0.8 | 1 |
| single_alpha1 | false | 1.0 | 1 |
| dual_alpha08 | true | 0.8 | 2 |
| dual_alpha1 | true | 1.0 | 2 |

All runs exited with code 0.  The case `controlDict` was restored after the
matrix run.

## Performance

| Variant | real [s] | execution@300 [s] | full evolve [s] | move+collide [s] | move [s] | buildOcc [s] | collision [s] | post [s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| single_alpha08 | 76.04 | 76.28 | 76.39094563 | 76.11197746 | 46.55428638 | 5.021724436 | 15.73732304 | 0.282117882 |
| single_alpha1 | 71.87 | 70.13 | 70.1149057 | 69.83898843 | 48.25538793 | 4.543092916 | 13.18486084 | 0.282188106 |
| dual_alpha08 | 76.66 | 74.74 | 75.68021614 | 75.41372781 | 47.06980561 | 5.556405264 | 8.493558745 | 0.282059611 |
| dual_alpha1 | 76.45 | 76.45 | 77.82488867 | 77.55081349 | 48.08009741 | 6.406645958 | 11.26681842 | 0.277000779 |

## DLB and Balance

| Variant | checks | rebalances | trigger steps | DLB wall max [s] | check max [s] | ParMETIS max [s] | DLB migration max [s] | migration wall max [s] | particles max/min | rank wall max/min |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| single_alpha08 | 300 | 5 | 100,150,200,250,300 | 29.79827948 | 28.92160611 | 0.200458618 | 0.730128052 | 4.41678672 | 1.387862198 | 1.497098646 |
| single_alpha1 | 300 | 6 | 50,100,150,200,250,300 | 23.15632431 | 22.19933865 | 0.223824617 | 0.804754436 | 3.850086869 | 1.088705342 | 1.36318408 |
| dual_alpha08 | 300 | 5 | 100,150,200,250,300 | 27.22620276 | 26.34259062 | 0.249000443 | 0.669413147 | 4.1863653 | 1.380494667 | 1.325046041 |
| dual_alpha1 | 300 | 6 | 50,100,150,200,250,300 | 26.03101864 | 24.96518877 | 0.346710019 | 0.777761428 | 4.235959699 | 1.115388119 | 1.268052282 |

## Correctness Gate

| Variant | iterations | final particles | final energy | cumulative global collisions | cumulative global candidates | stuck |
|---|---:|---:|---:|---:|---:|---:|
| single_alpha08 | 300 | 1958559 | 0.00195802461 | 41155540 | 75495600 | 0 |
| single_alpha1 | 300 | 1959502 | 0.001954727125 | 40365506 | 74833799 | 0 |
| dual_alpha08 | 300 | 1959047 | 0.001926009554 | 38165471 | 69976619 | 0 |
| dual_alpha1 | 300 | 1959958 | 0.001927133244 | 38141889 | 70587526 | 0 |

## Trigger Detail

| Variant | step | reason | load imbalance | cells changed |
|---|---:|---|---:|---:|
| single_alpha08 | 100 | threshold | 1.579621133 | 51027 |
| single_alpha08 | 150 | threshold | 1.892369649 | 42386 |
| single_alpha08 | 200 | threshold | 1.872885093 | 24022 |
| single_alpha08 | 250 | threshold | 1.99793543 | 44061 |
| single_alpha08 | 300 | threshold | 2.400333335 | 58138 |
| single_alpha1 | 50 | threshold | 1.593948863 | 54656 |
| single_alpha1 | 100 | SAR+threshold | 1.97131276 | 54427 |
| single_alpha1 | 150 | SAR+threshold | 1.94769748 | 52960 |
| single_alpha1 | 200 | SAR+threshold | 1.794343102 | 24730 |
| single_alpha1 | 250 | threshold | 2.040498237 | 33727 |
| single_alpha1 | 300 | SAR+threshold | 1.748246606 | 58577 |
| dual_alpha08 | 100 | threshold | 1.533757695 | 50509 |
| dual_alpha08 | 150 | SAR+threshold | 2.105280716 | 59620 |
| dual_alpha08 | 200 | SAR+threshold | 1.913087066 | 22 |
| dual_alpha08 | 250 | SAR+threshold | 2.024320449 | 59301 |
| dual_alpha08 | 300 | SAR+threshold | 1.8713309 | 16594 |
| dual_alpha1 | 50 | threshold | 1.549866987 | 47815 |
| dual_alpha1 | 100 | SAR+threshold | 1.667805848 | 41183 |
| dual_alpha1 | 150 | SAR+threshold | 1.748320515 | 9263 |
| dual_alpha1 | 200 | SAR+threshold | 1.948406569 | 13272 |
| dual_alpha1 | 250 | SAR+threshold | 1.587537611 | 16812 |
| dual_alpha1 | 300 | SAR+threshold | 2.165079218 | 33108 |

## Interpretation

1. Best end-to-end strategy in this matrix is `single_alpha1`.
   It has the best external real time (`71.87 s`) and best full evolve time
   (`70.11 s`).  It also has the best final particle-count balance
   (`1.089 max/min`) among the four variants.

2. Dual constraint does reduce collision phase on `zb`.
   `dual_alpha08` gives the lowest collision phase (`8.49 s`) versus
   `single_alpha1` (`13.18 s`) and `single_alpha08` (`15.74 s`).  This supports
   the idea that the second `N*(N-1)` constraint is physically meaningful for
   a collision-heavier case.

3. The dual-constraint collision gain does not pay back in total time here.
   `dual_alpha08` loses in move/build/DLB enough that full evolve stays
   `75.68 s`, about `5.57 s` slower than `single_alpha1`.  `dual_alpha1` has
   the best rank-wall ratio (`1.268`) but the worst full evolve time
   (`77.82 s`), so rank-wall ratio alone is not a sufficient objective.

4. For single constraint, `alpha=1` is clearly better than `alpha=0.8` on zb.
   It reduces real time by `4.17 s`, full evolve by `6.28 s`, collision by
   `2.55 s`, DLB wall by `6.64 s`, and particle max/min from `1.388` to
   `1.089`.

5. For dual constraint, `alpha=0.8` has a better internal profile than
   `alpha=1` even though external real time is close.  The full-evolve numbers
   are `75.68 s` vs `77.82 s`; the small real-time difference (`76.66 s` vs
   `76.45 s`) should not be over-interpreted from one run.

6. The `zb` result differs from `ourmesh`: here the cleanest automatic-DLB
   strategy is not `alpha=0.8`, but `single_alpha1`.  The likely explanation is
   structural: `zb` has enough collision work that a raw particle-count balance
   avoids underweighting dense cells, while the dual constraint is still too
   expensive as a hard second balance constraint.

## Comparison With Previous MPI8 Best

The previous formal repeated MPI8 result for `zb` is from
`repeat_nowrite_compute_20260610_parmetisfix_clean`, 3 runs:

| Source | Variant | real [s] | full evolve [s] | move [s] | buildOcc [s] | collision [s] | DLB rebalances |
|---|---|---:|---:|---:|---:|---:|---:|
| 2026-06-10 clean repeated best single run | `zb_MPI8_rep1` | 78.20 | 78.87242265 | 49.29770603 | 5.986395256 | 30.2890735 | 2 |
| 2026-06-10 clean repeated mean | `MPI8 mean` | 79.81 | 81.33 | 54.01 | 6.51 | 33.73 | 2.33 |
| 2026-06-12 current-source detail, no candidate-load gather | `MPI8 rank-detail` | 77.61 | 79.5623 | n/a | n/a | 31.7615 | n/a |
| 2026-06-13 matrix best | `single_alpha1` | 71.87 | 70.1149057 | 48.25538793 | 4.543092916 | 13.18486084 | 6 |

Against the previous clean repeated best single run (`78.20 s`), this matrix
best improves:

- real: `78.20 -> 71.87 s`, `-6.33 s` / `-8.1%`
- full evolve: `78.87 -> 70.11 s`, `-8.76 s` / `-11.1%`
- collision phase: `30.29 -> 13.18 s`, `-17.10 s` / `-56.5%`
- build occupancy: `5.99 -> 4.54 s`, `-1.44 s` / `-24.1%`

Against the 2026-06-12 current-source rank-detail run (`77.61 s`), this matrix
best improves:

- real: `77.61 -> 71.87 s`, `-5.74 s` / `-7.4%`
- full evolve: `79.56 -> 70.11 s`, `-9.45 s` / `-11.9%`
- collision phase: `31.76 -> 13.18 s`, `-18.58 s` / `-58.5%`

This is the current fastest observed `zb` MPI8 replicated-mesh single run.
It should still be treated as a single-run best until repeated, because the
historical MPI8 `zb` runs had nontrivial run-to-run spread.

## Artifacts

- Logs: `logs/`
- Control snapshots: `controlDicts/`
- Restored case controlDict: `controlDicts/zb_mpi8_controlDict_restored`
