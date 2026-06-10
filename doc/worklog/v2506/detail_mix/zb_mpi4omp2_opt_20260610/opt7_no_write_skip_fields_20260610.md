# opt7 no-write field skip record, 2026-06-10

## Formal run scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi4omp2`
- Command: `OMP_NUM_THREADS=2 mpirun -np 4 dsmcFoam+`
- Do not add `-parallel` for this replicated-mesh formal path.
- Steps: 300 (`deltaT 6.640048894e-08`, `endTime 1.9920146682e-05`)
- Output mode: formal no-write (`writeControl runTime`, `writeInterval 1.e-3`)
- Replicated mesh: `replicatedMesh true`
- Logs: `doc/worklog/v2506/detail_mix/zb_mpi4omp2_opt_20260610/logs/`

## Source decision

opt7 is kept on top of opt3 and opt4.

- opt3 keeps `writeRatesToTerminal no` truly silent while preserving per-step reaction counter reset.
- opt4 rebuilds cell occupancy after reactions only when this step actually reacted.
- opt7 skips measurement field calculation/write in `dsmcCloud::evolve_fields()` only when the whole run has no scheduled output and the binary collision model is not using the `inverseZvFormulation 2008` macroscopic temperature feedback path.

The opt7 guard is intentionally narrow:

- `runTime.outputTime()` disables the skip for actual output steps.
- Only `writeControl runTime` and `adjustableRunTime` are eligible.
- `writeInterval` must be larger than the full run span plus a half-`deltaT` margin.
- `inverseZvFormulation 2008` disables the skip because `dsmcCloud::relaxationCollision()` can read `fields().overallT(cellI)`.

Current `mpi4omp2` inputs satisfy the skip boundary:

- `controlDict`: `writeInterval 1.e-3`, larger than the 300-step run span.
- `dsmcProperties`: `BinaryCollisionModel LarsenBorgnakkeVariableHardSphere`, no explicit `inverseZvFormulation 2008`.
- `controllersDict`: both state and flux controller lists are empty.
- `fieldPropertiesDict`: six `dsmcVolFields` entries only (`N2`, `O2`, `NO`, `N`, `O`, `mixture`).
- Replicated-mesh `autoRebalance()` runs before `evolve_fields()` and does not consume `dsmcVolFields` output.

## Performance result

| state | real [s] | full evolve wall [s] | total profiled [s] | buildCellOccupancy [s] | collision phase [s] | post fields/output [s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 83.80 | 80.70886171 | 90.94736673 | 7.85489927 | 25.65323808 | 12.61317843 |
| opt4 | 77.20 | 73.54404033 | 77.48239900 | 3.469258241 | 22.17539249 | 11.81203185 |
| opt7 | 67.00 | 64.41471091 | 73.29035397 | 4.190947921 | 25.17229504 | 0.561783891 |

Relative to opt4, opt7 gives:

- `real`: 77.20 -> 67.00 s, -10.20 s (-13.2%).
- `full evolve wall`: 73.54404033 -> 64.41471091 s, -9.12932942 s (-12.4%).
- `post fields/output`: 11.81203185 -> 0.561783891 s, -11.250247959 s (-95.2%).
- `total profiled`: 77.482399 -> 73.29035397 s, -4.19204503 s (-5.4%).

`buildCellOccupancy` and `collision phase` are higher than the opt4 run, but the large `post fields/output` reduction still produces a clear end-to-end win. The replicated-mesh DLB trigger steps differ between formal runs, so move/collision subphase values are not bitwise-stable comparison targets.

## Correctness and run health

All listed runs reached `Total Iterations = 300`, `End main`, and the opt7 log ended with `exit 0`.

| state | final particles | stuck | collisions | candidates | total energy |
| --- | ---: | ---: | ---: | ---: | ---: |
| baseline | 1959546 | 0 | 233200 | 454408 | 0.001956581282 |
| opt4 | 1959186 | 0 | 221590 | 426826 | 0.001953886199 |
| opt7 | 1959580 | 0 | 226763 | 432748 | 0.001955012694 |

Against opt4, opt7 changes final particles by +394 (+0.020%) and total energy by about +0.058%. Collision and candidate counts also differ, but this is consistent with the current replicated-mesh DLB path not being bitwise invariant across formal runs:

- baseline DLB triggers: steps 120, 240
- opt3 DLB triggers: steps 170, 270
- opt4 DLB triggers: steps 140, 220
- opt7 DLB triggers: steps 100, 190

The current acceptance criterion is therefore not bitwise equality of stochastic/DLB counters, but successful completion, no stuck particles, same 300-step/no-write control state, comparable particle/energy totals, and a profile-supported reduction of the targeted field-statistics cost.

## Rejected changes

The earlier opt1, opt5, and opt6 attempts remain rejected and are not part of the kept source state. Do not reintroduce them when continuing from this record.

