# ourmesh mpi8 alpha=1 Dual=false auto DLB no-gate retest

Date: 2026-06-12

Purpose:

Retest the `alpha=1, DualConstraint=false` case with non-forced automatic DLB,
while preventing the current source's default particle gate from suppressing
automatic checks.

Case:

`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Temporary controls:

- `replicatedMeshDLBAlpha 1`
- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBForceSteps ()`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBImbalanceThreshold 1.5`
- `endTime 5e-05` / 500 steps

Trigger result:

| trigger # | step | reason | load imbalance |
|---:|---:|---|---:|
| 1 | 50 | threshold | 1.580024109 |
| 2 | 120 | SAR | 1.403008973 |
| 3 | 190 | SAR | 1.45007786 |
| 4 | 280 | SAR | 1.436682984 |
| 5 | 360 | SAR | 1.405377069 |
| 6 | 460 | SAR | 1.455940596 |

Summary confirms:

- `Phase C auto DLB checks = 500`
- `Phase C auto DLB rebalances = 6`
- `Phase C auto DLB triggered checks = 6`

Performance:

| metric | value |
|---|---:|
| real | 82.52 s |
| execution time at step 500 | 82.01 s |
| full evolve wall | 79.21398808 s |
| move+collide wall | 78.02151468 s |
| move only | 62.16768641 s |
| buildCellOccupancy | 9.398137272 s |
| collision phase | 4.314510839 s |
| migration wall max | 3.238523382 s |
| particles max/min | 1.152432608 |
| rank wall max/min | 1.254160637 |
| DLB wall max | 17.33056553 s |
| Phase C ParMETIS max | 0.34705024 s |
| Phase C migration max | 0.549230206 s |

Correctness gate:

| metric | value |
|---|---:|
| exit | 0 |
| Total Iterations | 500 |
| final particles | 2463707 |
| final total energy | 1.242442422 |
| cumulative global collisions | 9658420 |
| stuck particles | 0 |

Interpretation:

1. Non-forced auto DLB with particle gate disabled does not collapse to one
   trigger.  It triggered six times in 500 steps.

2. This still differs from strict forced8.  The report's eight-rebalance
   benchmark corresponds to the forced schedule, while this run follows the
   current legacy auto/SAR logic.

3. Compared with the one-trigger `alpha=1, Dual=false` run, this auto no-gate
   run is slightly faster in external real time (82.52 s vs 83.18 s), but slower
   in full evolve (79.21 s vs 78.24 s).  The important correction is not the
   small timing delta; it is that the previous one-trigger run was not a valid
   normal-DLB test.

4. Compared with forced8 `alpha=1, Dual=false`, auto no-gate is slower:
   `real 82.52 s` vs `77.85 s`, and `full evolve 79.21 s` vs `71.70 s`.

5. If the target is "about eight DLB actions in 500 steps" without force steps,
   the current legacy auto/SAR trigger policy is still too sparse.  It needs a
   real automatic trigger-policy change, not just disabling the particle gate.

Artifacts:

- Control snapshots: `controlDicts/`
- Logs: `logs/`
- Restored case controlDict snapshot: `controlDicts/ourmesh_mpi8_controlDict_restored`
