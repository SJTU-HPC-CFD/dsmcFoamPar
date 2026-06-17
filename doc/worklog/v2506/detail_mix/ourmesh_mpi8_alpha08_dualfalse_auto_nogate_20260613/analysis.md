# ourmesh mpi8 alpha=0.8 Dual=false auto DLB no-gate retest

Date: 2026-06-13

Purpose:

Retest the `alpha=0.8, DualConstraint=false` case with non-forced automatic
DLB and particle gate disabled, matching the previous `alpha=1` auto no-gate
test except for `replicatedMeshDLBAlpha`.

Case:

`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Temporary controls:

- `replicatedMeshDLBAlpha 0.8`
- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBForceSteps ()`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBImbalanceThreshold 1.5`
- `endTime 5e-05` / 500 steps

Trigger result:

| trigger # | step | reason | load imbalance | cells changed |
|---:|---:|---|---:|---:|
| 1 | 50 | threshold | 1.696056119 | 91908 |
| 2 | 110 | SAR | 1.353868675 | 163 |
| 3 | 200 | SAR | 1.371187748 | 156 |
| 4 | 300 | SAR | 1.441535897 | 85901 |
| 5 | 350 | threshold | 1.632395365 | 124 |
| 6 | 400 | threshold | 1.510806908 | 276 |
| 7 | 480 | SAR | 1.574635019 | 62100 |

Summary confirms:

- `Phase C auto DLB checks = 500`
- `Phase C auto DLB rebalances = 7`
- `Phase C auto DLB triggered checks = 7`

Performance:

| metric | value |
|---|---:|
| real | 79.49 s |
| execution time at step 500 | 79.84 s |
| full evolve wall | 74.16347718 s |
| move+collide wall | 72.90211629 s |
| move only | 58.24886678 s |
| buildCellOccupancy | 8.297763903 s |
| collision phase | 3.977779489 s |
| migration wall max | 3.269669417 s |
| particles max/min | 1.129305066 |
| rank wall max/min | 1.292758688 |
| DLB wall max | 17.64523032 s |
| Phase C ParMETIS max | 0.357175849 s |
| Phase C migration max | 0.682885786 s |

Correctness gate:

| metric | value |
|---|---:|
| exit | 0 |
| Total Iterations | 500 |
| final particles | 2463759 |
| final total energy | 1.24180869 |
| cumulative global collisions | 9574183 |
| cumulative global candidates | 16765931 |
| stuck particles | 0 |

Comparison:

| Variant | real [s] | full evolve [s] | move+collide [s] | move [s] | build occ [s] | collision [s] | DLB rebalances | particles max/min | rank wall max/min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| auto no-gate alpha=0.8 Dual=false | 79.49 | 74.16347718 | 72.90211629 | 58.24886678 | 8.297763903 | 3.977779489 | 7 | 1.129305066 | 1.292758688 |
| auto no-gate alpha=1.0 Dual=false | 82.52 | 79.21398808 | 78.02151468 | 62.16768641 | 9.398137272 | 4.314510839 | 6 | 1.152432608 | 1.254160637 |
| forced8 alpha=0.8 Dual=false | 77.79 | 75.00281689 | 73.80264439 | 59.22950381 | 8.836414037 | 4.943804734 | 8 | 1.063836386 | 1.291486858 |
| forced8 alpha=1.0 Dual=false | 77.85 | 71.69948578 | 70.42431176 | 57.15336762 | 6.064381571 | 4.093368093 | 8 | 1.076510267 | 1.306471816 |

Interpretation:

1. Under non-forced auto no-gate, `alpha=0.8` is clearly better than
   `alpha=1.0`: real time improves from 82.52 s to 79.49 s, and full evolve
   improves from 79.21 s to 74.16 s.  The largest improvements are in move and
   buildCellOccupancy.

2. The performance improvement is linked to different automatic trigger timing,
   not merely to lower DLB overhead.  `alpha=0.8` triggers seven times and
   performs large repartitions at steps 50, 300, and 480, while the previous
   `alpha=1.0` auto no-gate run triggered six times and mostly made small
   changes after step 50.

3. The run still does not beat the best current setting, `forced8 alpha=1.0
   Dual=false`, which remains faster in real time and full evolve.  The forced8
   alpha=1.0 run has much lower buildCellOccupancy and lower critical-path rank
   wall time.

4. Compared with `forced8 alpha=0.8 Dual=false`, this auto no-gate run has a
   slightly better internal full-evolve profile but worse external real time.
   Treat that as noise and trigger-policy sensitivity until repeated.

5. Current conclusion: for non-forced auto no-gate, `alpha=0.8` is better than
   `alpha=1.0`; for the global best measured state, forced8 alpha=1.0 is still
   the target to reproduce with a real automatic trigger policy.

Artifacts:

- Control snapshots: `controlDicts/`
- Logs: `logs/`
- Restored case controlDict snapshot: `controlDicts/ourmesh_mpi8_controlDict_restored`
