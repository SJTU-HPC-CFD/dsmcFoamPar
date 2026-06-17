# ourmesh mpi8 alpha=0.8 Dual=false auto no-gate adaptive-alpha retest

Date: 2026-06-13

Source state:

- Added `replicatedMeshDLBAdaptiveAlpha` support in
  `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.{H,C}`.
- Default behavior remains fixed alpha because `replicatedMeshDLBAdaptiveAlpha`
  defaults to `false`.
- The controller uses inter-DLB move/collision wall-time deltas to adjust
  `replicatedMeshDLBAlpha` after a DLB action; the updated alpha is used by the
  next ParMETIS repartition.

Case:

`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Temporary controls:

- `replicatedMeshDLBAlpha 0.8`
- `replicatedMeshDLBDualConstraint false`
- `replicatedMeshDLBForceSteps ()`
- `replicatedMeshDLBParticleGate false`
- `replicatedMeshDLBAdaptiveAlpha true`
- `replicatedMeshDLBAlphaMin 0.5`
- `replicatedMeshDLBAlphaMax 2.0`
- `replicatedMeshDLBAlphaGain 0.04`
- `replicatedMeshDLBAlphaMaxStep 0.08`
- `replicatedMeshDLBAlphaWorsenTol 0.05`
- `endTime 5e-05` / 500 steps

Trigger and alpha result:

| trigger # | step | reason | alpha used | load imbalance | cells changed | alpha after |
|---:|---:|---|---:|---:|---:|---:|
| 1 | 50 | threshold | 0.8 | 1.683007338 | 91396 | 0.8 |
| 2 | 110 | SAR | 0.8 | 1.357011839 | 72 | 0.7921385742 |
| 3 | 160 | SAR | 0.7921385742 | 1.437129799 | 54303 | 0.7860321889 |
| 4 | 220 | SAR | 0.7860321889 | 1.604594849 | 32 | 0.8022580925 |
| 5 | 270 | threshold | 0.8022580925 | 1.59778853 | 93433 | 0.7825165223 |
| 6 | 370 | threshold | 0.7825165223 | 1.551546613 | 243 | 0.7677191159 |
| 7 | 490 | SAR | 0.7677191159 | 1.487753508 | 25 | 0.754550874 |

Summary confirms:

- `Phase C auto DLB checks = 500`
- `Phase C auto DLB rebalances = 7`
- `Phase C auto DLB triggered checks = 7`
- `Phase C adaptive alpha final = 0.754550874`

Performance:

| metric | value |
|---|---:|
| real | 77.78 s |
| execution time at step 500 | 77.62 s |
| full evolve wall | 73.76475971 s |
| move+collide wall | 72.3988907 s |
| move only | 57.58577547 s |
| buildCellOccupancy | 9.096421846 s |
| collision phase | 3.871498427 s |
| migration wall max | 3.153045241 s |
| particles max/min | 1.139564912 |
| rank wall max/min | 1.209390382 |
| DLB wall max | 14.19554355 s |
| Phase C ParMETIS max | 0.391758379 s |
| Phase C migration max | 0.81956522 s |

Correctness gate:

| metric | value |
|---|---:|
| exit | 0 |
| Total Iterations | 500 |
| final particles | 2463767 |
| final total energy | 1.24170012 |
| cumulative global collisions | 9503318 |
| cumulative global candidates | 16785320 |
| stuck particles | 0 |

Comparison:

| Variant | real [s] | full evolve [s] | move+collide [s] | move [s] | build occ [s] | collision [s] | DLB rebalances | DLB wall max [s] | particles max/min | rank wall max/min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| auto no-gate alpha=0.8 adaptive | 77.78 | 73.76475971 | 72.3988907 | 57.58577547 | 9.096421846 | 3.871498427 | 7 | 14.19554355 | 1.139564912 | 1.209390382 |
| auto no-gate alpha=0.8 fixed | 79.49 | 74.16347718 | 72.90211629 | 58.24886678 | 8.297763903 | 3.977779489 | 7 | 17.64523032 | 1.129305066 | 1.292758688 |
| auto no-gate alpha=1.0 fixed | 82.52 | 79.21398808 | 78.02151468 | 62.16768641 | 9.398137272 | 4.314510839 | 6 | 17.33056553 | 1.152432608 | 1.254160637 |
| forced8 alpha=1.0 fixed | 77.85 | 71.69948578 | 70.42431176 | 57.15336762 | 6.064381571 | 4.093368093 | 8 | 17.64633934 | 1.076510267 | 1.306471816 |

Interpretation:

1. Adaptive alpha improves the non-forced auto no-gate case versus fixed
   `alpha=0.8`: real time improves by 1.71 s, execution time by 2.22 s, and
   full evolve by 0.40 s.

2. The improvement is not just the final alpha value.  The trigger path changes:
   adaptive alpha triggers at 50, 110, 160, 220, 270, 370, and 490, while fixed
   alpha=0.8 triggered at 50, 110, 200, 300, 350, 400, and 480.  The adaptive
   run moves large repartitions earlier, especially steps 160 and 270.

3. The adaptive run reaches the same external real-time band as the current
   forced8 alpha=1.0 best run (`77.78 s` vs `77.85 s`), but it does not beat it
   internally: full evolve is still 73.76 s versus 71.70 s.  Treat the real-time
   tie as encouraging but not final proof.

4. The rank-wall ratio improves strongly to 1.209, but particle-count balance is
   slightly worse than fixed alpha=0.8.  This supports the current direction:
   optimize actual rank work, not just final particle max/min.

Artifacts:

- Control snapshots: `controlDicts/`
- Logs: `logs/`
- Restored case controlDict snapshot: `controlDicts/ourmesh_mpi8_controlDict_restored`
