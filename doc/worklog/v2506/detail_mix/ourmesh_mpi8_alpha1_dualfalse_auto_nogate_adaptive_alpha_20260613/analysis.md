# ourmesh mpi8 alpha=1 Dual=false auto no-gate adaptive-alpha retest

Date: 2026-06-13

Purpose:

Test `replicatedMeshDLBAlpha 1` with `DualConstraint false`, automatic DLB,
particle gate disabled, and adaptive alpha enabled.

Case:

`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Source/control state:

- Current `dsmcFoam+` binary timestamp is newer than the adaptive-alpha source
  edit, so this run uses the adaptive-alpha implementation.
- `replicatedMeshDLBAdaptiveAlpha` defaults to `false`; this run explicitly
  enables it.
- The production case `controlDict` was restored after the run to the previous
  `replicatedMeshDLBAlpha 0.8` state.

Temporary controls:

- `replicatedMeshDLBAlpha 1`
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
| 1 | 50 | threshold | 1 | 1.62023377 | 92732 | 1 |
| 2 | 140 | SAR | 1 | 1.436387102 | 79 | 0.9855748902 |
| 3 | 200 | SAR | 0.9855748902 | 1.504136028 | 72523 | 1.004095329 |
| 4 | 250 | threshold | 1.004095329 | 1.73127502 | 150 | 0.980434722 |
| 5 | 300 | threshold | 0.980434722 | 1.807412559 | 38 | 1.006895056 |
| 6 | 350 | threshold | 1.006895056 | 1.714436198 | 16 | 0.9802192818 |
| 7 | 400 | threshold | 0.9802192818 | 1.763183819 | 93346 | 0.9517299959 |
| 8 | 450 | threshold | 0.9517299959 | 1.570077347 | 69 | 0.9407682276 |
| 9 | 500 | threshold | 0.9407682276 | 1.68288764 | 14 | 0.967545961 |

Summary confirms:

- `Phase C auto DLB checks = 500`
- `Phase C auto DLB rebalances = 9`
- `Phase C auto DLB triggered checks = 9`
- `Phase C adaptive alpha final = 0.967545961`

Performance:

| metric | value |
|---|---:|
| real | 86.63 s |
| execution time at step 500 | 85.60 s |
| full evolve wall | 81.08255063 s |
| move+collide wall | 79.90594767 s |
| move only | 63.31484159 s |
| buildCellOccupancy | 9.80492273 s |
| collision phase | 3.601251415 s |
| migration wall max | 3.896504854 s |
| particles max/min | 1.161694967 |
| rank wall max/min | 1.438447544 |
| DLB wall max | 25.49498549 s |
| Phase C check max | 24.38335062 s |
| Phase C ParMETIS max | 0.501525304 s |
| Phase C migration max | 0.868543655 s |
| Phase C post-diagnostic max | 0.001383885 s |

Correctness gate:

| metric | value |
|---|---:|
| exit | 0 |
| Total Iterations | 500 |
| final particles | 2463719 |
| final total energy | 1.241255987 |
| cumulative global collisions | 9456192 |
| cumulative global candidates | 16322398 |
| stuck particles | 0 |

Comparison:

| Variant | real [s] | full evolve [s] | move+collide [s] | move [s] | build occ [s] | collision [s] | DLB rebalances | DLB wall max [s] | particles max/min | rank wall max/min |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| auto no-gate alpha=1.0 adaptive | 86.63 | 81.08255063 | 79.90594767 | 63.31484159 | 9.80492273 | 3.601251415 | 9 | 25.49498549 | 1.161694967 | 1.438447544 |
| auto no-gate alpha=1.0 fixed | 82.52 | 79.21398808 | 78.02151468 | 62.16768641 | 9.398137272 | 4.314510839 | 6 | 17.33056553 | 1.152432608 | 1.254160637 |
| auto no-gate alpha=0.8 adaptive | 77.78 | 73.76475971 | 72.3988907 | 57.58577547 | 9.096421846 | 3.871498427 | 7 | 14.19554355 | 1.139564912 | 1.209390382 |
| forced8 alpha=1.0 fixed | 77.85 | 71.69948578 | 70.42431176 | 57.15336762 | 6.064381571 | 4.093368093 | 8 | 17.64633934 | 1.076510267 | 1.306471816 |

Interpretation:

1. `alpha=1` plus adaptive alpha is worse than fixed `alpha=1` in this run:
   real time regresses from 82.52 s to 86.63 s, and full evolve regresses from
   79.21 s to 81.08 s.

2. The regression is not due to ParMETIS compute itself.  ParMETIS max is only
   0.50 s.  The larger cost is the auto-DLB/check/rebalance path: DLB wall max
   rises to 25.49 s, and the summary reports 9 rebalances instead of the fixed
   run's 6.

3. The adaptive controller moved alpha only from 1.0 to 0.9675 by the end of
   500 steps and oscillated around 1.0.  Starting at alpha=1 therefore does not
   reach the better `alpha=0.8 adaptive` region within this run.

4. Rank-wall balance also got worse (`1.438`), despite the lower collision
   phase time.  The move/build/DLB path dominates the loss.

5. Current best among these auto no-gate adaptive tests remains
   `alpha=0.8, Dual=false, adaptiveAlpha=true`: real 77.78 s, full evolve
   73.76 s, and rank wall max/min 1.209.

Artifacts:

- Control snapshots: `controlDicts/`
- Logs: `logs/`
- Restored case controlDict snapshot: `controlDicts/ourmesh_mpi8_controlDict_restored`
