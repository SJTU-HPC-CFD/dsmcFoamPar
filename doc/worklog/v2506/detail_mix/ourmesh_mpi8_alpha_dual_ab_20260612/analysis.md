# ourmesh mpi8 alpha / DualConstraint A-B

Date: 2026-06-12

Case:
`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Build:
`source doc/scripts/build-dsmcFoam.sh` completed before the runs.

Important source semantics:

- Current source does not implement a runtime-adaptive alpha controller.
- `replicatedMeshDLBAlpha` is read once in `dsmcReplicatedMesh::initialize()`.
- `ParMETIS AdaptiveRepart` is the repartitioning algorithm name; it is not adaptive alpha.
- `DualConstraint false`: one vertex weight, `N^alpha`.
- `DualConstraint true`: two vertex weights, `N^alpha` and `N*(N-1)`.

Common run controls:

- `useOpenMP false`
- `openmpThreads 1`
- `endTime 5e-05` / 500 steps
- `profileSummary true`
- `profileDetail false`
- `replicatedMeshNoAlltoall true`
- `replicatedMeshFlatTransfer true`
- `replicatedMeshDLBForceSteps ()`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBCheckCollective allgather`
- `replicatedMeshDLBImbalanceThreshold 1.5`
- `replicatedMeshDLBVsizeExp 0`

Runs:

| Variant | alpha | DualConstraint | real [s] | full evolve [s] | move+collide [s] | move [s] | build occ [s] | collision [s] | migration max [s] | final particles max/min | rank wall max/min | DLB rebalances |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| alpha1_dualfalse | 1.0 | false | 83.18 | 78.23545799 | 77.30572329 | 62.94900865 | 9.520176916 | 4.313123272 | 3.047877293 | 1.3618536 | 1.316286826 | 1 |
| alpha08_dualfalse | 0.8 | false | 88.03 | 84.84078044 | 83.98431521 | 69.21700316 | 10.13078389 | 4.66468599 | 3.125185615 | 1.271053458 | 1.336800139 | 1 |
| alpha08_dualtrue | 0.8 | true | 94.37 | 91.04214661 | 90.17481276 | 72.54206034 | 10.86472974 | 4.63362736 | 3.825390303 | 1.463624307 | 1.421079604 | 1 |

Trigger confirmation:

- `alpha1_dualfalse`: step 50, `alpha=1 ncon=1`, changed 104127 / 104151 cells.
- `alpha08_dualfalse`: step 50, `alpha=0.8 ncon=1`, changed 94713 / 104151 cells.
- `alpha08_dualtrue`: step 50, `alpha=0.8 ncon=2`, changed 96754 / 104151 cells.

Correctness gate:

| Variant | final particles | final total energy | cumulative global collisions | cumulative global candidates | stuck particles |
|---|---:|---:|---:|---:|---:|
| alpha1_dualfalse | 2463766 | 1.242498636 | 9661571 | 17537068 | 0 |
| alpha08_dualfalse | 2463773 | 1.24264311 | 9666254 | 17561306 | 0 |
| alpha08_dualtrue | 2463891 | 1.242637847 | 9665727 | 17573758 | 0 |

Interpretation:

1. `alpha=0.8` is not beneficial in this current code/control state.  It improves the final particle-count proxy imbalance from 1.3619 to 1.2711, but the actual rank wall imbalance gets slightly worse, from 1.3163 to 1.3368.  The wall time also increases: real time +5.8%, full evolve +8.4%, move+collide +8.6%.

2. `DualConstraint true` is worse for this case.  With `alpha=0.8`, enabling the second collision proxy constraint makes real time +7.2% slower than `alpha08_dualfalse` and +13.5% slower than `alpha1_dualfalse`.  It also worsens final particles max/min to 1.4636 and rank wall max/min to 1.4211.

3. The DLB decision path is not the main overhead difference in these three runs.  All variants trigger exactly one repartition at step 50.  `Phase C ParMETIS max` is only 0.059-0.084 s and `Phase C auto DLB wall max` is about 1.31-1.44 s.  The end-to-end differences mostly appear after repartition in move/build/migration and rank waiting behavior.

4. For the current `ourmesh` mpi8 case, the cleanest DLB weight setting is still the simplest one: `DualConstraint false`, `replicatedMeshDLBAlpha 1.0`.  The compressed `N^0.8` weight reduces particle-count spread, but it apparently underweights heavy particle cells enough that the actual compute wall becomes worse.

Recommended next state:

- Keep `DualConstraint false` for this case.
- Use `replicatedMeshDLBAlpha 1.0` as the current best clean alpha setting unless a repeated run set contradicts it.
- Do not describe the current code as adaptive alpha; it is fixed alpha plus ParMETIS AdaptiveRepart.

Artifacts:

- Control snapshots: `controlDicts/`
- Logs: `logs/`
- Final restored controlDict snapshot: `controlDicts/ourmesh_mpi8_controlDict_restored`
