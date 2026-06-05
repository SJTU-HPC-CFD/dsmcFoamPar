# Stage 0 Baseline Source Audit - 2026-06-05

## Scope

- Working tree: `/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb`
- Read-only reference: `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
- Build log: `doc/worklog/v2506/detail/stage0_build_20260605.log`
- Baseline case root: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp`

## Build Status

`source doc/scripts/build-dsmcFoam.sh` completed successfully on 2026-06-05.

Built targets:

- `src/lagrangian/basic` -> `liblagrangian+`
- `src/lagrangian/molecularDynamics/general` -> `libgeneralMolecule`
- `src/lagrangian/dsmc` -> `libdsmcFoam+`
- `applications/solvers/discreteMethods/dsmc/dsmcFoam+`
- `applications/utilities/preProcessing/dsmc/dsmcInitialise+`

The build output still reports older binary mtimes for `dsmcFoam+` and `dsmcInitialise+`; the library rebuild did compile `dsmcVolFields.C` and relink `libdsmcFoam+.so`. Treat this as a successful compile check, not as proof that every solver object was rebuilt from source.

## Current Source Versus Reference

| Module | Current OFv1706 | Reference hyStrath_xcx | Stage 0 finding |
|---|---:|---:|---|
| `src/lagrangian/basic/Cloud/Cloud.C` | 591 lines | 2020 lines | Current tree lacks the reference OpenMP move path, replicated mesh guards, ordered parcel extraction, and move profiling helpers. |
| `src/lagrangian/dsmc/clouds/dsmcCloud.C` | 2829 lines | 5434 lines | Current tree lacks the reference replicated mesh lifecycle, profile controls, OpenMP collision/move controls, flat/moveOrdered occupancy support, and DLB reporting. |
| `src/lagrangian/dsmc/clouds/dsmcCloud.H` | 811 lines | 686 lines | API shape differs; direct copy is unsafe. Need field-by-field port and compile after each small slice. |
| `src/lagrangian/dsmc/clouds/dsmcCloudI.H` | 394 lines | 523 lines | Reference has replicated mesh and additional inline accessors; current tree has no matching replicated mesh accessors. |
| `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C` | 355 lines | 580 lines | Reference includes optimized/profiled collision path; current tree is baseline NTC. |
| `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.H` | 136 lines | 147 lines | Minor API drift; port must preserve OFv1706 constructor and runtime-selection compatibility. |
| `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C` | 2614 lines | 4901 lines | Current file already has local uncommitted edits; do not overwrite. Reference has replicated-output and flat occupancy/post-field optimizations. |
| `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H` | 262 lines | 295 lines | Current file already has local uncommitted edits; post-field port must be merged carefully. |
| `src/lagrangian/dsmc/collisions/derived/VariableHardSphere/VariableHardSphere.C` | 280 lines | 289 lines | Reference adds small collision model changes; verify `sigmaTcR()` API before porting cr2 lookup. |
| `src/lagrangian/dsmc/collisions/derived/VariableHardSphere/VariableHardSphere.H` | 159 lines | 168 lines | Header has small drift; do not copy blindly. |
| `src/lagrangian/dsmc/dynamicLoadBalancing/dsmcDynamicLoadBalancing.C` | 288 lines | 1848 lines | Current tree only has the older dynamic load balancing implementation. Reference contains large MPI/DLB work that must be staged after OpenMP basics. |
| `src/lagrangian/dsmc/dynamicLoadBalancing/dsmcDynamicLoadBalancing.H` | 134 lines | 227 lines | DLB API has expanded in reference. |
| `applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C` | 183 lines | 280 lines | Current solver lacks replicated mesh gather/write/migrate/cellOwner flow. |
| `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C` | missing | 2126 lines | Replicated mesh implementation has not been ported into current OFv1706 source tree. |
| `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H` | missing | 180 lines | Replicated mesh API has not been ported into current OFv1706 source tree. |

## Feature Probe

Current OFv1706 source has 112 matches for the broad feature probe, but these are mostly ordinary `sigmaTcR` symbols and baseline collision code. It does not contain the reference feature set:

- no `dsmcReplicatedMesh` source directory;
- no `replicatedMeshActive()` implementation in `dsmcCloud`;
- no solver-side replicated gather/write/migrate flow;
- no `profileSummary` / `profileDetail` controlDict support in current `dsmcCloud`;
- no `openmpMoveSchedule` / `openmpCollisionSchedule` controls;
- no `moveOrderedParcels_` data flow;
- no METIS/ParMETIS link line in current `src/lagrangian/dsmc/Make/options`.

Reference `src/lagrangian/dsmc/Make/files` includes:

- `replicatedMesh/dsmcReplicatedMesh.C`
- `replicatedMesh/dsmcLocalMesh.C`

Reference `Make/options` links:

- `-L$(PARMETIS_DIR)/lib -lparmetis -lmetis`
- solver link also includes `-lparmetis -lmetis -lGKlib`

Current OFv1706 `Make/options` only links decomposition methods and has no METIS/ParMETIS path. Stage 3 must first establish the link environment before adding replicated mesh files.

## Dirty Worktree Guard

The working tree already contains unrelated or previous edits:

- modified `doc/scripts/build-dsmcFoam.sh`
- modified `doc/worklog/v1706/tvib_min_parcels_threshold_20260528.md`
- modified `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`
- modified `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H`
- untracked v2506 worklog files

Do not overwrite the two `dsmcVolFields` files when porting post-field optimizations. Read and merge them explicitly.

## Stage 0 Gate

- Current source compiles.
- Baseline logs are present and parseable.
- The current source is not equivalent to the logged replicated/OpenMP-capable baseline; treat `ourmeshbkp` logs as historical reference results, and treat the current source as a clean compile base for staged porting.

Next implementation step: port only the low-disturbance profiling controls and timing sinks first, then rebuild and run a short smoke case before moving OpenMP data flow.
