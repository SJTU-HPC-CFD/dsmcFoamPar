# MPI Replicated Mesh DLB Port/Gap Analysis - 2026-06-06

Scope:

- Current tree: `/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb`
- Reference source: `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
- Question: relative to the reference source, what has the current v1706 tree implemented for MPI replicated mesh DLB, and what is still missing.

## Executive conclusion

The current v1706 tree is no longer at the original "missing replicatedMesh module" state. It has a real `dsmcReplicatedMesh` port, is compiled into the DSMC library, is connected into `dsmcCloud`, and has a replicated output path in `dsmcFoam+.C`. The core Phase A/B/C machinery exists: static owner decomposition, particle migration, manual rebalance, ParMETIS AdaptiveRepart auto-DLB, SAR/legacyWindow-style triggering, and replicated mesh profiling.

However, the port is not yet equivalent to the reference implementation. The largest remaining gaps are outside the core `dsmcReplicatedMesh` class:

1. `dsmcVolFields` only has a partial replicated reduce patch. It does not yet match the reference owned-cell / owned-boundary-face accumulation path.
2. `replicatedMeshDelayedReceive true` exists in the reference controlDict and the reference `dsmcCloud`, but the current v1706 `dsmcCloud` does not wire that async begin/finish path.
3. The current `ourmesh/mpi8replicatedmesh/system/controlDict` is not configured for replicated mesh DLB; the full DLB settings are visible in `ourmeshbkp`.
4. Existing replicated mesh logs are still under `doc/worklog/v2506/detail_omp/`, not the new `detail_mpi/` convention.
5. Current 500-step owner-filter run is correct enough on particle count, but still slower than the reference/history and has a large post-field/output cost.

## Implemented in the current v1706 tree

### 1. Source module and build integration

The current tree contains the replicated mesh source files:

- `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H`
- `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C`
- `src/lagrangian/dsmc/replicatedMesh/dsmcLocalMesh.H`
- `src/lagrangian/dsmc/replicatedMesh/dsmcLocalMesh.C`

They are included by `src/lagrangian/dsmc/Make/files:19-20`:

- `replicatedMesh/dsmcLocalMesh.C`
- `replicatedMesh/dsmcReplicatedMesh.C`

The `dsmcReplicatedMesh.H` public/private interface is effectively the reference interface: the current file has the same line count as the reference and `cmp` reports identical contents. The important API surface is visible at `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.H:45-174`:

- ownership state: `cellOwner_`, `myCells_`, `localMesh_`
- migration state: `migrateInterval_`, `migrateBegin()`, `migrateFinish()`, `migrateParticlesByCellOwner()`
- Phase C DLB state: `autoDLBEnabled_`, `sarSteps_`, `autoDLBTriggerMode_`, `autoRebalanceChecks_`, `autoRebalanceCount_`
- ParMETIS path: `reassignByParMetisAdaptiveRepart()`
- output/profile path: `writeCellOwner()`, `report()`, `gatherParcelsToRank0()`

### 2. Core replicated mesh / DLB implementation

The current `dsmcReplicatedMesh.C` implements the same major algorithmic stages as the reference:

- MPI setup, initial owner construction, local mesh view, controlDict parsing: `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:235-357`
- ParMETIS AdaptiveRepart DLB: `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:364-624`
- automatic DLB trigger/check path: `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:631-760`
- synchronous migration and flat transfer: `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:1178-1325`
- async migration implementation exists in the class: `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:1618-1922`
- end-of-run replicated mesh profile: `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:2103-2229`

Important v1706-specific differences from the reference:

- current uses `UPstream::init(argc, argv)` instead of the reference `UPstream::initNull()`;
- current uses `OStringStream` / a local `ISpanStream` compatibility wrapper instead of the reference `OCharStream` / `ISpanStream`;
- current default for `replicatedMeshFlatTransfer` is `true`, while the reference source default is `false`. The case controlDict can still override this.

These are compatibility differences, not missing DLB core logic.

### 3. `dsmcCloud` integration

Current `dsmcCloud` has the main replicated mesh hooks:

- member and accessors: `src/lagrangian/dsmc/clouds/dsmcCloud.H:229-230`, `src/lagrangian/dsmc/clouds/dsmcCloudI.H:230-260`
- constructor creates and initializes `dsmcReplicatedMesh` when `replicatedMesh true`: `src/lagrangian/dsmc/clouds/dsmcCloud.C:1331-1335`
- new parcel insertion is owner-filtered:
  - `addNewParcel`: `src/lagrangian/dsmc/clouds/dsmcCloud.C:913-922`
  - `addNewStuckParcel`: `src/lagrangian/dsmc/clouds/dsmcCloud.C:976-985`
- first-step distribution / migration:
  - initial distribution: `src/lagrangian/dsmc/clouds/dsmcCloud.C:1559-1571`
  - periodic owner migration: `src/lagrangian/dsmc/clouds/dsmcCloud.C:1612-1628`
  - manual rebalance steps: `src/lagrangian/dsmc/clouds/dsmcCloud.C:1631-1647`
  - post-collision auto-DLB call: `src/lagrangian/dsmc/clouds/dsmcCloud.C:1697-1710`
- replicated profile reporting is called from `printProfileSummary()`: `src/lagrangian/dsmc/clouds/dsmcCloud.C:2001-2003`

This means the current solver does call `autoRebalance()` once per evolve step when replicated mesh is active. As in the reference, "entered every step" is not the same as "full collective check every step" or "actual repartition every step"; the latter are separately gated by `sarSteps`, `replicatedMeshDLBMinGapSteps`, SAR/imbalance, and forced-step controls.

### 4. Solver-level output path

Current `applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C` has a replicated mesh output branch:

- after `dsmc.evolve()`, if replicated mesh is active, output-time gather is performed by `dsmc.replicatedMeshRef().gatherParcelsToRank0()`;
- only `dsmc.isOutputRank()` calls `runTime.write()`;
- after output, particles migrate back to owners and counts are updated.

Evidence: `applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C:153-194`.

This is directionally aligned with the reference solver path, which also gathers parcels to rank0, writes only from the output rank, then migrates particles back.

### 5. Partial replicated field reduce

Current `dsmcVolFields` includes a generic raw-MPI field reducer:

- helper: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:52-88`
- output-time reduce over many cell and boundary fields: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:2919-2979`
- final profile printing respects `cloud_.isOutputRank()`: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:3995-4005`

This is a partial port. It is not yet the same as the reference field/output implementation; see the missing section below.

## Partially implemented or not yet equivalent

### 1. `dsmcVolFields` replicated sampling/output is not reference-equivalent

The reference implementation uses `replicatedMesh().myCells()` to avoid each replicated rank accumulating all cells. Examples:

- shared sample cache loops over owned cells: reference `dsmcVolFields.C:728-738`, `984-994`, `1320-1327`
- field combine and cell reduction use owned cells: reference `dsmcVolFields.C:2634-2645`, `2889-2930`
- replicated output reduce is followed by non-output-rank skip and rank0-only post-processing: reference `dsmcVolFields.C:3388-3478`
- reference has precomputed owned boundary face lists: reference `dsmcVolFields.H:224-225`, `dsmcVolFields.C:2449-2468`

Current v1706 does not yet have these owned-cell/owned-face paths. It still has all-cell loops in key places:

- shared cache fallback touches all cells: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:609-613`
- shared cache build has no `myCells()` argument/path: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:2518-2535`
- field combine loop falls back to `dsmcNCum_.size()` instead of replicated `myCells()`: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:2565-2576`
- cell reduction still loops `forAll(dsmcNCum_, celli)`: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:2791-2797`
- boundary accumulation loops all faces on sampled patches: `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C:2856-2903`

So the current status is: reduce exists, but reduce-before-write is not enough by itself. The reference avoids duplicate accumulation by restricting what each rank samples before the final sum. This is the highest-priority correctness/performance gap.

### 2. `replicatedMeshDelayedReceive` is not wired in current `dsmcCloud`

The reference `dsmcCloud.C` has an explicit delayed receive path:

- reads `replicatedMeshDelayedReceive`: reference `dsmcCloud.C:3981-3984`
- calls `migrateFinish()` when a previous async migration is pending: reference `dsmcCloud.C:3984-3986`
- uses `migrateBegin()` / `migrateFinish()` around migration steps: reference `dsmcCloud.C:4008-4026`

The current `dsmcReplicatedMesh` class contains `migrateBegin()` / `migrateFinish()`, but current `dsmcCloud.C` only calls the synchronous `migrateParticlesByCellOwner()` path in the evolve loop: `src/lagrangian/dsmc/clouds/dsmcCloud.C:1612-1624`.

This means `replicatedMeshDelayedReceive true` in a controlDict currently gives a false sense of parity with the reference unless the cloud path is wired.

### 3. `Cloud.C` does not have reference-level replicated mesh awareness

The reference `src/lagrangian/basic/Cloud/Cloud.C` contains a generic `hasReplicatedMesh()` detection helper at reference `Cloud.C:147-155`. The current v1706 `src/lagrangian/basic/Cloud/Cloud.C` has OpenMP move-ordered helpers but no equivalent replicated mesh detector in the same area.

Current correctness relies on the outer `dsmcCloud` migration path rather than deeper `Cloud::move()` replicated-mesh awareness. That is acceptable for a first port, but it is not a full reference-equivalent move path.

### 4. Collision diagnostics/profile are less complete than the reference

Current `noTimeCounter.C` computes and stores `nCandidatesPerCell`, which the replicated DLB uses, but it does not have the reference's replicated rank-aware collision diagnostic and accumulated collision counters.

Reference examples:

- `cloud_.accumulateCollisionCounts(collisionCandidates, collisions)`: reference `noTimeCounter.C:492-495`
- replicated raw-MPI rank output when `!Pstream::parRun()`: reference `noTimeCounter.C:520-540`

Current output is more basic, for example `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C:754-770`.

This is mostly a diagnostics/comparison gap, but it matters for DLB analysis because rank-level candidates/collisions are needed to explain whether repartitioning actually improves the critical path.

### 5. Solver/evolve profile summary is not reference-equivalent under raw MPI

Current `dsmcCloud::printProfileSummary()` only does ordinary `Pstream::parRun()` reductions before printing the solver profile: `src/lagrangian/dsmc/clouds/dsmcCloud.C:1913-1951`. It does call `replicatedMesh_->report()` afterward, and that report performs MPI reductions for the replicated mesh summary: `src/lagrangian/dsmc/clouds/dsmcCloud.C:2001-2003`, `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:2103-2229`.

The reference has more explicit replicated raw-MPI final summary/profile handling:

- evolve profile printing is gated by `isOutputRank()`: reference `dsmcCloud.C:5188-5218`
- final replicated rank summary gathers parcels/collisions/candidates with MPI: reference `dsmcCloud.C:5388-5425`
- replicated mesh report is then printed when profiling is enabled: reference `dsmcCloud.C:5428-5431`

Therefore the current `DSMC solver profile summary` numbers in replicated raw-MPI logs should be treated cautiously; the `Replicated mesh profiling summary` is the stronger rank-reduced evidence.

### 6. Lagrangian write-on-proc API from the reference is not ported

The reference has `writeOnProc` plumbing in particle/IO write paths:

- reference `IOPosition.C:57-63`
- reference `particleTemplates.C:170-219`
- reference `dsmcParcelIO.C:282-290`

Current v1706 still uses unconditional field writes in these locations:

- `src/lagrangian/basic/IOPosition/IOPosition.C:51-60`
- `src/lagrangian/basic/particle/particleTemplates.C:160-184`
- `src/lagrangian/dsmc/parcels/dsmcParcelIO.C:406-428`

Because current solver only calls `runTime.write()` on the output rank in replicated mesh mode, this may not be an immediate correctness blocker for the current path. It is still not reference-equivalent and should be checked if output behavior changes or if empty-rank field writes appear.

### 7. Current working case is not configured for replicated mesh DLB

Current case file:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/system/controlDict`

Currently visible contents only include basic run controls through `functions {}` and do not include `replicatedMesh true` or any `replicatedMeshDLB*` controls.

The backup/reference case still has the intended settings:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/mpi8replicatedmesh/system/controlDict:36-59`

Key settings present there:

- `replicatedMesh true`
- `replicatedMeshDelayedReceive true`
- `replicatedMeshNoAlltoall true`
- `replicatedMeshFlatTransfer true`
- `replicatedMeshMigrateInterval 10`
- `replicatedMeshDLBDualConstraint true`
- `replicatedMeshAutoDLB true`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBCheckCollective allgather`

Before any new MPI replicated-mesh run, the current `ourmesh/mpi8replicatedmesh` controlDict must be restored or regenerated intentionally.

## Existing run evidence

### Smoke and 500-step evidence exists, but under old `detail_omp`

Relevant existing logs:

- `doc/worklog/v2506/detail_omp/stage2_runtime_smoke_mpi8replicatedmesh_np2_final_20260605.log`
- `doc/worklog/v2506/detail_omp/stage3_smoke_mpi8replicatedmesh_ownerfilter_10step_20260605.log`
- `doc/worklog/v2506/detail_omp/stage3_formal_ourmesh_mpi8replicatedmesh_500step_20260605.log`
- `doc/worklog/v2506/detail_omp/stage3_formal_ourmesh_mpi8replicatedmesh_ownerfilter_500step_20260605.log`
- summary: `doc/worklog/v2506/detail_omp/stage3_formal_ourmesh_500step_results_20260605.md`

These should be treated as historical evidence. New MPI/DLB work should write new logs and analysis under `doc/worklog/v2506/detail_mpi/`.

### Owner-filtered current run

From `doc/worklog/v2506/detail_omp/stage3_formal_ourmesh_mpi8replicatedmesh_ownerfilter_500step_20260605.log`:

- `real 184.13`
- `solver profile steps = 500`
- `move+collide wall [s] = 105.26`
- `move only [s] = 83.71`
- `buildCellOccupancy [s] = 6.78`
- `collision phase [s] = 14.25`
- `post fields/output [s] = 75.75`
- final particles = `2463897`
- collisions at final output = `34774`
- total energy = `1.239659286`
- migration calls = `57`
- particles per rank = min `206847`, max `411539`, max/min `1.989581671`
- `Phase C auto DLB checks = 500`
- `Phase C auto DLB rebalances = 6`
- `Phase C auto DLB wall max [s] = 1.291806841`
- `Phase C ParMETIS max [s] = 0.397185123`
- `Phase C migration max [s] = 0.913196513`

The summary file records that this fixed the particle explosion:

- unfiltered replicated run final particles = `8089972`
- owner-filtered run final particles = `2463897`
- `mpi8origin` current final particles = `2463565`

So owner filtering is a real correctness improvement.

### Reference/history comparison

Reference/history evidence in `run/.../ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610`:

- external wall seconds = `117.31`
- migration calls = `59`
- local particles final = `296197`
- particles per rank = min `228405`, max `405542`, max/min `1.775539064`
- `Phase C auto DLB checks = 500`
- `Phase C auto DLB rebalances = 8`
- `Phase C auto DLB wall max [s] = 25.67342357`
- `Phase C auto DLB check max [s] = 23.55898278`
- `Phase C ParMETIS max [s] = 0.320407026`
- `Phase C migration max [s] = 1.637910673`

Current owner-filter run is therefore functional but not performance-matched:

- current owner-filter `real 184.13 s`
- reference/history external wall `117.31 s`
- current post fields/output `75.75 s`, which is a major remaining target and aligns with the `dsmcVolFields` gap above.

## What remains to do

Priority order:

1. Fix `dsmcVolFields` replicated owned-cell and owned-boundary-face accumulation to match reference semantics before claiming output correctness.
2. Restore/regenerate `ourmesh/mpi8replicatedmesh/system/controlDict` with the intended replicated mesh DLB settings, and record the exact controlDict in `detail_mpi`.
3. Wire `replicatedMeshDelayedReceive` in current `dsmcCloud` or remove/disable that setting in the v1706 case until it is real.
4. Port or replace the missing reference rank-level collision/profile diagnostics so DLB decisions can be checked by per-rank candidates/collisions and wall time.
5. Re-run `mpi8replicatedmesh` 500 steps under `doc/worklog/v2506/detail_mpi/` and compare against `ourmeshbkp`, checking particle count, collision count, energy, wall time, and Phase C checks/rebalances separately.
6. Only after correctness/output parity, tune DLB cadence and migration options. Keep these three quantities separate in reports:
   - `autoRebalance()` entered per step;
   - collective load check cadence, e.g. `sarSteps`;
   - actual ParMETIS repartition count.
