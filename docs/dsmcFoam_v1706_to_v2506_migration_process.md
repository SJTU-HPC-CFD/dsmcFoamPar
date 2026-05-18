# dsmcFoam+ from OpenFOAM-v1706 to OpenFOAM-v2506

## 1. Purpose

This document records the migration process of `dsmcFoam+` from the
`hyStrath/OpenFOAM-v1706` code base to the `hyStrath_xcx/OpenFOAM-v2506`
environment.

It focuses on:

- migration requirements and scope
- architecture and API differences between `v1706` and `v2506`
- implementation strategy used in this project
- major technical difficulties
- validation and post-migration optimization direction

This is a process document for the actual migration path used in
`hyStrath_xcx`, not a generic OpenFOAM porting guide.


## 2. Background and Requirements

### 2.1 Original baseline

The original solver stack came from:

- `hyStrath`
- host platform: `OpenFOAM-v1706`
- solver: `dsmcFoam+`

The original code was not a small patch over upstream `dsmcFoam`. It included
its own customized lagrangian and DSMC stack, including:

- custom `lagrangian/basic`
- custom `lagrangian/dsmc`
- custom `dsmcCloud`
- custom `dsmcParcel`
- extended collision, reaction, initialization, boundary, and measurement logic

### 2.2 Target baseline

The target platform is:

- `OpenFOAM-v2506`
- working tree: `hyStrath_xcx`

### 2.3 Migration requirements

The migration target was not "full feature parity in one step". The practical
requirements were:

1. Build a usable `v2506`-based DSMC branch first.
2. Restore the minimum solver chain before restoring advanced features.
3. Avoid entangling the first migration phase with AMR and legacy dynamic load
   balancing.
4. Keep the result compatible with the `OpenFOAM-v2506` runtime and build
   system.
5. Preserve the original physical models as much as possible, unless adaptation
   was required by `v2506` API changes.


## 3. Migration Scope and Phasing

The actual migration was performed in stages.

### 3.1 Phase 1: minimum compilable baseline

Target:

- `liblagrangian+`
- `libdsmcFoam+`
- `dsmcFoam+`

At this stage the goal was:

- compile successfully
- run `dsmcFoam+ -help`
- establish a `v2506`-compatible core stack

Features intentionally deferred:

- AMR
- original MPI load balancing path
- some advanced coordinate/time-step models
- large parts of extended boundary, reaction, and measurement logic

### 3.2 Phase 2: minimum runnable solver chain

Target:

- `dsmcInitialise+`
- real case initialization
- time advancement in `dsmcFoam+`
- basic collision and macro field chain

This stage restored enough functionality to run representative non-reacting and
reacting cases.

### 3.3 Phase 3: verification and regression repair

Target:

- repair runtime regressions after compile success
- validate key output fields such as `Ma` and `Tvib`
- compare against the original `hyStrath` implementation

### 3.4 Phase 4: shared-memory performance work

After the `v2506` port became stable, a separate OpenMP optimization line was
introduced in `hyStrath_xcx`. This is post-migration work, but it depends
directly on the migrated architecture and is therefore part of the practical
history of this branch.


## 4. Architecture Differences: v1706 vs v2506

The core migration difficulty did not come from the solver main file alone.
The main problem was the lagrangian infrastructure underneath it.

### 4.1 Lagrangian base layer changed substantially

`hyStrath` in `v1706` relied on its own customized `lagrangian/basic` stack.
That stack was tightly coupled to the older OpenFOAM lagrangian API.

In `v2506`, the lagrangian base layer is different enough that the old
`lagrangian/basic` code cannot simply be copied over.

Typical incompatibilities included:

- particle tracking internals
- `Cloud` behavior and mapping hooks
- template instantiation and registration patterns
- changed assumptions in mesh/particle mapping support

### 4.2 Particle tracking model changed

The `v2506` particle implementation uses a newer tracking model than the
`v1706` branch expected. This affects:

- particle state representation
- face crossing / tracking logic
- mapping and mesh motion related interfaces

As a result, the old `hyStrath` lagrangian base could not be reused unchanged.

### 4.3 Cloud and mapping interfaces changed

Legacy code paths around:

- `Cloud`
- `autoMap()`
- mesh mapping support

were not directly portable. This is one of the reasons AMR and old dynamic load
balancing were explicitly excluded from the first migration stage.

### 4.4 Registration and linkage behavior changed

The port also had to deal with:

- missing or changed type registration patterns
- template linkage failures
- solver/library link-order issues

This was especially visible in early `liblagrangian+` and `libdsmcFoam+`
builds.


## 5. Migration Strategy Used in hyStrath_xcx

The key strategic decision was:

> do not force the entire `v1706` lagrangian base into `v2506`

Instead, the migration followed this rule:

- keep the `v2506` host lagrangian assumptions
- rebuild the DSMC layer on top of them
- restore `hyStrath`-specific logic incrementally

This led to the following implementation strategy.

### 5.1 New working tree

A separate working tree was created:

- `hyStrath_xcx`

This allowed the migration to proceed without destabilizing the original
`hyStrath` source tree.

### 5.2 Keep the new host, adapt the old solver

Instead of transplanting the entire old base layer, the migration reused the
`v2506` host environment and adapted:

- `dsmcCloud`
- `dsmcParcel`
- solver entry logic
- runtime selection registration
- build system wiring

### 5.3 Minimize hard dependencies early

The original `dsmcCloud` and related modules had hard dependencies on many
second-stage features. Early migration therefore required:

- slimming the cloud path
- adding temporary compatibility handling
- deferring or isolating optional functionality

### 5.4 Rebuild, then restore

The practical workflow was:

1. make the libraries build
2. make the solver link
3. make initialization run
4. make simple cases advance in time
5. reintroduce features and repair runtime regressions


## 6. Main Code Areas Involved

The migration concentrated on the following areas.

### 6.1 Core solver and cloud

- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
- `applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C`

### 6.2 Parcel layer

- `src/lagrangian/dsmc/parcels/dsmcParcel.H`
- `src/lagrangian/dsmc/parcels/dsmcParcel.C`
- `src/lagrangian/dsmc/parcels/dsmcParcelIO.C`

### 6.3 Build and registration layer

- `src/lagrangian/basic/Make/*`
- `src/lagrangian/dsmc/Make/*`
- `src/lagrangian/dsmc/dsmcTypes.C`

### 6.4 Initialization and collision chain

- `applications/utilities/preProcessing/dsmc/dsmcInitialise+/dsmcInitialise+.C`
- collision model classes
- partner selection classes

### 6.5 Macro fields and derived statistics

- `src/lagrangian/dsmc/macroscopicProperties/.../dsmcVolFields.C`


## 7. Main Technical Difficulties

### 7.1 Compile success did not imply runtime success

After the first successful build, runtime issues still appeared in:

- initialization
- particle tracking
- macro field generation
- parallel reduction paths

This required a second layer of migration work beyond API adaptation.

### 7.2 `Ma` field regression

One concrete regression during the port was a wrong Mach number field:

- `Ma` became uniformly zero in a non-reacting case

Root cause:

- a protective threshold used `molecularMass > SMALL`
- molecular mass is physically far smaller than `SMALL`
- the entire `Ma` path was skipped

Fix:

- change the relevant check to `VSMALL`

This was a migration-induced bug and had to be corrected in
`dsmcVolFields.C`.

### 7.3 `Tvib` analysis

The vibration-temperature path needed careful comparison against:

- original `hyStrath`
- `SPARTA`
- Bird 2013 definitions

The analysis showed:

- one suspicious multi-mode aggregation issue already existed in old
  `hyStrath`
- the low-excitation roughness in the tested single-mode `N2` case was mainly
  due to physics/statistics and threshold behavior, not the `v2506` port itself

### 7.4 Reactive initialization and runtime restoration

Reacting cases exposed additional fragility in:

- initialization chain
- reaction configuration
- field/property reads
- final runtime behavior under longer evolution

These issues were solved incrementally while validating both non-reacting and
reacting benchmark cases.

### 7.5 Parallel reduction and profiling correctness

During later OpenMP and MPI verification work, additional engineering issues
appeared:

- `Field`-based reductions were not portable as in the older implementation
- profiling output had to be restructured to measure main-loop wall time only
- final MPI profiling output originally contained a reduction deadlock bug,
  later fixed by ensuring all ranks participated in reductions


## 8. Verification Strategy

The migration was validated progressively.

### 8.1 Compile-level validation

- `liblagrangian+`
- `libdsmcFoam+`
- `dsmcFoam+`
- `dsmcInitialise+`

### 8.2 Minimum runtime validation

Representative cases were used to confirm:

- initialization works
- the solver advances
- serial, MPI, and OpenMP variants run correctly

### 8.3 Field-level validation

Specific checks included:

- `Ma` is non-zero and physically reasonable
- `Tvib` behavior is understood and compared to the legacy code and reference
  formulations

### 8.4 Reactive case validation

Reacting cylinder cases were used to compare:

- serial
- MPI
- OpenMP

against:

- total wall time
- profiled main-loop phases
- particle count trends
- collision statistics
- reaction-related trends


## 9. Post-Migration Shared-Memory Optimization

After the `v2506` migration stabilized, the code was extended with OpenMP-based
shared-memory optimization inspired by the load-decoupling paper and
`ParDSMC3D`.

### 9.1 What was adopted

- separate thinking for `Move/Index` and `Collision`
- OpenMP runtime controls
- thread-local RNG streams
- phase profiling inside the main loop
- internal batch-style views to improve shared-memory processing

### 9.2 What was not fully adopted

The current code still does not match `ParDSMC3D` completely because:

- the canonical particle container is still not replaced by a native
  cell-contiguous storage design
- `partition` scheduling does not always outperform `dynamic`
- some benefits of the paper depend on stronger data locality than the current
  `hyStrath_xcx` structure provides

### 9.3 Practical outcome

Despite that gap, the migrated `hyStrath_xcx` branch achieved:

- stable OpenMP support
- strong speedup on representative reacting cases
- in some longer reactive runs, OpenMP outperforming MPI because MPI collision
  load imbalance became severe


## 10. Lessons Learned

### 10.1 Port the solver stack, not just the solver file

For a code base like `dsmcFoam+`, migration cannot be done by editing only the
solver main program. The real work is in:

- cloud infrastructure
- parcel behavior
- initialization
- registration
- macro fields
- runtime control paths

### 10.2 Reuse the new host architecture where possible

Trying to force the full `v1706` lagrangian base into `v2506` would have been
more fragile and more expensive than rebuilding a compatible DSMC layer on top
of the new host assumptions.

### 10.3 Stage the work aggressively

Separating the work into:

- minimum compile
- minimum runtime
- correctness repair
- performance work

was necessary to keep the migration tractable.

### 10.4 Performance tuning is a separate phase

The OpenMP work was valuable, but it only made sense after the solver had
already become correct and stable on `v2506`.


## 11. Current Status

At the current stage, `hyStrath_xcx` represents:

- a successful port of `dsmcFoam+` from `OpenFOAM-v1706` to `OpenFOAM-v2506`
- with serial, MPI, and OpenMP execution paths available
- with major runtime regressions already repaired
- with a growing shared-memory optimization layer on top of the migrated base

The major remaining architectural gap relative to `ParDSMC3D` is still the data
layout:

- `hyStrath_xcx` has moved toward shared batch views
- but the canonical storage model is still tied to the OpenFOAM-style cloud
  container model

That is the main reason why some paper-style scheduling strategies are not yet
universally superior in the current code base.


## 12. Recommended Follow-up Documents

This process document should be complemented by:

1. a build-and-run guide for the current `hyStrath_xcx` branch
2. a validation report for non-reacting and reacting regression cases
3. a separate note for the OpenMP load-decoupling implementation in
   `hyStrath_xcx`
4. a file-level migration appendix generated from `git diff --name-only`

