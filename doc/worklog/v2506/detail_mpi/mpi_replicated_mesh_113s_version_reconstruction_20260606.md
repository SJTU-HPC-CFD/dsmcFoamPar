# MPI replicated-mesh 113s version reconstruction - 2026-06-06

## Purpose

This note reconstructs the code, case configuration, and build state behind the
best retained MPI replicated-mesh result:

- case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- log: `log.codex_mpi8_replicatedmesh_dlb_openmp2_500step_20260606`
- result: `real 113.48 s`

The reconstruction uses the Codex session history, current worklogs, and build
logs.  There is no separate source backup, git commit, or preserved old shared
library for the 113s binary.

## Reconstructed source manifest

At the 113s point, the retained source diff relative to `HEAD` was limited to
the MPI replicated-mesh/output path:

- `applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
- `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C`
- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`
- `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H`

The session `git status` immediately before the OpenMP2 500-step result showed
no active `Cloud.C`, `dsmcCloud.H`, `dsmcParcel.C`, or
`particleTemplates.C` diff.  Those files were part of the already-built OMP
baseline, not part of the MPI 113s delta.

Current re-check after the rollback attempts gives the same source manifest:

```text
git diff --name-only -- <113s-related source files>
applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C
src/lagrangian/dsmc/clouds/dsmcCloud.C
src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.H

git diff --stat -- <same files plus Cloud.C/dsmcCloud.H/dsmcParcel.C>
5 files changed, 299 insertions(+), 45 deletions(-)
```

`Cloud.C`, `dsmcCloud.H`, `dsmcParcel.C`, and `particleTemplates.C` still have
no source diff.  `git diff --check` on the five retained source files is clean.

The MPI delta consisted of:

- `dsmcVolFields`: owner-cell sampling/reduction for replicated mesh, owned
  boundary-face accumulation, and output-rank-only field computation after the
  collective reduce;
- `dsmcCloud`: `replicatedMeshDelayedReceive` wiring, async migration finish
  before manual/auto DLB and output, replicated raw-MPI profile reduction, and
  output-rank gating;
- `noTimeCounter`: replicated raw-MPI rank collision/candidate diagnostics and
  output-rank gating for ordinary logs;
- `dsmcFoam+.C`: finish pending async replicated migration before solver output.

The later `dsmcReplicatedMesh.C` candidate-weight change was not part of the
113s version.  It was tested after the 113s run and reverted.

## Reconstructed case configuration

The 113s run used the following full 500-step configuration:

```text
useOpenMP true;
openmpThreads 2;
openmpMoveSchedule static;
openmpMoveChunk 64;
openmpCollisionSchedule dynamic;
openmpCollisionChunk 8;

profileSummary true;
profileDetail false;
collisionFastRng true;

replicatedMesh true;
replicatedMeshDelayedReceive true;
replicatedMeshNoAlltoall true;
replicatedMeshFlatTransfer true;
replicatedMeshDecompMethod metis;
replicatedMeshMigrateInterval 10;
replicatedMeshDLBDualConstraint true;
replicatedMeshAutoDLB true;
replicatedMeshDLBSteps 50;
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBCheckCollective allgather;
replicatedMeshDLBSkipFixedKPostDiag false;
replicatedMeshDLBImbalanceThreshold 1.5;
replicatedMeshDLBFixedK 1;
replicatedMeshDLBInitialK 1024;
replicatedMeshDLBItr 1000;
replicatedMeshDLBUbvec 1.05;
replicatedMeshDLBUbvec1 1.5;
replicatedMeshDLBProfile false;
replicatedMeshDLBProfileSteps 5;
replicatedMeshDLBInitialAlpha 0.8;
replicatedMeshDLBAdaptiveKMode 0;
replicatedMeshDLBVsizeExp 0;
```

The current `controlDict` hash for this reconstructed configuration is:

```text
a2f2e53d3b3d94765420f33b5748c6330c8af0896ae1474eef339739a6cc4d17
```

## Build-state evidence

The source/config manifest above is recoverable.  The exact 113s binary state is
not fully recoverable from text logs alone.

The relevant build log is:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_replicated_mesh_dlb_build_20260606.log
```

That 05:52 build compiled only:

- `clouds/dsmcCloud.C`
- `collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C`
- `macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`
- `applications/.../dsmcFoam+.C`

It then linked `libdsmcFoam+.so` and `dsmcFoam+`.  The 113s run was launched
after changing only `controlDict` to `useOpenMP true; openmpThreads 2;`, without
a new rebuild.

The later 10:02 rollback build touched `Cloud.C` and `dsmcCloud.H` timestamps
while restoring them from `HEAD`.  That caused a broad rebuild of `dsmc`
objects, including `dsmcParcel.o`, `dsmcCloud.o`, `dsmcReplicatedMesh.o`,
`noTimeCounter.o`, `dsmcVolFields.o`, many boundary/controller objects, and a
fresh `libdsmcFoam+.so`.  After that rebuild, the formal run slowed to:

- `log.codex_mpi8_replicatedmesh_dlb_rollback113_500step_20260606`
- `real 127.81 s`
- `move only 78.42538533 s`
- `post fields/output 42.94822251 s`

The current artifact timestamps confirm that the old 05:52 shared-library state
has been overwritten:

```text
2026-06-06 05:52:12 platforms/linux64IccDPInt32Opt/bin/dsmcFoam+
2026-06-06 10:02:35 platforms/linux64IccDPInt32Opt/lib/libdsmcFoam+.so
2026-06-06 09:59:15 src/lagrangian/dsmc/Make/.../parcels/dsmcParcel.o
2026-06-06 09:59:20 src/lagrangian/dsmc/Make/.../clouds/dsmcCloud.o
2026-06-06 10:00:33 src/lagrangian/dsmc/Make/.../noTimeCounter.o
2026-06-06 10:02:15 src/lagrangian/dsmc/Make/.../dsmcVolFields.o
2026-06-06 05:52:12 applications/.../dsmcFoam+.o
```

The `dsmcParcel.C` / `particleTemplates.C` trail is therefore an object-state
trail, not a source-diff trail.  `dsmcParcel.o` was recompiled during the later
broad rebuild because the templated move path includes `Cloud.C` and cloud
headers.  `particleTemplates.C` appears through template instantiation notes in
that rebuild log; there is no standalone retained source edit for it.

## Current 500-step retest

The current artifact state was retested without rebuilding:

```text
log.codex_mpi8_replicatedmesh_current_retest_500step_20260606
controlDict sha256: a2f2e53d3b3d94765420f33b5748c6330c8af0896ae1474eef339739a6cc4d17
libdsmcFoam+.so sha256: 08dc33892973f46d0f21fb1f700a1019bfa49dbcefafad7e1db6317d0c3e5667
dsmcFoam+ sha256: 09d82cc54e60767a34e41cbe34dcfa185f2e1b348f0ee2144e25425860496ab1
```

Result:

```text
real 129.02 s
Total Iterations 500
final particles 2463789
final collisions 35332
final candidates 64389
final total energy 1.24050557
move+collide wall 94.19475275 s
move only 80.89923011 s
buildCellOccupancy 7.081801184 s
collision phase 5.660002116 s
post fields/output 43.42733056 s
full evolve wall 121.9367765 s
DLB checks 500
DLB rebalances 4
```

This retest confirms the current binary/object state is in the same performance
class as the rollback result (`127.81 s`), not the old 113s result.  The main
regression versus 113s remains move/post:

```text
113s log:     real 113.48, move 69.34624257, post 33.6372609
current log: real 129.02, move 80.89923011, post 43.42733056
```

## MPI-only 8-rank retest

The previous current retest used the restored `controlDict` with
`useOpenMP true; openmpThreads 2;`.  A separate MPI-only retest was run by
temporarily switching only:

```text
useOpenMP false;
openmpThreads 1;
```

The original `controlDict` was backed up as
`doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_mpionly_retest_20260606`
and restored after the run.  The restored file hash again matches:

```text
a2f2e53d3b3d94765420f33b5748c6330c8af0896ae1474eef339739a6cc4d17
```

Run log:

```text
log.codex_mpi8_replicatedmesh_mpionly_retest_500step_20260606
```

Result:

```text
real 138.51 s
Total Iterations 500
OpenMP enabled 0
OpenMP max threads 1
final particles 2463829
final collisions 35205
final candidates 65494
final total energy 1.240333986
move+collide wall 98.87990121 s
move only 77.41416822 s
buildCellOccupancy 7.276991174 s
collision phase 13.81711219 s
post fields/output 35.7837221 s
full evolve wall 128.8437596 s
DLB checks 500
DLB rebalances 4
```

Swap was checked before, during, and after this MPI-only retest and remained at
`0B` used.  This confirms the MPI-only current result is slower than the old
OpenMP2 113s result, but the cost distribution is different: MPI-only improves
post/output versus the current OMP2 retest, while collision becomes much slower.

## Manual post-113 source rollback

The rejected candidate-weight source attempt was documented as changing the
second ParMETIS constraint from the retained `N*(N-1)` collision proxy to a
measured candidate-count based weight.  Current code had already restored the
actual vertex weight to `N*(N-1)`, but still retained post-candidate dead code:

```text
nCandidatesPerCell / totalCandidates / collScale / local collDivisor
```

Those leftovers were removed manually from:

```text
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C
```

The resulting `buildPartitionParMetisAdaptive()` weight path is again only:

```text
constraint 0: N^alpha
constraint 1: N*(N-1)
```

Build:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_revert_post113_candidate_residual_build_20260606.log
```

The build recompiled `replicatedMesh/dsmcReplicatedMesh.C` and relinked
`libdsmcFoam+.so` / `dsmcFoam+`.

Validation log:

```text
log.codex_mpi8_replicatedmesh_post113_revert_500step_20260606
```

Result:

```text
real 115.88 s
Total Iterations 500
OpenMP enabled 1
OpenMP max threads 2
move+collide wall 80.69254188 s
move only 68.72489113 s
buildCellOccupancy 6.391809092 s
collision phase 5.408795021 s
post fields/output 39.84684517 s
full evolve wall 106.3524297 s
final particles 2463742
final collisions 34626
final candidates 64498
final total energy 1.2408346
DLB checks 500
DLB rebalances 4
particles max/min 1.526559674
```

This manual rollback restores the move path to the 113s performance class:

```text
113s log:          real 113.48, move 69.34624257, post 33.6372609
current bad retest: real 129.02, move 80.89923011, post 43.42733056
post-rollback:     real 115.88, move 68.72489113, post 39.84684517
```

The remaining gap to 113.48s is no longer the move kernel.  It is mostly in
post fields/output and run-to-run DLB/output timing.

## MPI-only retest after manual rollback

The user narrowed this round to MPI replicated-mesh DLB only and explicitly
requested that following tests do not enable OpenMP.  The case controlDict was
therefore left as:

```text
useOpenMP false;
openmpThreads 1;
```

The first sandboxed launch failed before solver startup because Intel MPI/Hydra
could not open a local listening socket:

```text
HYD_sock_listen_on_port: cannot open socket (Operation not permitted)
```

The same command was rerun outside the sandbox:

```text
source doc/scripts/env.sh
cd run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh
OMP_NUM_THREADS=1 /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_post113_revert_500step_20260606 2>&1
```

Validation log:

```text
log.codex_mpi8_replicatedmesh_mpionly_post113_revert_500step_20260606
```

Result:

```text
real 135.30 s
Total Iterations 500
OpenMP enabled 0
OpenMP max threads 1
move+collide wall 93.99630213 s
move only 72.59819434 s
buildCellOccupancy 7.432408543 s
collision phase 12.97955973 s
post fields/output 38.10104302 s
full evolve wall 125.5070652 s
final particles 2463986
final collisions 35211
final candidates 65708
final total energy 1.240619931
DLB checks 500
DLB rebalances 5
particles max/min 2.431223683
```

Swap was checked before, during, and after the run and remained at `0B` used.

Comparison against the previous MPI-only retest:

```text
pre-rollback MPI-only:  real 138.51, move 77.41416822, collision 13.81711219, post 35.7837221, rebalances 4
post-rollback MPI-only: real 135.30, move 72.59819434, collision 12.97955973, post 38.10104302, rebalances 5
```

For this MPI-only scope, the manual rollback improves end-to-end time by about
3.21 s versus the earlier MPI-only retest, mostly through the move path and a
smaller collision reduction.  This does not change the earlier conclusion that
the historical 113.48 s result was an OpenMP2 run, not a pure MPI run.

## Current performance comparison against ourmeshbkp

This comparison intentionally excludes the historical 113.48 s OpenMP2 run.
The only reference baseline considered here is the preserved
`ourmeshbkp/mpi8replicatedmesh` 500-step log.

Current version:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_mpionly_post113_revert_500step_20260606
8 MPI, OpenMP disabled, 500 steps
real 135.30 s
```

Reference baseline:

```text
run/.../ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610
500 steps
External wall seconds 117.31
```

Wall-time comparison:

| Metric | Current | ourmeshbkp | Delta | Relative |
|---|---:|---:|---:|---:|
| external wall | 135.30 s | 117.31 s | +17.99 s | +15.34% slower |
| full evolve wall | 125.5070652 s | 114.9207628 s | +10.5863024 s | +9.21% slower |

Comparable profile fields:

| Metric | Current | ourmeshbkp | Delta | Relative |
|---|---:|---:|---:|---:|
| move only | 72.59819434 s | 47.1652636 s | +25.43293074 s | +53.92% |
| buildCellOccupancy | 7.432408543 s | 10.7084157 s | -3.27600716 s | -30.59% |
| collision phase | 12.97955973 s | 10.89824478 s | +2.08131495 s | +19.10% |
| post fields/output | 38.10104302 s | 30.35563996 s | +7.74540306 s | +25.52% |

DLB and balance comparison:

| Metric | Current | ourmeshbkp | Delta |
|---|---:|---:|---:|
| DLB checks | 500 | 500 | 0 |
| DLB rebalances | 5 | 8 | -3 |
| migration calls | 56 | 59 | -3 |
| migration wall | 1.72624754 s | 4.406303137 s | -2.680055597 s |
| Phase C auto DLB wall max | 1.219253054 s | 25.67342357 s | -24.454170516 s |
| Phase C check max | 0.022978711 s | 23.55898278 s | -23.536004069 s |
| Phase C ParMETIS max | 0.342560905 s | 0.320407026 s | +0.022153879 s |
| Phase C migration max | 0.905531486 s | 1.637910673 s | -0.732379187 s |
| particles max/min | 2.431223683 | 1.775539064 | +0.655684619 |
| rank wall max/min | 1.077137871 | 1.080575306 | -0.003437435 |

Correctness-scale comparison:

| Metric | Current | ourmeshbkp | Delta |
|---|---:|---:|---:|
| final particles | 2463986 | 2466342 | -2356 (-0.096%) |
| last-step collisions | 35211 | 35158 | +53 (+0.151%) |
| last-step candidates | 65708 | 70116 | -4408 (-6.287%) |

Against the preserved `ourmeshbkp` baseline, the current run is 17.99 s slower
end-to-end.  The DLB overhead itself is much lower in the current run, but this
is outweighed by a much slower move phase and slower post fields/output.  The
rank-wall imbalance is essentially unchanged, so the particle max/min imbalance
increase does not directly translate into worse rank wall spread in this run.

Therefore the best current explanation is:

1. The source/config state of the 113s version can be reconstructed.
2. The exact 113s performance likely depended on the old 05:52 incremental
   object/library state.
3. That object state was overwritten by the 10:02 broad rebuild.
4. Without the old `.o`/`.so` files, the exact binary cannot be recreated only
   from chat/worklog text.

## Practical recovery boundary

The recoverable state is the source/config manifest above.  The unrecoverable
state is the old linked binary/object mixture.

For future performance checkpoints, keep both:

- the source/config hash list;
- a binary artifact snapshot of `platforms/linux64IccDPInt32Opt/lib/libdsmcFoam+.so`,
  `platforms/linux64IccDPInt32Opt/bin/dsmcFoam+`, and the relevant
  `src/lagrangian/*/Make/linux64IccDPInt32Opt/*.o` files before any diagnostic
  rebuild that touches shared headers.
