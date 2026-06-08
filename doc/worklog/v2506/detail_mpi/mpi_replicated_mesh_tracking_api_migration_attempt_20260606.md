# MPI replicated-mesh tracking API migration attempt - 2026-06-06

## Scope

Case:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Test policy:

- pure MPI only
- `useOpenMP false`
- `openmpThreads 1`
- 8 MPI ranks
- final comparison target remains `ourmeshbkp/mpi8replicatedmesh`, not the older OpenMP2 113s run

## Compatibility audit

The reference source uses the newer lower-level tracking API:

- `particle::trackToAndHitFace()`
- `particle::deviationFromMeshCentre()`
- barycentric particle coordinates and displacement/fraction-based tracking

The current OF-v1706 tree still uses the older particle core:

- particle position is stored as `position_`
- `dsmcParcel::move()` calls
  `trackToFace(position() + dt*Utracking, td, true)`
- `dsmcParcel::move()` owns the step-fraction update through `tEnd/dtCell`

A direct reference-style migration is therefore not a small `dsmcParcel.C`
change.  It would expand into the base particle storage model, constructors,
particle IO, parallel transfer, replicated-mesh flat transfer packing, cyclic
transform handling, and possibly newer support types such as barycentric
coordinates and vector/tensor transform helpers.

For this round, the tested minimum migration was a compatibility layer on top
of the existing OF-v1706 tracking core, not a wholesale barycentric particle
port.

## Source attempt

Temporary changes:

- added `particle::trackToAndHitFace(displacement, fraction, td, true)` as a
  wrapper around the old `trackToFace(position_ + displacement, td, true)` path;
- added `particle::deviationFromMeshCentre()`;
- changed the uniform-deltaT branch of `dsmcParcel::move()` to use the
  reference-style displacement/fraction call form;
- kept the non-uniform-deltaT branch on the old `tEnd/dtCell` path.

Build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_tracking_api_compat_build_20260606.log`

Build result:

- completed successfully;
- only the known OpenFOAM-v1706 template-instantiation warnings were present.

## 10-step smoke

Temporary case change:

- `endTime 1.e-06`
- `useOpenMP false`
- `openmpThreads 1`

Control backup:

- `doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_tracking_api_compat_smoke_20260606`

Log:

- `log.codex_mpi8_replicatedmesh_mpionly_tracking_api_compat_smoke_10step_20260606`

Result:

| metric | value |
|---|---:|
| real | 10.87 s |
| Total Iterations | 10 |
| OpenMP enabled | 0 |
| OpenMP max threads | 1 |
| move only | 1.218105333 s |
| buildCellOccupancy | 0.110084485 s |
| collision phase | 0.190928628 s |
| post fields/output | 0.85127116 s |
| full evolve wall | 2.498202913 s |
| particles | 2065849 |
| collisions | 2812 |
| collision candidates | 2910 |
| total energy | 1.150617816 |
| DLB checks | 10 |
| DLB rebalances | 0 |

No `Fatal`, `Segmentation`, `Floating`, `NaN`, or `nan` entries were found.

The 500-step controlDict was restored after the smoke.  Restored hash:

```text
31a78a69f94450ed6338861d2cd24cff5b0ace39d4536d2356f2cd81b47a9b87
```

## 500-step formal validation

Log:

- `log.codex_mpi8_replicatedmesh_mpionly_tracking_api_compat_500step_20260606`

Result:

| metric | pre-candidate current | compatibility API candidate | delta |
|---|---:|---:|---:|
| real | 135.30 s | 140.79 s | +5.49 s, 4.06% slower |
| move only | 72.59819434 s | 82.92498684 s | +10.32679250 s, 14.22% slower |
| buildCellOccupancy | 7.432408543 s | 8.25212966 s | +0.819721117 s |
| collision phase | 12.97955973 s | 15.08526616 s | +2.10570643 s |
| post fields/output | 38.10104302 s | 37.82146892 s | -0.27957410 s |
| full evolve wall | 125.5070652 s | 131.3125433 s | +5.8054781 s |
| DLB checks | 500 | 500 | unchanged |
| DLB rebalances | 5 | 4 | -1 |

Final physical summary:

| metric | pre-candidate current | compatibility API candidate |
|---|---:|---:|
| particles | 2463986 | 2463850 |
| stuck particles | 0 | 0 |
| collisions | 35211 | 34182 |
| collision candidates | 65708 | 60967 |
| total energy | 1.240619931 | 1.241065816 |

Against the preserved `ourmeshbkp` baseline:

| metric | ourmeshbkp | compatibility API candidate | gap |
|---|---:|---:|---:|
| external/real wall | 117.31 s | 140.79 s | +23.48 s, 20.02% slower |
| move only | 47.1652636 s | 82.92498684 s | +35.75972324 s |
| buildCellOccupancy | 10.7084157 s | 8.25212966 s | -2.45628604 s |
| collision phase | 10.89824478 s | 15.08526616 s | +4.18702138 s |
| post fields/output | 30.35563996 s | 37.82146892 s | +7.46582896 s |
| full evolve wall | 114.9207628 s | 131.3125433 s | +16.3917805 s |
| DLB rebalances | 8 | 4 | -4 |

## Decision

Rejected and reverted.

The compatibility API layer compiled and passed the 10-step smoke, but the
formal 500-step pure-MPI run made move materially worse.  This shows that merely
changing the outer call shape to `trackToAndHitFace` on top of the old
OF-v1706 Cartesian tracking core does not recover the reference performance.

The result also narrows the interpretation:

- the performance gap is not caused by the shape of the `dsmcParcel::move()`
  loop alone;
- a useful migration must change the lower-level tracking algorithm/data path,
  or find an equivalent fast path inside the existing `particle::trackToFace`
  implementation;
- a full barycentric particle port remains high-risk because it crosses IO and
  transfer semantics, not only move tracking.

## Revert state

The temporary source changes were manually reverted without using `git reset`,
`git checkout`, stash, or staging operations.

Revert build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_tracking_api_compat_revert_build_20260606.log`

Post-revert checks:

- `git diff -- src/lagrangian/basic/particle/particle.H src/lagrangian/basic/particle/particle.C src/lagrangian/basic/particle/particleTemplates.C src/lagrangian/dsmc/parcels/dsmcParcel.C` is empty;
- `rg "trackToAndHitFace|deviationFromMeshCentre|handleTrackedFaceTransition" src/lagrangian/basic/particle src/lagrangian/dsmc/parcels/dsmcParcel.C` has no matches;
- `system/controlDict` remains at the 500-step pure-MPI hash
  `31a78a69f94450ed6338861d2cd24cff5b0ace39d4536d2356f2cd81b47a9b87`.

