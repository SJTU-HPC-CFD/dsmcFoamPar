# MPI replicated-mesh move candidate rejections - 2026-06-06

## Kept baseline

Current kept formal baseline remains:

- case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- log: `log.codex_mpi8_replicatedmesh_dlb_openmp2_500step_20260606`
- configuration: `useOpenMP true`, `openmpThreads 2`, `openmpMoveSchedule static`, `openmpMoveChunk 64`
- formal result: `real = 113.48 s`, `move only = 69.34624257 s`, `full evolve wall = 104.73 s`, `post fields/output = 33.6372609 s`
- replicated summary: `rebalances = 4`, `migration calls = 55`, `particles max/min = 1.960998784`

The source and case configuration were restored after all rejected tests below.

## Rejected candidate: DLB alpha/gap/threshold retune

Temporary change:

- `endTime 2.e-05`
- `replicatedMeshDLBInitialAlpha 1.0`
- `replicatedMeshDLBMinGapSteps 25`
- `replicatedMeshDLBImbalanceThreshold 1.25`

Log:

- `log.codex_mpi8_replicatedmesh_dlb_move_alpha1_gap25_smoke_200step_20260606`

Result:

| metric | value |
|---|---:|
| real | 51.57 s |
| move only | 27.06028949 s |
| full evolve wall | 42.68083678 s |
| particles max/min | 3.475835014 |
| rebalances | 4 |

Decision: rejected.  The original 500-step kept log reached iteration 200 at about `ClockTime = 48 s`; this candidate reached `real = 51.57 s` and badly worsened particle balance.

## Rejected candidate: DLB gap/threshold retune with kept alpha

Temporary change:

- `endTime 2.e-05`
- kept `replicatedMeshDLBInitialAlpha 0.8`
- `replicatedMeshDLBMinGapSteps 25`
- `replicatedMeshDLBImbalanceThreshold 1.25`

Log:

- `log.codex_mpi8_replicatedmesh_dlb_alpha08_gap25_thr125_smoke_200step_20260606`

Result:

| metric | value |
|---|---:|
| real | 55.25 s |
| move only | 29.46732052 s |
| full evolve wall | 46.5771577 s |
| particles max/min | 1.589474703 |
| rebalances | 3 |

Decision: rejected.  Even though particle max/min improved versus the kept 500-step final summary, the 200-step wall and move time were much worse than the kept run's iteration-200 clock (`ClockTime = 48 s`).

## Rejected candidate: candidate-count collision DLB weight

Source attempt:

- changed the ParMETIS dual-constraint collision weight from `N*(N-1)` to measured `nCandidatesPerCell * collScale / collDivisor`, with the old `N*(N-1)` proxy as fallback.

Build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_move_candidate_candidateweight_build_20260606.log`

Smoke log:

- `log.codex_mpi8_replicatedmesh_dlb_candidateweight_smoke_200step_20260606`

Formal log:

- `log.codex_mpi8_replicatedmesh_dlb_candidateweight_500step_20260606`

Formal result:

| metric | kept baseline | candidate |
|---|---:|---:|
| real | 113.48 s | 123.50 s |
| move only | 69.34624257 s | 75.63816822 s |
| full evolve wall | 104.73 s | 117.0024191 s |
| post fields/output | 33.6372609 s | 42.0309142 s |
| particles max/min | 1.960998784 | 4.416719647 |
| rebalances | 4 | 4 |

Decision: rejected and reverted.  The measured-candidate weight made the final owner partition much worse for particle balance and increased move time.

Revert build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_move_candidate_candidateweight_revert_build_20260606.log`

## Rejected candidate: candidate-count weight plus `FixedK 128`

Temporary change on top of the source candidate:

- `endTime 2.e-05`
- `replicatedMeshDLBFixedK 128`

Log:

- `log.codex_mpi8_replicatedmesh_dlb_candidateweight_fixedk128_smoke_200step_20260606`

Result:

| metric | value |
|---|---:|
| real | 51.70 s |
| move only | 27.493807 s |
| full evolve wall | 43.50373681 s |
| particles max/min | 3.017444138 |
| rebalances | 2 |

Decision: rejected.  It was slower than both the source-only 200-step candidate and the kept baseline's iteration-200 clock.

## Rejected candidate: larger/smaller OpenMP move chunks

Temporary 10-step runs, with only `endTime 1.e-06` and `openmpMoveChunk` changed:

| chunk | log | real | move only | full evolve wall |
|---:|---|---:|---:|---:|
| 64 | `log.codex_mpi8_replicatedmesh_dlb_openmp2_smoke_10step_20260606` | 10.70 s | 1.288977197 s | 2.319679209 s |
| 4096 | `log.codex_mpi8_replicatedmesh_dlb_openmp2_chunk4096_smoke_10step_20260606` | 12.23 s | 1.456535383 s | 2.709025113 s |
| 128 | `log.codex_mpi8_replicatedmesh_dlb_openmp2_chunk128_smoke_10step_20260606` | 11.89 s | 1.610710126 s | 2.853151922 s |
| 32 | `log.codex_mpi8_replicatedmesh_dlb_openmp2_chunk32_smoke_10step_20260606` | 12.12 s | 1.556271252 s | 2.760580777 s |

Decision: keep `openmpMoveChunk 64`.  The tested larger and smaller chunks all degraded move and end-to-end smoke time.

## Diagnostic: OpenMP move stage profile

Source diagnostic retained:

- Added default-off `openmpMoveProfile` / `openmpMoveProfileInterval` controls.
- Added output-rank profile printing in both OpenMP move branches in `Cloud.C`.
- Added `dsmcCloud` accessors so the templated `Cloud.C` path reads the controls reliably.

Build logs:

- `doc/worklog/v2506/detail_mpi/stage_mpi_move_profile_build_20260606.log`
- `doc/worklog/v2506/detail_mpi/stage_mpi_move_profile_accessor_build_20260606.log`
- `doc/worklog/v2506/detail_mpi/stage_mpi_move_profile_pstream_build_20260606.log`

Temporary profile logs:

- `log.codex_mpi8_replicatedmesh_dlb_moveprofile_smoke_10step_20260606`
- `log.codex_mpi8_replicatedmesh_dlb_moveprofile_interval1_smoke_10step_20260606`

Key observation:

- The emitted line is `OpenMP move profile (Pstream)`, so this MPI replicated-mesh run enters the `Pstream::parRun()` OpenMP branch in `Cloud.C`.
- The no-candidate interval-1 profile run reports, on output rank over 10 steps:

| metric | value |
|---|---:|
| profile parcels | 2,524,883 |
| deleted | 0 |
| transfers | 0 |
| prepare wall | 0.099198002 s |
| kernel wall | 0.768926865 s |
| commit wall | 0.218797160 s |
| profile total wall | 1.086922027 s |
| solver move-only max | 1.389663182 s |
| real | 10.75 s |

Interpretation: move is still kernel dominated.  The output-rank local profile is lower than the solver's reduced move max, so it should be used for stage proportions and branch identification, not as the final performance metric.

Profile-off sanity log:

- `log.codex_mpi8_replicatedmesh_dlb_profileoff_smoke_10step_20260606`
- No `OpenMP move profile` lines are emitted when the controls are absent.
- The smoke completed normally with `Total Iterations = 10`.

## Rejected candidate: Pstream no-change fast path

Source attempt:

- In the Pstream OpenMP move branch, added a reduction over changed particles.
- If no rank had deletes or processor transfers, skipped the survivor/transfer/delete compaction path and stored the original particle order directly.

Build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_move_pstream_nochange_fastpath_build_20260606.log`

Smoke log:

- `log.codex_mpi8_replicatedmesh_dlb_pstream_nochange_fastpath_interval1_smoke_10step_20260606`

Result versus the no-candidate interval-1 profile run:

| metric | no candidate | no-change fast path |
|---|---:|---:|
| real | 10.75 s | 12.35 s |
| solver move-only max | 1.389663182 s | 1.498996041 s |
| profile kernel wall | 0.768926865 s | 0.828943878 s |
| profile commit wall | 0.218797160 s | 0.294153668 s |
| profile total wall | 1.086922027 s | 1.207208690 s |

Decision: rejected and reverted.  The added changed-particle reduction and branch did not pay for itself; both move and end-to-end smoke time degraded.

Revert build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_move_pstream_nochange_fastpath_revert_build_20260606.log`

## Final state

- `system/controlDict` restored to the pre-candidate hash:
  `a2f2e53d3b3d94765420f33b5748c6330c8af0896ae1474eef339739a6cc4d17`
- `src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C` restored; no residual diff from the rejected candidate remains.
- The Pstream no-change fast path was reverted.
- The default-off move profile instrumentation remains available for diagnostics.
- The best validated MPI replicated-mesh result remains `log.codex_mpi8_replicatedmesh_dlb_openmp2_500step_20260606`.

## Current profile-only 500-step validation

Formal run after reverting the Pstream no-change fast path and restoring the
500-step case configuration:

- Log: `log.codex_mpi8_replicatedmesh_dlb_profileonly_500step_20260606`
- Temporary controls: none; `openmpMoveProfile` is absent.
- Result: completed normally with `Total Iterations = 500`.

| metric | previous best | current profile-only |
|---|---:|---:|
| real | 113.48 s | 128.53 s |
| move only | 69.34624257 s | 78.51006665 s |
| buildCellOccupancy | 5.975890432 s | 6.727392902 s |
| collision phase | 5.523175253 s | 6.067793177 s |
| post fields/output | 33.63726090 s | 41.10214107 s |
| full evolve wall | 104.7277518 s | 121.9378848 s |
| particles max/min | 1.960998784 | 2.165081371 |
| DLB rebalances | 4 | 4 |

Final physical summary at `Time = 5e-05`:

- DSMC particles: `2463754`
- stuck particles: `0`
- total energy: `1.240411821`

Decision: current profile-only state does not improve the validated best.  Treat
`log.codex_mpi8_replicatedmesh_dlb_openmp2_500step_20260606` as the retained
performance baseline unless another candidate beats it in a fresh 500-step run.

## Pure-MPI move optimization attempts after scope narrowing

Scope update:

- Only MPI replicated-mesh DLB is considered in this round.
- OpenMP is disabled for all follow-up tests:
  `useOpenMP false`, `openmpThreads 1`, `OMP_NUM_THREADS=1`.
- The comparison target for this round is the preserved
  `ourmeshbkp/mpi8replicatedmesh` baseline, not the historical OpenMP2 113s run.

Pre-candidate current reference:

- Log: `log.codex_mpi8_replicatedmesh_mpionly_post113_revert_500step_20260606`

| metric | pre-candidate current |
|---|---:|
| real | 135.30 s |
| move only | 72.59819434 s |
| buildCellOccupancy | 7.432408543 s |
| collision phase | 12.97955973 s |
| post fields/output | 38.10104302 s |
| full evolve wall | 125.5070652 s |
| DLB rebalances | 5 |

ourmeshbkp reference:

- Log: `run/.../ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610`

| metric | ourmeshbkp |
|---|---:|
| external wall | 117.31 s |
| move only | 47.1652636 s |
| buildCellOccupancy | 10.7084157 s |
| collision phase | 10.89824478 s |
| post fields/output | 30.35563996 s |
| full evolve wall | 114.9207628 s |
| DLB rebalances | 8 |

### Rejected candidate: serial move pointer snapshot plus Cartesian-tracking cache

Source attempt:

- In `Cloud.C`, changed the OpenMP-off serial move branch to snapshot particle
  pointers into a `List<ParticleType*>` before moving, and moved the first-pass
  `stepFraction` reset into that loop.
- In `dsmcCloud`, cached whether Cartesian reduced-D tracking correction is
  required, replacing the per-parcel string/type check in `dsmcParcel::move()`.

Build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_move_opt_build_20260606.log`

Validation log:

- `log.codex_mpi8_replicatedmesh_mpionly_moveopt_500step_20260606`

Result:

| metric | pre-candidate | candidate |
|---|---:|---:|
| real | 135.30 s | 137.78 s |
| move only | 72.59819434 s | 75.15477516 s |
| buildCellOccupancy | 7.432408543 s | 9.276892873 s |
| collision phase | 12.97955973 s | 15.17461794 s |
| post fields/output | 38.10104302 s | 40.53209741 s |
| full evolve wall | 125.5070652 s | 127.90028 s |
| DLB rebalances | 5 | 4 |

Decision: rejected.  The pointer snapshot did not improve the pure-MPI move
path and worsened the full run.  This source change was manually reverted.

### Rejected candidate: Cartesian-tracking cache only

Source attempt:

- Reverted the `Cloud.C` pointer snapshot change.
- Kept only the `dsmcCloud` cached Cartesian reduced-D tracking flag and the
  corresponding `dsmcParcel::move()` lookup replacement.

Build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_move_opt_cacheonly_build_20260606.log`

Validation log:

- `log.codex_mpi8_replicatedmesh_mpionly_moveopt_cacheonly_500step_20260606`

Result:

| metric | pre-candidate | cache-only |
|---|---:|---:|
| real | 135.30 s | 139.87 s |
| move only | 72.59819434 s | 77.13016431 s |
| buildCellOccupancy | 7.432408543 s | 7.434989898 s |
| collision phase | 12.97955973 s | 14.44935205 s |
| post fields/output | 38.10104302 s | 36.80562033 s |
| full evolve wall | 125.5070652 s | 129.6942547 s |
| DLB rebalances | 5 | 5 |

Decision: rejected.  The cache-only change did not reduce move time in the
formal pure-MPI 500-step run.  This source change was manually reverted.

Revert build log:

- `doc/worklog/v2506/detail_mpi/stage_mpi_move_opt_revert_build_20260606.log`

Post-revert source state:

- The `Cloud.C` serial move branch is back to the pre-candidate direct
  `forAllIter` path.
- `dsmcCloud` no longer contains `constrainCartesianTracking_`.
- `dsmcParcel::move()` again computes the Cartesian/reduced-D condition locally.
- The binary was rebuilt after reverting the rejected candidates.

Technical conclusion:

The attempted micro-optimizations do not explain the gap to `ourmeshbkp`.
The meaningful reference-code difference is the lower-level particle tracking
API: the reference move path uses `particle::trackToAndHitFace()` and
`deviationFromMeshCentre()`, while the current v1706 path still uses the older
`trackToFace(position() + dt*Utracking, td, true)` loop.  Closing the move gap
therefore requires a deliberate particle-tracking API migration or an equivalent
current-API fast path, not another small wrapper-level change.

## Rejected candidate: compatibility `trackToAndHitFace` wrapper

Follow-up attempt:

- added a temporary OF-v1706 compatibility wrapper named
  `particle::trackToAndHitFace()` around the old
  `trackToFace(position_ + displacement, td, true)` implementation;
- added temporary `particle::deviationFromMeshCentre()`;
- changed the uniform-deltaT `dsmcParcel::move()` path to the reference-style
  displacement/fraction call shape.

Logs:

- build: `doc/worklog/v2506/detail_mpi/stage_mpi_tracking_api_compat_build_20260606.log`
- smoke: `log.codex_mpi8_replicatedmesh_mpionly_tracking_api_compat_smoke_10step_20260606`
- formal: `log.codex_mpi8_replicatedmesh_mpionly_tracking_api_compat_500step_20260606`
- revert build: `doc/worklog/v2506/detail_mpi/stage_mpi_tracking_api_compat_revert_build_20260606.log`

Formal result:

| metric | pre-candidate | compatibility API candidate |
|---|---:|---:|
| real | 135.30 s | 140.79 s |
| move only | 72.59819434 s | 82.92498684 s |
| buildCellOccupancy | 7.432408543 s | 8.25212966 s |
| collision phase | 12.97955973 s | 15.08526616 s |
| post fields/output | 38.10104302 s | 37.82146892 s |
| full evolve wall | 125.5070652 s | 131.3125433 s |
| DLB rebalances | 5 | 4 |

Decision: rejected and manually reverted.  The wrapper passed the 10-step pure
MPI smoke but worsened move by `10.32679250 s` in the 500-step run.  This
confirms that merely changing the `dsmcParcel::move()` call shape on top of the
old OF-v1706 tracking core is not enough; the remaining gap is in the lower
tracking algorithm/data path.  Detailed record:
`doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_tracking_api_migration_attempt_20260606.md`.
