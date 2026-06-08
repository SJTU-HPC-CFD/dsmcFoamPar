# MPI replicated mesh DLB implementation - 2026-06-06

## Scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Reference source: `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
- Build log: `doc/worklog/v2506/detail_mpi/stage_mpi_replicated_mesh_dlb_build_20260606.log`
- Formal run log: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_dlb_500step_20260606`

## Implemented changes

- `dsmcVolFields` now uses replicated owner cells for shared sample cache build, density-only OpenMP sampling, field combine, and cell reduction.
- `dsmcVolFields` now precomputes owned boundary faces and accumulates boundary fluxes only on owner ranks before the replicated output reduce.
- Replicated output-time field reduce still runs collectively, but only the output rank computes/writes output fields; all ranks still execute output reset.
- `replicatedMeshDelayedReceive` is wired in `dsmcCloud`: periodic replicated migration can use `migrateBegin()`/`migrateFinish()`, with pending async migration finished before manual DLB, auto DLB, and solver output.
- Raw/replicated profile printing now uses a rank-reduced max path and output-rank gating; collision diagnostics include candidates and acceptance rate.
- The fixed case `system/controlDict` was restored to the full replicated mesh DLB configuration from `ourmeshbkp`.

## Validation

Build:

- Command: `source doc/scripts/env.sh; source doc/scripts/build-dsmcFoam.sh`
- Result: success. The log contains existing OpenFOAM template warnings only; no compile/link errors.

Formal np8 500-step run:

- Command: `source doc/scripts/env.sh; cd .../ourmesh/mpi8replicatedmesh; /usr/bin/time -p mpirun -np 8 dsmcFoam+`
- Result: success, `real 152.07`.
- Final particles: `2,463,839`.
- Final collisions: `34,543`; candidates: `60,646`; acceptance rate: `0.569584144`.
- Final total energy: `1.240961863`.
- Solver profile: move+collide `107.1024787 s`, move only `82.58246361 s`, buildCellOccupancy `9.800232796 s`, collision `15.99929084 s`, post fields/output `43.16658299 s`.
- Replicated DLB profile: migration calls `55`; particles per rank min `223,516`, max `371,279`, max/min `1.661084665`; checks `500`; rebalances `4`; Phase C wall max `1.163136352 s`.

Output-path smoke:

- Temporary settings: `endTime 1.e-06`, `writeInterval 1.e-06`, then restored to the formal 500-step values.
- Field output was generated under `1e-06/` and includes mixture fields such as `p_mixture`, `Ttra_mixture`, `U_mixture`, `dsmcN_mixture`, wall heat/shear fields, and `dsmcSigmaTcRMax`.
- The run was manually terminated after the field output point because the later solver-level parcel gather/write stage was still running for more than three minutes. This produced the `BAD TERMINATION` lines in the smoke log. Treat this as a remaining parcel-output performance/operational risk, not a field-reduce failure.

## Comparison and residual risk

- Compared with the earlier owner-filter 500-step run (`real 184.13`, post fields/output `75.75 s`, final particles `2,463,897`), this run keeps the particle count consistent and reduces post fields/output to `43.17 s`.
- It is still slower than the historical reference target (`External wall seconds 117.31 s`), mainly due to move/build/post costs rather than Phase C DLB overhead.
- The fixed formal case does not normally trigger solver output because `writeInterval 1.e-3` is greater than `endTime 5.e-05`. Field output was separately smoke-tested; full parcel gather/write remains expensive and should be optimized or gated separately if frequent replicated output is required.

## 2026-06-06 OpenMP2 replicated-mesh DLB tuning

Reason:

- The first formal DLB run had `useOpenMP false`, so the hot replicated move, collision, and field sampling paths ran mostly serial inside each raw-MPI rank.
- Enabling `useOpenMP true` with `openmpThreads 8` was rejected. The 10-step smoke completed, but the 500-step candidate was terminated at iteration 380 after `ClockTime = 281 s`, clearly worse than the `152.07 s` baseline. The log is `log.codex_mpi8_replicatedmesh_dlb_openmp_500step_20260606` and contains the expected signal-15 termination lines.
- `openmpThreads 2` was then tested to avoid 8 MPI ranks oversubscribing the node with 64 OpenMP worker threads.

Smoke logs:

- `log.codex_mpi8_replicatedmesh_dlb_openmp_smoke_10step_20260606`: 8 threads, `real 18.04`, move `3.356057359 s`, collision `2.203699577 s`, post fields/output `1.661658776 s`.
- `log.codex_mpi8_replicatedmesh_dlb_openmp2_smoke_10step_20260606`: 2 threads, `real 10.70`, move `1.288977197 s`, collision `0.08681908 s`, post fields/output `0.914833782 s`.
- `log.codex_mpi8_replicatedmesh_dlb_openmp4_smoke_10step_20260606`: 4 threads, `real 13.41`, move `2.449368408 s`, collision `1.708530312 s`, post fields/output `1.188395968 s`; rejected.
- `log.codex_mpi8_replicatedmesh_dlb_openmp2_dynamicmove_smoke_10step_20260606`: 2 threads with dynamic move schedule, `real 10.82`, move `1.274089906 s`, collision `0.098398883 s`, post fields/output `0.951590381 s`; rejected because end-to-end smoke was slower than static.

Final controlDict choice:

- `useOpenMP true`
- `openmpThreads 2`
- `openmpMoveSchedule static`
- `openmpMoveChunk 64`
- `openmpCollisionSchedule dynamic`
- `openmpCollisionChunk 8`

Formal np8 500-step run:

- Log: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_dlb_openmp2_500step_20260606`
- Result: success, `real 113.48`.
- Final particles: `2,463,653`.
- Final total energy: `1.240928489`.
- Solver profile: move+collide `81.23246265 s`, move only `69.34624257 s`, buildCellOccupancy `5.975890432 s`, collision `5.523175253 s`, post fields/output `33.6372609 s`.
- OpenMP profile: enabled, max threads `2`, move enabled with `static, chunk 64`, collision `dynamic, chunk 8`.
- Replicated DLB profile: migration calls `55`; particles per rank min `226,916`, max `444,982`, max/min `1.960998784`; checks `500`; rebalances `4`; Phase C wall max `1.083275751 s`.

Comparison:

- Versus the immediate serial/OpenMP-off run: `152.07 -> 113.48 s`, a `25.38%` reduction in wall time.
- Versus `ourmeshbkp/mpi8origin` historical baseline: `160.95 -> 113.48 s`, a `29.49%` reduction.
- Versus `ourmeshbkp/mpi8replicatedmesh` reference: `117.31 -> 113.48 s`, `3.26%` faster end-to-end.
- Remaining gap is concentrated in move: current move `69.35 s` versus reference move `47.17 s`. Collision is now faster than the reference (`5.52 s` versus `10.90 s`), and post fields/output is close (`33.64 s` versus `30.36 s`).
