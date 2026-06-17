# MPI8 collision detail/load analysis, 2026-06-12

## Scope

- Case: `zb-cylinder-react`, 300 steps (`endTime 1.9920146682e-05`, `deltaT 6.640048894e-08`).
- Baseline source state: local dirty tree with owned collision-cell iterator and current noTimeCounter candidate-clear changes.
- Run mode for detailed analysis: replicated mesh raw MPI, `np=8`, `OMP_NUM_THREADS=1`, `profileSummary true`, `profileDetail true`, `replicatedMeshDLBProfile true`.
- The 300-step detail case was run from a `/tmp` hardlink copy of `run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`; durable artifacts are under this worklog directory.

## Instrumentation changes

1. `dsmcCloud::printProfileSummary()` now treats replicated raw MPI as `replicatedMeshActive() && replicatedMesh_->nProcs() > 1`, not `!Pstream::parRun()`.
   - Reason: `dsmcReplicatedMesh::initialize()` calls `UPstream::init()` when MPI was not initialized, which makes `Pstream::parRun()` true even though this is not a normal OpenFOAM `-parallel` case.
   - Effect: final rank profile gather now prints `parcels`, `ownedCollCells`, `finalCandidates`, `full`, `move+collide`, `move`, `build`, `collision`, `post`.
2. `noTimeCounter` now keeps per-rank cumulative local collisions/candidates and, under `profileDetail true`, gathers them at terminal-output intervals.
   - This is diagnostic only.
   - It prints current-window `localColl/localCand` and cumulative `cumColl/cumCand` by rank.

## Artifacts

- Serial/OMP/MPI baseline summary: `doc/worklog/v2506/detail_mix/zb_serial_baseline_20260612/zb_serial_collision_baseline_20260612.md`
- Rank profile smoke before fix: `logs/zb_MPI8_detail_smoke10_20260612.log`
- Rank profile smoke after raw-MPI fix: `logs/zb_MPI8_detail_smoke10_after_rawmpi_fix_20260612.log`
- 300-step rank profile after raw-MPI fix: `logs/zb_MPI8_detail_rank_after_rawmpi_fix_300step_20260612.log`
- 300-step candidate-load profile: `logs/zb_MPI8_collision_load_300step_20260612.log`
- Matching `.exit` files all contain `0`.

## Baseline speedups

The serial and OMP rows come from the same current-source 300-step baseline summary. The MPI candidate-load row includes extra diagnostic gather output, so use it for load explanation, not as a low-overhead production timing.

| mode | real [s] | full evolve [s] | collision [s] | real speedup vs serial | collision speedup vs serial |
| --- | ---: | ---: | ---: | ---: | ---: |
| Serial 1 core | 342.78 | 370.9618 | 49.2992 | 1.00x | 1.00x |
| OMP8 current | 74.05 | 72.01 | 11.39 | 4.63x | 4.33x |
| MPI8 rank-detail, no candidate-load gather | 77.61 | 79.5623 | 31.7615 | 4.42x | 1.55x |
| MPI8 candidate-load detail | 79.09 | 80.8120 | 35.0456 | 4.33x | 1.41x |

Key comparison: even after owned collision-cell traversal, MPI8 collision remains about `2.8-3.1x` slower than OMP8 collision in these current-source runs.

## MPI8 rank profile, 300-step candidate-load run

From `logs/zb_MPI8_collision_load_300step_20260612.log`:

| rank | final parcels | owned coll cells | final candidates | full [s] | move [s] | build [s] | collision [s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 217542 | 6159 | 26154 | 80.7902 | 42.4026 | 5.8625 | 26.9494 |
| 1 | 248726 | 7922 | 49102 | 80.7946 | 45.6135 | 5.0427 | 24.9989 |
| 2 | 216097 | 8423 | 62020 | 80.8002 | 46.3809 | 5.1856 | 24.2086 |
| 3 | 218125 | 6090 | 23305 | 80.7836 | 45.5581 | 4.6355 | 25.5087 |
| 4 | 298877 | 6197 | 105122 | 80.8120 | 51.8327 | 6.2974 | 18.1110 |
| 5 | 213043 | 7690 | 46706 | 80.7745 | 40.0356 | 3.6529 | 31.2918 |
| 6 | 290631 | 5890 | 47297 | 80.8094 | 37.4425 | 5.7620 | 31.6096 |
| 7 | 255445 | 6321 | 17904 | 80.7792 | 35.6970 | 3.8435 | 35.0456 |

- Full wall max/min: `1.000463931`.
- Collision wall max/min: `1.935047402`.
- Particles per rank max/min: `1.434536119`.
- Rank wall time from replicated mesh report: min `71.69`, max `71.87`, max/min `1.00251081`.

## Candidate-load evidence

Final terminal-output window at step 300:

| rank | localColl | localCand | cumColl | cumCand | cum acceptance |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 13680 | 26154 | 3849500 | 6599658 | 0.5833 |
| 1 | 22222 | 49102 | 4278565 | 8486185 | 0.5042 |
| 2 | 28330 | 62020 | 4458501 | 8769912 | 0.5084 |
| 3 | 11913 | 23305 | 2978781 | 5531274 | 0.5385 |
| 4 | 58184 | 105122 | 7866505 | 13951475 | 0.5638 |
| 5 | 21673 | 46706 | 7756903 | 15343425 | 0.5056 |
| 6 | 27555 | 47297 | 3135053 | 5210950 | 0.6016 |
| 7 | 9585 | 17904 | 4067411 | 7442975 | 0.5465 |

- Final-window local candidates max/min: `5.87142538`.
- Cumulative candidates max/min at step 300: `2.944458304`.
- Global collisions/candidates at step 300 output: `193142 / 377610`.

The inter-DLB actual load block after rebalance shows why this was hidden:

| rank | move [s] | coll [s] | total [s] |
| --- | ---: | ---: | ---: |
| 0 | 10.1884 | 9.9364 | 20.1248 |
| 1 | 13.0739 | 7.0176 | 20.0915 |
| 2 | 13.2227 | 6.7561 | 19.9787 |
| 3 | 14.8683 | 5.0630 | 19.9313 |
| 4 | 16.1476 | 3.7725 | 19.9201 |
| 5 | 10.0613 | 9.9097 | 19.9710 |
| 6 | 9.8103 | 9.3451 | 19.1555 |
| 7 | 9.2767 | 10.6803 | 19.9570 |

Total max/min is only `1.050603366`, while collision and candidate work are not balanced.

## Interpretation

MPI8 collision is not "doing no work"; it is not shortening the collision critical path like OMP8.

The measured reason is not gross total rank imbalance. DLB is balancing move+collision total time very well. The problem is that this balance hides an anti-correlation: ranks with larger move time often have smaller collision time and vice versa. The final full-evolve max/min is nearly 1.0, but collision-only max/min is about 1.94 and cumulative candidate max/min is about 2.94.

OMP8 gets a real shared-memory collision-kernel speedup because one process can dynamically schedule active collision cells across 8 threads. MPI8 replicated mesh gives each rank a static owned subset between DLB events; the slowest rank's single-thread collision time remains the critical path. Balancing total move+collision can make total wall time look good while leaving collision under-parallelized.

This explains why OMP8 collision cannot be fully reproduced by MPI8 with only owned-cell traversal and configuration switches. MPI8 needs a collision-aware partition/offload policy, not just total productive-time DLB.

## Next optimization targets

1. Add a collision-load constraint or score using cumulative local candidates/collision wall, not only particle/move-weight proxies.
2. Rebalance on collision-specific imbalance separately from total move+collision imbalance.
3. Consider a small intra-rank OpenMP collision kernel for MPI ranks, because MPI-only ownership cannot dynamically split a hot owned rank's cells within a step.
4. Keep the owned collision-cell iterator; it is necessary, but not sufficient.
