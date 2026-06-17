# Serial1 collision speedup check

Date: 2026-06-11

## Purpose

The question was whether the large OMP8 vs MPI8 replicated collision gap is caused by MPI load distribution rather than simply by an unoptimized MPI collision kernel.

The test adds a one-core serial reference and compares collision speedup:

- serial1: one process, `useOpenMP false`, `openmpThreads 1`, `replicatedMesh false`
- OMP8: one process, `useOpenMP true`, `openmpThreads 8`, `replicatedMesh false`
- MPI8 replicated: current cached-reference MPI8 result from `doc/worklog/v2506/detail_mix/mpi8_collision_cached_refs_20260611/summary.md`

All runs use the same 500-step no-write `ourmesh` setup and `collisionFastRng true`.

## Case setup

Serial case was copied from:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8`

Temporary case:

- `/tmp/hystrath_serial1_ompbase_20260611`

Control changes:

- `useOpenMP false`
- `openmpThreads 1`
- `replicatedMesh false`

## Results

| mode | exit | real | full evolve | move+collide | move | buildOcc | collision | post | particles | collisions | candidates | energy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| serial1 | 0 | 283.07 | 302.392635 | 301.756149 | 258.861539 | 27.853350 | 14.636285 | 0.633565 | 2463557 | 35746 | 68578 | 1.243308561 |
| OMP8 current | 0 | 59.21 | 57.779304 | 57.136494 | 47.208625 | 6.397944 | 3.003647 | 0.639841 | 2463637 | 35672 | 68655 | 1.243041817 |
| MPI8 replicated current mean | 0 | 92.69 | 92.260071 | 87.345999 | 70.295910 | 10.200031 | 25.957368 | 4.945385 | - | - | - | - |

Fatal-marker scan:

- no `FOAM FATAL`
- no `Segmentation`
- no `MPI_ABORT`
- no `Killed`
- no `ERROR` / `Error`

## Collision speedup

| comparison | collision time | speedup vs serial1 |
|---|---:|---:|
| serial1 | 14.636285 | 1.000000 |
| OMP8 current | 3.003647 | 4.872838 |
| MPI8 replicated current mean | 25.957368 | 0.563859 |

Equivalently:

- OMP8 collision is `4.87x` faster than serial1.
- MPI8 replicated collision is not faster than serial1; it is `1.77x` slower than serial1.

## Interpretation

This result strongly supports the load-distribution/synchronization explanation.

If MPI8 replicated collision were truly splitting the same collision workload across 8 ranks with a comparable single-rank kernel, its max-rank collision wall time should be well below the one-core serial reference. Instead, the MPI8 replicated collision phase is slower than serial1.

That means the remaining MPI8 collision gap is not primarily a missing micro-optimization in the serial kernel. The likely causes are:

- collision candidates are distributed very unevenly across owned cells/ranks;
- candidate work scales roughly with `nC*(nC - 1)`, so particle-count balance can hide collision imbalance;
- MPI replicated collision is bounded by the slowest rank and synchronization around the collision phase;
- OMP8 uses per-step dynamic scheduling over collision cells, while MPI8 uses coarse cell ownership.

## Next diagnostic

The next direct measurement should print per rank:

- `ownedCollCells.size()`
- sum of `nCandidatesPerCell` over owned cells
- max/top cell candidate counts
- local collision wall time

That will show whether one or a few ranks carry most of the collision candidate work.
