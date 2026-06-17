# MPI8 replicated collision kernel follow-up

Date: 2026-06-11

## Scope

Case:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Source state:

- Includes the earlier owned collision-cell iterator and serial fast collision kernel work.
- This follow-up keeps the real MPI8 replicated-mesh path (`openmpThreads 1`, `replicatedMesh true`) and does not use OpenMP to explain the MPI result.

## Code changes

Files:

- `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.H`
- `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C`

Changes:

- Removed the unused OMP collision `threadVelocities_` scratch list. The cached velocities were written but never read.
- Added `candidateCellsToClear_` plus per-thread `threadCandidateCells_` so `nCandidatesPerCell` clears only cells that had non-zero candidates on the previous call.
- Cleared stale candidate counts even when binary collision is inactive.
- Cached hot-path references outside the cell/candidate loops:
  - `constProps`
  - `binaryCollision`
  - `reactions`
  - `pairModelAddressing`
  - `nParticles`
  - `collisionSelectionRemainder`
  - `sigmaTcRMax`
  - `cellVolumes`
- Replaced repeated `cloud_.constProps(typeId)` calls with direct `constProps[typeId]` indexing after `typeId` has already been read from parcels.

## Build

Command:

```bash
source doc/scripts/build-dsmcFoam.sh
```

Result:

- build passed
- `dsmcFoam+` rebuilt successfully

## MPI8 candidate-clear intermediate

Directory:

- `doc/worklog/v2506/detail_mix/mpi_coll_candidate_clear_followup_20260611`

This was the first follow-up after adding precise candidate-count clearing and removing `threadVelocities_`, before cached hot-path references.

| rep | real | full evolve | move+collide | move | buildOcc | collision | post |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 94.57 | 94.355692 | 89.342542 | 71.093848 | 8.741130 | 30.723736 | 5.037581 |
| 2 | 91.53 | 89.745954 | 84.874939 | 66.089372 | 9.199431 | 27.842282 | 4.891650 |
| 3 | 94.41 | 92.466811 | 87.525419 | 69.064399 | 9.857113 | 33.212576 | 5.028937 |
| mean | 93.503333 | 92.189486 | 87.247633 | 68.749207 | 9.265891 | 30.592865 | 4.986056 |

Interpretation:

- Correctness passed, but collision mean did not improve over the previous serial-kernel mean.
- This confirmed the candidate clearing fix is useful for state hygiene, not enough by itself for collision speed.

## MPI8 cached-reference final run

Directory:

- `doc/worklog/v2506/detail_mix/mpi8_collision_cached_refs_20260611`

| rep | exit | real | full evolve | move+collide | move | buildOcc | collision | post | particles | collisions | candidates | energy | rebalances |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0 | 93.47 | 94.801022 | 89.906220 | 72.967099 | 10.628191 | 28.619710 | 4.934336 | 2463648 | 35032 | 67614 | 1.242479707 | 4 |
| 2 | 0 | 90.74 | 89.581171 | 84.665223 | 68.108751 | 9.984661 | 23.985004 | 4.929612 | 2463778 | 34481 | 65947 | 1.241685421 | 3 |
| 3 | 0 | 93.86 | 92.398019 | 87.466555 | 69.811881 | 9.987243 | 25.267391 | 4.972206 | 2463630 | 35287 | 67650 | 1.242255131 | 3 |
| mean | - | 92.690000 | 92.260071 | 87.345999 | 70.295910 | 10.200031 | 25.957368 | 4.945385 | - | - | - | - | - |

Fatal-marker scan:

- no `FOAM FATAL`
- no `Segmentation`
- no `MPI_ABORT`
- no `Killed`
- no `ERROR` / `Error`

## Comparison

Baseline means:

- clean MPI8 mean from `doc/worklog/v2506/detail_mix/mpi8_serial_collision_kernel_rngindex_20260611/summary.md`
- previous serial-kernel mean from the same summary

| metric | clean MPI8 mean | cached-ref mean | delta vs clean |
|---|---:|---:|---:|
| external real | 101.446667 | 92.690000 | -8.63% |
| full evolve wall | 98.629420 | 92.260071 | -6.46% |
| move+collide wall | 93.693971 | 87.345999 | -6.78% |
| move only | 69.070791 | 70.295910 | +1.77% |
| buildCellOccupancy | 10.613061 | 10.200031 | -3.89% |
| collision phase | 32.798558 | 25.957368 | -20.86% |
| post fields/output | 5.003583 | 4.945385 | -1.16% |

| metric | previous serial-kernel mean | cached-ref mean | delta |
|---|---:|---:|---:|
| external real | 96.056667 | 92.690000 | -3.50% |
| collision phase | 28.400884 | 25.957368 | -8.60% |

## Conclusion

This version satisfies the MPI8 replicated-mesh target more clearly than the previous serial-kernel step:

- real MPI8 total time improves from clean 101.45s to 92.69s
- MPI8 collision phase improves from clean 32.80s to 25.96s
- relative to the previous serial-kernel version, the cached-reference kernel still improves both total time and collision phase

The remaining limit is still Amdahl-dominated by move/build and replicated-mesh migration/check overhead. Collision is now no longer the only large direct target; the next practical target for end-to-end speed is the replicated move path and migration-side bookkeeping.
