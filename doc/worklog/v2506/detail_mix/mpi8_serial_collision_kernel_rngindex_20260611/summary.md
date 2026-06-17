# MPI8 replicated serial collision kernel optimization

Date: 2026-06-11

## Scope

Case:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Final test directory:

- `doc/worklog/v2506/detail_mix/mpi8_serial_collision_kernel_rngindex_20260611`

Baseline:

- `doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean/aggregate.csv`
- baseline row: `ourmesh,MPI8`, 3 runs, 500 steps, `useOpenMP false`

## Code change

Changed only the serial fallback in:

- `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C`

Main changes:

- Reused the existing thread-0 scratch lists for serial collision instead of allocating per-cell temporary lists.
- Traversed `cloud_.occupancyCollisionCells()` instead of scanning every mesh cell in the serial collision loop.
- Cached per-cell parcel pointers, `typeId`, and charge for the serial collision kernel.
- Used `pairModelAddressing()[typeIdP][typeIdQ]` instead of the `returnModelId(parcelP, parcelQ)` wrapper in the hot path.
- Replaced the serial fast-RNG hot path wrapper with `randomIndex(n)` for 0..n-1 indexing.

The OMP/mixed branch was left unchanged.

## Correctness

Final 3 runs:

| rep | exit | iterations | End main | particles | collisions | candidates | total energy |
|---:|---:|---:|---|---:|---:|---:|---:|
| 1 | 0 | 500 | true | 2463624 | 35726 | 64620 | 1.241986983 |
| 2 | 0 | 500 | true | 2463653 | 35391 | 67784 | 1.242416101 |
| 3 | 0 | 500 | true | 2463800 | 35318 | 67704 | 1.242474576 |

Fatal-marker scan:

- no `FOAM FATAL`
- no `Segmentation`
- no `MPI_ABORT`
- no `Killed`
- no `ERROR` / `Error`

## Performance

All values are 3-run means.

| metric | clean MPI8 mean [s] | serial kernel mean [s] | delta [s] | delta |
|---|---:|---:|---:|---:|
| external real | 101.446667 | 96.056667 | -5.390000 | -5.31% |
| full evolve wall | 98.629420 | 95.321816 | -3.307604 | -3.35% |
| move+collide wall | 93.693971 | 90.135724 | -3.558247 | -3.80% |
| move only | 69.070791 | 70.924566 | +1.853775 | +2.68% |
| buildCellOccupancy | 10.613061 | 10.291471 | -0.321590 | -3.03% |
| collision phase | 32.798558 | 28.400884 | -4.397674 | -13.41% |
| post fields/output | 5.003583 | 5.205413 | +0.201830 | +4.03% |

Per-run timings:

| rep | real [s] | full evolve [s] | move [s] | buildOcc [s] | collision [s] |
|---:|---:|---:|---:|---:|---:|
| 1 | 93.20 | 94.214825 | 68.185453 | 9.705546 | 29.738871 |
| 2 | 99.79 | 98.263434 | 74.932610 | 10.851112 | 30.729862 |
| 3 | 95.18 | 93.487189 | 69.655635 | 10.317755 | 24.733921 |

## Conclusion

The final `useOpenMP false` replicated MPI8 path now shows both collision and end-to-end speedup on the 500-step `ourmesh/mpi8replicatedmesh` benchmark. The dominant direct gain is in `collision phase` (-13.41%). End-to-end `real` improves by 5.31%, despite `move only` still being slightly slower than the clean mean.

The next larger target is still the replicated move path, because move remains the largest phase and absorbs part of the collision gain.
