# palphd3.3.1react-m8-fixoutput performance analysis - 2026-06-21

## Scope

Input logs:

- `results/palphd3.3.1react-m8-fixoutput/omp64_3628404.out`
- `results/palphd3.3.1react-m8-fixoutput/mpi16omp4_3632046.out`
- `results/palphd3.3.1react-m8-fixoutput/mpi64_3632049.out`

Reference older logs for same layouts:

- `results/palphd3.3.1react-m8/omp64_3628404.out`
- `results/palphd3.3.1react-m8/mpi16omp4_3629840.out`
- `results/palphd3.3.1react-m8/mpi64_3629848.out`

All three `fixoutput` runs complete `8750` iterations and end at `Time = 0.0035`.

## Final summary

| Case | Layout | ClockTime [s] | Full evolve wall [s] | Final particles | Final total energy | DLB rebalances | Rank wall max/min |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `omp64` | `1 x 64` | 20333 | 19393.22 | 9637612 | 46.94174131 | - | - |
| `mpi16omp4` | `16 x 4` | 1680 | 1399.62 | 9524016 | 46.9733164 | 30 | 1.1623 |
| `mpi64` | `64 x 1` | 1985 | 1727.96 | 9554809 | 47.16757125 | 69 | 1.2212 |

## Main conclusions

1. In the `fixoutput` series, `mpi16omp4` is now the fastest of the tested configurations.
   It beats `mpi64` on both end-to-end wall time and full-evolve wall time:
   `ClockTime 1680 s` vs `1985 s`,
   `full evolve wall 1399.62 s` vs `1727.96 s`.

2. Both mixed/MPI runs now finish with physically plausible final state metrics.
   The severe anomaly from the older `mpi64_3629848.out` run is gone.

3. Compared with `omp64`, the speedups are large:
   `mpi16omp4` gives about `12.10x` speedup in `ClockTime` and `13.86x` in `full evolve wall`;
   `mpi64` gives about `10.24x` and `11.22x` respectively.

## Phase breakdown

Fractions below are normalized by `full evolve wall`.

| Case | move only | buildCellOccupancy | collision phase | post fields/output |
| --- | ---: | ---: | ---: | ---: |
| `omp64` | 42.92% | 40.80% | 14.36% | 1.70% |
| `mpi16omp4` | 48.35% | 7.48% | 22.53% | 27.60% |
| `mpi64` | 47.06% | 18.47% | 20.93% | 18.40% |

Observations:

- `mpi16omp4` has a much smaller `buildCellOccupancy` share than both `omp64` and `mpi64`.
- `mpi16omp4` pays more in `post fields/output` share than `mpi64`, but still wins overall because move plus occupancy is much cheaper.

## CPU-timer interpretation

The `collision phase cpu [s]` values are not directly comparable between pure OpenMP and replicated raw-MPI layouts.

- `omp64`: `collision phase cpu = 177203.6`, `collision phase wall = 2785.44`, ratio about `63.6`, consistent with one process accumulating CPU across about 64 threads.
- `mpi16omp4`: `collision phase cpu = 1261.16`, `collision phase wall = 315.35`, ratio about `4.0`, consistent with one rank accumulating CPU across 4 threads.
- `mpi64`: `collision phase cpu = 360.28`, `collision phase wall = 361.64`, ratio about `1.0`, consistent with one single-thread rank.

For replicated raw-MPI runs, the profile summary reduces CPU timers with `MPI_MAX`, so the printed CPU value is the slowest-rank local CPU time, not the sum across all ranks.

## Comparison against the older `react-m8` series

### `mpi64`

Old run:

- `ClockTime = 1787 s`
- `full evolve wall = 1571.15 s`
- `Total energy = 31.99176119`
- `collision phase wall = 34.28 s`
- `collision phase cpu = 34.22 s`

New `fixoutput` run:

- `ClockTime = 1985 s`
- `full evolve wall = 1727.96 s`
- `Total energy = 47.16757125`
- `collision phase wall = 361.64 s`
- `collision phase cpu = 360.28 s`

Interpretation:

- The new `mpi64` is slower than the old one.
- But the old `mpi64` had an obviously wrong final thermodynamic state and an implausibly tiny collision phase.
- The new `mpi64` is the first `mpi64` result in these folders that looks comparable to the other layouts in final-state terms.

### `mpi16omp4`

Old run:

- `ClockTime = 2798 s`
- `full evolve wall = 2433.12 s`
- `Total energy = 47.13719112`
- `DLB rebalances = 143`
- `rank wall max/min = 1.5231`

New `fixoutput` run:

- `ClockTime = 1680 s`
- `full evolve wall = 1399.62 s`
- `Total energy = 46.9733164`
- `DLB rebalances = 30`
- `rank wall max/min = 1.1623`

Interpretation:

- `mpi16omp4` improves very strongly in the `fixoutput` series.
- Wall time is substantially better and the DLB behavior is much calmer.
- The final state stays in a plausible band close to `omp64`.

## Recommendation

For the currently available `fixoutput` results, the best-performing and plausibly-correct tested configuration is:

1. `mpi16omp4`
2. `mpi64`
3. `omp64`

## Caution

This ranking is only for the three tested layouts present in `results/palphd3.3.1react-m8-fixoutput`.
It does not prove that `mpi16omp4` would beat `mpi32omp2` if `mpi32omp2` were rerun under the same `fixoutput` source/control state.
