# palphd3.3.1react-m8 64-core mixed-parallel performance analysis - 2026-06-19

## Scope

Input logs:

- `results/palphd3.3.1react-m8/omp64_3628404.out`
- `results/palphd3.3.1react-m8/mpi4omp16_3629839.out`
- `results/palphd3.3.1react-m8/mpi8omp8_3629855.out`
- `results/palphd3.3.1react-m8/mpi16omp4_3629840.out`
- `results/palphd3.3.1react-m8/mpi32omp2_3629843.out`
- `results/palphd3.3.1react-m8/mpi64_3629848.out`

All six runs complete `8750` iterations and write final output at `Time = 0.0035`.
The mixed/MPI runs are replicated-mesh runs with Phase C auto DLB enabled.

This note extracts the final summary metrics and separates:

- end-to-end performance conclusions that are safe to compare
- diagnostics that look inconsistent and should not be treated as correctness proof

## Summary table

| Case | Layout | ClockTime [s] | Speedup vs `omp64` | Full evolve wall [s] | Speedup vs `omp64` full evolve | Final particles | Final total energy | DLB rebalances | Particles/rank max/min | Rank wall max/min |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `omp64` | `1 x 64` | 20333 | 1.00 | 19393.22 | 1.00 | 9637612 | 46.94174131 | - | - | - |
| `mpi4omp16` | `4 x 16` | 7251 | 2.80 | 6763.39 | 2.87 | 9621725 | 47.00506479 | 169 | 1.0903 | 2.0568 |
| `mpi8omp8` | `8 x 8` | 5322 | 3.82 | 4931.12 | 3.93 | 9582918 | 47.03540991 | 26 | 1.1358 | 1.1255 |
| `mpi16omp4` | `16 x 4` | 2798 | 7.27 | 2433.12 | 7.97 | 9633056 | 47.13719112 | 143 | 1.2124 | 1.5231 |
| `mpi32omp2` | `32 x 2` | 2087 | 9.74 | 1757.77 | 11.03 | 9471754 | 47.03835509 | 172 | 1.1639 | 1.4294 |
| `mpi64` | `64 x 1` | 1787 | 11.38 | 1571.15 | 12.34 | 9865853 | 31.99176119 | 91 | 1.5292 | 1.2220 |

## Phase breakdown

Fractions below are normalized by `full evolve wall`.

| Case | move only | buildCellOccupancy | collision phase | post fields/output |
| --- | ---: | ---: | ---: | ---: |
| `omp64` | 42.92% | 40.80% | 14.36% | 1.70% |
| `mpi4omp16` | 40.94% | 36.53% | 15.20% | 5.16% |
| `mpi8omp8` | 37.19% | 47.32% | 8.01% | 7.93% |
| `mpi16omp4` | 49.56% | 28.28% | 15.30% | 15.57% |
| `mpi32omp2` | 53.11% | 22.46% | 19.54% | 14.64% |
| `mpi64` | 46.11% | 16.74% | 2.18% | 23.93% |

## Main performance conclusions

1. At the level of total wall time, more MPI ranks clearly help on this machine for this case. The ordering is:
   `mpi64` < `mpi32omp2` < `mpi16omp4` < `mpi8omp8` < `mpi4omp16` << `omp64`.

2. `mpi32omp2` is the fastest configuration that still keeps final total energy in the same band as the other plausible runs.
   Its final total energy is `47.03835509`, close to `omp64 46.94174131`, `mpi4omp16 47.00506479`, `mpi8omp8 47.03540991`, and `mpi16omp4 47.13719112`.

3. `mpi8omp8` has the best rank-wall balance ratio (`1.1255`), but it is not the fastest run.
   This is a reminder that better balance does not automatically mean lower end-to-end time.

4. `mpi16omp4` and especially `mpi32omp2` cut `buildCellOccupancy` cost sharply relative to `omp64`.
   This is a major contributor to the speedup trend.

5. `omp64` is the clear worst performer.
   Its final `ClockTime` is `20333 s`, versus `2087 s` for `mpi32omp2`.

## DLB behavior

1. DLB frequency is highly configuration-dependent.
   Rebalances are `169` for `mpi4omp16`, `26` for `mpi8omp8`, `143` for `mpi16omp4`, `172` for `mpi32omp2`, and `91` for `mpi64`.

2. `mpi8omp8` reaches the best rank-wall balance with far fewer rebalances than `mpi16omp4` and `mpi32omp2`.
   That means the lower wall time of `mpi16omp4` and `mpi32omp2` is not simply “more DLB is better”; their decomposition is just a better overall fit for this workload.

3. `mpi64` still ends with a relatively poor particle balance (`1.5292`) and a final reported maximum imbalance of `21.77300838%`, but total time is still the lowest.
   That suggests the dominant cost for this run is not captured by particle-count balance alone.

## Data quality warnings

### 1. `mpi64` should not be treated as correctness-confirmed

`mpi64` is the fastest run, but its final thermodynamic state is far away from the other five runs:

- final total energy is `31.99176119`, while the others are all around `46.94 - 47.14`
- final average linear kinetic energy is `2.892352666e-19`, while the others are around `4.05e-19 - 4.14e-19`
- final average rotational energy is `~3.00e-20`, while the others are around `~5.13e-20`
- final average vibrational energy is `~5.03e-21`, while the others are around `~3.06e-20`

This is too large to explain as normal statistical drift. `mpi64` may be fast, but with the current evidence it is not a safe “best configuration” conclusion.

### 2. Some cumulative collision diagnostics are internally inconsistent

Two runs show clearly suspicious cumulative summaries:

- `mpi8omp8`: `collision cumulative global collisions = 973573119`, `collision cumulative global candidates = 52435904`, acceptance `18.56691779`
- `mpi4omp16`: `collision cumulative global collisions = 1173573770`, `collision cumulative global candidates = 1605173817`, acceptance `0.7311194324`

The `mpi8omp8` cumulative acceptance is physically impossible as printed because collisions exceed candidates by a large factor.
This means the cumulative collision diagnostics are not trustworthy for direct cross-run comparison in their current form.

Per-step local collision lines near output also differ strongly across layouts because replicated-mesh output is using rank-0 local diagnostics in these logs. Those lines should not be interpreted as global collision totals.

## Practical ranking

If the objective is pure wall time and correctness is ignored, the logs rank:

1. `mpi64`
2. `mpi32omp2`
3. `mpi16omp4`
4. `mpi8omp8`
5. `mpi4omp16`
6. `omp64`

If the objective is performance with at least basic cross-run plausibility of the final thermodynamic state, the safer ranking is:

1. `mpi32omp2`
2. `mpi16omp4`
3. `mpi8omp8`
4. `mpi4omp16`
5. `omp64`

`mpi64` should be excluded from the formal recommendation until the energy/state discrepancy is explained.

## Recommended next checks

1. Diff the runtime controls used by `mpi64` against `mpi32omp2` and `mpi16omp4`, especially chemistry, sampling, and any collision/output reduction switches.
2. Re-run one confirmation case for `mpi32omp2` and `mpi64` with the same source/binary state and preserved control snapshots.
3. For any formal report, use total wall time, full-evolve wall time, move/build/collision split, final particle count, and final energy as the primary comparison set.
4. Do not use the cumulative collision summary in these logs as correctness evidence until the reduction/output path is audited.
