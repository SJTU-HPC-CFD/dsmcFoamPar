# No-write compute performance repeat summary - 2026-06-10

## Scope

This report summarizes three repeat runs for each requested 8-core mode:

- `ourmesh`: 500 steps (`endTime 5.e-05`, `deltaT 1.e-07`)
- `zb-cylinder-react`: 300 steps (`endTime 1.9920146682e-05`, `deltaT 6.640048894e-08`)
- output writing disabled for these step ranges by keeping `writeControl runTime` and setting `writeInterval 1.e-3`, which is larger than both case end times
- `OMP8`, `MPI2xOMP4`, `MPI4xOMP2`, and `MPI8` use the replicated-mesh/raw-MPI path and do not use OpenFOAM `-parallel`
- `MPI8origin` is the standard OpenFOAM decomposed baseline, run as `mpirun -np 8 dsmcFoam+ -parallel`

Logs and control snapshots are under:

```text
doc/worklog/v2506/detail_mix/ourmesh_full_series_current_mpi8cfg_repeat3_20260618
```

## Correctness Status

All parsed runs passed the correctness gates: exit 0, expected iteration count, `End main`, no fatal/segmentation/MPI abort/NaN markers, nonzero final particles, zero stuck particles, and finite total energy.

## ourmesh Performance

Expected iterations: `500`.

| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 3/3 | 63.23 | 1.20 | 62.27-64.58 | 57.55 | 47.41 | 6.21 | 2.80 | 0.61 |
| MPI2xOMP4 | 3/3 | 67.43 | 0.31 | 67.15-67.77 | 61.79 | 50.32 | 7.45 | 3.23 | 0.67 |
| MPI4xOMP2 | 3/3 | 66.09 | 1.13 | 64.95-67.20 | 59.50 | 48.91 | 6.56 | 3.11 | 0.72 |
| MPI8 | 3/3 | 84.25 | 0.91 | 83.59-85.29 | 75.28 | 60.46 | 10.43 | 4.85 | 0.91 |
| MPI8origin | 3/3 | 126.63 | 0.97 | 125.73-127.65 | 126.50 | 109.11 | 15.21 | 11.49 | 0.46 |

Best mean external `real`: `OMP8` at `63.23 s`.

| mode | mean real delta vs OMP8 | percent |
|---|---:|---:|
| OMP8 | 0.00 s | 0.00% |
| MPI2xOMP4 | 4.20 s | 6.64% |
| MPI4xOMP2 | 2.86 s | 4.52% |
| MPI8 | 21.02 s | 33.24% |
| MPI8origin | 63.39 s | 100.25% |

### ourmesh Correctness Metrics

| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |
|---|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 2463661 | 2463653-2463667 | 35530 | 68713 | 1.2430686303 | 1.2429961450-1.2431566040 | 0 |
| MPI2xOMP4 | 2463668 | 2463586-2463725 | 7570 | 12972 | 1.2431412250 | 1.2431095240-1.2431945820 | 0 |
| MPI4xOMP2 | 2463677 | 2463625-2463724 | 3484 | 5934 | 1.2428320623 | 1.2427814920-1.2429191140 | 0 |
| MPI8 | 2463656 | 2463591-2463767 | 589 | 876 | 1.2426806073 | 1.2425475250-1.2427948950 | 0 |
| MPI8origin | 2463678 | 2463631-2463762 | 35460 | 68272 | 1.2432886157 | 1.2431735070-1.2433772370 | 0 |

### ourmesh Load-Balance Diagnostics

| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |
|---|---:|---:|---:|---:|---:|
| OMP8 | n/a | n/a | n/a | n/a | n/a |
| MPI2xOMP4 | 0.90 | 1.1647 | 1.2287 | 0 | 0 |
| MPI4xOMP2 | 0.93 | 1.1919 | 1.2199 | 0.67 | 0.67 |
| MPI8 | 2.61 | 1.3776 | 1.3492 | 1 | 1 |
| MPI8origin | n/a | n/a | n/a | n/a | n/a |

## zb Performance

Expected iterations: `300`.

| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|


### zb Correctness Metrics

| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |
|---|---:|---:|---:|---:|---:|---:|---:|

### zb Load-Balance Diagnostics

| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |
|---|---:|---:|---:|---:|---:|

## MPI8origin Baseline Comparison

`MPI8origin` is the original decomposed OpenFOAM `-parallel` baseline.  It is compared with both the single-process OpenMP baseline and the raw-MPI replicated-mesh `MPI8` mode for every case where the baseline was run.

| case | comparison | real mean | delta | full evolve mean | delta |
|---|---|---:|---:|---:|---:|
| ourmesh | MPI8origin vs OMP8 | 126.63 vs 63.23 | 63.39 s (+100.25%) | 126.50 vs 57.55 | 68.95 s (+119.80%) |
| ourmesh | MPI8origin vs replicated-mesh MPI8 | 126.63 vs 84.25 | 42.38 s (+50.30%) | 126.50 vs 75.28 | 51.22 s (+68.04%) |

## Per-Run Logs

| case | mode | rep | exit | real | iterations | correctness | log |
|---|---|---:|---:|---:|---:|---|---|
| ourmesh | OMP8 | 1 | 0 | 62.27 | 500 | True | `ourmesh_OMP8_rep1.log` |
| ourmesh | OMP8 | 2 | 0 | 64.58 | 500 | True | `ourmesh_OMP8_rep2.log` |
| ourmesh | OMP8 | 3 | 0 | 62.85 | 500 | True | `ourmesh_OMP8_rep3.log` |
| ourmesh | MPI2xOMP4 | 1 | 0 | 67.77 | 500 | True | `ourmesh_MPI2xOMP4_rep1.log` |
| ourmesh | MPI2xOMP4 | 2 | 0 | 67.15 | 500 | True | `ourmesh_MPI2xOMP4_rep2.log` |
| ourmesh | MPI2xOMP4 | 3 | 0 | 67.37 | 500 | True | `ourmesh_MPI2xOMP4_rep3.log` |
| ourmesh | MPI4xOMP2 | 1 | 0 | 66.12 | 500 | True | `ourmesh_MPI4xOMP2_rep1.log` |
| ourmesh | MPI4xOMP2 | 2 | 0 | 64.95 | 500 | True | `ourmesh_MPI4xOMP2_rep2.log` |
| ourmesh | MPI4xOMP2 | 3 | 0 | 67.20 | 500 | True | `ourmesh_MPI4xOMP2_rep3.log` |
| ourmesh | MPI8 | 1 | 0 | 85.29 | 500 | True | `ourmesh_MPI8_rep1.log` |
| ourmesh | MPI8 | 2 | 0 | 83.87 | 500 | True | `ourmesh_MPI8_rep2.log` |
| ourmesh | MPI8 | 3 | 0 | 83.59 | 500 | True | `ourmesh_MPI8_rep3.log` |
| ourmesh | MPI8origin | 1 | 0 | 125.73 | 500 | True | `ourmesh_MPI8origin_rep1.log` |
| ourmesh | MPI8origin | 2 | 0 | 126.50 | 500 | True | `ourmesh_MPI8origin_rep2.log` |
| ourmesh | MPI8origin | 3 | 0 | 127.65 | 500 | True | `ourmesh_MPI8origin_rep3.log` |
