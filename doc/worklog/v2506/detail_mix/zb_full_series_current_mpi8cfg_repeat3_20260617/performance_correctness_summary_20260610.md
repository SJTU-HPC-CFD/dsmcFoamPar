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
doc/worklog/v2506/detail_mix/zb_full_series_current_mpi8cfg_repeat3_20260617
```

## Correctness Status

All parsed runs passed the correctness gates: exit 0, expected iteration count, `End main`, no fatal/segmentation/MPI abort/NaN markers, nonzero final particles, zero stuck particles, and finite total energy.

## ourmesh Performance

Expected iterations: `500`.

| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|


### ourmesh Correctness Metrics

| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |
|---|---:|---:|---:|---:|---:|---:|---:|

### ourmesh Load-Balance Diagnostics

| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |
|---|---:|---:|---:|---:|---:|

## zb Performance

Expected iterations: `300`.

| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 3/3 | 55.41 | 0.33 | 55.07-55.73 | 50.85 | 37.06 | 3.92 | 8.68 | 0.17 |
| MPI2xOMP4 | 3/3 | 59.26 | 0.53 | 58.76-59.82 | 53.35 | 39.78 | 4.79 | 8.54 | 0.18 |
| MPI4xOMP2 | 3/3 | 60.48 | 0.96 | 59.48-61.39 | 54.02 | 40.46 | 3.83 | 9.61 | 0.19 |
| MPI8 | 3/3 | 69.32 | 2.26 | 66.82-71.22 | 60.67 | 46.39 | 4.98 | 10.75 | 0.19 |
| MPI8origin | 3/3 | 107.68 | 0.60 | 107.12-108.32 | 107.65 | 79.64 | 9.87 | 23.71 | 0.12 |

Best mean external `real`: `OMP8` at `55.41 s`.

| mode | mean real delta vs OMP8 | percent |
|---|---:|---:|
| OMP8 | 0.00 s | 0.00% |
| MPI2xOMP4 | 3.85 s | 6.95% |
| MPI4xOMP2 | 5.07 s | 9.16% |
| MPI8 | 13.91 s | 25.10% |
| MPI8origin | 52.27 s | 94.34% |

### zb Correctness Metrics

| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |
|---|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 1957932 | 1957750-1958123 | 242505 | 477999 | 0.0019611839 | 0.0019609183-0.0019614311 | 0 |
| MPI2xOMP4 | 1957737 | 1957530-1958110 | 84781 | 165642 | 0.0019609948 | 0.0019607680-0.0019612134 | 0 |
| MPI4xOMP2 | 1959961 | 1959699-1960120 | 7669 | 11483 | 0.0019600020 | 0.0019573627-0.0019617662 | 0 |
| MPI8 | 1959896 | 1959652-1960034 | 29252 | 53232 | 0.0019550757 | 0.0019523858-0.0019568601 | 0 |
| MPI8origin | 1958092 | 1957861-1958263 | 240806 | 474559 | 0.0019624003 | 0.0019622055-0.0019626217 | 0 |

### zb Load-Balance Diagnostics

| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |
|---|---:|---:|---:|---:|---:|
| OMP8 | n/a | n/a | n/a | n/a | n/a |
| MPI2xOMP4 | 2.08 | 1.2559 | 1.0643 | 1 | 1 |
| MPI4xOMP2 | 2.31 | 1.1721 | 1.2576 | 3 | 3 |
| MPI8 | 2.87 | 1.2299 | 1.3732 | 5 | 5 |
| MPI8origin | n/a | n/a | n/a | n/a | n/a |

## MPI8origin Baseline Comparison

`MPI8origin` is the original decomposed OpenFOAM `-parallel` baseline.  It is compared with both the single-process OpenMP baseline and the raw-MPI replicated-mesh `MPI8` mode for every case where the baseline was run.

| case | comparison | real mean | delta | full evolve mean | delta |
|---|---|---:|---:|---:|---:|
| zb | MPI8origin vs OMP8 | 107.68 vs 55.41 | 52.27 s (+94.34%) | 107.65 vs 50.85 | 56.80 s (+111.71%) |
| zb | MPI8origin vs replicated-mesh MPI8 | 107.68 vs 69.32 | 38.37 s (+55.35%) | 107.65 vs 60.67 | 46.98 s (+77.45%) |

## Per-Run Logs

| case | mode | rep | exit | real | iterations | correctness | log |
|---|---|---:|---:|---:|---:|---|---|
| zb | OMP8 | 1 | 0 | 55.43 | 300 | True | `zb_OMP8_rep1.log` |
| zb | OMP8 | 2 | 0 | 55.07 | 300 | True | `zb_OMP8_rep2.log` |
| zb | OMP8 | 3 | 0 | 55.73 | 300 | True | `zb_OMP8_rep3.log` |
| zb | MPI2xOMP4 | 1 | 0 | 58.76 | 300 | True | `zb_MPI2xOMP4_rep1.log` |
| zb | MPI2xOMP4 | 2 | 0 | 59.21 | 300 | True | `zb_MPI2xOMP4_rep2.log` |
| zb | MPI2xOMP4 | 3 | 0 | 59.82 | 300 | True | `zb_MPI2xOMP4_rep3.log` |
| zb | MPI4xOMP2 | 1 | 0 | 60.58 | 300 | True | `zb_MPI4xOMP2_rep1.log` |
| zb | MPI4xOMP2 | 2 | 0 | 61.39 | 300 | True | `zb_MPI4xOMP2_rep2.log` |
| zb | MPI4xOMP2 | 3 | 0 | 59.48 | 300 | True | `zb_MPI4xOMP2_rep3.log` |
| zb | MPI8 | 1 | 0 | 66.82 | 300 | True | `zb_MPI8_rep1.log` |
| zb | MPI8 | 2 | 0 | 71.22 | 300 | True | `zb_MPI8_rep2.log` |
| zb | MPI8 | 3 | 0 | 69.91 | 300 | True | `zb_MPI8_rep3.log` |
| zb | MPI8origin | 1 | 0 | 107.61 | 300 | True | `zb_MPI8origin_rep1.log` |
| zb | MPI8origin | 2 | 0 | 108.32 | 300 | True | `zb_MPI8origin_rep2.log` |
| zb | MPI8origin | 3 | 0 | 107.12 | 300 | True | `zb_MPI8origin_rep3.log` |
