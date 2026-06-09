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
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_rerun
```

## Correctness Status

All parsed runs passed the correctness gates: exit 0, expected iteration count, `End main`, no fatal/segmentation/MPI abort/NaN markers, nonzero final particles, zero stuck particles, and finite total energy.

## ourmesh Performance

Expected iterations: `500`.

| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 3/3 | 87.37 | 2.31 | 85.05-89.67 | 80.34 | 45.82 | 6.10 | 3.21 | 24.70 |
| MPI2xOMP4 | 3/3 | 91.06 | 0.62 | 90.40-91.64 | 83.43 | 46.03 | 6.70 | 14.20 | 22.70 |
| MPI4xOMP2 | 3/3 | 93.17 | 1.62 | 92.16-95.03 | 85.70 | 45.83 | 5.59 | 23.09 | 22.97 |
| MPI8 | 3/3 | 111.98 | 2.85 | 108.74-114.06 | 101.84 | 54.76 | 6.61 | 40.59 | 18.89 |
| MPI8origin | 3/3 | 147.94 | 0.63 | 147.45-148.65 | 146.11 | 111.13 | 14.52 | 18.22 | 21.97 |

Best mean external `real`: `OMP8` at `87.37 s`.

| mode | mean real delta vs OMP8 | percent |
|---|---:|---:|
| OMP8 | 0.00 s | 0.00% |
| MPI2xOMP4 | 3.69 s | 4.22% |
| MPI4xOMP2 | 5.80 s | 6.63% |
| MPI8 | 24.61 s | 28.17% |
| MPI8origin | 60.57 s | 69.33% |

### ourmesh Correctness Metrics

| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |
|---|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 2463582 | 2463551-2463606 | 35678 | 68771 | 1.2430186450 | 1.2429623320-1.2430972160 | 0 |
| MPI2xOMP4 | 2463822 | 2463794-2463874 | 35512 | 68229 | 1.2428026890 | 1.2424455950-1.2431445910 | 0 |
| MPI4xOMP2 | 2463878 | 2463827-2463927 | 35299 | 64612 | 1.2415806257 | 1.2413554050-1.2419720340 | 0 |
| MPI8 | 2463782 | 2463714-2463822 | 34290 | 60220 | 1.2406224307 | 1.2402538530-1.2409511330 | 0 |
| MPI8origin | 2463662 | 2463614-2463725 | 35677 | 68708 | 1.2432895513 | 1.2432503360-1.2433143520 | 0 |

### ourmesh Load-Balance Diagnostics

| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |
|---|---:|---:|---:|---:|---:|
| OMP8 | n/a | n/a | n/a | n/a | n/a |
| MPI2xOMP4 | 1.80 | 1.2712 | 1.0508 | 3.67 | 3.67 |
| MPI4xOMP2 | 2.27 | 1.2166 | 1.0837 | 3.33 | 3.33 |
| MPI8 | 2.37 | 2.5574 | 1.0698 | 4 | 4 |
| MPI8origin | n/a | n/a | n/a | n/a | n/a |

## zb Performance

Expected iterations: `300`.

| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 3/3 | 72.78 | 1.28 | 71.61-74.14 | 66.87 | 37.38 | 7.58 | 9.40 | 11.37 |
| MPI2xOMP4 | 3/3 | 78.90 | 0.59 | 78.38-79.54 | 72.68 | 41.03 | 10.66 | 14.80 | 11.34 |
| MPI4xOMP2 | 3/3 | 78.63 | 1.39 | 77.02-79.48 | 73.70 | 38.60 | 7.18 | 22.57 | 11.12 |
| MPI8 | 3/3 | 107.31 | 12.65 | 94.34-119.62 | 101.31 | 51.03 | 14.19 | 49.16 | 14.97 |
| MPI8origin | 3/3 | 129.74 | 0.61 | 129.36-130.44 | 130.11 | 83.07 | 19.44 | 25.44 | 14.36 |

Best mean external `real`: `OMP8` at `72.78 s`.

| mode | mean real delta vs OMP8 | percent |
|---|---:|---:|
| OMP8 | 0.00 s | 0.00% |
| MPI2xOMP4 | 6.12 s | 8.41% |
| MPI4xOMP2 | 5.85 s | 8.03% |
| MPI8 | 34.53 s | 47.45% |
| MPI8origin | 56.96 s | 78.26% |

### zb Correctness Metrics

| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |
|---|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 1957964 | 1957935-1957986 | 242244 | 477718 | 0.0019616803 | 0.0019614338-0.0019618318 | 0 |
| MPI2xOMP4 | 1958025 | 1957848-1958210 | 229159 | 445464 | 0.0019587546 | 0.0019582480-0.0019592069 | 0 |
| MPI4xOMP2 | 1959531 | 1959437-1959631 | 225657 | 429367 | 0.0019534190 | 0.0019516008-0.0019566291 | 0 |
| MPI8 | 1958609 | 1958500-1958785 | 221341 | 417574 | 0.0019569224 | 0.0019540346-0.0019585066 | 0 |
| MPI8origin | 1957751 | 1957650-1957864 | 240447 | 474010 | 0.0019622363 | 0.0019619917-0.0019626472 | 0 |

### zb Load-Balance Diagnostics

| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |
|---|---:|---:|---:|---:|---:|
| OMP8 | n/a | n/a | n/a | n/a | n/a |
| MPI2xOMP4 | 2.23 | 1.2515 | 1.0318 | 3.67 | 3.67 |
| MPI4xOMP2 | 1.79 | 1.6457 | 1.0256 | 2.33 | 2.33 |
| MPI8 | 2.45 | 1.7042 | 1.0342 | 2.33 | 2.33 |
| MPI8origin | n/a | n/a | n/a | n/a | n/a |

## MPI8origin Baseline Comparison

`MPI8origin` is the original decomposed OpenFOAM `-parallel` baseline.  It is compared with both the single-process OpenMP baseline and the raw-MPI replicated-mesh `MPI8` mode for every case where the baseline was run.

| case | comparison | real mean | delta | full evolve mean | delta |
|---|---|---:|---:|---:|---:|
| ourmesh | MPI8origin vs OMP8 | 147.94 vs 87.37 | 60.57 s (+69.33%) | 146.11 vs 80.34 | 65.77 s (+81.87%) |
| ourmesh | MPI8origin vs replicated-mesh MPI8 | 147.94 vs 111.98 | 35.96 s (+32.11%) | 146.11 vs 101.84 | 44.27 s (+43.47%) |
| zb | MPI8origin vs OMP8 | 129.74 vs 72.78 | 56.96 s (+78.26%) | 130.11 vs 66.87 | 63.24 s (+94.57%) |
| zb | MPI8origin vs replicated-mesh MPI8 | 129.74 vs 107.31 | 22.42 s (+20.90%) | 130.11 vs 101.31 | 28.80 s (+28.43%) |

## Per-Run Logs

| case | mode | rep | exit | real | iterations | correctness | log |
|---|---|---:|---:|---:|---:|---|---|
| ourmesh | OMP8 | 1 | 0 | 85.05 | 500 | True | `ourmesh_OMP8_rep1.log` |
| ourmesh | OMP8 | 2 | 0 | 89.67 | 500 | True | `ourmesh_OMP8_rep2.log` |
| ourmesh | OMP8 | 3 | 0 | 87.39 | 500 | True | `ourmesh_OMP8_rep3.log` |
| ourmesh | MPI2xOMP4 | 1 | 0 | 91.64 | 500 | True | `ourmesh_MPI2xOMP4_rep1.log` |
| ourmesh | MPI2xOMP4 | 2 | 0 | 90.40 | 500 | True | `ourmesh_MPI2xOMP4_rep2.log` |
| ourmesh | MPI2xOMP4 | 3 | 0 | 91.14 | 500 | True | `ourmesh_MPI2xOMP4_rep3.log` |
| ourmesh | MPI4xOMP2 | 1 | 0 | 92.31 | 500 | True | `ourmesh_MPI4xOMP2_rep1.log` |
| ourmesh | MPI4xOMP2 | 2 | 0 | 92.16 | 500 | True | `ourmesh_MPI4xOMP2_rep2.log` |
| ourmesh | MPI4xOMP2 | 3 | 0 | 95.03 | 500 | True | `ourmesh_MPI4xOMP2_rep3.log` |
| ourmesh | MPI8 | 1 | 0 | 108.74 | 500 | True | `ourmesh_MPI8_rep1.log` |
| ourmesh | MPI8 | 2 | 0 | 113.15 | 500 | True | `ourmesh_MPI8_rep2.log` |
| ourmesh | MPI8 | 3 | 0 | 114.06 | 500 | True | `ourmesh_MPI8_rep3.log` |
| ourmesh | MPI8origin | 1 | 0 | 148.65 | 500 | True | `ourmesh_MPI8origin_rep1.log` |
| ourmesh | MPI8origin | 2 | 0 | 147.45 | 500 | True | `ourmesh_MPI8origin_rep2.log` |
| ourmesh | MPI8origin | 3 | 0 | 147.73 | 500 | True | `ourmesh_MPI8origin_rep3.log` |
| zb | OMP8 | 1 | 0 | 71.61 | 300 | True | `zb_OMP8_rep1.log` |
| zb | OMP8 | 2 | 0 | 72.59 | 300 | True | `zb_OMP8_rep2.log` |
| zb | OMP8 | 3 | 0 | 74.14 | 300 | True | `zb_OMP8_rep3.log` |
| zb | MPI2xOMP4 | 1 | 0 | 78.38 | 300 | True | `zb_MPI2xOMP4_rep1.log` |
| zb | MPI2xOMP4 | 2 | 0 | 79.54 | 300 | True | `zb_MPI2xOMP4_rep2.log` |
| zb | MPI2xOMP4 | 3 | 0 | 78.78 | 300 | True | `zb_MPI2xOMP4_rep3.log` |
| zb | MPI4xOMP2 | 1 | 0 | 79.38 | 300 | True | `zb_MPI4xOMP2_rep1.log` |
| zb | MPI4xOMP2 | 2 | 0 | 77.02 | 300 | True | `zb_MPI4xOMP2_rep2.log` |
| zb | MPI4xOMP2 | 3 | 0 | 79.48 | 300 | True | `zb_MPI4xOMP2_rep3.log` |
| zb | MPI8 | 1 | 0 | 119.62 | 300 | True | `zb_MPI8_rep1.log` |
| zb | MPI8 | 2 | 0 | 107.98 | 300 | True | `zb_MPI8_rep2.log` |
| zb | MPI8 | 3 | 0 | 94.34 | 300 | True | `zb_MPI8_rep3.log` |
| zb | MPI8origin | 1 | 0 | 129.41 | 300 | True | `zb_MPI8origin_rep1.log` |
| zb | MPI8origin | 2 | 0 | 130.44 | 300 | True | `zb_MPI8origin_rep2.log` |
| zb | MPI8origin | 3 | 0 | 129.36 | 300 | True | `zb_MPI8origin_rep3.log` |
