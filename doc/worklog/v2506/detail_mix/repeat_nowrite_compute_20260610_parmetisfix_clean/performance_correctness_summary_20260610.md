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
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_parmetisfix_clean
```

## Correctness Status

All parsed runs passed the correctness gates: exit 0, expected iteration count, `End main`, no fatal/segmentation/MPI abort/NaN markers, nonzero final particles, zero stuck particles, and finite total energy.

## ourmesh Performance

Expected iterations: `500`.

| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 3/3 | 62.22 | 1.60 | 61.10-64.06 | 60.46 | 49.11 | 6.62 | 3.52 | 0.65 |
| MPI2xOMP4 | 3/3 | 70.30 | 1.48 | 68.60-71.18 | 68.12 | 54.41 | 7.54 | 16.87 | 1.07 |
| MPI4xOMP2 | 3/3 | 70.10 | 1.49 | 68.42-71.24 | 67.91 | 51.95 | 6.84 | 12.67 | 2.20 |
| MPI8 | 3/3 | 101.45 | 1.43 | 100.01-102.87 | 98.63 | 69.07 | 10.61 | 32.80 | 5.00 |
| MPI8origin | 3/3 | 135.07 | 2.31 | 132.44-136.78 | 144.27 | 117.47 | 16.31 | 20.56 | 0.48 |

Best mean external `real`: `OMP8` at `62.22 s`.

| mode | mean real delta vs OMP8 | percent |
|---|---:|---:|
| OMP8 | 0.00 s | 0.00% |
| MPI2xOMP4 | 8.08 s | 12.99% |
| MPI4xOMP2 | 7.88 s | 12.66% |
| MPI8 | 39.22 s | 63.04% |
| MPI8origin | 72.85 s | 117.07% |

### ourmesh Correctness Metrics

| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |
|---|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 2463611 | 2463545-2463681 | 35591 | 68786 | 1.2429850577 | 1.2428346060-1.2430942860 | 0 |
| MPI2xOMP4 | 2463726 | 2463659-2463786 | 35550 | 68586 | 1.2431112053 | 1.2430114520-1.2431917970 | 0 |
| MPI4xOMP2 | 2463860 | 2463788-2463938 | 35614 | 68353 | 1.2425777920 | 1.2423908420-1.2428044570 | 0 |
| MPI8 | 2463796 | 2463776-2463829 | 35101 | 67264 | 1.2421654850 | 1.2417871210-1.2424040990 | 0 |
| MPI8origin | 2463700 | 2463640-2463765 | 35517 | 68544 | 1.2432590103 | 1.2432013240-1.2433580940 | 0 |

### ourmesh Load-Balance Diagnostics

| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |
|---|---:|---:|---:|---:|---:|
| OMP8 | n/a | n/a | n/a | n/a | n/a |
| MPI2xOMP4 | 1.05 | 1.1623 | 1.0005 | 2.67 | 2.67 |
| MPI4xOMP2 | 1.53 | 1.2190 | 1.0018 | 2.33 | 2.33 |
| MPI8 | 2.96 | 1.3181 | 1.0030 | 2.67 | 2.67 |
| MPI8origin | n/a | n/a | n/a | n/a | n/a |

## zb Performance

Expected iterations: `300`.

| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 3/3 | 58.74 | 1.25 | 57.31-59.61 | 57.76 | 41.69 | 4.26 | 10.46 | 0.18 |
| MPI2xOMP4 | 3/3 | 63.66 | 1.22 | 62.49-64.92 | 63.59 | 42.89 | 5.11 | 13.43 | 0.27 |
| MPI4xOMP2 | 3/3 | 67.72 | 0.17 | 67.58-67.91 | 67.42 | 44.66 | 3.73 | 24.39 | 0.55 |
| MPI8 | 3/3 | 79.81 | 1.59 | 78.20-81.38 | 81.33 | 54.01 | 6.51 | 33.73 | 1.29 |
| MPI8origin | 3/3 | 114.80 | 1.20 | 113.42-115.60 | 124.21 | 90.73 | 11.23 | 28.69 | 0.13 |

Best mean external `real`: `OMP8` at `58.74 s`.

| mode | mean real delta vs OMP8 | percent |
|---|---:|---:|
| OMP8 | 0.00 s | 0.00% |
| MPI2xOMP4 | 4.92 s | 8.38% |
| MPI4xOMP2 | 8.98 s | 15.28% |
| MPI8 | 21.07 s | 35.86% |
| MPI8origin | 56.06 s | 95.44% |

### zb Correctness Metrics

| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |
|---|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 1957871 | 1957484-1958390 | 241823 | 477015 | 0.0019613657 | 0.0019609871-0.0019620457 | 0 |
| MPI2xOMP4 | 1958166 | 1958061-1958320 | 226886 | 441423 | 0.0019570060 | 0.0019510445-0.0019613185 | 0 |
| MPI4xOMP2 | 1959629 | 1959318-1960011 | 222185 | 429424 | 0.0019459436 | 0.0019434562-0.0019474198 | 0 |
| MPI8 | 1958375 | 1957998-1958771 | 206645 | 391782 | 0.0019360628 | 0.0019326233-0.0019393615 | 0 |
| MPI8origin | 1958023 | 1957683-1958417 | 240788 | 475115 | 0.0019624559 | 0.0019622125-0.0019627143 | 0 |

### zb Load-Balance Diagnostics

| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |
|---|---:|---:|---:|---:|---:|
| OMP8 | n/a | n/a | n/a | n/a | n/a |
| MPI2xOMP4 | 2.06 | 1.2048 | 1.0008 | 2.33 | 2.33 |
| MPI4xOMP2 | 2.47 | 1.2532 | 1.0013 | 2 | 2 |
| MPI8 | 2.39 | 1.3424 | 1.0020 | 2.33 | 2.33 |
| MPI8origin | n/a | n/a | n/a | n/a | n/a |

## MPI8origin Baseline Comparison

`MPI8origin` is the original decomposed OpenFOAM `-parallel` baseline.  It is compared with both the single-process OpenMP baseline and the raw-MPI replicated-mesh `MPI8` mode for every case where the baseline was run.

| case | comparison | real mean | delta | full evolve mean | delta |
|---|---|---:|---:|---:|---:|
| ourmesh | MPI8origin vs OMP8 | 135.07 vs 62.22 | 72.85 s (+117.07%) | 144.27 vs 60.46 | 83.81 s (+138.63%) |
| ourmesh | MPI8origin vs replicated-mesh MPI8 | 135.07 vs 101.45 | 33.62 s (+33.14%) | 144.27 vs 98.63 | 45.64 s (+46.28%) |
| zb | MPI8origin vs OMP8 | 114.80 vs 58.74 | 56.06 s (+95.44%) | 124.21 vs 57.76 | 66.45 s (+115.04%) |
| zb | MPI8origin vs replicated-mesh MPI8 | 114.80 vs 79.81 | 34.99 s (+43.85%) | 124.21 vs 81.33 | 42.88 s (+52.72%) |

## Per-Run Logs

| case | mode | rep | exit | real | iterations | correctness | log |
|---|---|---:|---:|---:|---:|---|---|
| ourmesh | OMP8 | 1 | 0 | 61.10 | 500 | True | `ourmesh_OMP8_rep1.log` |
| ourmesh | OMP8 | 2 | 0 | 61.51 | 500 | True | `ourmesh_OMP8_rep2.log` |
| ourmesh | OMP8 | 3 | 0 | 64.06 | 500 | True | `ourmesh_OMP8_rep3.log` |
| ourmesh | MPI2xOMP4 | 1 | 0 | 71.13 | 500 | True | `ourmesh_MPI2xOMP4_rep1.log` |
| ourmesh | MPI2xOMP4 | 2 | 0 | 71.18 | 500 | True | `ourmesh_MPI2xOMP4_rep2.log` |
| ourmesh | MPI2xOMP4 | 3 | 0 | 68.60 | 500 | True | `ourmesh_MPI2xOMP4_rep3.log` |
| ourmesh | MPI4xOMP2 | 1 | 0 | 70.64 | 500 | True | `ourmesh_MPI4xOMP2_rep1.log` |
| ourmesh | MPI4xOMP2 | 2 | 0 | 68.42 | 500 | True | `ourmesh_MPI4xOMP2_rep2.log` |
| ourmesh | MPI4xOMP2 | 3 | 0 | 71.24 | 500 | True | `ourmesh_MPI4xOMP2_rep3.log` |
| ourmesh | MPI8 | 1 | 0 | 100.01 | 500 | True | `ourmesh_MPI8_rep1.log` |
| ourmesh | MPI8 | 2 | 0 | 101.46 | 500 | True | `ourmesh_MPI8_rep2.log` |
| ourmesh | MPI8 | 3 | 0 | 102.87 | 500 | True | `ourmesh_MPI8_rep3.log` |
| ourmesh | MPI8origin | 1 | 0 | 136.78 | 500 | True | `ourmesh_MPI8origin_rep1.log` |
| ourmesh | MPI8origin | 2 | 0 | 132.44 | 500 | True | `ourmesh_MPI8origin_rep2.log` |
| ourmesh | MPI8origin | 3 | 0 | 135.99 | 500 | True | `ourmesh_MPI8origin_rep3.log` |
| zb | OMP8 | 1 | 0 | 57.31 | 300 | True | `zb_OMP8_rep1.log` |
| zb | OMP8 | 2 | 0 | 59.61 | 300 | True | `zb_OMP8_rep2.log` |
| zb | OMP8 | 3 | 0 | 59.30 | 300 | True | `zb_OMP8_rep3.log` |
| zb | MPI2xOMP4 | 1 | 0 | 63.57 | 300 | True | `zb_MPI2xOMP4_rep1.log` |
| zb | MPI2xOMP4 | 2 | 0 | 64.92 | 300 | True | `zb_MPI2xOMP4_rep2.log` |
| zb | MPI2xOMP4 | 3 | 0 | 62.49 | 300 | True | `zb_MPI2xOMP4_rep3.log` |
| zb | MPI4xOMP2 | 1 | 0 | 67.66 | 300 | True | `zb_MPI4xOMP2_rep1.log` |
| zb | MPI4xOMP2 | 2 | 0 | 67.58 | 300 | True | `zb_MPI4xOMP2_rep2.log` |
| zb | MPI4xOMP2 | 3 | 0 | 67.91 | 300 | True | `zb_MPI4xOMP2_rep3.log` |
| zb | MPI8 | 1 | 0 | 78.20 | 300 | True | `zb_MPI8_rep1.log` |
| zb | MPI8 | 2 | 0 | 81.38 | 300 | True | `zb_MPI8_rep2.log` |
| zb | MPI8 | 3 | 0 | 79.84 | 300 | True | `zb_MPI8_rep3.log` |
| zb | MPI8origin | 1 | 0 | 115.38 | 300 | True | `zb_MPI8origin_rep1.log` |
| zb | MPI8origin | 2 | 0 | 113.42 | 300 | True | `zb_MPI8origin_rep2.log` |
| zb | MPI8origin | 3 | 0 | 115.60 | 300 | True | `zb_MPI8origin_rep3.log` |
