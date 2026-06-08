# Stage21 Stability Rerun Results - 2026-06-06

## Scope

This note records three 500-step formal reruns of the current kept best OMP
configuration: stage21 `Cloud::move()` no-delete fast path.

No source change was made for this stability check.

Case:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8`

Logs:

- `log.codex_omp_stage21_stability_run1_bkpctrl_500step_20260606`
- `log.codex_omp_stage21_stability_run2_bkpctrl_500step_20260606`
- `log.codex_omp_stage21_stability_run3_bkpctrl_500step_20260606`

After the reruns, `ourmesh/omp8/system/controlDict` was verified as restored:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d  ourmesh/omp8/system/controlDict
```

## Per-Run Results

All three runs reached `End main` and `Total Iterations = 500`.

| metric | run1 | run2 | run3 |
|---|---:|---:|---:|
| real [s] | 99.84 | 96.01 | 99.30 |
| full evolve wall [s] | 99.04519771 | 96.50163862 | 98.96211664 |
| move only [s] | 54.36522341 | 52.42366711 | 54.42676369 |
| buildCellOccupancy [s] | 8.529617478 | 8.412417962 | 8.283163669 |
| collision phase [s] | 3.620491183 | 3.546368092 | 3.659187166 |
| post fields/output [s] | 31.95483556 | 31.55022843 | 32.00615842 |
| DSMC particles | 2463620 | 2463629 | 2463727 |
| stuck particles | 0 | 0 | 0 |
| collisions | 35726 | 35748 | 35678 |
| collision candidates | 68636 | 69204 | 68986 |
| total energy | 1.243124401 | 1.243145619 | 1.24312808 |

No `FOAM FATAL`, `NaN`, `Segmentation`, `MPI_ABORT`, or floating-point error
marker was found in the three logs.

## Stability Statistics

Statistics use the three reruns only.  `std` is sample standard deviation.

| metric | mean | min | max | range | std | CV |
|---|---:|---:|---:|---:|---:|---:|
| real [s] | 98.38333333 | 96.01 | 99.84 | 3.83 | 2.073025165 | 2.107090% |
| full evolve wall [s] | 98.16965099 | 96.50163862 | 99.04519771 | 2.54355909 | 1.445138251 | 1.472082% |
| move only [s] | 53.73855140 | 52.42366711 | 54.42676369 | 2.00309658 | 1.139138855 | 2.119780% |
| buildCellOccupancy [s] | 8.408399703 | 8.283163669 | 8.529617478 | 0.246453809 | 0.123276031 | 1.466106% |
| collision phase [s] | 3.608682147 | 3.546368092 | 3.659187166 | 0.112819074 | 0.057329101 | 1.588644% |
| post fields/output [s] | 31.83707414 | 31.55022843 | 32.00615842 | 0.45592999 | 0.249737569 | 0.784424% |
| DSMC particles | 2463658.67 | 2463620 | 2463727 | 107 | 59.34924880 | 0.002409% |
| collisions | 35717.33 | 35678 | 35748 | 70 | 35.79571669 | 0.100219% |
| collision candidates | 68942.00 | 68636 | 69204 | 568 | 286.5449354 | 0.415632% |
| total energy | 1.243132700 | 1.243124401 | 1.243145619 | 0.000021218 | 0.000011338 | 0.000912% |

## Comparison With First Stage21 Formal Run

The first accepted stage21 500-step formal result was:

| metric | first stage21 formal | stability rerun mean | mean - first |
|---|---:|---:|---:|
| real [s] | 100.65 | 98.38333333 | -2.26666667 |
| full evolve wall [s] | 93.28696662 | 98.16965099 | +4.88268437 |
| move only [s] | 51.17943521 | 53.73855140 | +2.55911619 |
| buildCellOccupancy [s] | 8.087915571 | 8.408399703 | +0.320484132 |
| collision phase [s] | 3.314541395 | 3.608682147 | +0.294140752 |
| post fields/output [s] | 30.15687469 | 31.83707414 | +1.68019945 |
| DSMC particles | 2463537 | 2463658.67 | +121.67 |
| collisions | 36005 | 35717.33 | -287.67 |
| collision candidates | 69029 | 68942.00 | -87.00 |
| total energy | 1.243069478 | 1.243132700 | +0.000063222 |

The repeated runs confirm correctness stability.  Particle count, collision
count, candidate count, and total energy remain in a narrow stochastic band, all
runs have zero stuck particles, and all runs complete the same 500 iterations.

The performance signal is stable enough to keep stage21 as the current best
configuration, but the first accepted stage21 run was faster in the profiled
`full evolve wall` path than the three reruns.  The rerun spread is mostly in
`move only`; `post fields/output` is tighter.  This should be treated as
run-to-run/environment noise unless a later controlled baseline shows a
repeatable regression.

## Conclusion

Stage21 remains the current kept best configuration after three formal 500-step
stability reruns.

Use the rerun mean as the stability envelope for later comparisons:

- `real = 98.38333333 s` mean, range `96.01-99.84 s`;
- `full evolve wall = 98.16965099 s` mean, range `96.50163862-99.04519771 s`;
- `move only = 53.73855140 s` mean, range `52.42366711-54.42676369 s`;
- `buildCellOccupancy = 8.408399703 s` mean, range
  `8.283163669-8.529617478 s`;
- `post fields/output = 31.83707414 s` mean, range
  `31.55022843-32.00615842 s`.
