# MPI+OpenMP mixed-parallel report pointer - 2026-06-10

This detail-level note is superseded by the complete top-level technical
report:

```text
doc/worklog/v2506/dsmcFoam_plus_mpi_omp_mixed_technical_report_20260610.md
```

The final no-write repeat dataset now contains:

```text
2 cases x 5 modes x 3 repeats = 30 runs
```

The complete parsed tables are in:

```text
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_rerun/performance_correctness_summary_20260610.md
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_rerun/results.csv
doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_rerun/aggregate.csv
```

The previous 27-run version lacked `zb-cylinder-react/MPI8origin`.  That case
was created from `zb-cylinder-react/omp8`, decomposed with `decomposePar` into
8 subdomains, and then added as three no-write `mpirun -np 8 dsmcFoam+ -parallel`
repeats.
