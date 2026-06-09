# zb-cylinder-react mixed parallel 8-core comparison - 2026-06-08

## Scope

Compare four 8-core modes on:

```text
run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react
```

The source case is `omp8`; the three MPI/mixed cases are copied directly from
that directory before editing `system/controlDict`.

| label | run form | case |
|---|---|---|
| OMP8 | `OMP_NUM_THREADS=8 dsmcFoam+` | `zb-cylinder-react/omp8` |
| MPI2xOMP4 | `OMP_NUM_THREADS=4 mpirun -np 2 dsmcFoam+` | `zb-cylinder-react/mpi2omp4` |
| MPI4xOMP2 | `OMP_NUM_THREADS=2 mpirun -np 4 dsmcFoam+` | `zb-cylinder-react/mpi4omp2` |
| MPI8 | `OMP_NUM_THREADS=1 mpirun -np 8 dsmcFoam+` | `zb-cylinder-react/mpi8` |

MPI modes use replicated mesh DLB and do not use OpenFOAM `-parallel`.

## Initial Control Policy

`omp8/system/controlDict` was inspected before case creation. The `zb` formal
step policy is retained for all four cases:

```text
endTime 1.9920146682e-05;
deltaT 6.640048894e-08;
profileSummary true;
profileDetail false;
```

This is a 300-step case by `endTime / deltaT`. The existing `omp8` directory
also contains older stage-26 schedule-sweep logs; those are template history,
not results from this comparison.

For the final timing run, `writeInterval` was set to `1.e-3` in all four
cases.  This keeps the four modes comparable and avoids raw-MPI replicated-mesh
final field/particle gather, which reached the target time in a 10-step
diagnostic but hung during output/redistribution cleanup.

## Source State

- same-tet area reuse remains active;
- barycentric tracker remains inactive;
- mixed boundary safety fix is active: patch boundary `controlParticle()` runs
  under `dsmcMoveBoundary` critical during OpenMP move;
- OpenMP move-detail aggregation instrumentation is active.
- Raw-MPI replicated mesh move fix is active: `Cloud::move()` uses OpenFOAM
  processor-patch transfer only for non-replicated OpenFOAM `-parallel` runs.
- Replicated-mesh initial particle distribution now rebuilds cell occupancy
  before the first move after deleting non-owned parcels.

## Preparation Log

- `20260608_115212`: copied `omp8` to `mpi2omp4`, `mpi4omp2`, and `mpi8`.
- Backed up the original `omp8/system/controlDict` in both the case directory
  and `doc/worklog/v2506/detail_mix`.
- Backed up each copied MPI/mixed case's original `controlDict` before editing.
- Saved the actual comparison `controlDict` files as:

```text
zb_omp8_controlDict_compare_20260608_115212
zb_mpi2omp4_controlDict_compare_20260608_115212
zb_mpi4omp2_controlDict_compare_20260608_115212
zb_mpi8_controlDict_compare_20260608_115212
```

Final no-write comparison controls were then saved as:

```text
zb_omp8_controlDict_compare_rawmpi_guard_nowrite_20260608_115212
zb_mpi2omp4_controlDict_compare_rawmpi_guard_nowrite_20260608_115212
zb_mpi4omp2_controlDict_compare_rawmpi_guard_nowrite_20260608_115212
zb_mpi8_controlDict_compare_rawmpi_guard_nowrite_20260608_115212
```

## Diagnostics Before Final Run

The first replicated-mesh MPI attempts failed after initial particle
distribution:

```text
void Foam::Cloud<Foam::dsmcParcel>::move<...>()
KILLED BY SIGNAL: 11 (Segmentation fault)
```

This reproduced in `MPI2xOMP4`, pure `MPI2`, and `MPI8`, so it was not a mixed
OpenMP-only failure.  A 10-step raw-MPI diagnostic after the `Cloud::move()`
guard fix reached `Time = 6.640048894e-07` and wrote fields, confirming the
first-move segfault was fixed.  That short run then hung during final
replicated-mesh output cleanup, so the formal comparison uses no-write controls.

## Logs

| label | log | exit |
|---|---|---:|
| OMP8 | `zb_compare_omp8_rawmpi_guard_nowrite_300step_20260608_115212.log` | 0 |
| MPI2xOMP4 | `zb_compare_mpi2omp4_replicated_rawmpi_guard_nowrite_300step_20260608_115212.log` | 0 |
| MPI4xOMP2 | `zb_compare_mpi4omp2_replicated_rawmpi_guard_nowrite_300step_20260608_115212.log` | 0 |
| MPI8 | `zb_compare_mpi8_replicated_rawmpi_guard_nowrite_300step_20260608_115212.log` | 0 |

No `FOAM FATAL`, `Fatal`, `Segmentation`, `SIGSEGV`, `MPI_ABORT`,
`BAD TERMINATION`, `Killed`, `nan`, or `NaN` marker was found in the four final
logs.

## Results

### Performance

| mode | real | full evolve | move+collide | move | buildOcc | collision | post | DLB rebalances |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 73.94 s | 69.05171913 s | 57.11835853 s | 38.42851716 s | 7.84149354 s | 9.697659374 s | 11.93026923 s | n/a |
| MPI2xOMP4 | 81.46 s | 78.55469086 s | 68.63900701 s | 43.38214653 s | 11.19479529 s | 15.93265544 s | 12.58936173 s | 4 |
| MPI4xOMP2 | 83.55 s | 81.26013632 s | 70.63147743 s | 40.08826542 s | 7.162368006 s | 26.48233093 s | 11.67984341 s | 2 |
| MPI8 | 116.09 s | 113.9744001 s | 99.90250354 s | 54.77817548 s | 16.43615009 s | 53.16270512 s | 17.551509 s | 2 |

### Correctness/Sanity

| mode | iterations | particles | stuck | collisions | candidates | total energy |
|---|---:|---:|---:|---:|---:|---:|
| OMP8 | 300 | 1957895 | 0 | 242643 | 478681 | 0.001961240578 |
| MPI2xOMP4 | 300 | 1957590 | 0 | 236438 | 454106 | 0.001960220523 |
| MPI4xOMP2 | 300 | 1959669 | 0 | 219624 | 412494 | 0.001952720265 |
| MPI8 | 300 | 1958433 | 0 | 213314 | 374248 | 0.001958114059 |

All four formal runs reached `Total Iterations = 300` and `End main`.
Final particle counts are within `2079` parcels across the four modes, and all
final stuck-particle counts are zero.

### Replicated Mesh DLB

| mode | migration calls | migration wall | particles max/min | rank wall max/min | DLB checks | triggered checks |
|---|---:|---:|---:|---:|---:|---:|
| MPI2xOMP4 | 35 | 2.221775225 s | 1.234741946 | 1.03590319 | 300 | 4 |
| MPI4xOMP2 | 33 | 2.278614359 s | 1.41406743 | 1.016788928 | 300 | 2 |
| MPI8 | 33 | 4.131586484 s | 1.607938539 | 1.029983705 | 300 | 2 |

## Interpretation

Single-run ranking by external `real`:

```text
OMP8 (73.94 s) < MPI2xOMP4 (81.46 s) < MPI4xOMP2 (83.55 s) < MPI8 (116.09 s)
```

Relative to OMP8:

| mode | real delta | full-evolve delta |
|---|---:|---:|
| MPI2xOMP4 | +7.52 s, +10.17% | +9.50297173 s, +13.76% |
| MPI4xOMP2 | +9.61 s, +13.00% | +12.20841719 s, +17.68% |
| MPI8 | +42.15 s, +57.01% | +44.92268097 s, +65.06% |

OMP8 remains the best single-run result for this `zb` case.  The mixed MPI+OMP
runs are close to OMP8, especially `MPI2xOMP4`, but replicated-mesh migration
and higher collision/build-occupancy costs offset the smaller rank-local working
sets.  Pure `MPI8` is much slower in this run because collision time grows to
`53.16270512 s` and move time also rises to `54.77817548 s`.

Between the two mixed splits, `MPI2xOMP4` is slightly faster in both `real` and
`full evolve`.  `MPI4xOMP2` has lower move and build-occupancy time, but its
collision time is much higher.  A repeated run would be needed before treating
the `2x4` versus `4x2` gap as stable, but this first formal `zb` run favors
`MPI2xOMP4`.

## Write-Fix Follow-up

After the no-write comparison above, the raw-MPI replicated-mesh output path was
continued and re-tested.

### Source State

- `Cloud::move()` still skips OpenFOAM processor-patch transfer in replicated
  raw-MPI mode.
- Solver output now gathers replicated-mesh parcels to rank0, temporarily
  disables the cloud's automatic `runTime.write()` path, writes non-cloud
  registered objects, then writes the gathered cloud and `cellOwner`.
- `writeGatheredCloudOnRank0()` now writes `cloudProperties` entries for every
  MPI rank, using the same conservative `particleCount` value for raw-MPI
  restart compatibility.
- The attempted full binary stream gather was rejected: it exited cleanly but
  wrote only the rank0-local parcel subset (`548696/1330276` for `np=2`,
  `203649/1330272` for `np=4`, `265748/1330322` for `np=8`).  The final path
  therefore keeps the already validated flat POD gather for output.

### Build

| log | exit |
|---|---:|
| `stage_write_flat_gather_cloudprops_build_20260608.log` | 0 |

### Final Write Diagnostics

All diagnostics use the same `zb-cylinder-react` 10-step output point:

```text
endTime 6.640048894e-07;
deltaT 6.640048894e-08;
writeInterval 6.640048894e-07;
```

| mode | log | exit | particles | gathered cloud | gather | rank0 write | migrate-back |
|---|---|---:|---:|---:|---:|---:|---:|
| MPI4xOMP2 | `zb_diag_mpi4omp2_flat_gather_cloudprops_write_10step_20260608.log` | 0 | 1330313 | 1330313 | 0.09 s | 3.07 s | 0.15 s |
| MPI8 | `zb_diag_mpi8_flat_gather_cloudprops_write_10step_20260608.log` | 0 | 1330305 | 1330305 | 0.07 s | 4.11 s | 0.14 s |

No `FOAM FATAL`, `Fatal`, `Segmentation`, `SIGSEGV`, `MPI_ABORT`,
`BAD TERMINATION`, `Killed`, `nan`, or `NaN` marker was found in either final
write diagnostic.  Both reached `Total Iterations = 10` and `End main`.

The final `mpi2omp4` output `cloudProperties` file includes `processor0` through
`processor7`, each with the same `particleCount 2395877`, so raw-MPI restarts no
longer depend on only `processor0` being present.

The temporary `mpi4omp2/system/controlDict` 10-step write edit was restored from
`zb_mpi4omp2_controlDict_before_write_current_10step_20260608`.
