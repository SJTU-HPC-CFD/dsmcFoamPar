# Mixed parallel 8-core comparison - 2026-06-08

## Scope

Compare four 500-step modes on the `ourmesh` cylinder-react case with the
current source tree:

| label | run form | case |
|---|---|---|
| OMP8 | `OMP_NUM_THREADS=8 dsmcFoam+` | `ourmesh/omp8` |
| MPI2xOMP4 | `OMP_NUM_THREADS=4 mpirun -np 2 dsmcFoam+` | `ourmesh/mix-mpi4omp2` |
| MPI4xOMP2 | `OMP_NUM_THREADS=2 mpirun -np 4 dsmcFoam+` | `ourmesh/mix-mpi4omp2` |
| MPI8 | `OMP_NUM_THREADS=1 mpirun -np 8 dsmcFoam+` | `ourmesh/mpi8replicatedmesh` |

MPI modes use replicated mesh DLB and do not use OpenFOAM `-parallel`.

Comparison timestamp:

```text
20260608_111426
```

## Source State

- same-tet area reuse remains active;
- barycentric tracker remains inactive;
- mixed boundary safety fix is active: patch boundary `controlParticle()` runs
  under `dsmcMoveBoundary` critical during OpenMP move;
- OpenMP move-detail aggregation instrumentation is active.

## Control Policy

Formal controls:

```text
endTime 5.e-05;
deltaT 1.e-07;
profileSummary true;
profileDetail false;
```

OMP8 temporarily uses the OMP/profile controlDict from
`ourmeshbkp/omp8/system/controlDict`, then restores its original
`ourmesh/omp8/system/controlDict`.

MPI2xOMP4 and MPI4xOMP2 share `ourmesh/mix-mpi4omp2`, with only
`openmpThreads` changed between the two runs.

MPI8 uses `ourmesh/mpi8replicatedmesh`.

## Logs

| label | log | exit |
|---|---|---:|
| OMP8 | `compare_omp8_500step_20260608_111426.log` | 0 |
| MPI2xOMP4 | `compare_mpi2omp4_replicated_500step_20260608_111426.log` | 0 |
| MPI4xOMP2 | `compare_mpi4omp2_replicated_500step_20260608_111426.log` | 0 |
| MPI8 | `compare_mpi8_replicated_500step_20260608_111426.log` | 0 |

All logs are under:

```text
doc/worklog/v2506/detail_mix
```

No `FOAM FATAL`, `Fatal`, `Segmentation`, `SIGSEGV`, `MPI_ABORT`,
`BAD TERMINATION`, `Killed`, `nan`, or `NaN` marker was found in the four logs.

## Results

### Performance

| mode | real | full evolve | move | buildOcc | collision | post | DLB rebalances |
|---|---:|---:|---:|---:|---:|---:|---:|
| OMP8 | 90.43 s | 86.13389632 s | 48.62578607 s | 6.476924612 s | 3.417942971 s | 27.08038795 s | n/a |
| MPI2xOMP4 | 110.16 s | 106.2217127 s | 69.75613658 s | 8.12312802 s | 4.594994329 s | 24.40489679 s | 3 |
| MPI4xOMP2 | 110.37 s | 106.6533562 s | 73.23105661 s | 7.712518004 s | 5.472380215 s | 26.18601688 s | 4 |
| MPI8 | 116.54 s | 108.9185741 s | 71.27684936 s | 8.088753547 s | 15.69056538 s | 19.52255321 s | 4 |

### Correctness/Sanity

| mode | iterations | particles | stuck | collisions | candidates | total energy |
|---|---:|---:|---:|---:|---:|---:|
| OMP8 | 500 | 2463775 | 0 | 35796 | 68628 | 1.243037974 |
| MPI2xOMP4 | 500 | 2463821 | 0 | 35550 | 68066 | 1.242551047 |
| MPI4xOMP2 | 500 | 2463867 | 0 | 35081 | 65010 | 1.241567594 |
| MPI8 | 500 | 2463813 | 0 | 34999 | 64832 | 1.240300762 |

Particle counts remain within `92` parcels across the four final states, and all
runs completed 500/500 iterations with zero stuck particles.

### Replicated Mesh DLB

| mode | migration calls | migration wall | particles max/min | rank wall max/min | DLB checks | triggered checks |
|---|---:|---:|---:|---:|---:|---:|
| MPI2xOMP4 | 54 | 1.488399117 s | 1.189804225 | 1.030502566 | 500 | 3 |
| MPI4xOMP2 | 55 | 1.703505911 s | 1.447793077 | 1.080703366 | 500 | 4 |
| MPI8 | 55 | 1.735631219 s | 1.728034676 | 1.059783885 | 500 | 4 |

## Interpretation

Single-run ranking by external `real`:

```text
OMP8 (90.43 s) < MPI2xOMP4 (110.16 s) ~= MPI4xOMP2 (110.37 s) < MPI8 (116.54 s)
```

Relative to OMP8:

| mode | real delta | full-evolve delta |
|---|---:|---:|
| MPI2xOMP4 | +19.73 s, +21.82% | +20.08781638 s, +23.32% |
| MPI4xOMP2 | +19.94 s, +22.05% | +20.51945988 s, +23.82% |
| MPI8 | +26.11 s, +28.88% | +22.78467778 s, +26.45% |

The main reason OMP8 wins is move time:

- OMP8 move is `48.62578607 s`;
- replicated MPI/mixed move remains `69.75613658-73.23105661 s`;
- mixed MPI+OMP reduces collision cost versus MPI8, but that is not enough to
  offset the move gap and replicated-mesh overhead.

MPI2xOMP4 and MPI4xOMP2 are effectively tied in this single run.  MPI2xOMP4 is
slightly faster in `real` and `full evolve`, while MPI4xOMP2 has slightly lower
`buildCellOccupancy`.  This difference is small enough that a 2-3 run repeat is
needed before treating one mixed split as definitively better.

MPI8 is slowest in this current auto-DLB comparison because collision time is
much larger (`15.69056538 s`) than in the mixed runs (`4.59-5.47 s`), while its
move time remains in the same high band.

## Restore

After the comparison, the touched control files were restored and copied to:

```text
omp8_controlDict_after_compare_restore_20260608_111426
mix_mpi4omp2_controlDict_after_compare_restore_20260608_111426
mpi8replicatedmesh_controlDict_after_compare_restore_20260608_111426
```
