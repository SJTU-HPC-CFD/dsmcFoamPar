# MPI replicated-mesh scope correction to pure MPI - 2026-06-07

## Correction

The previous OpenMP2 tracker atomic run was out of scope for this round.

The user clarified that this round must use only:

```text
MPI replicated-mesh DLB
useOpenMP false
openmpThreads 1
mpirun -np 8
```

The reference benchmark for this round is the pure MPI8 log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610
```

Known reference metrics from the existing worklog series:

| metric | pure MPI8 reference |
| --- | ---: |
| external wall | 117.31 s |
| move only | 47.1652636 s |
| buildCellOccupancy | 10.7084157 s |
| collision phase | 10.89824478 s |
| post fields/output | 30.35563996 s |
| full evolve wall | 114.9207628 s |
| DLB rebalances | 8 |

## Actions taken

1. Reverted the out-of-scope OpenMP tracker atomic source candidate:
   - `src/lagrangian/dsmc/faceTracker/dsmcFaceTracker.C`
   - `src/lagrangian/dsmc/parcels/dsmcParcel.C`
2. Restored the case control to pure MPI forced8:

```text
useOpenMP false;
openmpThreads 1;
endTime 5.e-05;
profileDetail false;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Restored controlDict hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

3. Rebuilt after the revert:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_scope_correction_rebuild_20260607.log
```

Build result: passed.

## Current valid comparison state

After reverting the out-of-scope OpenMP2 candidate and rebuilding, the pure
MPI8 forced8 case was rerun.

Current same-configuration log:

```text
log.codex_mpi8_replicatedmesh_purempi_scope_corrected_rebuild_500step_20260607
```

Result:

| metric | corrected current pure MPI forced8 | pure MPI8 reference | gap |
| --- | ---: | ---: | ---: |
| real | 131.16 s | 117.31 s | +16.85 s |
| move only | 70.58583509 s | 47.1652636 s | +24.42057149 s |
| buildCellOccupancy | 7.037232476 s | 10.7084157 s | -3.671183224 s |
| collision phase | 14.06214838 s | 10.89824478 s | +3.16490360 s |
| post fields/output | 41.44612826 s | 30.35563996 s | +11.09048830 s |
| full evolve wall | 129.1197774 s | 114.9207628 s | +14.1990146 s |
| DLB rebalances | 8 | 8 | 0 |

Correctness:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- final particles `2463800`
- stuck particles `0`
- final collisions `34542`
- final candidates `61800`
- final total energy `1.238986224`
- no `Fatal`, `Segmentation`, `Floating`, `NaN`, `nan`, or
  `BAD TERMINATION` match was found.

Previous same-round valid baseline, kept as historical context:

```text
log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607
```

| metric | current pure MPI forced8 | pure MPI8 reference | gap |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 117.31 s | +16.52 s |
| move only | 72.9567669 s | 47.1652636 s | +25.7915033 s |
| full evolve wall | 132.7056254 s | 114.9207628 s | +17.7848626 s |

The next candidate must be accepted only if it improves the current pure MPI
500-step baseline and is measured with `useOpenMP false`, `openmpThreads 1`,
and `mpirun -np 8`.

## Reference-control auto-DLB check

The reference `ourmeshbkp` controlDict does not set
`replicatedMeshDLBForceSteps`.  A same-source pure MPI8 run was therefore made
with the force-step list cleared, while keeping OpenMP disabled:

```text
useOpenMP false;
openmpThreads 1;
replicatedMeshDLBForceSteps ();
```

Backup before the temporary control change:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_purempi_reference_auto_dlb_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Log:

```text
log.codex_mpi8_replicatedmesh_purempi_reference_auto_dlb_500step_20260607
```

Result:

| metric | forced8 current | auto-DLB reference-control check |
| --- | ---: | ---: |
| real | 131.16 s | 132.49 s |
| move only | 70.58583509 s | 73.89721341 s |
| buildCellOccupancy | 7.037232476 s | 7.716358346 s |
| collision phase | 14.06214838 s | 14.55972334 s |
| post fields/output | 41.44612826 s | 41.32726588 s |
| full evolve wall | 129.1197774 s | 132.4919131 s |
| DLB rebalances | 8 | 5 |
| particles max/min | 4.082128483 | 2.956969229 |

Correctness:

- `Total Iterations = 500`
- `OpenMP enabled = 0`
- `OpenMP max threads = 1`
- final particles `2463914`
- stuck particles `0`
- final total energy `1.240376762`
- no `Fatal`, `Segmentation`, `Floating`, `NaN`, `nan`, or
  `BAD TERMINATION` match was found.

Decision: do not retain auto-DLB for the current source state.  It improves the
final particle max/min ratio but worsens `real`, `move only`, and
`full evolve wall`.  The case was restored to forced8 after the check:

```text
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```
