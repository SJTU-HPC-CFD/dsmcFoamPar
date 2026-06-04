# ourmesh three-method 500-step performance, 2026-06-05

Follow-up note for the later `omp8` active-mode source investigation:

- `omp8_activefix_followup_20260605.md`

Current effective-work and source-audit summary:

- `../effective_optimizations_current_code_audit_20260605.md`

## Scope

Case root:

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/pal-phd-ch3.3.1-cylinder-react/timecompare/ourmesh
```

The three prepared cases were run directly without cleaning or decomposing. `mpi8origin` used the prepared `processor0` through `processor7` directories and was run with `-parallel`. `mpi8replicatedmesh` was run in replicated mesh mode without `-parallel`.

Environment:

```bash
source /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/docs/env.sh
```

Commands:

```bash
cd case/pal-phd-ch3.3.1-cylinder-react/timecompare/ourmesh/mpi8origin
mpirun -np 8 dsmcFoam+ -parallel

cd ../mpi8replicatedmesh
mpirun -np 8 dsmcFoam+

cd ../omp8
dsmcFoam+
```

## Config notes

Run 1 used the case state as found:

- `mpi8origin`: `profileDetail true`
- `mpi8replicatedmesh`: `profileDetail true`
- `omp8`: `profileDetail false`

Run 2 was a repeat after setting all three methods to `profileDetail false`. `omp8` already had this setting. The two changed production-case files were backed up before editing:

- `mpi8origin/system/controlDict.profileDetailTrue.20260605.bak`
- `mpi8replicatedmesh/system/controlDict.profileDetailTrue.20260605.bak`

Current `system/controlDict` state after Run 2 is `profileSummary true` and `profileDetail false` in all three cases.

## Results, run 1

This run was not fully uniform because the two MPI methods had `profileDetail true`, while `omp8` already had `profileDetail false`.

| method | mode | profileDetail | external wall [s] | main loop [s] | full evolve [s] | move only [s] | move kernel [s] | collision [s] | auto DLB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `mpi8origin` | MPI 8, non-replicated, `-parallel` | true | 176.45 | 182.8824358 | 181.1585693 | 134.8334005 | 63.49502908 | 18.35668146 | off |
| `mpi8replicatedmesh` | MPI 8, replicated mesh | true | 141.97 | 140.2923869 | 139.2353594 | 61.64699802 | 60.18959779 | 11.25532484 | 8 rebalances, 22.04153597 s |
| `omp8` | OpenMP 8 | false | 124.32 | 125.8246969 | 122.3566544 | 71.4635492 | 60.47155209 | 11.82494929 | 0.000408767 s |

Final step summaries:

| method | particles | collisions | candidates |
|---|---:|---:|---:|
| `mpi8origin` | 2463429 | 35691 | 75917 |
| `mpi8replicatedmesh` | 2466014 | 35356 | 71537 |
| `omp8` | 2463573 | 35949 | 76170 |

`mpi8replicatedmesh` DLB trigger steps: `110, 160, 210, 260, 310, 360, 410, 460`.

## Results, run 2

This run used `profileDetail false` for all three methods.

| method | mode | profileDetail | external wall [s] | main loop [s] | full evolve [s] | move only [s] | move kernel [s] | collision [s] | auto DLB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `mpi8origin` | MPI 8, non-replicated, `-parallel` | false | 159.69 | 165.4740552 | 161.8667951 | 119.4785574 | 56.37815238 | 18.44549356 | off |
| `mpi8replicatedmesh` | MPI 8, replicated mesh | false | 127.37 | 126.2693867 | 125.199297 | 49.51627701 | 48.08979666 | 11.56662326 | 8 rebalances, 19.75266342 s |
| `omp8` | OpenMP 8 | false | 125.46 | 127.2241534 | 123.6637694 | 72.52985026 | 61.01169953 | 11.96272507 | 0.000353748 s |

Final step summaries:

| method | particles | collisions | candidates |
|---|---:|---:|---:|
| `mpi8origin` | 2463429 | 35691 | 75917 |
| `mpi8replicatedmesh` | 2466014 | 35356 | 71537 |
| `omp8` | 2463615 | 35625 | 75586 |

`mpi8replicatedmesh` final rank sums matched the replicated final summary exactly:

```text
particles:  2466014
collisions: 35356
candidates: 71537
```

## Results, run 3 confirmation

Run 3 was requested as a repeat confirmation after the OMP move profiling fix.
All three production case dictionaries were checked before running and were
already `profileSummary true` / `profileDetail false`.

Run id:

```text
confirm_pdFalse_20260605_024610
```

Logs:

```text
case/.../ourmesh/mpi8origin/log.confirm_pdFalse_20260605_024610
case/.../ourmesh/mpi8replicatedmesh/log.confirm_pdFalse_20260605_024610
case/.../ourmesh/omp8/log.confirm_pdFalse_20260605_024610
```

Commands:

```bash
cd case/pal-phd-ch3.3.1-cylinder-react/timecompare/ourmesh/mpi8origin
mpirun -np 8 dsmcFoam+ -parallel

cd ../mpi8replicatedmesh
mpirun -np 8 dsmcFoam+

cd ../omp8
OMP_NUM_THREADS=8 dsmcFoam+
```

Timing:

| method | mode | profileDetail | external wall [s] | main loop [s] | full evolve [s] | move only [s] | move kernel [s] | collision [s] | auto DLB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| `mpi8origin` | MPI 8, non-replicated, `-parallel` | false | 160.95 | 166.6245026 | 163.0536101 | 120.3658396 | 56.6983455 | 18.73261305 | off |
| `mpi8replicatedmesh` | MPI 8, replicated mesh | false | 117.31 | 118.3501321 | 114.9207628 | 47.1652636 | 45.76155029 | 10.89824478 | 8 rebalances, 25.67342357 s max wall |
| `omp8` | OpenMP 8 | false | 103.66 | 104.8437967 | 102.1225033 | 48.82365912 | 40.14260107 | 9.876225069 | 0 |

Final step summaries:

| method | particles | collisions | candidates | average total energy |
|---|---:|---:|---:|---:|
| `mpi8origin` | 2463429 | 35691 | 75917 | 1.031952906e-18 |
| `mpi8replicatedmesh` | 2466342 | 35158 | 70116 | 1.023733261e-18 |
| `omp8` | 2463391 | 35847 | 75838 | 1.032208342e-18 |

Run 3 external wall time ranks:

```text
omp8:               103.66 s
mpi8replicatedmesh: 117.31 s
mpi8origin:         160.95 s
```

Relative to Run 2, the current source keeps `mpi8origin` essentially unchanged,
improves replicated mesh, and restores `omp8` to the fastest method:

| method | run 2 wall [s] | run 3 wall [s] | delta [s] |
|---|---:|---:|---:|
| `mpi8origin` | 159.69 | 160.95 | +1.26 |
| `mpi8replicatedmesh` | 127.37 | 117.31 | -10.06 |
| `omp8` | 125.46 | 103.66 | -21.80 |

## Checks

All six Run 1/Run 2 logs reached `End`. For the Run 2 logs, no `FOAM FATAL`,
`MPI_ABORT`, `nan`, or segmentation failure was found.  The three Run 3
confirmation logs also reached `End` and had no `FOAM FATAL`, `MPI_ABORT`,
`Abort(`, `nan`, or segmentation failure markers.

The logs contain the existing non-fatal `IOstreamOption::compressionEnum`
warnings; `omp8` and replicated-mesh MPI runs also ended with the known
non-fatal `Finalizing MPI, but was initialized elsewhere` warning after `End`.

## Interpretation

Run 2 external wall time ranks:

```text
omp8:               125.46 s
mpi8replicatedmesh: 127.37 s
mpi8origin:         159.69 s
```

Changing `profileDetail` to `false` reduced the measured wall time for the two MPI methods in this repeat:

| method | run 1 wall [s] | run 2 wall [s] | delta [s] | delta |
|---|---:|---:|---:|---:|
| `mpi8origin` | 176.45 | 159.69 | -16.76 | -9.50% |
| `mpi8replicatedmesh` | 141.97 | 127.37 | -14.60 | -10.28% |
| `omp8` | 124.32 | 125.46 | +1.14 | +0.92% |

`omp8` was already `profileDetail false` in Run 1, so its Run 1 to Run 2 difference should be read as repeat-run variation, not a configuration effect. In the uniform `profileDetail false` run, `mpi8replicatedmesh` is close to `omp8` on external wall time and remains much faster than `mpi8origin`. The remaining replicated-mesh cost is still dominated by DLB/rebalance and post-field/output timing rather than collision work alone.

## Logs

Full logs copied into this directory:

```text
log.mpi8origin.perf_20260605
log.mpi8replicatedmesh.perf_20260605
log.omp8.perf_20260605
log.mpi8origin.perf_20260605.profileDetailFalse
log.mpi8replicatedmesh.perf_20260605.profileDetailFalse
log.omp8.perf_20260605.profileDetailFalse
log.mpi8origin.confirm_pdFalse_20260605_024610
log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610
log.omp8.confirm_pdFalse_20260605_024610
```

Run 3 confirmation bundle with logs and matching `controlDict` files:

```text
confirm_pdFalse_20260605_024610_bundle/mpi8origin.log
confirm_pdFalse_20260605_024610_bundle/mpi8origin.controlDict
confirm_pdFalse_20260605_024610_bundle/mpi8replicatedmesh.log
confirm_pdFalse_20260605_024610_bundle/mpi8replicatedmesh.controlDict
confirm_pdFalse_20260605_024610_bundle/omp8.log
confirm_pdFalse_20260605_024610_bundle/omp8.controlDict
```
