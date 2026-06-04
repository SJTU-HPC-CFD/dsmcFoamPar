# OMP8 active-mode follow-up, 2026-06-05

## Scope

This follow-up investigated why the `omp8` run from Run 2 was slower than the older
`omp8/log.dsmcFoam+.omp8` baseline.  No production case directory was cleaned or
rerun for this follow-up.  Runtime validation used:

```text
/tmp/dsmc_omp_activefix_20260605/omp8
```

The validation case is a copy of the `omp8` case with:

- `startTime 0`
- `endTime 5.e-05`
- `deltaT 1.e-07`
- `profileDetail false`
- `replicatedMesh true`
- OpenMP threads: 8

A recursive comparison between the validation copy and the current production
`omp8` case showed only extra log files in `/tmp`; the case dictionaries, mesh,
and initial fields matched.

## Source changes kept

The final kept source changes are:

- `src/lagrangian/basic/Cloud/Cloud.C`
  - After comparison with `src/gitbkp/Cloud.C`, removed the extra
    `moveTrackCallCount` profiling plumbing from the OpenMP move path.
  - The current OpenMP move branch now matches the `gitbkp` move algorithm for
    this case; remaining differences are outside the active OMP move kernel or
    are the older wall-time fields that also exist in `gitbkp`.
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
  - Replaced several single-rank-sensitive `replicatedMesh_.valid()` checks with
    `replicatedMeshActive()` in `evolve()` and `reportProfiling()`.
  - This prevents inactive single-rank replicated mesh objects from entering
    migration, rebalance, rank-time, and replicated report paths.
- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/parcels/dsmcParcel.H`
  - Removed the `moveTrackCallCount` storage, accumulation, and summary output
    chain that was absent from `src/gitbkp`.
- `src/lagrangian/dsmc/parcels/dsmcParcel.C`
  - Split `dsmcParcel::move()` so `profileDetail false` bypasses detailed move
    timing/counting work.
  - `moveItersPerCell_` and `moveItersPerCellCumulative_` updates now occur only
    in the detailed profiling path.

One additional experiment tightened the `hitPatch()` and `hitWallPatch()` false
path.  It did not improve the run and was reverted before the final rebuild.

## Build

Command:

```bash
cd /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
source docs/env.sh && wmake libso src/lagrangian/dsmc
```

Result: build passed.  The only warning observed was the pre-existing
`Cloud.C:847:16: warning: unused variable 'useParticlePartition'`.

## Validation logs

Sandboxed single-rank runs still fail during Intel MPI initialization with:

```text
OFI endpoint open failed
```

That failed sandbox log is:

```text
/tmp/dsmc_omp_activefix_20260605/omp8/log.activefix_fastpath_20260605
```

Valid timing runs were performed outside the sandbox:

```text
/tmp/dsmc_omp_activefix_20260605/omp8/log.activefix_profileguard_20260605
/tmp/dsmc_omp_activefix_20260605/omp8/log.activefix_dsmcparcel_split_20260605
/tmp/dsmc_omp_activefix_20260605/omp8/log.activefix_fastpath_unsandbox_20260605
/tmp/dsmc_omp_activefix_20260605/omp8/log.activefix_cloud_trackcount_unsandbox_20260605
/tmp/dsmc_omp_activefix_20260605/omp8/log.activefix_cloud_trackcount_nobind_unsandbox_20260605
```

All listed valid timing logs reached `End` and had no `FOAM FATAL`,
`MPI_ABORT`, `nan`, or segmentation fault markers.  They ended with the known
non-fatal warning:

```text
Finalizing MPI, but was initialized elsewhere
```

## Timing comparison

| log | external wall [s] | main loop [s] | full evolve [s] | move only [s] | move kernel [s] | move total [s] | build total [s] | collision [s] | post fields/output [s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| old baseline `case/.../omp8/log.dsmcFoam+.omp8` | n/a | 107.5181266 | n/a | 52.29264194 | 44.02423736 | 50.78470898 | 6.064044591 | 9.534489021 | 37.21961824 |
| Run 2 `log.omp8.perf_20260605.profileDetailFalse` | 125.46 | 127.2241534 | 123.6637694 | 72.52985026 | 61.01169953 | 70.42114138 | 7.190830959 | 11.96272507 | 31.97915458 |
| kept fix `log.activefix_profileguard_20260605` | 118.35 | 120.6044536 | 117.680763 | 69.50859184 | 58.98261391 | 67.55697343 | 7.045966163 | 11.64806653 | 29.48552793 |
| parcel false-path split `log.activefix_dsmcparcel_split_20260605` | 119.71 | 120.6232305 | 117.9381641 | 68.68627768 | 58.39925487 | 66.68204823 | 7.042852671 | 11.45707902 | 30.75795815 |
| Cloud track-count removed, bound `log.activefix_cloud_trackcount_unsandbox_20260605` | 104.88 | 104.559808 | 101.955195 | 48.53928657 | 39.9727554 | 46.90682001 | 6.48891632 | 9.834628195 | 37.10005607 |
| Cloud track-count removed, nobind `log.activefix_cloud_trackcount_nobind_unsandbox_20260605` | 103.77 | 103.4617452 | 100.914042 | 48.20069174 | 39.593751 | 46.69207251 | 6.434977963 | 9.831305983 | 36.45382145 |
| reverted experiment `log.activefix_fastpath_unsandbox_20260605` | 120.87 | 122.9996174 | 120.248936 | 70.34996158 | 59.72563745 | 68.37130077 | 7.083011571 | 11.48867526 | 31.33392443 |

Final-step consistency:

| log | particles | last-step collisions | average total energy |
|---|---:|---:|---:|
| old baseline | 2463509 | 35860 | 1.032168431e-18 |
| Run 2 slow log | 2463615 | 35625 | 1.03211336e-18 |
| kept fix | 2463398 | 35777 | 1.032158287e-18 |
| parcel false-path split | 2463358 | 35815 | 1.032281158e-18 |
| Cloud track-count removed, bound | 2463544 | 35876 | 1.032087876e-18 |
| Cloud track-count removed, nobind | 2463408 | 35652 | 1.032166391e-18 |
| reverted experiment | 2463316 | 35715 | 1.032218785e-18 |

## Conclusions

The `replicatedMeshActive()` and first `profileDetail false` move guard changes
recover part of the slowdown:

```text
external wall: 125.46 -> 118.35 s
main loop:     127.2241534 -> 120.6044536 s
move kernel:   61.01169953 -> 58.98261391 s
move only:     72.52985026 -> 69.50859184 s
```

The follow-up `dsmcParcel` false-path split only made a small move-kernel
difference and did not explain the remaining gap:

```text
move kernel: 58.98261391 -> 58.39925487 s
move only:   69.50859184 -> 68.68627768 s
```

The high-impact regression was the extra `moveTrackCallCount` profiling chain
found by comparing current `Cloud.C`/headers against `src/gitbkp`.  Removing it
restored the OMP8 move path:

```text
external wall: 118.35 -> 103.77 s
main loop:     120.6044536 -> 103.4617452 s
full evolve:   117.680763 -> 100.914042 s
move only:      69.50859184 -> 48.20069174 s
move kernel:    58.98261391 -> 39.593751 s
```

The nobind retest is important because an earlier bound run used
`OMP_PROC_BIND=close` and `OMP_PLACES=cores`.  The nobind run only set
`OMP_NUM_THREADS=8` and still reproduced the improvement, so the result is not
explained by OpenMP binding variables.

Compared with the older normal baseline, the current nobind result is now faster
on the main move metrics:

```text
main loop:   103.4617452 vs 107.5181266 s
move only:    48.20069174 vs 52.29264194 s
move kernel:  39.593751 vs 44.02423736 s
```

The old baseline was produced by an earlier code path that still reported
single-rank replicated mesh output as `[rank 0]` and printed a replicated mesh
final summary.  The current active-mode semantics do not match that older output
path.  That output-path difference should not be reintroduced just to match the
old log text; the performance regression itself is explained by the
`moveTrackCallCount` chain.

Additional source inspection after this run:

- `src/lagrangian/basic/Cloud/Cloud.C` no longer has the extra
  `moveTrackCallCount` OMP move plumbing that differed from `src/gitbkp`.
- The kept `profileDetail false` branch in `dsmcParcel::move()` now mirrors the
  old direct track/control flow closely.
- The remaining `dsmcParcel` detailed-profiling branches are not entered in this
  validated `profileDetail false` OMP8 path.
- More `hitPatch()` or `hitWallPatch()` fast paths are not supported by the
  timing evidence.

## Source comparison against `src/gitbkp`

The normal-performance OMP source backup was found at:

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/gitbkp
```

It contains the key files needed for comparison:

```text
Cloud.C
dsmcCloud.C
dsmcCloud.H
dsmcCloudI.H
dsmcParcel.C
dsmcParcel.H
```

High-signal diffs:

- `gitbkp/dsmcParcel.C` vs current `src/lagrangian/dsmc/parcels/dsmcParcel.C`
  - Current code adds `<chrono>`, `MoveProfileClock`, and `elapsedSeconds()`.
  - Current code adds `profileMoveDetail` and a detailed move path with
    `moveTrackWallTime`, `moveItersPerCell_`, face/cyclic hit counters, and
    boundary timing.
  - The kept current false path is now close to `gitbkp`, but the file still has
    extra profiling branches in `hitPatch()`, `hitProcessorPatch()`, and
    `hitWallPatch()`.
- `gitbkp/dsmcParcel.H` vs current `dsmcParcel.H`
  - `moveTrackCallCount` has now been removed from current `trackingData`.
  - Other move timing/counter fields already existed in `gitbkp`.
- `gitbkp/Cloud.C` vs current `src/lagrangian/basic/Cloud/Cloud.C`
  - The move algorithm is not materially different.
  - The `moveThreadTrackCallCounts` plumbing has now been removed.
  - Current `Cloud.C` still has serial-path inner-profile bookkeeping, but that
    is not entered by the validated OMP8 move kernel.
- `gitbkp/dsmcCloud.C/H/I` vs current cloud files
  - Current code adds `moveItersPerCell_`, `moveItersPerCellCumulative_`,
    extra evolve-profile fields, and `reportMoveCellHotspots()`.
  - Most hotspot/report work is guarded by `profileDetail true`.
  - Current code also includes the retained `replicatedMeshActive()` fix that
    recovered part of the slowdown and should not be mixed into the first
    profiling-regression isolation unless testing a pure historical baseline.

Isolation result:

1. `dsmcParcel` false-path cleanup was reasonable but not the main cause.
2. `Cloud.C`/header `moveTrackCallCount` cleanup was the main fix.
3. Full `dsmcCloud.C/H/I` replacement is not recommended for this regression;
   it would mix performance cleanup with the retained `replicatedMeshActive()`
   semantic fix.
