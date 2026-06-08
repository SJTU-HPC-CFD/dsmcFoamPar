# Stage9 Move Inline Reset and DSMC Stationary-Tet Fast Path Results - 2026-06-05

## Scope

Implemented and validated two move-focused changes:

- `src/lagrangian/basic/Cloud/Cloud.C`
  - keeps the serial `stepFraction = 0` sweep only for the non-OpenMP move path;
  - moves `stepFraction = 0` into the OpenMP move kernel;
  - in the MPI+OpenMP transfer loop, resets only during the first pass so received transfer particles are not reset repeatedly.
- `src/lagrangian/basic/particle/particleTemplates.C`
  - adds the stationary-mesh DSMC fast path in the DSMC `trackToFace(endPosition, td, true)` overload;
  - if the destination remains inside the current tet and no wall-impact-distance handling is required, sets `position_ = endPosition` and returns immediately;
  - keeps the ordinary `trackToFace(endPosition, td)` template out of this DSMC-only fast path;
  - keeps per-call `DynamicList<label> tris(4)` and cached mesh references in the tracking path to avoid shared `cloud.labels()` state under OpenMP and reduce repeated mesh lookups.

## Build

Log:

- `doc/worklog/v2506/detail/stage9_move_inline_reset_fastpath_build_20260605.log`

Result:

- `dsmcFoam+` rebuilt at `platforms/linux64IccDPInt32Opt/bin/dsmcFoam+`;
- `dsmcInitialise+` rebuilt at `platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+`;
- only the known OpenFOAM v1706 template-instantiation warnings were present;
- no compile or link errors were reported.

## 10-Step Pure OMP Smoke

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage9_inline_reset_fastpath_smoke_10step_20260605`

Temporary control settings:

- started from `ourmesh/omp8/system/controlDict`;
- changed `endTime` to `1.e-06`;
- added `profileSummary true`, `profileDetail false`;
- enabled `useOpenMP true`, `openmpThreads 8`, `openmpMove true`,
  `openmpMoveSchedule static`, `openmpMoveChunk 64`,
  `openmpCollisionSchedule dynamic`, `openmpCollisionChunk 8`, and
  `openmpMoveOrderedReuse true`;
- restored the original controlDict after the run.

Result:

| metric | value |
| --- | ---: |
| real | 7.25 s |
| user | 18.71 s |
| sys | 1.35 s |
| full evolve wall | 1.841310943 s |
| move+collide wall | 1.128834449 s |
| move only | 0.987154392 s |
| buildCellOccupancy | 0.112296824 s |
| collision phase | 0.01947797 s |
| post fields/output | 0.712365127 s |
| collisions | 2847 |
| particles | 2065844 |
| stuck particles | 0 |
| total energy | 1.150621639 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

No `Fatal`, `NaN`, `Floating point`, or `Segmentation` entries were reported.

Smoke comparison against stage8 pure OMP smoke:

| metric | stage8 | stage9 | delta |
| --- | ---: | ---: | ---: |
| move only | 1.040756554 s | 0.987154392 s | 0.053602162 s faster, 5.15% lower |
| particles | 2065850 | 2065844 | -6 (-0.000290%) |
| collisions | 2867 | 2847 | -20 (-0.697593%) |
| total energy | 1.150621191 | 1.150621639 | +4.48e-07 (+0.0000389%) |

## 10-Step MPI+OMP Smoke

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8origin/log.codex_mpi8_stage9_inline_reset_fastpath_smoke_10step_20260605`

Temporary control settings:

- started from `ourmesh/mpi8origin/system/controlDict`;
- changed `endTime` to `1.e-06`;
- added `profileSummary true`, `profileDetail false`;
- enabled `useOpenMP true`, `openmpThreads 2`, `openmpMove true`,
  `openmpMoveSchedule static`, `openmpMoveChunk 1`,
  `openmpCollisionSchedule dynamic`, `openmpCollisionChunk 1`, and
  `openmpMoveOrderedReuse true`;
- restored the original controlDict after the run.

The first sandboxed `mpirun` attempt failed before solver startup because Intel
MPI Hydra could not open a local bootstrap socket.  The rerun outside the
sandbox completed normally.

Result:

| metric | value |
| --- | ---: |
| real | 5.07 s |
| user | 49.25 s |
| sys | 6.62 s |
| full evolve wall | 3.153057971 s |
| move+collide wall | 2.815344517 s |
| move only | 2.635528082 s |
| buildCellOccupancy | 0.117927417 s |
| collision phase | 0.136523033 s |
| post fields/output | 0.575056466 s |
| collisions | 2912 |
| collision candidates | 3028 |
| particles | 2065860 |
| stuck particles | 0 |
| total energy | 1.150623507 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 2`
- `OpenMP move = 1 (static, chunk 1)`
- `OpenMP collision = dynamic, chunk 1`

No `Fatal`, `NaN`, `Floating point`, or `Segmentation` entries were reported.

Smoke comparison against stage8 MPI+OMP smoke:

| metric | stage8 | stage9 | delta |
| --- | ---: | ---: | ---: |
| move only | 2.156631315 s | 2.635528082 s | 0.478896767 s slower, 22.21% higher |
| particles | 2065866 | 2065860 | -6 (-0.000290%) |
| collisions | 2819 | 2912 | +93 (+3.299042%) |
| total energy | 1.150630334 | 1.150623507 | -6.827e-06 (-0.000593%) |

This smoke validates that the MPI+OpenMP path still starts, transfers, and
finishes cleanly.  It is not used as formal performance evidence.

## 500-Step Formal OMP Run

Formal comparison used the same `ourmeshbkp/omp8/system/controlDict` settings
as stage8.  For this run, `ourmeshbkp/omp8/system/controlDict` was temporarily
copied to `ourmesh/omp8/system/controlDict`, then the original
`ourmesh/omp8/system/controlDict` was restored.

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage9_inline_reset_fastpath_bkpctrl_500step_20260605`

Result:

| metric | value |
| --- | ---: |
| real | 127.02 s |
| user | 967.60 s |
| sys | 13.09 s |
| ClockTime | 127 s |
| full evolve wall | 128.4992481 s |
| move+collide wall | 76.50441151 s |
| move only | 64.61371202 s |
| buildCellOccupancy | 7.589475879 s |
| collision phase | 3.711549217 s |
| post fields/output | 51.98935525 s |
| move+collide cpu | 569.67 s |
| move only cpu | 480.91 s |
| buildCellOccupancy cpu | 56.01 s |
| collision phase cpu | 27.64 s |
| post fields/output cpu | 387.07 s |
| full evolve cpu | 956.76 s |
| collisions | 35672 |
| collision candidates | 68591 |
| collision acceptance rate | 0.5200682305 |
| particles | 2463607 |
| stuck particles | 0 |
| total energy | 1.243088811 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

No `Fatal`, `NaN`, `Floating point`, or `Segmentation` entries were reported.
The run completed 500 iterations and reached `End main`.

## Stage8 Comparison

| metric | stage8 move OMP complete | stage9 inline reset + fast path | delta |
| --- | ---: | ---: | ---: |
| real | 137.46 s | 127.02 s | 10.44 s faster, 7.59% lower |
| full evolve wall | 141.1479884 s | 128.4992481 s | 12.6487403 s faster, 8.96% lower |
| move+collide wall | 90.27564774 s | 76.50441151 s | 13.77123623 s faster, 15.25% lower |
| move only | 78.62462603 s | 64.61371202 s | 14.01091401 s faster, 17.82% lower |
| buildCellOccupancy | 7.449130794 s | 7.589475879 s | 0.140345085 s slower, 1.88% higher |
| collision phase | 3.631726566 s | 3.711549217 s | 0.079822651 s slower, 2.20% higher |
| post fields/output | 50.86713931 s | 51.98935525 s | 1.12221594 s slower, 2.21% higher |

The formal improvement is move-dominated.  The extra cost in
buildCellOccupancy, collision, and post fields/output is much smaller than the
move reduction.

## Correctness Check

Against stage8 move OMP complete:

| metric | stage8 | stage9 | difference |
| --- | ---: | ---: | ---: |
| particles | 2463681 | 2463607 | -74 (-0.003004%) |
| collisions | 35576 | 35672 | +96 (+0.269845%) |
| collision candidates | 68618 | 68591 | -27 (-0.039348%) |
| total energy | 1.243106701 | 1.243088811 | -1.789e-05 (-0.001439%) |
| stuck particles | 0 | 0 | unchanged |

The correctness drift relative to stage8 is small for this OMP regression
level.

## Historical Gap

Against historical `ourmeshbkp/omp8`:

| metric | historical bkp | stage9 inline reset + fast path | gap |
| --- | ---: | ---: | ---: |
| external/real wall | 103.66 s | 127.02 s | 23.36 s slower, 22.54% higher |
| full evolve wall | 102.1225033 s | 128.4992481 s | 26.3767448 s slower, 25.83% higher |
| move only | 48.82365912 s | 64.61371202 s | 15.7900529 s slower, 32.34% higher |
| buildCellOccupancy | 6.490896398 s | 7.589475879 s | 1.098579481 s slower, 16.92% higher |
| collision phase | 9.876225069 s | 3.711549217 s | 6.164675852 s faster, 62.42% lower |
| post fields/output | 36.92670834 s | 51.98935525 s | 15.06264691 s slower, 40.79% higher |
| particles | 2463391 | 2463607 | +216 (+0.008768%) |
| collisions | 35847 | 35672 | -175 (-0.488186%) |
| collision candidates | 75838 | 68591 | -7247 (-9.555895%) |

Stage9 closes a significant part of the stage8 move gap, but the historical
gap remains dominated by move and post fields/output.

## Cleanup

Restored control files:

- `ourmesh/omp8/system/controlDict` restored to
  `ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d`;
- `ourmesh/mpi8origin/system/controlDict` restored to
  `ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d`.

## Next Move Targets

The remaining OMP gap is not from the old serial `stepFraction` reset anymore.
The next move-only work should profile or port the retained-order move kernel
breakdown: parcel extraction, parallel tracking kernel, survivor rebuild, and
post-move ordered list maintenance.  Historical bkp has much lower move
kernel and post/output totals, so the next useful stage should add comparable
move sub-timers in this OFv1706 path before making another larger code change.
