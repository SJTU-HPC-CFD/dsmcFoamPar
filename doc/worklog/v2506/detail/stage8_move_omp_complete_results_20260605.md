# Stage8 OMP Move Complete Results - 2026-06-05

## Scope

This stage completed the retained move-OMP port for the current OFv1706 tree.
It still did not copy the newer v2506 Pstream API directly.  Instead, the
stage keeps the OFv1706 `Cloud::move(td, trackTime)` interface and adapts the
move OMP path to the old `PstreamBuffers`, `IDLList`, and `finishedSends`
transfer mechanism.

Changes applied in this stage:

- pre-move append capture now starts before
  `controllers_.controlBeforeMove()` and `boundaries_.controlBeforeMove()`, so
  parcels injected by pre-move controls are visible to ordered move reuse;
- `Cloud::move()` keeps its idempotent append-capture begin, so the early begin
  in `dsmcCloud::evolve_moveAndCollide()` remains balanced with the existing
  post-move end;
- the OpenMP move path is no longer restricted to `!Pstream::parRun()`;
- under `Pstream::parRun()`, the OMP kernel classifies keep, switch-processor,
  and delete parcels, then the old OFv1706 parallel transfer path sends and
  receives the switched parcels.

## Build

Command:

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
```

Log:

- `doc/worklog/v2506/detail/stage8_move_omp_complete_build_20260605.log`

Result:

- build completed;
- `dsmcFoam+` rebuilt at `platforms/linux64IccDPInt32Opt/bin/dsmcFoam+`;
- `dsmcInitialise+` rebuilt at `platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+`;
- only the known OpenFOAM v1706 template-instantiation warnings were present.

## 10-Step Pure OMP Smoke

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage8_move_omp_complete_smoke_10step_20260605`

Result:

| metric | value |
| --- | ---: |
| real | 4.60 s |
| user | 18.57 s |
| sys | 1.28 s |
| full evolve wall | 1.854459089 s |
| move+collide wall | 1.171918681 s |
| move only | 1.040756554 s |
| buildCellOccupancy | 0.103318462 s |
| collision phase | 0.017478417 s |
| post fields/output | 0.68242016 s |
| collisions | 2867 |
| particles | 2065850 |
| stuck particles | 0 |
| total energy | 1.150621191 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

No `Fatal`, `NaN`, `Floating point`, or `Segmentation` entries were reported.

## 10-Step MPI+OMP Smoke

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8origin/log.codex_mpi8_stage8_omp_move_smoke_10step_20260605`

Temporary control settings:

- started from `ourmesh/mpi8origin/system/controlDict`;
- changed `endTime` to `1.e-06`;
- added `profileSummary true`, `profileDetail false`;
- enabled `useOpenMP true`, `openmpThreads 2`, `openmpMove true`,
  `openmpMoveSchedule static`, `openmpMoveChunk 1`, and
  `openmpMoveOrderedReuse true`;
- restored the original controlDict after the run.

The first sandboxed `mpirun` attempt failed before solver startup because Intel
MPI Hydra could not open a local bootstrap socket.  The approved rerun outside
the sandbox completed normally.

Result:

| metric | value |
| --- | ---: |
| real | 4.54 s |
| user | 42.00 s |
| sys | 5.92 s |
| full evolve wall | 2.600737969 s |
| move+collide wall | 2.302683754 s |
| move only | 2.156631315 s |
| buildCellOccupancy | 0.116979306 s |
| collision phase | 0.089318975 s |
| post fields/output | 0.519829267 s |
| collisions | 2819 |
| collision candidates | 2933 |
| particles | 2065866 |
| stuck particles | 0 |
| total energy | 1.150630334 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 2`
- `OpenMP move = 1 (static, chunk 1)`
- `OpenMP collision = dynamic, chunk 1`

No `Fatal`, `NaN`, `Floating point`, or `Segmentation` entries were reported.
No new time directory was created by this smoke run.

## 500-Step Formal Run

Formal comparison used the same control settings as stage7.  For this run,
`ourmeshbkp/omp8/system/controlDict` was temporarily copied to
`ourmesh/omp8/system/controlDict`, then the original
`ourmesh/omp8/system/controlDict` was restored.

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage8_move_omp_complete_bkpctrl_500step_20260605`

Result:

| metric | value |
| --- | ---: |
| real | 137.46 s |
| user | 1059.81 s |
| sys | 15.71 s |
| ClockTime | 137 s |
| full evolve wall | 141.1479884 s |
| move+collide wall | 90.27564774 s |
| move only | 78.62462603 s |
| buildCellOccupancy | 7.449130794 s |
| collision phase | 3.631726566 s |
| post fields/output | 50.86713931 s |
| move+collide cpu | 672.86 s |
| move only cpu | 586.2 s |
| buildCellOccupancy cpu | 54.67 s |
| collision phase cpu | 26.98 s |
| post fields/output cpu | 379.58 s |
| full evolve cpu | 1052.44 s |
| collisions | 35576 |
| collision candidates | 68618 |
| collision acceptance rate | 0.5184645428 |
| particles | 2463681 |
| stuck particles | 0 |
| total energy | 1.243106701 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

No `Fatal`, `NaN`, `Floating point`, or `Segmentation` entries were reported.
The run completed 500 iterations and reached `End main`.

## Stage7 Comparison

| metric | stage7 move tracking | stage8 move OMP complete | delta |
| --- | ---: | ---: | ---: |
| real | 151.06 s | 137.46 s | 13.60 s faster, 9.00% lower |
| full evolve wall | 156.8234055 s | 141.1479884 s | 15.6754171 s faster, 10.00% lower |
| move+collide wall | 106.8374108 s | 90.27564774 s | 16.56176306 s faster, 15.50% lower |
| move only | 95.5824105 s | 78.62462603 s | 16.95778447 s faster, 17.74% lower |
| buildCellOccupancy | 7.194255854 s | 7.449130794 s | 0.25487494 s slower, 3.54% higher |
| collision phase | 3.516419771 s | 3.631726566 s | 0.115306795 s slower, 3.28% higher |
| post fields/output | 49.98115219 s | 50.86713931 s | 0.88598712 s slower, 1.77% higher |

The gain is dominated by move.  The small regressions in buildCellOccupancy,
collision, and post fields/output are much smaller than the move reduction.

## Correctness Check

Against stage7 move tracking:

| metric | stage7 | stage8 | difference |
| --- | ---: | ---: | ---: |
| particles | 2463537 | 2463681 | +144 (+0.005845%) |
| collisions | 35601 | 35576 | -25 (-0.070223%) |
| collision candidates | 68574 | 68618 | +44 (+0.064164%) |
| total energy | 1.242846706 | 1.243106701 | +0.000259995 (+0.020919%) |
| stuck particles | 0 | 0 | unchanged |

The correctness drift relative to stage7 is small for this OMP regression
level.

## Historical Gap

Against historical `ourmeshbkp/omp8`:

| metric | historical bkp | stage8 move OMP complete | gap |
| --- | ---: | ---: | ---: |
| external/real wall | 103.66 s | 137.46 s | 33.80 s slower |
| full evolve wall | 102.1225033 s | 141.1479884 s | 39.0254851 s slower |
| move only | 48.82365912 s | 78.62462603 s | 29.80096691 s slower |
| buildCellOccupancy | 6.490896398 s | 7.449130794 s | 0.958234396 s slower |
| collision phase | 9.876225069 s | 3.631726566 s | 6.244498503 s faster |
| post fields/output | 36.92670834 s | 50.86713931 s | 13.94043097 s slower |
| particles | 2463391 | 2463681 | +290 (+0.011773%) |
| collisions | 35847 | 35576 | -271 (-0.756004%) |
| collision candidates | 75838 | 68618 | -7220 (-9.520291%) |

The remaining historical gap is still dominated by move and post
fields/output.  Collision remains faster than the historical bkp run.

## Cleanup

`ourmesh/mpi8origin/system/controlDict` and
`ourmesh/omp8/system/controlDict` were restored after the smoke and formal
runs.

Restored hashes:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d  ourmesh/mpi8origin/system/controlDict
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d  ourmesh/omp8/system/controlDict
```

No `dsmcFoam+` or `mpirun` process was left running by the stage8 commands.
