# Stage 1 Minimal Profile Timer Port - 2026-06-05

## Scope

This is the first Stage 1 source change. It intentionally ports only the low-disturbance profiling controls needed to create fresh OFv1706 baselines.

Implemented:

- `profileSummary` controlDict switch in `dsmcCloud`.
- `profileDetail` controlDict switch, currently only reported as `basic-stage timers only`.
- Accumulated per-run timing for:
  - full `dsmcCloud::evolve()`;
  - `evolve_moveAndCollide()`;
  - `Cloud<dsmcParcel>::move()`;
  - `buildCellOccupancy()` calls that occur inside `evolve()`;
  - `collisions()`;
  - `evolve_fields()`.
- Solver-side `dsmc.printProfileSummary()` at stage end.
- Parallel summaries use max reduction over ranks for stage timings and print from rank 0.

Not implemented in this slice:

- OpenMP move/collision data flow.
- `moveOrderedParcels_`.
- flat occupancy storage.
- FastRng.
- replicated mesh.
- thread-level diagnostics.

## Files Changed

- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`
- `applications/solvers/discreteMethods/dsmc/dsmcFoam+/dsmcFoam+.C`

The timer uses `Time::elapsedCpuTime()` rather than `elapsedClockTime()` because the OFv1706 clock-time accessor is integer-second granularity on this build and produced zero-valued one-step smoke timings.

`buildCellOccupancy()` timing is gated by `profileTimingActive_`, so constructor-time initial occupancy creation is not counted in solver-step summaries.

## Build Verification

Successful builds:

- `doc/worklog/v2506/detail/stage1_profile_build_20260605.log`
- `doc/worklog/v2506/detail/stage1_profile_build_after_cpu_timer_20260605.log`
- `doc/worklog/v2506/detail/stage1_profile_build_after_timing_gate_20260605.log`

The build emits existing OpenFOAM template warnings but completes and refreshes `libdsmcFoam+.so`, `dsmcFoam+`, and `dsmcInitialise+`.

## Smoke Verification

Primary smoke:

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/omp8`
- Temporary controlDict changes:
  - `nTerminalOutputs 1`
  - `endTime 1.e-07`
  - `profileSummary true`
  - `profileDetail false`
- Command: `source doc/scripts/env.sh && OMP_NUM_THREADS=1 dsmcFoam+`
- Log: `doc/worklog/v2506/detail/stage1_profile_smoke_ourmeshbkp_omp8_timing_gate_20260605.log`

Result:

```text
DSMC solver profile summary
    solver profile steps          = 1
    move+collide wall [s]         = 0.48
    move only [s]                 = 0.39
    buildCellOccupancy [s]        = 0.03
    collision phase [s]           = 0.06
    post fields/output [s]        = 0.08
    total profiled [s]            = 0.56
    full evolve wall [s]          = 0.56
```

Correctness smoke values:

- Collisions: 2611
- DSMC particles: 2,058,569
- Total energy: 1.150482424
- No `FOAM FATAL`, `MPI_ABORT`, or `nan` markers in the smoke tail.

The temporary `controlDict` was restored after the run, and generated `boundaries` / `fieldMeasurements` directories were removed.

## Failed Smoke Note

`heatBath-5species` failed before entering `evolve()`:

- Log: `doc/worklog/v2506/detail/stage1_profile_smoke_heatBath_20260605.log`
- Failure: floating point exception in `dsmcVolFields::calculateWallUnitVectors()` during field creation.

This failure is outside the new profile code path and is not used as the Stage 1 verification case.

## Next Step

With fresh OFv1706 profiling now available, the next Stage 1 slice should port move-path structure in small compile-tested increments:

1. Extract/survivor rebuild accounting without OpenMP behavior change.
2. Add optional move profiling sub-stages.
3. Only then introduce OpenMP move scheduling and ordered parcel storage.
