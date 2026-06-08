# MPI replicated mesh move-detail profile - 2026-06-06

## Scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI, `useOpenMP false`, `openmpThreads 1`
- Retained formal configuration after the DLB sweep:
  `replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);`
- Purpose: locate the MPI-only move hot path before making more tracking changes.

This is a diagnostic run only.  The formal performance gate remains the 500-step
run against the retained forced8 configuration and the `ourmeshbkp`
replicated-mesh reference.

## Instrumentation

Default-off diagnostic counters/timers were added to the existing profiling
path:

- Control switch: `moveDetailProfile`
- Enabled only when both `profileDetail true` and `moveDetailProfile true`
- Counts are reduced by global sum:
  parcels, `trackToFace` calls, face hits, processor hits, cyclic hits, patch
  hits, and stuck hits
- Timing is reduced by rank max:
  total time around tracking calls, boundary handling time, and reserved tracker
  time

The production case was restored after the diagnostic run:

```text
profileDetail false;
useOpenMP false;
openmpThreads 1;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

## 10-step diagnostic run

Temporary controls:

```text
endTime 1.e-06;
profileDetail true;
moveDetailProfile true;
```

Run log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_mpionly_move_detail_profile_smoke_10step_20260606
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 10.27 s |
| move only | 1.17689206 s |
| move detail parcels | 20,630,641 |
| move detail track calls | 25,758,269 |
| move detail face hits | 5,135,781 |
| move detail processor hits | 0 |
| move detail cyclic hits | 0 |
| move detail patch hits | 16,452 |
| move detail stuck hits | 0 |
| move detail track max | 0.719672337 s |
| move detail boundary max | 0.001153038 s |
| OpenMP enabled | 0 |
| OpenMP max threads | 1 |

## Interpretation

The boundary path is not the move bottleneck in this case:

- Boundary handling is about `0.0012 s` over the 10-step diagnostic run.
- There are no processor/cyclic hits in this replicated-mesh owner-local path.
- Patch hits are present but small relative to total tracking volume.

The hot region is the old tracking core itself:

- The run executed about `25.76M` `trackToFace(..., true)` calls in 10 steps.
- The max-rank time around those calls is about `0.720 s`, while total move is
  about `1.177 s`.
- The remaining move time is therefore in the tracking call itself plus the
  per-segment geometry/constraint work around it, not in wall/cyclic/processor
  boundary model callbacks.

## Decision

Do not prioritize more DLB cadence tests, boundary-model fast paths, or wrapper
style tracking API migration for the next candidate.  The next source candidate
should be a local, easily reverted fast path inside the old
`particle::trackToFace(..., true)` core or the reduced-D constraint helpers it
calls.

Acceptance for any new source candidate remains strict:

- Build must pass.
- A short MPI-only smoke run may be used only for correctness and diagnostic
  direction.
- A retained candidate must beat the current forced8 500-step result on both:
  `real 132.29 s` and `move only 69.49945105 s`.
