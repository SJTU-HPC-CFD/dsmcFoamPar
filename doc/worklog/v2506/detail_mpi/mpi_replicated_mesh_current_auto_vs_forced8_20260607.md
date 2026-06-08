# MPI replicated mesh current auto-vs-forced8 state - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI, `useOpenMP false`, `openmpThreads 1`
- Run length: 500 steps (`endTime 5.e-05`, `deltaT 1.e-07`)
- Purpose: record the current same-source comparison after the retained
  tracking/DLB candidate sequence and before starting another source candidate.

## Current source and case state

The current case control is auto DLB, not forced-step DLB:

```text
useOpenMP false;
openmpThreads 1;
profileSummary true;
profileDetail false;
replicatedMeshDLBDualConstraint true;
replicatedMeshAutoDLB true;
replicatedMeshDLBForceSteps ();
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBMinGapSteps 50;
replicatedMeshDLBInitialAlpha 0.8;
```

Current `system/controlDict` sha256:

```text
a536a0e7c00df8009e857acc4b442eee5306537e1c2d6449cc9f578a0946c5b9
```

Important working-tree boundary:

- staged MPI replicated-mesh port/output/profile changes remain present;
- unstaged source changes include the retained same-tet non-normalised
  `particleTemplates.C` fast path plus default-off move-detail diagnostics in
  `dsmcCloud`/`dsmcParcel`;
- there are many untracked `detail_mpi` controlDict backups and build logs from
  the candidate sequence.  They are historical evidence and should not be
  cleaned while continuing this round.

## Latest same-source 500-step comparison

### Forced8 retained-control rerun

Log:

```text
log.codex_mpi8_replicatedmesh_retained_forced8_500step_20260607
```

Configuration delta:

```text
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Result:

| metric | value |
| --- | ---: |
| real | 138.82 s |
| move+collide wall | 99.42313411 s |
| move only | 75.29161688 s |
| buildCellOccupancy | 7.004017032 s |
| collision phase | 15.01691516 s |
| post fields/output | 42.1227724 s |
| full evolve wall | 135.4172667 s |
| DLB checks | 500 |
| DLB rebalances | 8 |
| migration calls | 59 |
| migration wall | 2.189168896 s |
| particles max/min | 2.270407632 |
| rank wall max/min | 1.068532931 |

Correctness:

- `Total Iterations = 500`
- final particles `2463823`, stuck `0`
- final collisions `34360`, candidates `59932`, acceptance `0.5733164253`
- final total energy `1.238928632`
- no Fatal/NaN/BAD TERMINATION match was found in the log.

### Auto-DLB rerun

Log:

```text
log.codex_mpi8_replicatedmesh_retained_auto_dlb_500step_20260607
```

Configuration:

```text
replicatedMeshDLBForceSteps ();
```

Result:

| metric | value |
| --- | ---: |
| real | 135.91 s |
| move+collide wall | 99.40221419 s |
| move only | 76.10286716 s |
| buildCellOccupancy | 8.563354644 s |
| collision phase | 14.747785 s |
| post fields/output | 44.67486854 s |
| full evolve wall | 134.4745938 s |
| DLB checks | 500 |
| DLB rebalances | 5 |
| migration calls | 56 |
| migration wall | 1.823290831 s |
| particles max/min | 1.544313227 |
| rank wall max/min | 1.106024384 |

Correctness:

- `Total Iterations = 500`
- final particles `2463933`, stuck `0`
- final collisions `34577`, candidates `62292`, acceptance `0.5550793039`
- final total energy `1.240443059`
- no Fatal/NaN/BAD TERMINATION match was found in the log.

## Comparison

| metric | forced8 | auto DLB | auto - forced8 |
| --- | ---: | ---: | ---: |
| real | 138.82 s | 135.91 s | -2.91 s |
| move only | 75.29161688 s | 76.10286716 s | +0.81125028 s |
| buildCellOccupancy | 7.004017032 s | 8.563354644 s | +1.559337612 s |
| collision phase | 15.01691516 s | 14.747785 s | -0.26893016 s |
| post fields/output | 42.1227724 s | 44.67486854 s | +2.55209614 s |
| full evolve wall | 135.4172667 s | 134.4745938 s | -0.9426729 s |
| DLB rebalances | 8 | 5 | -3 |
| migration wall | 2.189168896 s | 1.823290831 s | -0.365878065 s |
| particles max/min | 2.270407632 | 1.544313227 | -0.726094405 |
| rank wall max/min | 1.068532931 | 1.106024384 | +0.037491453 |

Against the preserved `ourmeshbkp/mpi8replicatedmesh` external wall reference
(`117.31 s`), the latest auto-DLB run is still `+18.60 s` slower.  The current
run is also far slower than the earlier historical same-tet pair
(`127.91-129.93 s`), so that pair should be treated as historical evidence, not
as the current same-environment baseline.

## Decision

For the current source/object/runtime state, auto DLB is the better immediate
end-to-end baseline than forced8:

- it improves `real` by `2.91 s` and `full evolve wall` by `0.94 s`;
- it reduces migration/rebalance work (`5` rebalances instead of `8`);
- it improves final particle max/min balance;
- it does not improve the move timer, so the main optimization direction remains
  the lower-level old tracking core.

Because the recent retained same-tet reruns showed significant run-to-run or
rebuild-state variation, the next step should be a same-configuration auto-DLB
repeat before retaining auto DLB as the formal current baseline for new source
candidates.

## Next work

1. Repeat the current auto-DLB 500-step case with unchanged source and
   `replicatedMeshDLBForceSteps ()`.
2. If the repeat stays in the same performance class, use the auto-DLB result as
   the current pure-MPI baseline for subsequent source candidates.
3. Continue source candidates only in the old tracking core around
   `particle::trackToFace(..., true)`, same-tet checks, and internal-tet walk.
4. Keep formal acceptance strict: a retained source candidate must pass build,
   complete 500 steps, keep particle/collision/energy consistency, and improve
   the current same-environment baseline without relying on 10-step smoke
   timings alone.

## Auto-DLB repeat

The current auto-DLB configuration was repeated with unchanged source and case
controls.

Sandbox note:

- the first sandboxed launch failed before solver startup with the known Intel
  MPI/Hydra socket error:
  `HYD_sock_listen_on_port: cannot open socket (Operation not permitted)`;
- the same command was rerun outside the sandbox and completed normally.

Repeat log:

```text
log.codex_mpi8_replicatedmesh_retained_auto_dlb_repeat_500step_20260607
```

Result:

| metric | previous auto | auto repeat | repeat - previous |
| --- | ---: | ---: | ---: |
| real | 135.91 s | 139.68 s | +3.77 s |
| move+collide wall | 99.40221419 s | 105.819344 s | +6.41712981 s |
| move only | 76.10286716 s | 81.44518205 s | +5.34231489 s |
| buildCellOccupancy | 8.563354644 s | 8.868006041 s | +0.304651397 s |
| collision phase | 14.747785 s | 16.09319453 s | +1.34540953 s |
| post fields/output | 44.67486854 s | 43.73627522 s | -0.93859332 s |
| full evolve wall | 134.4745938 s | 136.0861873 s | +1.6115935 s |
| DLB rebalances | 5 | 4 | -1 |
| particles max/min | 1.544313227 | 2.010086841 | +0.465773614 |
| rank wall max/min | 1.106024384 | 1.1503207 | +0.044296316 |

Correctness:

- `Total Iterations = 500`
- final particles `2463871`, stuck `0`
- final collisions `34031`, candidates `58389`, acceptance `0.5828323828`
- final total energy `1.240740027`
- no Fatal/NaN/BAD TERMINATION match was found in the log.

Decision update:

- the auto-DLB repeat did not reproduce the earlier `135.91 s` result;
- it is also slower than the latest forced8 run in both `real` and `move only`;
- do not replace forced8 with auto DLB as the retained pure-MPI baseline from
  this evidence;
- restore the case to forced8 before continuing source candidates so the
  acceptance gate stays aligned with the retained DLB configuration.

## Case restore after auto repeat

The auto-DLB controlDict was backed up before restoring forced8:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_restore_forced8_after_auto_repeat_20260607
sha256 a536a0e7c00df8009e857acc4b442eee5306537e1c2d6449cc9f578a0946c5b9
```

Current formal case controls are restored to:

```text
useOpenMP false;
openmpThreads 1;
profileDetail false;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Current `system/controlDict` sha256 after restore:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Retained control for the next source-candidate run is therefore forced8, not
auto DLB.

## Forced8 repeat after restore

After restoring forced8, one more 500-step repeat was run to establish a
same-environment retained baseline.

Repeat log:

```text
log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607
```

Result:

| metric | previous forced8 | forced8 repeat | repeat - previous |
| --- | ---: | ---: | ---: |
| real | 138.82 s | 133.83 s | -4.99 s |
| move+collide wall | 99.42313411 s | 97.10060854 s | -2.32252557 s |
| move only | 75.29161688 s | 72.9567669 s | -2.33484998 s |
| buildCellOccupancy | 7.004017032 s | 6.729516136 s | -0.274500896 s |
| collision phase | 15.01691516 s | 14.8398869 s | -0.17602826 s |
| post fields/output | 42.1227724 s | 40.37082178 s | -1.75195062 s |
| full evolve wall | 135.4172667 s | 132.7056254 s | -2.7116413 s |
| DLB rebalances | 8 | 8 | 0 |
| particles max/min | 2.270407632 | 1.412911358 | -0.857496274 |
| rank wall max/min | 1.068532931 | 1.053059781 | -0.01547315 |

Correctness:

- `Total Iterations = 500`
- final particles `2463567`, stuck `0`
- final collisions `34522`, candidates `59845`, acceptance `0.5768568803`
- final total energy `1.239014073`
- no Fatal/NaN/BAD TERMINATION match was found in the log.

Decision:

- forced8 remains the retained pure-MPI DLB configuration;
- current same-environment baseline range is:

```text
real 133.83-138.82 s
move only 72.9567669-75.29161688 s
full evolve wall 132.7056254-135.4172667 s
```

- for a new source candidate, the strict first-pass acceptance target should be
  the better repeat result: improve both `real 133.83 s` and
  `move only 72.9567669 s` in a 500-step pure-MPI run;
- if a candidate only improves by a small margin, repeat it before retaining.

## Current forced8 move-detail diagnostic

The restored forced8 controlDict was backed up before the diagnostic:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_current_move_detail_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Temporary controls:

```text
endTime 1.e-06;
profileDetail true;
moveDetailProfile true;
```

Diagnostic log:

```text
log.codex_mpi8_replicatedmesh_retained_forced8_move_detail_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 11.69 s |
| move only | 1.322484047 s |
| buildCellOccupancy | 0.11536334 s |
| collision phase | 0.195757522 s |
| post fields/output | 0.945544449 s |
| full evolve wall | 2.763975522 s |
| move detail parcels | 20,630,514 |
| move detail track calls | 25,758,151 |
| same-tet no-face | 11,139,106 |
| internal tet only | 9,483,253 |
| face hits | 5,135,792 |
| processor hits | 0 |
| cyclic hits | 0 |
| patch hits | 16,453 |
| stuck hits | 0 |
| track max | 0.896837629 s |
| boundary max | 0.001501033 s |
| DLB checks | 10 |
| DLB rebalances | 0 |

Correctness:

- final particles `2065851`, stuck `0`
- final collisions `2824`, candidates `2930`, acceptance `0.9638225256`
- final total energy `1.150623554`
- no Fatal/NaN/BAD TERMINATION match was found in the log.

The production controlDict was restored after the diagnostic.  The restored
hash again matches the pre-diagnostic forced8 backup:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Interpretation:

- the current diagnostic reproduces the earlier shape of the move hot path:
  same-tet no-face, internal tet-only, and face-hit counts are still the
  dominant tracking buckets;
- boundary, cyclic, processor, and stuck paths remain negligible for this case;
- there is no evidence that another DLB cadence test or boundary callback
  change is the right next source candidate;
- the next useful step should be a more granular old-tracking-core diagnostic
  that splits the `track max` time into `findTris`, second `tetLambda` scan,
  `tetNeighbour`, and `crossEdgeConnectedFace` time before attempting another
  algorithmic fast path.

Implementation boundary for that diagnostic:

- `particle::trackToFace(..., bool)` is not DSMC-only; the same template overload
  is also called from PIC and molecular-dynamics parcels;
- do not directly reference DSMC-specific `moveDetail*` fields inside
  `particleTemplates.C`, because that would risk breaking other template
  instantiations;
- a finer timing pass should first add a generic/trait-safe diagnostic interface
  or an equivalent local mechanism that compiles for non-DSMC `TrackData` too.
