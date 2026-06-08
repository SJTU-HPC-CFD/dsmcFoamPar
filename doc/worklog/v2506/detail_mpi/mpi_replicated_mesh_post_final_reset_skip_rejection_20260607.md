# MPI replicated mesh final-output reset-skip attempt - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI only, `useOpenMP false`, `openmpThreads 1`,
  `mpirun -np 8`
- Formal retained control after this attempt:
  `replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);`
- Accepted forced8 control hash:
  `8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2`

## Candidate

File:

```text
src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C
```

Temporary source change:

- Skip `resetFieldsAtOutput()` work when the output time is the final
  `endTime` output.
- Add wall-time accumulation for the existing `output reset` profile field.

Rationale:

- `fieldPropertiesDict` has `resetAtOutput on`;
- the final output reset has no later sampling window to serve;
- this does not change tracking, collisions, DLB, or sampled field math before
  the final write.

Important scope correction:

- The accepted formal control has `writeInterval 1.e-3` and `endTime 5.e-05`,
  so this case normally does not enter the `outputTime()` field-write/reset path.
- The 10-step post-detail smoke also reported `output compute = 0`,
  `field writes = 0`, and `output reset = 0`.
- Therefore the candidate did not target the active 500-step formal hotspot; the
  active cost is per-step field sampling/post work, not final field output.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_final_reset_skip_build_20260607.log
```

Build result: passed.

## 10-step smoke

Temporary controls:

```text
endTime 1.e-06;
profileDetail true;
profilePostDetail true;
useOpenMP false;
openmpThreads 1;
```

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_post_final_reset_skip_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_post_final_reset_skip_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 11.62 s |
| move only | 1.053147492 s |
| post fields/output | 0.834544903 s |
| full evolve wall | 2.492067129 s |
| OpenMP enabled | 0 |
| DLB rebalances | 0 |
| final particles | 2065872 |
| stuck particles | 0 |
| final collisions | 2826 |
| final candidates | 2926 |
| final total energy | 1.150639951 |

Smoke passed; no Fatal/NaN/BAD TERMINATION string was found.

The post timer was slightly lower than the recent same-control smoke, but this
is only a 10-step signal and the candidate can only remove final-output reset
work.  It is not enough to decide retention.

## 200-step signal

Temporary controls:

```text
endTime 2.e-05;
profileSummary true;
profileDetail false;
useOpenMP false;
openmpThreads 1;
```

Signal log:

```text
log.codex_mpi8_replicatedmesh_post_final_reset_skip_signal_200step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 200 |
| real | 54.36 s |
| move only | 26.41077164 s |
| move+collide wall | 34.8785451 s |
| buildCellOccupancy | 2.392489168 s |
| collision phase | 5.06123331 s |
| post fields/output | 16.33986782 s |
| full evolve wall | 48.47341764 s |
| migration wall time | 0.926998719 s |
| DLB checks | 200 |
| DLB rebalances | 2 |
| Phase C wall max | 0.798943981 s |
| Phase C migration max | 0.649998042 s |
| particles max/min | 2.609631215 |
| rank wall max/min | 1.092145529 |
| OpenMP enabled | 0 |
| final particles | 2220226 |
| stuck particles | 0 |
| final collisions | 15506 |
| final candidates | 25463 |
| final total energy | 1.169789004 |

Correctness:

- `Total Iterations = 200`;
- `OpenMP enabled = 0`, `OpenMP max threads = 1`;
- final stuck particles `0`;
- no Fatal/NaN/BAD TERMINATION string was found.

## Decision

Rejected and reverted.  Do not promote this candidate to 500-step formal.

Comparison notes:

- Current strict screening baseline remains
  `log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607`:
  `real 133.83 s`, `move only 72.9567669 s`,
  `post fields/output 40.37082178 s`.
- A simple 200/500 scale of that strict baseline gives about `53.53 s` real and
  `16.15 s` post/output.  The candidate signal was `54.36 s` real and
  `16.33986782 s` post/output, so the target timer did not show a clear win.
- Existing 200-step config-screening logs are not like-for-like source
  baselines, but they put this result in the same narrow range:
  `vsizeExp=0.25` had `real 55.00 s`, `post 17.03748328 s`;
  `DLBItr=100` had `real 54.93 s`, `post 17.1258495 s`.
- Since the formal run does not normally execute output reset at all, the
  plausible 500-step upside is effectively outside the current formal path and
  is not supported by the 200-step signal.

## Revert state

The source change was reverted, the case control was restored, and the binary
was rebuilt.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_final_reset_skip_revert_build_20260607.log
```

Restored control hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

There is no remaining unstaged diff in `dsmcVolFields.C` or the case
`system/controlDict` from this candidate.
