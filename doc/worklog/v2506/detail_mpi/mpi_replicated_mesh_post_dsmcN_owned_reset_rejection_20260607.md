# MPI replicated mesh instantaneous dsmcN owner-reset attempt - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI only, `useOpenMP false`, `openmpThreads 1`,
  `mpirun -np 8`
- Baseline control hash:
  `8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2`
- Source area:
  `src/lagrangian/dsmc/macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C`

This was a strict full-field sampling-body screen.  It did not touch OpenMP,
tracking, final output reset, no-fields mode, or shared-cache reuse across time
steps.

## Candidate

Temporary source change:

- in replicated mesh mode, replace the per-call full `dsmcN_ = 0.0` with an
  owner-cell-only clear for non-output sampling steps;
- keep a flag so that after an output-time MPI reduce of `dsmcN_`, the next
  sampling call performs a full clear before returning to owner-cell clears;
- add a temporary profile-detail line for `instantaneous reset`.

Rationale:

- `dsmcN_` is an instantaneous field and is reset every `calculateField()` call;
- the current replicated-mesh run samples only local owned cells before output
  reduction;
- if the full clear was material, owner-cell clearing could reduce per-step
  post work without changing full-field output semantics.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_dsmcN_owned_reset_build_20260607.log
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
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_post_dsmcN_owned_reset_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_post_dsmcN_owned_reset_smoke_10step_20260607
```

The first sandboxed launch failed before solver startup because Intel MPI Hydra
could not open its listening socket:

```text
HYD_sock_listen_on_port ... cannot open socket (Operation not permitted)
```

The same command was rerun outside the sandbox and completed.

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 9.58 s |
| move only | 1.192247914 s |
| buildCellOccupancy | 0.123001487 s |
| collision phase | 0.210958266 s |
| post fields/output | 0.865637994 s |
| full evolve wall | 2.492606307 s |
| dsmcVolFields instantaneous reset | 0.000559476 s |
| dsmcVolFields sample accumulation | 0.239334675 s |
| dsmcVolFields shared cache build | 0.167777043 s |
| dsmcVolFields cache allocate | 0.018111933 s |
| dsmcVolFields parcel accumulate | 0.148883679 s |
| dsmcVolFields field combine | 0.071534244 s |
| dsmcVolFields cell reduction | 0.005973127 s |
| dsmcVolFields boundary accumulation | 0.003791123 s |
| OpenMP enabled | 0 |
| OpenMP max threads | 1 |
| final particles | 2065869 |
| stuck particles | 0 |
| final collisions | 2948 |
| final candidates | 3060 |
| final total energy | 1.150624043 |

Smoke correctness passed; no Fatal/NaN/BAD TERMINATION string was found.

## Decision

Rejected and reverted.  Do not promote this candidate to 200-step signal.

The measured reset cost was negligible: `0.000559476 s` across 9 profiled
`dsmcVolFields` calls.  The overall `post fields/output` timer was
`0.865637994 s`, which is not better than the recent strict full-field detail
smoke (`0.834544903 s` in
`log.codex_mpi8_replicatedmesh_post_final_reset_skip_smoke_10step_20260607`).

This confirms that the active strict full-field cost is not the instantaneous
`dsmcN_` full clear.  Continuing this path would add state-machine risk around
output-time `dsmcN_` reduction without a measurable payoff.

## Revert state

The source change was reverted, the case control was restored, and the binary
was rebuilt.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_dsmcN_owned_reset_revert_build_20260607.log
```

Restored control hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

There is no remaining unstaged diff in `dsmcVolFields.C`,
`dsmcVolFields.H`, or the case `system/controlDict` from this candidate.
