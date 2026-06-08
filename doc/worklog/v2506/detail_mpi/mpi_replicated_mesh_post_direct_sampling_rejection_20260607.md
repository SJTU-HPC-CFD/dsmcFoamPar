# MPI replicated mesh direct dsmcVolFields sampling attempt - 2026-06-07

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
tracking, no-fields mode, output/reset, or shared-cache reuse across time steps.

## Candidate

Temporary source change:

- for replicated mesh, pure-MPI, non-density-only `dsmcVolFields` sampling, add
  a direct parcel-accumulation path when electronic, classification, and
  heat-flux/shear-stress measurements are disabled;
- accumulate the current single `mixture` field directly into the cumulative
  field arrays and species cumulative arrays;
- bypass the per-species shared sample cache and the subsequent field-combine
  pass for that restricted case;
- keep the existing shared-cache path as the fallback for other configurations;
- add a temporary profile-detail line for `direct parcel accum`.

Rationale:

- The current formal case has one `dsmcVolFields mixture` entry with all five
  species.
- A 10-step detail run showed the active field work was sampling-body work:
  `shared cache build 0.151396085 s` and `field combine 0.061793793 s`
  across 9 profiled calls.
- A direct single-field path could remove the combine pass without changing
  full-field output semantics.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_direct_sampling_build_20260607.log
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
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_post_direct_sampling_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_post_direct_sampling_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 9.22 s |
| move only | 1.274443693 s |
| buildCellOccupancy | 0.124734719 s |
| collision phase | 0.220984024 s |
| post fields/output | 0.836535452 s |
| full evolve wall | 2.49562147 s |
| dsmcVolFields sample accumulation | 0.199297762 s |
| dsmcVolFields direct parcel accum | 0.199266246 s |
| dsmcVolFields shared cache build | 0 s |
| dsmcVolFields cache allocate | 0 s |
| dsmcVolFields field combine | 0 s |
| dsmcVolFields cell reduction | 0.010563488 s |
| dsmcVolFields boundary accumulation | 0.00756991 s |
| OpenMP enabled | 0 |
| OpenMP max threads | 1 |
| final particles | 2065831 |
| stuck particles | 0 |
| final collisions | 2857 |
| final candidates | 2980 |
| final total energy | 1.15060919 |

Smoke correctness passed; no Fatal/NaN/BAD TERMINATION string was found.

## Decision

Rejected and reverted.  Do not promote this candidate to 200-step signal.

The field-internal `sample accumulation` improved only modestly versus the
recent strict full-field detail smoke:

| metric | strict detail smoke | direct sampling smoke | delta |
| --- | ---: | ---: | ---: |
| dsmcVolFields sample accumulation | 0.21323199 s | 0.199297762 s | -0.013934228 s |
| dsmcVolFields field combine | 0.061793793 s | 0 s | -0.061793793 s |
| dsmcVolFields shared cache build | 0.151396085 s | 0 s | -0.151396085 s |
| dsmcVolFields direct parcel accum | 0 s | 0.199266246 s | +0.199266246 s |
| post fields/output | 0.834544903 s | 0.836535452 s | +0.001990549 s |

The direct path removed the cache/combine stages, but the direct per-parcel
accumulation cost replaced nearly all of that work.  The total post timer did
not improve, so a 200-step signal run is not justified.

## Revert state

The source change was reverted, the case control was restored, and the binary
was rebuilt.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_post_direct_sampling_revert_build_20260607.log
```

Restored control hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

There is no remaining unstaged diff in `dsmcVolFields.C`,
`dsmcVolFields.H`, or the case `system/controlDict` from this candidate.
