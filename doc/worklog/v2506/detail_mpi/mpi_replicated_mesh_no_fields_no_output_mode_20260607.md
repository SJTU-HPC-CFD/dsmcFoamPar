# MPI replicated mesh no-fields/no-output mode - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI only, `useOpenMP false`, `openmpThreads 1`,
  `mpirun -np 8`
- Source changes: none.
- Temporary config candidate: clear `system/fieldPropertiesDict` to
  `dsmcFields ();`.

This is not a like-for-like full-field configuration.  It is valid only as a
no-field-measurement/no-field-output performance mode.

## Applicability

Current formal controls do not write field output:

```text
endTime 5.e-05;
writeInterval 1.e-3;
```

The current collision setup also does not require the macroscopic field
temperature path:

- `system/chemReactDict` has `reactions ();`;
- `constant/dsmcProperties` uses `LarsenBorgnakkeVariableHardSphere`;
- no `inverseZvFormulation` is specified, so the model default remains
  `invZvFormulation = 2`, which does not read `fields().overallT(cellI)`.

Therefore this mode is a legitimate performance screen when field sampling data
are not required.  It is not acceptable for runs that need `dsmcVolFields`
sampling/output or `inverseZvFormulation 2008`.

## Backups

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_no_fields_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2

doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_fieldPropertiesDict_before_no_fields_smoke_20260607
sha256 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
```

## 10-step smoke

Temporary config:

```text
endTime 1.e-06;
dsmcFields ();
profileSummary true;
profileDetail false;
useOpenMP false;
openmpThreads 1;
```

Log:

```text
log.codex_mpi8_replicatedmesh_no_fields_smoke_10step_20260607
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| real | 11.30 s |
| move only | 1.112012558 s |
| buildCellOccupancy | 0.108655299 s |
| collision phase | 0.199578 s |
| post fields/output | 0.635546126 s |
| full evolve wall | 2.302075802 s |
| OpenMP enabled | 0 |
| DLB rebalances | 0 |
| final particles | 2065837 |
| stuck particles | 0 |
| final collisions | 2829 |
| final candidates | 2925 |
| final total energy | 1.150613464 |

Smoke passed; no Fatal/NaN/BAD TERMINATION string was found.

## 200-step signal

Temporary config:

```text
endTime 2.e-05;
dsmcFields ();
profileSummary true;
profileDetail false;
useOpenMP false;
openmpThreads 1;
```

Log:

```text
log.codex_mpi8_replicatedmesh_no_fields_signal_200step_20260607
```

Result:

| metric | no-fields 200-step |
| --- | ---: |
| Total Iterations | 200 |
| real | 49.79 s |
| move only | 22.43553799 s |
| move+collide wall | 30.66654133 s |
| buildCellOccupancy | 2.532010612 s |
| collision phase | 5.041371943 s |
| post fields/output | 12.3229034 s |
| full evolve wall | 42.9284654 s |
| migration wall time | 0.666515448 s |
| DLB checks | 200 |
| DLB rebalances | 2 |
| particles max/min | 1.848781249 |
| rank wall max/min | 1.009292352 |
| OpenMP enabled | 0 |
| final particles | 2220259 |
| stuck particles | 0 |
| final collisions | 15248 |
| final candidates | 25438 |
| final total energy | 1.169688813 |

Comparison against the full-field 200-step signal from the final-reset-skip
attempt (`log.codex_mpi8_replicatedmesh_post_final_reset_skip_signal_200step_20260607`):

| metric | full fields | no fields | delta |
| --- | ---: | ---: | ---: |
| real | 54.36 s | 49.79 s | -4.57 s |
| move only | 26.41077164 s | 22.43553799 s | -3.97423365 s |
| post fields/output | 16.33986782 s | 12.3229034 s | -4.01696442 s |
| full evolve wall | 48.47341764 s | 42.9284654 s | -5.54495224 s |

The 200-step signal is strong enough for 500-step formal validation.

## 500-step formal

Temporary config:

```text
endTime 5.e-05;
dsmcFields ();
profileSummary true;
profileDetail false;
useOpenMP false;
openmpThreads 1;
```

Log:

```text
log.codex_mpi8_replicatedmesh_no_fields_500step_20260607
```

Result against the current strict full-field baseline
`log.codex_mpi8_replicatedmesh_retained_forced8_repeat_500step_20260607`:

| metric | strict full fields | no fields | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 119.69 s | -14.14 s |
| move only | 72.9567669 s | 64.72452806 s | -8.23223884 s |
| move+collide wall | 97.10060854 s | 88.04482464 s | -9.05578390 s |
| buildCellOccupancy | 6.729516136 s | 6.855826684 s | +0.126310548 s |
| collision phase | 14.8398869 s | 14.1775748 s | -0.66231210 s |
| post fields/output | 40.37082178 s | 30.60678571 s | -9.76403607 s |
| full evolve wall | 132.7056254 s | 118.6090831 s | -14.09654230 s |
| migration wall time | 2.481249298 s | 2.055102923 s | -0.426146375 s |
| DLB checks | 500 | 500 | 0 |
| DLB rebalances | 8 | 8 | 0 |
| Phase C wall max | 2.292849742 s | 1.840267165 s | -0.452582577 s |
| Phase C migration max | 1.687160342 s | 1.240092133 s | -0.447068209 s |
| particles max/min | 1.412911358 | 2.676803069 | +1.263891711 |
| rank wall max/min | 1.053059781 | 1.003844246 | -0.049215535 |

Correctness:

- `Total Iterations = 500`;
- `OpenMP enabled = 0`, `OpenMP max threads = 1`;
- final particles `2463718`, stuck `0`;
- final collisions `34446`, candidates `61940`, acceptance `0.5561188247`;
- final total energy `1.2393756`;
- no Fatal/NaN/BAD TERMINATION string was found.

Historical-reference note:

- The preserved historical reference reports `External wall seconds 117.31`.
- This no-fields run reaches `real 119.69`, but the comparison is not
  equivalent because field sampling is disabled here.

## Decision

Retain this as a documented no-output performance mode, not as the default
strict full-field configuration.

- If field measurements/output are not needed for a timing campaign, clearing
  `dsmcFields` saves about `14.14 s` versus the current same-tree strict
  full-field forced8 baseline.
- If field output or macroscopic field temperature coupling is needed, keep the
  original `dsmcVolFields` configuration.
- The default case files were restored after the test.

Restored hashes:

```text
controlDict         8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
fieldPropertiesDict 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
```
