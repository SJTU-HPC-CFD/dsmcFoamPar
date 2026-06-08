# MPI replicated mesh flat-transfer-false stream path rejection - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI only, `useOpenMP false`, `openmpThreads 1`,
  `mpirun -np 8`
- Baseline control hash:
  `8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2`
- Baseline fieldProperties hash:
  `73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9`

This was a configuration-only screen to test whether the delayed-receive /
async migration path can be activated without changing source code.  It did not
touch OpenMP, tracking, post/output, no-fields mode, or sampling semantics.

## Candidate

Temporary control change:

```text
replicatedMeshFlatTransfer false;
endTime 1.e-06;
```

Rationale:

- current accepted controls include `replicatedMeshDelayedReceive true`,
  `replicatedMeshNoAlltoall true`, and `replicatedMeshFlatTransfer true`;
- source review showed that `dsmcReplicatedMesh::migrateBegin()` immediately
  falls back to the flat owner-cell transfer path when `useFlatTransfer_` is
  true, so the delayed receive setting does not create a real overlap path in
  the accepted configuration;
- setting `replicatedMeshFlatTransfer false` routes into the stream fallback,
  which is the only existing configuration-level way to exercise the async
  begin/finish migration path.

Control backup before the candidate:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_flattransfer_false_smoke_20260607
sha256 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

## 10-step smoke

Smoke log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_flattransfer_false_smoke_10step_20260607
```

The run reached replicated-mesh initialization and did activate the intended
non-flat path:

```text
Replicated mesh: no-Alltoall async migration enabled
Replicated mesh: flatTransfer=0 found=1
```

It then entered initial particle distribution but did not reach the first time
iteration.  The last progress lines before termination were:

```text
Starting time loop

Replicated mesh: initial particle distribution

Replicated mesh: rank 0 kept 251365 parcels, deleted 1806661 non-owned parcels
```

Final log state:

| item | value |
| --- | --- |
| log size | 22212 bytes |
| created | 2026-06-07 12:43:57 CST |
| last modified | 2026-06-07 21:18:39 CST |
| final recheck | 2026-06-07 21:20:16 CST |
| solver completion markers | absent |
| termination markers | `BAD TERMINATION`, ranks 0-7 `KILLED BY SIGNAL: 9` |
| `/usr/bin/time` `real` line | `real 30882.88` |
| `/usr/bin/time` `user` line | `user 167245.59` |
| `/usr/bin/time` `sys` line | `sys 79811.70` |
| first timestep output | absent |

## Decision

Rejected.  Do not promote this candidate to 200-step or 500-step runs.

This candidate is not a slow-but-valid performance result.  It stalled before
time-step execution for more than eight hours while using the stream transfer
fallback, then terminated with SIGKILL on all ranks.  This matches the high
risk identified in the configuration screening note for the OFv1706 port.  The
current accepted flat-transfer path remains the only validated replicated-mesh
transfer path for this strict full-fields MPI8 case.

Implication for future work:

- do not spend more runs on `replicatedMeshFlatTransfer false` as a pure
  configuration switch;
- if delayed receive is still needed, implement a narrower source-level path
  that preserves the validated flat-transfer serialization and gives
  `migrateBegin()` / `migrateFinish()` a real overlap window, instead of
  falling back to the stream migration path.

## Revert state

The case control was restored to the accepted strict full-fields configuration:

```text
system/controlDict sha256
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2

system/fieldPropertiesDict sha256
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
```

Restored key controls:

```text
endTime 5.e-05;
replicatedMeshDelayedReceive true;
replicatedMeshFlatTransfer true;
replicatedMeshAutoDLB true;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
replicatedMeshDLBProfile false;
```

There is no source diff associated with this candidate.
