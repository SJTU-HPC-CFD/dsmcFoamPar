# MPI8 replicated-mesh subphase profiling, 2026-06-12

## Scope

Implemented low-perturbation subphase profiling for the current MPI8
replicated-mesh path.  The instrumentation only accumulates local wall-time
with `std::chrono` during the timestep.  MPI aggregation is done in final
summary/report output, so it does not add per-step collectives.

## Source changes

- `noTimeCounter`: collision subphase timers are now active when
  `profileSummary true`, not only when `profileDetail true`.
- `collisionPartnerSelection`: added virtual subphase timing accessors so
  `dsmcCloud::printProfileSummary()` can report them without depending on the
  concrete model.
- `dsmcReplicatedMesh`: split the active flat POD migration path into:
  pack, local prep, size exchange, request post, wait, deserialize, post, and
  candidate gather.
- `dsmcReplicatedMesh::updateParticleCounts()` now reports its own final max
  wall time.
- `dsmcReplicatedMesh::report()` now does one final `MPI_Allreduce` for the
  migration/update-count timing maxima.

## Build

- `wmake libso src/lagrangian/dsmc`: pass.
- `wmake applications/solvers/discreteMethods/dsmc/dsmcFoam+`: pass.

## Smoke run

Case copy:

```text
/tmp/hystrath_zb_mpi8_subphase_tJ9NCe
```

Run:

```text
OMP_NUM_THREADS=1 mpirun -np 8 dsmcFoam+
```

Log files:

- `logs/zb_mpi8_subphase_smoke2_20260612.log`
- `logs/zb_mpi8_subphase_smoke2_20260612.exit`
- `logs/zb_mpi8_subphase_smoke2_20260612.time`

Exit code: `0`.

Key 2-step smoke metrics:

```text
collision phase [s]           = 0.06820428
collision localLoop max [s]   = 0.017023999
collision reduce max [s]      = 0.061179396
collision sigmaBC max [s]     = 0.002096955
collision accounted max [s]   = 0.068177392

migration wall max [s]        = 0.057154865
migration size exchange max [s] = 0.050879434
migration wait max [s]        = 0.000699068
migration candidate gather max [s] = 0.000827159
updateParticleCounts max [s]  = 0.02856863
```

## Immediate read

This smoke is not a formal performance result, but it confirms the suspected
shape: in MPI8 replicated mesh, the collision stage is not dominated by the
local collision kernel on this short run.  The collision reduce path dominates
the reported collision subphases.  Separately, regular replicated-mesh
migration exposes a large size-exchange synchronization cost, and
`updateParticleCounts()` is also visible as a standalone collective cost.

The next useful measurement is the normal 500-step `ourmesh/mpi8replicatedmesh`
case with these summary timers enabled and `profileDetail false`, so the final
report can separate collision local work from reduction/sigmaBC/migration/count
collective time without per-step diagnostic distortion.
