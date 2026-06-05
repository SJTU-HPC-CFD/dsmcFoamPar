# Stage 0 Full Performance Profile - 2026-06-05

## Inputs

Baseline logs parsed from:

- `ourmeshbkp/omp8/log.omp8.confirm_pdFalse_20260605_024610`
- `ourmeshbkp/mpi8origin/log.mpi8origin.confirm_pdFalse_20260605_024610`
- `ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610`

All three logs are 500-step runs to `endTime 5.e-05`, `deltaT 1.e-07`. The case `controlDict` files in `ourmeshbkp` do not currently show `profileSummary/profileDetail`; the logs nevertheless contain solver profile summaries, so they came from a profiling-capable source revision and should be treated as historical benchmark logs.

## Clean Baseline Table

| Mode | External final ClockTime | Full evolve wall | Move only | buildCellOccupancy | Collision phase | Post fields/output | Auto DLB/rebalance | Final particles / total particles | Last-step collisions | Energy check |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `omp8` | 104 s | 102.1225033 s | 48.82365912 s | 6.490896398 s | 9.876225069 s | 36.92670834 s | 0 s | 2,463,391 | not reported | avg total energy `1.032208342e-18` |
| `mpi8origin` | 160 s | 163.0536101 s | 120.3658396 s | 6.924980537 s | 18.73261305 s | 17.02521609 s | 0 s | 2,463,429 | not reported | avg total energy `1.031952906e-18` |
| `mpi8replicatedmesh` | 117 s | 114.9207628 s | 47.1652636 s | 10.7084157 s | 10.89824478 s | 30.35563996 s | 15.78659015 s | total 2,466,342 | 35,158 | per-rank average energy ranges `8.338e-19` to `1.140e-18` |

No `FOAM FATAL`, `MPI_ABORT`, or lowercase `nan` markers were found in the parsed tails.

## Time Share By Mode

Percentages below use `full evolve wall` as denominator.

| Mode | Move | buildCellOccupancy | Collision | Post fields/output | Auto DLB |
|---|---:|---:|---:|---:|---:|
| `omp8` | 47.81% | 6.36% | 9.67% | 36.16% | 0.00% |
| `mpi8origin` | 73.82% | 4.25% | 11.49% | 10.44% | 0.00% |
| `mpi8replicatedmesh` | 41.04% | 9.32% | 9.48% | 26.42% | 13.74% |

## Hotspot Interpretation

### OMP8

`omp8` is move-heavy and post-field-heavy:

- move only: 48.82 s
- post fields/output: 36.93 s
- collision: 9.88 s
- build occupancy: 6.49 s

The immediate Amdahl ceiling for collision-only work is low because collision is under 10% of full evolve. Stage 1 should prioritize move data flow and post-field/occupancy reuse before deeper collision tuning.

### MPI Origin

`mpi8origin` is dominated by ordinary MPI move transfer/finalize:

- move only: 120.37 s
- move transfer/delete finalize: 63.64 s inside move profiling
- move parallel kernel wall: 56.70 s

This validates replicated mesh / ownership migration as the main MPI direction. Collision tuning alone cannot fix the main critical path.

### MPI Replicated Mesh

`mpi8replicatedmesh` removes the ordinary transfer/finalize hotspot:

- move transfer/delete finalize: 0 s
- move only: 47.17 s
- auto DLB/rebalance: 15.79 s
- Phase C auto DLB rebalances: 8
- particles per rank final: min 228,405, max 405,542, max/min 1.7755

The replicated path is faster than `mpi8origin` by full evolve wall, but it pays a visible DLB tax. Later Stage 3 work should keep DLB trigger/gap tuning separate from movement and output correctness.

## Correctness And Consistency Notes

- `omp8` and `mpi8origin` final particle counts agree closely: 2,463,391 versus 2,463,429.
- `mpi8replicatedmesh` reports total particles 2,466,342, within about 0.12% of `mpi8origin`; this is close enough for a baseline comparison but must be rechecked after every port.
- `omp8` and `mpi8origin` final average total energy are close: `1.032208342e-18` versus `1.031952906e-18`.
- Replicated mesh energy is printed per-rank in the parsed tail, not as one global average; Stage 3 output fixes must verify global reduced fields and derived temperatures.

## Stage 0 Gate

Stage 0 passes for analysis:

- compile check succeeded;
- historical baseline logs are parseable;
- performance bottlenecks map to concrete stage priorities.

Stage 0 does not yet pass for implementation equivalence:

- current source does not contain the profiling/OpenMP/replicated mesh implementation that produced the historical logs;
- the next source edit should start with minimal profiling instrumentation so fresh OFv1706 runs can be compared against the historical baseline.
