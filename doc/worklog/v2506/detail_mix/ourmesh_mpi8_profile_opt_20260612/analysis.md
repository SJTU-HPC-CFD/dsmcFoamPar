# ourmesh MPI8 replicated-mesh collision/migration communication optimization

Date: 2026-06-12

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Formal length: 500 steps, `endTime 5.e-05`, `deltaT 1.e-07`, no field writes
  in this interval.
- Run mode: `OMP_NUM_THREADS=1 mpirun -np 8 dsmcFoam+ -case <case>`,
  raw MPI replicated mesh, no `-parallel`.
- Control state: `replicatedMesh true`, `replicatedMeshNoAlltoall true`,
  `replicatedMeshFlatTransfer true`, `replicatedMeshAutoDLB true`,
  `replicatedMeshDLBProfile false`.

## Code changes

1. Flat POD migration no longer uses the per-step size `MPI_Alltoall` when
   `replicatedMeshNoAlltoall true`.
   Each peer sends exactly one message, including 0-byte messages; receivers use
   `MPI_Probe`/`MPI_Get_count` to discover the actual byte count and then
   receive into per-peer buffers. This removes the fixed size-exchange collective
   from the flat transfer path.

2. Replicated-mesh `noTimeCounter` now defaults to avoiding global collision
   count reductions during normal terminal diagnostics:
   `collisionReduceOnlyOnOutput` defaults to true and
   `collisionOutputGlobalReduce` defaults to false for replicated mesh.
   The terminal `Collisions` / `Collision candidates` lines are explicitly marked
   as rank-0 local diagnostics when this path is active.

3. `printProfileSummary()` now does one final global sum for cumulative
   collisions/candidates. This keeps a comparable final correctness/profile
   counter without reintroducing the repeated output-time collision reduction.

4. The profile summary now reports collision subphases and migration subphases:
   local collision loop, collision reduce, sigma boundary update, migration
   size exchange, request post, wait, deserialize, candidate gather, and DLB
   check/repartition buckets.

## Build

```bash
source doc/scripts/build-dsmcFoam.sh
```

Build passed after the final summary-counter change.

## 500-step results

Historical reference from
`doc/worklog/v2506/detail_mix/mpi8_collision_cached_refs_20260611/summary.md`:

| state | real [s] | full evolve [s] | move [s] | build [s] | coll [s] | post [s] |
|---|---:|---:|---:|---:|---:|---:|
| cached-ref mean | 92.69 | 92.260071 | 70.295910 | 10.200031 | 25.957368 | 4.945385 |

This worklog:

| state | real [s] | full evolve [s] | move [s] | build [s] | coll [s] | post [s] | coll reduce max [s] | mig sizeX max [s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| before opt | 100.26 | 89.482993 | 67.365519 | 10.362814 | 26.988595 | 4.763223 | 25.456105 | 2.747550 |
| flat no-Alltoall | 96.33 | 93.838338 | 71.939852 | 10.192066 | 24.341288 | 5.106053 | 22.669874 | 0 |
| no-Alltoall + reduce still active | 91.10 | 80.618083 | 62.016005 | 10.060068 | 26.621681 | 1.151625 | 25.559251 | 0 |
| output collision reduce skipped | 87.34 | 83.307570 | 64.395523 | 9.972819 | 4.783962 | 1.279978 | 0 | 0 |
| final source + final global summary | 90.79 | 85.164136 | 67.539543 | 10.527126 | 5.409990 | 0.956470 | 0 | 0 |

Final source log:

- log:
  `logs/ourmesh_mpi8_profile500_final_summary_20260612.log`
- time:
  `logs/ourmesh_mpi8_profile500_final_summary_20260612.time`
- exit: `0`
- fatal-marker scan: no `FOAM FATAL`, `MPI_ABORT`, `Segmentation`, `Killed`,
  `ERROR`, or `Error`.
- final particles: `2,463,713`
- stuck particles: `0`
- total energy: `1.242486044`
- cumulative global collisions/candidates:
  `9,672,441 / 17,592,071`

## Interpretation

The collision metric is now genuinely accelerated:

- cached-ref collision phase: `25.957368 s`
- final-source collision phase: `5.409990 s`
- improvement: about `79.2%`

The end-to-end 500-step wall time also improves, but only modestly:

- cached-ref mean external real: `92.69 s`
- final-source external real: `90.79 s`
- improvement: about `2.0%`

The fastest observed source state in this directory was
`after_output_reduce_skip` at `87.34 s`, but the final rerun landed at `90.79 s`
because move/DLB stochastic balance changed. Both runs have `collision reduce max
= 0` and `migration size exchange max = 0`; the remaining variation is not in
the collision kernel.

The flat migration change did remove the intended fixed collective:

- before opt migration size exchange max: `2.747550 s`
- final-source migration size exchange max: `0 s`

However, the same migration step still has arrival-skew wait:

- final-source migration wait max: `2.519002 s`
- final-source migration candidate gather max: `0.136097 s`

So this optimization removes one fixed collective, but not the phase skew created
by unequal move/build work before ranks reach migration.

## Why the total speedup is limited

After removing collision reductions, the local collision loop is no longer the
main problem. In the final 500-step run:

- collision local loop max: `5.260910 s`
- collision reduce max: `0 s`
- move only max: `67.539543 s`
- buildCellOccupancy max: `10.527126 s`
- rank evolve wall max/min: `1.422620`

That is the key reason MPI8 cannot reproduce OMP8's collision-dominated speedup.
OMP has one process and shared-memory scheduling over active cells; MPI still has
phase boundaries where faster ranks wait for slower ranks after move/build and
before collective or point-to-point synchronization. Removing a collective often
moves the wait to the next synchronization site unless the underlying per-rank
move/build skew is reduced.

## Current status

This source state satisfies the immediate communication target:

- flat POD migration no longer performs the size `MPI_Alltoall`;
- per-output collision-count global reductions are removed from the default
  replicated-mesh path;
- final global collision/candidate totals are retained with one summary-time
  reduction;
- 500-step `mpi8 replicatedMesh` shows both collision-phase speedup and a real
  total-time improvement versus the cached-ref mean.

The next profitable target is not more collision-kernel work. It is reducing
move/build imbalance and the arrival-skew wait that appears in migration/DLB.
