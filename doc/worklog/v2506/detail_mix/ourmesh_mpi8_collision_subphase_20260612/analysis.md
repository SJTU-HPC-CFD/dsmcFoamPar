# ourmesh MPI8 collision subphase / synchronization analysis

Date: 2026-06-12

## Scope

- Case: `pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`.
- Formal case length: 500 steps, `endTime 5.e-05`, `deltaT 1.e-07`, no writes.
- Source state: dirty local tree with owned collision-cell iterator, cached serial
  `noTimeCounter` hot-path work, and this worklog's extra collision subphase
  instrumentation.
- Default path after this note: `collisionReduceOnlyOnOutput false`, so normal
  production behavior still performs the per-step global collision count reductions.

## Instrumentation added

Files:

- `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.H`
- `src/lagrangian/dsmc/collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C`

Added `profileDetail true` diagnostics for `noTimeCounter::collide()`:

- cumulative local active-cell/kernel loop wall time;
- cumulative global collision-count reduction wall time;
- cumulative `sigmaTcRMax.correctBoundaryConditions()` wall time;
- cumulative pre-output collision wall time by rank.

Also added an experimental control:

```text
collisionReduceOnlyOnOutput false;
```

When set true, global `Collisions` / `Collision candidates` reductions are done
only at terminal-output steps. This was tested below and kept default-off because
it improves the collision metric but not total time.

## Existing reference numbers

From `serial1_coll_speedup_20260611/summary.md` and
`mpi8_collision_cached_refs_20260611/summary.md`:

| mode | real [s] | full evolve [s] | move [s] | build [s] | collision [s] | post [s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Serial1 | 283.07 | 302.392635 | 258.861539 | 27.853350 | 14.636285 | 0.633565 |
| OMP8 current | 59.21 | 57.779304 | 47.208625 | 6.397944 | 3.003647 | 0.639841 |
| MPI8 cached-ref mean | 92.69 | 92.260071 | 70.295910 | 10.200031 | 25.957368 | 4.945385 |

This is the symptom: serial collision is `14.64 s`, while MPI8 reported
collision is `25.96 s`.

## Detail run, DLB profile off

Run:

- log: `logs/ourmesh_MPI8_collision_subphase_nodlbprofile_500step_20260612.log`
- time: `logs/ourmesh_MPI8_collision_subphase_nodlbprofile_500step_20260612.time`
- exit: `0`
- config: `profileDetail true`, `nTerminalOutputs 500`,
  `replicatedMeshDLBProfile false`.

Summary:

| metric | value |
| --- | ---: |
| external real [s] | 107.60 |
| full evolve wall [s] | 109.399549 |
| move+collide wall [s] | 103.936834 |
| move only [s] | 80.998946 |
| buildCellOccupancy [s] | 12.079373 |
| collision phase [s] | 34.279284 |
| post fields/output [s] | 5.496350 |
| collisions / candidates | 33289 / 56371 |
| stuck particles | 0 |
| total energy | 1.242001602 |

Rank/subphase evidence:

| rank | final parcels | final candidates | move [s] | build [s] | collision [s] | localLoop [s] | reduce [s] | sigmaBC [s] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 323108 | 729 | 55.039608 | 8.832424 | 34.279284 | 1.213284 | 32.396364 | 0.662951 |
| 1 | 325620 | 6256 | 75.112118 | 10.215834 | 14.853831 | 2.949352 | 11.255324 | 0.643505 |
| 2 | 310510 | 776 | 64.844325 | 9.370426 | 24.964066 | 2.094887 | 22.198543 | 0.665100 |
| 3 | 299875 | 3615 | 61.039195 | 9.613316 | 28.238396 | 3.717052 | 23.822533 | 0.693292 |
| 4 | 284368 | 12337 | 63.137356 | 10.279052 | 25.605254 | 1.341336 | 23.600391 | 0.656730 |
| 5 | 298907 | 548 | 80.998946 | 12.079373 | 7.611305 | 2.316645 | 4.610417 | 0.679141 |
| 6 | 309339 | 3375 | 72.838439 | 10.484291 | 16.708375 | 1.790760 | 14.238271 | 0.674084 |
| 7 | 312016 | 28735 | 80.805912 | 10.057366 | 9.754863 | 5.285824 | 3.811512 | 0.652033 |

Key ratios:

- max `localLoop` = `5.285824 s`
- max `reduce` = `32.396364 s`
- max `collision phase` = `34.279284 s`
- final-window candidates max/min = `52.44`
- cumulative candidates max/min = `15.93`
- rank collision max/min = `4.50`
- rank full max/min = `1.00007`

Interpretation:

The MPI8 collision phase is not mostly collision kernel time. The slowest
measured local active-cell loop is only `5.29 s`, already faster than the
single-core serial collision time `14.64 s`. The reported MPI collision phase is
large because fast ranks enter `noTimeCounter` reductions early and wait for
ranks that are still in move/build/local collision work. That wait is charged to
`collision phase`.

## Experimental defer-reduce run

This run used an intermediate source state equivalent to:

```text
collisionReduceOnlyOnOutput true;
```

Run:

- log: `logs/ourmesh_MPI8_deferred_reduce_500step_20260612.log`
- time: `logs/ourmesh_MPI8_deferred_reduce_500step_20260612.time`
- exit: `0`

Result:

| metric | cached-ref mean | defer-reduce run |
| --- | ---: | ---: |
| real [s] | 92.69 | 107.01 |
| full evolve [s] | 92.260071 | 109.027011 |
| collision phase [s] | 25.957368 | 5.300807 |
| post fields/output [s] | 4.945385 | 1.474260 |
| Phase C DLB check max [s] | ~0.036 | 31.307800 |
| Phase C DLB rebalances | 3-4 | 8 |

This proves the direct source of the large collision metric: removing the
per-step collision reductions collapses `collision phase` to about `5.3 s`.
However, total time does not improve. The same synchronization wait moves to the
next collective DLB check/allgather, so this cannot be kept as a default
optimization.

## Experimental defer-reduce plus no-DLB run

Run:

- log: `logs/ourmesh_MPI8_deferred_reduce_nodlb_500step_20260612.log`
- time: `logs/ourmesh_MPI8_deferred_reduce_nodlb_500step_20260612.time`
- exit: `0`
- config: equivalent to `collisionReduceOnlyOnOutput true` plus
  `replicatedMeshAutoDLB false`.

Result:

| metric | value |
| --- | ---: |
| real [s] | 104.71 |
| full evolve wall [s] | 103.808698 |
| move+collide wall [s] | 102.801386 |
| move only [s] | 83.100710 |
| buildCellOccupancy [s] | 14.932237 |
| collision phase [s] | 5.369201 |
| post fields/output [s] | 1.063522 |
| migration wall [s] | 17.776834 |

Disabling DLB avoids the huge DLB-check bucket, but total time is still not
better than the cached-ref MPI8 mean. Waiting and imbalance move into
move/build/migration instead.

## Conclusion

The user's observation is correct: on `ourmesh`, "serial coll is small while MPI
coll is large" means extra overhead is being introduced. The detailed profile
shows that the extra overhead is not extra useful collision computation. It is
blocking synchronization/arrival skew being attributed to collision.

More precisely:

1. MPI8 local collision kernel work is about `5.3 s` max-rank in the detail run,
   not `25-35 s`.
2. The large MPI8 `collision phase` comes from collectives inside
   `noTimeCounter::collide()`, mainly global collision/candidate reductions.
3. Removing those reductions from every step improves the collision metric but
   not total time, because the wait transfers to DLB or migration collectives.
4. OMP8 avoids this failure mode because it runs one process with dynamic
   shared-memory scheduling over active collision cells; it does not have
   inter-rank phase barriers inside collision.

Therefore the remaining MPI8 problem is not just "optimize the serial collision
kernel". That part is already much smaller than the reported collision phase.
The real target is to reduce or overlap global synchronization and to make the
owned-cell partition phase-aware:

- keep the active owned-cell iterator and cached serial kernel;
- use `profileDetail true` to compare `localLoop` vs `reduce` before trusting
  any collision-phase speedup;
- do not default-enable `collisionReduceOnlyOnOutput` until a total-time run
  improves;
- next code target should be DLB/migration synchronization: either reduce the
  frequency of blocking collective checks without moving the wait elsewhere, or
  redesign the DLB score to balance move/build/collision phase timing and
  collision candidates separately;
- hybrid MPI+OMP remains the practical way to reproduce OMP's dynamic
  collision-cell scheduling inside each owned MPI partition.

## Verification

- Build: `source doc/scripts/env.sh && source doc/scripts/build-dsmcFoam.sh`
  passed after adding the instrumentation and default-off switch.
- All three runs in this worklog exited `0`.
- Fatal-marker scans found no `FOAM FATAL`, `Segmentation`, `MPI_ABORT`,
  `Killed`, `ERROR`, or `Error`.
- Final source-state smoke:
  `logs/ourmesh_MPI8_final_default_detail_smoke10_20260612.log`, 10 steps,
  `profileDetail true`, default `collisionReduceOnlyOnOutput false`, exit `0`.
  It printed the new collision subphase block and ended with `0` stuck
  particles.
- The original
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/system/controlDict`
  was restored and verified byte-identical to
  `controlDicts/ourmesh_mpi8_subphase_controlDict_before_20260612` after a
  hardlink-copy test temporarily changed the source case controls.
