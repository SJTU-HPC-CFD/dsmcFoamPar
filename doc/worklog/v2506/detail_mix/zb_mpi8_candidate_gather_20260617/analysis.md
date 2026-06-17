# zb mpi8 allProcCandidates_ candidate-gather test, 2026-06-17

## Source Finding

`allProcCandidates_` is a cached per-rank candidate-count list in
`dsmcReplicatedMesh`.

Search result in the current tree:

- Written in `migrateParticlesByCellOwner()` after migration.
- Written in the non-flat migration path after migration.
- Exposed by `allProcCandidates()`.
- No current solver or replicated-mesh DLB path reads the accessor.

Therefore the per-migration `MPI_Allgather(localCands)` is not required by the
current `zb` replicated-mesh DLB path.

## Source Change

Added an opt-in control:

```text
replicatedMeshGatherCandidates false;
```

Default is `false`.

When disabled, the code skips the post-migration candidate-count loop and
`MPI_Allgather`.  The old behavior can be restored by setting:

```text
replicatedMeshGatherCandidates true;
```

The retained DLB controls, cell ownership, migration data exchange, particle
counts, and collision logic are unchanged.

## Validation Run

Case:

`run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/mpi8`

Base configuration:

`alltoall_guard50`

Additional explicit control:

```text
replicatedMeshGatherCandidates false;
```

Build result: passed.

Run result: exit `0`.

## Results

| Run | real [s] | execution@300 [s] | full evolve [s] | migration max [s] | sizeX max [s] | wait max [s] | candidate gather max [s] | DLB wall max [s] | DLB check max [s] | particles max/min | final energy | stuck |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| old alltoall_guard50 | 70.86 | 69.54 | 62.49474379 | 3.720826666 | 2.79109105 | 0.070676603 | 0.227600363 | 4.506035192 | 3.778185699 | 1.200627186 | 0.001959286903 | 0 |
| old alltoall_guard50_rep2 | 70.93 | 69.70 | 62.42188458 | 3.360279869 | 2.250693878 | 0.078226887 | 0.190361815 | 5.43824471 | 4.659641272 | 1.222443149 | 0.001953648664 | 0 |
| gatherCandidatesFalse | 76.51 | 75.44 | 61.84406319 | 3.702970187 | 2.79401046 | 0.072973595 | 0.000005455 | 6.298159618 | 5.550338799 | 1.241792915 | 0.001957306167 | 0 |

## Interpretation

The code change works mechanically: `migration candidate gather max` drops from
about `0.19-0.23 s` to approximately zero.

However, this single validation run did not improve end-to-end `real` time.
External time worsened to `76.51 s`, while `full evolve` improved to `61.84 s`.
The regression is therefore not explained by the candidate gather itself; this
sample had worse DLB/check timing and a larger `updateParticleCounts max`
(`0.261 s`) than the previous alltoall_guard50 repeats.

Conclusion:

- Keep `replicatedMeshGatherCandidates false` as a safe cleanup/overhead
  reduction for the current code path, because `allProcCandidates_` is unused.
- Do not claim an end-to-end speedup from this change based on one run.
- If this becomes the default for final numbers, repeat `alltoall_guard50` with
  the new binary at least twice more.


## Retest After Machine Fluctuation

A first rep2 attempt failed before solver startup with Intel MPI/Hydra:

```text
HYD_sock_listen_on_port: cannot open socket (Operation not permitted)
exit = 255
```

That failed run was caused by sandboxed MPI bootstrap restrictions and is not a
solver result.  It is excluded from performance comparison.

The same control was rerun outside the sandbox as
`alltoall_guard50_gatherCandidatesFalse_rep2_rerun` and completed with exit `0`.

| Run | real [s] | execution@300 [s] | full evolve [s] | candidate gather max [s] | updateParticleCounts max [s] | DLB wall max [s] | DLB check max [s] | particles max/min | final energy | stuck |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gatherCandidatesFalse rep1 | 76.51 | 75.44 | 61.84406319 | 0.000005455 | 0.261166333 | 6.298159618 | 5.550338799 | 1.241792915 | 0.001957306167 | 0 |
| gatherCandidatesFalse rep2 rerun | 68.72 | 67.73 | 55.5965783 | 0.000004354 | 0.233794993 | 5.497012967 | 4.794927666 | 1.234144307 | 0.001951684872 | 0 |

Updated interpretation:

- The candidate gather removal is mechanically stable: both valid runs reduce
  `migration candidate gather max` to approximately zero.
- The first valid run was affected by machine/runtime fluctuation; the retest is
  much faster and reaches `real=68.72 s`.
- Compared with old alltoall_guard50 repeats (`70.86/70.93 s`), the retest shows
  a plausible positive end-to-end result, but the valid two-run spread
  (`76.51` vs `68.72 s`) is too wide to claim a stable mean speedup yet.
- Keep the source change because `allProcCandidates_` is unused in the current
  path and the removed collective is semantically unnecessary.  For final claims,
  use more repeats under stable machine conditions.

## Artifacts

- Log: `logs/alltoall_guard50_gatherCandidatesFalse.log`
- Control snapshot:
  `controlDicts/alltoall_guard50_gatherCandidatesFalse_controlDict`
