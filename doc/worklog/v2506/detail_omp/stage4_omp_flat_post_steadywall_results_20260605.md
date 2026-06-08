# Stage 4 OMP Flat Post and Steady-Wall Profile Results

Date: 2026-06-05

Scope:
- OMP-only line.
- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8`
- Reference config/log source: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp/omp8`
- Reference code remained read-only: `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`

## Code Changes

1. Fixed the `dsmcVolFields.C` brace/scope error introduced while adding the flat occupancy sampling path.
2. Kept `dsmcVolFields::calculateField()` on a light OMP sampling path:
   - if flat occupancy ordered parcels are available, sample by cell in parallel;
   - each cell is written by one thread;
   - no atomics were introduced.
   - controlled by `openmpFieldSampling`, default `true`.
3. Replaced DSMC profile wall timers in `dsmcCloud.C` with `std::chrono::steady_clock`.
   - Reason: OFv1706 `clockTime` uses `gettimeofday()` and one 10-step run produced non-monotonic wall accounting (`post fields/output [s] < 0` and `full evolve wall < move+collide wall`).
   - CPU profile timers were left unchanged.

## Build

Command:

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
```

Status:
- build passed after the `dsmcVolFields.C` scope fix;
- build passed again after switching profile wall timers to `steady_clock`;
- build passed after adding the `openmpFieldSampling` control.

## Smoke/Profile Logs

All runs temporarily copied `ourmeshbkp/omp8/system/controlDict` into `ourmesh/omp8/system/controlDict` and restored the original file with an absolute-path trap.

Logs:
- `ourmesh/omp8/log.codex_omp_flat_post_steadywall_1step_20260605`
- `ourmesh/omp8/log.codex_omp_flat_post_steadywall_10step_20260605`
- `ourmesh/omp8/log.codex_omp_flat_post_steadywall_50step_20260605`
- `ourmesh/omp8/log.codex_omp_flat_post_steadywall_500step_20260605`
- `ourmesh/omp8/log.codex_omp_flat_post_defaultoff_steadywall_50step_20260605`
- `ourmesh/omp8/log.codex_omp_flat_post_defaultoff_steadywall_500step_20260605`

Smoke results:

| steps | real [s] | full evolve wall [s] | move [s] | buildCellOccupancy [s] | collision [s] | post [s] |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 5.06 | 0.289183826 | 0.149866546 | 0.010563336 | 0.069453426 | 0.058480048 |
| 10 | 5.95 | 3.023313479 | 1.540721655 | 0.095197040 | 0.739597923 | 0.637757900 |
| 50 | 21.56 | 16.77391922 | 8.297439705 | 0.531471265 | 4.336554933 | 3.556509818 |

The steady-wall profile is internally consistent in these runs: `full evolve wall` is greater than the directly accounted phase sum, and `post fields/output` is non-negative.

## Formal 500-Step Result

Current stage4 OMP8:

| metric | value |
| --- | ---: |
| real | 265.13 s |
| user | 2044.67 s |
| sys | 40.99 s |
| ClockTime | 265 s |
| full evolve wall | 272.2832057 s |
| move only | 144.8259732 s |
| buildCellOccupancy | 7.89583818 s |
| collision phase | 68.18518398 s |
| post fields/output | 50.78234158 s |
| collisions | 35753 |
| particles | 2463724 |
| total energy | 1.24318294 |

Diagnostic run with `openmpFieldSampling false`:

| metric | value |
| --- | ---: |
| real | 277.45 s |
| user | 1512.14 s |
| sys | 29.58 s |
| ClockTime | 277 s |
| full evolve wall | 285.4980434 s |
| move only | 136.1531591 s |
| buildCellOccupancy | 8.376696021 s |
| collision phase | 69.16212251 s |
| post fields/output | 71.3218158 s |
| collisions | 35290 |
| particles | 2463691 |
| total energy | 1.243336011 |

This shows that flat field sampling is a local improvement in the current code:
- `post fields/output`: `71.32 s -> 50.78 s`;
- `real`: `277.45 s -> 265.13 s`.

Previous stage3 current-tree OMP8 formal result:

| metric | value |
| --- | ---: |
| real | 248.80 s |
| ClockTime | 249 s |
| collisions | 35283 |
| particles | 2463725 |
| total energy | 1.243229348 |

Note: stage3 phase labels used the older timing path and are not reliable as wall-profile values. Use `real`/`ClockTime` for end-to-end comparison.

Historical `ourmeshbkp/omp8` result:

| metric | value |
| --- | ---: |
| external wall | 103.66 s |
| ClockTime | 104 s |
| full evolve wall | 102.1225033 s |
| move only | 48.82365912 s |
| buildCellOccupancy | 6.490896398 s |
| collision phase | 9.876225069 s |
| post fields/output | 36.92670834 s |
| collisions | 35847 |
| particles | 2463391 |

## Assessment

Correctness is close enough for this stage:
- final particle count matches the previous current-tree OMP run within 1 parcel;
- final total energy differs from stage3 by about `4.64e-05`;
- final collision count is closer to the historical OMP8 count than stage3 was, but still stochastic.

Performance is not acceptable as an improvement:
- stage4 default `real = 265.13 s`, slower than stage3 `real = 248.80 s`;
- historical OMP8 remains much faster at about `103.66 s`;
- `buildCellOccupancy` is now near the historical scale, so it is no longer the main OMP gap;
- the remaining formal 500-step gaps are dominated by `move`, `collision`, and then `post`.

Conclusion:
- Keep the steady-wall profile fix.
- Keep flat field sampling enabled by default, because it is faster than the serial field sampling path in the current tree.
- Do not claim overall OMP success yet: the current default is still slower than stage3 and much slower than `ourmeshbkp/omp8`.
- The next OMP work should focus on why current `move only` and `collision phase` are far slower than `ourmeshbkp/omp8`, before attempting broader collision/post refactors.
