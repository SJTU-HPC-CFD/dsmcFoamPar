# Stage6 OMP Full Move Port Results - 2026-06-05

## Scope

This stage continued only on the OMP line. It ported the retained reference move outer-control layer that is compatible with the OFv1706 particle-tracking API:

- `Cloud<ParticleType>::move()` now has a reference-style OMP extract/kernel/commit path;
- the OMP move path can reuse `moveOrderedParcels()` plus appended parcels instead of always scanning the linked list;
- survivor/delete commit is counted and assembled per thread before serial deletion;
- `dsmcCloud` now has move append capture support:
  - `beginMoveAppendCapture()`;
  - `endMoveAppendCapture()`;
  - `recordMoveAppendedParcel()`;
  - `moveAppendedParcels()`;
- `addNewParcel()` and `addNewStuckParcel()` record appended parcels while capture is active;
- `rebuildMoveOrderedParcels()` and `setMoveOrderedParcels()` now generate `ompNumThreads + 1` offsets instead of a single 2-entry range.

The reference `useMoveParticlePartition` branch is still hardcoded false in the reference `Cloud.C`, so this stage did not turn it into a new default kernel. The v2506 particle core uses a different `trackToAndHitFace()` API and was not copied over the OFv1706 `trackToFace()` implementation.

## Build

Command:

```bash
source /home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb/doc/scripts/build-dsmcFoam.sh
```

Result:

- build completed;
- `dsmcFoam+` rebuilt at `platforms/linux64IccDPInt32Opt/bin/dsmcFoam+`;
- `dsmcInitialise+` rebuilt at `platforms/linux64IccDPInt32Opt/bin/dsmcInitialise+`.

## 10-step OMP Smoke

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_move_fullport_smoke_10step_omp_on_20260605`

Temporary control settings:

- copied OMP/profiling settings into `system/controlDict`;
- changed `endTime` to `1.e-06`;
- restored the original `controlDict` after the run.

Result:

| metric | value |
| --- | ---: |
| iterations | 10 |
| move only | 1.10950086 s |
| buildCellOccupancy | 0.10550111 s |
| collision phase | 0.019319499 s |
| full evolve wall | 1.932411347 s |
| collisions | 2860 |
| collision candidates | 2980 |
| particles | 2065869 |
| stuck particles | 0 |
| total energy | 1.150622572 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 1`

No `Fatal`, `NaN`, or `Segmentation` entries were reported.

## 500-step Formal Run

Formal comparison must use the same control settings as stage5. For this run, `ourmeshbkp/omp8/system/controlDict` was temporarily copied to `ourmesh/omp8/system/controlDict`, then the original `ourmesh/omp8/system/controlDict` was restored.

Log:

- `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_move_fullport_bkpctrl_500step_20260605`

Result:

| metric | value |
| --- | ---: |
| real | 161.79 s |
| user | 1250.08 s |
| sys | 20.79 s |
| ClockTime | 161 s |
| full evolve wall | 167.0823348 s |
| move+collide wall | 115.3128968 s |
| move only | 103.6495219 s |
| buildCellOccupancy | 7.483010573 s |
| collision phase | 3.613524884 s |
| post fields/output | 51.76373504 s |
| move+collide cpu | 860.04 s |
| move only cpu | 772.92 s |
| buildCellOccupancy cpu | 55.05 s |
| collision phase cpu | 27.09 s |
| post fields/output cpu | 387.05 s |
| full evolve cpu | 1247.10 s |
| collisions | 35898 |
| collision candidates | 68794 |
| collision acceptance rate | 0.5218187633 |
| particles | 2463609 |
| stuck particles | 0 |
| total energy | 1.242916724 |

OpenMP settings reported by the log:

- `OpenMP enabled = 1`
- `OpenMP max threads = 8`
- `OpenMP move = 1 (static, chunk 64)`
- `OpenMP collision = dynamic, chunk 8`

No `Fatal`, `NaN`, or `Segmentation` entries were reported.

## Stage5 Comparison

| metric | stage5 retained OMP | stage6 full move port | delta |
| --- | ---: | ---: | ---: |
| real | 163.13 s | 161.79 s | 1.34 s faster, 0.82% lower |
| full evolve wall | 167.3012128 s | 167.0823348 s | 0.2188780 s faster, 0.13% lower |
| move only | 104.1127577 s | 103.6495219 s | 0.4632358 s faster, 0.45% lower |
| buildCellOccupancy | 7.730683037 s | 7.483010573 s | 0.247672464 s faster, 3.20% lower |
| collision phase | 3.632583082 s | 3.613524884 s | 0.019058198 s faster, 0.52% lower |
| post fields/output | 51.2608737 s | 51.76373504 s | 0.50286134 s slower |

This is a valid but small improvement. The full move outer-control port does not explain the remaining gap to the historical `ourmeshbkp` OMP move time.

## Correctness Check

Against stage5 retained OMP:

| metric | stage5 | stage6 | difference |
| --- | ---: | ---: | ---: |
| particles | 2463658 | 2463609 | -49 (-0.001989%) |
| collisions | 35466 | 35898 | +432 (+1.218067%) |
| collision candidates | 68727 | 68794 | +67 (+0.097488%) |
| total energy | 1.242958948 | 1.242916724 | -0.000042224 (-0.003397%) |
| stuck particles | 0 | 0 | unchanged |

The result is correct for the current-tree OMP regression level. Collision count changed by about 1.22%, but candidates, particle count, and total energy remain very close, and the run completed without numerical/runtime errors.

## Historical Gap

Against historical `ourmeshbkp/omp8`:

| metric | historical bkp | stage6 full move port | gap |
| --- | ---: | ---: | ---: |
| external/real wall | 103.66 s | 161.79 s | 58.13 s slower |
| full evolve wall | 102.1225033 s | 167.0823348 s | 64.9598315 s slower |
| move only | 48.82365912 s | 103.6495219 s | 54.82586278 s slower |
| buildCellOccupancy | 6.490896398 s | 7.483010573 s | 0.992114175 s slower |
| collision phase | 9.876225069 s | 3.613524884 s | 6.262700185 s faster |
| post fields/output | 36.92670834 s | 51.76373504 s | 14.8370267 s slower |
| particles | 2463391 | 2463609 | +218 |
| collisions | 35847 | 35898 | +51 |

The main remaining gap is still move. Collision remains faster than the historical bkp run. Since the retained move outer layer only reduced formal move time by about 0.45%, the remaining move gap is likely below the `Cloud::move()` outer loop, in the OFv1706 tracking core, mesh/tet tracking path, or other version-specific particle-base behavior.

## Notes

A separate non-comparable 500-step run was also made with `openmpCollisionChunk 1`:

- log: `log.codex_omp_move_fullport_steadywall_500step_20260605`;
- `real = 159.19 s`;
- `move only = 98.8591279 s`;
- `collision phase = 3.692947083 s`.

This run is not used as the formal stage comparison because stage5 and `ourmeshbkp/omp8/system/controlDict` use `openmpCollisionChunk 8`.

## Cleanup

`ourmesh/omp8/system/controlDict` was restored after all runs.

Restored hash:

```text
ff1129081c467f7dd3a4854d94fcc26e2fa5cdf94daee11a9682b91d7a1b5b6d
```

No `dsmcFoam+` process was left running.
