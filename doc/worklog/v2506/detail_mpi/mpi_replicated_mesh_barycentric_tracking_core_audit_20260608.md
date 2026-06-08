# MPI replicated mesh barycentric tracking-core audit - 2026-06-08

## Scope

- Mode: pure MPI8 replicated mesh DLB only.  No OpenMP candidate is considered.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Reference source:
  `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
- Purpose: decide whether the reference barycentric tracking core has a
  narrow, separable migration path for the current OFv1706 Cartesian
  tracking implementation.

Strict full-fields formal baseline for the current tree remains:

| metric | value |
| --- | ---: |
| external real | 133.83 s |
| move only | 72.9567669 s |
| post fields/output | 40.37082178 s |

Formal case hashes before the 500-step candidate run and after cleanup:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

## Source Audit

The current OFv1706 particle core stores a Cartesian position:

```text
src/lagrangian/basic/particle/particle.H:129-151
```

The DSMC move path still calls the old endpoint API:

```text
src/lagrangian/dsmc/parcels/dsmcParcel.C:235-265
dt *= trackToFace(position() + dt*Utracking, td, true);
```

The current non-same-tet path in the DSMC `trackToFace(..., true)` overload
still does:

```text
src/lagrangian/basic/particle/particleTemplates.C:981-1068
```

That block:

1. builds `tetAreas`;
2. normalises the four area vectors;
3. builds `tetPlaneBasePtIs`;
4. runs `findTris(endPosition, ...)`;
5. optionally checks wall-impact faces;
6. runs outer `tetLambda(position_, endPosition, ...)` over the candidates.

The reference tree is structurally different.  It stores particle coordinates
as barycentric state:

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/particle/particle.H:136-167
```

Its tracking core uses persistent barycentric coordinates and reverse tet
transforms:

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/particle/particle.H:180-234
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/particle/particle.C:701-828
```

The reference `trackToStationaryTri()` computes local displacement in
barycentric space, then chooses the first barycentric coordinate that reaches
zero.  This replaces the current `findTris()` screening plus outer
`tetLambda()` selection.  However, it also updates `coordinates_`, not only
`position_`.

The reference topology transitions are coupled to the barycentric state:

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/particle/particle.C:172-244
```

`changeTet()` / `changeFace()` / `changeCell()` apply `reflect()` and
`rotate()` to the barycentric coordinates.  This is not equivalent to the
current `tetNeighbour()` path, which only updates `tetFaceI_`, `tetPtI_`, and
eventually `cellI_` because the actual position remains Cartesian.

The reference also changes IO and transfer semantics:

```text
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/particle/particle.H:791-795
/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/basic/particle/particle.C:1124-1177
```

Therefore a full barycentric particle port is not a narrow move-path patch.  It
crosses particle storage, constructors, IO compatibility, processor/cyclic
transfer, patch-hit handling, and replicated flat transfer packing.

## Candidate Tested

A narrow partial candidate was still worth testing:

- keep current Cartesian `position_` storage;
- keep current patch and transfer semantics;
- for static mesh and cells that do not require wall-impact-distance handling,
  compute a transient barycentric reverse transform for the current tet;
- select the first crossed tet triangle directly from the local barycentric
  displacement;
- fall back to the old `findTris()` / `tetLambda()` path for degenerate,
  near-face, moving-mesh, and wall-impact cases.

This candidate tried to change the lower-level tet-hit selection without
repeating the earlier wrapper-level `trackToAndHitFace()` compatibility layer.

Temporary source file:

```text
src/lagrangian/basic/particle/particleTemplates.C
```

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_tet_hit_build_20260608.log
```

Build result: passed, with only the known OpenFOAM-v1706 template-instantiation
warnings.

## 10-step Smoke

Temporary controls:

```text
endTime 1.e-06;
profileDetail true;
moveDetailProfile true;
```

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_barycentric_tet_hit_smoke_20260608
```

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_tet_hit_smoke_10step_20260608
```

Result:

| metric | value |
| --- | ---: |
| real | 12.65 s |
| Total Iterations | 10 |
| move only | 1.193627235 s |
| post fields/output | 0.941428578 s |
| full evolve wall | 2.622754324 s |
| move detail track calls | 25,758,324 |
| same-tet no-face | 11,139,482 |
| internal tet only | 9,482,928 |
| face hits | 5,135,914 |
| track max | 0.755557751 s |
| particles | 2,065,850 |
| stuck particles | 0 |
| collisions | 2,891 |
| collision candidates | 2,988 |
| total energy | 1.150620381 |
| DLB rebalances | 0 |
| OpenMP enabled | 0 |
| OpenMP max threads | 1 |

No `Fatal`, `Segmentation`, `Floating`, `NaN`, or `nan` entries were found.

## 200-step Signal

Temporary controls:

```text
endTime 2.e-05;
profileDetail true;
moveDetailProfile true;
```

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_barycentric_tet_hit_signal_20260608
```

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_tet_hit_signal_200step_20260608
```

Signal comparison:

| metric | previous sampled signal | transient barycentric candidate |
| --- | ---: | ---: |
| real | 58.53 s | 56.82 s |
| move only | 30.93683107 s | 29.18178768 s |
| DLB rebalances | 2 | 2 |

Candidate 200-step details:

| metric | value |
| --- | ---: |
| post fields/output | 16.86202402 s |
| full evolve wall | 52.10364208 s |
| move detail track calls | 533,647,392 |
| same-tet no-face | 232,587,947 |
| internal tet only | 195,288,939 |
| face hits | 105,770,506 |
| track max | 14.92046008 s |
| particles | 2,220,244 |
| stuck particles | 0 |
| collisions | 15,240 |
| collision candidates | 25,108 |
| total energy | 1.169676904 |
| OpenMP enabled | 0 |
| OpenMP max threads | 1 |

The 200-step signal was positive, so the candidate was advanced to a strict
500-step full-fields formal run.

## 500-step Formal

Formal controls were restored before the run:

```text
endTime 5.e-05;
profileDetail false;
```

Control and field hashes:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_barycentric_tet_hit_formal_20260608
```

Log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_barycentric_tet_hit_formal_500step_20260608
```

Formal comparison against the current strict full-fields baseline:

| metric | strict baseline | transient barycentric candidate | delta |
| --- | ---: | ---: | ---: |
| external real | 133.83 s | 138.37 s | +4.54 s |
| move only | 72.9567669 s | 76.32640482 s | +3.36963792 s |
| post fields/output | 40.37082178 s | 41.49157059 s | +1.12074881 s |

Candidate formal details:

| metric | value |
| --- | ---: |
| move+collide wall | 101.5629632 s |
| buildCellOccupancy | 7.199321263 s |
| collision phase | 15.58302029 s |
| full evolve wall | 137.4397216 s |
| particles | 2,463,669 |
| stuck particles | 0 |
| collisions | 34,213 |
| collision candidates | 61,017 |
| total energy | 1.238646238 |
| DLB checks | 500 |
| DLB rebalances | 8 |
| OpenMP enabled | 0 |
| OpenMP max threads | 1 |

No hard runtime errors were found, but the formal performance gate failed:
both external wall time and move-only time were worse than the strict
full-fields baseline.

## Decision

Rejected and reverted.

The candidate showed a useful 200-step diagnostic signal, but it did not
survive the 500-step full-fields formal test.  The formal regression means this
transient local barycentric hit-selection path should not be retained.

Interpretation:

- A wrapper-level `trackToAndHitFace()` compatibility layer was already rejected
  earlier.
- A transient barycentric selection helper avoids the full storage/IO port, but
  it still adds per-tet reverse-transform work and changes tolerance/selection
  details enough that formal performance does not improve.
- The reference barycentric core gets its value as a coherent particle storage
  and topology-update model, not as a small helper injected into the current
  Cartesian `trackToFace()` loop.
- A full barycentric particle port remains too broad for the current narrow
  pure-MPI move optimisation round because it crosses particle IO, processor
  transfer, cyclic handling, replicated flat transfer, and boundary semantics.

No further narrow barycentric/tracking-core candidate is recommended in this
round.

## Cleanup State

Temporary candidate code was removed from:

```text
src/lagrangian/basic/particle/particleTemplates.C
```

Cleanup rebuild log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_barycentric_tet_hit_revert_build_20260608.log
```

Cleanup rebuild result: passed, with only the known OpenFOAM-v1706
template-instantiation warnings.

Post-cleanup checks:

- no `dsmcTrackStaticTetFirstCross`, `selectedTetTri`, or `reverseTransform`
  candidate symbols remain in `src/lagrangian/basic/particle` or
  `src/lagrangian/dsmc/parcels`;
- `git diff -- src/lagrangian/basic/particle/particleTemplates.C
  src/lagrangian/basic/particle/particle.H
  src/lagrangian/basic/particle/particleI.H
  src/lagrangian/dsmc/parcels/dsmcParcel.C
  src/lagrangian/dsmc/parcels/dsmcParcel.H` is empty;
- formal case hashes are restored to:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

The retained same-tet non-normalised helper in `particleTemplates.C` was not
rolled back.
