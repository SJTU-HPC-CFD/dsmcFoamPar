# MPI replicated mesh move transition-cache feasibility - 2026-06-07

## Scope

- Mode: pure MPI8 replicated mesh DLB, no OpenMP evidence.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Current strict full-fields formal baseline:
  - `real = 133.83 s`
  - `move only = 72.9567669 s`
  - `post fields/output = 40.37082178 s`
- Historical reference / older contract remains separate:
  - `ourmeshbkp` external wall `117.31 s`
  - `move only = 47.1652636 s`
- `no-fields` remains diagnostic/optional only and is not mixed into this
  full-fields gate.

## Baseline Runtime State Check

Current formal control hashes:

```text
system/controlDict        8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
system/fieldPropertiesDict 73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9
```

The current control keeps:

```text
profileSummary true;
profileDetail false;
replicatedMesh true;
replicatedMeshFlatTransfer true;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

There is no active `moveDetailProfile` setting in the formal control.  The
retained same-tet non-normalised helper in `particleTemplates.C` remains part
of the current source state and must not be removed by this candidate.

## Mesh Scale

Read-only counts from the replicated mesh case:

| item | count |
| --- | ---: |
| points | 210000 |
| faces | 417452 |
| internal faces | 207454 |
| boundary faces | 209998 |
| face arity | all `4(...)` |
| owner-side tets | 834904 |
| internal neighbour-side tets | 414908 |
| total cell-tet states | 1249812 |

For this all-quad mesh, each face side contributes `4 - 2 = 2` tets.  The
state count is therefore:

```text
2 * faces + 2 * internalFaces = 2*417452 + 2*207454 = 1249812
```

## Memory Budget

Narrow topology-transition table:

```text
tri1/tri2/tri3 transitions per state = 3
payload per transition = next tetFace + next tetPt = 2 labels
label size in linux64IccDPInt32Opt = 4 bytes
1249812 * 3 * 2 * 4 = 29995488 bytes = 28.61 MiB
```

Offset tables for mapping `(cell, tetFace, tetPt)` to state:

```text
owner-side start per face        = 417452 labels
neighbour-side start per face    = 417452 labels
two labelLists total             = 834904 * 4 = 3.18 MiB
```

Approximate narrow table total:

```text
28.61 MiB + 3.18 MiB = 31.79 MiB
```

Static geometry cache is materially heavier.  Caching four `vector` face
normals per state would cost:

```text
1249812 * 4 * 3 scalars * 8 bytes = 119982336 bytes = 114.42 MiB
```

That excludes plane-base IDs and other metadata, so the first implementation
should not cache normals.

## Candidate Design

Implement only a topology table for non-moving DSMC tracking:

```text
(cellI, tetFaceI, tetPtI, triI=1/2/3) -> (next tetFaceI, next tetPtI)
```

Rules:

- Generate the table from the current `tetNeighbour()` logic, including the
  `crossEdgeConnectedFace()` search and duplicate-face guard.
- Do not cache or alter tri0.  tri0 remains the existing face-hit path and the
  later owner/neighbour or boundary dispatch keeps changing `cellI_`.
- Do not change `findTris()`, `tetLambda()`, same-tet geometry, wall handling,
  tracker handling, boundary handling, or replicated-mesh migration.
- Gate the runtime path with `moveTetTransitionCache true`.
- Set the cache pointers only for `!mesh_.moving()`; moving mesh falls back to
  the current `tetNeighbour()` path.

## Validation Gate

Run order:

1. 10-step smoke: build/run only, check no addressing failure and no correctness
   anomaly.
2. 200-step signal: must include forced DLB steps 120 and 170.
3. 500-step formal: compare against the strict full-fields baseline above.

Accept only if the 500-step full-fields run improves both:

- end-to-end `real`;
- `move only`.

If the formal run regresses, reject and revert the source candidate.  A positive
10-step or 200-step signal is not enough to retain this code.
