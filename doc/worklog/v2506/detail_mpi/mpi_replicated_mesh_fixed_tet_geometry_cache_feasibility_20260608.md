# MPI replicated mesh fixed tet-geometry cache feasibility - 2026-06-08

## Scope

- Mode: pure MPI8 replicated mesh DLB only.  No OpenMP evidence is used here.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Purpose: read-only feasibility, memory bill, and call-chain audit for a
  fixed-mesh static tet-geometry cache.
- No source, binary, or case configuration was changed for this note.

Current strict full-fields formal baseline remains:

| metric | value |
| --- | ---: |
| external real | 133.83 s |
| move only | 72.9567669 s |
| post fields/output | 40.37082178 s |

Current case hashes:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

Post-revert symbol check:

```text
no moveTetTransition / tetTransition / cachedTetNeighbour / moveTransition / TransitionCache symbols remain
```

The current source state is not clean upstream.  It intentionally retains the
same-tet non-normalised helper in
`src/lagrangian/basic/particle/particleTemplates.C`, and that helper must not
be rolled back.

## Call Chain

The active DSMC path is:

```text
dsmcCloud::evolve_moveAndCollide()
  -> dsmcParcel::move(td, trackTime)
    -> particle::trackToFace(position() + dt*Utracking, td, true)
```

The fixed-mesh tracking hot path pulls global mesh views at the top of
`particle::trackToFace(..., true)`:

```text
particleTemplates.C:789-799
faces, points, cellCentres, faceOwner, faceNeighbour, tetBasePtIs,
cellVolumes, boundaryMesh, nInternalFaces, movingMesh, wall-distance flag
```

The retained same-tet fast path builds the current tet and calls
`dsmcTrackInsideTetNoNormalise()`:

```text
particleTemplates.C:803-831
```

The non-fast loop then repeats the static tet reconstruction on every loop
iteration:

```text
particleTemplates.C:910-1000
```

Per tet-walk loop, current code:

- updates tet ownership with `tetNeighbour(triI)` when an internal tet tri was
  crossed;
- recomputes `own`, `tetBasePtI`, `basePtI`, `facePtI`, `otherFacePtI`,
  `fPtAI`, and `fPtBI`;
- constructs a `tetPointRef` from `cellCentres[cellI]` and three face points;
- computes `tet.Sa()`, `tet.Sb()`, `tet.Sc()`, and `tet.Sd()`;
- normalises all four area vectors with four `mag()` calls;
- rebuilds the four `tetPlaneBasePtIs`;
- calls `findTris()`;
- evaluates outer `tetLambda(position_, endPosition, ...)` for candidate tris.

`findTris()` itself always evaluates four lambdas from tet centre to
`endPosition`:

```text
particleI.H:31-65
```

The outer selection evaluates lambdas from the parcel's actual current segment:

```text
particleTemplates.C:1042-1067
```

The two lambda streams are not direct duplicates.  `findTris()` uses
`Ct -> endPosition`, while the outer pass uses `position_ -> endPosition`.

The area vectors are not precomputed in OpenFOAM's `tetPointRef`.  They are
recreated by triangle normals:

```text
tetrahedronI.H:136-158
Sa = normal(b,c,d)
Sb = normal(a,d,c)
Sc = normal(a,b,d)
Sd = normal(a,c,b)
```

Internal tet transitions are still topology-heavy when they call
`crossEdgeConnectedFace()`:

```text
particleI.H:376-479  tetNeighbour(triI)
particleI.H:482-601  crossEdgeConnectedFace()
```

This is important because the rejected transition-table candidate already
removed part of this topology cost and still failed formal validation.

## Existing Measurements

The useful move subprofile from 2026-06-07 showed the remaining shape:

| metric, 200-step signal | value |
| --- | ---: |
| track calls | 533,640,621 |
| same-tet fast | 232,603,167 |
| internal tet only | 195,273,674 |
| face hits | 105,763,780 |
| tet-walk loops | 574,273,748 |
| findTris candidates | 417,477,672 |
| lambda tests | 417,477,672 |
| internal tet transitions | 273,236,294 |
| cell-face crossings | 105,763,780 |
| wall checks | 0 |

Derived:

```text
non-fast calls = 301,037,454
tet-walk loops / non-fast call = 1.907649
findTris candidates / tet-walk loop = 0.726966
track timer / move only = 53.669019%
```

Approximate lambda-style evaluations for that 200-step diagnostic:

```text
4*574,273,748 + 417,477,672 = 2,714,572,664
```

A geometry cache can reduce repeated normal construction and normalisation.  It
does not remove `findTris()` or the outer lambda dot-product work.

## Mesh Size And State Count

Current mesh header counts:

| item | count |
| --- | ---: |
| points | 210,000 |
| faces | 417,452 |
| owner entries | 417,452 |
| neighbour entries / internal faces | 207,454 |

All observed face records are quads in this case.  For quad faces, each side of
a face contributes two decomposed tets.  The state count is therefore:

```text
owner-side states    = 417,452 * 2 =   834,904
neighbour-side states = 207,454 * 2 =   414,908
total tet states                         1,249,812
```

This matches the previous transition-cache runtime report:

```text
states=1249812
transitions=3749436
table=28.60592651 MiB
offsets=3.184906006 MiB
```

The environment uses DP scalars and 32-bit labels:

```text
WM_PRECISION_OPTION=DP
WM_LABEL_SIZE=32
```

## Memory Bill

For `1,249,812` tet states:

| cache item | formula | per-rank size |
| --- | --- | ---: |
| four normalised area vectors | `states * 4 * 3 * 8` | 114.423706 MiB |
| four plane-base labels | `states * 4 * 4` | 19.070618 MiB |
| combined normals + plane bases | above two rows | 133.494324 MiB |
| one label offset per state | `(states + 1) * 4` | 4.767658 MiB |

If the cache also stores the four tet vertex point labels, add another
`19.070618 MiB` per rank.  If it stores both normalised and unnormalised
normals, add another `114.423706 MiB` per rank.

This is materially larger than the rejected transition table.  That candidate
added about `31.79 MiB` per rank and regressed the strict 500-step formal run:

| metric | current forced8 baseline | transition cache | delta |
| --- | ---: | ---: | ---: |
| real | 133.83 s | 137.97 s | +4.14 s |
| move only | 72.9567669 s | 75.98126057 s | +3.02449367 s |
| post fields/output | 40.37082178 s | 43.02448373 s | +2.65366195 s |

The geometry cache would be roughly four times larger than that rejected
candidate before storing any optional vertex labels.

## Feasibility Assessment

### Full per-state normal cache

Technically feasible:

- build a table for `(cell side, tetFaceI, tetPtI)`;
- store four normalised area vectors and four plane-base labels;
- use it only when `!mesh_.moving()` and DSMC tracking is active;
- rebuild only after mesh topology or point motion invalidates it.

Risk:

- high memory footprint: about `133.49 MiB/rank` for normals plus plane bases;
- every tet-walk loop replaces arithmetic with large-table reads;
- the case already showed that a smaller `31.79 MiB/rank` transition table
  worsened formal runtime;
- the cache cannot remove the 2.7B lambda-style evaluations seen in the
  200-step diagnostic;
- it would touch generic `particle` tracking plus `dsmcCloud` control plumbing,
  so the blast radius is larger than the likely benefit.

Verdict: do not implement as the next candidate without a narrower timing
diagnostic proving that normal construction and normalisation alone are a large
part of track time.

### Metadata-only tet-state cache

Technically feasible:

- cache `tetBasePtI`, `basePtI`, `fPtAI`, `fPtBI`, and plane-base labels;
- leave normal construction in the hot loop.

Risk:

- much smaller memory bill, but much smaller benefit;
- most expensive arithmetic remains: four triangle normals, four magnitudes,
  `findTris()`, and outer `tetLambda()`;
- unlikely to deliver the required `5-10 s` move improvement.

Verdict: not a good formal candidate unless a diagnostic shows face/base-point
index reconstruction is unexpectedly expensive.

### Lazy or LRU geometry cache

Not attractive for this path:

- track calls are extremely frequent and distributed across many tet states;
- lookup overhead and branch pressure are added to every non-fast path;
- miss handling still performs the original reconstruction;
- it is harder to reason about after DLB and migration.

Verdict: reject for this round.

## Recommended Next Step

Do not jump directly to a full static geometry cache candidate.

The useful next diagnostic is a split inside the non-fast tracking loop:

- time/count tet metadata reconstruction;
- time/count `tet.Sa/Sb/Sc/Sd` plus four `mag()` normalisations;
- time/count `findTris()`;
- time/count outer candidate `tetLambda()` selection;
- keep the diagnostic under `profileDetail true` and a new explicit switch;
- run only `10-step smoke -> 200-step signal`, then remove instrumentation.

Acceptance before any geometry-cache candidate:

1. 10-step category counts must match the existing same-tet/internal/face-hit
   shape.
2. 200-step signal must show normal construction and normalisation as a
   meaningful fraction of track time.
3. Only then implement the smallest cache variant supported by the diagnostic.
4. Any real candidate still must pass the strict full-fields 500-step gate:
   improve both `real=133.83 s` and `move only=72.9567669 s`.

## Decision

No geometry-cache source candidate was implemented on 2026-06-08.

Reason:

- the full per-state geometry cache has a large per-rank memory bill
  (`133.49 MiB` before optional vertex labels);
- the previous smaller transition table failed formal validation;
- current evidence shows many remaining lambda evaluations that a geometry
  cache cannot remove;
- without a narrower timing split, a full cache would be another high-cost
  formal candidate with weak evidence.

The next implementation work, if continuing this branch, should be
diagnostic-only timing of the non-fast tet-walk loop, not a production cache.
