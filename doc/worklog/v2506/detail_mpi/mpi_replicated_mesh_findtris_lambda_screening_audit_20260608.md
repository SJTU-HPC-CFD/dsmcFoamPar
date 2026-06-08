# MPI replicated mesh findTris / lambda screening audit - 2026-06-08

## Scope

- Mode: pure MPI8 replicated mesh DLB only.  No OpenMP evidence is used here.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Purpose: decide whether a narrow current-tree `findTris()` / `tetLambda()`
  source candidate remains after the 2026-06-08 tet-walk sampled profile.
- This is a read-only audit.  It does not change source, binary, or case
  configuration.

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

## Current Code Path

In the DSMC `particle::trackToFace(..., true)` overload, the non-fast path
does this on every tet-walk loop:

```text
particleTemplates.C:910-1070
```

The key sequence is:

1. apply `tetNeighbour(triI)` after an internal tet crossing;
2. rebuild face/base-point metadata and construct `tetPointRef`;
3. build and normalise four tet area vectors;
4. set four plane-base point labels;
5. call `findTris(endPosition, tris, tet, tetAreas, tetPlaneBasePtIs, tol)`;
6. run the outer `tetLambda(position_, endPosition, ...)` selection over the
   `findTris()` candidates.

`findTris()` itself is:

```text
particleI.H:31-65
```

It clears a four-entry `DynamicList`, computes `Ct = tet.centre()`, then calls
`tetLambda(Ct, endPosition, ...)` four times.  The outer selection later calls
`tetLambda(position_, endPosition, ...)` only for candidates.

These lambda streams are not duplicates:

```text
findTris line: Ct -> endPosition
outer line:    position_ -> endPosition
```

Therefore the `findTris()` lambda values cannot simply be reused for the outer
nearest-hit selection.

## Latest Sampled Evidence

The 2026-06-08 sampled profile reported the following 200-step split:

| section | share of sampled split |
| --- | ---: |
| `findTris()` | 39.618008% |
| `tetNeighbour()` / edge-connected topology | 24.881215% |
| area normals + plane bases | 16.289417% |
| metadata reconstruction | 9.471115% |
| outer lambda selection | 9.740245% |

The same 200-step run reported:

| metric | value |
| --- | ---: |
| track calls | 533,643,543 |
| same-tet no-face | 232,602,591 |
| internal tet only | 195,276,408 |
| face hits | 105,764,544 |
| tet-walk loops | 574,273,170 |
| tet-walk samples | 560,016 |

This means a narrow candidate must materially change `findTris()` itself or the
topology transition path to have enough upside.  Optimising only plane-base
labels, normal construction, or outer candidate selection is too small.

## Already Rejected Or Excluded Routes

### Fixed tris array / inline findTris

Already tested and rejected in:

```text
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_tracking_tet_walk_followup_20260606.md
```

Candidate:

- replace `DynamicList<label> tris(4)` with `FixedList<label, 4>` plus `nTris`;
- inline the old `findTris()` loop for this DSMC path.

Result:

| metric | tet-walk diagnostic | fixed tris array |
| --- | ---: | ---: |
| move only | 1.268614198 s | 1.329774664 s |
| track max | 0.783745186 s | 0.818885698 s |

Decision was reject and revert.  Do not repeat this as a "container cleanup"
candidate.

### Static tetLambda helper

Already tested and formally rejected.

Candidate:

- keep the existing geometry and tolerance logic;
- for static mesh, call a local helper that avoids the `mesh_.moving()` branch
  inside `tetLambda()`.

The 10-step smoke was positive, but the 500-step formal failed:

| metric | forced8 baseline | same-tet retained range | static lambda |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 127.91-129.93 s | 135.56 s |
| move only | 69.49945105 s | 67.56621687-69.11276280 s | 71.57134787 s |
| full evolve wall | 123.0267073 s | 118.7495988-120.5815986 s | 130.1452013 s |

Decision was reject and revert.  Do not repeat a local static-lambda helper.

### Outside-plane tri prefilter

Already tested and formally rejected.

Candidate:

- replace `findTris()` centre-ray screening with direct outside-plane tests;
- leave outer `tetLambda(position_, endPosition, ...)` selection unchanged.

It had a strong 10-step move signal, but formal failed:

| metric | forced8 baseline | same-tet retained range | outside-plane |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 127.91-129.93 s | 143.30 s |
| move only | 69.49945105 s | 67.56621687-69.11276280 s | 74.50898682 s |
| full evolve wall | 123.0267073 s | 118.7495988-120.5815986 s | 137.3230093 s |

Decision was reject and revert.  Do not re-open outside-plane screening.

### Direct lambda reuse

Invalid for this code path.

`findTris()` uses `Ct -> endPosition`; the outer selection uses
`position_ -> endPosition`.  Reusing the centre-ray lambda for the parcel
segment changes the tracking decision and is not a valid optimisation.

### Metadata-only or geometry-only cache

The sampled profile shows metadata is under `10%` and area normals plus plane
bases are about `16-17%` of the sampled tet-walk split.  The full geometry
cache feasibility note also estimated about `133.49 MiB/rank` for normals plus
plane-base labels.  This is too large for the likely gain and smaller than the
dominant `findTris()` and topology-transition costs.

## Remaining Narrow Candidate Space

No low-risk narrow current-tree candidate remains in `findTris()` /
`tetLambda()` that is not already rejected, invalid, or too small:

- replacing `DynamicList` was tried and worsened the tracking timer;
- bypassing the static `tetLambda()` branch was tried and failed formal;
- replacing centre-ray screening with outside-plane tests was tried and failed
  formal;
- direct lambda reuse is semantically wrong;
- geometry-only caching targets a smaller bucket and has a high memory bill;
- candidate-order or early-exit changes would alter tie/tolerance behaviour and
  are high-risk without a new tracking formulation.

The sampled split does not justify another current-tree micro-cache or
arithmetic patch.

## Decision

Do not implement a new `findTris()` / `tetLambda()` narrow source candidate in
this round.

The next meaningful direction is a larger tracking-core audit:

1. Compare the reference barycentric/displacement tracking path against the
   current OFv1706 Cartesian `trackToFace()` path at the lower-level algorithm
   and data-layout level.
2. Identify whether a partial port can change candidate screening and topology
   transition together, not only wrap the old `trackToFace()` call shape.
3. Keep the strict acceptance gate: a retained candidate must improve both the
   current full-fields 500-step `real=133.83 s` and
   `move only=72.9567669 s`.

The earlier compatibility wrapper for `trackToAndHitFace()` already failed; the
remaining barycentric/tracking API work must be a deeper migration audit, not a
wrapper-level retry.
