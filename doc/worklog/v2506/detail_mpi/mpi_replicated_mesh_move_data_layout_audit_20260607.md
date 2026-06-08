# MPI replicated mesh move data-layout / mesh-invariant audit - 2026-06-07

## Scope

- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI8 replicated mesh DLB only, `useOpenMP false`,
  `openmpThreads 1`.
- This audit is read-only.  It does not change source, binary, or case
  configuration.

Carry-over boundaries:

- Keep the retained same-tet non-normalised helper in
  `src/lagrangian/basic/particle/particleTemplates.C`.
- Do not re-open the rejected tracking micro-patches: same-face increment,
  direct same-tet normal calculation, outside-plane tri prefilter,
  `tetNeighbour(0)` early return, same-tet sqrt-free tolerance, direct
  `Utracking`, `hasWallImpactDistance=false`, or trackingData cache.
- Do not mix OpenMP evidence into this MPI-only round.

## Current Move Path Shape

The current tree remains on the v1706 Cartesian-position tracking path:

```text
src/lagrangian/dsmc/parcels/dsmcParcel.C:158-396
src/lagrangian/basic/particle/particleTemplates.C:776-1325
```

Important mechanics:

- `dsmcParcel::move()` computes `dtCell`, `tEnd`, and `dt`, then calls
  `trackToFace(position() + dt*Utracking, td, true)`
  (`dsmcParcel.C:200-265`).
- `particle::trackToFace(..., true)` pulls static mesh data from
  `mesh_.faces()`, `mesh_.points()`, `mesh_.cellCentres()`,
  `mesh_.faceOwner()`, `mesh_.faceNeighbour()`, `mesh_.tetBasePtIs()`,
  `mesh_.cellVolumes()`, and `mesh_.boundaryMesh()` at the top of each call
  (`particleTemplates.C:789-799`).
- The retained same-tet fast path builds the current tet from face/base-point
  state and tests whether the end position stays inside it
  (`particleTemplates.C:803-831`).
- If that fails, the main tracking loop rebuilds the tet-local geometry,
  normalises four tet face-area vectors, builds plane base-point IDs, calls
  `findTris()`, then calls `tetLambda()` for candidate triangle hits
  (`particleTemplates.C:910-1067`).
- Internal tet transitions go through `tetNeighbour()` and, for edge-connected
  transitions, `crossEdgeConnectedFace()`, which scans the current cell faces
  and calls `edgeDirection()` on candidate faces
  (`particleI.H:339-525`).
- Face hits dispatch through internal-face owner/neighbour updates or through
  `boundaryMesh.whichPatch(faceI_)` and an `isA<>` patch-type chain
  (`particleTemplates.C:1151-1256`).

This is not the same structure as the historical reference.  The reference
source stores barycentric coordinates and uses displacement/fraction tracking
through `trackToAndHitFace()`; the already-written reference audit records why a
wrapper-level migration did not reproduce the reference move time:

```text
doc/worklog/v2506/detail_mpi/mpi_replicated_mesh_reference_move_diff_audit_20260607
```

## Existing Move Split

The 10-step tet-walk detail smoke remains the most useful local split:

```text
log.codex_mpi8_replicatedmesh_mpionly_tet_walk_detail_smoke_10step_20260606
```

| metric | value | share of track calls |
| --- | ---: | ---: |
| move detail track calls | 25,758,199 | 100.00% |
| same-tet no-face | 11,139,176 | 43.25% |
| internal tet only | 9,483,256 | 36.82% |
| face hits | 5,135,767 | 19.94% |
| processor hits | 0 | 0.00% |
| cyclic hits | 0 | 0.00% |
| patch hits | 16,439 | 0.06% of track calls |
| tracker max | 0 s | n/a |
| boundary max | 0.00097108 s | n/a |

The later retained forced8 move-detail smoke is consistent:

```text
log.codex_mpi8_replicatedmesh_retained_forced8_move_detail_smoke_10step_20260607
```

It reports processor hits `0`, cyclic hits `0`, tracker max `0`, patch hits
`16,453`, and boundary max `0.001501033 s`.

Interpretation:

- Boundary and tracker handling are not the move bottleneck for this case.
- Same-tet is already a large bucket and already has the retained helper.
- The next structural target is the internal-tet walk and repeated static mesh
  geometry/topology reconstruction, not another patch-dispatch or wrapper tweak.

## Static Mesh Invariants In The Hot Loop

For this fixed mesh formal case, the following are stable between mesh rebuilds:

- `faces`, `points`, `cellCentres`, `faceOwner`, `faceNeighbour`,
  `tetBasePtIs`, `cellVolumes`, `nInternalFaces`;
- per-face `basePtI`, face point count, face point successor/predecessor
  relations;
- per `(cell, tetFace, tetPt)` orientation: owner-side flag, `fPtAI/fPtBI`,
  the four tet plane base points, and the non-moving tet geometry;
- per internal face owner/neighbour target;
- per boundary face patch ID and patch kind;
- per boundary patch model ID arrays already materialised by
  `dsmcBoundaries::patchToModelIds()` and
  `dsmcBoundaries::cyclicBoundaryToModelIds()`;
- replicated `cellOwner_` and `myCells_` are stable between DLB repartitions,
  but they are not the same class of invariant as mesh topology because DLB
  changes them.

Current code already has a `dsmcLocalMesh` view with owned plus 1-ring halo
cells (`dsmcLocalMesh.C:36-101`), but the move tracking path still uses global
mesh indices and global mesh storage.  A local-indexed move path would be a
larger ownership/packing change, not a low-risk local tweak.

## What Not To Target Next

The following are poor next candidates for this case:

- Patch dispatch and boundary model calls.  Patch hits are only about
  `0.06%` of track calls and measured boundary time is around `0.001 s` in
  10-step detail runs.
- Face tracker calls.  `trackerActive_` is false unless a `dsmcFluxSurface`
  field is configured (`dsmcCloud.C:615-634`), and current move-detail logs
  show tracker max `0`.
- Direct `Utracking` or trackingData per-step cache.  Both had formal or gate
  failures and were reverted.
- Forcing `hasWallImpactDistance=false`.  The no-wall override improved a
  smoke move timer but regressed the 500-step formal gate.
- Another hand-written special case inside `tetNeighbour()` or the same-tet
  test.  Several variants passed small smoke signals but failed to survive
  formal comparison.

## Plausible Structural Direction

The next credible direction is not another arithmetic micro-patch.  It is a
diagnostic-first, table-driven treatment of static mesh invariants.

Candidate class A: exact topology-transition table.

- Build, at startup or after mesh/topology rebuild, an exact table for the
  static mesh:
  `(cellI, tetFaceI, tetPtI, triI) -> next tetFace/tetPt or face-crossing`.
- Generate it by reusing the current logic, not by writing a new heuristic.
- Rebuild or invalidate it whenever mesh topology changes.
- Use it first only under `!mesh_.moving()` and DSMC tracking.
- Validate with a 10-step detail run that the counts for same-tet, internal
  tet-only, and face hits match the existing baseline before interpreting time.

Candidate class B: static tet-geometry cache.

- Cache per-tet `basePtI`, `fPtAI`, `fPtBI`, `tetPlaneBasePtIs`, and possibly
  unnormalised/normalised tet face-area vectors for the fixed mesh.
- The retained same-tet helper already reduced one geometry path, so this
  should only be considered if a diagnostic shows that the internal-tet walk is
  still spending time in repeated geometry construction rather than in
  `tetLambda()` itself.
- Memory use can be material on this mesh, so this needs a measured table-size
  estimate before implementation.

Candidate class C: boundary face patch classification.

- Precompute `faceToPatch`, `patchKind`, and patch model IDs for boundary faces.
- This is safe conceptually, but it is low priority for the current case because
  patch hits are tiny and processor/cyclic hits are zero in the move-detail
  smoke.

## Required Diagnostic Before Implementation

Before implementing any table/cache candidate, add diagnostic-only counters
around the current path:

- count calls to `crossEdgeConnectedFace()`;
- count internal tet transitions by `triI`;
- count face sizes and cell face counts encountered by tracked parcels;
- count `boundaryMesh.whichPatch()` calls and patch kinds hit;
- count `hasWallImpactDistance && cellHasWallFaces[cellI]` guard hits;
- optionally time the internal-tet walk block separately from same-tet and
  face-hit dispatch.

Run order:

1. 10-step pure MPI8 smoke with detail counters.
2. If counts are stable and overhead is acceptable, 200-step signal.
3. Remove instrumentation, rebuild, then only implement the narrow candidate
   supported by the diagnostic.

Acceptance remains the current strict full-fields formal gate:

```text
real 133.83 s
move only 72.9567669 s
```

The candidate must improve both in a 500-step pure-MPI run.

## Conclusion

The useful move-path work now sits at the data-layout boundary:

- current high-frequency cost is inside legacy static-mesh tet tracking;
- boundary/tracker/patch model handling is too small to prioritize;
- the historical reference move path is a larger barycentric stack, not a
  one-function migration;
- a future current-tree candidate should be driven by explicit internal-tet
  transition counts and generated static mesh tables, not by more isolated
  tracking micro-edits.
