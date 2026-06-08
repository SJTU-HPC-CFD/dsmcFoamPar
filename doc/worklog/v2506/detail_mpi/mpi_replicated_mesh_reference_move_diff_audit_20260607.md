# MPI replicated-mesh reference move-path diff audit - 2026-06-07

## Scope

- Current tree: `/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb`
- Reference source: `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx`
- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI8 replicated mesh only; OpenMP is out of scope for this audit.

This audit is read-only with respect to the reference source.

## Baseline labels

Keep these baselines separate:

| label | source/control meaning | real/external wall | move only |
| --- | --- | ---: | ---: |
| current strict full-fields baseline | current retained source, formal full-fields case | 133.83 s | 72.9567669 s |
| historical reference / older contract | preserved `ourmeshbkp` pure MPI8 reference log | 117.31 s | 47.1652636 s |
| no-fields optional mode | output fields disabled; not full-fields equivalent | 119.69 s | not a strict full-fields baseline |

The historical reference is useful for locating the large move gap, but its
profile-output code path is older.  Current candidate decisions still use the
current strict full-fields baseline as the acceptance gate.

## Current move path

Current `dsmcParcel::move()` is an OF-v1706 Cartesian-position tracking loop.
Evidence:

- `src/lagrangian/dsmc/parcels/dsmcParcel.C:200-203` computes `dtCell` and
  `tEnd`;
- `src/lagrangian/dsmc/parcels/dsmcParcel.C:212-229` loops while `tEnd` remains;
- `src/lagrangian/dsmc/parcels/dsmcParcel.C:241` and `:264` call
  `trackToFace(position() + dt*Utracking, td, true)`;
- `src/lagrangian/dsmc/parcels/dsmcParcel.C:268-275` manually updates
  `tEnd` and `stepFraction()`;
- `src/lagrangian/basic/particle/particleTemplates.C:263-357` begins the
  legacy `trackToFace(endPosition, td)` implementation;
- `src/lagrangian/basic/particle/particleTemplates.C:1000-1149` shows the DSMC
  tracking loop repeatedly building tet-local geometry, finding candidate tris,
  computing `tetLambda`, and changing `faceI_` / tet indices;
- `src/lagrangian/basic/particle/particleTemplates.C:1151-1256` then handles
  face crossing and patch dispatch.

The retained same-tet non-normalised helper lives inside this current
`trackToFace(..., true)` path and must remain retained; this audit does not
roll it back.

## Reference move path

The reference source uses a different particle representation and a different
lower-level tracking API.  Evidence:

- `/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc/parcels/dsmcParcel.C:156-163`
  loops on `stepFraction() < 1`, computes `deviationFromMeshCentre()`, and calls
  `trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td)`;
- the profiled branch at reference `dsmcParcel.C:237-250` uses the same
  displacement/fraction form and records track-wall time around
  `trackToAndHitFace`;
- reference `particle.H:137-155` stores barycentric particle state
  (`coordinates_`, `tetFacei_`, `tetPti_`, `stepFraction_`);
- reference `particle.H:651-662` declares `trackToAndHitFace()` and
  `deviationFromMeshCentre()`;
- reference `particle.C:632-692` contains the barycentric `track()` loop with
  fallback tracking and step-fraction update;
- reference `particle.C:1101` defines `deviationFromMeshCentre()`;
- reference `particleIO.C` reads/writes barycentric coordinates and tet indices,
  not only Cartesian `position`.

The reference path therefore is not a wrapper around the current
`trackToFace(position()+dt*Utracking, ...)` call shape.  It depends on the
particle storage model, IO, transfer semantics, tet state, and barycentric
transform helpers.

## Structural differences

| area | current tree | reference source | implication |
| --- | --- | --- | --- |
| particle position state | Cartesian `position_` plus legacy tet indices | barycentric `coordinates_` with `tetFacei_`, `tetPti_`, `stepFraction_` | not a local `dsmcParcel.C` change |
| move loop progress | explicit `tEnd/dtCell` countdown | `stepFraction() < 1` with displacement/fraction tracking | different fraction accounting |
| track target | absolute `position() + dt*Utracking` | relative displacement adjusted by `deviationFromMeshCentre()` | different lower-level geometry math |
| tracking API | `trackToFace(endPosition, td, true)` | `trackToAndHitFace(displacement, fraction, cloud, td)` | wrapper alone does not reproduce reference |
| IO/transfer surface | current v1706 particle IO and replicated flat transfer | barycentric state appears in IO and transfer paths | full port crosses correctness boundaries |
| patch dispatch | legacy `trackToFace` handles boundary patch classes after face hit | reference `hitFace` / `hitBoundaryFace` are integrated with barycentric tracking | patch behavior must be audited together |

## Evidence from rejected minimum migration

The compatibility `trackToAndHitFace` wrapper attempt compiled and passed smoke,
but failed the formal 500-step run:

| metric | pre-candidate current | compatibility API candidate |
| --- | ---: | ---: |
| real | 135.30 s | 140.79 s |
| move only | 72.59819434 s | 82.92498684 s |
| buildCellOccupancy | 7.432408543 s | 8.25212966 s |
| collision phase | 12.97955973 s | 15.08526616 s |
| post fields/output | 38.10104302 s | 37.82146892 s |

This result supports the structural reading: changing only the outer call shape
does not recover the reference move performance.

## Conclusion

The move gap to the historical reference is not explained by one small
`dsmcParcel::move()` detail, a `trackingData` cache, or a wrapper named
`trackToAndHitFace`.  The reference uses a broader barycentric particle/tracking
stack.  Closing the reference move gap would require one of two larger paths:

1. port the barycentric particle and tracking stack with IO, transfer, cyclic,
   and replicated-mesh packing semantics audited together; or
2. derive a current-tree equivalent fast path inside legacy
   `particle::trackToFace(..., true)` based on measured invariants, without
   changing the particle storage model.

The second path is lower risk for this v1706 branch, but it should be driven by
data-layout and mesh-invariant analysis rather than additional one-off tracking
micro-patches from the rejected list.
