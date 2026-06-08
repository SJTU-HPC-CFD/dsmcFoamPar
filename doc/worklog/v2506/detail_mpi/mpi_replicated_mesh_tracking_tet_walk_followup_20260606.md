# MPI replicated mesh tracking tet-walk follow-up - 2026-06-06

## Scope

- Case: `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Mode: pure MPI, `useOpenMP false`, `openmpThreads 1`
- Retained DLB configuration:
  `replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);`
- Current retained source candidate:
  `src/lagrangian/basic/particle/particleTemplates.C`
  same-tet non-normalised inside check.

This note records the finer move-detail split and the rejected
tracking fast-path attempts after the retained same-tet candidate.

## Refined move-detail split

Diagnostic controls:

```text
endTime 1.e-06;
profileDetail true;
moveDetailProfile true;
```

Diagnostic log:

```text
log.codex_mpi8_replicatedmesh_mpionly_tet_walk_detail_smoke_10step_20260606
```

Result:

| metric | value |
| --- | ---: |
| Total Iterations | 10 |
| move only | 1.268614198 s |
| move detail track calls | 25,758,199 |
| same-tet no-face | 11,139,176 |
| internal tet only | 9,483,256 |
| face hits | 5,135,767 |
| track max | 0.783745186 s |
| boundary max | 0.00097108 s |
| OpenMP max threads | 1 |

Interpretation:

- The same-tet path is a large piece of the tracking workload, and the retained
  non-normalised inside check already targets that path.
- The next large bucket is not boundary handling.  About `9.48M` of `25.76M`
  10-step `trackToFace(..., true)` calls stay within the same cell but walk to a
  different internal tet.
- Further move candidates should focus on the old tracking core around
  `tetNeighbour()` and internal-tet transitions, not DLB cadence, boundary
  callbacks, reduced-D constraints, or wrapper-style API migration.

## Rejected candidate: quad edge fast path

Candidate:

- File: `src/lagrangian/basic/particle/particleI.H`
- Function: `particle::crossEdgeConnectedFace()`
- Idea: for the common hexahedral/quad-face case, find the neighbouring face by
  direct four-point edge checks before falling back to the generic
  `face::edgeDirection()` logic.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_quad_edge_fastpath_build_20260606.log
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_mpionly_quad_edge_fastpath_smoke_10step_20260606
```

Comparison:

| metric | tet-walk diagnostic | quad edge fast path |
| --- | ---: | ---: |
| move only | 1.268614198 s | 1.327664556 s |
| move only cpu | 1.18 s | 1.26 s |
| track calls | 25,758,199 | 25,758,347 |
| same-tet no-face | 11,139,176 | 11,139,132 |
| internal tet only | 9,483,256 | 9,483,333 |
| face hits | 5,135,767 | 5,135,882 |
| track max | 0.783745186 s | 0.857916341 s |
| boundary max | 0.00097108 s | 0.001503829 s |

Decision: reject and do not run a 500-step validation.

The candidate made the targeted tracking timer worse in the 10-step diagnostic.
The extra specialised face scan did not beat the generic path on this mesh, so
the change was removed from `particleI.H`.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_quad_edge_fastpath_revert_build_20260606.log
```

Revert build result: passed, with only the existing OpenFOAM template
instantiation warnings.  After the revert, `particleI.H` has no remaining diff
from this rejected candidate.

## Rejected candidate: fixed tris array

Candidate:

- File: `src/lagrangian/basic/particle/particleTemplates.C`
- Function: DSMC `particle::trackToFace(..., true)` overload
- Idea: replace the per-call `DynamicList<label> tris(4)` with a fixed
  four-slot `FixedList<label, 4>` plus `nTris`, and inline the old
  `findTris()` loop for this DSMC path only.

Reasoning:

- The diagnostic run executes about `25.76M` tracking calls over 10 steps.
- The candidate tried to remove a small heap-backed dynamic container from the
  repeated internal-tet path without changing tracking decisions.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_fixed_tris_build_20260606.log
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_mpionly_fixed_tris_smoke_10step_20260606
```

Comparison:

| metric | tet-walk diagnostic | fixed tris array |
| --- | ---: | ---: |
| real | 11.63 s | 11.36 s |
| move only | 1.268614198 s | 1.329774664 s |
| move only cpu | 1.18 s | 1.24 s |
| track calls | 25,758,199 | 25,758,383 |
| same-tet no-face | 11,139,176 | 11,139,678 |
| internal tet only | 9,483,256 | 9,482,760 |
| face hits | 5,135,767 | 5,135,945 |
| track max | 0.783745186 s | 0.818885698 s |
| boundary max | 0.00097108 s | 0.001518258 s |

Decision: reject and do not run a 500-step validation.

Although this removed the dynamic list construction, the local inlined loop made
the tracking timer worse in the 10-step diagnostic.  The change was reverted.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_fixed_tris_revert_build_20260606.log
```

Revert build result: passed, with only the existing OpenFOAM template
instantiation warnings.  After the revert, `particleTemplates.C` again only has
the retained same-tet non-normalised inside-check diff.

## Temporary diagnostic: tri split

Purpose:

- Split the internal tracking loop by selected tet triangle to locate which
  `tetNeighbour()` branches are most active.
- This diagnostic was default-off but inserted checks inside the tracking loop,
  so it was used only for one 10-step run and then removed before formal
  candidate validation.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_move_detail_tri_split_build_20260606.log
```

Smoke log:

```text
log.codex_mpi8_replicatedmesh_mpionly_tri_split_smoke_10step_20260606
```

Result:

| metric | value |
| --- | ---: |
| track calls | 25,758,346 |
| same-tet no-face | 11,139,638 |
| internal tet only | 9,482,842 |
| face hits | 5,135,866 |
| tet loop iterations | 27,916,564 |
| tri0 hits | 5,135,866 |
| tri1 hits | 4,844,653 |
| tri2 hits | 5,689,524 |
| tri3 hits | 2,763,679 |

Interpretation:

- `tri0` equals the real face-hit count.
- Internal tet transitions are split across `tri1/2/3`, with `tri2` the largest
  bucket and `tri1` close behind.
- This justified one narrow `tetNeighbour()` same-face indexing candidate, but
  the diagnostic itself was removed after collecting the data.

## Rejected candidate: same-face tet index increment

Candidate:

- File: `src/lagrangian/basic/particle/particleI.H`
- Function: `particle::tetNeighbour(label triI)`
- Idea: in guarded same-face branches, replace `f.fcIndex(tetPtI_)` with
  `++tetPtI_` and `f.rcIndex(tetPtI_)` with `--tetPtI_`.

Reasoning:

- The existing guards prevent wraparound in those branches, so the replacement
  is topologically equivalent.
- The candidate targets frequent `tri2/tri3` same-face internal tet moves while
  leaving `crossEdgeConnectedFace()` untouched.

Build logs:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_same_face_inc_build_20260606.log
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_same_face_inc_notri_build_20260606.log
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_same_face_inc_formal_build_20260606.log
```

Initial detail smoke, compared with the temporary tri-split diagnostic:

| metric | tri-split baseline | same-face increment |
| --- | ---: | ---: |
| move only | 1.475285141 s | 1.417757571 s |
| move only cpu | 1.35 s | 1.29 s |
| track max | 0.88837908 s | 0.899226443 s |

Because `track max` did not improve and the tri-split diagnostic added
hot-loop overhead, a no-detail A/B was run after removing the tri-split
diagnostic.

No-detail 10-step A/B:

| metric | same-tet baseline | same-face increment |
| --- | ---: | ---: |
| real | 10.86 s | 8.63 s |
| move+collide wall | 1.775343523 s | 1.697345013 s |
| move only | 1.223100367 s | 1.21153991 s |
| full evolve wall | 2.433862566 s | 2.181271882 s |

The no-detail smoke was mildly positive, so the candidate was allowed one
500-step formal run.

Formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_same_face_inc_500step_20260606
```

Formal comparison:

| metric | forced8 baseline | same-tet retained range | same-face increment |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 127.91-129.93 s | 124.57 s |
| move only | 69.49945105 s | 67.56621687-69.11276280 s | 72.55420743 s |
| move+collide wall | 90.87479808 s | 88.73592885-90.75870920 s | 94.51620716 s |
| full evolve wall | 123.0267073 s | 118.7495988-120.5815986 s | 122.918634 s |
| particles max/min | 2.95295203 | 1.585453037-1.854821906 | 2.289947846 |
| DLB rebalances | 8 | 8 | 8 |

Decision: reject and revert.

Although external `real` improved in this single run, the required move metric
regressed badly against both the forced8 baseline and the retained same-tet
candidate.  This candidate does not meet the move-performance objective.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_same_face_inc_revert_build_20260606.log
```

Revert build result: passed, with only the existing OpenFOAM template
instantiation warnings.  After the revert, `particleI.H` has no remaining diff.

## Rejected candidate: static tetLambda helper

Candidate:

- File: `src/lagrangian/basic/particle/particleTemplates.C`
- Function: DSMC `particle::trackToFace(..., true)` overload
- Idea: keep the normalised tet geometry and old tolerance logic, but for the
  static-mesh path call a local `dsmcTrackTetLambdaStatic()` helper instead of
  the generic `tetLambda()` function, avoiding the repeated `mesh_.moving()`
  branch and generic moving-mesh parameter path inside the hot loop.
- Moving-mesh tracking still falls back to the original `findTris()` and
  `tetLambda()` paths.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_static_lambda_build_20260606.log
```

No-detail 10-step smoke:

```text
log.codex_mpi8_replicatedmesh_mpionly_static_lambda_nodetail_smoke_10step_20260606
```

Comparison:

| metric | same-tet baseline | static lambda |
| --- | ---: | ---: |
| real | 10.86 s | 10.71 s |
| move+collide wall | 1.775343523 s | 1.647035455 s |
| move only | 1.223100367 s | 1.121454305 s |
| full evolve wall | 2.433862566 s | 2.336798233 s |

10-step smoke decision: promising enough for formal validation.

The required 500-step formal run was attempted with:

```text
bash -lc "source doc/scripts/env.sh; cd run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh; /usr/bin/time -p mpirun -np 8 dsmcFoam+ > log.codex_mpi8_replicatedmesh_mpionly_static_lambda_500step_20260606 2>&1"
```

The command initially did not start because the escalation approval-review layer
returned a 503 service error.  After retrying, the run completed on 2026-06-07.

Formal 500-step log:

```text
log.codex_mpi8_replicatedmesh_mpionly_static_lambda_500step_20260606
```

Log note: this file name was later accidentally reused after the source had
already been rebuilt back to the same-tet state.  That later run was copied to
`log.codex_mpi8_replicatedmesh_mpionly_sametet_current_rerun_500step_20260607`.
The static-lambda metrics below are retained from the earlier extraction made
before that overwrite.

Formal comparison:

| metric | forced8 baseline | same-tet retained range | static lambda |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 127.91-129.93 s | 135.56 s |
| move only | 69.49945105 s | 67.56621687-69.11276280 s | 71.57134787 s |
| move+collide wall | 90.87479808 s | 88.73592885-90.75870920 s | 95.16644694 s |
| full evolve wall | 123.0267073 s | 118.7495988-120.5815986 s | 130.1452013 s |
| particles max/min | 2.95295203 | 1.585453037-1.854821906 | 1.430501722 |
| DLB rebalances | 8 | 8 | 8 |

Decision: reject and revert.

The 10-step signal did not hold at 500 steps.  `move only` regressed against
both the forced8 baseline and the retained same-tet candidate, so the source
change was removed.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_revert_static_lambda_build_20260607.log
```

Revert build result: passed, with only the existing OpenFOAM template
instantiation warnings.  After the revert, `particleTemplates.C` retains only
the same-tet non-normalised inside check from this candidate sequence.

## Rejected candidate: direct same-tet normal calculation

Candidate:

- File: `src/lagrangian/basic/particle/particleTemplates.C`
- Function: DSMC `particle::trackToFace(..., true)` same-tet fast path
- Idea: keep the retained non-normalised half-space test, but avoid constructing
  `tetPointRef` / `triangle` temporaries for the same-tet check by computing the
  four face normals directly from the tet points.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_direct_same_tet_build_20260607.log
```

No-detail 10-step smoke:

```text
log.codex_mpi8_replicatedmesh_mpionly_direct_same_tet_smoke_10step_20260607
```

Comparison:

| metric | same-tet no-detail baseline | direct same-tet |
| --- | ---: | ---: |
| real | 10.86 s | 12.25 s |
| move+collide wall | 1.775343523 s | 1.886675629 s |
| move only | 1.223100367 s | 1.282788205 s |
| full evolve wall | 2.433862566 s | 2.750851197 s |

Decision: reject and revert.

The direct-normal form degraded the smoke move metric, so it was not promoted to
500-step validation.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_direct_same_tet_revert_build_20260607.log
```

Revert build result: passed, with only the existing OpenFOAM template
instantiation warnings.  The formal case controls were restored to 500-step
pure MPI after the smoke.

## Rejected candidate: outside-plane tri prefilter

Candidate:

- File: `src/lagrangian/basic/particle/particleTemplates.C`
- Function: DSMC `particle::trackToFace(..., true)` internal tet-walk loop
- Idea: on static meshes, replace the `findTris()` ray test from tet centre to
  `endPosition` with a direct outside-plane test for the four tet planes.  The
  later `tetLambda(position_, endPosition, ...)` selection of the nearest
  crossing was left unchanged, and moving meshes still used the original
  `findTris()` path.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_outside_plane_build_20260607.log
```

No-detail 10-step smoke:

```text
log.codex_mpi8_replicatedmesh_mpionly_outside_plane_smoke_10step_20260607
```

Smoke comparison:

| metric | same-tet no-detail baseline | outside-plane |
| --- | ---: | ---: |
| real | 10.86 s | 10.81 s |
| move+collide wall | 1.775343523 s | 1.442562229 s |
| move only | 1.223100367 s | 0.851679179 s |
| full evolve wall | 2.433862566 s | 2.141245599 s |

The smoke result was strongly positive, so the candidate was promoted to one
500-step formal run.

Formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_outside_plane_500step_20260607
```

Formal comparison:

| metric | forced8 baseline | same-tet retained range | outside-plane |
| --- | ---: | ---: | ---: |
| real | 132.29 s | 127.91-129.93 s | 143.30 s |
| move only | 69.49945105 s | 67.56621687-69.11276280 s | 74.50898682 s |
| move+collide wall | 90.87479808 s | 88.73592885-90.75870920 s | 100.9467873 s |
| full evolve wall | 123.0267073 s | 118.7495988-120.5815986 s | 137.3230093 s |
| particles max/min | 2.95295203 | 1.585453037-1.854821906 | 1.911010006 |
| DLB rebalances | 8 | 8 | 8 |

Because same-day runs were slower than the historical retained range, the
outside-plane candidate was also compared against a clean same-tet confirm run
after reverting the candidate:

```text
log.codex_mpi8_replicatedmesh_mpionly_sametet_revert_confirm_500step_20260607
```

| metric | same-tet confirm | outside-plane |
| --- | ---: | ---: |
| real | 139.66 s | 143.30 s |
| move only | 75.14025444 s | 74.50898682 s |
| move+collide wall | 99.01058420 s | 100.9467873 s |
| full evolve wall | 134.0601301 s | 137.3230093 s |
| particles max/min | 2.265839756 | 1.911010006 |
| DLB rebalances | 8 | 8 |

Decision: reject and revert.

Despite the large 10-step move improvement, the 500-step formal run regressed
against both the forced8 baseline and the retained same-tet range.  Against the
near-time same-tet confirm it only shaved about `0.63 s` from `move only`, while
making `real`, `move+collide wall`, and `full evolve wall` worse.  The source
change was removed.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_outside_plane_revert_build_20260607.log
```

Revert build result: passed, with only the existing OpenFOAM template
instantiation warnings.  After the revert, `particleTemplates.C` again only
has the retained same-tet non-normalised inside-check diff from this candidate
sequence.

## Rejected candidate: `tetNeighbour(0)` early return

Candidate:

- File: `src/lagrangian/basic/particle/particleI.H`
- Function: `particle::tetNeighbour(label triI)`
- Idea: `triI == 0` represents crossing the real cell face.  The existing
  `case 0` intentionally does not modify `tetFaceI_` or `tetPtI_`, but the
  function still loaded face ownership, face geometry, base point, and adjacent
  point indices before entering the switch.  The candidate returned immediately
  for `triI == 0`.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_tetneigh_tri0_build_20260607.log
```

No-detail 10-step smoke:

```text
log.codex_mpi8_replicatedmesh_mpionly_tetneigh_tri0_smoke_10step_20260607
```

Smoke comparison:

| metric | same-tet no-detail baseline | tri0 early-return |
| --- | ---: | ---: |
| real | 10.86 s | 11.76 s |
| move+collide wall | 1.775343523 s | 1.831595762 s |
| move only | 1.223100367 s | 1.145058446 s |
| full evolve wall | 2.433862566 s | 2.616137885 s |

The smoke `move only` timer improved, but `real`, `move+collide wall`, and
`full evolve wall` regressed.  Because the code change was semantically narrow,
it was allowed one 500-step formal run.

Formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_tetneigh_tri0_500step_20260607
```

Formal comparison:

| metric | forced8 baseline | same-tet retained range | same-tet same-day confirm | tri0 early-return |
| --- | ---: | ---: | ---: | ---: |
| real | 132.29 s | 127.91-129.93 s | 139.66 s | 141.81 s |
| move only | 69.49945105 s | 67.56621687-69.11276280 s | 75.14025444 s | 79.32223967 s |
| move+collide wall | 90.87479808 s | 88.73592885-90.75870920 s | 99.01058420 s | 104.0407785 s |
| full evolve wall | 123.0267073 s | 118.7495988-120.5815986 s | 134.0601301 s | 137.5466661 s |
| particles max/min | 2.95295203 | 1.585453037-1.854821906 | 2.265839756 | 1.778463094 |
| DLB rebalances | 8 | 8 | 8 | 8 |

Decision: reject and revert.

The 500-step formal result regressed in both required metrics, including
against the same-day same-tet confirm run.  The source change was removed.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_tetneigh_tri0_revert_build_20260607.log
```

Revert build result: passed, with only the existing OpenFOAM template
instantiation warnings.  After the revert, `particleI.H` has no remaining diff
from this candidate.

## Rejected candidate: same-tet sqrt-free tolerance

Candidate:

- File: `src/lagrangian/basic/particle/particleTemplates.C`
- Function: retained same-tet non-normalised inside check
- Idea: keep the same non-normalised half-space check, but replace
  `d > SMALL*(mag(n) + VSMALL)` with a squared comparison to avoid the `sqrt`
  inside `mag(n)` on outside-plane failures.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_sametet_sqrtfree_build_20260607.log
```

No-detail 10-step smoke:

```text
log.codex_mpi8_replicatedmesh_mpionly_sametet_sqrtfree_smoke_10step_20260607
```

Smoke comparison:

| metric | same-tet no-detail baseline | same-tet sqrt-free |
| --- | ---: | ---: |
| real | 10.86 s | 14.46 s |
| move+collide wall | 1.775343523 s | 2.489695381 s |
| move only | 1.223100367 s | 1.506014525 s |
| full evolve wall | 2.433862566 s | 3.345659358 s |

Decision: reject and revert without 500-step validation.

The 10-step smoke regressed in every relevant timer, including `move only`, so
the candidate was not promoted to a formal run.

Revert build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_tracking_sametet_sqrtfree_revert_build_20260607.log
```

Revert build result: passed, with only the existing OpenFOAM template
instantiation warnings.  After the revert, `particleTemplates.C` again only has
the retained same-tet non-normalised inside-check diff.

## Current same-tet retained retest

After all rejected candidates above were reverted and rebuilt, the retained
same-tet source was retested without further source changes.

Formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_sametet_retest_500step_20260607
```

Result:

| metric | historical same-tet run 1 | historical same-tet confirm | same-day confirm | retained retest |
| --- | ---: | ---: | ---: | ---: |
| real | 127.91 s | 129.93 s | 139.66 s | 139.34 s |
| move only | 67.56621687 s | 69.11276280 s | 75.14025444 s | 76.67553366 s |
| move+collide wall | 88.73592885 s | 90.75870920 s | 99.01058420 s | 100.2404160 s |
| full evolve wall | 118.7495988 s | 120.5815986 s | 134.0601301 s | 135.0886023 s |
| post fields/output | 36.72935704 s | 36.30228583 s | 40.38098607 s | 41.91755684 s |
| particles max/min | 1.585453037 | 1.854821906 | 2.265839756 | 2.759808264 |
| DLB rebalances | 8 | 8 | 8 | 8 |

Interpretation: the clean retained retest confirms that the current runtime
environment / run outcome is slower than the earlier retained same-tet formal
pair.  The source state is still the retained same-tet version, but the latest
two same-day retained runs are both around `139 s` real with higher move,
post/output, and final particle imbalance than the historical `127.91-129.93 s`
pair.

## Clean rebuild retest after reverting local same-tet diff

The local `same-tet non-normalised` diff was later removed, leaving
`particleTemplates.C` on the original `tet.inside(endPosition)` fast path.  The
clean source was rebuilt and tested once with the same 500-step pure-MPI case.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_clean_retained_rebuild_20260607.log
```

Formal log:

```text
log.codex_mpi8_replicatedmesh_mpionly_clean_rebuild_retest_500step_20260607
```

Result:

| metric | clean rebuild retest | same-tet same-day retest | historical same-tet best |
| --- | ---: | ---: | ---: |
| real | 145.79 s | 139.34 s | 127.91 s |
| move only | 83.03607914 s | 76.67553366 s | 67.56621687 s |
| move+collide wall | 107.2331797 s | 100.2404160 s | 88.73592885 s |
| full evolve wall | 140.4947836 s | 135.0886023 s | 118.7495988 s |
| post fields/output | 41.59040581 s | 41.91755684 s | 36.72935704 s |
| particles max/min | 1.514010367 | 2.759808264 | 1.585453037 |
| DLB rebalances | 8 | 8 | 8 |

Interpretation: clean source is slower than the same-day same-tet retest by
about `6.45 s` real and `6.36 s` move-only.  This confirms that the retained
same-tet non-normalised fast path was still beneficial, even though the same-day
runtime environment was slower than the earlier historical retained runs.

A follow-up control-confirm rerun used the same clean source and forced8
control:

```text
log.codex_mpi8_replicatedmesh_mpionly_control_confirm_rerun_500step_20260607
```

| metric | clean rebuild retest | control-confirm rerun |
| --- | ---: | ---: |
| real | 145.79 s | 144.89 s |
| move only | 83.03607914 s | 83.37830902 s |
| move+collide wall | 107.2331797 s | 108.1305460 s |
| full evolve wall | 140.4947836 s | 142.6358027 s |
| post fields/output | 41.59040581 s | 44.51808941 s |
| particles max/min | 1.514010367 | 1.447143374 |
| DLB rebalances | 8 | 8 |

## Reference comparison caveat

The preserved `ourmeshbkp` log remains the external wall-time reference, but its
internal profile fields are not emitted with the same code path as the current
logs.

Evidence from
`run/.../ourmeshbkp/mpi8replicatedmesh/log.mpi8replicatedmesh.confirm_pdFalse_20260605_024610`:

- the profile output is the older detailed rank-format block, including
  `Move profiling summary`, `Evolve profiling summary`, and
  `Cloud name: dsmc [rank 2]`;
- `move only [s] = 47.1652636`, but the same block also reports
  `auto DLB/rebalance [s] = 15.78659015` and `full evolve wall [s] =
  114.9207628`;
- the replicated summary reports rank evolve max `99.13380862 s`, matching the
  detailed rank block's `total profiled [s] = 99.12887323`;
- the preserved controlDict has `useOpenMP false` but `openmpThreads 8`.

Current pure-MPI logs print the newer `DSMC solver profile summary` after a
raw-MPI max reduction of stage timers in `dsmcCloud::printProfileSummary()`.
Therefore the external wall time remains comparable, but the `move only` field
is not guaranteed to be a like-for-like timer against `ourmeshbkp`.  Candidate
decisions in this round still use the current 500-step move metric to avoid
retaining changes that merely improve `real` while making the current move path
worse.

## Current state after this note

Formal case controls are restored:

```text
endTime 5.e-05;
profileDetail false;
useOpenMP false;
openmpThreads 1;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
replicatedMeshDLBDualConstraint true;
replicatedMeshDLBInitialAlpha 0.8;
```

Retain:

- forced8 DLB configuration
- `particleTemplates.C` same-tet non-normalised inside check
- default-off move-detail diagnostic instrumentation

Reject:

- `particleI.H` quad/hexa edge fast path in `crossEdgeConnectedFace()`
- `particleTemplates.C` fixed tris array in DSMC `trackToFace(..., true)`
- temporary tri-split tracking-loop diagnostic is not retained
- `particleI.H` same-face `++tetPtI_` / `--tetPtI_` increment candidate
- `particleTemplates.C` static `tetLambda` helper candidate
- `particleTemplates.C` direct same-tet normal calculation candidate
- `particleTemplates.C` outside-plane tri prefilter candidate
- `particleI.H` `tetNeighbour(0)` early-return candidate
- `particleTemplates.C` same-tet sqrt-free tolerance candidate
