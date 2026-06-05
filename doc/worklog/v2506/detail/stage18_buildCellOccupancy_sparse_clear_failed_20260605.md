# Stage 18: buildCellOccupancy sparse clear attempt, rejected

Date: 2026-06-05

## Scope

This was the first `buildCellOccupancy()` follow-up after the kept stage17 post
field optimization.

Attempted changes:

- Reuse a member `occupancyTotalCounts_` instead of constructing
  `labelList totalCounts(nCells, 0)` per call.
- Clear per-thread cell counts only for the previous call's active cells.
- Remove the full `nextCellOffsets(occupancyCellOffsets_)` copy by assigning
  per-thread write offsets from active cells.

Touched files during the attempt:

- `src/lagrangian/dsmc/clouds/dsmcCloud.H`
- `src/lagrangian/dsmc/clouds/dsmcCloud.C`

## Build and smoke

Build log:

```text
doc/worklog/v2506/detail/stage18_buildCellOccupancy_sparse_clear_build_20260605.log
```

Smoke log:

```text
run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/omp8/log.codex_omp_stage18_buildCellOccupancy_sparse_clear_smoke_10step_20260605
```

10-step smoke comparison:

| metric | stage17 kept smoke | stage18 smoke | delta |
|---|---:|---:|---:|
| real | 6.72 s | 6.89 s | 0.17 s slower |
| move only | 0.826445053 s | 0.846352675 s | 0.019907622 s slower |
| buildCellOccupancy | 0.106092180 s | 0.113830229 s | 0.007738049 s slower |
| collision phase | 0.019019326 s | 0.019111045 s | 0.000091719 s slower |
| post fields/output | 0.632757959 s | 0.627169476 s | 0.005588483 s faster |

Smoke correctness tail:

```text
Collisions                = 2827
Collision candidates      = 2942
Number of DSMC particles  = 2065855
Total energy              = 1.15062639
End main
```

## Decision

Rejected.  The target metric `buildCellOccupancy` regressed from
`0.106092180 s` to `0.113830229 s` in the 10-step smoke, and end-to-end smoke
also regressed.

The code was reverted to the kept stage17 state and rebuilt.

Revert build log:

```text
doc/worklog/v2506/detail/stage18_buildCellOccupancy_sparse_clear_revert_build_20260605.log
```

No 500-step run was performed for this rejected stage.
