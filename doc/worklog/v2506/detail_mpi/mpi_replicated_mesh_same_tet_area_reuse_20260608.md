# MPI replicated mesh same-tet area reuse - 2026-06-08

## Scope

- Mode: pure MPI8 replicated mesh DLB only, no OpenMP.
- Case:
  `run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`
- Source state before this candidate: accepted sparse boundary clean plus the
  retained same-tet non-normalised helper.
- Candidate file:
  `src/lagrangian/basic/particle/particleTemplates.C`

Strict full-fields hashes used for formal runs:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2  system/controlDict
73a93c7487ec2a5e54c6a89f30d886207432bdc5ba9db1a5348e072423b55be9  system/fieldPropertiesDict
```

## Rationale

The retained DSMC same-tet fast path computes the current tet's four
unnormalised area vectors to test whether `endPosition` remains inside the
same tet.  When that fast path fails, the first iteration of the legacy
tet-walk loop is still on the same tet and previously recomputed the same four
area vectors before normalising them for `findTris()` and `tetLambda()`.

This candidate reuses the failed same-tet fast-path area vectors for that first
non-fast tet-walk iteration.  It is not a cross-step cache, transition table,
geometry cache, outside-plane screening replacement, or barycentric partial
tracking change.  It preserves the existing `findTris()` / outer
`tetLambda()` selection semantics.

## Source Change

- Split `dsmcTrackInsideTetNoNormalise()` so it can consume precomputed
  unnormalised tet area vectors.
- In the static, no-wall-impact same-tet fast path, compute
  `initialTetAreas[0..3]` once.
- If the same-tet check succeeds, return as before.
- If it fails, mark those areas reusable.
- In the first non-fast tet-walk loop, reuse `initialTetAreas` when
  `triI == -1`; later loops still compute fresh areas after topology changes.

Build log:

```text
doc/worklog/v2506/detail_mpi/stage_mpi_same_tet_area_reuse_build_20260608.log
```

Build result: passed, with only the known OpenFOAM-v1706 template warnings.

## 10-step Smoke

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_same_tet_area_reuse_smoke_20260608
```

Temporary controls:

```text
endTime 1.e-06;
profileSummary true;
profileDetail true;
moveDetailProfile true;
```

Log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_same_tet_area_reuse_smoke_10step_20260608
```

| metric | value |
| --- | ---: |
| real | 10.83 s |
| move only | 1.13728321 s |
| post fields/output | 0.333113415 s |
| full evolve wall | 1.850278409 s |
| move detail track calls | 25,758,229 |
| same-tet no-face | 11,139,519 |
| internal tet only | 9,482,988 |
| face hits | 5,135,722 |
| track max | 0.765096637 s |
| boundary max | 0.001273451 s |

Correctness:

```text
Total Iterations = 10
Number of DSMC particles = 2065865
Number of stuck particles = 0
Collisions = 2853
Total energy = 1.150631021
OpenMP enabled = 0
```

No `Fatal`, `NaN`, or `BAD TERMINATION` string was found.

## 200-step Signal

Control backup:

```text
doc/worklog/v2506/detail_mpi/mpi8replicatedmesh_controlDict_before_same_tet_area_reuse_signal_20260608
```

Temporary controls:

```text
endTime 2.e-05;
profileSummary true;
profileDetail true;
moveDetailProfile true;
```

Log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_same_tet_area_reuse_signal_200step_20260608
```

| metric | post-sparse structure diagnostic | same-tet area reuse |
| --- | ---: | ---: |
| real | 49.77 s | 47.47 s |
| move only | 28.04566287 s | 27.50361695 s |
| buildCellOccupancy | 2.348702429 s | 2.361648577 s |
| collision phase | 4.791171992 s | 4.901541412 s |
| post fields/output | 6.769591856 s | 6.663389341 s |
| full evolve wall | 41.58454539 s | 40.73366843 s |
| DLB rebalances | 2 | 2 |

Detail metrics:

| metric | value |
| --- | ---: |
| move detail track calls | 533,638,673 |
| same-tet no-face | 232,602,894 |
| internal tet only | 195,275,119 |
| face hits | 105,760,660 |
| patch hits | 439,068 |
| track max | 15.82798699 s |
| boundary max | 0.043582848 s |

Correctness:

```text
Total Iterations = 200
Number of DSMC particles = 2220288
Number of stuck particles = 0
Collisions = 15259
Total energy = 1.169563047
OpenMP enabled = 0
```

The signal was positive in aggregate `move only` and end-to-end wall time, but
the detail track timer was noisy.  The candidate therefore advanced to strict
500-step formal validation rather than being accepted from signal evidence.

## 500-step Formal

Formal controls were restored before both runs:

```text
endTime 5.e-05;
profileSummary true;
profileDetail false;
```

Formal log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_same_tet_area_reuse_formal_500step_20260608
```

Formal repeat log:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_same_tet_area_reuse_formal_repeat_500step_20260608
```

Comparison against the current strict full-fields baseline and the accepted
sparse-clean repeat:

| metric | strict baseline | sparse-clean repeat | area reuse formal | area reuse repeat |
| --- | ---: | ---: | ---: | ---: |
| real | 133.83 s | 113.24 s | 104.32 s | 107.45 s |
| move only | 72.9567669 s | 70.66798424 s | 67.43927759 s | 67.84144046 s |
| buildCellOccupancy | 6.729516136 s | 7.264959245 s | 6.65048193 s | 6.342166926 s |
| collision phase | 14.8398869 s | 14.50986947 s | 14.09394481 s | 14.86660294 s |
| post fields/output | 40.37082178 s | 18.95231379 s | 17.64348673 s | 17.2290821 s |
| full evolve wall | 132.7056254 s | 110.5389487 s | 104.4499775 s | 105.2572367 s |
| DLB rebalances | 8 | 8 | 8 | 8 |

Formal correctness:

```text
Total Iterations = 500
Number of DSMC particles = 2463520
Number of stuck particles = 0
Collisions = 33374
Total energy = 1.238430851
OpenMP enabled = 0
```

Formal repeat correctness:

```text
Total Iterations = 500
Number of DSMC particles = 2463734
Number of stuck particles = 0
Collisions = 34451
Total energy = 1.239092344
OpenMP enabled = 0
```

No `Fatal`, `NaN`, or `BAD TERMINATION` string was found in either formal log.

## Decision

Retain the same-tet area reuse candidate.

Reason:

- Both strict 500-step full-fields formal runs beat the current strict
  baseline in `real`, `move only`, and `full evolve wall`.
- Both formal runs also beat the accepted sparse-clean repeat in `real`,
  `move only`, `post fields/output`, and `full evolve wall`.
- The change is local to the static/no-wall first tet-walk iteration after a
  same-tet fast-path failure and does not alter `findTris()` or outer
  `tetLambda()` selection.

Post-candidate strict full-fields observed range:

```text
real 104.32-107.45 s
move only 67.43927759-67.84144046 s
post fields/output 17.2290821-17.64348673 s
full evolve wall 104.4499775-105.2572367 s
```

This remains separate from the historical `ourmeshbkp` reference
(`real 117.31 s`, `move only 47.1652636 s`) and from no-fields optional mode.

## Post-run Checks

- Case hashes are restored to the strict full-fields values shown at the top.
- `git diff --check` over `particleTemplates.C` and the case `controlDict`
  passed.
- The full-repo `git diff --check` is still expected to report the unrelated
  pre-existing EOF blank-line issue in `doc/worklog/v2506/dsmcFoam++tips.md`;
  that file was not touched by this candidate.

## Rollback Rebuild Retest

After the later barycentric-tracking experiment was rolled back to the
area-reuse source state, the first rebuild attempt failed because generated
`lnInclude` links under `src/lagrangian/basic/lnInclude` still pointed at
`particle/bkp` barycentric files.  The generated links for
`particle.C`, `particleI.H`, and `particleTemplates.C` were restored to the
main `particle/` directory, the stale `Make/linux64IccDPInt32Opt` directories
for `liblagrangian+`, `libdsmcFoam+`, and `dsmcFoam+` were removed, and the
solver was rebuilt successfully.

Build log:

```text
doc/worklog/v2506/detail_mpi/build_area_reuse_retest_20260608.log
```

Retest controls:

```text
controlDict sha256 = 8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
endTime 5.e-05;
deltaT 1.e-07;
profileSummary true;
profileDetail false;
useOpenMP false;
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

Retest logs:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_area_reuse_rebuild_retest_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_area_reuse_rebuild_retest_repeat_500step_20260608
```

Comparison including the rollback rebuild retest:

| metric | area reuse formal | area reuse repeat | rebuild retest | rebuild retest repeat |
| --- | ---: | ---: | ---: | ---: |
| real | 104.32 s | 107.45 s | 108.69 s | 116.04 s |
| move+collide wall | 89.86885435 s | 91.06343871 s | 84.97968617 s | 94.10253623 s |
| move only | 67.43927759 s | 67.84144046 s | 62.51511822 s | 70.26858009 s |
| buildCellOccupancy | 6.65048193 s | 6.342166926 s | 6.482688895 s | 7.243594039 s |
| collision phase | 14.09394481 s | 14.86660294 s | 13.66598807 s | 15.21113079 s |
| post fields/output | 17.64348673 s | 17.2290821 s | 17.10465238 s | 18.14279997 s |
| full evolve wall | 104.4499775 s | 105.2572367 s | 100.2458092 s | 108.3720354 s |
| real - full evolve | -0.1299775 s | 2.1927633 s | 8.4441908 s | 7.6679646 s |
| migration wall | 2.084649607 s | 2.156514815 s | 2.325213275 s | 2.33008606 s |
| particles per rank max/min | 1.691773706 | 1.609575878 | 2.047642043 | 2.092819679 |
| rank wall max/min | 1.03494898 | 1.035808851 | 1.024034227 | 1.045664207 |
| DLB rebalances | 8 | 8 | 8 | 8 |

Rollback retest correctness:

```text
Total Iterations = 500
Number of DSMC particles = 2463891
Number of stuck particles = 0
Collisions = 34150
Total energy = 1.238664777
OpenMP enabled = 0
```

Rollback repeat correctness:

```text
Total Iterations = 500
Number of DSMC particles = 2463653
Number of stuck particles = 0
Collisions = 34748
Total energy = 1.239284822
OpenMP enabled = 0
```

No `Fatal`, `NaN`, `BAD TERMINATION`, `Segmentation`, or `abort` string was
found in either rollback retest log.

Retest conclusion:

- The rollback source state still builds and runs the strict 500-step case
  correctly after stale generated `lnInclude` links are repaired.
- The first rollback retest had a faster profiled core than the original
  area-reuse formal run (`full evolve 100.2458092 s` versus `104.4499775 s`,
  `move only 62.51511822 s` versus `67.43927759 s`), but external `real`
  was worse because `real - full evolve` grew to `8.4441908 s`.
- The repeat run was slower in both `real` and profiled core
  (`real 116.04 s`, `full evolve 108.3720354 s`).
- Therefore this rollback rebuild retest does not show stable reproduction of
  the earlier `real 104.32-107.45 s` range, even though the first retest still
  shows the area-reuse move kernel can be faster than the earlier formal logs.

## Backup Suffix Hardening

The remaining backup files under `src/lagrangian/basic/particle/bkp` were
renamed so that no file ends in `.C` or `.H`.  This prevents
`wmake -j lnInclude` from re-linking backup barycentric or backup particle
implementation files into `src/lagrangian/basic/lnInclude`.

The only name collision was `particleTemplates.C`, because
`particleTemplates.Cbkp` already existed and the two files were different.
The `.C` file was therefore preserved as
`particleTemplates.barycentric.Cbkp`.

Verification:

```text
find src/lagrangian/basic/particle/bkp -maxdepth 1 -type f \( -name '*.C' -o -name '*.H' \)
```

returned no files.

After rerunning `wmake -j lnInclude` in `src/lagrangian/basic`, no generated
link points to `../particle/bkp/*`; the active particle links remain:

```text
lnInclude/particle.C -> ../particle/particle.C
lnInclude/particle.H -> ../particle/particle.H
lnInclude/particleI.H -> ../particle/particleI.H
lnInclude/particleTemplates.C -> ../particle/particleTemplates.C
```

## Particle Max/Min Interpretation

The rollback retests changed `particles per rank max/min` even though the
global final particle count stayed in the same range:

| run | final DSMC particles | final collisions | rank 0 final local particles | particles per rank max/min | rank wall max/min |
| --- | ---: | ---: | ---: | ---: | ---: |
| area reuse formal | 2463520 | 33374 | 217522 | 1.691773706 | 1.03494898 |
| area reuse repeat | 2463734 | 34451 | 356560 | 1.609575878 | 1.035808851 |
| rebuild retest | 2463891 | 34150 | 398594 | 2.047642043 | 1.024034227 |
| rebuild retest repeat | 2463653 | 34748 | 453288 | 2.092819679 | 1.045664207 |

This is not evidence that the rollback changed the global particle population;
the total particle counts differ only by normal stochastic/run-to-run amounts.
The larger max/min comes from the owner partition selected by the runtime DLB.

Relevant implementation details:

- `dsmcReplicatedMesh::report()` prints `particles per rank` from
  `allParticleCounts_`, which is last filled by `updateParticleCounts()` using
  each rank's current `cloud_.size()`.
- `updateParticleCounts()` is called during migration and after DLB migration;
  the report does not recompute a new balance objective by itself.
- `reassignByParMetisAdaptiveRepart()` builds ParMETIS vertex weights from the
  current `cellOccupancy()` at the trigger step, not from the final particle
  distribution at step 500.
- With the current controls, DLB uses dual constraints:
  `N^alpha` for move balance and `N*(N-1)` as a collision proxy, with
  `alpha=0.8`, `ubvec=1.05`, and the second constraint relaxed to
  `ubvec1=1.5`.
- The objective is therefore closer to move/collision wall-time balance than
  strict final particle-count equality.  This matches the logs: the two
  rollback retests have worse particle-count max/min (`2.05-2.09`) while rank
  wall-time max/min remains tight (`1.024-1.046`).
- The owner maps themselves differ between runs.  At the final forced DLB
  step 470, rank 0 owned `9708`, `12724`, `14479`, and `15731` cells in the
  four logs respectively, so the final local-particle extrema are not being
  produced by a fixed post-rollback partition.

Conclusion: the changed particle max/min is primarily a runtime DLB/ParMETIS
partition outcome driven by stochastic per-cell occupancy and the dual
move/collision weighting, not a direct source-state difference introduced by
the rollback.  It should be evaluated together with rank wall-time balance and
`full evolve wall`, not as a standalone correctness regression.

## Current Forced vs Auto-Trigger Retest

The suffix-hardened current area-reuse binary was retested with the strict
500-step full-fields case in two control modes:

1. Current forced-DLB mode:

```text
replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470);
```

2. Temporary auto-trigger mode:

```text
replicatedMeshDLBForceSteps ();
replicatedMeshDLBTriggerMode legacyWindow;
replicatedMeshDLBImbalanceThreshold 1.5;
```

The forced `controlDict` was backed up as:

```text
run/.../ourmesh/mpi8replicatedmesh/system/controlDict.codex_forced_backup_20260608_0905
```

After the auto-trigger runs, `controlDict` was restored to the original forced
hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Logs:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_current_forced_retest1_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_current_forced_retest2_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_current_forced_retest3_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_current_auto_dlb_retest1_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_current_auto_dlb_retest2_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_current_auto_dlb_retest3_500step_20260608
```

Per-run metrics:

| mode | run | real | full evolve | move | buildOcc | collision | post | migration wall | DLB | particles max/min | rank wall max/min | trigger steps |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| forced | 1 | 110.93 s | 102.3943661 s | 64.88674413 s | 7.227137593 s | 14.06330151 s | 18.863859 s | 2.361587943 s | 8 | 2.607125043 | 1.054816859 | 120,170,220,270,320,370,420,470 |
| forced | 2 | 109.48 s | 101.7065025 s | 64.48090461 s | 6.946939193 s | 13.30541235 s | 18.63083526 s | 2.176022909 s | 8 | 2.933420366 | 1.045288255 | 120,170,220,270,320,370,420,470 |
| forced | 3 | 109.26 s | 101.4888815 s | 64.17474708 s | 6.484722268 s | 13.7652318 s | 17.35575066 s | 2.186805419 s | 8 | 2.08765907 | 1.029473946 | 120,170,220,270,320,370,420,470 |
| auto | 1 | 121.83 s | 114.5789268 s | 75.69324276 s | 9.934031429 s | 16.39261879 s | 22.3487361 s | 2.038778871 s | 4 | 1.543833199 | 1.086112047 | 90,200,320,490 |
| auto | 2 | 109.93 s | 103.0170939 s | 66.11991646 s | 7.764739976 s | 14.30610954 s | 19.09327943 s | 1.556616604 s | 4 | 2.150588263 | 1.054087429 | 100,170,280,410 |
| auto | 3 | 115.69 s | 108.3855002 s | 71.16215016 s | 8.170085713 s | 15.48269779 s | 19.69242345 s | 1.768112806 s | 4 | 1.971021618 | 1.063735236 | 80,140,290,480 |

Three-run summary:

| metric | forced mean | forced range | auto mean | auto range | auto - forced mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| real | 109.89 s | 109.26-110.93 s | 115.816667 s | 109.93-121.83 s | +5.926667 s |
| full evolve | 101.86325 s | 101.4888815-102.3943661 s | 108.660507 s | 103.0170939-114.5789268 s | +6.797257 s |
| move | 64.5141319 s | 64.17474708-64.88674413 s | 70.9917698 s | 66.11991646-75.69324276 s | +6.477638 s |
| buildOcc | 6.88626635 s | 6.484722268-7.227137593 s | 8.62295237 s | 7.764739976-9.934031429 s | +1.736686 s |
| collision | 13.7113152 s | 13.30541235-14.06330151 s | 15.3938087 s | 14.30610954-16.39261879 s | +1.682494 s |
| post | 18.2834816 s | 17.35575066-18.863859 s | 20.3781463 s | 19.09327943-22.3487361 s | +2.094665 s |
| migration wall | 2.24147209 s | 2.176022909-2.361587943 s | 1.78783609 s | 1.556616604-2.038778871 s | -0.453636 s |
| DLB rebalances | 8 | 8-8 | 4 | 4-4 | -4 |
| particles max/min | 2.54273483 | 2.08765907-2.933420366 | 1.88848103 | 1.543833199-2.150588263 | -0.654254 |
| rank wall max/min | 1.04319302 | 1.029473946-1.054816859 | 1.06797824 | 1.054087429-1.086112047 | +0.024785 |

Correctness:

```text
All six logs completed Total Iterations = 500.
All six logs had Number of stuck particles = 0 at the final report.
No Fatal, NaN, BAD TERMINATION, Segmentation, or abort string was found.
```

Interpretation:

- The current forced-step policy is clearly faster in this 3-run A/B:
  `real` mean improves by `5.93 s` and `full evolve` mean improves by
  `6.80 s` versus auto-trigger.
- Auto-trigger cuts actual rebalances from 8 to 4 and saves about
  `0.45 s` of migration wall time on average, but that saving is much smaller
  than the added move/build/collision/post cost.
- Auto-trigger often reports better particle-count max/min, but worse
  rank-wall max/min.  For performance, the rank-wall and `full evolve`
  metrics are the more relevant indicators here.
- The auto-trigger steps are SAR-driven and run-to-run variable
  (`90/200/320/490`, `100/170/280/410`, `80/140/290/480`).  With this
  threshold/SAR configuration it rebalances too sparsely for this case.

Decision: keep the forced-step schedule as the current formal comparison
control.  Treat the unforced auto-trigger mode as slower for this case unless
the SAR/threshold policy is retuned.

## Threshold-Only Optimized Trigger Retest

The current area-reuse binary was then retested with the new runtime trigger
controls added to `dsmcReplicatedMesh::autoRebalance()`.  The new controls keep
the existing behavior unless explicitly set in `controlDict`; this test used:

```text
replicatedMeshDLBForceSteps ();
replicatedMeshDLBMinStartStep 120;
replicatedMeshDLBRemainingGuardSteps 20;
replicatedMeshDLBUseSAR false;
replicatedMeshDLBImbalanceThreshold 1.08;
```

Build log:

```text
doc/worklog/v2506/detail_mpi/build_dlb_trigger_controls_20260608.log
```

The forced `controlDict` was backed up before the temporary edit as:

```text
run/.../ourmesh/mpi8replicatedmesh/system/controlDict.codex_forced_backup_20260608_0922
```

After the optimized-trigger runs, `controlDict` was restored to the forced
schedule.  Current `controlDict` and the 09:22 backup have the same hash:

```text
8863b96bc6a149404bb5a70eebcfb1eaafb859984a83cc19816735a1460d09a2
```

Logs:

```text
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_optimized_trigger_retest1_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_optimized_trigger_retest2_500step_20260608
run/.../ourmesh/mpi8replicatedmesh/log.codex_mpi8_replicatedmesh_optimized_trigger_retest3_500step_20260608
```

Per-run metrics:

| mode | run | real | full evolve | move | buildOcc | collision | post | migration wall | DLB | particles max/min | rank wall max/min | trigger steps |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| optimized trigger | 1 | 110.58 s | 103.2169454 s | 66.33302837 s | 6.677452896 s | 14.02316395 s | 17.86551863 s | 2.206863526 s | 7 | 2.932479725 | 1.046480255 | 170,220,270,320,370,420,470 |
| optimized trigger | 2 | 107.76 s | 100.2996596 s | 63.31827823 s | 6.998684833 s | 13.30115936 s | 18.24323326 s | 1.848968246 s | 5 | 1.900177096 | 1.038339469 | 170,220,270,320,370 |
| optimized trigger | 3 | 109.92 s | 102.8316708 s | 65.76631816 s | 6.54054088 s | 13.64635924 s | 17.74926116 s | 2.188250273 s | 7 | 1.861853701 | 1.036753823 | 170,220,270,320,370,420,470 |

Three-mode summary:

| metric | forced mean | auto mean | optimized mean | optimized - forced | optimized - auto |
| --- | ---: | ---: | ---: | ---: | ---: |
| real | 109.89 s | 115.816667 s | 109.42 s | -0.47 s | -6.396667 s |
| full evolve | 101.86325 s | 108.660507 s | 102.116092 s | +0.252842 s | -6.544415 s |
| move | 64.5141319 s | 70.9917698 s | 65.1392083 s | +0.625076 s | -5.852562 s |
| buildOcc | 6.88626635 s | 8.62295237 s | 6.73889287 s | -0.147373 s | -1.884059 s |
| collision | 13.7113152 s | 15.3938087 s | 13.6568942 s | -0.054421 s | -1.736915 s |
| post | 18.2834816 s | 20.3781463 s | 17.952671 s | -0.330811 s | -2.425475 s |
| migration wall | 2.24147209 s | 1.78783609 s | 2.08136068 s | -0.160111 s | +0.293525 s |
| DLB rebalances | 8 | 4 | 6.333333 | -1.666667 | +2.333333 |
| particles max/min | 2.54273483 | 1.88848103 | 2.23150351 | -0.311231 | +0.343022 |
| rank wall max/min | 1.04319302 | 1.06797824 | 1.04052452 | -0.002669 | -0.027454 |

Optimized-trigger ranges:

| metric | range |
| --- | ---: |
| real | 107.76-110.58 s |
| full evolve | 100.2996596-103.2169454 s |
| move | 63.31827823-66.33302837 s |
| DLB rebalances | 5-7 |
| particles max/min | 1.861853701-2.932479725 |
| rank wall max/min | 1.036753823-1.046480255 |

Correctness:

```text
All three optimized-trigger logs completed Total Iterations = 500.
All three optimized-trigger logs had Number of stuck particles = 0 at the final report.
No Fatal, NaN, BAD TERMINATION, Segmentation, or abort string was found.
```

Interpretation:

- The optimized threshold-only trigger is clearly better than the old unforced
  auto/SAR policy for this case: mean `real` improves by `6.40 s`, mean
  `full evolve` improves by `6.54 s`, and rank-wall balance is tighter.
- Relative to forced8 it is very close.  Mean external `real` is `0.47 s`
  lower, but mean `full evolve` is `0.25 s` higher and mean `move` is `0.63 s`
  higher.  The result therefore does not clearly beat forced8 under the stricter
  solver-profile criterion.
- The optimized trigger consistently skips step 120.  It then follows the
  forced cadence while the threshold is exceeded; one run also skipped 420 and
  470, giving the best single result (`real 107.76 s`, `full evolve
  100.2996596 s`).

Decision: use the new trigger controls as a useful tuning path, but keep the
forced8 schedule as the current formal comparison control until a tuned policy
beats forced8 in both `real` and solver-profile metrics.

## Trigger-Control Rollback To Area-Reuse Auto Config

After reviewing the long-run DSMC implications, the threshold-only trigger
experiment was rolled back.  The reason is that `threshold=1.08`,
`minGapSteps=50`, and fixed 500-step guards are too case-specific for production
DSMC runs that may cover thousands of steps.

Rollback actions:

```text
src/lagrangian/dsmc/replicatedMesh/dsmcReplicatedMesh.C:
  removed replicatedMeshDLBMinStartStep / RemainingGuardSteps / UseSAR runtime controls
  restored the original currentStep < 30 guard
  restored SAR triggering to sar_ > 0.0

run/.../ourmesh/mpi8replicatedmesh/system/controlDict:
  replicatedMeshDLBForceSteps ();
  replicatedMeshDLBTriggerMode legacyWindow;
  replicatedMeshDLBMinGapSteps 50;
  replicatedMeshDLBCheckCollective allgather;
  replicatedMeshDLBImbalanceThreshold 1.5;
```

The forced8 `controlDict` was backed up before this rollback as:

```text
run/.../ourmesh/mpi8replicatedmesh/system/controlDict.codex_forced_backup_before_trigger_revert_20260608_rollback
```

Build/verification:

```text
doc/worklog/v2506/detail_mpi/build_dlb_trigger_revert_20260608.log
```

- Build completed successfully.
- `libdsmcFoam+.so` was rebuilt at `2026-06-08 10:08`.
- `dsmcReplicatedMesh.C` has no remaining diff from the trigger-control
  experiment.
- Current `controlDict` contains no `replicatedMeshDLBMinStartStep`,
  `replicatedMeshDLBRemainingGuardSteps`, or `replicatedMeshDLBUseSAR` entries.

Current state: area-reuse source plus non-forced legacy auto/SAR DLB trigger
configuration.  This is not a new benchmark result; no 500-step retest was run
after the rollback.
