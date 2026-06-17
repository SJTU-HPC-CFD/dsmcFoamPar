# ourmesh mpi8 alpha / DualConstraint forced8 A-B

Date: 2026-06-12

Reason for rerun:

The previous `ourmesh_mpi8_alpha_dual_ab_20260612` comparison only triggered
one DLB event in 500 steps.  That is not the strict DLB benchmark state described
in `doc/worklog/v2506/dsmcFoam_plus_mpi_replicated_mesh_dlb_technical_report_20260608.md`,
where the retained formal comparison uses about eight DLB rebalances in 500
steps, normally via the forced8 schedule.

Current discrepancy:

- Current case before this rerun had `replicatedMeshDLBForceSteps ();`.
- Current source also has `replicatedMeshDLBParticleGate` defaulting to `true`.
- In non-forced mode, particle gate can skip expensive global checks once final
  particle max/min is below the threshold, and with `profileDetail false` this is
  silent in the normal log.
- Therefore the previous one-rebalance alpha comparison was a non-forced,
  particle-gated DLB test, not the strict forced8 DLB test.

Case:

`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh`

Common controls:

- `useOpenMP false`
- `openmpThreads 1`
- `endTime 5e-05` / 500 steps
- `profileSummary true`
- `profileDetail false`
- `replicatedMeshNoAlltoall true`
- `replicatedMeshFlatTransfer true`
- `replicatedMeshDLBForceSteps (120 170 220 270 320 370 420 470)`
- `replicatedMeshDLBTriggerMode legacyWindow`
- `replicatedMeshDLBMinGapSteps 50`
- `replicatedMeshDLBCheckCollective allgather`
- `replicatedMeshDLBImbalanceThreshold 1.5`
- `replicatedMeshDLBVsizeExp 0`

Trigger verification:

| Variant | alpha | DualConstraint | forced trigger count | summary rebalances | trigger steps |
|---|---:|---|---:|---:|---|
| forced8_alpha1_dualfalse | 1.0 | false | 8 | 8 | 120, 170, 220, 270, 320, 370, 420, 470 |
| forced8_alpha08_dualfalse | 0.8 | false | 8 | 8 | 120, 170, 220, 270, 320, 370, 420, 470 |
| forced8_alpha08_dualtrue | 0.8 | true | 8 | 8 | 120, 170, 220, 270, 320, 370, 420, 470 |

Performance:

| Variant | real [s] | full evolve [s] | move+collide [s] | move [s] | build occ [s] | collision [s] | migration max [s] | particles max/min | rank wall max/min | DLB wall max [s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| forced8_alpha1_dualfalse | 77.85 | 71.69948578 | 70.42431176 | 57.15336762 | 6.064381571 | 4.093368093 | 3.384457675 | 1.076510267 | 1.306471816 | 17.64633934 |
| forced8_alpha08_dualfalse | 77.79 | 75.00281689 | 73.80264439 | 59.22950381 | 8.836414037 | 4.943804734 | 3.128942142 | 1.063836386 | 1.291486858 | 17.92789238 |
| forced8_alpha08_dualtrue | 85.27 | 82.69961835 | 81.52390108 | 65.12548298 | 10.37490256 | 4.409779295 | 3.209079645 | 1.221412599 | 1.22995283 | 16.17028558 |

Correctness gate:

| Variant | final particles | final total energy | cumulative global collisions | stuck particles | exit |
|---|---:|---:|---:|---:|---:|
| forced8_alpha1_dualfalse | 2463765 | 1.24055354 | 9401450 | 0 | 0 |
| forced8_alpha08_dualfalse | 2463895 | 1.242222512 | 9588033 | 0 | 0 |
| forced8_alpha08_dualtrue | 2463575 | 1.241877649 | 9584130 | 0 | 0 |

Interpretation:

1. The report expectation is correct for the strict forced8 benchmark.  With
   forced steps restored, all three variants perform eight DLB rebalances.

2. The previous one-rebalance run should be treated as invalid for judging
   "normal DLB" alpha/dual behavior.

3. Under forced8, `DualConstraint true` is still slower in this current source:
   it is +7.4 s real versus `alpha=0.8, Dual=false` and +7.4 s versus
   `alpha=1, Dual=false`.  It improves rank-wall max/min, but this does not
   pay back in full evolve or real time.

4. `alpha=1, Dual=false` and `alpha=0.8, Dual=false` are close in external
   real time, but internal profile favors `alpha=1`: full evolve is 71.70 s
   versus 75.00 s, and move/build/collision are all lower.

5. For formal 500-step DLB comparisons, use forced8 or explicitly disable the
   particle gate.  Otherwise non-forced auto/SAR plus default particle gate can
   under-trigger and make DLB look artificially cheap.

Artifacts:

- Control snapshots: `controlDicts/`
- Logs: `logs/`
- Restored case controlDict snapshot: `controlDicts/ourmesh_mpi8_controlDict_restored`
