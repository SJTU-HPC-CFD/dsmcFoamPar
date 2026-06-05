# Stage 2 replicated mesh runtime port - 2026-06-05

## Scope

This note records the first runtime integration slice for the OFv1706
replicated mesh path.  The earlier stage-2 build work had already made
`dsmcReplicatedMesh` compile and link; this slice connects it to `dsmcCloud`
and `dsmcFoam+`.

## Source changes

- Added `autoPtr<dsmcReplicatedMesh> replicatedMesh_` to `dsmcCloud`.
- Added `replicatedMeshActive()`, `isOutputRank()`, and replicated mesh accessors.
- Construct `dsmcReplicatedMesh` when `controlDict` has `replicatedMesh true`.
- Integrated the runtime order in `dsmcCloud::evolve_moveAndCollide()`:
  - initial distribution before the first move;
  - post-move migration by `cellOwner`;
  - build occupancy after migration;
  - collision;
  - replicated auto-DLB check after collision controls.
- Updated `dsmcFoam+` write path:
  - gather parcels to rank 0 at output time;
  - only output rank writes;
  - migrate parcels back to owner after output;
  - call `Pstream::exit(0)` in `main()` when the replicated path initialized
    MPI without the normal `-parallel` `ParRunControl`.
- Made replicated flat transfer the default because the stream fallback fails
  on OFv1706 `IStringStream` binary parcel reads.
- Added flat `TransferData` gather for output-time parcel collection.
- Added replicated output-time MPI allreduce for the main `dsmcVolFields`
  cumulative volume/boundary arrays before derived fields are recomputed.

## Validation

Builds:

- `stage2_replicated_mesh_runtime_build_20260605.log`: passed.
- `stage2_replicated_mesh_finalize_build_20260605.log`: passed.
- `stage2_replicated_mesh_fields_reduce_mpi_build_20260605.log`: passed.
- `stage2_replicated_mesh_output_timing_build_20260605.log`: passed.

Smoke tests:

- `stage2_runtime_smoke_omp8_20260605.log`: passed, one serial step.
- `stage2_runtime_smoke_mpi8replicatedmesh_np2_finalize_20260605.log`: passed,
  `mpirun -np 2 dsmcFoam+` without `-parallel`.
- `stage2_runtime_smoke_mpi8replicatedmesh_np2_after_fieldreduce_20260605.log`:
  passed after the field-reduce changes.
- `stage2_runtime_smoke_mpi8replicatedmesh_np2_final_20260605.log`: passed
  after adding output-segment timing.

Important observed output from the latest replicated smoke:

- `Replicated mesh: flatTransfer=1 found=0`
- `migration calls = 2`
- `particles per rank = min 579067 max 1479501 max/min 2.554973777`
- `rank wall time (evolve, 1 steps): min 0.43 max 0.43 max/min 1`
- `Total Iterations = 1`
- `End main`

## Fixed failures

- The initial stream-transfer smoke failed with:
  - `Expected a '(' while reading binaryBlock`
  - log: `stage2_runtime_smoke_mpi8replicatedmesh_np2_20260605.log`
- Enabling flat transfer fixed the parcel migration failure.
- After the first flat-transfer success, mpirun still ended with `SIGKILL`
  after the solver summary because the manually initialized MPI path did not
  finalize through OpenFOAM.  Adding the `Pstream::exit(0)` branch in `main()`
  fixed this; `stage2_runtime_smoke_mpi8replicatedmesh_np2_finalize_20260605.log`
  exits with status 0.

## Remaining limitation

The forced output-time stress smoke with `writeInterval 1.e-07` did not finish
within 60 seconds on the full 2.06M-parcel case:

- `stage2_runtime_smoke_mpi8replicatedmesh_np2_output_mpi_20260605.log`
- status `124` from the timeout wrapper.

The log reached the post-evolve `dsmc.info()` output and then spent the timeout
window in the replicated output/write segment.  This is not a correctness crash,
but output-time gather/write for the full parcel cloud still needs a smaller
targeted case or a longer dedicated validation run.
After this timeout, the solver output path was instrumented to print:

- `gather parcels [s]`
- `rank0 write [s]`
- `migrate-back [s]`

## Next steps

1. Run a smaller replicated output case, or let the full output stress run
   complete without a 60 second timeout, to validate gather/write/cellOwner
   output end to end.
2. Add timing around replicated output gather, rank-0 write, and migrate-back.
3. Continue with move-side OpenMP data flow and flat occupancy only after the
   replicated output path has one completed output-time validation.
