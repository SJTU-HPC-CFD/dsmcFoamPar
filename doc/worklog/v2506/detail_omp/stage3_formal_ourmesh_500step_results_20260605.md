# v2506 500-Step Formal Results

Case root:
`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh`

Reference:
`run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmeshbkp`

## OMP8

Current:
- `real` 248.80 s
- `ClockTime` 249 s
- `full evolve wall` 1384.77 s
- `move only` 547.21 s
- `buildCellOccupancy` 194.51 s
- `collision phase` 473.59 s
- `post fields/output` 168.87 s
- collisions 35283
- particles 2463725
- total energy 1.243229348

History:
- `ClockTime` 104 s
- `full evolve wall` 102.1225033 s
- collisions 35847

Status:
- correct enough for 500-step run
- slower than history

## MPI8 Origin

Current:
- `real` 177.89 s
- `ClockTime` 177 s
- `full evolve wall` 174.26 s
- `move only` 132.82 s
- `buildCellOccupancy` 14.26 s
- `collision phase` 19.46 s
- `post fields/output` 31.64 s
- collisions 35616
- particles 2463565
- total energy 1.243249082

History:
- `full evolve wall` 163.0536101 s
- collisions 35691

Status:
- correctness close
- slower than history

## MPI8 Replicated Mesh

Current owner-filtered run:
- `real` 184.13 s
- `ClockTime` 184 s
- `full evolve wall` 174.31 s
- `move only` 83.71 s
- `buildCellOccupancy` 6.78 s
- `collision phase` 14.25 s
- `post fields/output` 75.75 s
- collisions 34774
- particles 2463897
- total energy 1.239659286
- local particles final 411539
- migration calls 57
- Phase C auto DLB rebalances 6

History:
- `full evolve wall` 114.9207628 s
- `migration calls` 59
- `local particles final` 296197
- `Phase C auto DLB rebalances` 8

Status:
- particle count no longer explodes
- still slower than history
- replicated DLB behavior still not fully matched

## Implementation Notes

- `Cloud::move()` now has a local OpenMP move path for non-MPI runs.
- `particle::trackToFace()` no longer uses shared cloud-level tri buffers.
- `dsmcCloud` now respects `openmpThreads` from `controlDict`.
- replicated mesh now filters new parcel insertion by owning rank.
