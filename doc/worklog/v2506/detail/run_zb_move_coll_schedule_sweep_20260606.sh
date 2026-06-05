#!/usr/bin/env bash

ROOT="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb"
CASE_DIR="$ROOT/run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react/omp8"
CONTROL="$CASE_DIR/system/controlDict"
RUN_DIR="$CASE_DIR/schedule_sweep_stage26_20260606"
CHUNKS=(8 16 32 64 128 256)
SCHEDULES=(static dynamic)

source "$ROOT/doc/scripts/env.sh"
set -u

if [ -e "$RUN_DIR/manifest.tsv" ]; then
    echo "Refusing to overwrite existing sweep manifest: $RUN_DIR/manifest.tsv" >&2
    exit 2
fi

mkdir -p "$RUN_DIR"

BACKUP="$CONTROL.codex_bak_20260606_before_move_coll_schedule_sweep"
if [ -e "$BACKUP" ]; then
    BACKUP="$CONTROL.codex_bak_20260606_before_move_coll_schedule_sweep_$(date +%H%M%S)"
fi

cp "$CONTROL" "$BACKUP"
cp "$CONTROL" "$RUN_DIR/controlDict.before_sweep"

restore_control_dict()
{
    cp "$BACKUP" "$CONTROL"
}
trap restore_control_dict EXIT

printf 'group\tmoveSchedule\tmoveChunk\tcollisionSchedule\tcollisionChunk\texitCode\telapsedSeconds\tlog\tcontrolDict\n' > "$RUN_DIR/manifest.tsv"

set_control()
{
    local move_sched="$1"
    local move_chunk="$2"
    local coll_sched="$3"
    local coll_chunk="$4"

    foamDictionary "$CONTROL" -entry openmpMoveSchedule -set "$move_sched" >/dev/null
    foamDictionary "$CONTROL" -entry openmpMoveChunk -set "$move_chunk" >/dev/null
    foamDictionary "$CONTROL" -entry openmpCollisionSchedule -set "$coll_sched" >/dev/null
    foamDictionary "$CONTROL" -entry openmpCollisionChunk -set "$coll_chunk" >/dev/null
}

run_case()
{
    local group="$1"
    local move_sched="$2"
    local move_chunk="$3"
    local coll_sched="$4"
    local coll_chunk="$5"
    local label="$group"_move-"$move_sched"-"$move_chunk"_coll-"$coll_sched"-"$coll_chunk"
    local log="$RUN_DIR/log.codex_zb_${label}_300step_20260606"
    local ctl="$RUN_DIR/controlDict.${label}"
    local start_stamp
    local elapsed
    local status

    set_control "$move_sched" "$move_chunk" "$coll_sched" "$coll_chunk"
    cp "$CONTROL" "$ctl"

    start_stamp="$(date '+%F %T')"
    echo "START $start_stamp $label"

    SECONDS=0
    (
        cd "$CASE_DIR"
        OMP_NUM_THREADS=8 /usr/bin/time -p dsmcFoam+
    ) > "$log" 2>&1
    status=$?
    elapsed=$SECONDS

    echo "END $(date '+%F %T') $label exit=$status elapsed=${elapsed}s"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$group" "$move_sched" "$move_chunk" "$coll_sched" "$coll_chunk" \
        "$status" "$elapsed" "$log" "$ctl" >> "$RUN_DIR/manifest.tsv"
}

# Move schedule sweep: isolate move scheduling while keeping collision at the
# current stage25 baseline, dynamic chunk 8.
for sched in "${SCHEDULES[@]}"; do
    for chunk in "${CHUNKS[@]}"; do
        run_case "moveSweep" "$sched" "$chunk" dynamic 8
    done
done

# Collision schedule sweep: isolate collision scheduling while keeping move at
# the current stage25 baseline, static chunk 64.
for sched in "${SCHEDULES[@]}"; do
    for chunk in "${CHUNKS[@]}"; do
        run_case "collisionSweep" static 64 "$sched" "$chunk"
    done
done

restore_control_dict
trap - EXIT
cp "$CONTROL" "$RUN_DIR/controlDict.after_restore"
echo "Sweep complete: $RUN_DIR"
