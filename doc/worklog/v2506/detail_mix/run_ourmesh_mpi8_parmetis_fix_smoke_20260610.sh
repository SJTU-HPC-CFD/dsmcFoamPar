#!/usr/bin/env bash

set -u

ROOT="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb"
DETAIL="$ROOT/doc/worklog/v2506/detail_mix/mpi8_parmetis_fix_20260610"
CASE_DIR="$ROOT/run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8replicatedmesh"
CTRL="$CASE_DIR/system/controlDict"
BACKUP="$DETAIL/controlDict_before_ourmesh_mpi8_20260610"
CTRL_COPY="$DETAIL/controlDict_ourmesh_mpi8_320step_nowrite_20260610"
LOG="$DETAIL/ourmesh_MPI8_320step_parmetis_fix_20260610.log"
EXIT_FILE="$DETAIL/ourmesh_MPI8_320step_parmetis_fix_20260610.exit"

mkdir -p "$DETAIL"
cp "$CTRL" "$BACKUP"

restore_control()
{
    if [ -f "$BACKUP" ]
    then
        cp "$BACKUP" "$CTRL"
    fi
}

on_exit()
{
    local status=$?
    restore_control
    exit "$status"
}
trap on_exit EXIT

set_entry()
{
    local file="$1"
    local key="$2"
    local value="$3"

    if grep -Eq "^[[:space:]]*${key}[[:space:]]" "$file"
    then
        sed -i -E "s|^[[:space:]]*${key}[[:space:]].*;|${key} ${value};|" "$file"
    else
        printf "\n%s %s;\n" "$key" "$value" >> "$file"
    fi
}

set_entry "$CTRL" "startFrom" "startTime"
set_entry "$CTRL" "startTime" "0"
set_entry "$CTRL" "stopAt" "endTime"
set_entry "$CTRL" "endTime" "3.2e-05"
set_entry "$CTRL" "deltaT" "1.e-07"
set_entry "$CTRL" "writeControl" "runTime"
set_entry "$CTRL" "writeInterval" "1.e-3"
set_entry "$CTRL" "profileSummary" "true"
set_entry "$CTRL" "profileDetail" "false"
set_entry "$CTRL" "runTimeModifiable" "no"
cp "$CTRL" "$CTRL_COPY"

set +u
source "$ROOT/doc/scripts/env.sh"
set -u

(
    cd "$CASE_DIR"
    export OMP_NUM_THREADS=1
    echo "case=ourmesh"
    echo "mode=MPI8"
    echo "steps=320"
    echo "controlDict=$CTRL_COPY"
    echo "start=$(date -Iseconds)"
    /usr/bin/time -p /usr/bin/timeout 600 mpirun -np 8 dsmcFoam+
) > "$LOG" 2>&1
status=$?
printf "%s\n" "$status" > "$EXIT_FILE"

restore_control
trap - EXIT
exit "$status"
