#!/usr/bin/env bash

set -u

ROOT="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb"
DETAIL="$ROOT/doc/worklog/v2506/detail_mix"
OUTDIR="$DETAIL/repeat_nowrite_compute_20260610_rerun"
LOGDIR="$OUTDIR/logs"
CTRLDIR="$OUTDIR/controlDicts"
MANIFEST="$OUTDIR/run_manifest.tsv"

CASE_DIR="$ROOT/run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh/mpi8origin"
CTRL="$CASE_DIR/system/controlDict"
BACKUP="$CTRLDIR/ourmesh_MPI8origin_controlDict_before_20260610"

mkdir -p "$LOGDIR" "$CTRLDIR"

set +u
source "$ROOT/doc/scripts/env.sh"
set -u

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

configure_nowrite_control()
{
    local file="$1"

    set_entry "$file" "startFrom" "startTime"
    set_entry "$file" "startTime" "0"
    set_entry "$file" "stopAt" "endTime"
    set_entry "$file" "endTime" "5.e-05"
    set_entry "$file" "deltaT" "1.e-07"
    set_entry "$file" "writeControl" "runTime"
    set_entry "$file" "writeInterval" "1.e-3"
    set_entry "$file" "profileSummary" "true"
    set_entry "$file" "profileDetail" "false"
    set_entry "$file" "runTimeModifiable" "no"
}

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

if [ ! -f "$MANIFEST" ]
then
    printf "case\tmode\trep\tnp\tomp_threads\tstart\tend\texit\tlog\tcontrolDict\n" > "$MANIFEST"
fi

cp "$CTRL" "$BACKUP"

for rep in 1 2 3
do
    log="$LOGDIR/ourmesh_MPI8origin_rep${rep}.log"
    exit_file="$LOGDIR/ourmesh_MPI8origin_rep${rep}.exit"
    ctrl_copy="$CTRLDIR/ourmesh_MPI8origin_rep${rep}_controlDict"

    cp "$BACKUP" "$CTRL"
    configure_nowrite_control "$CTRL"
    cp "$CTRL" "$ctrl_copy"

    start_ts="$(date -Iseconds)"
    printf "ourmesh\tMPI8origin\t%s\t8\t1\t%s\t" "$rep" "$start_ts" >> "$MANIFEST"

    (
        cd "$CASE_DIR"
        export OMP_NUM_THREADS=1
        echo "case=ourmesh"
        echo "mode=MPI8origin"
        echo "rep=$rep"
        echo "np=8"
        echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
        echo "case_dir=$CASE_DIR"
        echo "controlDict=$ctrl_copy"
        echo "start=$start_ts"
        /usr/bin/time -p /usr/bin/timeout 2400 mpirun -np 8 dsmcFoam+ -parallel
    ) > "$log" 2>&1
    status=$?

    end_ts="$(date -Iseconds)"
    printf "%s\n" "$status" > "$exit_file"
    printf "%s\t%s\t%s\t%s\n" "$end_ts" "$status" "$log" "$ctrl_copy" >> "$MANIFEST"
done

restore_control
trap - EXIT
