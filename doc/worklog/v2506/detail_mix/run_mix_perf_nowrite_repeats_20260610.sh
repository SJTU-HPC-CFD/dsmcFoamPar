#!/usr/bin/env bash

set -u

ROOT="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb"
DETAIL="$ROOT/doc/worklog/v2506/detail_mix"
OUTDIR="$DETAIL/repeat_nowrite_compute_20260610_rerun"
LOGDIR="$OUTDIR/logs"
CTRLDIR="$OUTDIR/controlDicts"
MANIFEST="$OUTDIR/run_manifest.tsv"

OURMESH_ROOT="$ROOT/run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh"
ZB_ROOT="$ROOT/run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react"

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
    local case_name="$2"

    set_entry "$file" "startFrom" "startTime"
    set_entry "$file" "startTime" "0"
    set_entry "$file" "stopAt" "endTime"

    if [ "$case_name" = "ourmesh" ]
    then
        set_entry "$file" "endTime" "5.e-05"
        set_entry "$file" "deltaT" "1.e-07"
    else
        set_entry "$file" "endTime" "1.9920146682e-05"
        set_entry "$file" "deltaT" "6.640048894e-08"
    fi

    set_entry "$file" "writeControl" "runTime"
    set_entry "$file" "writeInterval" "1.e-3"
    set_entry "$file" "profileSummary" "true"
    set_entry "$file" "profileDetail" "false"
    set_entry "$file" "runTimeModifiable" "no"
}

declare -a BACKUP_CASES=(
    "ourmesh_omp8|$OURMESH_ROOT/omp8"
    "ourmesh_mix|$OURMESH_ROOT/mix-mpi4omp2"
    "ourmesh_mpi8|$OURMESH_ROOT/mpi8replicatedmesh"
    "zb_omp8|$ZB_ROOT/omp8"
    "zb_mpi2omp4|$ZB_ROOT/mpi2omp4"
    "zb_mpi4omp2|$ZB_ROOT/mpi4omp2"
    "zb_mpi8|$ZB_ROOT/mpi8"
)

backup_controls()
{
    local entry label dir
    for entry in "${BACKUP_CASES[@]}"
    do
        IFS='|' read -r label dir <<< "$entry"
        cp "$dir/system/controlDict" "$CTRLDIR/${label}_controlDict_before_20260610"
    done
}

restore_controls()
{
    local entry label dir backup
    for entry in "${BACKUP_CASES[@]}"
    do
        IFS='|' read -r label dir <<< "$entry"
        backup="$CTRLDIR/${label}_controlDict_before_20260610"
        if [ -f "$backup" ]
        then
            cp "$backup" "$dir/system/controlDict"
        fi
    done
}

on_exit()
{
    local status=$?
    restore_controls
    exit "$status"
}

trap on_exit EXIT

run_one()
{
    local case_name="$1"
    local mode="$2"
    local rep="$3"
    local case_dir="$4"
    local template="$5"
    local np="$6"
    local omp_threads="$7"
    local timeout_s="$8"

    local ctrl="$case_dir/system/controlDict"
    local log="$LOGDIR/${case_name}_${mode}_rep${rep}.log"
    local exit_file="$LOGDIR/${case_name}_${mode}_rep${rep}.exit"
    local ctrl_copy="$CTRLDIR/${case_name}_${mode}_rep${rep}_controlDict"
    local start_ts end_ts status

    cp "$template" "$ctrl"
    configure_nowrite_control "$ctrl" "$case_name"
    cp "$ctrl" "$ctrl_copy"

    start_ts="$(date -Iseconds)"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t" \
        "$case_name" "$mode" "$rep" "$np" "$omp_threads" "$start_ts" >> "$MANIFEST"

    (
        cd "$case_dir"
        export OMP_NUM_THREADS="$omp_threads"
        echo "case=$case_name"
        echo "mode=$mode"
        echo "rep=$rep"
        echo "np=$np"
        echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
        echo "case_dir=$case_dir"
        echo "controlDict=$ctrl_copy"
        echo "start=$start_ts"
        if [ "$np" -eq 1 ]
        then
            /usr/bin/time -p /usr/bin/timeout "$timeout_s" dsmcFoam+
        else
            /usr/bin/time -p /usr/bin/timeout "$timeout_s" mpirun -np "$np" dsmcFoam+
        fi
    ) > "$log" 2>&1
    status=$?

    end_ts="$(date -Iseconds)"
    printf "%s\n" "$status" > "$exit_file"
    printf "%s\t%s\t%s\t%s\n" "$end_ts" "$status" "$log" "$ctrl_copy" >> "$MANIFEST"
}

backup_controls

printf "case\tmode\trep\tnp\tomp_threads\tstart\tend\texit\tlog\tcontrolDict\n" > "$MANIFEST"

for rep in 1 2 3
do
    run_one "ourmesh" "OMP8" "$rep" \
        "$OURMESH_ROOT/omp8" \
        "$DETAIL/omp8_controlDict_compare_20260608_111426" \
        1 8 900

    run_one "ourmesh" "MPI2xOMP4" "$rep" \
        "$OURMESH_ROOT/mix-mpi4omp2" \
        "$DETAIL/mpi2omp4_controlDict_compare_20260608_111426" \
        2 4 1200

    run_one "ourmesh" "MPI4xOMP2" "$rep" \
        "$OURMESH_ROOT/mix-mpi4omp2" \
        "$DETAIL/mpi4omp2_controlDict_compare_20260608_111426" \
        4 2 1200

    run_one "ourmesh" "MPI8" "$rep" \
        "$OURMESH_ROOT/mpi8replicatedmesh" \
        "$DETAIL/mpi8replicatedmesh_controlDict_compare_20260608_111426" \
        8 1 1500

    run_one "zb" "OMP8" "$rep" \
        "$ZB_ROOT/omp8" \
        "$DETAIL/zb_omp8_controlDict_compare_rawmpi_guard_nowrite_20260608_115212" \
        1 8 900

    run_one "zb" "MPI2xOMP4" "$rep" \
        "$ZB_ROOT/mpi2omp4" \
        "$DETAIL/zb_mpi2omp4_controlDict_compare_rawmpi_guard_nowrite_20260608_115212" \
        2 4 1200

    run_one "zb" "MPI4xOMP2" "$rep" \
        "$ZB_ROOT/mpi4omp2" \
        "$DETAIL/zb_mpi4omp2_controlDict_compare_rawmpi_guard_nowrite_20260608_115212" \
        4 2 1200

    run_one "zb" "MPI8" "$rep" \
        "$ZB_ROOT/mpi8" \
        "$DETAIL/zb_mpi8_controlDict_compare_rawmpi_guard_nowrite_20260608_115212" \
        8 1 1500
done

restore_controls
trap - EXIT
