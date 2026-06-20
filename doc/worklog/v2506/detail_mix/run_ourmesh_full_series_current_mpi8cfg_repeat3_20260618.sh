#!/usr/bin/env bash

set -u

ROOT="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb"
DETAIL="$ROOT/doc/worklog/v2506/detail_mix"
OUTDIR="${OUTDIR:-$DETAIL/ourmesh_full_series_current_mpi8cfg_repeat3_20260618}"
LOGDIR="$OUTDIR/logs"
CTRLDIR="$OUTDIR/controlDicts"
MANIFEST="$OUTDIR/run_manifest.tsv"

OUR_ROOT="$ROOT/run/hyStrath/dsmcFoam+/xcx_test/pal-phd-3.3.1-cylinder-react/timecompare/ourmesh"

OMP8_CASE="$OUR_ROOT/omp8"
MIX_CASE="$OUR_ROOT/mix-mpi4omp2"
MPI8_CASE="$OUR_ROOT/mpi8replicatedmesh"
MPI8ORIGIN_CASE="$OUR_ROOT/mpi8origin"

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

configure_common_nowrite()
{
    local file="$1"
    local use_openmp="$2"
    local omp_threads="$3"

    set_entry "$file" "nTerminalOutputs" "10"
    set_entry "$file" "startFrom" "startTime"
    set_entry "$file" "startTime" "0"
    set_entry "$file" "stopAt" "endTime"
    set_entry "$file" "endTime" "5.e-05"
    set_entry "$file" "deltaT" "1.e-07"
    set_entry "$file" "writeControl" "runTime"
    set_entry "$file" "writeInterval" "1.e-3"
    set_entry "$file" "runTimeModifiable" "no"
    set_entry "$file" "profileSummary" "true"
    set_entry "$file" "profileDetail" "false"
    set_entry "$file" "collisionFastRng" "true"
    set_entry "$file" "useOpenMP" "$use_openmp"
    set_entry "$file" "openmpThreads" "$omp_threads"
    set_entry "$file" "openmpMoveSchedule" "static"
    set_entry "$file" "openmpMoveChunk" "64"
    set_entry "$file" "openmpCollisionSchedule" "dynamic"
    set_entry "$file" "openmpCollisionChunk" "8"
}

configure_replicated_formal()
{
    local file="$1"
    local omp_threads="$2"

    configure_common_nowrite "$file" "true" "$omp_threads"
    if [ "$omp_threads" = "1" ]
    then
        set_entry "$file" "useOpenMP" "false"
    fi

    set_entry "$file" "replicatedMesh" "true"
    set_entry "$file" "replicatedMeshDelayedReceive" "true"
    set_entry "$file" "replicatedMeshNoAlltoall" "true"
    set_entry "$file" "replicatedMeshFlatTransfer" "true"
    set_entry "$file" "replicatedMeshDecompMethod" "metis"
    set_entry "$file" "replicatedMeshMigrateInterval" "10"
    set_entry "$file" "replicatedMeshDLBDualConstraint" "false"
    set_entry "$file" "replicatedMeshAutoDLB" "true"
    set_entry "$file" "replicatedMeshDLBSteps" "50"
    set_entry "$file" "replicatedMeshDLBForceSteps" "( )"
    set_entry "$file" "replicatedMeshDLBTriggerMode" "legacyWindow"
    set_entry "$file" "replicatedMeshDLBMinGapSteps" "50"
    set_entry "$file" "replicatedMeshDLBCheckCollective" "allgather"
    set_entry "$file" "replicatedMeshDLBSkipPostDiag" "false"
    set_entry "$file" "replicatedMeshDLBImbalanceThreshold" "1.5"
    set_entry "$file" "replicatedMeshDLBItr" "1000"
    set_entry "$file" "replicatedMeshDLBUbvec" "1.05"
    set_entry "$file" "replicatedMeshDLBProfile" "false"
    set_entry "$file" "replicatedMeshDLBProfileSteps" "5"
    set_entry "$file" "replicatedMeshDLBAlpha" "0.8"
    set_entry "$file" "replicatedMeshDLBVsizeExp" "0"
}

configure_origin_nowrite()
{
    local file="$1"
    configure_common_nowrite "$file" "false" "1"
}

declare -a BACKUPS=(
    "omp8|$OMP8_CASE/system/controlDict"
    "mix|$MIX_CASE/system/controlDict"
    "mpi8|$MPI8_CASE/system/controlDict"
    "mpi8origin|$MPI8ORIGIN_CASE/system/controlDict"
)

backup_controls()
{
    local entry label file
    for entry in "${BACKUPS[@]}"
    do
        IFS='|' read -r label file <<< "$entry"
        cp "$file" "$CTRLDIR/${label}_controlDict_before_20260618"
    done
}

restore_controls()
{
    local entry label file backup
    for entry in "${BACKUPS[@]}"
    do
        IFS='|' read -r label file <<< "$entry"
        backup="$CTRLDIR/${label}_controlDict_before_20260618"
        if [ -f "$backup" ]
        then
            cp "$backup" "$file"
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
    local mode="$1"
    local rep="$2"
    local case_dir="$3"
    local np="$4"
    local omp_threads="$5"
    local timeout_s="$6"
    local ctrl="$case_dir/system/controlDict"
    local ctrl_copy="$CTRLDIR/ourmesh_${mode}_rep${rep}_controlDict"
    local log="$LOGDIR/ourmesh_${mode}_rep${rep}.log"
    local exit_file="$LOGDIR/ourmesh_${mode}_rep${rep}.exit"
    local start_ts end_ts status

    case "$mode" in
        OMP8)
            configure_common_nowrite "$ctrl" "true" "8"
            set_entry "$ctrl" "replicatedMesh" "false"
            set_entry "$ctrl" "replicatedMeshAutoDLB" "false"
            ;;
        MPI2xOMP4)
            configure_replicated_formal "$ctrl" "4"
            ;;
        MPI4xOMP2)
            configure_replicated_formal "$ctrl" "2"
            ;;
        MPI8)
            configure_replicated_formal "$ctrl" "1"
            ;;
        MPI8origin)
            configure_origin_nowrite "$ctrl"
            ;;
        *)
            echo "unknown mode: $mode" >&2
            return 2
            ;;
    esac

    cp "$ctrl" "$ctrl_copy"

    start_ts="$(date -Iseconds)"
    printf "ourmesh\t%s\t%s\t%s\t%s\t%s\t" \
        "$mode" "$rep" "$np" "$omp_threads" "$start_ts" >> "$MANIFEST"

    (
        cd "$case_dir"
        export OMP_NUM_THREADS="$omp_threads"
        echo "case=ourmesh"
        echo "mode=$mode"
        echo "rep=$rep"
        echo "np=$np"
        echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
        echo "case_dir=$case_dir"
        echo "controlDict=$ctrl_copy"
        echo "start=$start_ts"
        if [ "$mode" = "MPI8origin" ]
        then
            /usr/bin/time -p /usr/bin/timeout "$timeout_s" mpirun -np "$np" dsmcFoam+ -parallel
        elif [ "$np" -eq 1 ]
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
    run_one "OMP8" "$rep" "$OMP8_CASE" 1 8 900
    run_one "MPI2xOMP4" "$rep" "$MIX_CASE" 2 4 1200
    run_one "MPI4xOMP2" "$rep" "$MIX_CASE" 4 2 1200
    run_one "MPI8" "$rep" "$MPI8_CASE" 8 1 1500
    run_one "MPI8origin" "$rep" "$MPI8ORIGIN_CASE" 8 1 1800
done

restore_controls
trap - EXIT
