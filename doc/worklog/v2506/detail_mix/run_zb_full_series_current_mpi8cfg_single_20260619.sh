#!/usr/bin/env bash

set -u

ROOT="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb"
DETAIL="$ROOT/doc/worklog/v2506/detail_mix"
OUTDIR="${OUTDIR:-$DETAIL/zb_full_series_current_mpi8cfg_single_20260619}"
LOGDIR="$OUTDIR/logs"
CTRLDIR="$OUTDIR/controlDicts"
MANIFEST="$OUTDIR/run_manifest.tsv"

ZB_ROOT="$ROOT/run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react"

OMP8_CASE="$ZB_ROOT/omp8"
MPI2OMP4_CASE="$ZB_ROOT/mpi2omp4"
MPI4OMP2_CASE="$ZB_ROOT/mpi4omp2"
MPI8_CASE="$ZB_ROOT/mpi8"
MPI8ORIGIN_CASE="$ZB_ROOT/mpi8origin"

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
    set_entry "$file" "endTime" "1.9920146682e-05"
    set_entry "$file" "deltaT" "6.640048894e-08"
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

configure_replicated_best()
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
    set_entry "$file" "replicatedMeshNoAlltoall" "false"
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
    set_entry "$file" "replicatedMeshDLBAlpha" "1"
    set_entry "$file" "replicatedMeshDLBVsizeExp" "0"
    set_entry "$file" "replicatedMeshDLBParticleGate" "false"
    set_entry "$file" "replicatedMeshDLBAdaptiveAlpha" "false"
    set_entry "$file" "replicatedMeshSARSteps" "50"
    set_entry "$file" "replicatedMeshDLBMinRemainingSteps" "50"
    set_entry "$file" "replicatedMeshGatherCandidates" "false"
    set_entry "$file" "replicatedMeshOverlapSizeExchange" "false"
}

configure_origin_nowrite()
{
    local file="$1"
    configure_common_nowrite "$file" "false" "1"
}

declare -a BACKUPS=(
    "omp8|$OMP8_CASE/system/controlDict"
    "mpi2omp4|$MPI2OMP4_CASE/system/controlDict"
    "mpi4omp2|$MPI4OMP2_CASE/system/controlDict"
    "mpi8|$MPI8_CASE/system/controlDict"
    "mpi8origin|$MPI8ORIGIN_CASE/system/controlDict"
)

backup_controls()
{
    local entry label file
    for entry in "${BACKUPS[@]}"
    do
        IFS='|' read -r label file <<< "$entry"
        cp "$file" "$CTRLDIR/${label}_controlDict_before_20260619"
    done
}

restore_controls()
{
    local entry label file backup
    for entry in "${BACKUPS[@]}"
    do
        IFS='|' read -r label file <<< "$entry"
        backup="$CTRLDIR/${label}_controlDict_before_20260619"
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

prepare_origin_case()
{
    local decomp="$MPI8ORIGIN_CASE/system/decomposeParDict"
    set_entry "$decomp" "numberOfSubdomains" "8"
    if [ ! -d "$MPI8ORIGIN_CASE/processor0" ]
    then
        (
            cd "$MPI8ORIGIN_CASE"
            decomposePar
        ) > "$LOGDIR/zb_MPI8origin_decomposePar_20260619.log" 2>&1
    fi
}

run_one()
{
    local mode="$1"
    local case_dir="$2"
    local np="$3"
    local omp_threads="$4"
    local timeout_s="$5"
    local ctrl="$case_dir/system/controlDict"
    local ctrl_copy="$CTRLDIR/zb_${mode}_controlDict"
    local log="$LOGDIR/zb_${mode}.log"
    local exit_file="$LOGDIR/zb_${mode}.exit"
    local start_ts end_ts status

    case "$mode" in
        OMP8)
            configure_common_nowrite "$ctrl" "true" "8"
            ;;
        MPI2xOMP4)
            configure_replicated_best "$ctrl" "4"
            ;;
        MPI4xOMP2)
            configure_replicated_best "$ctrl" "2"
            ;;
        MPI8)
            configure_replicated_best "$ctrl" "1"
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
    printf "zb\t%s\t%s\t%s\t%s\t" \
        "$mode" "$np" "$omp_threads" "$start_ts" >> "$MANIFEST"

    (
        cd "$case_dir"
        export OMP_NUM_THREADS="$omp_threads"
        echo "case=zb"
        echo "mode=$mode"
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
prepare_origin_case

printf "case\tmode\tnp\tomp_threads\tstart\tend\texit\tlog\tcontrolDict\n" > "$MANIFEST"

run_one "OMP8" "$OMP8_CASE" 1 8 900
run_one "MPI2xOMP4" "$MPI2OMP4_CASE" 2 4 1200
run_one "MPI4xOMP2" "$MPI4OMP2_CASE" 4 2 1200
run_one "MPI8" "$MPI8_CASE" 8 1 1200
run_one "MPI8origin" "$MPI8ORIGIN_CASE" 8 1 1800

restore_controls
trap - EXIT
