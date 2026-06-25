#!/usr/bin/env bash

set -u

ROOT="/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb"
OUTDIR="$ROOT/doc/worklog/v2506/detail_mix/zb_4mode_300step_nowrite_migrate1_retest_20260622"
CASEDIR="$OUTDIR/cases"
LOGDIR="$OUTDIR/logs"
CTRLDIR="$OUTDIR/controlDicts"
MANIFEST="$OUTDIR/run_manifest.tsv"

ZB_ROOT="$ROOT/run/hyStrath/dsmcFoam+/xcx_test/zb-cylinder-react"

mkdir -p "$CASEDIR" "$LOGDIR" "$CTRLDIR"

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

prepare_case_copy()
{
    local src="$1"
    local dst="$2"

    rm -rf "$dst"
    mkdir -p "$dst"

    cp -a "$src/0" "$src/constant" "$src/system" "$src/boundaries" "$src/fieldMeasurements" "$dst"/
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

configure_replicated_current()
{
    local file="$1"
    local omp_threads="$2"

    configure_common_nowrite "$file" "true" "$omp_threads"
    if [ "$omp_threads" = "1" ]
    then
        set_entry "$file" "useOpenMP" "false"
    fi

    set_entry "$file" "replicatedMesh" "true"
    set_entry "$file" "replicatedMeshWriteMode" "processor"
    set_entry "$file" "replicatedMeshProcessorWriteTimeMesh" "false"
    set_entry "$file" "replicatedMeshDelayedReceive" "true"
    set_entry "$file" "replicatedMeshNoAlltoall" "false"
    set_entry "$file" "replicatedMeshFlatTransfer" "true"
    set_entry "$file" "replicatedMeshDecompMethod" "metis"
    set_entry "$file" "replicatedMeshMigrateInterval" "1"
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

configure_decompose()
{
    local file="$1"
    local n="$2"

    set_entry "$file" "numberOfSubdomains" "$n"
    set_entry "$file" "method" "scotch"
}

run_one()
{
    local mode="$1"
    local src_case="$2"
    local np="$3"
    local omp_threads="$4"
    local timeout_s="$5"

    local case_dir="$CASEDIR/$mode"
    local ctrl="$case_dir/system/controlDict"
    local decomp="$case_dir/system/decomposeParDict"
    local ctrl_copy="$CTRLDIR/zb_${mode}_controlDict"
    local decomp_copy="$CTRLDIR/zb_${mode}_decomposeParDict"
    local log="$LOGDIR/zb_${mode}.log"
    local exit_file="$LOGDIR/zb_${mode}.exit"
    local start_ts end_ts status

    prepare_case_copy "$src_case" "$case_dir"

    case "$mode" in
        OMP8)
            configure_common_nowrite "$ctrl" "true" "8"
            configure_decompose "$decomp" "1"
            ;;
        MPI2xOMP4)
            configure_replicated_current "$ctrl" "4"
            configure_decompose "$decomp" "2"
            ;;
        MPI4xOMP2)
            configure_replicated_current "$ctrl" "2"
            configure_decompose "$decomp" "4"
            ;;
        MPI8)
            configure_replicated_current "$ctrl" "1"
            configure_decompose "$decomp" "8"
            ;;
        MPI8origin)
            configure_origin_nowrite "$ctrl"
            configure_decompose "$decomp" "8"
            (
                cd "$case_dir"
                decomposePar
            ) > "$LOGDIR/zb_MPI8origin_decomposePar.log" 2>&1
            ;;
        *)
            echo "unknown mode: $mode" >&2
            return 2
            ;;
    esac

    cp "$ctrl" "$ctrl_copy"
    cp "$decomp" "$decomp_copy"

    start_ts="$(date -Iseconds)"
    printf "zb\t%s\t%s\t%s\t%s\t" \
        "$mode" "$np" "$omp_threads" "$start_ts" >> "$MANIFEST"

    (
        cd "$case_dir"
        export OMP_NUM_THREADS="$omp_threads"
        export OMP_PROC_BIND=close
        export OMP_PLACES=cores
        export I_MPI_PIN=on
        export I_MPI_PIN_DOMAIN="$omp_threads"
        export I_MPI_PIN_ORDER=compact
        echo "case=zb"
        echo "mode=$mode"
        echo "np=$np"
        echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"
        echo "case_dir=$case_dir"
        echo "controlDict=$ctrl_copy"
        echo "decomposeParDict=$decomp_copy"
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

printf "case\tmode\tnp\tomp_threads\tstart\tend\texit\tlog\tcontrolDict\n" > "$MANIFEST"

run_one "OMP8" "$ZB_ROOT/omp8" 1 8 1200
run_one "MPI2xOMP4" "$ZB_ROOT/mpi2omp4" 2 4 1800
run_one "MPI4xOMP2" "$ZB_ROOT/mpi4omp2" 4 2 1800
run_one "MPI8" "$ZB_ROOT/mpi8" 8 1 1800

