#!/bin/bash
# Run heatBath-5species test case with hyStrath_dlb dsmcFoam+
# Usage: source run-heatBath-5species.sh [nProcs]
# Example: source run-heatBath-5species.sh 4

set -e

HYSTRATH_DLB=/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb
CASE_DIR=$HYSTRATH_DLB/run/hyStrath/dsmcFoam+/heatBath-5species
nProcs=${1:-1}

# Source environment
source ~/intel/oneapi/setvars.sh --force 2>/dev/null
source ~/code/OpenFoam/OF-1706/OpenFOAM-v1706/etc/bashrc

# Override paths to use hyStrath_dlb binaries/libs
export WM_PROJECT_USER_DIR=$HYSTRATH_DLB
export FOAM_USER_LIBBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/lib
export FOAM_USER_APPBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/bin
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH

echo "=== Environment ==="
echo "  APPBIN: $FOAM_USER_APPBIN"
echo "  LIBBIN: $FOAM_USER_LIBBIN"
echo "  CASE:   $CASE_DIR"
echo "  nProcs: $nProcs"
echo ""

# Build dsmcInitialise+ if not present
if [ ! -f "$FOAM_USER_APPBIN/dsmcInitialise+" ]; then
    echo "=== Building dsmcInitialise+ ==="
    cd $HYSTRATH_DLB/applications/utilities/preProcessing/dsmc/dsmcInitialise+
    wmake -j
    echo ""
fi

# Clean case
echo "=== Cleaning case ==="
cd $CASE_DIR
rm -rf 0 boundaries fieldMeasurements processor* log.* > /dev/null 2>&1
rm -rf constant/polyMesh > /dev/null 2>&1

# Run
echo "=== Running blockMesh ==="
blockMesh > log.blockMesh 2>&1
echo "  Done."

echo "=== Running dsmcInitialise+ ==="
dsmcInitialise+ > log.dsmcInitialise+ 2>&1
echo "  Done."

if [ $nProcs -eq 1 ]; then
    echo "=== Running dsmcFoam+ (serial) ==="
    dsmcFoam+ > log.dsmcFoam+ 2>&1
    echo "  Done."
else
    # Update decomposeParDict
    sed -i "s/numberOfSubdomains.*/numberOfSubdomains $nProcs;/" system/decomposeParDict
    echo "=== Running decomposePar ==="
    decomposePar -latestTime > log.decomposePar 2>&1
    echo "  Done."
    echo "=== Running dsmcFoam+ (parallel, $nProcs procs) ==="
    mpirun -np $nProcs dsmcFoam+ -parallel > log.dsmcFoam+ 2>&1
    echo "  Done."
    echo "=== Running reconstructPar ==="
    reconstructPar -latestTime > log.reconstructPar 2>&1
    echo "  Done."
fi

echo ""
echo "=== Finished ==="
echo "Logs in: $CASE_DIR/log.*"
echo "Check Tvib: grep 'Tvib' $CASE_DIR/log.dsmcFoam+ | tail -5"
