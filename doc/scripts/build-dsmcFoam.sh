#!/bin/bash
# Build dsmcFoam+ from hyStrath_dlb into local directory
# Usage: source build-dsmcFoam.sh

set -e

HYSTRATH_DLB=/home/superxcx/code/OpenFoam/OF-1706/work/v1706-sparta

# Source environment
source ~/intel/oneapi/setvars.sh --force 2>/dev/null
source ~/code/OpenFoam/OF-1706/OpenFOAM-v1706/etc/bashrc

# Override paths to keep everything local
export WM_PROJECT_USER_DIR=$HYSTRATH_DLB
export FOAM_USER_LIBBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/lib
export FOAM_USER_APPBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/bin
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH

echo "=== Build targets ==="
echo "  LIBBIN: $FOAM_USER_LIBBIN"
echo "  APPBIN: $FOAM_USER_APPBIN"
echo "  USER_DIR: $WM_PROJECT_USER_DIR"
echo ""

# 1. Build liblagrangian+ (basic lagrangian library)
echo "=== Building liblagrangian+ ==="
cd $HYSTRATH_DLB/src/lagrangian/basic
wmake -j lnInclude 2>/dev/null || true
wmake -j libso

# 2. Build libgeneralMolecule
echo "=== Building libgeneralMolecule ==="
cd $HYSTRATH_DLB/src/lagrangian/molecularDynamics/general
wmake -j lnInclude 2>/dev/null || true
wmake -j libso

# 3. Build libdsmcFoam+
echo "=== Building libdsmcFoam+ ==="
cd $HYSTRATH_DLB/src/lagrangian/dsmc
wmake -j lnInclude 2>/dev/null || true
wmake -j libso

# 4. Build dsmcFoam+ solver
echo "=== Building dsmcFoam+ ==="
cd $HYSTRATH_DLB/applications/solvers/discreteMethods/dsmc/dsmcFoam+
wmake -j

# 5. Build dsmcInitialise+ utility
echo "=== Building dsmcInitialise+ ==="
cd $HYSTRATH_DLB/applications/utilities/preProcessing/dsmc/dsmcInitialise+
wmake -j

echo ""
echo "=== Done ==="
echo "dsmcFoam+ binary: $FOAM_USER_APPBIN/dsmcFoam+"
echo "dsmcInitialise+ binary: $FOAM_USER_APPBIN/dsmcInitialise+"

ls -la $FOAM_USER_APPBIN/dsmcFoam+ 2>/dev/null
ls -la $FOAM_USER_APPBIN/dsmcInitialise+ 2>/dev/null
