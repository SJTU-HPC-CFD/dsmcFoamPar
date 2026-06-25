#!/bin/bash
# Build dsmcFoam+ from hyStrath_dlb into local directory
# Usage: source build-dsmcFoam.sh

set -e

HYSTRATH_DLB=/publicfs01/fs1-m8/home/m8s000774/users/xiao_chen_xiang/611/dsmcFoamPar/v1706/base/dsmcFoamPar-1706-dlb

# Source environment
source /publicfs01/fs1-m8/home/m8s000774/users/xiao_chen_xiang/package/mbbq/OpenFOAM-v1706/etc/rebashrc

# Override paths to keep everything local
export WM_PROJECT_USER_DIR=$HYSTRATH_DLB
export FOAM_USER_LIBBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/lib
export FOAM_USER_APPBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/bin
export FOAM_PROJECT_LIBBIN=$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib
export FOAM_PROJECT_APPBIN=$WM_PROJECT_DIR/platforms/$WM_OPTIONS/bin
export FOAM_PROJECT_MPI_LIBBIN=$FOAM_PROJECT_LIBBIN/$FOAM_MPI
export PATH=$FOAM_USER_APPBIN:$FOAM_PROJECT_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$FOAM_PROJECT_LIBBIN:$FOAM_PROJECT_MPI_LIBBIN:$LD_LIBRARY_PATH

echo "=== Build targets ==="
echo "  LIBBIN: $FOAM_USER_LIBBIN"
echo "  APPBIN: $FOAM_USER_APPBIN"
echo "  USER_DIR: $WM_PROJECT_USER_DIR"
echo "  BASEBIN: $FOAM_PROJECT_APPBIN"
echo ""

# 1. Build liblagrangian+ (basic lagrangian library)
echo "=== Building liblagrangian+ ==="
cd $HYSTRATH_DLB/src/lagrangian/basic
wmake -j lnInclude 2>/dev/null || true
wmake -j libso

# 2. Build libdecompose
echo "=== Building libdecompose ==="
cd $HYSTRATH_DLB/src/parallel/decompose/decompose
wmake -j lnInclude 2>/dev/null || true
wmake -j libso

# 3. Build libgeneralMolecule
echo "=== Building libgeneralMolecule ==="
cd $HYSTRATH_DLB/src/lagrangian/molecularDynamics/general
wmake -j lnInclude 2>/dev/null || true
wmake -j libso

# 4. Build libdsmcFoam+
echo "=== Building libdsmcFoam+ ==="
cd $HYSTRATH_DLB/src/lagrangian/dsmc
wmake -j lnInclude 2>/dev/null || true
wmake -j libso

# 5. Build dsmcFoam+ solver
echo "=== Building dsmcFoam+ ==="
cd $HYSTRATH_DLB/applications/solvers/discreteMethods/dsmc/dsmcFoam+
wmake -j

# 6. Build dsmcInitialise+ utility
echo "=== Building dsmcInitialise+ ==="
cd $HYSTRATH_DLB/applications/utilities/preProcessing/dsmc/dsmcInitialise+
wmake -j

echo ""
echo "=== Done ==="
echo "dsmcFoam+ binary: $FOAM_USER_APPBIN/dsmcFoam+"
echo "dsmcInitialise+ binary: $FOAM_USER_APPBIN/dsmcInitialise+"
ls -la $FOAM_USER_APPBIN/dsmcFoam+ 2>/dev/null
ls -la $FOAM_USER_APPBIN/dsmcInitialise+ 2>/dev/null
cd $HYSTRATH_DLB
