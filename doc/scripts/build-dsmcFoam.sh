#!/bin/bash
# Build dsmcFoam+ from hyStrath_dlb into local directory
# Usage: source build-dsmcFoam.sh

set -e

HYSTRATH_DLB=/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb

# Source environment
source ~/intel/oneapi/setvars.sh --force 2>/dev/null
source ~/code/OpenFoam/OF-1706/OpenFOAM-v1706/etc/bashrc

# Override paths to keep everything local
export WM_PROJECT_USER_DIR=$HYSTRATH_DLB
export FOAM_USER_LIBBIN=$HYSTRATH_DLB/platforms/linux64IccDPInt32Opt/lib
export FOAM_USER_APPBIN=$HYSTRATH_DLB/platforms/linux64IccDPInt32Opt/bin

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

echo ""
echo "=== Done ==="
echo "dsmcFoam+ binary: $FOAM_USER_APPBIN/dsmcFoam+"
ls -la $FOAM_USER_APPBIN/dsmcFoam+ 2>/dev/null
