source /home/superxcx/intel/oneapi/setvars.sh --force 2>/dev/null
source /home/superxcx/code/OpenFoam/OF-2506/OpenFOAM-v2506/etc/bashrc

export WM_PROJECT_USER_DIR=/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:/home/superxcx/code/DSMC/dsmcFoam++/parmetis-install/lib:$LD_LIBRARY_PATH

cd $WM_PROJECT_USER_DIR

echo "=== Building lagrangian/basic ==="
wmake -j src/lagrangian/basic

echo "=== Building lagrangian/molecularDynamics/general ==="
wmake -j src/lagrangian/molecularDynamics/general

echo "=== Building lagrangian/dsmc ==="
wmake -j src/lagrangian/dsmc

echo "=== Building dsmcFoam+ solver ==="
wmake -j applications/solvers/discreteMethods/dsmc/dsmcFoam+

echo "=== Building dsmcInitialise+ ==="
wmake -j applications/utilities/preProcessing/dsmc/dsmcInitialise+

echo "=== Build complete ==="
echo "Library: $FOAM_USER_LIBBIN/libdsmcFoam+.so"
echo "Binary:  $FOAM_USER_APPBIN/dsmcFoam+"
echo "Binary:  $FOAM_USER_APPBIN/dsmcInitialise+"
