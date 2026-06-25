HYSTRATH_DLB=/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb

# Source environment
source ~/intel/oneapi/setvars.sh --force 2>/dev/null
source ~/code/OpenFoam/OF-1706/OpenFOAM-v1706/etc/bashrc

# Override paths to keep everything local
export WM_PROJECT_USER_DIR=$HYSTRATH_DLB
export FOAM_USER_LIBBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/lib
export FOAM_USER_APPBIN=$HYSTRATH_DLB/platforms/$WM_OPTIONS/bin
export FOAM_PROJECT_LIBBIN=$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib
export FOAM_PROJECT_APPBIN=$WM_PROJECT_DIR/platforms/$WM_OPTIONS/bin
export FOAM_PROJECT_MPI_LIBBIN=$FOAM_PROJECT_LIBBIN/$FOAM_MPI

# Prefer this branch's applications, then the matching base OpenFOAM tools.
# This avoids picking stale ~/OpenFOAM user binaries such as reconstructPar.
export PATH=$FOAM_USER_APPBIN:$FOAM_PROJECT_APPBIN:$PATH
export PARMETIS_DIR=/home/superxcx/code/DSMC/dsmcFoam++/parmetis-install
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$PARMETIS_DIR/lib:$FOAM_PROJECT_LIBBIN:$FOAM_PROJECT_MPI_LIBBIN:$LD_LIBRARY_PATH
