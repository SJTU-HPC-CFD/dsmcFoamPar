source ~/intel/oneapi/setvars.sh >/dev/null 2>&1
source /home/superxcx/code/OpenFoam/OF-2506/OpenFOAM-v2506/etc/bashrc >/dev/null 2>&1
export WM_PROJECT_USER_DIR=/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export LD_LIBRARY_PATH=${FOAM_USER_LIBBIN}:$LD_LIBRARY_PATH
wmake libso /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/src/lagrangian/dsmc