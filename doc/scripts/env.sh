HYSTRATH_DLB=/home/superxcx/code/OpenFoam/OF-1706/work/v1706-sparta

# Source environment
source ~/intel/oneapi/setvars.sh --force 2>/dev/null
source ~/code/OpenFoam/OF-1706/OpenFOAM-v1706/etc/bashrc

# Override paths to keep everything local
export WM_PROJECT_USER_DIR=$HYSTRATH_DLB
export FOAM_USER_LIBBIN=$HYSTRATH_DLB/platforms/linux64IccDPInt32Opt/lib
export FOAM_USER_APPBIN=$HYSTRATH_DLB/platforms/linux64IccDPInt32Opt/bin
export PATH=$FOAM_USER_APPBIN:$PATH
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH