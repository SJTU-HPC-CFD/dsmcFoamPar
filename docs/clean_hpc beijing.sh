

# HPC build script for hyStrath_xcx on OpenFOAM-v2506
# Preferred usage:
#   source build_dsmc_hpc.sh
#   source build_dsmc_hpc.sh clean}"
source ~/users/xiao_chen_xiang/software/OpenFOAM/OF-v2506/env.sh
    export DSMC_USER_DIR="/publicfs01/fs1-m8/home/m8s000774/users/xiao_chen_xiang/611/dsmcFoamPar/v2506/dlb-openmp"

    export WM_PROJECT_USER_DIR="${DSMC_USER_DIR}"
    export FOAM_USER_APPBIN="${WM_PROJECT_USER_DIR}/platforms/${WM_OPTIONS}/bin"
    export FOAM_USER_LIBBIN="${WM_PROJECT_USER_DIR}/platforms/${WM_OPTIONS}/lib"
    export LD_LIBRARY_PATH="${FOAM_USER_LIBBIN}:${LD_LIBRARY_PATH:-}"

    wclean libso "${WM_PROJECT_USER_DIR}/src/lagrangian/basic" 
    wclean libso "${WM_PROJECT_USER_DIR}/src/lagrangian/molecularDynamics/general" 
    wclean libso "${WM_PROJECT_USER_DIR}/src/lagrangian/dsmc" 
    wclean "${WM_PROJECT_USER_DIR}/applications/solvers/discreteMethods/dsmc/dsmcFoam+" 
    wclean "${WM_PROJECT_USER_DIR}/applications/utilities/preProcessing/dsmc/dsmcInitialise+"
    find "${WM_PROJECT_USER_DIR}" -name '*.dep' -delete



