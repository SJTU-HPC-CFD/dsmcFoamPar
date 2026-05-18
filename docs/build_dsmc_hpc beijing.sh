

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
wmakeLnIncludeAll "${WM_PROJECT_USER_DIR}/src/lagrangian/basic" 
wmake -j32 libso "${WM_PROJECT_USER_DIR}/src/lagrangian/basic" 
wmakeLnIncludeAll "${WM_PROJECT_USER_DIR}/src/lagrangian/molecularDynamics/general" 
wmake -j32 libso "${WM_PROJECT_USER_DIR}/src/lagrangian/molecularDynamics/general" 
wmakeLnIncludeAll "${WM_PROJECT_USER_DIR}/src/lagrangian/dsmc" 
wmake -j32 libso "${WM_PROJECT_USER_DIR}/src/lagrangian/dsmc" 
wmake -j32 "${WM_PROJECT_USER_DIR}/applications/utilities/preProcessing/dsmc/dsmcInitialise+" 
wmake -j32 "${WM_PROJECT_USER_DIR}/applications/solvers/discreteMethods/dsmc/dsmcFoam+"

source /public1/soft/modules/module.sh
module purge
module load gcc/11.2-para mpi/oneAPI/2022.1
source /public4/home/a0s001239/users/xiao_chen_xiang/software/OF-v2506/OpenFOAM-v2506/etc/bashrc

source /public4/home/a0s001239/users/xiao_chen_xiang/software/OF-v2506/env.sh

