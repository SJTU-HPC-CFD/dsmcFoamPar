

# HPC build script for hyStrath_xcx on OpenFOAM-v2506
# Preferred usage:
#   source build_dsmc_hpc.sh

module purge 
module load compiler/gcc/11.4.0 
module load compiler/intel/oneapi/2021.3.0

export OF2506_DIR="/work/home/acpovf0xmc/downloads/openfoam/v2506/OpenFOAM-v2506"
export DSMC_USER_DIR="/work/home/acpovf0xmc/users/xiao_chen_xiang/DSMC/dsmcFoamPar/base"
source "${OF2506_DIR}/etc/bashrc"
export WM_PROJECT_USER_DIR="${DSMC_USER_DIR}"
export FOAM_USER_APPBIN="${WM_PROJECT_USER_DIR}/platforms/${WM_OPTIONS}/bin"
export FOAM_USER_LIBBIN="${WM_PROJECT_USER_DIR}/platforms/${WM_OPTIONS}/lib"
export LD_LIBRARY_PATH="${FOAM_USER_LIBBIN}:${LD_LIBRARY_PATH:-}"

wmakeLnIncludeAll "${WM_PROJECT_USER_DIR}/src/lagrangian/basic" 
wmake -j libso "${WM_PROJECT_USER_DIR}/src/lagrangian/basic" 
wmakeLnIncludeAll "${WM_PROJECT_USER_DIR}/src/lagrangian/molecularDynamics/general" 
wmake -j libso "${WM_PROJECT_USER_DIR}/src/lagrangian/molecularDynamics/general" 
wmakeLnIncludeAll "${WM_PROJECT_USER_DIR}/src/lagrangian/dsmc" 
wmake -j libso "${WM_PROJECT_USER_DIR}/src/lagrangian/dsmc" 
wmake -j "${WM_PROJECT_USER_DIR}/applications/utilities/preProcessing/dsmc/dsmcInitialise+" 
wmake -j "${WM_PROJECT_USER_DIR}/applications/solvers/discreteMethods/dsmc/dsmcFoam+"

# run case:
# OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp8 2>&1

# 使用说明：
# controDict里修改openmpThreads，目前只用openmp并行，不用mpi并行