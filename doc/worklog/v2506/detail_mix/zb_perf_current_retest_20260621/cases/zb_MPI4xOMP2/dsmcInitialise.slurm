#!/bin/bash

#SBATCH --job-name=DSMC-test_n2noreact
#SBATCH -p amd_m8_768
#SBATCH -N 1
#SBATCH --ntasks-per-node=192
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --exclusive

ulimit -s unlimited
ulimit -l unlimited

source /public1/home/m8s000774/users/xiao_chen_xiang/package/mbbq/OpenFOAM-v1706/etc/rebashrc

sh clean.sh

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

runApplication dsmcInitialise+

