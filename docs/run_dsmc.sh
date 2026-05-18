source ~/intel/oneapi/setvars.sh >/dev/null 2>&1
source /home/superxcx/code/OpenFoam/OF-2506/OpenFOAM-v2506/etc/bashrc >/dev/null 2>&1
export WM_PROJECT_USER_DIR=/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export LD_LIBRARY_PATH=${FOAM_USER_LIBBIN}:$LD_LIBRARY_PATH

${FOAM_USER_APPBIN}/dsmcInitialise+ > log.dsmcInitialise+
sh updateDt.sh 
OMP_NUM_THREADS=1 mpirun -np 8 ${FOAM_USER_APPBIN}/dsmcFoam+ -parallel > log.dsmcFoam+
OMP_NUM_THREADS=8 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+
OMP_NUM_THREADS=1 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+

OMP_NUM_THREADS=1 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.series-old 2>&1 &
OMP_NUM_THREADS=1 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.series 2>&1 &
OMP_NUM_THREADS=2 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp2 2>&1 &
OMP_NUM_THREADS=4 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp4 2>&1 &
OMP_NUM_THREADS=8 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp8 2>&1 &
OMP_NUM_THREADS=16 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp16 2>&1 &
OMP_NUM_THREADS=32 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp32 2>&1 &
OMP_NUM_THREADS=64 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp64 2>&1 &

OMP_NUM_THREADS=64 OMP_PROC_BIND=close OMP_PLACES=cores nohup srun --cpu-bind=cores ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp64.bind 2>&1 &
OMP_NUM_THREADS=64 OMP_PROC_BIND=close OMP_PLACES=cores nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp64.bind 2>&1 &

OMP_NUM_THREADS=2 OMP_PROC_BIND=close OMP_PLACES=cores  nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp2.bind.occ 2>&1 &
OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores  nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp4.bind.occ 2>&1 &
OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores  nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp8.bind.occ 2>&1 &
OMP_NUM_THREADS=16 OMP_PROC_BIND=close OMP_PLACES=cores nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp16.bind.occ 2>&1 &
OMP_NUM_THREADS=32 OMP_PROC_BIND=close OMP_PLACES=cores nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp32.bind.occ 2>&1 &
OMP_NUM_THREADS=64 OMP_PROC_BIND=close OMP_PLACES=cores nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp64.bind.occ 2>&1 &

OMP_NUM_THREADS=1 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.series-old 2>&1 &
OMP_NUM_THREADS=1 nohup ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.series 2>&1 &
OMP_NUM_THREADS=1 nohup mpirun -np 2 ${FOAM_USER_APPBIN}/dsmcFoam+ -parallel > log.dsmcFoam+.mpi2 2>&1 &
OMP_NUM_THREADS=1 nohup mpirun -np 4 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.mpi4 2>&1 &
OMP_NUM_THREADS=1 nohup mpirun -np 8 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.mpi8 2>&1 &
OMP_NUM_THREADS=1 nohup mpirun -np 16 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.mpi16 2>&1 &
OMP_NUM_THREADS=1 nohup mpirun -np 32 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.mpi32 2>&1 &
OMP_NUM_THREADS=1 nohup mpirun -np 64 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.mpi64 2>&1 &

OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores mpirun -n 2 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+.omp4.mpi2

reconstructPar -latestTime 
foamToEnsight  -latestTime 


source /home/superxcx/intel/oneapi/setvars.sh --force 2>/dev/null; source /home/superxcx/code/OpenFoam/OF-2506/OpenFOAM-v2506/etc/bashrc; export WM_PROJECT_USER_DIR=/home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx; export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin; export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib; export PATH=$FOAM_USER_APPBIN:$PATH; export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:/home/superxcx/code/DSMC/dsmcFoam++/parmetis-install/lib:/home/superxcx/miniconda3/envs/cap/lib:$LD_LIBRARY_PATH; cd /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinder_react/mixparallel/allmesh/omp4_mpi2_mpiopt_noreact; rm -rf processor*/0/lagrangian processor*/0/uniform; for i in 0 1; do cp -r 0/lagrangian processor$i/0/; done; export OMP_NUM_THREADS=4; /usr/bin/time -v mpirun -np 2 dsmcFoam+ -parallel > log.dsmcFoam+.fix_all 2>&1