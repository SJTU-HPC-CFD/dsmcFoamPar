${FOAM_USER_APPBIN}/dsmcInitialise+ > log.dsmcInitialise+
sh updateDt.sh 
OMP_NUM_THREADS=1 mpirun -np 8 ${FOAM_USER_APPBIN}/dsmcFoam+ -parallel > log.dsmcFoam+
OMP_NUM_THREADS=8 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+
OMP_NUM_THREADS=1 ${FOAM_USER_APPBIN}/dsmcFoam+ > log.dsmcFoam+

reconstructPar -latestTime 
foamToEnsight  -latestTime 
