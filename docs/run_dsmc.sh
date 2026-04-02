${FOAM_USER_APPBIN}/dsmcInitialise+ > log.dsmcInitialise+
sh update.sh
mpirun -np 8 ${FOAM_USER_APPBIN}/dsmcFoam+ -parallel > log.dsmcFoam+
reconstructPar -latestTime 
foamToEnsight  -latestTime 