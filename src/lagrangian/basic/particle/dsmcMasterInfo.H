/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\      /   O peration   |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

Description
    One-copy diagnostics helper for replicated-mesh runs: without the
    -parallel flag every rank is an independent OpenFOAM "serial" process
    (Pstream::master() is true for all ranks, even under mpirun), so rank
    gating must query MPI / the launcher environment directly.  Safe
    before and after MPI_Init.

SourceFiles
    dsmcMasterInfo.H

\*---------------------------------------------------------------------------*/

#ifndef dsmcMasterInfo_H
#define dsmcMasterInfo_H

#include "mpi.h"
#include <cstdlib>
#include <cstdio>

namespace Foam
{

//- True if this process should print one-copy diagnostics (master rank,
//  or any rank when running genuinely serial)
inline bool dsmcIsPrintingRank()
{
    // 1. MPI initialised: query the communicator
    int init = 0;
    MPI_Initialized(&init);
    if (init)
    {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank == 0;
    }

    // 2. mpirun-launched but MPI not yet initialised (replicated-mesh
    //    construction phase): the launcher exports PMI_RANK
    const char* pmi = std::getenv("PMI_RANK");
    if (pmi)
    {
        return std::atoi(pmi) == 0;
    }

    // 3. genuinely serial
    return true;
}


}

#endif

// ************************************************************************* //
