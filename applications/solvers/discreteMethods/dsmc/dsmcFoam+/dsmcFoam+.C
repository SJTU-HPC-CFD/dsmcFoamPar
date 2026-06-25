/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    dsmcFoam+
    dsmcFoam+ -AMR

Description
    Direct simulation Monte Carlo (DSMC) solver for 3D, transient, multi-
    species flows

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dsmcCloud.H"
#include <chrono>

namespace
{
    typedef std::chrono::steady_clock steadyWallClock;

    scalar wallElapsed(const steadyWallClock::time_point& start)
    {
        return std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            steadyWallClock::now() - start
        ).count();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool run
(
    argList& args,
    scalar& currentIterationTime,
    scalar& previousIterationTime,
    label& noRestart,
    label& noIteration
);


int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "AMR",
        "activate Adaptive Mesh Refinement"
    );

    #include "setRootCase.H"

    scalar currentIterationTime = 0.0;
    scalar previousIterationTime = 1.0;

    label noRestart = 0;
    label noIteration = 1;
    label totNoIteration = 0;
    bool restart = false;

    do
    {
        noRestart += 1;

        restart = run
        (
            args,
            currentIterationTime,
            previousIterationTime,
            noRestart,
            noIteration
        );

        totNoIteration += noIteration - 1;
        noIteration = 1;

    } while(restart);

    Info<< "Total Iterations = " << totNoIteration << "\n"
        << "End main\n" << endl;

    if (Pstream::parRun())
    {
        Pstream::exit(0);
    }

    return 0;
}


bool run
(
    argList& args,
    scalar& currentIterationTime,
    scalar& previousIterationTime,
    label& noRestart,
    label& noIteration
)
{
    #include "createTime.H"

    const bool activateAMR = args.optionFound("AMR");

    #include "createDynamicFvMesh.H"

    Info<< nl << "Constructing dsmcCloud " << endl;

    dsmcCloud dsmc(runTime, "dsmc", mesh);

    Info<< "\nStarting time loop\n" << endl;

    label infoCounter = 0;
    label noRefinement = 0;
    const bool profileSummary = dsmc.profileSummaryEnabled();
    scalar loopWall = 0.0;
    scalar amrWall = 0.0;
    scalar timeBannerWall = 0.0;
    scalar evolveWall = 0.0;
    scalar infoWall = 0.0;
    scalar asyncFinishWall = 0.0;
    scalar outputOwnerAlignWall = 0.0;
    scalar outputGatherWall = 0.0;
    scalar outputWriteWall = 0.0;
    scalar outputCellOwnerWall = 0.0;
    scalar outputMigrateBackWall = 0.0;
    scalar serialWriteWall = 0.0;
    scalar iterationAccountingWall = 0.0;
    scalar stagePrintWall = 0.0;
    scalar loadBalanceCheckWall = 0.0;

    while (runTime.loop())
    {
        const steadyWallClock::time_point loopWallStart =
            steadyWallClock::now();

        if (activateAMR)
        {
            const steadyWallClock::time_point amrWallStart =
                steadyWallClock::now();
            scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

            mesh.update();

            if (mesh.changing())
            {
                 ++noRefinement;

                 Info<< "Execution time for mesh.update() = "
                     << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                     << " s" << endl;
            }
            if (profileSummary)
            {
                amrWall += wallElapsed(amrWallStart);
            }
        }

        infoCounter++;

        if (infoCounter >= dsmc.nTerminalOutputs())
        {
            const steadyWallClock::time_point timeBannerWallStart =
                steadyWallClock::now();
            Info<< "Time = " << runTime.timeName() << nl << endl;
            if (profileSummary)
            {
                timeBannerWall += wallElapsed(timeBannerWallStart);
            }
        }

        const steadyWallClock::time_point evolveWallStart =
            steadyWallClock::now();
        dsmc.evolve();
        if (profileSummary)
        {
            evolveWall += wallElapsed(evolveWallStart);
        }

        if (infoCounter >= dsmc.nTerminalOutputs())
        {
            const steadyWallClock::time_point infoWallStart =
                steadyWallClock::now();
            dsmc.info();
            if (profileSummary)
            {
                infoWall += wallElapsed(infoWallStart);
            }
        }

        if (dsmc.replicatedMeshActive())
        {
            if (dsmc.replicatedMeshRef().asyncMigrationPending())
            {
                const steadyWallClock::time_point asyncFinishWallStart =
                    steadyWallClock::now();
                dsmc.replicatedMeshRef().migrateFinish();
                dsmc.replicatedMeshRef().updateParticleCounts();
                dsmc.clearMoveOrderedParcels();
                if (profileSummary)
                {
                    asyncFinishWall += wallElapsed(asyncFinishWallStart);
                }
            }

            const bool isOutput = runTime.outputTime();
            scalar gatherWall = 0.0;
            scalar writeWall = 0.0;
            scalar migrateBackWall = 0.0;

            if (dsmc.replicatedMeshRef().processorWriteEnabled())
            {
                if (isOutput)
                {
                    const steadyWallClock::time_point tOwnerAlign0 =
                        steadyWallClock::now();
                    dsmc.replicatedMeshRef().migrateParticlesByCellOwner();
                    dsmc.replicatedMeshRef().updateParticleCounts();
                    gatherWall = wallElapsed(tOwnerAlign0);
                    outputOwnerAlignWall += gatherWall;

                    const steadyWallClock::time_point tWrite0 =
                        steadyWallClock::now();
                    dsmc.replicatedMeshRef().writeProcessorOutput();
                    writeWall = wallElapsed(tWrite0);
                    outputWriteWall += writeWall;
                }

                if (isOutput && dsmc.isOutputRank())
                {
                    Info<< "Replicated mesh processor output timing:" << nl
                        << "    owner-align migrate [s] = " << gatherWall << nl
                        << "    processor fields/cloud write [s] = "
                        << writeWall
                        << nl << endl;
                }
            }
            else
            {
                if (isOutput)
                {
                    const steadyWallClock::time_point tGather0 =
                        steadyWallClock::now();
                    dsmc.replicatedMeshRef().gatherParcelsToRank0();
                    gatherWall = wallElapsed(tGather0);
                    outputGatherWall += gatherWall;
                }

                if (dsmc.isOutputRank())
                {
                    IOobject::writeOption oldCloudWriteOpt = dsmc.writeOpt();
                    dsmc.writeOpt() = IOobject::NO_WRITE;

                    const steadyWallClock::time_point tWrite0 =
                        steadyWallClock::now();
                    runTime.write();
                    dsmc.writeOpt() = oldCloudWriteOpt;

                    if (isOutput)
                    {
                        dsmc.replicatedMeshRef().writeGatheredCloudOnRank0();
                    }

                    writeWall = wallElapsed(tWrite0);
                    outputWriteWall += writeWall;

                    if (isOutput)
                    {
                        const steadyWallClock::time_point tCellOwner0 =
                            steadyWallClock::now();
                        dsmc.replicatedMeshRef().writeCellOwner();
                        outputCellOwnerWall += wallElapsed(tCellOwner0);
                    }
                }

                if (isOutput)
                {
                    const steadyWallClock::time_point tMigrateBack0 =
                        steadyWallClock::now();
                    dsmc.replicatedMeshRef().migrateParticlesByCellOwner();
                    dsmc.replicatedMeshRef().updateParticleCounts();
                    migrateBackWall = wallElapsed(tMigrateBack0);
                    outputMigrateBackWall += migrateBackWall;
                }

                if (isOutput && dsmc.isOutputRank())
                {
                    Info<< "Replicated mesh output timing:" << nl
                        << "    gather parcels [s]       = " << gatherWall << nl
                        << "    rank0 write [s]           = " << writeWall << nl
                        << "    migrate-back [s]          = " << migrateBackWall
                        << nl << endl;
                }
            }
        }
        else
        {
            const steadyWallClock::time_point serialWriteWallStart =
                steadyWallClock::now();
            runTime.write();
            if (profileSummary)
            {
                serialWriteWall += wallElapsed(serialWriteWallStart);
            }
        }

        const steadyWallClock::time_point iterationAccountingWallStart =
            steadyWallClock::now();
        previousIterationTime =
            max(runTime.elapsedCpuTime() - currentIterationTime, 1e-3);
        if (profileSummary)
        {
            iterationAccountingWall += wallElapsed(iterationAccountingWallStart);
        }

        if (infoCounter >= dsmc.nTerminalOutputs())
        {
            const steadyWallClock::time_point stagePrintWallStart =
                steadyWallClock::now();
            Info<< nl << "Stage " << noRestart << "." << noRefinement
                << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << "  Iteration " << noIteration << " ("
                << previousIterationTime << " s)"
                << nl << endl;

            infoCounter = 0;
            if (profileSummary)
            {
                stagePrintWall += wallElapsed(stagePrintWallStart);
            }
        }

        currentIterationTime = runTime.elapsedCpuTime();
        noIteration += 1;

        const steadyWallClock::time_point loadBalanceCheckWallStart =
            steadyWallClock::now();
        dsmc.loadBalanceCheck();
        if (profileSummary)
        {
            loadBalanceCheckWall += wallElapsed(loadBalanceCheckWallStart);
            loopWall += wallElapsed(loopWallStart);
        }
    }

    Info<< "End stage " << noRestart << "\n" << endl;

    if (profileSummary)
    {
        const scalar localAccounted =
            amrWall + timeBannerWall + evolveWall + infoWall
          + asyncFinishWall + outputOwnerAlignWall + outputGatherWall
          + outputWriteWall + outputCellOwnerWall + outputMigrateBackWall
          + serialWriteWall + iterationAccountingWall + stagePrintWall
          + loadBalanceCheckWall;

        scalar values[] =
        {
            loopWall,
            amrWall,
            timeBannerWall,
            evolveWall,
            infoWall,
            asyncFinishWall,
            outputOwnerAlignWall,
            outputGatherWall,
            outputWriteWall,
            outputCellOwnerWall,
            outputMigrateBackWall,
            serialWriteWall,
            iterationAccountingWall,
            stagePrintWall,
            loadBalanceCheckWall,
            localAccounted,
            loopWall - localAccounted
        };

        if (dsmc.replicatedMeshActive() && dsmc.replicatedMeshRef().nProcs() > 1)
        {
            int mpiInit = 0;
            MPI_Initialized(&mpiInit);
            if (mpiInit)
            {
                MPI_Allreduce
                (
                    MPI_IN_PLACE,
                    values,
                    17,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD
                );
            }
        }
        else if (Pstream::parRun())
        {
            for (label i = 0; i < 17; ++i)
            {
                reduce(values[i], maxOp<scalar>());
            }
        }

        const bool printProfile =
            (dsmc.replicatedMeshActive() && dsmc.replicatedMeshRef().myRank() == 0)
         || (!dsmc.replicatedMeshActive() && Pstream::master());

        if (printProfile)
        {
            const scalar componentMaxSum =
                values[1] + values[2] + values[3] + values[4] + values[5]
              + values[6] + values[7] + values[8] + values[9] + values[10]
              + values[11] + values[12] + values[13] + values[14];

            Info<< nl
                << "DSMC solver loop profile summary" << nl
                << "    loop wall max [s]             = " << values[0] << nl
                << "    AMR max [s]                   = " << values[1] << nl
                << "    time banner max [s]           = " << values[2] << nl
                << "    dsmc.evolve max [s]           = " << values[3] << nl
                << "    dsmc.info max [s]             = " << values[4] << nl
                << "    async finish max [s]          = " << values[5] << nl
                << "    output owner-align max [s]    = " << values[6] << nl
                << "    output gather max [s]         = " << values[7] << nl
                << "    output write max [s]          = " << values[8] << nl
                << "    output cellOwner max [s]      = " << values[9] << nl
                << "    output migrate-back max [s]   = " << values[10] << nl
                << "    serial runTime.write max [s]  = " << values[11] << nl
                << "    iteration accounting max [s]  = " << values[12] << nl
                << "    stage print max [s]           = " << values[13] << nl
                << "    loadBalanceCheck max [s]      = " << values[14] << nl
                << "    loop accounted local max [s]  = " << values[15] << nl
                << "    loop residual local max [s]   = " << values[16] << nl
                << "    loop component max sum [s]    = "
                << componentMaxSum << nl
                << endl;
        }
    }

    dsmc.printProfileSummary();

    if (dsmc.dynamicLoadBalancing().performBalance())
    {
        dsmc.loadBalance(noRefinement);
        return true;
    }

    return false;
}


// ************************************************************************* //
