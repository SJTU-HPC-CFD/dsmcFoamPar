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

    while (runTime.loop())
    {
        if (activateAMR)
        {
            scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

            mesh.update();

            if (mesh.changing())
            {
                 ++noRefinement;

                 Info<< "Execution time for mesh.update() = "
                     << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                     << " s" << endl;
            }
        }

        infoCounter++;

        if (infoCounter >= dsmc.nTerminalOutputs())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        dsmc.evolve();

        if (infoCounter >= dsmc.nTerminalOutputs())
        {
            dsmc.info();
        }

        if (dsmc.replicatedMeshActive())
        {
            if (dsmc.replicatedMeshRef().asyncMigrationPending())
            {
                dsmc.replicatedMeshRef().migrateFinish();
                dsmc.replicatedMeshRef().updateParticleCounts();
                dsmc.clearMoveOrderedParcels();
            }

            const bool isOutput = runTime.outputTime();
            scalar gatherWall = 0.0;
            scalar writeWall = 0.0;
            scalar migrateBackWall = 0.0;

            if (dsmc.replicatedMeshRef().processorWriteEnabled())
            {
                if (isOutput)
                {
                    const scalar tWrite0 = runTime.elapsedCpuTime();
                    dsmc.replicatedMeshRef().writeProcessorOutput();
                    writeWall = runTime.elapsedCpuTime() - tWrite0;
                }

                if (isOutput && dsmc.isOutputRank())
                {
                    Info<< "Replicated mesh processor output timing:" << nl
                        << "    processor write [s]    = " << writeWall
                        << nl << endl;
                }
            }
            else
            {
                if (isOutput)
                {
                    const scalar tGather0 = runTime.elapsedCpuTime();
                    dsmc.replicatedMeshRef().gatherParcelsToRank0();
                    gatherWall = runTime.elapsedCpuTime() - tGather0;
                }

                if (dsmc.isOutputRank())
                {
                    IOobject::writeOption oldCloudWriteOpt = dsmc.writeOpt();
                    dsmc.writeOpt() = IOobject::NO_WRITE;

                    const scalar tWrite0 = runTime.elapsedCpuTime();
                    runTime.write();
                    dsmc.writeOpt() = oldCloudWriteOpt;

                    if (isOutput)
                    {
                        dsmc.replicatedMeshRef().writeGatheredCloudOnRank0();
                    }

                    writeWall = runTime.elapsedCpuTime() - tWrite0;

                    if (isOutput)
                    {
                        dsmc.replicatedMeshRef().writeCellOwner();
                    }
                }

                if (isOutput)
                {
                    const scalar tMigrateBack0 = runTime.elapsedCpuTime();
                    dsmc.replicatedMeshRef().migrateParticlesByCellOwner();
                    dsmc.replicatedMeshRef().updateParticleCounts();
                    migrateBackWall = runTime.elapsedCpuTime() - tMigrateBack0;
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
            runTime.write();
        }

        previousIterationTime =
            max(runTime.elapsedCpuTime() - currentIterationTime, 1e-3);

        if (infoCounter >= dsmc.nTerminalOutputs())
        {
            Info<< nl << "Stage " << noRestart << "." << noRefinement
                << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << "  Iteration " << noIteration << " ("
                << previousIterationTime << " s)"
                << nl << endl;

            infoCounter = 0;
        }

        currentIterationTime = runTime.elapsedCpuTime();
        noIteration += 1;

        dsmc.loadBalanceCheck();
    }

    Info<< "End stage " << noRestart << "\n" << endl;

    dsmc.printProfileSummary();

    if (dsmc.dynamicLoadBalancing().performBalance())
    {
        dsmc.loadBalance(noRefinement);
        return true;
    }

    return false;
}


// ************************************************************************* //
