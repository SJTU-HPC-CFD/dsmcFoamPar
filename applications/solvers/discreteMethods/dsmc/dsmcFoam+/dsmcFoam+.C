/*---------------------------------------------------------------------------*\\
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

Description
    Stage-1 core DSMC solver for the v2506 migration.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dsmcCloud.H"
#include "dsmcDynamicLoadBalancing.H"

#include <cstdlib>
#include <chrono>

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< nl << "Constructing dsmcCloud" << endl;
    dsmcCloud dsmc(runTime, "dsmc", mesh);
    dsmcDynamicLoadBalancing loadBalancer(runTime, mesh, dsmc);

    Info<< "\nStarting time loop\n" << endl;

    label infoCounter = 0;
    const bool solverStageProbe =
        runTime.controlDict().lookupOrDefault<bool>("solverStageProbe", false);
    const auto loopStart = std::chrono::steady_clock::now();

    while (runTime.loop())
    {
        ++infoCounter;
        const bool emitStepDiagnostics = (infoCounter >= dsmc.nTerminalOutputs());
        dsmc.setStepDiagnosticOutput(emitStepDiagnostics);

        if (emitStepDiagnostics && dsmc.isOutputRank())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": before evolve" << nl << endl;
        }
        dsmc.evolve();

        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": after evolve" << nl << endl;
        }
        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": before loadBalancer.update" << nl << endl;
        }
        loadBalancer.update();
        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": after loadBalancer.update" << nl << endl;
        }
        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": before loadBalancer.perform" << nl << endl;
        }
        loadBalancer.perform();
        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": after loadBalancer.perform" << nl << endl;
        }

        if (emitStepDiagnostics)
        {
            dsmc.info();
            infoCounter = 0;
        }

        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": before runTime.write" << nl << endl;
        }
        if (dsmc.replicatedMeshActive())
        {
            const bool isOutput = runTime.outputTime();
            if (isOutput)
            {
                dsmc.replicatedMeshRef().gatherParcelsToRank0();
            }
            if (dsmc.isOutputRank())
            {
                runTime.write();
            }
            if (isOutput)
            {
                dsmc.replicatedMeshRef().migrateParticlesByCellOwner();
            }
        }
        else
        {
            runTime.write();
        }
        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": after runTime.write" << nl << endl;
        }
        if (dsmc.isOutputRank())
        {
            runTime.printExecutionTime(Info);
        }
    }

    const scalar mainLoopWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            std::chrono::steady_clock::now() - loopStart
        ).count();

    if (dsmc.isOutputRank())
    {
        Info<< nl
            << "Main loop profiling summary:" << nl
            << "    main loop wall time [s]      = "
            << mainLoopWallTime << nl
            << endl;
    }

    dsmc.reportProfiling();

    if (dsmc.isOutputRank())
    {
        Info<< "End\n" << endl;
    }

    if (dsmc.replicatedMeshActive() && !Pstream::parRun())
    {
        MPI_Finalize();
    }

    UPstream::exit(0);
}

// ************************************************************************* //


