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
    using clock_type = std::chrono::steady_clock;
    const auto loopStart = clock_type::now();

    scalar solverRunTimeLoopWallTime = 0.0;
    scalar solverHeaderWallTime = 0.0;
    scalar solverEvolveWallTime = 0.0;
    scalar solverLoadUpdateWallTime = 0.0;
    scalar solverLoadPerformWallTime = 0.0;
    scalar solverInfoWallTime = 0.0;
    scalar solverWriteWallTime = 0.0;
    scalar solverStepResidualWallTime = 0.0;
    label solverProfileSteps = 0;

    while (true)
    {
        const auto tRunTimeLoop0 = clock_type::now();
        if (!runTime.loop())
        {
            solverRunTimeLoopWallTime +=
                std::chrono::duration<scalar>(clock_type::now() - tRunTimeLoop0).count();
            break;
        }
        const auto tStepStart = clock_type::now();
        solverRunTimeLoopWallTime +=
            std::chrono::duration<scalar>(tStepStart - tRunTimeLoop0).count();

        ++infoCounter;
        const bool emitStepDiagnostics = (infoCounter >= dsmc.nTerminalOutputs());
        dsmc.setStepDiagnosticOutput(emitStepDiagnostics);

        if (emitStepDiagnostics && dsmc.isOutputRank())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }
        const auto tAfterHeader = clock_type::now();
        solverHeaderWallTime +=
            std::chrono::duration<scalar>(tAfterHeader - tStepStart).count();

        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": before evolve" << nl << endl;
        }
        const auto tEvolve0 = clock_type::now();
        dsmc.evolve();
        const auto tEvolve1 = clock_type::now();
        solverEvolveWallTime +=
            std::chrono::duration<scalar>(tEvolve1 - tEvolve0).count();

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
        const auto tLoadUpdate0 = clock_type::now();
        loadBalancer.update();
        const auto tLoadUpdate1 = clock_type::now();
        solverLoadUpdateWallTime +=
            std::chrono::duration<scalar>(tLoadUpdate1 - tLoadUpdate0).count();
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
        const auto tLoadPerform0 = clock_type::now();
        loadBalancer.perform();
        const auto tLoadPerform1 = clock_type::now();
        solverLoadPerformWallTime +=
            std::chrono::duration<scalar>(tLoadPerform1 - tLoadPerform0).count();
        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": after loadBalancer.perform" << nl << endl;
        }

        scalar stepInfoWallTime = 0.0;
        if (emitStepDiagnostics)
        {
            const auto tInfo0 = clock_type::now();
            dsmc.info();
            const auto tInfo1 = clock_type::now();
            stepInfoWallTime =
                std::chrono::duration<scalar>(tInfo1 - tInfo0).count();
            solverInfoWallTime += stepInfoWallTime;
            infoCounter = 0;
        }

        if (solverStageProbe)
        {
            Pout<< "Solver stage rank " << Pstream::myProcNo()
                << " timeIndex " << runTime.timeIndex()
                << " time " << runTime.timeName()
                << ": before runTime.write" << nl << endl;
        }
        const auto tWrite0 = clock_type::now();
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
                if (isOutput)
                {
                    dsmc.replicatedMeshRef().writeCellOwner();
                }
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
        const auto tStepEnd = clock_type::now();
        solverWriteWallTime +=
            std::chrono::duration<scalar>(tStepEnd - tWrite0).count();

        const scalar accountedStep =
            std::chrono::duration<scalar>(tAfterHeader - tStepStart).count()
          + std::chrono::duration<scalar>(tEvolve1 - tEvolve0).count()
          + std::chrono::duration<scalar>(tLoadUpdate1 - tLoadUpdate0).count()
          + std::chrono::duration<scalar>(tLoadPerform1 - tLoadPerform0).count()
          + stepInfoWallTime
          + std::chrono::duration<scalar>(tStepEnd - tWrite0).count();
        solverStepResidualWallTime +=
            std::chrono::duration<scalar>(tStepEnd - tStepStart).count()
          - accountedStep;
        ++solverProfileSteps;
    }

    const scalar mainLoopWallTime =
        std::chrono::duration_cast<std::chrono::duration<scalar>>
        (
            clock_type::now() - loopStart
        ).count();

    if (dsmc.isOutputRank())
    {
        const scalar solverAccountedWallTime =
            solverRunTimeLoopWallTime
          + solverHeaderWallTime
          + solverEvolveWallTime
          + solverLoadUpdateWallTime
          + solverLoadPerformWallTime
          + solverInfoWallTime
          + solverWriteWallTime
          + solverStepResidualWallTime;

        Info<< nl
            << "Main loop profiling summary:" << nl
            << "    main loop wall time [s]      = "
            << mainLoopWallTime << nl
            << "    solver profile steps          = "
            << solverProfileSteps << nl
            << "    runTime.loop [s]              = "
            << solverRunTimeLoopWallTime << nl
            << "    loop header/output [s]        = "
            << solverHeaderWallTime << nl
            << "    dsmc.evolve call [s]          = "
            << solverEvolveWallTime << nl
            << "    loadBalancer.update [s]       = "
            << solverLoadUpdateWallTime << nl
            << "    loadBalancer.perform [s]      = "
            << solverLoadPerformWallTime << nl
            << "    dsmc.info [s]                 = "
            << solverInfoWallTime << nl
            << "    runTime.write/print [s]       = "
            << solverWriteWallTime << nl
            << "    step residual [s]             = "
            << solverStepResidualWallTime << nl
            << "    solver accounted [s]          = "
            << solverAccountedWallTime << nl
            << "    main-accounted residual [s]   = "
            << mainLoopWallTime - solverAccountedWallTime << nl
            << endl;
    }

    dsmc.reportProfiling();

    if (dsmc.isOutputRank())
    {
        Info<< "End\n" << endl;
    }

    UPstream::exit(0);
}

// ************************************************************************* //
