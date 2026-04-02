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

#include <cstdlib>

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< nl << "Constructing dsmcCloud" << endl;
    dsmcCloud dsmc(runTime, "dsmc", mesh);

    Info<< "\nStarting time loop\n" << endl;

    label infoCounter = 0;

    while (runTime.loop())
    {
        ++infoCounter;

        if (infoCounter >= dsmc.nTerminalOutputs())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
        }

        dsmc.evolve();

        if (infoCounter >= dsmc.nTerminalOutputs())
        {
            dsmc.info();
            infoCounter = 0;
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    // The migrated lagrangian stack still has a destructor-time cleanup issue.
    // Exit directly after a successful run to avoid tearing down stale state.
    UPstream::exit(0);
}

// ************************************************************************* //


