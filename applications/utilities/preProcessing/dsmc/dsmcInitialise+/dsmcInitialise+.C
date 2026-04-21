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
    dsmcInitialise+

Description
    Minimal initialiser for the v2506 dsmcFoam+ migration.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dsmcCloud.H"
#include "OFstream.H"

namespace
{

void ensureEmptyCloudBootstrap
(
    const Foam::Time& runTime,
    const Foam::word& cloudName
)
{
    using namespace Foam;

    const fileName timeDir(runTime.path()/runTime.timeName());
    const fileName lagrangianDir(timeDir/cloud::prefix/cloudName);
    const fileName uniformDir(timeDir/"uniform");
    const fileName uniformLagrangianDir(uniformDir/cloud::prefix/cloudName);

    mkDir(timeDir);
    mkDir(uniformDir);
    mkDir(uniformDir/cloud::prefix);
    mkDir(uniformLagrangianDir);
    mkDir(timeDir/cloud::prefix);
    mkDir(lagrangianDir);

    const fileName cloudPropertiesFile(uniformLagrangianDir/"cloudProperties");
    if (!isFile(cloudPropertiesFile))
    {
        OFstream os(cloudPropertiesFile);

        os  << "FoamFile\n"
            << "{\n"
            << "    version     2.0;\n"
            << "    format      ascii;\n"
            << "    class       dictionary;\n"
            << "    location    \"" << runTime.timeName()
            << "/uniform/" << cloud::prefix << '/' << cloudName << "\";\n"
            << "    object      cloudProperties;\n"
            << "}\n\n"
            << "geometry        coordinates;\n\n"
            << "processor0\n"
            << "{\n"
            << "    particleCount   0;\n"
            << "}\n";
    }

    const fileName positionsFile(lagrangianDir/"positions");
    if (!isFile(positionsFile))
    {
        OFstream os(positionsFile);

        os  << "FoamFile\n"
            << "{\n"
            << "    version     2.0;\n"
            << "    format      ascii;\n"
            << "    class       Cloud;\n"
            << "    location    \"" << runTime.timeName()
            << '/' << cloud::prefix << '/' << cloudName << "\";\n"
            << "    object      positions;\n"
            << "}\n\n"
            << "0\n"
            << "(\n"
            << ")\n";
    }
}

}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dsmcInitialiseDict
    (
        IOobject
        (
            "dsmcInitialiseDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Initialising dsmc for Time = " << runTime.timeName() << nl << endl;

    ensureEmptyCloudBootstrap(runTime, "dsmc");

    dsmcCloud dsmc(runTime, "dsmc", mesh, false);
    dsmc.initialiseFromDict(dsmcInitialiseDict);

    Info<< "Initialised DSMC parcels = " << dsmc.nParcels() << nl << endl;


    if (!runTime.writeNow())
    {
        FatalErrorIn(args.executable())
            << "Failed writing initialised DSMC fields." << nl
            << exit(FatalError);
    }

    Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //


