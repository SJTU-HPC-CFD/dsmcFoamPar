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

Class
    dsmcVariableTimeStepModel

Description

\*----------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "dsmcVariableTimeStepModel.H"
#include "dsmcCloud.H"
#include <mpi.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcVariableTimeStepModel, 0);

    addToRunTimeSelectionTable
    (
        dsmcTimeStepModel,
        dsmcVariableTimeStepModel,
        fvMesh
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dsmcVariableTimeStepModel::findRefCell()
{
    // Find the cell with minimum volume
    // This cell will serve as a reference cell to set the
    // nParticle/timeStep ratio

    const scalarField& volumeCells = mesh_.V();
    scalar minVolume = gMin(volumeCells);

    if (Pstream::parRun())
    {
        reduce(minVolume, minOp<scalar>());
    }

    forAll(nParticles_, celli)
    {
        if (mag(volumeCells[celli] - minVolume) < SMALL)
        {
            refCell_ = celli;
            break;
        }
    }
}

void dsmcVariableTimeStepModel::findRefCellValue()
{
    const scalarField& volumeCells = mesh_.V();
    
    // 找到本地最小体积单元及其参数
    scalar localMinVolume = GREAT;
    scalar localNParticle = 0;
    scalar localDeltaT = 0;
    label localRefCell = -1;
    
    forAll(volumeCells, celli)
    {
        if (volumeCells[celli] < localMinVolume)
        {
            localMinVolume = volumeCells[celli];
            localNParticle = nParticles_[celli];
            localDeltaT = deltaT_[celli];
            localRefCell = celli;
        }
    }

    // 在并行情况下，获取全局最小值和对应的参数
    if (Pstream::parRun())
    {
        // 收集所有进程的最小体积和对应参数
        List<scalar> allMinVolumes(Pstream::nProcs());
        List<scalar> allNParticles(Pstream::nProcs());
        List<scalar> allDeltaTs(Pstream::nProcs());
        
        allMinVolumes[Pstream::myProcNo()] = localMinVolume;
        allNParticles[Pstream::myProcNo()] = localNParticle;
        allDeltaTs[Pstream::myProcNo()] = localDeltaT;
        
        // 找到全局最小体积
        reduce(allMinVolumes, minOp<List<scalar>>());
        
        // 找到拥有全局最小体积的进程
        label minProc = findMin(allMinVolumes);
        
        // 准备要广播的数据
        scalar broadcastData[3];
        if (minProc == Pstream::myProcNo())
        {
            broadcastData[0] = localMinVolume;
            broadcastData[1] = localNParticle;
            broadcastData[2] = localDeltaT;
            
            Info<< "Variable time-step model info from ref cell:" << nl
                << "- Reference cell found on processor " << minProc << nl
                << "- Minimum cell volume: " << localMinVolume << nl
                << "- Reference nParticle: " << localNParticle << nl
                << "- Reference deltaT: " << localDeltaT << nl
                << endl;
        }
        
        // 使用MPI进行广播
        MPI_Bcast
        (
            broadcastData,
            3,
            MPI_DOUBLE,
            minProc,
            MPI_COMM_WORLD
        );
        
        // 所有进程更新数据
        minVolumeRef_ = broadcastData[0];
        nParticleRef_ = broadcastData[1];
        deltaTRef_ = broadcastData[2];
        
        // 其他进程输出信息
        if (minProc != Pstream::myProcNo())
        {
            Info<< "Variable time-step model info from other cells:" << nl
                << "- Reference cell found on processor " << minProc << nl
                << "- Minimum cell volume: " << minVolumeRef_ << nl
                << "- Reference nParticle: " << nParticleRef_ << nl
                << "- Reference deltaT: " << deltaTRef_ << nl
                << endl;
        }
    }
    else
    {
        nParticleRef_ = localNParticle;
        deltaTRef_ = localDeltaT;
        minVolumeRef_ = localMinVolume;
        
        Info<< "Variable time-step model info:" << nl
            << "- Reference cell: " << localRefCell << nl
            << "- Minimum cell volume: " << localMinVolume << nl
            << "- Reference nParticle: " << localNParticle << nl
            << "- Reference deltaT: " << localDeltaT << nl
            << endl;
    }
}


void dsmcVariableTimeStepModel::updatenParticles()
{
    findRefCellValue();

    const scalarField& volumeCells = mesh_.V();
    const scalar minVolume = minVolumeRef_;

    const scalar nParticleRef = nParticleRef_;

    forAll(nParticles_, celli)
    {
        nParticles_[celli] = nParticleRef*volumeCells[celli]/minVolume;
    }

    forAll(nParticles_.boundaryField(), patchi)
    {
        fvPatchScalarField& pnParticles =
            nParticles_.boundaryFieldRef()[patchi];

        forAll(pnParticles, facei)
        {
            pnParticles[facei] =
                nParticles_[mesh_.boundaryMesh()[patchi].faceCells()[facei]];
        }
    }
}


void dsmcVariableTimeStepModel::updateTimeStep()
{
    const scalar nParticleTimeStepRatio =
        nParticleRef_/deltaTRef_;

    forAll(deltaT_, celli)
    {
        deltaT_[celli] = nParticles_[celli]/nParticleTimeStepRatio;
    }
}


void dsmcVariableTimeStepModel::updateVariableTimeStepMethod()
{
    updatenParticles();

    updateTimeStep();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcVariableTimeStepModel::dsmcVariableTimeStepModel
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    dsmcTimeStepModel(t, mesh, cloud),
    cloud_(cloud),
    refCell_(-1),
    deltaT_
    (
        IOobject
        (
            "deltaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("deltaT", dimTime, deltaTValueOrg())
    )
{
    nParticles_.writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVariableTimeStepModel::~dsmcVariableTimeStepModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcVariableTimeStepModel::checkTimeStepModelInputs()
{
    const bool nParticlesFromFile =
        cloud_.particleProperties().lookupOrDefault<bool>
        (
            "nEquivalentParticlesFromFile",
            false
        );

    if(nParticlesFromFile)
    {
        nParticles_.regIOobject::read();

        findRefCellValue();
    }
    else
    {
        updatenParticles();
    }

    updateTimeStep();

    writeTimeStepModelInfo();
}


void dsmcVariableTimeStepModel::update()
{
    updateVariableTimeStepMethod();
}


void dsmcVariableTimeStepModel::writeTimeStepModelInfo() const
{
    Info<< "Variable time-step model:" << nl
        << "- Initial time-step [sec]" << tab
        << dsmcTimeStepModel::deltaTValue(0) << nl
        << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
