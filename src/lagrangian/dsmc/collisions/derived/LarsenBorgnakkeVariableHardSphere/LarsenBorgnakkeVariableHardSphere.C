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

\*---------------------------------------------------------------------------*/

#include "LarsenBorgnakkeVariableHardSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include <chrono>

namespace
{
    inline std::chrono::steady_clock::time_point wallClockNowLBVHS()
    {
        return std::chrono::steady_clock::now();
    }

    inline Foam::scalar wallSecondsLBVHS
    (
        const std::chrono::steady_clock::time_point& start,
        const std::chrono::steady_clock::time_point& stop
    )
    {
        return std::chrono::duration<Foam::scalar>(stop - start).count();
    }
}

namespace Foam
{
    defineTypeNameAndDebug(LarsenBorgnakkeVariableHardSphere, 0);
    addToRunTimeSelectionTable(BinaryCollisionModel, LarsenBorgnakkeVariableHardSphere, dictionary);
}

Foam::LarsenBorgnakkeVariableHardSphere::LarsenBorgnakkeVariableHardSphere
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    VariableHardSphere(dict, cloud),
    coeffDictLB_(dict.isDict(typeName + "Coeffs") ? dict.subDict(typeName + "Coeffs") : dictionary()),
    rotationalRelaxationCollisionNumber_(coeffDictLB_.lookupOrDefault<scalar>("rotationalRelaxationCollisionNumber", 5.0)),
    vibrationalRelaxationCollisionNumber_(coeffDictLB_.lookupOrDefault<scalar>("vibrationalRelaxationCollisionNumber", 0.0)),
    invZvFormulation_(2),
    electronicRelaxationCollisionNumber_(coeffDictLB_.lookupOrDefault<scalar>("electronicRelaxationCollisionNumber", 500.0)),
    detailRedistributePWallTime_(0.0),
    detailRedistributeQWallTime_(0.0),
    detailScatterWallTime_(0.0),
    detailCollideCalls_(0)
{
    const word inverseZvFormulationVersion(coeffDictLB_.lookupOrDefault<word>("inverseZvFormulation", word::null));

    if (inverseZvFormulationVersion == "pre-2008")
    {
        invZvFormulation_ = 0;
    }
    else if (inverseZvFormulationVersion == "2008")
    {
        invZvFormulation_ = 1;
    }
}

Foam::LarsenBorgnakkeVariableHardSphere::~LarsenBorgnakkeVariableHardSphere()
{
    if (cloud_.profilingDetailEnabled() && detailCollideCalls_ > 0)
    {
        Info<< nl
            << "LarsenBorgnakkeVariableHardSphere detail profile:" << nl
            << "    collide calls             = " << detailCollideCalls_ << nl
            << "    redistribute P wall [s]  = " << detailRedistributePWallTime_ << nl
            << "    redistribute Q wall [s]  = " << detailRedistributeQWallTime_ << nl
            << "    scatter wall [s]         = " << detailScatterWallTime_ << nl
            << endl;
    }
}

void Foam::LarsenBorgnakkeVariableHardSphere::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    scalar cR
)
{
    const bool profileDetail = cloud_.profilingDetailEnabled();
    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();
    const scalar mR = mP*mQ/(mP + mQ);
    scalar translationalEnergy = 0.5*mR*magSqr(UP - UQ);
    const scalar omegaPQ = 0.5*(cloud_.constProps(typeIdP).omega() + cloud_.constProps(typeIdQ).omega());

    if (profileDetail)
    {
        const auto redistributePStart = wallClockNowLBVHS();
        redistribute(pP, translationalEnergy, omegaPQ);
        detailRedistributePWallTime_ +=
            wallSecondsLBVHS(redistributePStart, wallClockNowLBVHS());

        const auto redistributeQStart = wallClockNowLBVHS();
        redistribute(pQ, translationalEnergy, omegaPQ);
        detailRedistributeQWallTime_ +=
            wallSecondsLBVHS(redistributeQStart, wallClockNowLBVHS());
    }
    else
    {
        redistribute(pP, translationalEnergy, omegaPQ);
        redistribute(pQ, translationalEnergy, omegaPQ);
    }

    translationalEnergy = max(translationalEnergy, 0.0);
    cR = translationalEnergy > VSMALL ? sqrt(2.0*translationalEnergy/mR) : 0.0;

    if (profileDetail)
    {
        const auto scatterStart = wallClockNowLBVHS();
        VariableHardSphere::scatter(pP, pQ, cellI, cR);
        detailScatterWallTime_ += wallSecondsLBVHS(scatterStart, wallClockNowLBVHS());
        ++detailCollideCalls_;
    }
    else
    {
        VariableHardSphere::scatter(pP, pQ, cellI, cR);
    }
}

void Foam::LarsenBorgnakkeVariableHardSphere::redistribute
(
    dsmcParcel& p,
    scalar& translationalEnergy,
    const scalar omegaPQ,
    const bool postReaction
)
{
    const label typeId = p.typeId();
    const dsmcParcel::constantProperties& cP = cloud_.constProps(typeId);

    if (cP.type() == 0)
    {
        return;
    }

    const scalar inverseRotationalCollisionNumber = 1.0/rotationalRelaxationCollisionNumber_;
    const scalar inverseElectronicCollisionNumber = 1.0/electronicRelaxationCollisionNumber_;

    scalar& ERot = p.ERot();
    label& ELevel = p.ELevel();

    if (inverseElectronicCollisionNumber > cloud_.rndGen().sample01<scalar>())
    {
        const scalar preCollisionEEle = cP.electronicEnergyList()[ELevel];
        const scalar Ec = translationalEnergy + preCollisionEEle;

        ELevel = cloud_.postCollisionElectronicEnergyLevel(Ec, cP.nElectronicLevels(), omegaPQ, cP.electronicEnergyList(), cP.electronicDegeneracyList());
        translationalEnergy = max(Ec - cP.electronicEnergyList()[ELevel], 0.0);
    }

    if (cP.nVibrationalModes() > 0)
    {
        const scalarList preCollisionEVib = cP.eVib(p.vibLevel());

        forAll(cP.thetaV(), i)
        {
            const scalar Ec = translationalEnergy + preCollisionEVib[i];
            const label iMax = Ec/(constant::physicoChemical::k.value()*cP.thetaV_m(i));

            if (iMax > 0)
            {
                p.vibLevel()[i] = cloud_.postCollisionVibrationalEnergyLevel(postReaction, p.vibLevel()[i], iMax, cP.thetaV_m(i), cP.thetaD(), cP.TrefZv_m(i), omegaPQ, cP.Zref_m(i), Ec, vibrationalRelaxationCollisionNumber_, invZvFormulation_, p.cell());
                translationalEnergy = max(Ec - cP.eVib_m(i, p.vibLevel()[i]), 0.0);
            }
        }
    }

    const scalar rotationalDof = cP.rotationalDegreesOfFreedom();

    if (rotationalDof > 0 && inverseRotationalCollisionNumber > cloud_.rndGen().sample01<scalar>())
    {
        const scalar Ec = translationalEnergy + ERot;
        const scalar ChiB = 2.5 - omegaPQ;
        ERot = cloud_.postCollisionRotationalEnergy(rotationalDof, ChiB)*Ec;
        translationalEnergy = max(Ec - ERot, 0.0);
    }
}

const Foam::dictionary& Foam::LarsenBorgnakkeVariableHardSphere::coeffDict() const
{
    return coeffDictLB_;
}

// ************************************************************************* //
