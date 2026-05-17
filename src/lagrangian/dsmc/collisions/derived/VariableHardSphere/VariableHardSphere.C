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

#include "VariableHardSphere.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

namespace Foam
{
    defineTypeNameAndDebug(VariableHardSphere, 0);
    addToRunTimeSelectionTable
    (
        BinaryCollisionModel,
        VariableHardSphere,
        dictionary
    );
};

Foam::VariableHardSphere::VariableHardSphere
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    BinaryCollisionModel(dict, cloud),
    coeffDict_
    (
        dict.isDict(typeName + "Coeffs")
      ? dict.subDict(typeName + "Coeffs")
      : dictionary()
    ),
    Tref_(coeffDict_.lookupOrDefault<scalar>("Tref", 273.0)),
    pairTablesReady_(false),
    nTypes_(0)
{}


Foam::VariableHardSphere::~VariableHardSphere()
{}


bool Foam::VariableHardSphere::active() const
{
    return true;
}


Foam::scalar Foam::VariableHardSphere::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
    const label typeIdP = pP.typeId();
    const label typeIdQ = pQ.typeId();

    // cr2 = |Up-Uq|^2 — avoids sqrt in mag()
    const vector& Up = pP.U();
    const vector& Uq = pQ.U();
    const scalar cr2 =
        (Up.x() - Uq.x())*(Up.x() - Uq.x())
      + (Up.y() - Uq.y())*(Up.y() - Uq.y())
      + (Up.z() - Uq.z())*(Up.z() - Uq.z());

    if (cr2 < VSMALL)
    {
        return 0.0;
    }

    // Lazy-init per-species-pair sigmaTcR lookup tables (thread-safe, once)
    if (!pairTablesReady_)
    {
        #pragma omp critical(VhsPrecompute)
        {
            if (!pairTablesReady_)
            {
                const label nt = cloud_.constProps().size();
                nTypes_ = nt;
                const label K = cr2TableSize_;
                pairCr2Table_.setSize(nt*nt*K);

                const scalar k = physicoChemical::k.value();
                const scalar cr2Min = 1e-2;
                const scalar cr2Max = 1e8;
                cr2LogMin_ = log(cr2Min);
                cr2InvDLog_ = (K - 1) / (log(cr2Max) - cr2LogMin_);

                for (label i = 0; i < nt; ++i)
                {
                    const scalar di = cloud_.constProps(i).d();
                    const scalar omegai = cloud_.constProps(i).omega();
                    const scalar mi = cloud_.constProps(i).mass();

                    for (label j = 0; j < nt; ++j)
                    {
                        const scalar dj = cloud_.constProps(j).d();
                        const scalar omegaj = cloud_.constProps(j).omega();
                        const scalar mj = cloud_.constProps(j).mass();

                        const scalar dAvg = 0.5*(di + dj);
                        const scalar omegaAvg = 0.5*(omegai + omegaj);
                        const scalar mR = mi*mj/(mi + mj);

                        const scalar base = 2.0*k*Tref_/mR;
                        const scalar denom = exp(Foam::lgamma(2.5 - omegaAvg));

                        const scalar A = pi*sqr(dAvg)*pow(base, omegaAvg - 0.5)/denom;
                        const scalar Bhalf = 1.0 - omegaAvg;
                        const label baseIdx = (i*nt + j)*K;

                        for (label m = 0; m < K; ++m)
                        {
                            const scalar cr2m =
                                cr2Min * exp((scalar(m) / (K - 1))
                                * (log(cr2Max) - cr2LogMin_));
                            pairCr2Table_[baseIdx + m] = A * pow(cr2m, Bhalf);
                        }
                    }
                }
                pairTablesReady_ = true;
            }
        }
    }

    // Linear interpolation in log-spaced cr2 table
    const scalar logCr2 = log(cr2);
    const scalar pos = (logCr2 - cr2LogMin_) * cr2InvDLog_;
    const label K = cr2TableSize_;
    label k = label(pos);
    if (k < 0) k = 0;
    if (k >= K - 1) k = K - 2;
    const scalar t = pos - scalar(k);

    const label baseIdx = (typeIdP*nTypes_ + typeIdQ)*K;
    return pairCr2Table_[baseIdx + k]*(1.0 - t)
         + pairCr2Table_[baseIdx + k + 1]*t;
}


void Foam::VariableHardSphere::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label cellI,
    scalar cR
)
{
    scatter(pP, pQ, cellI, cR);
}


void Foam::VariableHardSphere::scatter
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    const label,
    scalar cR
)
{
    postCollisionVelocities
    (
        pP.typeId(),
        pQ.typeId(),
        pP.U(),
        pQ.U(),
        cR
    );

    const label classificationP = pP.classification();
    const label classificationQ = pQ.classification();

    if (classificationP == 0 && (classificationQ == 1 || classificationQ == 2))
    {
        pP.classification() = 2;
    }

    if (classificationQ == 0 && (classificationP == 1 || classificationP == 2))
    {
        pQ.classification() = 2;
    }
}


void Foam::VariableHardSphere::postCollisionVelocities
(
    const label typeIdP,
    const label typeIdQ,
    vector& UP,
    vector& UQ,
    scalar cR
)
{
    if (cR == -1)
    {
        cR = mag(UP - UQ);
    }

    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();
    const vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

    const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
    const scalar sinTheta = sqrt(1.0 - sqr(cosTheta));
    const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    const vector postCollisionRelativeU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    UP = Ucm + postCollisionRelativeU*mQ/(mP + mQ);
    UQ = Ucm - postCollisionRelativeU*mP/(mP + mQ);
}


void Foam::VariableHardSphere::postReactionVelocities
(
    const label typeIdP,
    const label typeIdQ,
    vector& UP,
    vector& UQ,
    scalar cR
)
{
    const scalar mP = cloud_.constProps(typeIdP).mass();
    const scalar mQ = cloud_.constProps(typeIdQ).mass();

    const scalar cosTheta = 2.0*cloud_.rndGen().sample01<scalar>() - 1.0;
    const scalar sinTheta = sqrt(1.0 - sqr(cosTheta));
    const scalar phi = twoPi*cloud_.rndGen().sample01<scalar>();

    const vector postCollisionRelativeU =
        cR
       *vector
        (
            cosTheta,
            sinTheta*cos(phi),
            sinTheta*sin(phi)
        );

    UQ = UP - postCollisionRelativeU*mP/(mP + mQ);
    UP += postCollisionRelativeU*mQ/(mP + mQ);
}


void Foam::VariableHardSphere::redistribute
(
    dsmcParcel&,
    scalar&,
    const scalar,
    const bool
)
{}


const Foam::dictionary& Foam::VariableHardSphere::coeffDict() const
{
    return coeffDict_;
}

// ************************************************************************* //
