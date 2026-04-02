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

#include "dsmcParcel.H"
#include "dsmcCloud.H"
#include "meshTools.H"

bool Foam::dsmcParcel::move
(
    dsmcCloud& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    if (isFree())
    {
        if (newParcel_ != -1)
        {
            stepFraction() = cloud.rndGen().sample01<scalar>();
            newParcel_ = -1;
        }

        vector Utracking = U_;

        while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
        {
            Utracking = U_;
            meshTools::constrainDirection(mesh(), mesh().solutionD(), Utracking);

            const vector d = deviationFromMeshCentre();
            const scalar f = 1 - stepFraction();
            trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);

            if (face() != -1)
            {
                cloud.tracker().trackParcelFaceTransition(*this);

                forAll(cloud.boundaries().cyclicBoundaryModels(), c)
                {
                    const labelList& faces = cloud.boundaries().cyclicBoundaryModels()[c]->allFaces();

                    if (Foam::hyCompat::indexOf(faces, face()) != -1)
                    {
                        cloud.boundaries().cyclicBoundaryModels()[c]->controlMol(*this, td);
                    }
                }
            }
        }
    }
    else
    {
        forAll(cloud.boundaries().patchBoundaryModels(), c)
        {
            if (cloud.boundaries().patchBoundaryModels()[c]->patchId() == stuck().wallTemperature()[1])
            {
                cloud.boundaries().patchBoundaryModels()[c]->controlParticle(*this, td);
                break;
            }
        }
    }

    return td.keepParticle;
}

bool Foam::dsmcParcel::hitPatch(dsmcCloud& cloud, trackingData& td)
{
    const label patchIndex = patch();

    if (patchIndex >= 0 && patchIndex < cloud.boundaries().patchToModelIds().size())
    {
        const label patchModelId = cloud.boundaries().patchToModelIds()[patchIndex];

        if (patchModelId >= 0)
        {
            cloud.boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
        }
    }

    return false;
}

void Foam::dsmcParcel::hitProcessorPatch(dsmcCloud&, trackingData& td)
{
    td.switchProcessor = true;
}

void Foam::dsmcParcel::hitWallPatch(dsmcCloud& cloud, trackingData& td)
{
    const label patchIndex = patch();

    if (patchIndex >= 0 && patchIndex < cloud.boundaries().patchToModelIds().size())
    {
        const label patchModelId = cloud.boundaries().patchToModelIds()[patchIndex];

        if (patchModelId >= 0)
        {
            cloud.boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
            return;
        }
    }

    cloud.handleWallInteraction(*this, td);
}

void Foam::dsmcParcel::transformProperties(const tensor& T)
{
    particle::transformProperties(T);
    U_ = transform(T, U_);
}

void Foam::dsmcParcel::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}

#include "dsmcParcelIO.C"

// ************************************************************************* //

Foam::label Foam::dsmcParcel::TrackedParcel::nDELETED = 0;



