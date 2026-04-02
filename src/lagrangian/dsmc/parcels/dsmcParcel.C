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

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace
{
inline bool useOpenMPMoveCriticals(const Foam::dsmcCloud& cloud)
{
    #ifdef _OPENMP
    return cloud.openmpMoveEnabled() && omp_in_parallel();
    #else
    return false;
    #endif
}

inline void trackParcelFaceTransitionThreadSafe
(
    Foam::dsmcCloud& cloud,
    const Foam::dsmcParcel& p
)
{
    cloud.tracker().trackParcelFaceTransition(p);
}

inline void controlCyclicBoundaryThreadSafe
(
    Foam::dsmcCloud& cloud,
    const Foam::label boundaryI,
    Foam::dsmcParcel& p,
    Foam::dsmcParcel::trackingData& td
)
{
    if (useOpenMPMoveCriticals(cloud))
    {
        #pragma omp critical(dsmcMoveBoundary)
        {
            cloud.boundaries().cyclicBoundaryModels()[boundaryI]->controlMol(p, td);
        }
    }
    else
    {
        cloud.boundaries().cyclicBoundaryModels()[boundaryI]->controlMol(p, td);
    }
}

inline void controlPatchBoundaryThreadSafe
(
    Foam::dsmcCloud& cloud,
    const Foam::label boundaryI,
    Foam::dsmcParcel& p,
    Foam::dsmcParcel::trackingData& td
)
{
    Foam::dsmcPatchBoundary& model =
        *cloud.boundaries().patchBoundaryModels()[boundaryI];
    const Foam::word& modelType = model.type();
    const bool threadSafePatchModel =
        modelType == "dsmcDiffuseWallPatch"
     || modelType == "dsmcSpecularWallPatch";

    if (useOpenMPMoveCriticals(cloud))
    {
        if (threadSafePatchModel)
        {
            model.controlParticle(p, td);
        }
        else
        {
            #pragma omp critical(dsmcMoveBoundary)
            {
                model.controlParticle(p, td);
            }
        }
    }
    else
    {
        model.controlParticle(p, td);
    }
}
}

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
                if (cloud.trackerActive())
                {
                    trackParcelFaceTransitionThreadSafe(cloud, *this);
                }

                forAll(cloud.boundaries().cyclicBoundaryModels(), c)
                {
                    const labelList& faces = cloud.boundaries().cyclicBoundaryModels()[c]->allFaces();

                    if (Foam::hyCompat::indexOf(faces, face()) != -1)
                    {
                        controlCyclicBoundaryThreadSafe(cloud, c, *this, td);
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
                controlPatchBoundaryThreadSafe(cloud, c, *this, td);
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
            controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
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
            controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
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



