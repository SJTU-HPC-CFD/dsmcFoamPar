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
#include "Pstream.H"

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <chrono>

namespace
{
using clock_type = std::chrono::steady_clock;

inline bool useOpenMPMoveCriticals(const Foam::dsmcCloud& cloud)
{
    #ifdef _OPENMP
    return cloud.openmpMoveEnabled() && omp_in_parallel();
    #else
    return false;
    #endif
}

inline bool useOpenMPMoveTrackCriticals
(
    const Foam::dsmcCloud& cloud,
    const Foam::label celli
)
{
    #ifdef _OPENMP
    return
        cloud.openmpMoveEnabled()
     && omp_in_parallel()
     && Foam::Pstream::parRun()
     && cloud.openmpMoveGuardCell(celli);
    #else
    (void)cloud;
    (void)celli;
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
    const bool recordMoveDetail = cloud.profilingDetailEnabled();

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
            auto moveTrackStep = [&]()
            {
                Utracking = U_;
                meshTools::constrainDirection(mesh(), mesh().solutionD(), Utracking);

                const vector d = deviationFromMeshCentre();
                const scalar f = 1 - stepFraction();
                if (recordMoveDetail)
                {
                    const auto tTrackBegin = clock_type::now();
                    trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);
                    td.moveTrackWallTime +=
                        std::chrono::duration<scalar>(clock_type::now() - tTrackBegin).count();
                }
                else
                {
                    trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);
                }

                if (face() != -1)
                {
                    if (recordMoveDetail)
                    {
                        ++td.moveFaceHitCount;
                    }

                    if (cloud.trackerActive())
                    {
                        if (recordMoveDetail)
                        {
                            const auto tTrackerBegin = clock_type::now();
                            trackParcelFaceTransitionThreadSafe(cloud, *this);
                            td.moveTrackerWallTime +=
                                std::chrono::duration<scalar>(clock_type::now() - tTrackerBegin).count();
                        }
                        else
                        {
                            trackParcelFaceTransitionThreadSafe(cloud, *this);
                        }
                    }

                    const label patchIndex = patch();

                    if
                    (
                        patchIndex >= 0
                     && patchIndex
                      < cloud.boundaries().cyclicBoundaryToModelIds().size()
                    )
                    {
                        const label cyclicModelId =
                            cloud.boundaries().cyclicBoundaryToModelIds()[patchIndex];

                        if (cyclicModelId >= 0)
                        {
                            if (recordMoveDetail)
                            {
                                const auto tBoundaryBegin = clock_type::now();
                                controlCyclicBoundaryThreadSafe
                                (
                                    cloud,
                                    cyclicModelId,
                                    *this,
                                    td
                                );
                                td.moveBoundaryWallTime +=
                                    std::chrono::duration<scalar>(clock_type::now() - tBoundaryBegin).count();
                                ++td.moveCyclicHitCount;
                            }
                            else
                            {
                                controlCyclicBoundaryThreadSafe
                                (
                                    cloud,
                                    cyclicModelId,
                                    *this,
                                    td
                                );
                            }
                        }
                    }
                }
            };

            if (useOpenMPMoveTrackCriticals(cloud, cell()))
            {
                #pragma omp critical(dsmcMoveTrack)
                {
                    moveTrackStep();
                }
            }
            else
            {
                moveTrackStep();
            }
        }
    }
    else
    {
        const label patchIndex = stuck().wallTemperature()[1];

        if
        (
            patchIndex >= 0
         && patchIndex < cloud.boundaries().patchToModelIds().size()
        )
        {
            const label patchModelId =
                cloud.boundaries().patchToModelIds()[patchIndex];

            if (patchModelId >= 0)
            {
                if (recordMoveDetail)
                {
                    const auto tBoundaryBegin = clock_type::now();
                    controlPatchBoundaryThreadSafe
                    (
                        cloud,
                        patchModelId,
                        *this,
                        td
                    );
                    td.moveBoundaryWallTime +=
                        std::chrono::duration<scalar>(clock_type::now() - tBoundaryBegin).count();
                    ++td.moveStuckHitCount;
                }
                else
                {
                    controlPatchBoundaryThreadSafe
                    (
                        cloud,
                        patchModelId,
                        *this,
                        td
                    );
                }
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
            if (cloud.profilingDetailEnabled())
            {
                const auto tBoundaryBegin = clock_type::now();
                controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
                td.moveBoundaryWallTime +=
                    std::chrono::duration<scalar>(clock_type::now() - tBoundaryBegin).count();
                ++td.movePatchHitCount;
            }
            else
            {
                controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
            }
        }
    }

    return false;
}

void Foam::dsmcParcel::hitProcessorPatch(dsmcCloud& cloud, trackingData& td)
{
    td.switchProcessor = true;
    if (cloud.profilingDetailEnabled())
    {
        ++td.moveProcessorHitCount;
    }
}

void Foam::dsmcParcel::hitWallPatch(dsmcCloud& cloud, trackingData& td)
{
    const label patchIndex = patch();

    if (patchIndex >= 0 && patchIndex < cloud.boundaries().patchToModelIds().size())
    {
        const label patchModelId = cloud.boundaries().patchToModelIds()[patchIndex];

        if (patchModelId >= 0)
        {
            if (cloud.profilingDetailEnabled())
            {
                const auto tBoundaryBegin = clock_type::now();
                controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
                td.moveBoundaryWallTime +=
                    std::chrono::duration<scalar>(clock_type::now() - tBoundaryBegin).count();
                ++td.movePatchHitCount;
            }
            else
            {
                controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
            }
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



