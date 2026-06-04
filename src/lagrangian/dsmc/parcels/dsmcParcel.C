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

#include <chrono>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace
{

using MoveProfileClock = std::chrono::steady_clock;

inline Foam::scalar elapsedSeconds
(
    const MoveProfileClock::time_point& start
)
{
    return std::chrono::duration<Foam::scalar>
    (
        MoveProfileClock::now() - start
    ).count();
}

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
    const bool profileMoveDetail = cloud.profilingDetailEnabled();

    if (cloud.replicatedMeshActive() && cell() >= 0
        && cell() < cloud.mesh().nCells())
    {
        localCellI_ = cloud.replicatedMesh().localMesh().toLocal(cell());
    }

    if (!profileMoveDetail)
    {
        if (isFree())
        {
            if (newParcel_ != -1)
            {
                stepFraction() = td.moveRng.sample01();
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

        return td.keepParticle;
    }

    if (isFree())
    {
        if (newParcel_ != -1)
        {
            stepFraction() = td.moveRng.sample01();
            newParcel_ = -1;
        }

        vector Utracking = U_;

        while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
        {
            Utracking = U_;
            meshTools::constrainDirection(mesh(), mesh().solutionD(), Utracking);

            const vector d = deviationFromMeshCentre();
            const scalar f = 1 - stepFraction();

            const label cellBefore = cell();
            const vector trackDisplacement = f*trackTime*Utracking - d;

            const auto tTrack0 = MoveProfileClock::now();
            trackToAndHitFace(trackDisplacement, f, cloud, td);
            td.moveTrackWallTime += elapsedSeconds(tTrack0);
            if
            (
                cellBefore >= 0
             && cellBefore < cloud.moveItersPerCell().size()
            )
            {
                ++cloud.moveItersPerCell()[cellBefore];
                ++cloud.moveItersPerCellCumulative()[cellBefore];
            }

            if (face() != -1)
            {
                ++td.moveFaceHitCount;

                if (cloud.trackerActive())
                {
                    const auto tTracker0 = MoveProfileClock::now();
                    trackParcelFaceTransitionThreadSafe(cloud, *this);
                    td.moveTrackerWallTime += elapsedSeconds(tTracker0);
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
                        ++td.moveCyclicHitCount;
                        const auto tBoundary0 = MoveProfileClock::now();
                        controlCyclicBoundaryThreadSafe
                        (
                            cloud,
                            cyclicModelId,
                            *this,
                            td
                        );
                        td.moveBoundaryWallTime += elapsedSeconds(tBoundary0);
                    }
                }
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
                ++td.moveStuckHitCount;
                const auto tBoundary0 = MoveProfileClock::now();
                controlPatchBoundaryThreadSafe
                (
                    cloud,
                    patchModelId,
                    *this,
                    td
                );
                td.moveBoundaryWallTime += elapsedSeconds(tBoundary0);
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
            const bool profileMoveDetail = cloud.profilingDetailEnabled();
            if (!profileMoveDetail)
            {
                controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
                return false;
            }

            ++td.movePatchHitCount;
            const auto tBoundary0 = MoveProfileClock::now();
            controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
            td.moveBoundaryWallTime += elapsedSeconds(tBoundary0);
        }
    }

    return false;
}

void Foam::dsmcParcel::hitProcessorPatch(dsmcCloud& cloud, trackingData& td)
{
    if (cloud.profilingDetailEnabled())
    {
        ++td.moveProcessorHitCount;
    }
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
            const bool profileMoveDetail = cloud.profilingDetailEnabled();
            if (!profileMoveDetail)
            {
                controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
                return;
            }

            ++td.movePatchHitCount;
            const auto tBoundary0 = MoveProfileClock::now();
            controlPatchBoundaryThreadSafe(cloud, patchModelId, *this, td);
            td.moveBoundaryWallTime += elapsedSeconds(tBoundary0);
            return;
        }
    }

    const bool profileMoveDetail = cloud.profilingDetailEnabled();
    if (!profileMoveDetail)
    {
        cloud.handleWallInteraction(*this, td);
        return;
    }

    const auto tBoundary0 = MoveProfileClock::now();
    cloud.handleWallInteraction(*this, td);
    td.moveBoundaryWallTime += elapsedSeconds(tBoundary0);
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

// Fast binary write for migration — bypasses the virtual operator<< chain
// for the dsmcParcel-specific data.  Base particle data uses the standard
// particle::operator<< which writes raw bytes in binary mode.
// List/Field types are written as raw size+data to avoid operator<< overhead.
void Foam::dsmcParcel::writeBinaryFast(Ostream& os) const
{
    // Base particle data: use existing operator<< (binary = raw write)
    os << static_cast<const particle&>(*this);

    // dsmcParcel fixed fields — raw writes, same format as operator<< binary path
    os.write(reinterpret_cast<const char*>(&U_), sizeof(vector));
    os.write(reinterpret_cast<const char*>(&RWF_), sizeof(scalar));
    os.write(reinterpret_cast<const char*>(&ERot_), sizeof(scalar));
    os.write(reinterpret_cast<const char*>(&ELevel_), sizeof(label));
    os.write(reinterpret_cast<const char*>(&typeId_), sizeof(label));
    os.write(reinterpret_cast<const char*>(&newParcel_), sizeof(label));
    os.write(reinterpret_cast<const char*>(&classification_), sizeof(label));

    // vibLevel_ — use operator<< to ensure token/format compatibility
    os << vibLevel_;

    // Stuck data — use operator<< for correct binary format
    os << label(stuck_ != nullptr ? 1 : 0);
    if (stuck_)
    {
        os << stuck_->wallTemperature() << stuck_->wallVectors();
    }

    os.check("dsmcParcel::writeBinaryFast");
}


void Foam::dsmcParcel::packTransfer(TransferData& td) const
{
    const barycentric& c = coordinates();
    td.coordinates[0] = c[0];
    td.coordinates[1] = c[1];
    td.coordinates[2] = c[2];
    td.coordinates[3] = c[3];
    td.celli = cell();
    td.tetFacei = tetFace();
    td.tetPti = tetPt();
    td.facei = face();
    td.stepFraction = stepFraction();
    td.origProc = origProc();
    td.origId = origId();

    td.U[0] = U_.x();
    td.U[1] = U_.y();
    td.U[2] = U_.z();
    td.RWF = RWF_;
    td.ERot = ERot_;
    td.ELevel = ELevel_;
    td.typeId = typeId_;
    td.newParcel = newParcel_;
    td.classification = classification_;

    td.nVibModes = min(vibLevel_.size(), label(maxVibModes));
    for (label i = 0; i < td.nVibModes; ++i)
        td.vibLevel[i] = vibLevel_[i];
    for (label i = td.nVibModes; i < maxVibModes; ++i)
        td.vibLevel[i] = 0;
}


Foam::dsmcParcel* Foam::dsmcParcel::unpackTransfer
(
    const polyMesh& mesh,
    const TransferData& td
)
{
    labelList vib(td.nVibModes);
    for (label i = 0; i < td.nVibModes; ++i)
        vib[i] = td.vibLevel[i];

    auto* p = new dsmcParcel
    (
        mesh,
        barycentric(td.coordinates[0], td.coordinates[1],
                    td.coordinates[2], td.coordinates[3]),
        td.celli,
        td.tetFacei,
        td.tetPti,
        vector(td.U[0], td.U[1], td.U[2]),
        td.RWF,
        td.ERot,
        td.ELevel,
        td.typeId,
        td.newParcel,
        td.classification,
        vib
    );

    p->face() = td.facei;
    p->stepFraction() = td.stepFraction;
    p->origProc() = td.origProc;
    p->origId() = td.origId;

    return p;
}


#include "dsmcParcelIO.C"

// ************************************************************************* //

Foam::label Foam::dsmcParcel::TrackedParcel::nDELETED = 0;
