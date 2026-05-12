/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017, 2020 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "Cloud.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamBuffers.H"
#include "mapPolyMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include <chrono>
#include <unordered_set>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace Foam
{
namespace cloudOpenMP
{
template<class TrackCloudType>
inline auto moveEnabled(const TrackCloudType& cloud, int)
-> decltype(cloud.openmpMoveEnabled(), bool())
{
    return cloud.openmpMoveEnabled();
}

template<class TrackCloudType>
inline bool moveEnabled(const TrackCloudType&, long)
{
    return false;
}

template<class TrackCloudType>
inline auto moveThreads(const TrackCloudType& cloud, int)
-> decltype(cloud.ompNumThreads(), label())
{
    return cloud.ompNumThreads();
}

template<class TrackCloudType>
inline label moveThreads(const TrackCloudType&, long)
{
    return 1;
}

template<class TrackCloudType>
inline auto moveSchedule(const TrackCloudType& cloud, int)
-> decltype(cloud.openmpMoveSchedule(), word())
{
    return cloud.openmpMoveSchedule();
}

template<class TrackCloudType>
inline word moveSchedule(const TrackCloudType&, long)
{
    return "static";
}

template<class TrackCloudType>
inline auto moveChunk(const TrackCloudType& cloud, int)
-> decltype(cloud.openmpMoveChunk(), label())
{
    return cloud.openmpMoveChunk();
}

template<class TrackCloudType>
inline label moveChunk(const TrackCloudType&, long)
{
    return 64;
}

template<class TrackCloudType>
inline auto moveStageProbeEnabled(const TrackCloudType& cloud, int)
-> decltype(cloud.moveStageProbeEnabled(), bool())
{
    return cloud.moveStageProbeEnabled();
}

template<class TrackCloudType>
inline bool moveStageProbeEnabled(const TrackCloudType&, long)
{
    return false;
}

template<class TrackCloudType>
inline auto hasParticlePartition(const TrackCloudType& cloud, int)
-> decltype
(
    cloud.particleLoadStart(),
    cloud.particleLoadEnd(),
    cloud.cellOccupancy(),
    bool()
)
{
    return true;
}

template<class TrackCloudType>
inline bool hasParticlePartition(const TrackCloudType&, long)
{
    return false;
}

template<class TrackCloudType>
inline auto particleLoadStart(const TrackCloudType& cloud, int)
-> decltype(cloud.particleLoadStart())
{
    return cloud.particleLoadStart();
}

template<class TrackCloudType>
inline auto particleLoadEnd(const TrackCloudType& cloud, int)
-> decltype(cloud.particleLoadEnd())
{
    return cloud.particleLoadEnd();
}

template<class TrackCloudType>
inline auto cellOccupancy(const TrackCloudType& cloud, int)
-> decltype(cloud.cellOccupancy())
{
    return cloud.cellOccupancy();
}

template<class TrackCloudType, class ParticleType>
inline auto storeMoveOrderedParcels
(
    TrackCloudType& cloud,
    const List<ParticleType*>& parcels,
    const labelList& threadOffsets,
    int
)
-> decltype(cloud.storeMoveOrderedParcels(parcels, threadOffsets), void())
{
    cloud.storeMoveOrderedParcels(parcels, threadOffsets);
}

template<class TrackCloudType, class ParticleType>
inline void storeMoveOrderedParcels
(
    TrackCloudType&,
    const List<ParticleType*>&,
    const labelList&,
    long
)
{}

template<class TrackCloudType>
inline auto clearMoveOrderedParcels(TrackCloudType& cloud, int)
-> decltype(cloud.clearMoveOrderedParcels(), void())
{
    cloud.clearMoveOrderedParcels();
}

template<class TrackCloudType>
inline void clearMoveOrderedParcels(TrackCloudType&, long)
{}

template<class TrackCloudType>
inline auto hasMoveOrderedParcels(const TrackCloudType& cloud, int)
-> decltype(cloud.hasMoveOrderedParcels(), bool())
{
    return cloud.hasMoveOrderedParcels();
}

template<class TrackCloudType>
inline bool hasMoveOrderedParcels(const TrackCloudType&, long)
{
    return false;
}

template<class TrackCloudType, class ParticleType>
inline auto moveOrderedParcels(const TrackCloudType& cloud, int)
-> decltype(cloud.moveOrderedParcels())
{
    return cloud.moveOrderedParcels();
}

template<class TrackCloudType>
inline auto moveOrderedThreadOffsets(const TrackCloudType& cloud, int)
-> decltype(cloud.moveOrderedThreadOffsets())
{
    return cloud.moveOrderedThreadOffsets();
}

template<class TrackCloudType>
inline labelList moveOrderedThreadOffsets(const TrackCloudType&, long)
{
    return labelList();
}

template<class TrackCloudType, class ParticleType>
inline List<ParticleType*> moveOrderedParcels(const TrackCloudType&, long)
{
    return List<ParticleType*>();
}

template<class TrackCloudType>
inline auto beginMoveAppendCapture(TrackCloudType& cloud, int)
-> decltype(cloud.beginMoveAppendCapture(), void())
{
    cloud.beginMoveAppendCapture();
}

template<class TrackCloudType>
inline void beginMoveAppendCapture(TrackCloudType&, long)
{}

template<class TrackCloudType>
inline auto endMoveAppendCapture(TrackCloudType& cloud, int)
-> decltype(cloud.endMoveAppendCapture(), void())
{
    cloud.endMoveAppendCapture();
}

template<class TrackCloudType>
inline void endMoveAppendCapture(TrackCloudType&, long)
{}

template<class TrackCloudType>
inline auto beginMoveDeferredAppendStage(TrackCloudType& cloud, int)
-> decltype(cloud.beginMoveDeferredAppendStage(), void())
{
    cloud.beginMoveDeferredAppendStage();
}

template<class TrackCloudType>
inline void beginMoveDeferredAppendStage(TrackCloudType&, long)
{}

template<class TrackCloudType, class ParticleType>
inline auto pendingMoveParcels(const TrackCloudType& cloud, int)
-> decltype(cloud.pendingMoveParcels())
{
    return cloud.pendingMoveParcels();
}

template<class TrackCloudType, class ParticleType>
inline DynamicList<ParticleType*> pendingMoveParcels(const TrackCloudType&, long)
{
    return DynamicList<ParticleType*>();
}

template<class TrackCloudType, class ParticleType>
inline auto moveAppendedParcels(const TrackCloudType& cloud, int)
-> decltype(cloud.moveAppendedParcels())
{
    return cloud.moveAppendedParcels();
}

template<class TrackCloudType, class ParticleType>
inline DynamicList<ParticleType*> moveAppendedParcels(const TrackCloudType&, long)
{
    return DynamicList<ParticleType*>();
}

template<class TrackCloudType>
inline auto clearPendingMoveParcels(TrackCloudType& cloud, int)
-> decltype(cloud.clearPendingMoveParcels(), void())
{
    cloud.clearPendingMoveParcels();
}

template<class TrackCloudType>
inline void clearPendingMoveParcels(TrackCloudType&, long)
{}

template<class TrackCloudType>
inline auto clearMoveAppendedParcels(TrackCloudType& cloud, int)
-> decltype(cloud.clearMoveAppendedParcels(), void())
{
    cloud.clearMoveAppendedParcels();
}

template<class TrackCloudType>
inline void clearMoveAppendedParcels(TrackCloudType&, long)
{}

template<class TrackCloudType, class ParticleType>
inline auto recordMoveAppendedParcel(TrackCloudType& cloud, ParticleType* pPtr, int)
-> decltype(cloud.recordMoveAppendedParcel(pPtr), void())
{
    cloud.recordMoveAppendedParcel(pPtr);
}

template<class TrackCloudType, class ParticleType>
inline void recordMoveAppendedParcel(TrackCloudType&, ParticleType*, long)
{}

template<class TrackCloudType, class ParticleType>
inline auto recordPendingMoveParcel(TrackCloudType& cloud, ParticleType* pPtr, int)
-> decltype(cloud.recordPendingMoveParcel(pPtr), void())
{
    cloud.recordPendingMoveParcel(pPtr);
}

template<class TrackCloudType, class ParticleType>
inline void recordPendingMoveParcel(TrackCloudType&, ParticleType*, long)
{}

template<class TrackCloudType>
inline auto recordMoveThreadCounts
(
    TrackCloudType& cloud,
    const labelList& counts,
    int
)
-> decltype(cloud.recordMoveThreadCounts(counts), void())
{
    cloud.recordMoveThreadCounts(counts);
}

template<class TrackCloudType>
inline void recordMoveThreadCounts(TrackCloudType&, const labelList&, long)
{}

template<class TrackCloudType>
inline auto recordMoveLoopPasses
(
    TrackCloudType& cloud,
    const label nPasses,
    int
)
-> decltype(cloud.recordMoveLoopPasses(nPasses), void())
{
    cloud.recordMoveLoopPasses(nPasses);
}

template<class TrackCloudType>
inline void recordMoveLoopPasses(TrackCloudType&, const label, long)
{}

template<class TrackCloudType>
inline auto recordMoveDeferredParcels
(
    TrackCloudType& cloud,
    const label nDeferred,
    int
)
-> decltype(cloud.recordMoveDeferredParcels(nDeferred), void())
{
    cloud.recordMoveDeferredParcels(nDeferred);
}

template<class TrackCloudType>
inline void recordMoveDeferredParcels(TrackCloudType&, const label, long)
{}

template<class TrackCloudType>
inline auto recordMoveReceivedParcels
(
    TrackCloudType& cloud,
    const label nReceived,
    int
)
-> decltype(cloud.recordMoveReceivedParcels(nReceived), void())
{
    cloud.recordMoveReceivedParcels(nReceived);
}

template<class TrackCloudType>
inline void recordMoveReceivedParcels(TrackCloudType&, const label, long)
{}

template<class TrackCloudType>
inline auto recordMoveExtractDetail
(
    TrackCloudType& cloud,
    const scalar deferredWallTime,
    const scalar orderedReuseWallTime,
    const scalar fullScanWallTime,
    const label deferredParcels,
    const label orderedReuseParcels,
    const label fullScanParcels,
    int
)
-> decltype
(
    cloud.recordMoveExtractDetail
    (
        deferredWallTime,
        orderedReuseWallTime,
        fullScanWallTime,
        deferredParcels,
        orderedReuseParcels,
        fullScanParcels
    ),
    void()
)
{
    cloud.recordMoveExtractDetail
    (
        deferredWallTime,
        orderedReuseWallTime,
        fullScanWallTime,
        deferredParcels,
        orderedReuseParcels,
        fullScanParcels
    );
}

template<class TrackCloudType>
inline void recordMoveExtractDetail
(
    TrackCloudType&,
    const scalar,
    const scalar,
    const scalar,
    const label,
    const label,
    const label,
    long
)
{}

template<class TrackCloudType>
inline auto recordMoveFirstPassReuseCheck
(
    TrackCloudType& cloud,
    const bool hasOrdered,
    const bool offsetsOk,
    const label currentSize,
    const label priorSize,
    const label pendingSize,
    int
)
-> decltype
(
    cloud.recordMoveFirstPassReuseCheck
    (
        hasOrdered,
        offsetsOk,
        currentSize,
        priorSize,
        pendingSize
    ),
    void()
)
{
    cloud.recordMoveFirstPassReuseCheck
    (
        hasOrdered,
        offsetsOk,
        currentSize,
        priorSize,
        pendingSize
    );
}

template<class TrackCloudType>
inline void recordMoveFirstPassReuseCheck
(
    TrackCloudType&,
    const bool,
    const bool,
    const label,
    const label,
    const label,
    long
)
{}

template<class TrackCloudType>
inline auto recordMoveThreadProfile
(
    TrackCloudType& cloud,
    const labelList& counts,
    const scalarField& wallTimes,
    int
)
-> decltype(cloud.recordMoveThreadProfile(counts, wallTimes), void())
{
    cloud.recordMoveThreadProfile(counts, wallTimes);
}

template<class TrackCloudType>
inline void recordMoveThreadProfile
(
    TrackCloudType& cloud,
    const labelList& counts,
    const scalarField&,
    long
)
{
    cloudOpenMP::recordMoveThreadCounts(cloud, counts, 0);
}

template<class TrackCloudType>
inline auto recordMoveInnerProfile
(
    TrackCloudType& cloud,
    const scalarField& trackWallTimes,
    const scalarField& trackerWallTimes,
    const scalarField& boundaryWallTimes,
    const labelList& faceHitCounts,
    const labelList& cyclicHitCounts,
    const labelList& stuckHitCounts,
    const labelList& patchHitCounts,
    const labelList& processorHitCounts,
    int
)
-> decltype
(
    cloud.recordMoveInnerProfile
    (
        trackWallTimes,
        trackerWallTimes,
        boundaryWallTimes,
        faceHitCounts,
        cyclicHitCounts,
        stuckHitCounts,
        patchHitCounts,
        processorHitCounts
    ),
    void()
)
{
    cloud.recordMoveInnerProfile
    (
        trackWallTimes,
        trackerWallTimes,
        boundaryWallTimes,
        faceHitCounts,
        cyclicHitCounts,
        stuckHitCounts,
        patchHitCounts,
        processorHitCounts
    );
}

template<class TrackCloudType>
inline void recordMoveInnerProfile
(
    TrackCloudType&,
    const scalarField&,
    const scalarField&,
    const scalarField&,
    const labelList&,
    const labelList&,
    const labelList&,
    const labelList&,
    const labelList&,
    long
)
{}

template<class TrackCloudType>
inline auto recordMovePhaseProfile
(
    TrackCloudType& cloud,
    const scalar preControlWallTime,
    const scalar resetSetupWallTime,
    const scalar extractWallTime,
    const scalar kernelWallTime,
    const scalar commitWallTime,
    const scalar transferFinalizeWallTime,
    int
)
-> decltype
(
    cloud.recordMovePhaseProfile
    (
        preControlWallTime,
        resetSetupWallTime,
        extractWallTime,
        kernelWallTime,
        commitWallTime,
        transferFinalizeWallTime
    ),
    void()
)
{
    cloud.recordMovePhaseProfile
    (
        preControlWallTime,
        resetSetupWallTime,
        extractWallTime,
        kernelWallTime,
        commitWallTime,
        transferFinalizeWallTime
    );
}

template<class TrackCloudType>
inline void recordMovePhaseProfile
(
    TrackCloudType&,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    long
)
{}

template<class TrackCloudType>
inline auto recordMoveCommitCounts
(
    TrackCloudType& cloud,
    const label extracted,
    const label survivors,
    const label transferred,
    const label deleted,
    int
)
-> decltype(cloud.recordMoveCommitCounts(extracted, survivors, transferred, deleted), void())
{
    cloud.recordMoveCommitCounts(extracted, survivors, transferred, deleted);
}

template<class TrackCloudType>
inline void recordMoveCommitCounts
(
    TrackCloudType&,
    const label,
    const label,
    const label,
    const label,
    long
)
{}
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::checkPatches() const
{
    for (const polyPatch& pp : polyMesh_.boundaryMesh())
    {
        const auto* camipp = isA<cyclicAMIPolyPatch>(pp);

        if (camipp && camipp->owner() && camipp->AMI().distributed())
        {
            FatalErrorInFunction
                << "Particle tracking across AMI patches is only currently "
                << "supported for cases where the AMI patches reside on a "
                << "single processor" << abort(FatalError);
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const Foam::zero,
    const word& cloudName
)
:
    cloud(pMesh, cloudName),
    polyMesh_(pMesh),
    geometryType_(cloud::geometryType::COORDINATES)
{
    checkPatches();

    (void)polyMesh_.solutionD();
    (void)polyMesh_.tetBasePtIs();
    (void)polyMesh_.oldCellCentres();
}


template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles
)
:
    Cloud<ParticleType>(pMesh, Foam::zero{}, cloudName)
{
    if (particles.size())
    {
        IDLList<ParticleType>::operator=(particles);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::addParticle(ParticleType* pPtr)
{
    this->append(pPtr);
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteParticle(ParticleType& p)
{
    delete(this->remove(&p));
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteLostParticles()
{
    for (ParticleType& p : *this)
    {
        if (p.cell() == -1)
        {
            WarningInFunction
                << "deleting lost particle at position " << p.position()
                << endl;

            deleteParticle(p);
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::cloudReset(const Cloud<ParticleType>& c)
{
    // Reset particle count and particles only
    // - not changing the cloud object registry or reference to the polyMesh
    ParticleType::particleCount_ = 0;
    IDLList<ParticleType>::operator=(c);
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::Cloud<ParticleType>::move
(
    TrackCloudType& cloud,
    typename ParticleType::trackingData& td,
    const scalar trackTime
)
{
    using clock_type = std::chrono::steady_clock;
    const auto tResetSetupBegin = clock_type::now();

    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();
    const globalMeshData& pData = polyMesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Indexing of equivalent patch on neighbour processor into the
    // procPatches list on the neighbour
    const labelList& procPatchNeighbours = pData.processorPatchNeighbours();

    // Which processors this processor is connected to
    const labelList& neighbourProcs = pData.topology().procNeighbours();

    // Clear the global positions as these are about to change
    globalPositionsPtr_.clear();


    // For v2112 and earlier: pre-assembled lists of particles
    // to be transferred and target patch on a per processor basis.
    // Apart from memory overhead of assembling the lists this adds
    // allocations/de-allocation when building linked-lists.

    // Now stream particle transfer tuples directly into PstreamBuffers.
    // Use a local cache of UOPstream wrappers for the formatters
    // (since there are potentially many particles being shifted about).


    // Allocate transfer buffers,
    // automatic clearStorage when UIPstream closes is disabled.
    PstreamBuffers pBufs;
    pBufs.allowClearRecv(false);

    // Cache of opened UOPstream wrappers
    PtrList<UOPstream> UOPstreamPtrs(Pstream::nProcs());

    #ifdef _OPENMP
    const bool useOpenMPMove =
        cloudOpenMP::moveEnabled(cloud, 0)
     && cloudOpenMP::moveThreads(cloud, 0) > 1
     && this->size() > 1;
    const label moveThreads = cloudOpenMP::moveThreads(cloud, 0);
    const bool moveStageProbe = cloudOpenMP::moveStageProbeEnabled(cloud, 0);
    const word moveSchedule = cloudOpenMP::moveSchedule(cloud, 0);
    const label moveChunk = max(cloudOpenMP::moveChunk(cloud, 0), label(1));
    const bool moveOrderedReuse =
        cloud.mesh().time().controlDict().template lookupOrDefault<bool>
        (
            "openmpMoveOrderedReuse",
            true
        );
    const bool useParticlePartition =
        useOpenMPMove && cloudOpenMP::hasParticlePartition(cloud, 0);
    const bool useMoveParticlePartition =
        useParticlePartition && !Pstream::parRun();
    #else
    const bool useOpenMPMove = false;
    const label moveThreads = 1;
    const bool moveStageProbe = false;
    const word moveSchedule = "static";
    const label moveChunk = 64;
    const bool moveOrderedReuse = false;
    const bool useParticlePartition = false;
    const bool useMoveParticlePartition = false;
    #endif

    scalar moveExtractWallTime = 0.0;
    scalar moveKernelWallTime = 0.0;
    scalar moveCommitWallTime = 0.0;
    scalar moveTransferFinalizeWallTime = 0.0;
    scalar moveResetSetupWallTime = 0.0;
    bool resetPending = true;
    label moveLoopPasses = 0;
    const bool accumulateMixedMoveOrdered =
        moveOrderedReuse
     && useOpenMPMove
     && Pstream::parRun()
     && !useMoveParticlePartition;
    DynamicList<ParticleType*> accumulatedMoveOrdered;

    if (accumulateMixedMoveOrdered)
    {
        accumulatedMoveOrdered.setCapacity(this->size());
    }

    if (useOpenMPMove && Pstream::parRun())
    {
        // Pre-build lazy particle-tracking mesh data before entering the
        // threaded move kernel on decomposed meshes.
        (void)polyMesh_.tetBasePtIs();
        (void)polyMesh_.oldCellCentres();
    }

    moveResetSetupWallTime +=
        std::chrono::duration<scalar>(clock_type::now() - tResetSetupBegin).count();

    // While there are particles to transfer
    while (true)
    {
        ++moveLoopPasses;
        const auto tLoopSetupBegin = clock_type::now();

        if (moveOrderedReuse && useOpenMPMove && !useMoveParticlePartition)
        {
            cloudOpenMP::beginMoveAppendCapture(cloud, 0);
        }

        // Reset transfer buffers
        pBufs.clear();

        // Rewind existing streams
        forAll(UOPstreamPtrs, proci)
        {
            auto* osptr = UOPstreamPtrs.get(proci);
            if (osptr)
            {
                osptr->rewind();
            }
        }

        moveResetSetupWallTime +=
            std::chrono::duration<scalar>(clock_type::now() - tLoopSetupBegin).count();

        if (useOpenMPMove)
        {
            List<ParticleType*> particles;
            labelList threadOffsets(moveThreads + 1, 0);
            const auto rawDeferredParcels =
                cloudOpenMP::pendingMoveParcels<TrackCloudType, ParticleType>(cloud, 0);
            DynamicList<ParticleType*> filteredDeferredParcels;
            const DynamicList<ParticleType*>* deferredParcelsPtr = &rawDeferredParcels;
            const auto appendedParcels =
                cloudOpenMP::moveAppendedParcels<TrackCloudType, ParticleType>(cloud, 0);

            if (rawDeferredParcels.size() > 1)
            {
                std::unordered_set<ParticleType*> seenDeferred;
                seenDeferred.reserve(rawDeferredParcels.size()*2);
                filteredDeferredParcels.setCapacity(rawDeferredParcels.size());

                forAll(rawDeferredParcels, i)
                {
                    ParticleType* pPtr = rawDeferredParcels[i];

                    if (pPtr && seenDeferred.insert(pPtr).second)
                    {
                        filteredDeferredParcels.append(pPtr);
                    }
                }

                if (filteredDeferredParcels.size() != rawDeferredParcels.size())
                {
                    deferredParcelsPtr = &filteredDeferredParcels;
                }
            }

            const DynamicList<ParticleType*>& deferredParcels = *deferredParcelsPtr;
            cloudOpenMP::recordMoveDeferredParcels(cloud, deferredParcels.size(), 0);
            const auto tExtractBegin = clock_type::now();
            scalar deferredExtractWallTime = 0.0;
            scalar orderedReuseExtractWallTime = 0.0;
            scalar fullScanExtractWallTime = 0.0;
            label deferredExtractParcels = 0;
            label orderedReuseExtractParcels = 0;
            label fullScanExtractParcels = 0;

            if (useMoveParticlePartition)
            {
                const auto& particleLoadStart = cloudOpenMP::particleLoadStart(cloud, 0);
                const auto& particleLoadEnd = cloudOpenMP::particleLoadEnd(cloud, 0);
                labelList threadCounts(moveThreads, 0);
                labelList occupancyThreadCounts(moveThreads, 0);
                labelList deferredThreadCounts(moveThreads, 0);

                auto ownerThreadForCell =
                [&](const label celli, const label fallbackIndex)
                {
                    if (celli >= 0)
                    {
                        for (label threadI = 0; threadI < moveThreads; ++threadI)
                        {
                            if
                            (
                                celli >= particleLoadStart[threadI]
                             && celli < particleLoadEnd[threadI]
                            )
                            {
                                return threadI;
                            }
                        }
                    }

                    return min(moveThreads - 1, fallbackIndex*moveThreads/max(deferredParcels.size(), label(1)));
                };

                #ifdef _OPENMP
                #pragma omp parallel for num_threads(moveThreads) schedule(static)
                #endif
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    label localCount = 0;

                    for
                    (
                        label celli = particleLoadStart[threadI];
                        celli < particleLoadEnd[threadI];
                        ++celli
                    )
                    {
                        localCount += cloud.occupancyCount(celli);
                    }

                    threadCounts[threadI] = localCount;
                    occupancyThreadCounts[threadI] = localCount;
                }

                forAll(deferredParcels, i)
                {
                    const label threadI = ownerThreadForCell(deferredParcels[i]->cell(), i);
                    ++deferredThreadCounts[threadI];
                    ++threadCounts[threadI];
                }

                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    threadOffsets[threadI + 1] =
                        threadOffsets[threadI] + threadCounts[threadI];
                }

                particles.setSize(threadOffsets.last());

                #ifdef _OPENMP
                #pragma omp parallel for num_threads(moveThreads) schedule(static)
                #endif
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    label particlei = threadOffsets[threadI];

                    for
                    (
                        label celli = particleLoadStart[threadI];
                        celli < particleLoadEnd[threadI];
                        ++celli
                    )
                    {
                        for
                        (
                            label occI = cloud.occupancyStart(celli);
                            occI < cloud.occupancyEnd(celli);
                            ++occI
                        )
                        {
                            particles[particlei++] = cloud.occupancyParcel(occI);
                        }
                    }
                }

                labelList deferredWriteOffsets(moveThreads, 0);

                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    deferredWriteOffsets[threadI] =
                        threadOffsets[threadI] + occupancyThreadCounts[threadI];
                }

                forAll(deferredParcels, i)
                {
                    const label threadI = ownerThreadForCell(deferredParcels[i]->cell(), i);
                    particles[deferredWriteOffsets[threadI]++] = deferredParcels[i];
                }

            }
            else
            {
                const bool canReuseMoveOrdered =
                    [&]()
                    {
                        const bool firstPass =
                            accumulateMixedMoveOrdered && moveLoopPasses == 1;

                        if (!firstPass)
                        {
                            return false;
                        }

                        const bool hasOrdered =
                            cloudOpenMP::hasMoveOrderedParcels(cloud, 0);
                        const label priorSize =
                            hasOrdered
                          ? cloudOpenMP::moveOrderedParcels<TrackCloudType, ParticleType>(cloud, 0).size()
                          : 0;
                        const label appendedSize = appendedParcels.size();
                        const bool offsetsOk =
                            cloudOpenMP::moveOrderedThreadOffsets(cloud, 0).size()
                         == moveThreads + 1;

                        cloudOpenMP::recordMoveFirstPassReuseCheck
                        (
                            cloud,
                            hasOrdered,
                            offsetsOk,
                            this->size(),
                            priorSize,
                            appendedSize,
                            0
                        );

                        return
                            hasOrdered
                         && offsetsOk
                         && priorSize + appendedSize == this->size();
                    }();

                if (accumulateMixedMoveOrdered && moveLoopPasses > 1)
                {
                    const auto tBranchBegin = clock_type::now();
                    particles.setSize(deferredParcels.size());

                    forAll(deferredParcels, i)
                    {
                        particles[i] = deferredParcels[i];
                    }

                    deferredExtractWallTime =
                        std::chrono::duration<scalar>(clock_type::now() - tBranchBegin).count();
                    deferredExtractParcels = particles.size();
                }
                else if (canReuseMoveOrdered)
                {
                    const auto tBranchBegin = clock_type::now();
                    const auto& priorMoveOrdered =
                        cloudOpenMP::moveOrderedParcels<TrackCloudType, ParticleType>(cloud, 0);
                    particles.setSize(priorMoveOrdered.size() + appendedParcels.size());

                    forAll(priorMoveOrdered, i)
                    {
                        particles[i] = priorMoveOrdered[i];
                    }

                    forAll(appendedParcels, i)
                    {
                        particles[priorMoveOrdered.size() + i] = appendedParcels[i];
                    }

                    for (label threadI = 0; threadI < moveThreads; ++threadI)
                    {
                        threadOffsets[threadI] = threadI*particles.size()/moveThreads;
                    }

                    threadOffsets[moveThreads] = particles.size();
                    orderedReuseExtractWallTime =
                        std::chrono::duration<scalar>(clock_type::now() - tBranchBegin).count();
                    orderedReuseExtractParcels = particles.size();
                }
                else
                {
                    const auto tBranchBegin = clock_type::now();
                    particles.setSize(this->size());
                    label particlei = 0;

                    for (ParticleType& p : *this)
                    {
                        particles[particlei++] = &p;
                    }

                    fullScanExtractWallTime =
                        std::chrono::duration<scalar>(clock_type::now() - tBranchBegin).count();
                    fullScanExtractParcels = particles.size();
                }

                if (!canReuseMoveOrdered)
                {
                    for (label threadI = 0; threadI < moveThreads; ++threadI)
                    {
                        threadOffsets[threadI] = threadI*particles.size()/moveThreads;
                    }

                    threadOffsets[moveThreads] = particles.size();
                }

                if (accumulateMixedMoveOrdered && moveLoopPasses == 1)
                {
                    cloudOpenMP::beginMoveDeferredAppendStage(cloud, 0);
                }
            }

            cloudOpenMP::clearPendingMoveParcels(cloud, 0);

            moveExtractWallTime +=
                std::chrono::duration<scalar>(clock_type::now() - tExtractBegin).count();

            if (moveStageProbe)
            {
                label invalidStartCount = 0;
                label sampleInvalidCell = -1;

                forAll(particles, i)
                {
                    const label celli = particles[i]->cell();

                    if (celli < 0 || celli >= polyMesh_.nCells())
                    {
                        ++invalidStartCount;

                        if (sampleInvalidCell == -1)
                        {
                            sampleInvalidCell = celli;
                        }
                    }
                }

                Pout<< "Move stage rank " << Pstream::myProcNo()
                    << " pass " << moveLoopPasses
                    << ": after extract particles=" << particles.size()
                    << " deferred=" << deferredParcels.size()
                    << " appended=" << appendedParcels.size()
                    << " invalidStart=" << invalidStartCount;

                if (sampleInvalidCell != -1)
                {
                    Pout<< " sampleInvalidCell=" << sampleInvalidCell;
                }

                Pout<< nl << endl;
            }

            cloudOpenMP::recordMoveExtractDetail
            (
                cloud,
                deferredExtractWallTime,
                orderedReuseExtractWallTime,
                fullScanExtractWallTime,
                deferredExtractParcels,
                orderedReuseExtractParcels,
                fullScanExtractParcels,
                0
            );

            const bool inlineReset = resetPending;
            resetPending = false;

            List<label> keepParticleFlags(particles.size(), 1);
            List<label> switchProcessorFlags(particles.size(), 0);
            labelList processedParticleCounts(moveThreads, 0);
            scalarField moveThreadWallTimes(moveThreads, 0.0);
            scalarField moveThreadTrackWallTimes(moveThreads, 0.0);
            scalarField moveThreadTrackerWallTimes(moveThreads, 0.0);
            scalarField moveThreadBoundaryWallTimes(moveThreads, 0.0);
            labelList moveThreadFaceHitCounts(moveThreads, 0);
            labelList moveThreadCyclicHitCounts(moveThreads, 0);
            labelList moveThreadStuckHitCounts(moveThreads, 0);
            labelList moveThreadPatchHitCounts(moveThreads, 0);
            labelList moveThreadProcessorHitCounts(moveThreads, 0);

            #ifdef _OPENMP
            omp_sched_t sched = omp_sched_static;

            if (moveSchedule == "dynamic")
            {
                sched = omp_sched_dynamic;
            }
            else if (moveSchedule == "guided")
            {
                sched = omp_sched_guided;
            }

            omp_set_schedule(sched, int(moveChunk));
            #endif

            const auto tKernelBegin = clock_type::now();
            if (moveStageProbe)
            {
                Pout<< "Move stage rank " << Pstream::myProcNo()
                    << " pass " << moveLoopPasses
                    << ": before move kernel particles=" << particles.size()
                    << nl << endl;
            }

            #pragma omp parallel num_threads(moveThreads)
            {
                typename ParticleType::trackingData localTd(cloud);
                const label threadI =
                #ifdef _OPENMP
                    omp_get_thread_num();
                #else
                    0;
                #endif

                label localProcessedCount = 0;
                const auto tBegin = clock_type::now();

                if (useMoveParticlePartition)
                {
                    for (label i = threadOffsets[threadI]; i < threadOffsets[threadI + 1]; ++i)
                    {
                        ParticleType& p = *particles[i];

                        if (inlineReset)
                        {
                            p.reset();
                        }

                        localTd.switchProcessor = false;
                        localTd.keepParticle = true;

                        keepParticleFlags[i] = p.move(cloud, localTd, trackTime) ? 1 : 0;
                        switchProcessorFlags[i] = localTd.switchProcessor ? 1 : 0;
                        ++localProcessedCount;
                    }
                }
                else
                {
                    #pragma omp for schedule(runtime)
                    for (label i = 0; i < particles.size(); ++i)
                    {
                        ParticleType& p = *particles[i];

                        if (moveStageProbe && (i % 25000 == 0))
                        {
                            #pragma omp critical(dsmcMoveProbe)
                            {
                                Pout<< "Move stage rank " << Pstream::myProcNo()
                                    << " pass " << moveLoopPasses
                                    << " thread " << threadI
                                    << ": moving i=" << i
                                    << " cell=" << p.cell()
                                    << " face=" << p.face()
                                    << " stepFraction=" << p.stepFraction()
                                    << " orig=" << p.origProc()
                                    << ':' << p.origId()
                                    << nl << endl;
                            }
                        }

                        if (inlineReset)
                        {
                            p.reset();
                        }

                        localTd.switchProcessor = false;
                        localTd.keepParticle = true;

                        keepParticleFlags[i] = p.move(cloud, localTd, trackTime) ? 1 : 0;
                        switchProcessorFlags[i] = localTd.switchProcessor ? 1 : 0;
                        ++localProcessedCount;
                    }
                }

                processedParticleCounts[threadI] = localProcessedCount;
                moveThreadWallTimes[threadI] =
                    std::chrono::duration<scalar>(clock_type::now() - tBegin).count();
                moveThreadTrackWallTimes[threadI] = localTd.moveTrackWallTime;
                moveThreadTrackerWallTimes[threadI] = localTd.moveTrackerWallTime;
                moveThreadBoundaryWallTimes[threadI] = localTd.moveBoundaryWallTime;
                moveThreadFaceHitCounts[threadI] = localTd.moveFaceHitCount;
                moveThreadCyclicHitCounts[threadI] = localTd.moveCyclicHitCount;
                moveThreadStuckHitCounts[threadI] = localTd.moveStuckHitCount;
                moveThreadPatchHitCounts[threadI] = localTd.movePatchHitCount;
                moveThreadProcessorHitCounts[threadI] = localTd.moveProcessorHitCount;
            }
            moveKernelWallTime +=
                std::chrono::duration<scalar>(clock_type::now() - tKernelBegin).count();

            if (moveStageProbe)
            {
                label processedTotal = 0;
                label processorHitTotal = 0;
                label invalidKeepTotal = 0;
                label sampleInvalidKeepCell = -1;

                forAll(processedParticleCounts, threadI)
                {
                    processedTotal += processedParticleCounts[threadI];
                    processorHitTotal += moveThreadProcessorHitCounts[threadI];
                }

                forAll(particles, i)
                {
                    if (!keepParticleFlags[i])
                    {
                        continue;
                    }

                    const label celli = particles[i]->cell();

                    if (celli < 0 || celli >= polyMesh_.nCells())
                    {
                        ++invalidKeepTotal;

                        if (sampleInvalidKeepCell == -1)
                        {
                            sampleInvalidKeepCell = celli;
                        }
                    }
                }

                Pout<< "Move stage rank " << Pstream::myProcNo()
                    << " pass " << moveLoopPasses
                    << ": after move kernel processed=" << processedTotal
                    << " processorHits=" << processorHitTotal
                    << " invalidKeep=" << invalidKeepTotal;

                if (sampleInvalidKeepCell != -1)
                {
                    Pout<< " sampleInvalidKeepCell=" << sampleInvalidKeepCell;
                }

                Pout<< nl << endl;
            }

            cloudOpenMP::recordMoveThreadProfile
            (
                cloud,
                processedParticleCounts,
                moveThreadWallTimes,
                0
            );
            cloudOpenMP::recordMoveInnerProfile
            (
                cloud,
                moveThreadTrackWallTimes,
                moveThreadTrackerWallTimes,
                moveThreadBoundaryWallTimes,
                moveThreadFaceHitCounts,
                moveThreadCyclicHitCounts,
                moveThreadStuckHitCounts,
                moveThreadPatchHitCounts,
                moveThreadProcessorHitCounts,
                0
            );

            const auto tCommitBegin = clock_type::now();
            labelList commitStart(moveThreads + 1, 0);

            if (useMoveParticlePartition)
            {
                commitStart = threadOffsets;
            }
            else
            {
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    commitStart[threadI] = threadI*particles.size()/moveThreads;
                }

                commitStart[moveThreads] = particles.size();
            }

            labelList survivingCounts(moveThreads, 0);
            labelList transferCounts(moveThreads, 0);
            labelList deleteCounts(moveThreads, 0);

            #ifdef _OPENMP
            #pragma omp parallel for num_threads(moveThreads) schedule(static)
            #endif
            for (label threadI = 0; threadI < moveThreads; ++threadI)
            {
                label localSurvivors = 0;
                label localTransfers = 0;
                label localDeletes = 0;

                for
                (
                    label i = commitStart[threadI];
                    i < commitStart[threadI + 1];
                    ++i
                )
                {
                    if (keepParticleFlags[i])
                    {
                        if (switchProcessorFlags[i])
                        {
                            ++localTransfers;
                        }
                        else
                        {
                            ++localSurvivors;
                        }
                    }
                    else
                    {
                        ++localDeletes;
                    }
                }

                survivingCounts[threadI] = localSurvivors;

                transferCounts[threadI] = localTransfers;
                deleteCounts[threadI] = localDeletes;
            }

            labelList survivingOffsets(moveThreads + 1, 0);
            labelList transferOffsets(moveThreads + 1, 0);
            labelList deleteOffsets(moveThreads + 1, 0);

            for (label threadI = 0; threadI < moveThreads; ++threadI)
            {
                survivingOffsets[threadI + 1] =
                    survivingOffsets[threadI] + survivingCounts[threadI];
                transferOffsets[threadI + 1] =
                    transferOffsets[threadI] + transferCounts[threadI];
                deleteOffsets[threadI + 1] =
                    deleteOffsets[threadI] + deleteCounts[threadI];
            }

            List<ParticleType*> survivingList(survivingOffsets.last());
            List<ParticleType*> transferList(transferOffsets.last());
            List<ParticleType*> deleteList(deleteOffsets.last());

            #ifdef _OPENMP
            #pragma omp parallel for num_threads(moveThreads) schedule(static)
            #endif
            for (label threadI = 0; threadI < moveThreads; ++threadI)
            {
                label survivingI = survivingOffsets[threadI];
                label transferI = transferOffsets[threadI];
                label deleteI = deleteOffsets[threadI];

                for
                (
                    label i = commitStart[threadI];
                    i < commitStart[threadI + 1];
                    ++i
                )
                {
                    if (keepParticleFlags[i])
                    {
                        if (switchProcessorFlags[i])
                        {
                            transferList[transferI++] = particles[i];
                        }
                        else
                        {
                            survivingList[survivingI++] = particles[i];
                        }
                    }
                    else
                    {
                        deleteList[deleteI++] = particles[i];
                    }
                }
            }

            forAll(transferList, i)
            {
                ParticleType& p = *transferList[i];

                #ifdef FULLDEBUG
                if
                (
                    !Pstream::parRun()
                 || !p.onBoundaryFace()
                 || procPatchNeighbours[p.patch()] < 0
                )
                {
                    FatalErrorInFunction
                        << "Switch processor flag is true when no parallel "
                        << "transfer is possible. This is a bug."
                        << exit(FatalError);
                }
                #endif

                const label patchi = p.patch();

                const label toProci =
                (
                    refCast<const processorPolyPatch>(pbm[patchi])
                    .neighbProcNo()
                );

                auto* osptr = UOPstreamPtrs.get(toProci);
                if (!osptr)
                {
                    osptr = new UOPstream(toProci, pBufs);
                    UOPstreamPtrs.set(toProci, osptr);
                }

                p.prepareForParallelTransfer();
                (*osptr) << procPatchNeighbours[patchi] << p;
                deleteParticle(p);
            }

            forAll(deleteList, i)
            {
                deleteParticle(*deleteList[i]);
            }

            if (accumulateMixedMoveOrdered)
            {
                forAll(survivingList, i)
                {
                    accumulatedMoveOrdered.append(survivingList[i]);
                }
            }
            else
            {
                cloudOpenMP::storeMoveOrderedParcels
                (
                    cloud,
                    survivingList,
                    survivingOffsets,
                    0
                );
            }

            cloudOpenMP::recordMoveCommitCounts
            (
                cloud,
                particles.size(),
                survivingList.size(),
                transferList.size(),
                deleteList.size(),
                0
            );
            moveCommitWallTime +=
                std::chrono::duration<scalar>(clock_type::now() - tCommitBegin).count();

            if (moveStageProbe)
            {
                label invalidSurvivorCount = 0;
                label invalidTransferCount = 0;
                label sampleInvalidSurvivorCell = -1;
                label sampleInvalidTransferCell = -1;

                forAll(survivingList, i)
                {
                    const label celli = survivingList[i]->cell();

                    if (celli < 0 || celli >= polyMesh_.nCells())
                    {
                        ++invalidSurvivorCount;

                        if (sampleInvalidSurvivorCell == -1)
                        {
                            sampleInvalidSurvivorCell = celli;
                        }
                    }
                }

                forAll(transferList, i)
                {
                    const label celli = transferList[i]->cell();

                    if (celli < 0 || celli >= polyMesh_.nCells())
                    {
                        ++invalidTransferCount;

                        if (sampleInvalidTransferCell == -1)
                        {
                            sampleInvalidTransferCell = celli;
                        }
                    }
                }

                Pout<< "Move stage rank " << Pstream::myProcNo()
                    << " pass " << moveLoopPasses
                    << ": after commit survivors=" << survivingList.size()
                    << " transfers=" << transferList.size()
                    << " deletes=" << deleteList.size()
                    << " invalidSurvivors=" << invalidSurvivorCount
                    << " invalidTransfers=" << invalidTransferCount;

                if (sampleInvalidSurvivorCell != -1)
                {
                    Pout<< " sampleInvalidSurvivorCell="
                        << sampleInvalidSurvivorCell;
                }

                if (sampleInvalidTransferCell != -1)
                {
                    Pout<< " sampleInvalidTransferCell="
                        << sampleInvalidTransferCell;
                }

                Pout<< nl << endl;
            }
        }
        else
        {
            const bool inlineReset = resetPending;
            resetPending = false;

            const auto tBegin = clock_type::now();
            label processedParticleCount = 0;
            List<ParticleType*> particles(this->size());
            label particlei = 0;

            for (ParticleType& p : *this)
            {
                particles[particlei++] = &p;
            }

            forAll(particles, i)
            {
                ParticleType& p = *particles[i];

                if (inlineReset)
                {
                    p.reset();
                }

                // Move the particle
                bool keepParticle = p.move(cloud, td, trackTime);
                ++processedParticleCount;

                // If the particle is to be kept
                // (i.e. it hasn't passed through an inlet or outlet)
                if (keepParticle)
                {
                    // If the particle is going to switch processors, stream it
                    // into transfer buffers
                    if (td.switchProcessor)
                    {
                        #ifdef FULLDEBUG
                        if
                        (
                            !Pstream::parRun()
                         || !p.onBoundaryFace()
                         || procPatchNeighbours[p.patch()] < 0
                        )
                        {
                            FatalErrorInFunction
                                << "Switch processor flag is true when no parallel "
                                << "transfer is possible. This is a bug."
                                << exit(FatalError);
                        }
                        #endif

                        const label patchi = p.patch();

                        const label toProci =
                        (
                            refCast<const processorPolyPatch>(pbm[patchi])
                            .neighbProcNo()
                        );

                        // Get/create output stream
                        auto* osptr = UOPstreamPtrs.get(toProci);
                        if (!osptr)
                        {
                            osptr = new UOPstream(toProci, pBufs);
                            UOPstreamPtrs.set(toProci, osptr);
                        }

                        p.prepareForParallelTransfer();

                        // Tuple: (patchi particle)
                        (*osptr) << procPatchNeighbours[patchi] << p;

                        // Can now remove from my list
                        deleteParticle(p);
                    }
                }
                else
                {
                    deleteParticle(p);
                }
            }

            labelList processedCounts(1, processedParticleCount);
            scalarField moveWallTimes
            (
                1,
                std::chrono::duration<scalar>(clock_type::now() - tBegin).count()
            );
            moveKernelWallTime += moveWallTimes[0];
            cloudOpenMP::recordMoveThreadProfile(cloud, processedCounts, moveWallTimes, 0);
        }

        if (!Pstream::parRun())
        {
            break;
        }

        const auto tTransferBegin = clock_type::now();
        if (moveStageProbe)
        {
            Pout<< "Move stage rank " << Pstream::myProcNo()
                << " pass " << moveLoopPasses
                << ": before finishedNeighbourSends" << nl << endl;
        }
        pBufs.finishedNeighbourSends(neighbourProcs);
        if (moveStageProbe)
        {
            Pout<< "Move stage rank " << Pstream::myProcNo()
                << " pass " << moveLoopPasses
                << ": after finishedNeighbourSends" << nl << endl;
        }

        if (moveStageProbe)
        {
            Pout<< "Move stage rank " << Pstream::myProcNo()
                << " pass " << moveLoopPasses
                << ": before hasRecvData reduce" << nl << endl;
        }
        if (!returnReduceOr(pBufs.hasRecvData()))
        {
            // No parcels to transfer
            moveTransferFinalizeWallTime +=
                std::chrono::duration<scalar>(clock_type::now() - tTransferBegin).count();
            break;
        }
        if (moveStageProbe)
        {
            Pout<< "Move stage rank " << Pstream::myProcNo()
                << " pass " << moveLoopPasses
                << ": after hasRecvData reduce" << nl << endl;
        }

        // Retrieve from receive buffers
        label receivedCount = 0;
        label invalidReceivedCount = 0;
        label sampleInvalidReceivedCell = -1;
        if (moveStageProbe)
        {
            Pout<< "Move stage rank " << Pstream::myProcNo()
                << " pass " << moveLoopPasses
                << ": before recv loop" << nl << endl;
        }
        for (const label proci : neighbourProcs)
        {
            if (pBufs.recvDataCount(proci))
            {
                UIPstream is(proci, pBufs);

                // Read out each (patchi particle) tuple
                while (!is.eof())
                {
                    label patchi = pTraits<label>(is);
                    auto* newp = new ParticleType(polyMesh_, is);

                    // The real patch index
                    patchi = procPatches[patchi];

                    (*newp).correctAfterParallelTransfer(patchi, td);

                    const label correctedCell = newp->cell();

                    if (correctedCell < 0 || correctedCell >= polyMesh_.nCells())
                    {
                        ++invalidReceivedCount;

                        if (sampleInvalidReceivedCell == -1)
                        {
                            sampleInvalidReceivedCell = correctedCell;
                        }
                    }

                    addParticle(newp);
                    cloudOpenMP::recordPendingMoveParcel(cloud, newp, 0);
                    ++receivedCount;
                }
            }
        }
        if (moveStageProbe)
        {
            Pout<< "Move stage rank " << Pstream::myProcNo()
                << " pass " << moveLoopPasses
                << ": after recv loop received=" << receivedCount
                << " invalidReceived=" << invalidReceivedCount;

            if (sampleInvalidReceivedCell != -1)
            {
                Pout<< " sampleInvalidReceivedCell=" << sampleInvalidReceivedCell;
            }

            Pout<< nl << endl;
        }

        cloudOpenMP::recordMoveReceivedParcels(cloud, receivedCount, 0);

        moveTransferFinalizeWallTime +=
            std::chrono::duration<scalar>(clock_type::now() - tTransferBegin).count();

    }

    if (accumulateMixedMoveOrdered)
    {
        labelList finalThreadOffsets(moveThreads + 1, 0);

        for (label threadI = 0; threadI < moveThreads; ++threadI)
        {
            finalThreadOffsets[threadI] =
                threadI*accumulatedMoveOrdered.size()/moveThreads;
        }

        finalThreadOffsets[moveThreads] = accumulatedMoveOrdered.size();

        List<ParticleType*> finalMoveOrdered(accumulatedMoveOrdered.size());

        forAll(accumulatedMoveOrdered, i)
        {
            finalMoveOrdered[i] = accumulatedMoveOrdered[i];
        }

        cloudOpenMP::storeMoveOrderedParcels
        (
            cloud,
            finalMoveOrdered,
            finalThreadOffsets,
            0
        );
    }

    cloudOpenMP::recordMoveLoopPasses(cloud, moveLoopPasses, 0);
    cloudOpenMP::recordMovePhaseProfile
    (
        cloud,
        0.0,
        moveResetSetupWallTime,
        moveExtractWallTime,
        moveKernelWallTime,
        moveCommitWallTime,
        moveTransferFinalizeWallTime,
        0
    );
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::autoMap(const mapPolyMesh& mapper)
{
    if (!globalPositionsPtr_)
    {
        FatalErrorInFunction
            << "Global positions are not available. "
            << "Cloud::storeGlobalPositions has not been called."
            << exit(FatalError);
    }

    // Reset stored data that relies on the mesh
    //    polyMesh_.clearCellTree();
    cellWallFacesPtr_.clear();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    (void)polyMesh_.tetBasePtIs();
    (void)polyMesh_.oldCellCentres();

    const vectorField& positions = globalPositionsPtr_();

    label i = 0;
    for (ParticleType& p : *this)
    {
        p.autoMap(positions[i], mapper);
        ++i;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writePositions() const
{
    OFstream os
    (
        this->db().time().path()/this->name() + "_positions.obj"
    );

    for (const ParticleType& p : *this)
    {
        const point position(p.position());
        os  << "v "
            << position.x() << ' '
            << position.y() << ' '
            << position.z() << nl;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::storeGlobalPositions() const
{
    // Store the global positions for later use by autoMap. It would be
    // preferable not to need this. If the mapPolyMesh object passed to autoMap
    // had a copy of the old mesh then the global positions could be recovered
    // within autoMap, and this pre-processing would not be necessary.

    globalPositionsPtr_.reset(new vectorField(this->size()));
    vectorField& positions = globalPositionsPtr_();

    label i = 0;
    for (const ParticleType& p : *this)
    {
        positions[i] = p.position();
        ++i;
    }
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudIO.C"

// ************************************************************************* //
