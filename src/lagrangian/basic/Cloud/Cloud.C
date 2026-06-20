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

\*---------------------------------------------------------------------------*/

#include "Cloud.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "mapPolyMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"
#include "cyclicAMIPolyPatch.H"

#ifdef _OPENMP
    #include <omp.h>
#endif
#include <cstdint>
#include <type_traits>

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
inline auto moveOrderedReuseEnabled(const TrackCloudType& cloud, int)
-> decltype
(
    cloud.controlDict().template lookupOrDefault<bool>
    (
        "openmpMoveOrderedReuse",
        true
    ),
    bool()
)
{
    return cloud.controlDict().template lookupOrDefault<bool>
    (
        "openmpMoveOrderedReuse",
        true
    );
}

template<class TrackCloudType>
inline bool moveOrderedReuseEnabled(const TrackCloudType&, long)
{
    return false;
}

template<class TrackData>
inline auto setMoveSeed(TrackData& td, const label threadI, int)
-> decltype(td.moveRng = td.moveRng, void())
{
    typedef typename std::remove_reference
    <
        decltype(td.moveRng)
    >::type RngType;

    td.moveRng = RngType
    (
        uint64_t(threadI)*10000ULL + 42ULL
    );
}

template<class TrackData>
inline void setMoveSeed(TrackData&, const label, ...)
{}

template<class TrackCloudType>
inline auto replicatedMeshActive(const TrackCloudType& cloud, int)
-> decltype(cloud.replicatedMeshActive(), bool())
{
    return cloud.replicatedMeshActive();
}

template<class TrackCloudType>
inline bool replicatedMeshActive(const TrackCloudType&, long)
{
    return false;
}

template<class TrackData>
inline auto copyMoveDetailProfile(TrackData& dst, const TrackData& src, int)
-> decltype(dst.moveDetailProfile = src.moveDetailProfile, void())
{
    dst.moveDetailProfile = src.moveDetailProfile;
}

template<class TrackData>
inline void copyMoveDetailProfile(TrackData&, const TrackData&, long)
{}

template<class TrackData>
inline auto accumulateMoveDetailProfile(TrackData& dst, const TrackData& src, int)
-> decltype
(
    dst.moveParcels += src.moveParcels,
    dst.moveTrackCalls += src.moveTrackCalls,
    dst.moveSameTetNoFaceHits += src.moveSameTetNoFaceHits,
    dst.moveInternalTetNoFaceHits += src.moveInternalTetNoFaceHits,
    dst.moveFaceHits += src.moveFaceHits,
    dst.moveProcessorHits += src.moveProcessorHits,
    dst.moveCyclicHits += src.moveCyclicHits,
    dst.movePatchHits += src.movePatchHits,
    dst.moveStuckHits += src.moveStuckHits,
    dst.moveTrackWallTime += src.moveTrackWallTime,
    dst.moveTrackerWallTime += src.moveTrackerWallTime,
    dst.moveBoundaryWallTime += src.moveBoundaryWallTime,
    void()
)
{
    dst.moveParcels += src.moveParcels;
    dst.moveTrackCalls += src.moveTrackCalls;
    dst.moveSameTetNoFaceHits += src.moveSameTetNoFaceHits;
    dst.moveInternalTetNoFaceHits += src.moveInternalTetNoFaceHits;
    dst.moveFaceHits += src.moveFaceHits;
    dst.moveProcessorHits += src.moveProcessorHits;
    dst.moveCyclicHits += src.moveCyclicHits;
    dst.movePatchHits += src.movePatchHits;
    dst.moveStuckHits += src.moveStuckHits;
    dst.moveTrackWallTime += src.moveTrackWallTime;
    dst.moveTrackerWallTime += src.moveTrackerWallTime;
    dst.moveBoundaryWallTime += src.moveBoundaryWallTime;
}

template<class TrackData>
inline void accumulateMoveDetailProfile(TrackData&, const TrackData&, long)
{}

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

template<class TrackCloudType, class ParticleType>
inline auto transferMoveOrderedParcels
(
    TrackCloudType& cloud,
    List<ParticleType*>& parcels,
    const labelList& threadOffsets,
    int
)
-> decltype(cloud.transferMoveOrderedParcels(parcels, threadOffsets), void())
{
    cloud.transferMoveOrderedParcels(parcels, threadOffsets);
}

template<class TrackCloudType, class ParticleType>
inline void transferMoveOrderedParcels
(
    TrackCloudType&,
    List<ParticleType*>&,
    const labelList&,
    long
)
{}

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

template<class TrackCloudType, class ParticleType>
inline List<ParticleType*> moveOrderedParcels(const TrackCloudType&, long)
{
    return List<ParticleType*>();
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

template<class TrackCloudType, class ParticleListType>
inline auto appendMoveOrderedParcels
(
    TrackCloudType& cloud,
    const ParticleListType& parcels,
    int
)
-> decltype(cloud.appendBatchToMoveOrdered(parcels), void())
{
    cloud.appendBatchToMoveOrdered(parcels);
}

template<class TrackCloudType, class ParticleListType>
inline void appendMoveOrderedParcels
(
    TrackCloudType&,
    const ParticleListType&,
    long
)
{}
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::checkPatches() const
{
    const polyBoundaryMesh& pbm = polyMesh_.boundaryMesh();
    bool ok = true;
    forAll(pbm, patchI)
    {
        if (isA<cyclicAMIPolyPatch>(pbm[patchI]))
        {
            const cyclicAMIPolyPatch& cami =
                refCast<const cyclicAMIPolyPatch>(pbm[patchI]);

            if (cami.owner())
            {
                ok = ok && (cami.AMI().singlePatchProc() != -1);
            }
        }
    }

    if (!ok)
    {
        FatalErrorIn("void Foam::Cloud<ParticleType>::initCloud(const bool)")
            << "Particle tracking across AMI patches is only currently "
            << "supported for cases where the AMI patches reside on a "
            << "single processor" << abort(FatalError);
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::calcCellWallFaces() const
{
    cellWallFacesPtr_.reset(new PackedBoolList(pMesh().nCells(), false));

    PackedBoolList& cellWallFaces = cellWallFacesPtr_();

    const polyBoundaryMesh& patches = polyMesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<wallPolyPatch>(patches[patchI]))
        {
            const polyPatch& patch = patches[patchI];

            const labelList& pFaceCells = patch.faceCells();

            forAll(pFaceCells, pFCI)
            {
                cellWallFaces[pFaceCells[pFCI]] = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    labels_(),
    nTrackingRescues_(),
    cellWallFacesPtr_()
{
    checkPatches();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    IDLList<ParticleType>::operator=(particles);
}


template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh, cloudName),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    labels_(),
    nTrackingRescues_(),
    cellWallFacesPtr_()
{
    checkPatches();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    IDLList<ParticleType>::operator=(particles);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
const Foam::PackedBoolList& Foam::Cloud<ParticleType>::cellHasWallFaces()
const
{
    if (!cellWallFacesPtr_.valid())
    {
        calcCellWallFaces();
    }

    return cellWallFacesPtr_();
}


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
    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        ParticleType& p = pIter();
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
    // Reset particle cound and particles only
    // - not changing the cloud object registry or reference to the polyMesh
    ParticleType::particleCount_ = 0;
    IDLList<ParticleType>::operator=(c);
}


template<class ParticleType>
template<class TrackData>
void Foam::Cloud<ParticleType>::move(TrackData& td, const scalar trackTime)
{
    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();
    const bool useProcessorPatchTransfer =
        Pstream::parRun()
     && !cloudOpenMP::replicatedMeshActive(td.cloud(), 0);

    const labelList emptyLabelList;
    const labelList* procPatchesPtr = &emptyLabelList;
    const labelList* procPatchIndicesPtr = &emptyLabelList;
    const labelList* procPatchNeighboursPtr = &emptyLabelList;
    const labelList* neighbourProcsPtr = &emptyLabelList;

    if (useProcessorPatchTransfer)
    {
        const globalMeshData& pData = polyMesh_.globalData();

        procPatchesPtr = &pData.processorPatches();
        procPatchIndicesPtr = &pData.processorPatchIndices();
        procPatchNeighboursPtr = &pData.processorPatchNeighbours();
        neighbourProcsPtr = &pData[Pstream::myProcNo()];
    }

    // Which patches are processor patches
    const labelList& procPatches = *procPatchesPtr;

    // Indexing of patches into the procPatches list
    const labelList& procPatchIndices = *procPatchIndicesPtr;

    // Indexing of equivalent patch on neighbour processor into the
    // procPatches list on the neighbour
    const labelList& procPatchNeighbours = *procPatchNeighboursPtr;

    // Which processors this processor is connected to
    const labelList& neighbourProcs = *neighbourProcsPtr;

    // Indexing from the processor number into the neighbourProcs list
    labelList neighbourProcIndices(Pstream::nProcs(), -1);

    forAll(neighbourProcs, i)
    {
        neighbourProcIndices[neighbourProcs[i]] = i;
    }

    #ifdef _OPENMP
    const bool useOpenMPMove =
        cloudOpenMP::moveEnabled(td.cloud(), 0)
     && cloudOpenMP::moveThreads(td.cloud(), 0) > 1
     && this->size() > 1;
    #else
    const bool useOpenMPMove = false;
    #endif

    if (!useOpenMPMove)
    {
        // Initialise the stepFraction moved for the particles
        forAllIter(typename Cloud<ParticleType>, *this, pIter)
        {
            pIter().stepFraction() = 0;
        }
    }

    // Reset nTrackingRescues
    nTrackingRescues_ = 0;

    if (useOpenMPMove)
    {
        #ifdef _OPENMP
        // Build demand-driven mesh data before particle tracking enters
        // the OpenMP region.  Several tracking accessors are lazily
        // initialised and are not safe to construct concurrently.
        (void)polyMesh_.solutionD();
        (void)polyMesh_.tetBasePtIs();
        (void)polyMesh_.cells();
        (void)polyMesh_.oldPoints();
        (void)polyMesh_.cellVolumes();
        (void)polyMesh_.cellCentres();
        (void)polyMesh_.faceAreas();
        (void)polyMesh_.faceCentres();
        const polyBoundaryMesh& boundaryMesh = polyMesh_.boundaryMesh();
        (void)boundaryMesh.patchID();
        forAll(boundaryMesh, patchI)
        {
            (void)boundaryMesh[patchI].faceCells();
        }
        (void)this->cellHasWallFaces();

        const label moveThreads =
            max(label(1), cloudOpenMP::moveThreads(td.cloud(), 0));
        const word moveSchedule = cloudOpenMP::moveSchedule(td.cloud(), 0);
        const label moveChunk =
            max(label(1), cloudOpenMP::moveChunk(td.cloud(), 0));
        const bool moveOrderedReuse =
            cloudOpenMP::moveOrderedReuseEnabled(td.cloud(), 0);

        if (useProcessorPatchTransfer)
        {
            if (moveOrderedReuse)
            {
                cloudOpenMP::beginMoveAppendCapture(td.cloud(), 0);
            }

            List<IDLList<ParticleType> > particleTransferLists
            (
                neighbourProcs.size()
            );
            List<DynamicList<label> > patchIndexTransferLists
            (
                neighbourProcs.size()
            );
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
            bool resetPending = true;

            while (true)
            {
                particleTransferLists = IDLList<ParticleType>();

                forAll(patchIndexTransferLists, i)
                {
                    patchIndexTransferLists[i].clear();
                }

                List<ParticleType*> particles(this->size());
                label particlei = 0;

                forAllIter(typename Cloud<ParticleType>, *this, pIter)
                {
                    particles[particlei++] = &pIter();
                }

                particles.setSize(particlei);

                labelList threadOffsets(moveThreads + 1, 0);
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    threadOffsets[threadI] =
                        threadI*particles.size()/moveThreads;
                }

                threadOffsets[moveThreads] = particles.size();

                List<unsigned char> keepParticleFlags
                (
                    particles.size(),
                    static_cast<unsigned char>(1)
                );
                List<unsigned char> switchProcessorFlags
                (
                    particles.size(),
                    static_cast<unsigned char>(0)
                );
                const bool inlineReset = resetPending;
                resetPending = false;

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

                #pragma omp parallel num_threads(moveThreads)
                {
                    TrackData localTd(td.cloud());
                    cloudOpenMP::copyMoveDetailProfile(localTd, td, 0);
                    cloudOpenMP::setMoveSeed
                    (
                        localTd,
                        omp_get_thread_num(),
                        0
                    );

                    #pragma omp for schedule(runtime)
                    for (label i = 0; i < particles.size(); ++i)
                    {
                        if (inlineReset)
                        {
                            particles[i]->stepFraction() = 0;
                        }

                        localTd.switchProcessor = false;
                        localTd.keepParticle = true;
                        keepParticleFlags[i] =
                            particles[i]->move(localTd, trackTime) ? 1 : 0;
                        switchProcessorFlags[i] =
                            localTd.switchProcessor ? 1 : 0;
                    }

                    #pragma omp critical(dsmcMoveDetailProfile)
                    {
                        cloudOpenMP::accumulateMoveDetailProfile(td, localTd, 0);
                    }
                }

                labelList survivorCounts(moveThreads, 0);
                labelList transferCounts(moveThreads, 0);
                labelList deleteCounts(moveThreads, 0);

                #pragma omp parallel for num_threads(moveThreads) schedule(static)
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    label localSurvivors = 0;
                    label localTransfers = 0;
                    label localDeletes = 0;

                    for
                    (
                        label i = threadOffsets[threadI];
                        i < threadOffsets[threadI + 1];
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

                    survivorCounts[threadI] = localSurvivors;
                    transferCounts[threadI] = localTransfers;
                    deleteCounts[threadI] = localDeletes;
                }

                labelList survivorOffsets(moveThreads + 1, 0);
                labelList transferOffsets(moveThreads + 1, 0);
                labelList deleteOffsets(moveThreads + 1, 0);

                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    survivorOffsets[threadI + 1] =
                        survivorOffsets[threadI] + survivorCounts[threadI];
                    transferOffsets[threadI + 1] =
                        transferOffsets[threadI] + transferCounts[threadI];
                    deleteOffsets[threadI + 1] =
                        deleteOffsets[threadI] + deleteCounts[threadI];
                }

                List<ParticleType*> survivingParticles(survivorOffsets.last());
                List<ParticleType*> transferParticles(transferOffsets.last());
                List<ParticleType*> deletedParticles(deleteOffsets.last());

                #pragma omp parallel for num_threads(moveThreads) schedule(static)
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    label survivorI = survivorOffsets[threadI];
                    label transferI = transferOffsets[threadI];
                    label deleteI = deleteOffsets[threadI];

                    for
                    (
                        label i = threadOffsets[threadI];
                        i < threadOffsets[threadI + 1];
                        ++i
                    )
                    {
                        if (keepParticleFlags[i])
                        {
                            if (switchProcessorFlags[i])
                            {
                                transferParticles[transferI++] = particles[i];
                            }
                            else
                            {
                                survivingParticles[survivorI++] = particles[i];
                            }
                        }
                        else
                        {
                            deletedParticles[deleteI++] = particles[i];
                        }
                    }
                }

                forAll(transferParticles, i)
                {
                    ParticleType& p = *transferParticles[i];

                    if (p.face() >= pMesh().nInternalFaces())
                    {
                        const label patchI = pbm.whichPatch(p.face());

                        if (procPatchIndices[patchI] != -1)
                        {
                            const label n = neighbourProcIndices
                            [
                                refCast<const processorPolyPatch>
                                (
                                    pbm[patchI]
                                ).neighbProcNo()
                            ];

                            p.prepareForParallelTransfer(patchI, td);

                            particleTransferLists[n].append(this->remove(&p));

                            patchIndexTransferLists[n].append
                            (
                                procPatchNeighbours[patchI]
                            );
                        }
                    }
                }

                forAll(deletedParticles, i)
                {
                    deleteParticle(*deletedParticles[i]);
                }

                cloudOpenMP::storeMoveOrderedParcels
                (
                    td.cloud(),
                    survivingParticles,
                    survivorOffsets,
                    0
                );

                pBufs.clear();

                forAll(particleTransferLists, i)
                {
                    if (particleTransferLists[i].size())
                    {
                        UOPstream particleStream(neighbourProcs[i], pBufs);

                        particleStream
                            << patchIndexTransferLists[i]
                            << particleTransferLists[i];
                    }
                }

                labelList allNTrans(Pstream::nProcs());
                pBufs.finishedSends(allNTrans);

                bool transfered = false;

                forAll(allNTrans, i)
                {
                    if (allNTrans[i])
                    {
                        transfered = true;
                        break;
                    }
                }
                reduce(transfered, orOp<bool>());

                if (!transfered)
                {
                    break;
                }

                forAll(neighbourProcs, i)
                {
                    const label neighbProci = neighbourProcs[i];
                    const label nRec = allNTrans[neighbProci];

                    if (nRec)
                    {
                        UIPstream particleStream(neighbProci, pBufs);

                        labelList receivePatchIndex(particleStream);

                        IDLList<ParticleType> newParticles
                        (
                            particleStream,
                            typename ParticleType::iNew(polyMesh_)
                        );

                        label pI = 0;

                        forAllIter
                        (
                            typename Cloud<ParticleType>,
                            newParticles,
                            newpIter
                        )
                        {
                            ParticleType& newp = newpIter();

                            const label patchI =
                                procPatches[receivePatchIndex[pI++]];

                            newp.correctAfterParallelTransfer(patchI, td);

                            addParticle(newParticles.remove(&newp));
                        }
                    }
                }
            }

            if (cloud::debug)
            {
                reduce(nTrackingRescues_, sumOp<label>());

                if (nTrackingRescues_ > 0)
                {
                    Info<< nTrackingRescues_
                        << " tracking rescue corrections" << endl;
                }
            }

            return;
        }

        if (moveOrderedReuse)
        {
            cloudOpenMP::beginMoveAppendCapture(td.cloud(), 0);
        }

        List<ParticleType*> particlesStorage;
        const List<ParticleType*>* particlesPtr = nullptr;
        bool usingCloudOrdered = false;
        labelList threadOffsets(moveThreads + 1, 0);
        const auto& appendedParcels =
            cloudOpenMP::moveAppendedParcels
            <
                typename std::remove_reference<decltype(td.cloud())>::type,
                ParticleType
            >(td.cloud(), 0);

        const bool hasOrdered =
            cloudOpenMP::hasMoveOrderedParcels(td.cloud(), 0);
        const labelList orderedThreadOffsets =
            cloudOpenMP::moveOrderedThreadOffsets(td.cloud(), 0);
        const bool offsetsOk = orderedThreadOffsets.size() == moveThreads + 1;
        const label priorSize =
            hasOrdered
          ? cloudOpenMP::moveOrderedParcels
            <
                typename std::remove_reference<decltype(td.cloud())>::type,
                ParticleType
            >(td.cloud(), 0).size()
          : 0;

        const bool canReuseMoveOrdered =
            moveOrderedReuse
         && hasOrdered
         && offsetsOk
         && priorSize + appendedParcels.size() == this->size();

        if (canReuseMoveOrdered)
        {
            if (appendedParcels.size())
            {
                cloudOpenMP::appendMoveOrderedParcels
                (
                    td.cloud(),
                    appendedParcels,
                    0
                );
            }

            const auto& moveOrdered =
                cloudOpenMP::moveOrderedParcels
                <
                    typename std::remove_reference<decltype(td.cloud())>::type,
                    ParticleType
                >(td.cloud(), 0);

            particlesPtr = &moveOrdered;
            usingCloudOrdered = true;

            const labelList updatedThreadOffsets =
                cloudOpenMP::moveOrderedThreadOffsets(td.cloud(), 0);

            if (updatedThreadOffsets.size() == moveThreads + 1)
            {
                threadOffsets = updatedThreadOffsets;
            }
        }
        else
        {
            particlesStorage.setSize(this->size());
            label particlei = 0;

            forAllIter(typename Cloud<ParticleType>, *this, pIter)
            {
                particlesStorage[particlei++] = &pIter();
            }

            particlesStorage.setSize(particlei);
            particlesPtr = &particlesStorage;
        }

        const List<ParticleType*>& particles = *particlesPtr;

        if (threadOffsets.last() != particles.size())
        {
            for (label threadI = 0; threadI < moveThreads; ++threadI)
            {
                threadOffsets[threadI] = threadI*particles.size()/moveThreads;
            }

            threadOffsets[moveThreads] = particles.size();
        }

        List<unsigned char> keepParticleFlags
        (
            particles.size(),
            static_cast<unsigned char>(1)
        );
        const bool inlineReset = true;
        label deletedParticleCount = 0;

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

        #pragma omp parallel num_threads(moveThreads) reduction(+:deletedParticleCount)
        {
            TrackData localTd(td.cloud());
            cloudOpenMP::copyMoveDetailProfile(localTd, td, 0);
            cloudOpenMP::setMoveSeed(localTd, omp_get_thread_num(), 0);

            #pragma omp for schedule(runtime)
            for (label i = 0; i < particles.size(); ++i)
            {
                if (inlineReset)
                {
                    particles[i]->stepFraction() = 0;
                }

                localTd.switchProcessor = false;
                localTd.keepParticle = true;
                const bool keepParticle =
                    particles[i]->move(localTd, trackTime);
                keepParticleFlags[i] = keepParticle ? 1 : 0;

                if (!keepParticle)
                {
                    ++deletedParticleCount;
                }
            }

            #pragma omp critical(dsmcMoveDetailProfile)
            {
                cloudOpenMP::accumulateMoveDetailProfile(td, localTd, 0);
            }
        }

        if (deletedParticleCount == 0)
        {
            if (!usingCloudOrdered)
            {
                cloudOpenMP::storeMoveOrderedParcels
                (
                    td.cloud(),
                    particles,
                    threadOffsets,
                    0
                );
            }

            if (cloud::debug)
            {
                reduce(nTrackingRescues_, sumOp<label>());

                if (nTrackingRescues_ > 0)
                {
                    Info<< nTrackingRescues_
                        << " tracking rescue corrections" << endl;
                }
            }

            return;
        }

        labelList survivorCounts(moveThreads, 0);
        labelList deleteCounts(moveThreads, 0);

        #pragma omp parallel for num_threads(moveThreads) schedule(static)
        for (label threadI = 0; threadI < moveThreads; ++threadI)
        {
            label localSurvivors = 0;
            label localDeletes = 0;

            for
            (
                label i = threadOffsets[threadI];
                        i < threadOffsets[threadI + 1];
                        ++i
                    )
                    {
                        if (keepParticleFlags[i])
                        {
                            ++localSurvivors;
                        }
                else
                {
                    ++localDeletes;
                }
            }

            survivorCounts[threadI] = localSurvivors;
            deleteCounts[threadI] = localDeletes;
        }

        labelList survivorOffsets(moveThreads + 1, 0);
        labelList deleteOffsets(moveThreads + 1, 0);

        for (label threadI = 0; threadI < moveThreads; ++threadI)
        {
            survivorOffsets[threadI + 1] =
                survivorOffsets[threadI] + survivorCounts[threadI];
            deleteOffsets[threadI + 1] =
                deleteOffsets[threadI] + deleteCounts[threadI];
        }

        if (deleteOffsets.last() == 0)
        {
            if (!usingCloudOrdered)
            {
                cloudOpenMP::storeMoveOrderedParcels
                (
                    td.cloud(),
                    particles,
                    threadOffsets,
                    0
                );
            }

            if (cloud::debug)
            {
                reduce(nTrackingRescues_, sumOp<label>());

                if (nTrackingRescues_ > 0)
                {
                    Info<< nTrackingRescues_
                        << " tracking rescue corrections" << endl;
                }
            }

            return;
        }

        List<ParticleType*> survivingParticles(survivorOffsets.last());
        List<ParticleType*> deletedParticles(deleteOffsets.last());

        #pragma omp parallel for num_threads(moveThreads) schedule(static)
        for (label threadI = 0; threadI < moveThreads; ++threadI)
        {
            label survivorI = survivorOffsets[threadI];
            label deleteI = deleteOffsets[threadI];

            for
            (
                label i = threadOffsets[threadI];
                        i < threadOffsets[threadI + 1];
                        ++i
                    )
                    {
                        if (keepParticleFlags[i])
                        {
                            survivingParticles[survivorI++] = particles[i];
                        }
                else
                {
                    deletedParticles[deleteI++] = particles[i];
                }
            }
        }

        forAll(deletedParticles, i)
        {
            deleteParticle(*deletedParticles[i]);
        }

        cloudOpenMP::transferMoveOrderedParcels
        (
            td.cloud(),
            survivingParticles,
            survivorOffsets,
            0
        );
        #endif

        if (cloud::debug)
        {
            reduce(nTrackingRescues_, sumOp<label>());

            if (nTrackingRescues_ > 0)
            {
                Info<< nTrackingRescues_
                    << " tracking rescue corrections" << endl;
            }
        }

        return;
    }


    // List of lists of particles to be transfered for all of the
    // neighbour processors
    List<IDLList<ParticleType> > particleTransferLists
    (
        neighbourProcs.size()
    );

    // List of destination processorPatches indices for all of the
    // neighbour processors
    List<DynamicList<label> > patchIndexTransferLists
    (
        neighbourProcs.size()
    );

    // Allocate transfer buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // While there are particles to transfer
    while (true)
    {
        particleTransferLists = IDLList<ParticleType>();

        forAll(patchIndexTransferLists, i)
        {
            patchIndexTransferLists[i].clear();
        }

        // Loop over all particles
        forAllIter(typename Cloud<ParticleType>, *this, pIter)
        {
            ParticleType& p = pIter();

            // Move the particle
            bool keepParticle = p.move(td, trackTime);

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                // If we are running in parallel and the particle is on a
                // boundary face
                if
                (
                    useProcessorPatchTransfer
                 && p.face() >= pMesh().nInternalFaces()
                )
                {
                    label patchI = pbm.whichPatch(p.face());

                    // ... and the face is on a processor patch
                    // prepare it for transfer
                    if (procPatchIndices[patchI] != -1)
                    {
                        label n = neighbourProcIndices
                        [
                            refCast<const processorPolyPatch>
                            (
                                pbm[patchI]
                            ).neighbProcNo()
                        ];

                        p.prepareForParallelTransfer(patchI, td);

                        particleTransferLists[n].append(this->remove(&p));

                        patchIndexTransferLists[n].append
                        (
                            procPatchNeighbours[patchI]
                        );
                    }
                }
            }
            else
            {
                deleteParticle(p);
            }
        }

        if (!useProcessorPatchTransfer)
        {
            break;
        }


        // Clear transfer buffers
        pBufs.clear();

        // Stream into send buffers
        forAll(particleTransferLists, i)
        {
            if (particleTransferLists[i].size())
            {
                UOPstream particleStream
                (
                    neighbourProcs[i],
                    pBufs
                );

                particleStream
                    << patchIndexTransferLists[i]
                    << particleTransferLists[i];
            }
        }

        /*// Start sending. Sets number of bytes transferred
        labelListList allNTrans(Pstream::nProcs());
        pBufs.finishedSends(allNTrans);

        bool transfered = false;

        forAll(allNTrans, i)
        {
            forAll(allNTrans[i], j)
            {
                if (allNTrans[i][j])
                {
                    transfered = true;
                    break;
                }
            }
        }

        if (!transfered)
        {
            break;
        }

        // Retrieve from receive buffers
        forAll(neighbourProcs, i)
        {
            label neighbProci = neighbourProcs[i];

            label nRec = allNTrans[neighbProci][Pstream::myProcNo()];

            if (nRec)
            {
                UIPstream particleStream(neighbProci, pBufs);

                labelList receivePatchIndex(particleStream);

                IDLList<ParticleType> newParticles
                (
                    particleStream,
                    typename ParticleType::iNew(polyMesh_)
                );

                label pI = 0;

                forAllIter(typename Cloud<ParticleType>, newParticles, newpIter)
                {
                    ParticleType& newp = newpIter();

                    label patchI = procPatches[receivePatchIndex[pI++]];

                    newp.correctAfterParallelTransfer(patchI, td);

                    addParticle(newParticles.remove(&newp));
                }
            }
        }*/ // DELETED VINCENT

        // NEW VINCENT
        // Start sending. Sets number of bytes transferred
        labelList allNTrans(Pstream::nProcs());
        pBufs.finishedSends(allNTrans);


        bool transfered = false;

        forAll(allNTrans, i)
        {
            if (allNTrans[i])
            {
                transfered = true;
                break;
            }
        }
        reduce(transfered, orOp<bool>());

        if (!transfered)
        {
            break;
        }

        // Retrieve from receive buffers
        forAll(neighbourProcs, i)
        {
            label neighbProci = neighbourProcs[i];

            label nRec = allNTrans[neighbProci];

            if (nRec)
            {
                UIPstream particleStream(neighbProci, pBufs);

                labelList receivePatchIndex(particleStream);

                IDLList<ParticleType> newParticles
                (
                    particleStream,
                    typename ParticleType::iNew(polyMesh_)
                );

                label pI = 0;

                forAllIter(typename Cloud<ParticleType>, newParticles, newpIter)
                {
                    ParticleType& newp = newpIter();

                    label patchI = procPatches[receivePatchIndex[pI++]];

                    newp.correctAfterParallelTransfer(patchI, td);

                    addParticle(newParticles.remove(&newp));
                    //addParticle(new ParticleType(newp));
                    //delete(newParticles.remove(&newp));
                }
            }
        }
    }

    if (cloud::debug)
    {
        reduce(nTrackingRescues_, sumOp<label>());

        if (nTrackingRescues_ > 0)
        {
            Info<< nTrackingRescues_ << " tracking rescue corrections" << endl;
        }
    }
}


template<class ParticleType>
template<class TrackData>
void Foam::Cloud<ParticleType>::autoMap
(
    TrackData& td,
    const mapPolyMesh& mapper
)
{
    if (cloud::debug)
    {
        Info<< "Cloud<ParticleType>::autoMap(TrackData&, const mapPolyMesh&) "
            << "for lagrangian cloud " << cloud::name() << endl;
    }

    const labelList& reverseCellMap = mapper.reverseCellMap();
    const labelList& reverseFaceMap = mapper.reverseFaceMap();

    // Reset stored data that relies on the mesh
//    polyMesh_.clearCellTree();
    cellWallFacesPtr_.clear();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    scalar lostParticles = 0; // NEW VINCENT
    scalar totParticles = 0; // NEW VINCENT

    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        ParticleType& p = pIter();
        ++totParticles;

        if (reverseCellMap[p.cell()] >= 0)
        {
            p.cell() = reverseCellMap[p.cell()];

            if (p.face() >= 0 && reverseFaceMap[p.face()] >= 0)
            {
                p.face() = reverseFaceMap[p.face()];
            }
            else
            {
                p.face() = -1;
            }

            //p.initCellFacePt();
            // Mass will not be conserved but might be OK if sampling is not
            // yet enabled (because refinement steps are still in progress)
            if (p.initCellFacePtOrDeleteLostParticle()) // NEW VINCENT
            {
                ++lostParticles;
                deleteParticle(p);
            }
        }
        else
        {
            label trackStartCell = mapper.mergedCell(p.cell());

            if (trackStartCell < 0)
            {
                trackStartCell = 0;
                p.cell() = 0;
            }
            else
            {
                p.cell() = trackStartCell;
            }

            vector pos = p.position();

            const_cast<vector&>(p.position()) =
                polyMesh_.cellCentres()[trackStartCell];

            p.stepFraction() = 0;

            //p.initCellFacePt();
            // Mass will not be conserved but might be OK if sampling is not
            // yet enabled (because refinement steps are still in progress)
            if (p.initCellFacePtOrDeleteLostParticle()) // NEW VINCENT
            {
                ++lostParticles;
                deleteParticle(p);
                break;
            }

            p.track(pos, td);
        }
    }

    if (lostParticles > 0) // NEW VINCENT
    {
        Info<< "Lost particles deleted due to change in topology:" << tab
            << lostParticles << "/" << totParticles << endl;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/this->name() + "_positions.obj"
    );

    forAllConstIter(typename Cloud<ParticleType>, *this, pIter)
    {
        const ParticleType& p = pIter();
        pObj<< "v " << p.position().x() << " " << p.position().y() << " "
            << p.position().z() << nl;
    }

    pObj.flush();
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudIO.C"

// ************************************************************************* //
