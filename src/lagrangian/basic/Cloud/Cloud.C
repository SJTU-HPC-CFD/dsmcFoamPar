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
    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();
    const globalMeshData& pData = polyMesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Indexing of equivalent patch on neighbour processor into the
    // procPatches list on the neighbour
    const labelList& procPatchNeighbours = pData.processorPatchNeighbours();

    // Which processors this processor is connected to
    const labelList& neighbourProcs = pData.topology().procNeighbours();

    // Initialise the stepFraction moved for the particles
    for (ParticleType& p : *this)
    {
        p.reset();
    }

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
    const word moveSchedule = cloudOpenMP::moveSchedule(cloud, 0);
    const label moveChunk = max(cloudOpenMP::moveChunk(cloud, 0), label(1));
    const bool useParticlePartition =
        useOpenMPMove && cloudOpenMP::hasParticlePartition(cloud, 0);
    #else
    const bool useOpenMPMove = false;
    const label moveThreads = 1;
    const word moveSchedule = "static";
    const label moveChunk = 64;
    const bool useParticlePartition = false;
    #endif

    // While there are particles to transfer
    while (true)
    {
        cloudOpenMP::clearMoveOrderedParcels(cloud, 0);

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

        if (useOpenMPMove)
        {
            List<ParticleType*> particles;
            labelList threadOffsets(moveThreads + 1, 0);

            if (useParticlePartition)
            {
                const auto& cellOccupancy = cloudOpenMP::cellOccupancy(cloud, 0);
                const auto& particleLoadStart = cloudOpenMP::particleLoadStart(cloud, 0);
                const auto& particleLoadEnd = cloudOpenMP::particleLoadEnd(cloud, 0);

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
                        localCount += cellOccupancy[celli].size();
                    }

                    threadOffsets[threadI + 1] = threadOffsets[threadI] + localCount;
                }

                particles.setSize(threadOffsets.last());

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
                        const auto& cellParcels = cellOccupancy[celli];

                        forAll(cellParcels, i)
                        {
                            particles[particlei++] = cellParcels[i];
                        }
                    }
                }
            }
            else
            {
                particles.setSize(this->size());
                label particlei = 0;

                for (ParticleType& p : *this)
                {
                    particles[particlei++] = &p;
                }

                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    threadOffsets[threadI] = threadI*particles.size()/moveThreads;
                }

                threadOffsets[moveThreads] = particles.size();
            }

            List<label> keepParticleFlags(particles.size(), 1);
            List<label> switchProcessorFlags(particles.size(), 0);
            labelList processedParticleCounts(moveThreads, 0);

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

                if (useParticlePartition)
                {
                    for (label i = threadOffsets[threadI]; i < threadOffsets[threadI + 1]; ++i)
                    {
                        ParticleType& p = *particles[i];

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

                        localTd.switchProcessor = false;
                        localTd.keepParticle = true;

                        keepParticleFlags[i] = p.move(cloud, localTd, trackTime) ? 1 : 0;
                        switchProcessorFlags[i] = localTd.switchProcessor ? 1 : 0;
                        ++localProcessedCount;
                    }
                }

                processedParticleCounts[threadI] = localProcessedCount;
            }

            cloudOpenMP::recordMoveThreadCounts(cloud, processedParticleCounts, 0);

            DynamicList<ParticleType*> survivingParticles(particles.size());
            labelList survivingThreadCounts(moveThreads, 0);

            forAll(particles, i)
            {
                ParticleType& p = *particles[i];
                label ownerThread = 0;

                if (useParticlePartition)
                {
                    while (ownerThread + 1 < threadOffsets.size() && i >= threadOffsets[ownerThread + 1])
                    {
                        ++ownerThread;
                    }
                }

                if (keepParticleFlags[i])
                {
                    if (switchProcessorFlags[i])
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
                    else
                    {
                        survivingParticles.append(&p);

                        if (useParticlePartition && ownerThread < survivingThreadCounts.size())
                        {
                            ++survivingThreadCounts[ownerThread];
                        }
                    }
                }
                else
                {
                    deleteParticle(p);
                }
            }

            if (useParticlePartition)
            {
                List<ParticleType*> survivingList(survivingParticles.size());

                forAll(survivingParticles, i)
                {
                    survivingList[i] = survivingParticles[i];
                }

                labelList survivingOffsets(moveThreads + 1, 0);

                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    survivingOffsets[threadI + 1] =
                        survivingOffsets[threadI] + survivingThreadCounts[threadI];
                }

                cloudOpenMP::storeMoveOrderedParcels(cloud, survivingList, survivingOffsets, 0);
            }
        }
        else
        {
            label processedParticleCount = 0;
            for (ParticleType& p : *this)
            {
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
            cloudOpenMP::recordMoveThreadCounts(cloud, processedCounts, 0);
        }

        if (!Pstream::parRun())
        {
            break;
        }

        pBufs.finishedNeighbourSends(neighbourProcs);

        if (!returnReduceOr(pBufs.hasRecvData()))
        {
            // No parcels to transfer
            break;
        }

        // Retrieve from receive buffers
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
                    addParticle(newp);
                }
            }
        }
    }
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
