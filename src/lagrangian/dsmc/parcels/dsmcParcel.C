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

#include "dsmcParcel.H"
#include "dsmcCloud.H"
#include "meshTools.H"

#include <cstdlib>
#include <cstring>
#include <atomic>
#include <mutex>
#include <new>

// * * * * * * * * * * * * H6 transfer round-trip check * * * * * * * * * * //

// Field-by-field verification of the flat TransferData record (default
// off; enabled by replicatedMeshVerifyTransfer).  The check runs on the
// migration pack path, which is single-threaded per rank, so plain file
// statics are enough.
static bool  dsmcVerifyTransferOn_ = false;
static Foam::label dsmcVerifyTransferLimit_ = 0;
static Foam::label dsmcVerifyTransferCount_ = 0;
static Foam::label dsmcVerifyTransferMismatch_ = 0;

// * * * * * * * * * * * * * * * M6 slab arena  * * * * * * * * * * * * * * //

// §25: per-thread allocation magazine (bounded cache over the global
// arena).  Declared before the arena class since operator new/delete
// below use it.
struct dsmcParcelMagazine
{
    static constexpr int capacity = 64;
    void* slots[capacity];
    int n = 0;
    int tid = -1;      // cached OMP thread index, set at refill (§36.4)
};

// §36.4: per-thread allocation counters.  §36.3 put these inside the arena
// (i.e. under the global spinlock).  With ~48 threads queueing on that
// lock, the added critical-section work is amplified by roughly N/2, and
// the measurement showed it: unpack went 249 -> 1325 ns/parcel and the
// loop +20% (job 4916140 vs the pre-§36 baseline, same test).  The
// counters now live in a padded per-thread array, updated once per
// 64-slot magazine refill, i.e. OUTSIDE the lock and off the shared cache
// line.  Summing the array gives the process total; the per-thread slots
// also avoid false sharing.
struct dsmcAllocCounter
{
    long long v;
    char pad[56];
};
constexpr int dsmcMaxAllocThreads = 256;
static dsmcAllocCounter dsmcAllocSlots_[dsmcMaxAllocThreads];
static dsmcAllocCounter dsmcAllocRefills_[dsmcMaxAllocThreads];

// §36.3: real allocation counters for the unpack/rebuild path.  The
// magazine batches arena traffic (one global-arena lock per 64 slots), so
// the actual allocator work was invisible in the P0 report and the parcel
// count ("append") was often mistaken for an allocation count.  The arena
// counters are incremented under the arena spinlock, so they are exact
// without atomics; unpackTransfer runs in the serial receive loop, so the
// vibrational-scratch counters need no atomics either.
//   allocCalls — global-arena slot hand-outs (the real allocation count;
//                equals parcel constructions rounded up to a magazine)
//   blockAllocs — arena block mallocs (the actual ::malloc traffic)
//   peakInUse  — high-water mark of live slots (arena memory watermark)
//   vibReallocs — DynamicList vibScratch growths inside unpackTransfer
//   vibMaxModes — largest vibrational-mode vector seen in unpackTransfer
static long long dsmcVibReallocCount_ = 0;
static Foam::label dsmcVibMaxModes_ = 0;

// Fixed-size slot allocator for dsmcParcel.  Blocks of `blockParcels_`
// parcels are carved from ::malloc'd chunks; freed slots are recycled
// through an intrusive free list stored in the slot itself (the first
// pointer-sized word, always overwritten by the next construction).
// Thread safety: a per-thread arena avoids locking on the hot path; the
// OMP thread number keys the arena.  MPI ranks never share memory.
namespace Foam
{

class dsmcParcelArena
{
public:

    static const size_t slotSize_ = sizeof(dsmcParcel);

private:

    // 4096 parcels per block ~ 8-10 MB; keeps block count small while
    // bounding waste for shrinking clouds.
    static const label blockParcels_ = 4096;

    struct FreeSlot
    {
        FreeSlot* next_;
    };

    FreeSlot* freeList_;
    DynamicList<void*> blocks_;
    label inUse_;

    // §36.4: blockAllocs_ is incremented only in allocateBlock() (one per
    // 4096 parcels, under the arena lock) — too rare to matter.  The
    // per-slot counters are NOT here: see dsmcAllocSlots_ above.
    label blockAllocs_;

    void allocateBlock()
    {
        const size_t blockBytes = slotSize_*blockParcels_;
        void* mem = std::malloc(blockBytes);
        if (!mem)
        {
            FatalErrorInFunction
                << "dsmcParcel arena: failed to allocate "
                << blockBytes << " bytes" << exit(FatalError);
        }
        blocks_.append(mem);
        ++blockAllocs_;

        // Chain the new slots (ascending addresses for locality).
        char* base = reinterpret_cast<char*>(mem);
        for (label i = blockParcels_ - 1; i >= 0; --i)
        {
            FreeSlot* s = reinterpret_cast<FreeSlot*>(base + i*slotSize_);
            s->next_ = freeList_;
            freeList_ = s;
        }
    }

public:

    dsmcParcelArena()
    :
        freeList_(nullptr),
        blocks_(),
        inUse_(0),
        blockAllocs_(0)
    {}

    ~dsmcParcelArena()
    {
        forAll(blocks_, i)
        {
            std::free(blocks_[i]);
        }
    }

    void* allocate()
    {
        if (!freeList_)
        {
            allocateBlock();
        }
        FreeSlot* s = freeList_;
        freeList_ = s->next_;
        ++inUse_;
        return s;
    }

    void deallocate(void* ptr)
    {
        if (!ptr)
        {
            return;
        }
        FreeSlot* s = reinterpret_cast<FreeSlot*>(ptr);
        s->next_ = freeList_;
        freeList_ = s;
        --inUse_;
    }

    // §36.4 accessors.  blockAllocs_/inUse_ are read at report time, when
    // no thread is inside the arena (the writers are serialized by the
    // arena lock, so the values are consistent).
    label blockAllocs() const { return blockAllocs_; }
    label inUse() const { return inUse_; }

    // M6-ratchet fix: a single global arena.  Per-thread arenas ratchet
    // when parcels are created and destroyed on different threads
    // (boundary injection/outflow, DLB repacking): each thread's arena
    // grows to its own peak and never shrinks, so total memory tracks
    // the SUM of per-thread peaks, which drifts upward without bound.
    // A global arena recycles slots across threads; the spinlock costs
    // ~100ns per op against a churn of ~1% of the cloud per step.
    static dsmcParcelArena& globalArena()
    {
        static dsmcParcelArena arena;
        return arena;
    }

    static std::atomic_flag lock_;

    static void lockArena()
    {
        while (lock_.test_and_set(std::memory_order_acquire))
        {
        }
    }

    static void unlockArena()
    {
        lock_.clear(std::memory_order_release);
    }
};


} // End namespace Foam

#include <chrono>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace
{
using MoveDetailClock = std::chrono::steady_clock;

inline Foam::scalar elapsedMoveDetailSeconds
(
    const MoveDetailClock::time_point& start
)
{
    return std::chrono::duration<Foam::scalar>
    (
        MoveDetailClock::now() - start
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

inline Foam::scalar moveSample01
(
    Foam::dsmcCloud& cloud,
    Foam::dsmcParcel::trackingData& td
)
{
    if (useOpenMPMoveCriticals(cloud) && td.fastRng)
    {
        return td.moveRng.sample01();
    }

    #ifdef _OPENMP
    if (useOpenMPMoveCriticals(cloud))
    {
        Foam::scalar value = 0;
        #pragma omp critical(dsmcMoveRandom)
        {
            value = cloud.rndGen().sample01<Foam::scalar>();
        }
        return value;
    }
    #endif

    return cloud.rndGen().sample01<Foam::scalar>();
}

inline void trackParcelFaceTransitionThreadSafe
(
    Foam::dsmcCloud& cloud,
    const Foam::dsmcParcel& p
)
{
    if (useOpenMPMoveCriticals(cloud))
    {
        #pragma omp critical(dsmcMoveTracker)
        {
            cloud.tracker().trackParcelFaceTransition(p);
        }
    }
    else
    {
        cloud.tracker().trackParcelFaceTransition(p);
    }
}

inline void controlCyclicBoundaryThreadSafe
(
    Foam::dsmcCloud& cloud,
    const Foam::label modelI,
    Foam::dsmcParcel& p,
    Foam::dsmcParcel::trackingData& td
)
{
    if (useOpenMPMoveCriticals(cloud))
    {
        #pragma omp critical(dsmcMoveBoundary)
        {
            cloud.boundaries().cyclicBoundaryModels()[modelI]->controlMol(p, td);
        }
    }
    else
    {
        cloud.boundaries().cyclicBoundaryModels()[modelI]->controlMol(p, td);
    }
}

inline void controlPatchBoundaryThreadSafe
(
    Foam::dsmcCloud& cloud,
    const Foam::label modelI,
    Foam::dsmcParcel& p,
    Foam::dsmcParcel::trackingData& td
)
{
    Foam::dsmcPatchBoundary& model =
        cloud.boundaries().patchBoundaryModels()[modelI]();

    if (useOpenMPMoveCriticals(cloud))
    {
        // Patch models may update shared wall measurements and use cloud RNG.
        #pragma omp critical(dsmcMoveBoundary)
        {
            model.controlParticle(p, td);
        }
    }
    else
    {
        model.controlParticle(p, td);
    }
}
}


// * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * * //

Foam::label Foam::dsmcParcel::TrackedParcel::nDELETED = 0;

const Foam::label Foam::dsmcParcel::trackingData::offloadNBuckets = 12;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcParcel::trackingData::offloadProbeInit(bool observeOn)
{
    offloadProbeOn = observeOn;
    offloadCurBucket = -1;
    offloadCurSafe = false;
    for (label b = 0; b < offloadNBuckets; ++b)
    {
        offHistP[b] = 0;
        offHistT[b] = 0;
        offHistB[b] = 0;
    }
    offSafeP = 0;
    offSafeT = 0;
    offSafeB = 0;
    offMaxMagSqrU = 0.0;
}


void Foam::dsmcParcel::trackingData::offloadBeginParticle
(
    label celli,
    scalar magSqrU,
    scalar dtCell
)
{
    if (!offloadProbeOn)
    {
        offloadCurBucket = -1;
        offloadCurSafe = false;
        return;
    }

    if (magSqrU > offMaxMagSqrU)
    {
        offMaxMagSqrU = magSqrU;
    }

    const label b = cloud().offloadDistBucket(celli);
    offloadCurBucket = (b >= 0 && b < offloadNBuckets) ? b : 0;
    ++offHistP[offloadCurBucket];

    // dlb_offload.md §3: exact per-particle screen.  safeDist_[c] is a
    // conservative lower bound (Dijkstra over cell centres minus the cell
    // circumradius, sources at 0) on the clearance of any point of the
    // cell from any boundary face; a particle is safe when that clearance
    // exceeds twice its step displacement (margin 2).
    const scalar clearance = cloud().offloadCellClearance(celli);
    offloadCurSafe =
        clearance > 0
     && clearance*clearance > 4*magSqrU*dtCell*dtCell;
    if (offloadCurSafe)
    {
        ++offSafeP;
    }
}


void Foam::dsmcParcel::trackingData::offloadProbeAccumulate
(
    const trackingData& other
)
{
    if (!other.offloadProbeOn)
    {
        return;
    }
    for (label b = 0; b < offloadNBuckets; ++b)
    {
        offHistP[b] += other.offHistP[b];
        offHistT[b] += other.offHistT[b];
        offHistB[b] += other.offHistB[b];
    }
    offSafeP += other.offSafeP;
    offSafeT += other.offSafeT;
    offSafeB += other.offSafeB;
    offMaxMagSqrU = max(offMaxMagSqrU, other.offMaxMagSqrU);
}


bool Foam::dsmcParcel::move
(
    dsmcParcel::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    dsmcCloud& cloud = td.cloud();
    const dsmcBoundaries& boundaries = cloud.boundaries();
    const List<label>& cyclicBoundaryToModelIds =
        boundaries.cyclicBoundaryToModelIds();
    const List<label>& patchToModelIds = boundaries.patchToModelIds();
    const bool moveDetailProfile = td.moveDetailProfile;

    if (moveDetailProfile)
    {
        ++td.moveParcels;
    }

    if (isFree())
    {
        const polyMesh& mesh = cloud.pMesh();
        const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
        const bool cartesianTracking =
            cloud.coordSystem().type() == "dsmcCartesian";
        const bool constrainCartesianTracking =
            cartesianTracking
         && (mesh.nGeometricD() < 3 || mesh.nSolutionD() < 3);
        const bool trackerActive = cloud.trackerActive();
        const bool uniformDeltaT = cloud.uniformDeltaT();

        if (newParcel() != -1)
        {
            // note: this justifies that freshly inserted parcels should be
            // tracked as if they had passed the boundary face on which they
            // have been inserted in the time step in which they are inserted.
            stepFraction() = moveSample01(cloud, td);
            newParcel() = -1;
        }

        //scalar tEnd = (1.0 - stepFraction())*trackTime; // OLD FORMULATION
        label orgCell = cell(); // NEW VINCENT
        scalar dtCell = uniformDeltaT ? trackTime : cloud.deltaTValue(orgCell);
        scalar tEnd = (1.0 - stepFraction())*dtCell; // NEW VINCENT

        // dlb_offload.md P1: charge this particle's whole move work to its
        // starting cell (start-cell attribution, not per-crossing cell).
        td.offloadBeginParticle(orgCell, magSqr(U_), dtCell);
        //const scalar dtMax = tEnd; // OLD FORMULATION

        // For reduced-D cases, the velocity used to track needs to be
        // constrained, but the actual U_ of the parcel must not be
        // altered or used, as it is altered by patch interactions one
        // needs to retain its 3D value for collision purposes.
        vector Utracking = U_;

        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            Utracking = U_;

            if (constrainCartesianTracking)
            {
                // Apply correction to position for reduced-D cases,
                // but not for axisymmetric cases
                meshTools::constrainToMeshCentre(mesh, position());

                // Apply correction to velocity to constrain tracking for
                // reduced-D cases,  but not for axisymmetric cases
                meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);
            }

            //- Set the Lagrangian time-step
            //scalar dt = min(dtMax, tEnd); // OLD FORMULATION
            scalar dt = tEnd; // NEW VINCENT

            orgCell = cell(); // NEW VINCENT
            label tetFaceBefore = -1;
            label tetPtBefore = -1;

            td.offloadTrackEvent();

            if (moveDetailProfile)
            {
                ++td.moveTrackCalls;
                tetFaceBefore = tetFace();
                tetPtBefore = tetPt();
                // Flag on the starting cell: cell() is stale after
                // trackToFace.
                const bool nearWallStart =
                    td.nearWallCells
                  ? (*td.nearWallCells)[cell()]
                  : false;
                const auto tTrack0 = MoveDetailClock::now();
                dt *= trackToFace(position() + dt*Utracking, td, true);
                const scalar tTrackElapsed =
                    elapsedMoveDetailSeconds(tTrack0);
                td.moveTrackWallTime += tTrackElapsed;

                // §11.6: split by nearWall vs interior cell.
                if (td.nearWallCells)
                {
                    if (nearWallStart)
                    {
                        ++td.moveWallCellTrackCalls;
                        td.moveWallCellTrackWallTime += tTrackElapsed;
                    }
                    else
                    {
                        ++td.moveInteriorTrackCalls;
                        td.moveInteriorTrackWallTime += tTrackElapsed;
                    }
                }

                if (face() != -1)
                {
                    ++td.moveFaceHits;
                }
                else if
                (
                    cell() == orgCell
                 && tetFace() == tetFaceBefore
                 && tetPt() == tetPtBefore
                )
                {
                    ++td.moveSameTetNoFaceHits;
                }
                else
                {
                    ++td.moveInternalTetNoFaceHits;
                }
            }
            else
            {
                dt *= trackToFace(position() + dt*Utracking, td, true);
            }
            const label destCell = cell(); // NEW VINCENT

            tEnd -= dt;

            if (!uniformDeltaT)
            {
                dtCell = cloud.deltaTValue(orgCell);
            }

            stepFraction() = 1.0 - tEnd/dtCell; // NEW VINCENT
            //stepFraction() = 1.0 - tEnd/trackTime; // OLD FORMULATION

            /*if (destCell != orgCell)
            {
                tEnd *= cloud.deltaTValue(destCell)
                    /cloud.deltaTValue(orgCell);
            } // NEW VINCENT*/

            //- face tracking info
            if (face() != -1 && trackerActive)
            {
                //- measure flux properties
                if (moveDetailProfile)
                {
                    const auto tTracker0 = MoveDetailClock::now();
                    trackParcelFaceTransitionThreadSafe(cloud, *this);
                    td.moveTrackerWallTime += elapsedMoveDetailSeconds(tTracker0);
                }
                else
                {
                    trackParcelFaceTransitionThreadSafe(cloud, *this);
                }
            }

            if (onBoundary() && td.keepParticle)
            {
                const label patchIndex = patch(face());

                if (isA<processorPolyPatch>(pbMesh[patchIndex]))
                {
                    if (moveDetailProfile)
                    {
                        ++td.moveProcessorHits;
                    }
                    td.switchProcessor = true;
                }

                if
                (
                    patchIndex >= 0
                 && patchIndex < cyclicBoundaryToModelIds.size()
                )
                {
                    const label cyclicModelId = cyclicBoundaryToModelIds[patchIndex];

                    if (cyclicModelId >= 0)
                    {
                        td.offloadBoundaryEvent();
                        if (moveDetailProfile)
                        {
                            ++td.moveCyclicHits;
                            const auto tBoundary0 = MoveDetailClock::now();
                            controlCyclicBoundaryThreadSafe
                            (
                                cloud,
                                cyclicModelId,
                                *this,
                                td
                            );
                            td.moveBoundaryWallTime +=
                                elapsedMoveDetailSeconds(tBoundary0);
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
        }
    }
    else
    {
        //- The stuck particle is considered for desorption
        const label patchIndex = stuck().wallTemperature()[1];

        if
        (
            patchIndex >= 0
         && patchIndex < patchToModelIds.size()
        )
        {
            const label patchModelId = patchToModelIds[patchIndex];

            if (patchModelId >= 0)
            {
                // then this patch is the "dsmc*Sticking*WallPatch" the particle is stuck on
                if (moveDetailProfile)
                {
                    ++td.moveStuckHits;
                    const auto tBoundary0 = MoveDetailClock::now();
                    controlPatchBoundaryThreadSafe
                    (
                        cloud,
                        patchModelId,
                        *this,
                        td
                    );
                    td.moveBoundaryWallTime +=
                        elapsedMoveDetailSeconds(tBoundary0);
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


bool Foam::dsmcParcel::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::dsmcParcel::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::dsmcParcel::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    //- find which patch has been hit
    label patchIndex = wpp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];

    //- apply a boundary model when a molecule collides with this poly patch
    td.offloadBoundaryEvent();
    if (td.moveDetailProfile)
    {
        ++td.movePatchHits;
        const auto tBoundary0 = MoveDetailClock::now();
        controlPatchBoundaryThreadSafe(td.cloud(), patchModelId, *this, td);
        td.moveBoundaryWallTime += elapsedMoveDetailSeconds(tBoundary0);
    }
    else
    {
        controlPatchBoundaryThreadSafe(td.cloud(), patchModelId, *this, td);
    }
}


void Foam::dsmcParcel::hitPatch
(
    const polyPatch& pp,
    trackingData& td
)
{
    //- find which patch has been hit
    label patchIndex = pp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];

    //- apply a boundary model when a molecule collides with this poly patch
    td.offloadBoundaryEvent();
    if (td.moveDetailProfile)
    {
        ++td.movePatchHits;
        const auto tBoundary0 = MoveDetailClock::now();
        controlPatchBoundaryThreadSafe(td.cloud(), patchModelId, *this, td);
        td.moveBoundaryWallTime += elapsedMoveDetailSeconds(tBoundary0);
    }
    else
    {
        controlPatchBoundaryThreadSafe(td.cloud(), patchModelId, *this, td);
    }
}


void Foam::dsmcParcel::transformProperties
(
    const tensor& T
)
{
   particle::transformProperties(T);
   U_ = transform(T, U_);
}


void Foam::dsmcParcel::transformProperties
(
    const vector& separation
)
{
    particle::transformProperties(separation);
}


// * * * * * * * * * * * non-flat transfer record  * * * * * * * * * * * * //

namespace
{

//- Version of the explicit non-flat migration record.  A reader rejects
//  any other version: records from a different layout must not be mixed.
const Foam::label dsmcTransferRecordVersion = 1;

//- Trailing marker of a complete record (truncation/misalignment check).
const Foam::label dsmcTransferRecordEnd = 0x5A17C0DE;

template<class T>
inline void dsmcWriteRaw(Foam::Ostream& os, const T& v)
{
    os.write
    (
        reinterpret_cast<const char*>(&v),
        std::streamsize(sizeof(T))
    );
}

template<class T>
inline void dsmcReadRaw(Foam::Istream& is, T& v)
{
    is.read(reinterpret_cast<char*>(&v), std::streamsize(sizeof(T)));
}

//- Fixed part of the record.  The variable-length data (stuck wall
//  points, vibrational levels) follows it, then the end marker.
struct DsmcTransferRecordFixed
{
    Foam::label version;

    Foam::scalar position[3];
    Foam::scalar stepFraction;
    Foam::label celli;
    Foam::label facei;
    Foam::label tetFacei;
    Foam::label tetPti;
    Foam::label origProc;
    Foam::label origId;

    Foam::scalar U[3];
    Foam::scalar RWF;
    Foam::scalar ERot;
    Foam::label ELevel;
    Foam::label typeId;
    Foam::label newParcel;
    Foam::label classification;

    Foam::label tracked;
    Foam::label inPatchId;
    Foam::label storePositions;
    Foam::scalar t0;
    Foam::scalar initialPosition[3];
    Foam::scalar currentPosition[3];
    Foam::scalar distanceTravelled[3];
    Foam::scalar msd[3];

    Foam::label nWallPoints;
    Foam::label nVibModes;
};

} // End anonymous namespace


void Foam::dsmcParcel::writeBinaryFast(Ostream& os) const
{
    // clean N1: explicit, versioned record.  The previous implementation
    // streamed the parcel through operator<< (ASCII tokens plus framed raw
    // blocks, including the base-particle serializer).  That made the
    // record layout implicit and length-dependent on the token reader, and
    // a mismatch only surfaced later as a stream error during migration.
    // Here every field has an explicit layout and the record ends with a
    // marker; see readTransferRecord().
    DsmcTransferRecordFixed r;
    std::memset(&r, 0, sizeof(r));

    r.version = dsmcTransferRecordVersion;

    r.position[0] = position_.x();
    r.position[1] = position_.y();
    r.position[2] = position_.z();
    r.stepFraction = stepFraction_;
    r.celli = cellI_;
    r.facei = faceI_;
    r.tetFacei = tetFaceI_;
    r.tetPti = tetPtI_;
    r.origProc = origProc_;
    r.origId = origId_;

    r.U[0] = U_.x();
    r.U[1] = U_.y();
    r.U[2] = U_.z();
    r.RWF = RWF_;
    r.ERot = ERot_;
    r.ELevel = ELevel_;
    r.typeId = typeId_;
    r.newParcel = newParcel_;
    r.classification = classification_;

    r.tracked = tracked_.tracked() ? 1 : 0;
    r.inPatchId = tracked_.inPatchId();
    r.storePositions = tracked_.storePositions() ? 1 : 0;
    r.t0 = tracked_.initialTime();
    const vector& tIP = tracked_.initialPosition();
    const vector& tCP = tracked_.currentPosition();
    const vector& tDT = tracked_.distanceTravelledVector();
    const vector& tMSD = tracked_.meanSquaredDisplacementVector();
    r.initialPosition[0] = tIP.x();
    r.initialPosition[1] = tIP.y();
    r.initialPosition[2] = tIP.z();
    r.currentPosition[0] = tCP.x();
    r.currentPosition[1] = tCP.y();
    r.currentPosition[2] = tCP.z();
    r.distanceTravelled[0] = tDT.x();
    r.distanceTravelled[1] = tDT.y();
    r.distanceTravelled[2] = tDT.z();
    r.msd[0] = tMSD.x();
    r.msd[1] = tMSD.y();
    r.msd[2] = tMSD.z();

    r.nWallPoints = stuck_ ? stuck_->wallTemperature().size() : 0;
    r.nVibModes = vibLevel_.size();

    dsmcWriteRaw(os, r);

    if (r.nWallPoints > 0)
    {
        const scalarField& wT = stuck_->wallTemperature();
        const vectorField& wV = stuck_->wallVectors();
        os.write
        (
            reinterpret_cast<const char*>(wT.cdata()),
            std::streamsize(wT.size()*sizeof(scalar))
        );
        os.write
        (
            reinterpret_cast<const char*>(wV.cdata()),
            std::streamsize(wV.size()*sizeof(vector))
        );
    }

    if (r.nVibModes > 0)
    {
        os.write
        (
            reinterpret_cast<const char*>(vibLevel_.cdata()),
            std::streamsize(r.nVibModes*sizeof(label))
        );
    }

    const label endMarker = dsmcTransferRecordEnd;
    dsmcWriteRaw(os, endMarker);

    os.check("dsmcParcel::writeBinaryFast");
}


Foam::dsmcParcel* Foam::dsmcParcel::readTransferRecord
(
    const polyMesh& mesh,
    Istream& is
)
{
    DsmcTransferRecordFixed r;
    dsmcReadRaw(is, r);

    if (r.version != dsmcTransferRecordVersion)
    {
        FatalErrorInFunction
            << "Unsupported non-flat parcel record version " << r.version
            << " (this build writes/reads version "
            << dsmcTransferRecordVersion << ")"
            << abort(FatalError);
    }
    if
    (
        r.nWallPoints < 0
     || r.nWallPoints > 1024
     || r.nVibModes < 0
     || r.nVibModes > label(maxVibModes)
    )
    {
        FatalErrorInFunction
            << "Invalid non-flat parcel record: nWallPoints "
            << r.nWallPoints << " nVibModes " << r.nVibModes
            << abort(FatalError);
    }

    scalarField wallT;
    vectorField wallV;
    if (r.nWallPoints > 0)
    {
        wallT.setSize(r.nWallPoints);
        wallV.setSize(r.nWallPoints);
        is.read
        (
            reinterpret_cast<char*>(wallT.data()),
            std::streamsize(r.nWallPoints*sizeof(scalar))
        );
        is.read
        (
            reinterpret_cast<char*>(wallV.data()),
            std::streamsize(r.nWallPoints*sizeof(vector))
        );
    }

    labelList vib(r.nVibModes);
    if (r.nVibModes > 0)
    {
        is.read
        (
            reinterpret_cast<char*>(vib.data()),
            std::streamsize(r.nVibModes*sizeof(label))
        );
    }

    label endMarker = 0;
    dsmcReadRaw(is, endMarker);
    if (endMarker != dsmcTransferRecordEnd)
    {
        FatalErrorInFunction
            << "Truncated or misaligned non-flat parcel record: end marker "
            << endMarker << " (expected " << dsmcTransferRecordEnd << ")"
            << abort(FatalError);
    }

    auto* p = new dsmcParcel
    (
        mesh,
        vector(r.position[0], r.position[1], r.position[2]),
        vector(r.U[0], r.U[1], r.U[2]),
        r.RWF,
        r.ERot,
        r.ELevel,
        r.celli,
        r.tetFacei,
        r.tetPti,
        r.typeId,
        r.newParcel,
        r.classification,
        vib
    );

    p->face() = r.facei;
    p->stepFraction() = r.stepFraction;
    p->origProc() = r.origProc;
    p->origId() = r.origId;

    p->tracked_.setMigrationState
    (
        r.tracked != 0,
        r.inPatchId,
        r.storePositions != 0,
        r.t0,
        vector
        (
            r.initialPosition[0], r.initialPosition[1], r.initialPosition[2]
        ),
        vector
        (
            r.currentPosition[0], r.currentPosition[1], r.currentPosition[2]
        ),
        vector
        (
            r.distanceTravelled[0],
            r.distanceTravelled[1],
            r.distanceTravelled[2]
        ),
        vector(r.msd[0], r.msd[1], r.msd[2])
    );

    if (r.nWallPoints > 0)
    {
        p->stuck_ = new StuckParcel(wallT, wallV);
    }

    return p;
}


void Foam::dsmcParcel::packTransfer(TransferData& td) const
{
    // State contract: TransferData only carries the move input below.
    // Stuck (wall temperature/momentum) and tracked (tracer) states are not
    // representable, and the vibrational capacity is fixed — refuse instead
    // of silently resetting or truncating them.
    if (isStuck())
    {
        FatalErrorInFunction
            << "Cannot transfer a stuck parcel: the StuckParcel wall "
            << "temperature/momentum state is not part of TransferData. "
            << "Keep sticking-wall states out of replicated/DLB migration "
            << "or extend the transfer schema."
            << abort(FatalError);
    }
    if (tracked().tracked())
    {
        FatalErrorInFunction
            << "Cannot transfer a tracked parcel: the TrackedParcel tracer "
            << "state is not part of TransferData."
            << abort(FatalError);
    }
    if (vibLevel_.size() > label(maxVibModes))
    {
        FatalErrorInFunction
            << "Cannot transfer a parcel with " << vibLevel_.size()
            << " vibrational modes (TransferData capacity "
            << int(maxVibModes) << "): silent truncation is not allowed."
            << abort(FatalError);
    }

    td.position[0] = position().x();
    td.position[1] = position().y();
    td.position[2] = position().z();
    td.celli = cell();
    td.tetFacei = tetFace();
    td.tetPti = tetPt();

    td.U[0] = U_.x();
    td.U[1] = U_.y();
    td.U[2] = U_.z();
    td.RWF = RWF_;
    td.ERot = ERot_;
    td.ELevel = ELevel_;
    td.typeId = typeId_;
    td.newParcel = newParcel_;
    td.classification = classification_;

    td.nVibModes = vibLevel_.size();
    for (label i = 0; i < td.nVibModes; ++i)
    {
        td.vibLevel[i] = vibLevel_[i];
    }
    for (label i = td.nVibModes; i < maxVibModes; ++i)
    {
        td.vibLevel[i] = 0;
    }

    // H6: compare the record with a real unpack, field by field.  The
    // values are copied bit-for-bit through the POD, so any difference is
    // a schema/offset defect, not rounding.  Identity fields carried by
    // the base particle (origProc_/origId_) are not part of TransferData
    // and are re-assigned on insertion; that limitation is explicitly
    // recorded in the state contract.
    if (dsmcVerifyTransferOn_ && dsmcVerifyTransferCount_ < dsmcVerifyTransferLimit_)
    {
        ++dsmcVerifyTransferCount_;

        dsmcParcel* back = dsmcParcel::unpackTransfer(mesh(), td);

        const char* firstName = nullptr;
        scalar firstA = 0;
        scalar firstB = 0;
        const label nDiff =
            compareTransferFields(*this, *back, firstName, firstA, firstB);
        delete back;

        if (nDiff != 0)
        {
            ++dsmcVerifyTransferMismatch_;
            FatalErrorInFunction
                << "H6 flat transfer round-trip mismatch: " << nDiff
                << " field(s) differ, first " << (firstName ? firstName : "?")
                << " unpacked " << firstA << " != packed " << firstB
                << " (cell " << cell() << " type " << typeId_ << ")"
                << abort(FatalError);
        }
    }
}


void Foam::dsmcParcel::verifyTransferRoundTrip
(
    const bool enabled,
    const label limit
)
{
    dsmcVerifyTransferOn_ = enabled;
    dsmcVerifyTransferLimit_ = max(label(0), limit);
    dsmcVerifyTransferCount_ = 0;
    dsmcVerifyTransferMismatch_ = 0;
}


void Foam::dsmcParcel::transferRoundTripStats
(
    label& nVerified,
    label& nMismatch
)
{
    nVerified = dsmcVerifyTransferCount_;
    nMismatch = dsmcVerifyTransferMismatch_;
}


Foam::label Foam::dsmcParcel::compareTransferFields
(
    const dsmcParcel& a,
    const dsmcParcel& b,
    const char*& firstName,
    scalar& firstA,
    scalar& firstB
)
{
    label nDiff = 0;
    firstName = nullptr;
    firstA = 0;
    firstB = 0;

    auto cmp = [&](const char* name, const scalar va, const scalar vb)
    {
        if (va != vb)
        {
            ++nDiff;
            if (!firstName)
            {
                firstName = name;
                firstA = va;
                firstB = vb;
            }
        }
    };

    cmp("position.x", a.position().x(), b.position().x());
    cmp("position.y", a.position().y(), b.position().y());
    cmp("position.z", a.position().z(), b.position().z());
    cmp("cell", scalar(a.cell()), scalar(b.cell()));
    cmp("tetFace", scalar(a.tetFace()), scalar(b.tetFace()));
    cmp("tetPt", scalar(a.tetPt()), scalar(b.tetPt()));
    cmp("U.x", a.U_.x(), b.U_.x());
    cmp("U.y", a.U_.y(), b.U_.y());
    cmp("U.z", a.U_.z(), b.U_.z());
    cmp("RWF", a.RWF_, b.RWF_);
    cmp("ERot", a.ERot_, b.ERot_);
    cmp("ELevel", scalar(a.ELevel_), scalar(b.ELevel_));
    cmp("typeId", scalar(a.typeId_), scalar(b.typeId_));
    cmp("newParcel", scalar(a.newParcel_), scalar(b.newParcel_));
    cmp("classification", scalar(a.classification_), scalar(b.classification_));
    cmp("nVibModes", scalar(a.vibLevel_.size()), scalar(b.vibLevel_.size()));
    if (a.vibLevel_.size() == b.vibLevel_.size())
    {
        forAll(a.vibLevel_, i)
        {
            cmp("vibLevel", scalar(a.vibLevel_[i]), scalar(b.vibLevel_[i]));
        }
    }

    // Stuck / tracked state (the non-flat serializer carries it; the flat
    // record refuses it, so both sides are absent there and compare equal).
    cmp("isStuck", a.isStuck() ? 1 : 0, b.isStuck() ? 1 : 0);
    if (a.isStuck() && b.isStuck())
    {
        const scalarField& ta = a.stuck().wallTemperature();
        const scalarField& tb = b.stuck().wallTemperature();
        cmp("stuckWT.size", scalar(ta.size()), scalar(tb.size()));
        if (ta.size() == tb.size())
        {
            forAll(ta, i)
            {
                cmp("stuckWT", ta[i], tb[i]);
            }
        }
        const vectorField& va = a.stuck().wallVectors();
        const vectorField& vb = b.stuck().wallVectors();
        cmp("stuckWV.size", scalar(va.size()), scalar(vb.size()));
        if (va.size() == vb.size())
        {
            forAll(va, i)
            {
                cmp("stuckWV.x", va[i].x(), vb[i].x());
                cmp("stuckWV.y", va[i].y(), vb[i].y());
                cmp("stuckWV.z", va[i].z(), vb[i].z());
            }
        }
    }

    cmp("tracked", a.tracked().tracked() ? 1 : 0, b.tracked().tracked() ? 1 : 0);
    cmp
    (
        "inPatchId",
        scalar(a.tracked().inPatchId()),
        scalar(b.tracked().inPatchId())
    );
    cmp
    (
        "storePositions",
        a.tracked().storePositions() ? 1 : 0,
        b.tracked().storePositions() ? 1 : 0
    );
    cmp("tracked.t0", a.tracked().initialTime(), b.tracked().initialTime());
    cmp
    (
        "tracked.initialPosition.x",
        a.tracked().initialPosition().x(),
        b.tracked().initialPosition().x()
    );
    cmp
    (
        "tracked.initialPosition.y",
        a.tracked().initialPosition().y(),
        b.tracked().initialPosition().y()
    );
    cmp
    (
        "tracked.initialPosition.z",
        a.tracked().initialPosition().z(),
        b.tracked().initialPosition().z()
    );
    cmp
    (
        "tracked.currentPosition.x",
        a.tracked().currentPosition().x(),
        b.tracked().currentPosition().x()
    );
    cmp
    (
        "tracked.currentPosition.y",
        a.tracked().currentPosition().y(),
        b.tracked().currentPosition().y()
    );
    cmp
    (
        "tracked.currentPosition.z",
        a.tracked().currentPosition().z(),
        b.tracked().currentPosition().z()
    );
    cmp
    (
        "tracked.distance.x",
        a.tracked().distanceTravelledVector().x(),
        b.tracked().distanceTravelledVector().x()
    );
    cmp
    (
        "tracked.distance.y",
        a.tracked().distanceTravelledVector().y(),
        b.tracked().distanceTravelledVector().y()
    );
    cmp
    (
        "tracked.distance.z",
        a.tracked().distanceTravelledVector().z(),
        b.tracked().distanceTravelledVector().z()
    );
    cmp("tracked.msd", a.tracked().meanSquaredDisplacementVector().x(),
        b.tracked().meanSquaredDisplacementVector().x());
    cmp("tracked.msd.y", a.tracked().meanSquaredDisplacementVector().y(),
        b.tracked().meanSquaredDisplacementVector().y());
    cmp("tracked.msd.z", a.tracked().meanSquaredDisplacementVector().z(),
        b.tracked().meanSquaredDisplacementVector().z());

    return nDiff;
}


Foam::dsmcParcel* Foam::dsmcParcel::unpackTransfer
(
    const polyMesh& mesh,
    const TransferData& td
){
    // §34.3: OF-v1706 List::setSize() reallocates whenever the size
    // changes (List.C:283 — no capacity retention; verified), so the
    // previous labelList scratch still allocated for every parcel whose
    // vibrational-mode count alternated between species — which made the
    // §26 "no benefit" measurement uninformative about allocation cost.
    // A DynamicList keeps its storage across clear().
    static thread_local DynamicList<label> vibScratch;
    const label vibCapBefore = vibScratch.capacity();
    vibScratch.clear();
    for (label i = 0; i < td.nVibModes; ++i)
    {
        vibScratch.append(td.vibLevel[i]);
    }

    // §36.3: a capacity change is a real reallocation; counting them shows
    // whether the DynamicList actually removed the per-parcel allocation
    // (a single growth covers every later parcel of the same/smaller size).
    if (vibScratch.capacity() != vibCapBefore)
    {
        ++dsmcVibReallocCount_;
    }
    if (td.nVibModes > dsmcVibMaxModes_)
    {
        dsmcVibMaxModes_ = td.nVibModes;
    }

    return new dsmcParcel
    (
        mesh,
        vector(td.position[0], td.position[1], td.position[2]),
        vector(td.U[0], td.U[1], td.U[2]),
        td.RWF,
        td.ERot,
        td.ELevel,
        td.celli,
        td.tetFacei,
        td.tetPti,
        td.typeId,
        td.newParcel,
        td.classification,
        vibScratch
    );
}

bool Foam::dsmcParcel::relocateStuckParcel
(
    const polyMesh& mesh
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // find the closest patch and patch face indices, as this is the patch/face
    // on which this parcel is stuck.
    scalar closestFaceDistance = GREAT;
    label closestPatchi = -1;
    label closestPatchFacei = -1;

    forAll(mesh.cells()[cell()], i)
    {
        const label facei = mesh.cells()[cell()][i];

        // find corresponding boundary patch
        const label patchi = bMesh.whichPatch(facei);

        if (patchi != -1)
        {
            // this is a boundary patch, i.e. a potential candidate for the
            // patch to which parcel is stuck.
            const polyPatch& wpp = bMesh[patchi];

            // patch face index:
            const label patchFacei = wpp.whichFace(facei);

            // calculate distance between parcel position and this patch face:
            pointHit pHit
            (
                wpp[patchFacei].nearestPoint(position(), wpp.points())
            );

            if (pHit.hit())
            {
                const scalar distance = pHit.distance();
                if (distance < closestFaceDistance)
                {
                    // found new closest patch face
                    closestFaceDistance = distance;
                    closestPatchi = patchi;
                    closestPatchFacei = patchFacei;
                }
            }
        }
    }

    // did we find a closest patch face?
    if (closestPatchi != -1)
    {
        // yes, reset stuck parcel information to this patch face
        stuck().wallTemperature()[1] = closestPatchi;
        stuck().wallTemperature()[2] = closestPatchFacei;
        return true;
    }
    // no, indicates fatal error
    return false;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Operators  * * * * * * * * * * * * * * * * //

std::atomic_flag Foam::dsmcParcelArena::lock_ = ATOMIC_FLAG_INIT;

// §36.4: unpack-path allocation counters, read by the replicated-mesh P0
// report.  No arena lock is needed: the per-slot counters above are
// per-thread and only written while the owning thread is inside a parcel
// allocation; the report runs after every parallel region has been left,
// so the values are stable.  blockAllocs_/inUse_ are serialized by the
// arena lock anyway.
Foam::label Foam::dsmcParcel::unpackAllocCount()
{
    long long total = 0;
    for (int i = 0; i < dsmcMaxAllocThreads; ++i)
    {
        total += dsmcAllocSlots_[i].v;
    }
    return label(total);
}

Foam::label Foam::dsmcParcel::unpackMagazineRefillCount()
{
    long long total = 0;
    for (int i = 0; i < dsmcMaxAllocThreads; ++i)
    {
        total += dsmcAllocRefills_[i].v;
    }
    return label(total);
}

Foam::label Foam::dsmcParcel::unpackArenaBlockCount()
{
    return dsmcParcelArena::globalArena().blockAllocs();
}

Foam::label Foam::dsmcParcel::unpackArenaInUse()
{
    return dsmcParcelArena::globalArena().inUse();
}

Foam::label Foam::dsmcParcel::unpackVibReallocCount()
{
    return label(dsmcVibReallocCount_);
}

Foam::label Foam::dsmcParcel::unpackVibMaxModes()
{
    return dsmcVibMaxModes_;
}

// M7 phase 1: dedicated arena for StuckParcel cold data.  Stuck parcels are
// a small fraction of the cloud, so a single shared arena (no per-thread
// split) is sufficient; a mutex guards the rare concurrent wall-stick.
namespace
{

Foam::dsmcParcelArena& stuckParcelArena()
{
    static Foam::dsmcParcelArena arena;
    return arena;
}

} // anonymous namespace


void* Foam::dsmcParcel::StuckParcel::operator new(size_t sz)
{
    // StuckParcel is small and allocation is rare: the general arena slot
    // size is larger than needed, so use a simple mutex-protected pool of
    // exact-size slots via the shared parcel arena of the master thread.
    if (sz != sizeof(Foam::dsmcParcel::StuckParcel))
    {
        return std::malloc(sz);
    }
    static std::mutex m;
    std::lock_guard<std::mutex> lock(m);
    return stuckParcelArena().allocate();
}


void Foam::dsmcParcel::StuckParcel::operator delete(void* ptr) noexcept
{
    if (!ptr)
    {
        return;
    }
    static std::mutex m;
    std::lock_guard<std::mutex> lock(m);
    stuckParcelArena().deallocate(ptr);
}


void* Foam::dsmcParcel::operator new(size_t sz)
{
    if (sz != dsmcParcelArena::slotSize_)
    {
        // Derived-class or padded request: fall back to the heap.
        static bool warnedOnce = false;
        if (!warnedOnce)
        {
            warnedOnce = true;
            Info<< "dsmcParcel arena FALLBACK: sz=" << label(sz)
                << " slotSize_=" << label(dsmcParcelArena::slotSize_)
                << " (this message prints once)" << endl;
        }
        return std::malloc(sz);
    }

    // §25: per-thread magazine.  The global arena spinlock costs ~100ns
    // per op at the design's spread churn, but a migration is a 48-thread
    // allocation BURST: the deser probe (job 4844354) measured 2.36us per
    // allocation in the parallel unpack loop while the serial append
    // (no allocation) ran at 83ns.  Each thread now takes a fixed-size
    // batch under one lock and pops locally, so the shared cache line is
    // touched once per 64 allocations.  The magazine is a bounded cache
    // over the SINGLE global arena — the M6 ratchet fix (no per-thread
    // pools that grow without bound) is preserved.
    static thread_local dsmcParcelMagazine magazine;
    if (magazine.n == 0)
    {
        #ifdef _OPENMP
        const int ompTid = omp_in_parallel() ? omp_get_thread_num() : 0;
        #else
        const int ompTid = 0;
        #endif
        magazine.tid =
        (
            ompTid >= 0 && ompTid < dsmcMaxAllocThreads ? ompTid : 0
        );

        dsmcParcelArena& arena = dsmcParcelArena::globalArena();
        dsmcParcelArena::lockArena();
        for (label i = 0; i < magazine.capacity; ++i)
        {
            magazine.slots[i] = arena.allocate();
        }
        dsmcParcelArena::unlockArena();
        magazine.n = magazine.capacity;

        // §36.4: count the batch here — per-thread, outside the lock.
        dsmcAllocSlots_[magazine.tid].v += magazine.capacity;
        ++dsmcAllocRefills_[magazine.tid].v;
    }
    return magazine.slots[--magazine.n];
}


void* Foam::dsmcParcel::operator new(size_t, void* ptr) noexcept
{
    return ptr;
}


void Foam::dsmcParcel::operator delete(void* ptr) noexcept
{
    if (!ptr)
    {
        return;
    }
    // §25: return to the thread's magazine; flush the whole batch to the
    // global arena only when the magazine is full (one lock per 64 frees).
    // Slots from the fallback heap path cannot be distinguished from
    // arena slots by value; the arena only ever receives sizeof(dsmcParcel)
    // requests, so any delete of a base-size parcel is an arena slot.
    static thread_local dsmcParcelMagazine magazine;
    if (magazine.n < magazine.capacity)
    {
        magazine.slots[magazine.n++] = ptr;
        return;
    }

    dsmcParcelArena& arena = dsmcParcelArena::globalArena();
    dsmcParcelArena::lockArena();
    for (label i = 0; i < magazine.n; ++i)
    {
        arena.deallocate(magazine.slots[i]);
    }
    arena.deallocate(ptr);
    dsmcParcelArena::unlockArena();
    magazine.n = 0;
}


void Foam::dsmcParcel::operator delete(void* ptr, size_t sz) noexcept
{
    Foam::dsmcParcel::operator delete(ptr);
}


#include "dsmcParcelIO.C"


// ************************************************************************* //
