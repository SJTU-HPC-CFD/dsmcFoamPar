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
#include <mutex>
#include <new>

// * * * * * * * * * * * * * * * M6 slab arena  * * * * * * * * * * * * * * //

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
        inUse_(0)
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

    static dsmcParcelArena& threadArena()
    {
        #ifdef _OPENMP
        static const label nMax = 256;
        label tid = omp_get_thread_num();
        if (tid >= nMax)
        {
            tid = 0;  // fall back to the master arena
        }
        #else
        const label tid = 0;
        #endif
        static dsmcParcelArena arenas[256];
        return arenas[tid];
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


void Foam::dsmcParcel::writeBinaryFast(Ostream& os) const
{
    os << *this;
    os.check("dsmcParcel::writeBinaryFast");
}


void Foam::dsmcParcel::packTransfer(TransferData& td) const
{
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

    td.nVibModes = min(vibLevel_.size(), label(maxVibModes));
    for (label i = 0; i < td.nVibModes; ++i)
    {
        td.vibLevel[i] = vibLevel_[i];
    }
    for (label i = td.nVibModes; i < maxVibModes; ++i)
    {
        td.vibLevel[i] = 0;
    }
}


Foam::dsmcParcel* Foam::dsmcParcel::unpackTransfer
(
    const polyMesh& mesh,
    const TransferData& td
)
{
    labelList vib(td.nVibModes);
    for (label i = 0; i < td.nVibModes; ++i)
    {
        vib[i] = td.vibLevel[i];
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
        vib
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
        return std::malloc(sz);
    }
    return dsmcParcelArena::threadArena().allocate();
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
    // Slots from the fallback heap path cannot be distinguished from arena
    // slots by value; the arena only ever receives sizeof(dsmcParcel)
    // requests, so any delete of a base-size parcel is an arena slot.
    dsmcParcelArena::threadArena().deallocate(ptr);
}


void Foam::dsmcParcel::operator delete(void* ptr, size_t sz) noexcept
{
    Foam::dsmcParcel::operator delete(ptr);
}


#include "dsmcParcelIO.C"


// ************************************************************************* //
