/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.
    hyStrath is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    hyStrath is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
\*---------------------------------------------------------------------------*/

#include "dsmcReplicatedMesh.H"
#include "dsmcCloud.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "IOdictionary.H"
#include "IOField.H"
#include "IOPosition.H"
#include "labelIOList.H"
#include "HashSet.H"
#include "passiveParticleCloud.H"
#include "volFields.H"
#include "scotchDecomp.H"
#include "domainDecomposition.H"
#include "fvFieldDecomposer.H"
#include "parmetis.h"
#include <mpi.h>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>

namespace Foam
{

class ISpanStream
:
    public IStringStream
{
public:
    ISpanStream
    (
        const char* buffer,
        const label size,
        streamFormat format = IOstream::BINARY
    )
    :
        IStringStream(string(buffer, size), format)
    {}
};


template<class GeoField>
label writeProcessorVolumeFields
(
    const fvMesh& completeMesh,
    const fvFieldDecomposer& decomposer
)
{
    HashTable<const GeoField*> fields = completeMesh.lookupClass<GeoField>();
    label nWritten = 0;

    for
    (
        typename HashTable<const GeoField*>::const_iterator iter =
            fields.cbegin();
        iter != fields.cend();
        ++iter
    )
    {
        const GeoField& field = *iter();
        if (field.writeOpt() != IOobject::AUTO_WRITE)
        {
            continue;
        }
        decomposer.decomposeField(field)().write();
        ++nWritten;
    }

    return nWritten;
}

// ============================================================================
// Constructor / Destructor
// ============================================================================

void dsmcReplicatedMesh::addMigrationProfileSample
(
    const scalar wallTime,
    const scalar commTime
)
{
    migrationCommWallTime_ += commTime;

    switch (migrationProfileContext_)
    {
        case migrationProfileRegular:
            regularMigrationWallTime_ += wallTime;
            regularMigrationCommWallTime_ += commTime;
            break;
        case migrationProfileDLB:
            dlbMigrationWallTime_ += wallTime;
            dlbMigrationCommWallTime_ += commTime;
            break;
        case migrationProfileOutput:
            outputMigrationWallTime_ += wallTime;
            outputMigrationCommWallTime_ += commTime;
            break;
        case migrationProfileManual:
            manualMigrationWallTime_ += wallTime;
            manualMigrationCommWallTime_ += commTime;
            break;
        case migrationProfileInitial:
            initialMigrationWallTime_ += wallTime;
            initialMigrationCommWallTime_ += commTime;
            break;
        default:
            otherMigrationWallTime_ += wallTime;
            otherMigrationCommWallTime_ += commTime;
            break;
    }
}


void dsmcReplicatedMesh::addAutoRebalanceCommTime(const scalar commTime)
{
    autoRebalanceExplicitCommWallTime_ += commTime;
}


dsmcReplicatedMesh::dsmcReplicatedMesh(dsmcCloud& cloud, const fvMesh& mesh)
:
    cloud_(cloud), mesh_(mesh),
    cellOwner_(mesh.nCells(), -1), localMesh_(mesh), myCells_(0),
    localParticleCount_(0), allParticleCounts_(0),
    migrationWallTime_(0.0),
    migrationPackWallTime_(0.0),
    migrationLocalPrepWallTime_(0.0),
    migrationSizeExchangeWallTime_(0.0),
    migrationRequestPostWallTime_(0.0),
    migrationWaitWallTime_(0.0),
    migrationDeserializeWallTime_(0.0),
    migrationCandidateGatherWallTime_(0.0),
    migrationPostWallTime_(0.0),
    updateParticleCountsWallTime_(0.0),
    migrationCalls_(0),
    evolveStepTime_(0.0), evolveTimeSteps_(0),
    migrationProfileContext_(migrationProfileOther),
    migrationCommWallTime_(0.0),
    regularMigrationWallTime_(0.0),
    regularMigrationCommWallTime_(0.0),
    dlbMigrationWallTime_(0.0),
    dlbMigrationCommWallTime_(0.0),
    outputMigrationWallTime_(0.0),
    outputMigrationCommWallTime_(0.0),
    outputGatherCommWallTime_(0.0),
    manualMigrationWallTime_(0.0),
    manualMigrationCommWallTime_(0.0),
    initialMigrationWallTime_(0.0),
    initialMigrationCommWallTime_(0.0),
    otherMigrationWallTime_(0.0),
    otherMigrationCommWallTime_(0.0),
    active_(false), nProcs_(1), myRank_(0),
    migrateInterval_(1), stepCounter_(0),
    rebalanceCount_(0), totalCellsChanged_(0), totalParcelsMigrated_(0),
    // Phase C defaults
    autoDLBEnabled_(false),
    imbalanceThreshold_(1.15),
    dlbSteps_(50),
    sarSteps_(10),
    autoDLBTriggerMode_("cumulative"),
    forcedDLBSteps_(),
    lastAutoRebalanceStep_(-1000),
    autoRebalanceCount_(0),
    // ParDSMC3D trigger state
    productiveTime_(0.0),
    lastEvolveTime_(0.0),
    lastMigrationTime_(0.0),
    tidl_(0.0),
    tdecps_(0.0),
    w1_(0.0), w2_(0.0),
    ndecps_(0),
    sarEvalSteps_(0),
    nEvalPeriods_(0),
    sar_(0.0),
    postDLBSnapshotCountdown_(0),
    autoRebalanceChecks_(0),
    autoRebalanceTriggeredChecks_(0),
    autoRebalanceWallTime_(0.0),
    autoRebalanceCheckWallTime_(0.0),
    autoRebalanceRepartWallTime_(0.0),
    autoRebalanceMigrationWallTime_(0.0),
    autoRebalanceWriteWallTime_(0.0),
    autoRebalancePostDiagWallTime_(0.0),
    autoRebalanceExplicitCommWallTime_(0.0),
    dlbAlpha_(1.0),
    adaptiveAlpha_(false),
    adaptiveAlphaMin_(0.5),
    adaptiveAlphaMax_(2.0),
    adaptiveAlphaGain_(0.04),
    adaptiveAlphaMaxStep_(0.08),
    adaptiveAlphaWorsenTol_(0.05),
    adaptiveAlphaLastImbalance_(GREAT),
    adaptiveAlphaLastStep_(0.0),
    asyncMigrationPending_(false),
    asyncRecvSize_(0),
    useNoAlltoall_(false),
    useFlatTransfer_(false),
    transferChunkBytes_(256*1024*1024),
    gatherCandidates_(false),
    overlapSizeExchange_(false),
    writeMode_("gathered"),
    // Wall-time trigger metric (§5.6)
    triggerWallTime_(false),
    lastStepWallTime_(0.0),
    stepFullWallCum_(0.0),
    lastStepWait_(0.0),
    lastSizeExchangeWall_(0.0),
    lastWaitWall_(0.0),
    lastUpcWall_(0.0),
    // Futile-rebalance back-off (§5.3)
    backOffEnabled_(false),
    backOffChangedRatio_(1e-3),
    backOffTol_(0.05),
    backOffFutileRuns_(2),
    backOffSkipChecks_(6),
    backOffEscalate_(1.5),
    backOffRemaining_(0),
    futileCount_(0),
    lastExecChangedRatio_(0.0),
    lastExecImbalance_(0.0),
    processorWriteDecompVersion_(-1),
    processorWriteMeshInstance_(word::null)
{}


dsmcReplicatedMesh::~dsmcReplicatedMesh()
{}


// ============================================================================
// computeCellOwnerScotch — Phase A
// ============================================================================

void dsmcReplicatedMesh::computeCellOwnerScotch()
{
    const label nCells = mesh_.nCells();
    const word decompMethod =
        mesh_.time().controlDict().lookupOrDefault<word>
        ("replicatedMeshDecompMethod", "block");

    if (decompMethod == "metis")
    {
        // Serial METIS on rank 0, broadcast result
        const label nIntFaces = mesh_.nInternalFaces();
        const labelList& faceOwner = mesh_.faceOwner();
        const labelList& faceNei = mesh_.faceNeighbour();

        // Build CSR adjacency
        labelList cellDeg(nCells, 0);
        for (label fi = 0; fi < nIntFaces; ++fi)
        {
            ++cellDeg[faceOwner[fi]];
            ++cellDeg[faceNei[fi]];
        }
        List<idx_t> xadj(nCells + 1);
        xadj[0] = 0;
        for (label i = 0; i < nCells; ++i)
            xadj[i+1] = xadj[i] + cellDeg[i];

        List<idx_t> adjncy(xadj[nCells]);
        labelList off(nCells, 0);
        for (label i = 0; i < nCells; ++i) off[i] = xadj[i];
        for (label fi = 0; fi < nIntFaces; ++fi)
        {
            const label own = faceOwner[fi];
            const label nei = faceNei[fi];
            adjncy[off[own]++] = nei;
            adjncy[off[nei]++] = own;
        }

        List<idx_t> part(nCells, 0);
        if (myRank_ == 0)
        {
            idx_t nvtxs = nCells;
            idx_t nparts = nProcs_;
            idx_t ncon = 1;
            idx_t edgecut = 0;
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions(options);
            options[METIS_OPTION_SEED] = 42;

            METIS_PartGraphKway
            (
                &nvtxs, &ncon, xadj.data(), adjncy.data(),
                nullptr, nullptr, nullptr,
                &nparts, nullptr, nullptr, options,
                &edgecut, part.data()
            );

            Info<< "Replicated mesh: METIS decomposition complete, "
                << nCells << " cells -> " << nProcs_ << " ranks" << endl;
        }
        MPI_Bcast(part.data(), nCells, MPI_INT, 0, MPI_COMM_WORLD);
        forAll(cellOwner_, i) cellOwner_[i] = label(part[i]);
    }
    else if (decompMethod == "scotch")
    {
        pointField cellCentres(nCells);
        forAll(mesh_.cells(), cellI) cellCentres[cellI] = mesh_.cellCentres()[cellI];
        dictionary decompDict;
        decompDict.set("method", "scotch");
        decompDict.set("numberOfSubdomains", nProcs_);
        dictionary scotchCoeffs;
        decompDict.set("scotchCoeffs", scotchCoeffs);
        scotchDecomp decomposer(decompDict);
        cellOwner_ = decomposer.decompose(mesh_, cellCentres);
        Info<< "Replicated mesh: Scotch decomposition complete, "
            << nCells << " cells -> " << nProcs_ << " ranks" << endl;
    }
    else if (decompMethod == "skew")
    {
        // Skewed: rank 0 gets 50% of cells; rest distributed among others
        const label nHalf = nCells / 2;
        for (label cellI = 0; cellI < nHalf; ++cellI)
            cellOwner_[cellI] = 0;
        for (label cellI = nHalf; cellI < nCells; ++cellI)
        {
            cellOwner_[cellI] = 1 +
                label((scalar(cellI - nHalf)/scalar(nCells - nHalf))*scalar(nProcs_ - 1));
            if (cellOwner_[cellI] >= nProcs_) cellOwner_[cellI] = nProcs_ - 1;
        }
        Info<< "Replicated mesh: skewed decomposition complete, "
            << nCells << " cells -> " << nProcs_ << " ranks"
            << " (rank 0 heavy)" << endl;
    }
    else
    {
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            cellOwner_[cellI] = label((scalar(cellI)/scalar(nCells))*scalar(nProcs_));
            if (cellOwner_[cellI] >= nProcs_) cellOwner_[cellI] = nProcs_ - 1;
        }
        Info<< "Replicated mesh: block decomposition complete, "
            << nCells << " cells -> " << nProcs_ << " ranks" << endl;
    }
}


// ============================================================================
// rebuildMyCells — Phase A
// ============================================================================

void dsmcReplicatedMesh::rebuildMyCells()
{
    const label nCells = mesh_.nCells();
    myCells_.clear();
    myCells_.setCapacity(nCells / nProcs_ + 1);

    label nT = 1;
    #ifdef _OPENMP
    if (cloud_.openmpEnabled())
    {
        nT = max(cloud_.ompNumThreads(), label(1));
    }
    #endif

    if (nT > 1 && nCells >= 2*nT)
    {
        // Two-pass fill: parallel counts, serial prefix, parallel fill.
        // Preserves ascending global cell order of the serial path.
        List<label> threadBegin(nT + 1);
        for (label tI = 0; tI <= nT; ++tI)
        {
            threadBegin[tI] = label
            (
                (static_cast<long long>(tI)
               * static_cast<long long>(nCells))
              / static_cast<long long>(nT)
            );
        }

        List<label> threadCounts(nT, 0);
        #pragma omp parallel num_threads(nT)
        {
            #ifdef _OPENMP
            const label tid = omp_get_thread_num();
            #else
            const label tid = 0;
            #endif
            label count = 0;
            for (label cellI = threadBegin[tid]; cellI < threadBegin[tid + 1]; ++cellI)
            {
                if (cellOwner_[cellI] == myRank_)
                {
                    ++count;
                }
            }
            threadCounts[tid] = count;
        }

        List<label> threadStart(nT + 1, 0);
        for (label tI = 0; tI < nT; ++tI)
        {
            threadStart[tI + 1] = threadStart[tI] + threadCounts[tI];
        }
        myCells_.setSize(threadStart[nT]);

        #pragma omp parallel num_threads(nT)
        {
            #ifdef _OPENMP
            const label tid = omp_get_thread_num();
            #else
            const label tid = 0;
            #endif
            label idx = threadStart[tid];
            for (label cellI = threadBegin[tid]; cellI < threadBegin[tid + 1]; ++cellI)
            {
                if (cellOwner_[cellI] == myRank_)
                {
                    myCells_[idx++] = cellI;
                }
            }
        }
    }
    else
    {
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            if (cellOwner_[cellI] == myRank_) myCells_.append(cellI);
        }
    }

    Info<< "Replicated mesh: rank " << myRank_ << " owns "
        << myCells_.size() << " / " << nCells << " cells" << endl;
}


// ============================================================================
// validateCellOwnerMap -- validate a DLB owner-map update before migration
// ============================================================================

void dsmcReplicatedMesh::validateCellOwnerMap(const char* context) const
{
    const label nCells = mesh_.nCells();

    if (cellOwner_.size() != nCells)
    {
        FatalErrorInFunction
            << context << ": cellOwner size " << cellOwner_.size()
            << " does not match mesh cell count " << nCells
            << exit(FatalError);
    }

    label nOwnedByMap = 0;
    forAll(cellOwner_, cellI)
    {
        const label owner = cellOwner_[cellI];
        if (owner < 0 || owner >= nProcs_)
        {
            FatalErrorInFunction
                << context << ": invalid owner " << owner
                << " for global cell " << cellI
                << " (valid range [0, " << nProcs_ - 1 << "])"
                << exit(FatalError);
        }

        if (owner == myRank_)
        {
            ++nOwnedByMap;
        }
    }

    // Duplicate detection by membership bitmap.  The previous HashSet<label>
    // relied on Foam::Hash<label> (identity hash) with power-of-two bucketing
    // (HashTableCore::canonicalSize): for rank-dependent ownership
    // distributions this collapsed many keys into few buckets, turning
    // insert/found into O(myCells^2) chain walks -- the dominant DLB cost on
    // the 13M-cell mesh (measured: one rank grinding 15+ min inside the
    // chain-walk loop, bjm8 jobs 4746365/4746420).  A 13 MB bitmap gives
    // guaranteed O(1) membership regardless of key distribution; the transient
    // is negligible against a ~90 GB rank footprint.
    List<bool> listedOwned(nCells, false);
    forAll(myCells_, localI)
    {
        const label cellI = myCells_[localI];
        if (cellI < 0 || cellI >= nCells)
        {
            FatalErrorInFunction
                << context << ": myCells contains invalid global cell "
                << cellI << exit(FatalError);
        }
        if (cellOwner_[cellI] != myRank_)
        {
            FatalErrorInFunction
                << context << ": myCells contains global cell " << cellI
                << " owned by rank " << cellOwner_[cellI]
                << ", not local rank " << myRank_
                << exit(FatalError);
        }
        if (listedOwned[cellI])
        {
            FatalErrorInFunction
                << context << ": myCells contains duplicate global cell "
                << cellI << exit(FatalError);
        }
        listedOwned[cellI] = true;
    }

    if (myCells_.size() != nOwnedByMap)
    {
        FatalErrorInFunction
            << context << ": myCells has " << myCells_.size()
            << " cells, but cellOwner assigns " << nOwnedByMap
            << " cells to rank " << myRank_
            << exit(FatalError);
    }

    forAll(cellOwner_, cellI)
    {
        if (cellOwner_[cellI] == myRank_ && !listedOwned[cellI])
        {
            FatalErrorInFunction
                << context << ": owned global cell " << cellI
                << " is absent from myCells" << exit(FatalError);
        }
    }

    label totalOwned = nOwnedByMap;
    MPI_Allreduce(MPI_IN_PLACE, &totalOwned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (totalOwned != nCells)
    {
        FatalErrorInFunction
            << context << ": owner-map coverage is " << totalOwned
            << " cells across all ranks, expected " << nCells
            << exit(FatalError);
    }

    // Prove that every rank has exactly the same replicated owner map.
    // A 64-bit wrapping checksum per rank is compared via two MIN/MAX
    // allreduces.  Wrapping addition is associative and commutative, so the
    // checksum is rank-order independent.  This replaces the previous
    // per-chunk MIN/MAX allreduce over the full owner map (2 * nCells/256K
    // collectives; ~100 on a 13M-cell mesh), which was a measurable DLB
    // cost and an observed intermittent hcoll/UCX hang point on bjm8
    // (job 4738842 DLB #2, job 4744935 DLB #1: all ranks spinning inside
    // the collective, gdb stacks in validateCellOwnerMap).
    const std::uint64_t checksumMix = 0x9E3779B97F4A7C15ULL;
    std::uint64_t localChecksum = 0;
    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        localChecksum +=
            static_cast<std::uint64_t>(cellOwner_[cellI] + 1)
          * (checksumMix ^ static_cast<std::uint64_t>(cellI));
    }

    unsigned long long minChecksum = localChecksum;
    unsigned long long maxChecksum = localChecksum;
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        &minChecksum,
        1,
        MPI_UNSIGNED_LONG_LONG,
        MPI_MIN,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        &maxChecksum,
        1,
        MPI_UNSIGNED_LONG_LONG,
        MPI_MAX,
        MPI_COMM_WORLD
    );

    if (minChecksum != maxChecksum)
    {
        FatalErrorInFunction
            << context << ": replicated owner map differs between MPI ranks"
            // Ostream has no unsigned long long overload in OFv1706; the
            // long cast preserves the 64-bit pattern for diagnostics.
            << " (checksum min=" << static_cast<long>(minChecksum)
            << ", max=" << static_cast<long>(maxChecksum) << ")"
            << exit(FatalError);
    }

    if (myRank_ == 0)
    {
        Info<< "Phase C DLB owner map validated: " << nCells
            << " cells, " << nProcs_ << " ranks" << endl;
    }
}


// ============================================================================
// initialize — Phase A
// ============================================================================

void dsmcReplicatedMesh::initialize()
{
    int mpiInit = 0;
    MPI_Initialized(&mpiInit);
    if (!mpiInit)
    {
        int argc = 0;
        char** argv = nullptr;
        UPstream::init(argc, argv);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs_);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
    allParticleCounts_.setSize(nProcs_, 0);

    dlbAlpha_ = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBAlpha", 1.0);
    adaptiveAlpha_ = mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshDLBAdaptiveAlpha", false);
    adaptiveAlphaMin_ = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBAlphaMin", 0.5);
    adaptiveAlphaMax_ = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBAlphaMax", 2.0);
    adaptiveAlphaGain_ = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBAlphaGain", 0.04);
    adaptiveAlphaMaxStep_ = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBAlphaMaxStep", 0.08);
    adaptiveAlphaWorsenTol_ = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBAlphaWorsenTol", 0.05);
    if (adaptiveAlphaMin_ > adaptiveAlphaMax_)
    {
        const scalar tmp = adaptiveAlphaMin_;
        adaptiveAlphaMin_ = adaptiveAlphaMax_;
        adaptiveAlphaMax_ = tmp;
    }
    dlbAlpha_ = max(adaptiveAlphaMin_, min(adaptiveAlphaMax_, dlbAlpha_));


    if (nProcs_ < 2)
    {
        cellOwner_ = 0;
        myCells_.setSize(mesh_.nCells());
        forAll(mesh_.cells(), cellI) myCells_.append(cellI);
        Info<< "Replicated mesh: single-rank run, all "
            << mesh_.nCells() << " cells owned" << endl;
        active_ = false;
        return;
    }

    computeCellOwnerScotch();
    rebuildMyCells();

    // Build local mesh view (owned cells + 1-ring halo)
    localMesh_.build(myCells_);

    // Fixed to 1 — every-step migration is required for correctness.
    // Larger intervals cause particles to drift into non-owned cells
    // where they are invisible to collision (up to 66% collision loss
    // at interval=10).  ParDSMC3D and standard OpenFOAM MPI both
    // migrate every step.
    migrateInterval_ = mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshMigrateInterval", 1);
    stepCounter_ = 0;

    // Phase B: manual rebalance steps
    rebalanceSteps_.clear();
    if (mesh_.time().controlDict().found("replicatedMeshRebalanceSteps"))
    {
        rebalanceSteps_ = labelList
        (
            mesh_.time().controlDict().lookup("replicatedMeshRebalanceSteps")
        );
    }

    // ---- Phase C config ----
    autoDLBEnabled_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshAutoDLB", false);

    // ---- Wall-time trigger metric (§5.6) ----
    triggerWallTime_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshDLBTriggerWallTime", false);

    // ---- Futile-rebalance back-off (§5.3) ----
    backOffEnabled_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshDLBBackOff", false);
    backOffChangedRatio_ =
        mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBBackOffChangedRatio", 1e-3);
    backOffTol_ =
        mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBBackOffTol", 0.05);
    backOffFutileRuns_ =
        max(label(1), mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBBackOffFutileRuns", 2));
    backOffSkipChecks_ =
        max(label(1), mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBBackOffSkipChecks", 6));
    backOffEscalate_ =
        mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBBackOffEscalate", 1.5);
    if (backOffEnabled_)
    {
        Info<< "Replicated mesh: DLB back-off enabled (changedRatio<"
            << backOffChangedRatio_ << " tol=" << backOffTol_
            << " futileRuns=" << backOffFutileRuns_
            << " skipChecks=" << backOffSkipChecks_
            << " escalate=" << backOffEscalate_ << "x)" << nl
            << "Replicated mesh: DLB wall trigger metric = "
            << (triggerWallTime_ ? "per-step full wall" : "CPU (legacy)")
            << endl;
    }

    // ---- No-Alltoall async migration ----
    useNoAlltoall_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshNoAlltoall", false);
    if (useNoAlltoall_)
    {
        prevRecvSizes_.setSize(nProcs_, 0);
        perPeerRecvCapacity_.setSize(nProcs_, 0);
        perPeerRecvBuf_.setSize(nProcs_);
        Info<< "Replicated mesh: no-Alltoall async migration enabled" << endl;
    }

    useFlatTransfer_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshFlatTransfer", true);
    Info<< "Replicated mesh: flatTransfer=" << useFlatTransfer_
        << " found=" << mesh_.time().controlDict().found("replicatedMeshFlatTransfer")
        << endl;
    if (useFlatTransfer_)
    {
        Info<< "Replicated mesh: flat POD transfer enabled" << endl;
    }

    const label requestedTransferChunkMB =
        mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshTransferChunkMB", 256);
    const label transferChunkMB =
        max(label(1), min(requestedTransferChunkMB, label(512)));
    if (transferChunkMB != requestedTransferChunkMB)
    {
        WarningInFunction
            << "replicatedMeshTransferChunkMB="
            << requestedTransferChunkMB
            << " is outside [1, 512]; using "
            << transferChunkMB << " MiB" << endl;
    }
    transferChunkBytes_ = transferChunkMB*1024*1024;
    Info<< "Replicated mesh: transfer chunk limit="
        << transferChunkBytes_/(1024*1024) << " MiB"
        << " (" << transferChunkBytes_ << " bytes)" << endl;

    gatherCandidates_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshGatherCandidates", false);
    Info<< "Replicated mesh: gatherCandidates=" << gatherCandidates_
        << endl;

    overlapSizeExchange_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshOverlapSizeExchange", false);
    Info<< "Replicated mesh: overlapSizeExchange=" << overlapSizeExchange_
        << endl;

    writeMode_ =
        mesh_.time().controlDict().lookupOrDefault<word>
        ("replicatedMeshWriteMode", "gathered");
    if (writeMode_ != "gathered" && writeMode_ != "processor")
    {
        FatalErrorInFunction
            << "Unknown replicatedMeshWriteMode " << writeMode_
            << ". Valid values are gathered and processor."
            << exit(FatalError);
    }
    Info<< "Replicated mesh: writeMode=" << writeMode_ << endl;

    if (autoDLBEnabled_)
    {
        imbalanceThreshold_ = mesh_.time().controlDict()
            .lookupOrDefault<scalar>("replicatedMeshDLBImbalanceThreshold", 1.15);
        dlbSteps_ = mesh_.time().controlDict()
            .lookupOrDefault<label>("replicatedMeshDLBSteps", 50);
        if (dlbSteps_ < 10) dlbSteps_ = 10;
        sarSteps_ = mesh_.time().controlDict()
            .lookupOrDefault<label>("replicatedMeshSARSteps", 10);
        if (sarSteps_ < 1) sarSteps_ = 1;
        autoDLBTriggerMode_ = mesh_.time().controlDict()
            .lookupOrDefault<word>("replicatedMeshDLBTriggerMode", "cumulative");
        forcedDLBSteps_ = mesh_.time().controlDict()
            .lookupOrDefault<labelList>("replicatedMeshDLBForceSteps", labelList());
        const word checkCollective = mesh_.time().controlDict()
            .lookupOrDefault<word>("replicatedMeshDLBCheckCollective", "allgather");
        const bool skipPostDiag = mesh_.time().controlDict()
            .lookupOrDefault<bool>("replicatedMeshDLBSkipPostDiag", false);
        const label minGapSteps = mesh_.time().controlDict()
            .lookupOrDefault<label>("replicatedMeshDLBMinGapSteps", 0);
        const bool particleGate = mesh_.time().controlDict()
            .lookupOrDefault<bool>("replicatedMeshDLBParticleGate", true);
        const scalar particleGateThreshold = mesh_.time().controlDict()
            .lookupOrDefault<scalar>
            (
                "replicatedMeshDLBParticleGateThreshold",
                imbalanceThreshold_
            );

        Info<< "Phase C auto DLB: enabled (ParMETIS AdaptiveRepart)"
            << ", imbalanceThreshold=" << imbalanceThreshold_
            << ", dlbSteps=" << dlbSteps_
            << ", sarSteps=" << sarSteps_
            << ", triggerMode=" << autoDLBTriggerMode_
            << ", minGapSteps=" << minGapSteps
            << ", checkCollective=" << checkCollective
            << ", particleGate=" << particleGate
            << ", particleGateThreshold=" << particleGateThreshold
            << ", skipPostDiag=" << skipPostDiag
            << ", adaptiveAlpha=" << adaptiveAlpha_;
        if (adaptiveAlpha_)
        {
            Info<< " alpha0=" << dlbAlpha_
                << " alphaRange=[" << adaptiveAlphaMin_ << ","
                << adaptiveAlphaMax_ << "]"
                << " alphaGain=" << adaptiveAlphaGain_
                << " alphaMaxStep=" << adaptiveAlphaMaxStep_;
        }
        Info<< endl;
        if (forcedDLBSteps_.size())
        {
            Info<< "Phase C auto DLB: forced trigger steps "
                << forcedDLBSteps_ << endl;
        }
    }

    if (writeMode_ == "processor")
    {
        restoreProcessorCheckpoint();
    }

    active_ = true;
    Info<< "Replicated mesh: initialized with " << nProcs_ << " MPI ranks"
        << ", migrate interval " << migrateInterval_ << endl;
}


// ============================================================================
// Restore a processor-write checkpoint selected before cloud evolution.
// ============================================================================

bool dsmcReplicatedMesh::restoreProcessorCheckpoint()
{
    if (writeMode_ != "processor")
    {
        return false;
    }

    const Time& runTime = mesh_.time();
    if
    (
        runTime.controlDict().lookupOrDefault<word>
        ("startFrom", "latestTime") != "latestTime"
    )
    {
        return false;
    }

    const word timeName = runTime.timeName();
    const fileName processorRoot = runTime.path();

    label nCheckpointProcs = 0;
    while
    (
        isDir
        (
            processorRoot
           /fileName
            (
                word("processor") + Foam::name(nCheckpointProcs)
            )
        )
    )
    {
        ++nCheckpointProcs;
    }

    if (!nCheckpointProcs)
    {
        return false;
    }

    if (nCheckpointProcs != nProcs_)
    {
        FatalErrorInFunction
            << "Processor checkpoint at time " << timeName
            << " contains " << nCheckpointProcs << " processor directories,"
            << " but the current run has " << nProcs_ << " MPI ranks."
            << " Processor-write restart requires the same MPI rank count."
            << exit(FatalError);
    }

    const fileName processorDir =
        processorRoot
       /fileName(word("processor") + Foam::name(myRank_));
    const fileName timePolyMeshDir =
        processorDir/timeName/polyMesh::meshSubDir;
    const fileName constantPolyMeshDir =
        processorDir/runTime.constant()/polyMesh::meshSubDir;

    word procMeshInstance;
    if (isFile(timePolyMeshDir/"cellProcAddressing"))
    {
        procMeshInstance = timeName;
    }
    else if (isFile(constantPolyMeshDir/"cellProcAddressing"))
    {
        procMeshInstance = runTime.constant();
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find cellProcAddressing for processor" << myRank_
            << " at checkpoint time " << timeName << exit(FatalError);
    }

    fileName processorCasePath
    (
        runTime.caseName()
       /fileName(word("processor") + Foam::name(myRank_))
    );
    Time processorDb
    (
        Time::controlDictName,
        runTime.rootPath(),
        processorCasePath,
        word("system"),
        word("constant")
    );
    processorDb.setTime(runTime);

    fvMesh procMesh
    (
        IOobject
        (
            mesh_.name(),
            procMeshInstance,
            processorDb,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellProcAddressing",
            procMesh.facesInstance(),
            procMesh.meshSubDir,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (!cellProcAddressing.size())
    {
        FatalErrorInFunction
            << "Processor" << myRank_ << " checkpoint has no cells"
            << " in cellProcAddressing at time " << timeName
            << exit(FatalError);
    }

    const label nCells = mesh_.nCells();
    labelList localOwner(nCells, -1);
    labelList localOwnerCount(nCells, 0);

    forAll(cellProcAddressing, procCellI)
    {
        const label globalCellI = cellProcAddressing[procCellI];
        if (globalCellI < 0 || globalCellI >= nCells)
        {
            FatalErrorInFunction
                << "Invalid global cell " << globalCellI
                << " in processor" << myRank_ << " cellProcAddressing"
                << exit(FatalError);
        }

        localOwner[globalCellI] = myRank_;
        localOwnerCount[globalCellI] = 1;
    }

    labelList ownerCount(nCells, 0);
    cellOwner_.setSize(nCells, -1);
    MPI_Allreduce
    (
        localOwner.data(),
        cellOwner_.data(),
        nCells,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        localOwnerCount.data(),
        ownerCount.data(),
        nCells,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD
    );

    forAll(ownerCount, cellI)
    {
        if (ownerCount[cellI] != 1 || cellOwner_[cellI] < 0)
        {
            FatalErrorInFunction
                << "Processor checkpoint does not provide exactly one owner"
                << " for global cell " << cellI
                << " (owner count=" << ownerCount[cellI]
                << ", owner=" << cellOwner_[cellI] << ")"
                << exit(FatalError);
        }
    }

    rebuildMyCells();
    localMesh_.build(myCells_);

    volScalarField processorSigma
    (
        IOobject
        (
            cloud_.name() + "SigmaTcRMax",
            timeName,
            processorDb,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        procMesh
    );

    if (processorSigma.size() != cellProcAddressing.size())
    {
        FatalErrorInFunction
            << "Processor" << myRank_ << " sigmaTcRMax size "
            << processorSigma.size() << " does not match local cell count "
            << cellProcAddressing.size() << exit(FatalError);
    }

    scalarField localSigma(nCells, 0.0);
    scalarField globalSigma(nCells, 0.0);
    forAll(cellProcAddressing, procCellI)
    {
        localSigma[cellProcAddressing[procCellI]] = processorSigma[procCellI];
    }

    MPI_Allreduce
    (
        localSigma.data(),
        globalSigma.data(),
        nCells,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );

    forAll(cloud_.sigmaTcRMax(), cellI)
    {
        cloud_.sigmaTcRMax()[cellI] = globalSigma[cellI];
    }
    cloud_.sigmaTcRMax().correctBoundaryConditions();

    passiveParticleCloud positions
    (
        procMesh,
        cloud_.name(),
        IDLList<passiveParticle>()
    );
    const fileName positionsPath =
        processorDir/timeName/"lagrangian"/cloud_.name()/"positions";
    if (isFile(positionsPath))
    {
        IOPosition<Cloud<passiveParticle>> positionIO(positions);
        positionIO.readData(positions, false);
        positionIO.close();
    }

    const label nParcels = positions.size();
    if (nParcels)
    {
        IOField<vector> U
        (
            positions.fieldIOobject("U", IOobject::MUST_READ)
        );
        IOField<scalar> RWF
        (
            positions.fieldIOobject("radialWeight", IOobject::READ_IF_PRESENT),
            scalarField(nParcels, 1.0)
        );
        IOField<scalar> ERot
        (
            positions.fieldIOobject("ERot", IOobject::READ_IF_PRESENT),
            scalarField(nParcels, 0.0)
        );
        IOField<label> ELevel
        (
            positions.fieldIOobject("ELevel", IOobject::READ_IF_PRESENT),
            labelField(nParcels, 0)
        );
        IOField<label> typeId
        (
            positions.fieldIOobject("typeId", IOobject::MUST_READ)
        );
        IOField<label> newParcel
        (
            positions.fieldIOobject("newParcel", IOobject::MUST_READ)
        );
        IOField<label> classification
        (
            positions.fieldIOobject("classification", IOobject::MUST_READ)
        );

        IOField<label> stuckToWall
        (
            positions.fieldIOobject("stuckToWall", IOobject::READ_IF_PRESENT),
            labelField(nParcels, 0)
        );
        IOField<scalarField> wallTemperature
        (
            positions.fieldIOobject
            (
                "wallTemperature",
                IOobject::READ_IF_PRESENT
            )
        );
        if (wallTemperature.size() != nParcels)
        {
            wallTemperature.setSize(nParcels);
            forAll(wallTemperature, i)
            {
                wallTemperature[i] = scalarField(4, 0.0);
            }
        }

        IOField<vectorField> wallVectors
        (
            positions.fieldIOobject
            (
                "wallVectors",
                IOobject::READ_IF_PRESENT
            )
        );
        if (wallVectors.size() != nParcels)
        {
            wallVectors.setSize(nParcels);
            forAll(wallVectors, i)
            {
                wallVectors[i] = vectorField(4, vector::zero);
            }
        }

        IOField<label> isTracked
        (
            positions.fieldIOobject("isTracked", IOobject::READ_IF_PRESENT),
            labelField(nParcels, 0)
        );
        IOField<label> inPatchId
        (
            positions.fieldIOobject("inPatchId", IOobject::READ_IF_PRESENT),
            labelField(nParcels, -1)
        );
        IOField<scalar> tracerInitialTime
        (
            positions.fieldIOobject
            (
                "tracerInitialTime",
                IOobject::READ_IF_PRESENT
            ),
            scalarField(nParcels, 0.0)
        );
        IOField<vector> tracerInitialPosition
        (
            positions.fieldIOobject
            (
                "tracerInitialPosition",
                IOobject::READ_IF_PRESENT
            ),
            vectorField(nParcels, vector::zero)
        );
        IOField<vector> tracerCurrentPosition
        (
            positions.fieldIOobject
            (
                "tracerCurrentPosition",
                IOobject::READ_IF_PRESENT
            ),
            vectorField(nParcels, vector::zero)
        );
        IOField<vector> tracerDistanceTravelled
        (
            positions.fieldIOobject
            (
                "tracerDistanceTravelled",
                IOobject::READ_IF_PRESENT
            ),
            vectorField(nParcels, vector::zero)
        );

        IOField<labelField> vibLevel
        (
            positions.fieldIOobject("vibLevel", IOobject::READ_IF_PRESENT)
        );
        if (vibLevel.size() != nParcels)
        {
            vibLevel.setSize(nParcels);
            forAll(vibLevel, i)
            {
                vibLevel[i].setSize(0);
            }
        }

        labelField defaultOrigProc(nParcels, myRank_);
        labelField defaultOrigId(nParcels, 0);
        forAll(defaultOrigId, i)
        {
            defaultOrigId[i] = i;
        }
        IOField<label> origProcId
        (
            positions.fieldIOobject("origProcId", IOobject::READ_IF_PRESENT),
            defaultOrigProc
        );
        IOField<label> origId
        (
            positions.fieldIOobject("origId", IOobject::READ_IF_PRESENT),
            defaultOrigId
        );

        if
        (
            U.size() != nParcels
         || RWF.size() != nParcels
         || ERot.size() != nParcels
         || ELevel.size() != nParcels
         || typeId.size() != nParcels
         || newParcel.size() != nParcels
         || classification.size() != nParcels
         || stuckToWall.size() != nParcels
         || wallTemperature.size() != nParcels
         || wallVectors.size() != nParcels
         || isTracked.size() != nParcels
         || inPatchId.size() != nParcels
         || tracerInitialTime.size() != nParcels
         || tracerInitialPosition.size() != nParcels
         || tracerCurrentPosition.size() != nParcels
         || tracerDistanceTravelled.size() != nParcels
         || vibLevel.size() != nParcels
         || origProcId.size() != nParcels
         || origId.size() != nParcels
        )
        {
            FatalErrorInFunction
                << "Processor" << myRank_
                << " checkpoint field sizes do not match positions size "
                << nParcels << exit(FatalError);
        }

        if (cloud_.size())
        {
            cloud_.clear();
        }

        label parcelI = 0;
        forAllConstIter(Cloud<passiveParticle>, positions, iter)
        {
            const passiveParticle& position = iter();
            const label localCellI = position.cell();
            if (localCellI < 0 || localCellI >= cellProcAddressing.size())
            {
                FatalErrorInFunction
                    << "Invalid local cell " << localCellI
                    << " in processor" << myRank_ << " positions"
                    << exit(FatalError);
            }

            const label globalCellI = cellProcAddressing[localCellI];
            label tetFaceI = -1;
            label tetPtI = -1;
            mesh_.findTetFacePt
            (
                globalCellI,
                position.position(),
                tetFaceI,
                tetPtI
            );
            if (tetFaceI < 0 || tetPtI < 0)
            {
                FatalErrorInFunction
                    << "Cannot locate tetrahedron for restored parcel "
                    << parcelI << " in global cell " << globalCellI
                    << " at position " << position.position()
                    << exit(FatalError);
            }

            dsmcParcel* parcel = new dsmcParcel
            (
                mesh_,
                position.position(),
                U[parcelI],
                RWF[parcelI],
                ERot[parcelI],
                ELevel[parcelI],
                globalCellI,
                tetFaceI,
                tetPtI,
                typeId[parcelI],
                newParcel[parcelI],
                classification[parcelI],
                vibLevel[parcelI]
            );

            parcel->origProc() = origProcId[parcelI];
            parcel->origId() = origId[parcelI];

            if (stuckToWall[parcelI])
            {
                parcel->setStuck
                (
                    wallTemperature[parcelI],
                    wallVectors[parcelI]
                );
            }

            if (isTracked[parcelI])
            {
                parcel->setTracked
                (
                    true,
                    inPatchId[parcelI],
                    tracerInitialTime[parcelI],
                    tracerInitialPosition[parcelI],
                    tracerDistanceTravelled[parcelI]
                );
                parcel->tracked().updateCurrentPosition
                (
                    tracerCurrentPosition[parcelI]
                );
            }

            cloud_.addParticle(parcel);
            ++parcelI;
        }
    }
    else if (cloud_.size())
    {
        cloud_.clear();
    }

    IOdictionary uniformPropsDict
    (
        IOobject
        (
            Cloud<dsmcParcel>::cloudPropertiesName,
            timeName,
            "uniform"/cloud::prefix/cloud_.name(),
            processorDb,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    const word procName("processor" + Foam::name(myRank_));
    if
    (
        !uniformPropsDict.found(procName)
     || !uniformPropsDict.subDict(procName).found("particleCount")
    )
    {
        FatalErrorInFunction
            << "Processor checkpoint cloudProperties has no particleCount"
            << " for " << procName << exit(FatalError);
    }
    uniformPropsDict.subDict(procName).lookup("particleCount")
        >> dsmcParcel::particleCount_;

    cloud_.reBuildCellOccupancy();
    localParticleCount_ = cloud_.size();
    MPI_Allgather
    (
        &localParticleCount_,
        1,
        MPI_INT,
        allParticleCounts_.data(),
        1,
        MPI_INT,
        MPI_COMM_WORLD
    );

    Info<< "Replicated mesh: rank " << myRank_
        << " restored processor checkpoint " << timeName
        << " with " << localParticleCount_ << " parcels and "
        << myCells_.size() << " owned cells" << endl;

    return true;
}


// ============================================================================
// reassignByParMetisAdaptiveRepart — Phase C: ParMETIS incremental repartition
// ============================================================================

label dsmcReplicatedMesh::reassignByParMetisAdaptiveRepart()
{
    const label nCells = mesh_.nCells();
    const label nIntFaces = mesh_.nInternalFaces();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNei = mesh_.faceNeighbour();
    const bool dlbProfile = mesh_.time().controlDict().lookupOrDefault<bool>
    (
        "replicatedMeshDLBProfile",
        false
    );

    // OMP thread count for DLB mesh-side bookkeeping loops.  Mirrors the
    // migration-path convention: controlDict openmpThreads, degrading to
    // serial when OpenMP is disabled or compiled out.  The if(nT > 1)
    // clause keeps a single loop body for both paths.
    label nT = 1;
    #ifdef _OPENMP
    if (cloud_.openmpEnabled())
    {
        nT = max(cloud_.ompNumThreads(), label(1));
    }
    #endif

    // ---- Distributed graph: each rank holds nCells/nProcs_ vertices --------
    const label localN = nCells / nProcs_;
    const label myStart = myRank_ * localN;
    const label myEnd = (myRank_ == nProcs_ - 1) ? nCells : (myRank_ + 1) * localN;
    const label myN = myEnd - myStart;

    // vtxdist: vertex distribution
    List<idx_t> vtxdist(nProcs_ + 1);
    for (label i = 0; i <= nProcs_; ++i)
        vtxdist[i] = (i < nProcs_) ? i * localN : nCells;

    // Build LOCAL xadj/adjncy for my portion of the graph
    // First count local degrees
    labelList localDeg(myN, 0);
    for (label fi = 0; fi < nIntFaces; ++fi)
    {
        const label own = faceOwner[fi];
        const label nei = faceNei[fi];
        if (own >= myStart && own < myEnd) ++localDeg[own - myStart];
        if (nei >= myStart && nei < myEnd) ++localDeg[nei - myStart];
    }

    List<idx_t> xadj(myN + 1);
    xadj[0] = 0;
    for (label i = 0; i < myN; ++i)
        xadj[i+1] = xadj[i] + localDeg[i];

    List<idx_t> adjncy(xadj[myN]);
    labelList off(myN, 0);
    for (label i = 0; i < myN; ++i) off[i] = xadj[i];

    for (label fi = 0; fi < nIntFaces; ++fi)
    {
        const label own = faceOwner[fi];
        const label nei = faceNei[fi];
        if (own >= myStart && own < myEnd)
            adjncy[off[own - myStart]++] = nei;
        if (nei >= myStart && nei < myEnd)
            adjncy[off[nei - myStart]++] = own;
    }

    // ParMETIS owns vertices by vtxdist, not by the current DSMC cellOwner_.
    // After a DLB step those two layouts differ, so local cellOccupancy() is
    // not valid for all vertices in this rank's vtxdist range.  Build a
    // replicated global particle-count field before constructing weights.
    List<idx_t> localCellParticles(nCells, 0);
    List<idx_t> globalCellParticles(nCells, 0);
    // cellOccupancy() lazily materialises: resolve the reference serially
    // before entering any parallel region.
    const DynamicList<DynamicList<dsmcParcel*>>& occCells =
        cloud_.cellOccupancy();
    const label nOcc = min(label(occCells.size()), nCells);
    #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static)
    for (label cellI = 0; cellI < nOcc; ++cellI)
    {
        localCellParticles[cellI] = idx_t(occCells[cellI].size());
    }
    if (dlbProfile && myRank_ == 0)
    {
        Info<< "Phase C ParMETIS profile: entering particle-count allreduce"
            << endl;
    }
    const auto tParticleAllreduce0 = std::chrono::steady_clock::now();
    MPI_Allreduce
    (
        localCellParticles.data(),
        globalCellParticles.data(),
        nCells,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD
    );
    addAutoRebalanceCommTime
    (
        std::chrono::duration<scalar>
        (
            std::chrono::steady_clock::now() - tParticleAllreduce0
        ).count()
    );
    if (dlbProfile && myRank_ == 0)
    {
        Info<< "Phase C ParMETIS profile: particle-count allreduce complete"
            << endl;
    }

    const idx_t maxParMetisWeight = std::numeric_limits<idx_t>::max()/4;
    auto positiveWeight = [&](const scalar value) -> idx_t
    {
        if (!(value > scalar(1)))
        {
            return idx_t(1);
        }
        if (value > scalar(maxParMetisWeight))
        {
            return maxParMetisWeight;
        }
        return max(idx_t(1), idx_t(value + 0.5));
    };

    if (myRank_ == 0)
    {
        idx_t maxCellParticles = 0;
        label activeCells = 0;
        label totalParticles = 0;
        #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static) \
            reduction(+:activeCells) reduction(+:totalParticles) \
            reduction(max:maxCellParticles)
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            const idx_t nPart = globalCellParticles[cellI];
            if (nPart > 0)
            {
                ++activeCells;
                totalParticles += label(nPart);
                if (nPart > maxCellParticles)
                {
                    maxCellParticles = nPart;
                }
            }
        }
        Info<< "Phase C ParMETIS weights: particles=" << totalParticles
            << " activeCells=" << activeCells
            << " maxCellParticles=" << maxCellParticles << endl;
    }

    // ---- Vertex weights: N^alpha (particle count, compressed) --------------
    const scalar moveTime = cloud_.evolveMoveWallTime();
    const scalar collTime = cloud_.evolveCollisionWallTime();

    idx_t ncon = 1;
    idx_t wgtflag = 2;  // vertex weights only

    List<idx_t> vwgt(myN * ncon, 1);
    #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static)
    for (label i = 0; i < myN; ++i)
    {
        const label gi = myStart + i;
        const idx_t nPart = globalCellParticles[gi];

        const scalar w = (nPart > 1)
            ? std::pow(scalar(nPart), dlbAlpha_) : scalar(nPart);
        vwgt[i] = positiveWeight(w);
    }

    // Weight = sum of particle counts of the two cells sharing the face.
    List<idx_t> adjwgt(xadj[myN], 1);
    {
        labelList localOff(myN, 0);
        for (label i = 0; i < myN; ++i) localOff[i] = xadj[i];

        for (label fi = 0; fi < nIntFaces; ++fi)
        {
            const label own = faceOwner[fi];
            const label nei = faceNei[fi];
            // Particle count on both sides of this face
            const idx_t nOwn = globalCellParticles[own];
            const idx_t nNei = globalCellParticles[nei];
            const idx_t ew = positiveWeight(scalar(nOwn) + scalar(nNei));

            if (own >= myStart && own < myEnd)
                adjwgt[localOff[own - myStart]++] = ew;
            if (nei >= myStart && nei < myEnd)
                adjwgt[localOff[nei - myStart]++] = ew;
        }
    }

    // ---- Current partition (LOCAL) -----------------------------------------
    List<idx_t> part(myN);
    #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static)
    for (label i = 0; i < myN; ++i)
        part[i] = cellOwner_[myStart + i];

    // ---- ParMETIS parameters -----------------------------------------------
    idx_t numflag = 0;  // C-style numbering
    idx_t nparts = nProcs_;
    List<real_t> tpwgts(ncon * nparts, real_t(1.0) / real_t(nparts));
    const real_t ubvecVal = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBUbvec", 1.05);
    List<real_t> ubvec(ncon, ubvecVal);

    Info<< "Phase C ParMETIS: alpha=" << dlbAlpha_ << " ncon=" << ncon
        << " ubvec=[";
    for (label c = 0; c < ncon; ++c)
        Info<< (c > 0 ? "," : "") << ubvec[c];
    Info<< "] moveRatio=" << (moveTime / max(moveTime + collTime, SMALL))
        << " (moveT=" << moveTime << "s, collT=" << collTime << "s)" << endl;

    real_t itr = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBItr", 100.0);
    idx_t options[4] =
    {
        1,                          // use custom options
        0,                          // no ParMETIS debug output
        42,                         // random seed
        PARMETIS_PSR_UNCOUPLED     // use the current part[] as the input partition
    };
    idx_t edgecut = 0;
    // vsize: 2-rank uses N (conservative), 4+ rank uses N^vsExp (vsExp=0 means free migration)
    const scalar vsExp = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBVsizeExp", 0.5);
    List<idx_t> vsize(myN, 1);
    if (nProcs_ <= 2)
    {
        #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static)
        for (label i = 0; i < myN; ++i)
        {
            const label gi = myStart + i;
            const idx_t nPart = globalCellParticles[gi];
            vsize[i] = positiveWeight(scalar(nPart));
        }
    }
    else if (vsExp > SMALL)
    {
        #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static)
        for (label i = 0; i < myN; ++i)
        {
            const label gi = myStart + i;
            const idx_t nPart = globalCellParticles[gi];
            if (nPart > 1)
            {
                vsize[i] = positiveWeight(std::pow(scalar(nPart), vsExp));
            }
        }
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    Info<< "Phase C ParMETIS: AdaptiveRepart (" << nCells << " cells, "
        << ncon << " constraints, ubvec=" << ubvec[0] << ")" << endl;

    if (dlbProfile && myRank_ == 0)
    {
        Info<< "Phase C ParMETIS profile: entering AdaptiveRepart" << endl;
    }
    ParMETIS_V3_AdaptiveRepart
    (
        vtxdist.data(), xadj.data(), adjncy.data(),
        vwgt.data(), vsize.data(), nullptr,
        &wgtflag, &numflag, &ncon, &nparts,
        tpwgts.data(), ubvec.data(), &itr,
        options, &edgecut, part.data(), &comm
    );
    if (dlbProfile && myRank_ == 0)
    {
        Info<< "Phase C ParMETIS profile: AdaptiveRepart complete" << endl;
    }

    // ---- Gather new partition from all ranks --------------------------------
    // Each rank has its LOCAL portion of the new partition.
    // Allgather to reconstruct the full cellOwner_.
    List<idx_t> fullPart(nCells, 0);
    List<int> recvCounts(nProcs_);
    List<int> displs(nProcs_);
    for (label i = 0; i < nProcs_; ++i)
    {
        displs[i] = (i < nProcs_) ? i * localN : 0;
        recvCounts[i] = (i == nProcs_ - 1) ? (nCells - i * localN) : localN;
    }
    if (dlbProfile && myRank_ == 0)
    {
        Info<< "Phase C ParMETIS profile: entering partition allgatherv"
            << endl;
    }
    const auto tPartitionAllgather0 = std::chrono::steady_clock::now();
    MPI_Allgatherv(part.data(), myN, MPI_INT,
                   fullPart.data(), recvCounts.data(), displs.data(),
                   MPI_INT, MPI_COMM_WORLD);
    addAutoRebalanceCommTime
    (
        std::chrono::duration<scalar>
        (
            std::chrono::steady_clock::now() - tPartitionAllgather0
        ).count()
    );
    if (dlbProfile && myRank_ == 0)
    {
        Info<< "Phase C ParMETIS profile: partition allgatherv complete"
            << endl;
    }

    const bool remapEnabled = mesh_.time().controlDict().lookupOrDefault<bool>
    (
        "replicatedMeshDLBRemap",
        false
    );
    const word remapMode = mesh_.time().controlDict().lookupOrDefault<word>
    (
        "replicatedMeshDLBRemapMode",
        "greedyOverlap"
    );

    label changedBeforeRemap = 0;
    #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static) \
        reduction(+:changedBeforeRemap)
    for (label i = 0; i < nCells; ++i)
    {
        if (label(fullPart[i]) != cellOwner_[i])
        {
            ++changedBeforeRemap;
        }
    }

    label remapSavedCells = 0;

    if (remapEnabled && remapMode == "greedyOverlap")
    {
        // Per-thread flattened nProcs x nProcs overlap matrices, merged
        // serially (M1 two-pass pattern) to avoid write conflicts on the
        // shared matrix entries.
        List<labelList> threadOverlap(nT);
        for (label tI = 0; tI < nT; ++tI)
        {
            threadOverlap[tI].setSize(nProcs_*nProcs_, 0);
        }

        #pragma omp parallel num_threads(nT) if (nT > 1)
        {
            #ifdef _OPENMP
            const label tid = omp_get_thread_num();
            #else
            const label tid = 0;
            #endif
            labelList& myOverlap = threadOverlap[tid];

            #pragma omp for schedule(static)
            for (label i = 0; i < nCells; ++i)
            {
                const label newPart = label(fullPart[i]);
                const label oldRank = cellOwner_[i];
                if
                (
                    newPart >= 0 && newPart < nProcs_
                 && oldRank >= 0 && oldRank < nProcs_
                )
                {
                    ++myOverlap[newPart*nProcs_ + oldRank];
                }
            }
        }

        List<labelList> overlap(nProcs_);
        for (label newPart = 0; newPart < nProcs_; ++newPart)
        {
            overlap[newPart].setSize(nProcs_, 0);
        }
        for (label tI = 0; tI < nT; ++tI)
        {
            const labelList& myOverlap = threadOverlap[tI];
            for (label newPart = 0; newPart < nProcs_; ++newPart)
            {
                for (label oldRank = 0; oldRank < nProcs_; ++oldRank)
                {
                    overlap[newPart][oldRank] +=
                        myOverlap[newPart*nProcs_ + oldRank];
                }
            }
        }

        labelList remap(nProcs_, -1);
        labelList newPartUsed(nProcs_, 0);
        labelList oldRankUsed(nProcs_, 0);

        for (label matchI = 0; matchI < nProcs_; ++matchI)
        {
            label bestNewPart = -1;
            label bestOldRank = -1;
            label bestOverlap = -1;

            for (label newPart = 0; newPart < nProcs_; ++newPart)
            {
                if (newPartUsed[newPart]) continue;

                for (label oldRank = 0; oldRank < nProcs_; ++oldRank)
                {
                    if (oldRankUsed[oldRank]) continue;

                    const label value = overlap[newPart][oldRank];
                    if
                    (
                        value > bestOverlap
                     || (
                            value == bestOverlap
                         && (
                                bestNewPart < 0
                             || newPart < bestNewPart
                             || (
                                    newPart == bestNewPart
                                 && oldRank < bestOldRank
                                )
                            )
                        )
                    )
                    {
                        bestOverlap = value;
                        bestNewPart = newPart;
                        bestOldRank = oldRank;
                    }
                }
            }

            if (bestNewPart >= 0 && bestOldRank >= 0)
            {
                remap[bestNewPart] = bestOldRank;
                newPartUsed[bestNewPart] = 1;
                oldRankUsed[bestOldRank] = 1;
            }
        }

        for (label newPart = 0; newPart < nProcs_; ++newPart)
        {
            if (remap[newPart] >= 0) continue;

            for (label oldRank = 0; oldRank < nProcs_; ++oldRank)
            {
                if (!oldRankUsed[oldRank])
                {
                    remap[newPart] = oldRank;
                    oldRankUsed[oldRank] = 1;
                    break;
                }
            }
        }

        label changedAfterRemap = 0;
        label remapBadIndex = -1;
        #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static) \
            reduction(+:changedAfterRemap)
        for (label i = 0; i < nCells; ++i)
        {
            const label newPart = label(fullPart[i]);
            if
            (
                newPart >= 0 && newPart < nProcs_
             && remap[newPart] >= 0
            )
            {
                const label newOwner = remap[newPart];
                fullPart[i] = idx_t(newOwner);
                if (newOwner != cellOwner_[i])
                {
                    ++changedAfterRemap;
                }
            }
            else
            {
                #pragma omp critical(dsmcReplicatedMeshRemapBadIndex)
                if (remapBadIndex < 0)
                {
                    remapBadIndex = i;
                }
            }
        }

        if (remapBadIndex >= 0)
        {
            const label newPart = label(fullPart[remapBadIndex]);
            FatalErrorInFunction
                << "Invalid replicated mesh DLB remap: newPart="
                << newPart << " nProcs=" << nProcs_
                << abort(FatalError);
        }

        remapSavedCells = changedBeforeRemap - changedAfterRemap;

        Info<< "Phase C ParMETIS remap: mode=" << remapMode
            << " changedBefore=" << changedBeforeRemap
            << " changedAfter=" << changedAfterRemap
            << " saved=" << remapSavedCells
            << " remap=(";
        forAll(remap, i)
        {
            Info<< (i ? " " : "") << i << "->" << remap[i];
        }
        Info<< ")" << endl;
    }
    else if (remapEnabled)
    {
        WarningInFunction
            << "Unsupported replicatedMeshDLBRemapMode '" << remapMode
            << "'. Falling back to ParMETIS part labels." << endl;
    }

    // ---- Update cellOwner_ ------------------------------------------------
    label nChanged = 0;
    #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static) \
        reduction(+:nChanged)
    for (label i = 0; i < nCells; ++i)
    {
        const label newOwner = label(fullPart[i]);
        if (newOwner != cellOwner_[i])
        {
            cellOwner_[i] = newOwner;
            ++nChanged;
        }
    }

    // ---- Rebuild local structures only when ownership changed --------------
    // ParMETIS can return the current partition (especially after a remap).
    // Avoid invalidating local state and rebuilding the local mesh in that case.
    if (nChanged > 0)
    {
        rebuildMyCells();
        localMesh_.build(myCells_);
    }

    // The ParMETIS output is globally reconstructed above. Validate the
    // resulting replicated map before particles are migrated to it.
    validateCellOwnerMap("ParMETIS AdaptiveRepart");

    totalCellsChanged_ += nChanged;

    Info<< "Phase C ParMETIS AdaptiveRepart: " << nChanged
        << " / " << nCells << " cells changed ("
        << scalar(nChanged) / scalar(nCells) << ")";
    if (remapEnabled)
    {
        Info<< " remapBefore=" << changedBeforeRemap
            << " remapSaved=" << remapSavedCells;
    }
    Info<< nl;

    return nChanged;
}


// ============================================================================
// autoRebalance — Phase C: automatic DLB trigger + execution
// ============================================================================

void dsmcReplicatedMesh::autoRebalance()
{
    if (!autoDLBEnabled_ || !active_ || nProcs_ < 2) return;

    const auto tAuto0 = std::chrono::steady_clock::now();
    ++autoRebalanceChecks_;

    const label currentStep = stepCounter_;
    const bool dlbProfile = mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshDLBProfile", false);
    const label minGapSteps = max
    (
        label(0),
        mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBMinGapSteps", 0)
    );
    const label stepsSinceRebalance = currentStep - lastAutoRebalanceStep_;
    const bool minGapSatisfied =
        (minGapSteps == 0)
     || (autoRebalanceCount_ == 0)
     || (lastAutoRebalanceStep_ < 0)
     || (stepsSinceRebalance >= minGapSteps);
    const word checkCollective = mesh_.time().controlDict().lookupOrDefault<word>
        ("replicatedMeshDLBCheckCollective", "allgather");
    const bool useAllreduceCheck = (checkCollective == "allreduce");

    auto gatherLoadExtrema =
        [&](const scalar localT, scalar& maxT, scalar& minT)
        {
            const auto tComm0 = std::chrono::steady_clock::now();
            if (useAllreduceCheck)
            {
                MPI_Allreduce
                (
                    &localT,
                    &maxT,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD
                );
                MPI_Allreduce
                (
                    &localT,
                    &minT,
                    1,
                    MPI_DOUBLE,
                    MPI_MIN,
                    MPI_COMM_WORLD
                );
            }
            else
            {
                scalarList allT(nProcs_, 0.0);
                MPI_Allgather
                (
                    &localT,
                    1,
                    MPI_DOUBLE,
                    allT.data(),
                    1,
                    MPI_DOUBLE,
                    MPI_COMM_WORLD
                );
                maxT = max(allT);
                minT = min(allT);
            }
            addAutoRebalanceCommTime
            (
                std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now() - tComm0
                ).count()
            );
        };

    // Post-DLB snapshot: per-rank load for 5 steps after DLB
    if (dlbProfile && postDLBSnapshotCountdown_ > 0)
    {
        const scalar myMoveT = cloud_.evolveMoveWallTime();
        const scalar myCollT = cloud_.evolveCollisionWallTime();
        scalarList snapMove(nProcs_, 0.0), snapColl(nProcs_, 0.0);
        const auto tPostSnapshotComm0 = std::chrono::steady_clock::now();
        MPI_Allgather(&myMoveT, 1, MPI_DOUBLE,
                      snapMove.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(&myCollT, 1, MPI_DOUBLE,
                      snapColl.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
        addAutoRebalanceCommTime
        (
            std::chrono::duration<scalar>
            (
                std::chrono::steady_clock::now() - tPostSnapshotComm0
            ).count()
        );
        if (lastDLBMoveT_.size() == nProcs_)
        {
            const label nSnap = mesh_.time().controlDict().lookupOrDefault<label>
                ("replicatedMeshDLBProfileSteps", 5);
            const label step = nSnap + 1 - postDLBSnapshotCountdown_;
            scalar maxLoad = 0.0, minLoad = GREAT;
            Info<< "  [post-DLB step " << step << "] per-rank load:";
            for (label i = 0; i < nProcs_; ++i)
            {
                const scalar load = (snapMove[i] - lastDLBMoveT_[i])
                    + (snapColl[i] - lastDLBCollT_[i]);
                Info<< " r" << i << "=" << load;
                if (load > maxLoad) maxLoad = load;
                if (load < minLoad) minLoad = load;
            }
            Info<< " | max/min=" << maxLoad / max(minLoad, SMALL) << nl;
        }
        --postDLBSnapshotCountdown_;
    }

    // ---- Accumulate productive time ----
    // Legacy metric: CPU time of move+collision minus migration wall.
    // Wall metric (§5.6): per-iteration full wall reported by the solver
    // main loop minus that iteration's collective waits.  Both inputs lag
    // one step, so the pair refers to the same iteration: lastStepWallTime_
    // is set at the end of iteration N-1, and lastStepWait_ holds the wait
    // delta captured at the previous autoRebalance (waits of iteration N-1).
    scalar stepProductive = 0.0;
    if (triggerWallTime_)
    {
        stepProductive = max(lastStepWallTime_ - lastStepWait_, scalar(0));
        stepFullWallCum_ += lastStepWallTime_;
    }
    else
    {
        const scalar stepEvolve = evolveStepTime_ - lastEvolveTime_;
        const scalar stepMig = migrationWallTime_ - lastMigrationTime_;
        stepProductive = stepEvolve - stepMig;
    }
    // Collective-wait delta of the current step, consumed one step later.
    const scalar sizeExchangeNow = migrationSizeExchangeWallTime_;
    const scalar waitNow = migrationWaitWallTime_;
    const scalar upcNow = updateParticleCountsWallTime_;
    lastStepWait_ =
        (sizeExchangeNow - lastSizeExchangeWall_)
      + (waitNow - lastWaitWall_)
      + (upcNow - lastUpcWall_);
    lastSizeExchangeWall_ = sizeExchangeNow;
    lastWaitWall_ = waitNow;
    lastUpcWall_ = upcNow;
    productiveTime_ += stepProductive;
    lastEvolveTime_ = evolveStepTime_;
    lastMigrationTime_ = migrationWallTime_;
    ++ndecps_;
    ++sarEvalSteps_;

    // Minimum step guard (ParMETIS needs stable particle distribution)
    if (currentStep < 30) return;

    bool triggered = false;
    bool forcedTriggered = false;
    bool globalTriggerDecisionComputed = false;
    // Global imbalance measured at this check (>0 when a check ran); used
    // by the back-off gate and recorded for the executed rebalance.
    scalar checkImbalance = -1.0;

    if (forcedDLBSteps_.size())
    {
        forAll(forcedDLBSteps_, i)
        {
            if (forcedDLBSteps_[i] == currentStep)
            {
                forcedTriggered = true;
                break;
            }
        }
    }

    const bool particleGateEnabled = mesh_.time().controlDict()
        .lookupOrDefault<bool>("replicatedMeshDLBParticleGate", true);
    const scalar particleGateThreshold = mesh_.time().controlDict()
        .lookupOrDefault<scalar>
        (
            "replicatedMeshDLBParticleGateThreshold",
            imbalanceThreshold_
        );
    bool particleGateAllowsGlobalCheck = true;
    scalar particleGateImbalance = 1.0;
    if
    (
        particleGateEnabled
     && !forcedTriggered
     && allParticleCounts_.size() == nProcs_
    )
    {
        label maxParticles = 0;
        label minParticles = labelMax;
        forAll(allParticleCounts_, i)
        {
            maxParticles = max(maxParticles, allParticleCounts_[i]);
            minParticles = min(minParticles, allParticleCounts_[i]);
        }

        if (maxParticles > 0 && minParticles < labelMax)
        {
            particleGateImbalance =
                scalar(maxParticles) / max(scalar(minParticles), SMALL);
        }
        particleGateAllowsGlobalCheck =
            particleGateImbalance > particleGateThreshold;
    }

    if (forcedDLBSteps_.size())
    {
        if (sarEvalSteps_ >= sarSteps_)
        {
            const auto tForcedBarrier0 = std::chrono::steady_clock::now();
            MPI_Barrier(MPI_COMM_WORLD);
            addAutoRebalanceCommTime
            (
                std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now() - tForcedBarrier0
                ).count()
            );
            sarEvalSteps_ = 0;
        }

        if (forcedTriggered)
        {
            triggered = true;
            Info<< "\nPhase C auto DLB triggered (forced)"
                << " at step " << currentStep << nl;
        }
    }
    else if (autoDLBTriggerMode_ == "legacyWindow")
    {
        if (ndecps_ % sarSteps_ == 0)
        {
            if (!particleGateAllowsGlobalCheck)
            {
                if (dlbProfile && myRank_ == 0)
                {
                    Info<< "Phase C auto DLB check skipped by particle gate"
                        << " at step " << currentStep
                        << ": particle max/min=" << particleGateImbalance
                        << " threshold=" << particleGateThreshold << nl;
                }
            }
            else
            {
                const scalar localT = productiveTime_;
                scalar maxT = 0.0;
                scalar minT = GREAT;
                gatherLoadExtrema(localT, maxT, minT);

                if (maxT > SMALL)
                {
                    const scalar loadImbalance = maxT / max(minT, SMALL);
                    checkImbalance = loadImbalance;
                    printPerRankStepWall();

                    bool sarTriggered = false;

                    tidl_ += maxT - minT;
                    ++nEvalPeriods_;
                    w2_ = (tidl_ + tdecps_) / scalar(nEvalPeriods_);

                    if (nEvalPeriods_ == 1)
                    {
                        w1_ = w2_;
                    }
                    else
                    {
                        sar_ = w2_ - w1_;
                        w1_ = w2_;
                        sarTriggered = (sar_ > 0.0);
                    }

                    const bool thresholdTriggered =
                        (ndecps_ >= dlbSteps_)
                     && (loadImbalance > imbalanceThreshold_);

                    triggered =
                        minGapSatisfied && (sarTriggered || thresholdTriggered);
                    globalTriggerDecisionComputed = true;
                    if (triggered)
                    {
                        const char* reason =
                            sarTriggered && thresholdTriggered
                          ? "SAR+threshold"
                          : (sarTriggered ? "SAR" : "threshold");

                        Info<< "\nPhase C auto DLB triggered (" << reason
                            << ") at step " << currentStep
                            << ": sar=" << sar_
                            << " loadImbalance=" << loadImbalance
                            << " threshold=" << imbalanceThreshold_
                            << " particleMaxMin=" << particleGateImbalance
                            << nl;
                    }
                }
            }

            productiveTime_ = 0.0;
            if (ndecps_ >= dlbSteps_) ndecps_ = 0;
        }
    }
    else if (sarEvalSteps_ >= sarSteps_)
    {
        const bool thresholdEligible =
            (autoRebalanceCount_ > 0) && (ndecps_ >= dlbSteps_);

        if (!particleGateAllowsGlobalCheck)
        {
            if (dlbProfile && myRank_ == 0)
            {
                Info<< "Phase C auto DLB check skipped by particle gate"
                    << " at step " << currentStep
                    << ": particle max/min=" << particleGateImbalance
                    << " threshold=" << particleGateThreshold << nl;
            }
        }
        else
        {
            const scalar localT = productiveTime_;
            scalar maxT = 0.0;
            scalar minT = GREAT;
            gatherLoadExtrema(localT, maxT, minT);

            if (maxT > SMALL)
            {
                const scalar loadImbalance = maxT / max(minT, SMALL);
                checkImbalance = loadImbalance;
                printPerRankStepWall();

                // SAR trend update
                tidl_ += maxT - minT;
                ++nEvalPeriods_;
                w2_ = (tidl_ + tdecps_) / scalar(nEvalPeriods_);

                bool sarTriggered = false;
                if (nEvalPeriods_ == 1)
                {
                    w1_ = w2_;
                }
                else
                {
                    sar_ = w2_ - w1_;
                    w1_ = w2_;
                    sarTriggered = (sar_ > 0.0);
                }

                const bool thresholdTriggered =
                    thresholdEligible && (loadImbalance > imbalanceThreshold_);

                triggered =
                    minGapSatisfied && (sarTriggered || thresholdTriggered);
                globalTriggerDecisionComputed = true;
                if (triggered)
                {
                    const char* reason =
                        sarTriggered && thresholdTriggered
                      ? "SAR+threshold"
                      : (sarTriggered ? "SAR" : "threshold");

                    Info<< "\nPhase C auto DLB triggered (" << reason
                        << ") at step " << currentStep
                        << ": sar=" << sar_
                        << " loadImbalance=" << loadImbalance
                        << " threshold=" << imbalanceThreshold_
                        << " particleMaxMin=" << particleGateImbalance
                        << nl;
                }
            }
        }

        // Current SAR sampling window has been consumed, but productiveTime_
        // remains cumulative until an actual DLB trigger completes.
        sarEvalSteps_ = 0;
    }

    if (triggered)
    {
        const label minRemainingSteps = max
        (
            label(0),
            mesh_.time().controlDict().lookupOrDefault<label>
            (
                "replicatedMeshDLBMinRemainingSteps",
                0
            )
        );

        if (minRemainingSteps > 0)
        {
            const scalar startTime = mesh_.time().controlDict()
                .lookupOrDefault<scalar>("startTime", 0.0);
            const scalar endTime = mesh_.time().controlDict()
                .lookupOrDefault<scalar>("endTime", mesh_.time().endTime().value());
            const scalar deltaT = max(mesh_.time().deltaT().value(), SMALL);
            const label totalSteps = max
            (
                label(0),
                label((endTime - startTime)/deltaT + 0.5)
            );
            const label remainingSteps = totalSteps - currentStep;

            if (remainingSteps < minRemainingSteps)
            {
                if (myRank_ == 0)
                {
                    Info<< "Phase C auto DLB skipped near end at step "
                        << currentStep
                        << ": remainingSteps=" << remainingSteps
                        << " minRemainingSteps=" << minRemainingSteps
                        << nl;
                }
                triggered = false;
            }
        }
    }

    // ---- Futile-rebalance back-off (§5.3) ----
    // Applies to auto triggers only (checkImbalance > 0; forced triggers
    // carry checkImbalance = -1).  All inputs are globally identical on
    // every rank, so the bookkeeping stays collective-consistent.
    if (triggered && backOffEnabled_ && checkImbalance > 0)
    {
        const bool futile =
            lastExecChangedRatio_ < backOffChangedRatio_
         && checkImbalance >= lastExecImbalance_ - backOffTol_;
        futileCount_ = futile ? (futileCount_ + 1) : 0;
        if (futileCount_ >= backOffFutileRuns_ && backOffRemaining_ <= 0)
        {
            backOffRemaining_ = backOffSkipChecks_;
        }

        const scalar escalateImbalance =
            imbalanceThreshold_ * backOffEscalate_;
        if (backOffRemaining_ > 0 && checkImbalance < escalateImbalance)
        {
            if (myRank_ == 0)
            {
                Info<< "Phase C auto DLB skipped by back-off at step "
                    << currentStep
                    << ": remaining=" << backOffRemaining_
                    << " imbalance=" << checkImbalance
                    << " lastChangedRatio=" << lastExecChangedRatio_
                    << nl;
            }
            --backOffRemaining_;
            triggered = false;
        }
        else if (backOffRemaining_ > 0)
        {
            // Safety valve: imbalance well above threshold — execute now.
            if (myRank_ == 0)
            {
                Info<< "Phase C auto DLB back-off overridden at step "
                    << currentStep
                    << ": imbalance=" << checkImbalance
                    << " >= escalate=" << escalateImbalance << nl;
            }
            backOffRemaining_ = 0;
        }
    }

    if (globalTriggerDecisionComputed)
    {
        const int localTriggered = triggered ? 1 : 0;
        int triggeredSum = 0;
        const auto tTriggerComm0 = std::chrono::steady_clock::now();
        MPI_Allreduce
        (
            &localTriggered,
            &triggeredSum,
            1,
            MPI_INT,
            MPI_SUM,
            MPI_COMM_WORLD
        );
        addAutoRebalanceCommTime
        (
            std::chrono::duration<scalar>
            (
                std::chrono::steady_clock::now() - tTriggerComm0
            ).count()
        );

        if
        (
            triggeredSum != 0
         && triggeredSum != nProcs_
         && myRank_ == 0
        )
        {
            Info<< "Phase C auto DLB trigger mismatch across ranks at step "
                << currentStep
                << ": local decisions were inconsistent, promoting to global OR"
                << nl;
        }

        triggered = (triggeredSum != 0);
    }

    const auto tCheckEnd = std::chrono::steady_clock::now();
    autoRebalanceCheckWallTime_ +=
        std::chrono::duration<scalar>(tCheckEnd - tAuto0).count();

    if (!triggered)
    {
        autoRebalanceWallTime_ +=
            std::chrono::duration<scalar>(tCheckEnd - tAuto0).count();
        return;
    }

    ++autoRebalanceTriggeredChecks_;

    // Time the decomposition for cost accounting
    const auto tDecStart = std::chrono::steady_clock::now();

    // Reassign by ParMETIS AdaptiveRepart
    const label nChanged = reassignByParMetisAdaptiveRepart();

    // Record the executed rebalance for back-off bookkeeping (§5.3) and
    // reset the back-off state: a fresh partition needs re-evaluation.
    lastExecChangedRatio_ =
        scalar(nChanged)/max(scalar(mesh_.nCells()), scalar(1));
    lastExecImbalance_ = max(checkImbalance, scalar(1));
    backOffRemaining_ = 0;
    // §11.8: reset the futile counter only when the executed rebalance was
    // actually effective.  The previous unconditional reset (worklog
    // §11.8) pinned futileCount_ at 1 in a futile steady state (remap
    // rejecting the ParMETIS output), so the back-off could never reach
    // futileRuns and every check re-entered ParMETIS (500wcell-1bparticle
    // job 4776404: 159 executions, each changing ~50 of 5.82M cells).
    if (lastExecChangedRatio_ >= backOffChangedRatio_)
    {
        futileCount_ = 0;
    }
    const auto tRepartEnd = std::chrono::steady_clock::now();
    autoRebalanceRepartWallTime_ +=
        std::chrono::duration<scalar>(tRepartEnd - tDecStart).count();

    // Skip migration if partition unchanged
    if (nChanged > 0)
    {
        const auto tMigrate0 = std::chrono::steady_clock::now();
        const label oldMigrationContext = migrationProfileContext_;
        setMigrationProfileContext(migrationProfileDLB);
        migrateParticlesByCellOwner();
        updateParticleCounts();
        setMigrationProfileContext(oldMigrationContext);
        const auto tMigrate1 = std::chrono::steady_clock::now();
        autoRebalanceMigrationWallTime_ +=
            std::chrono::duration<scalar>(tMigrate1 - tMigrate0).count();
    }

    const auto tDecEnd = std::chrono::steady_clock::now();
    const scalar localDecWall =
        std::chrono::duration<scalar>(tDecEnd - tDecStart).count();
    scalar globalDecWall = localDecWall;
    const auto tDecWallComm0 = std::chrono::steady_clock::now();
    MPI_Allreduce
    (
        &localDecWall,
        &globalDecWall,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    addAutoRebalanceCommTime
    (
        std::chrono::duration<scalar>
        (
            std::chrono::steady_clock::now() - tDecWallComm0
        ).count()
    );
    tdecps_ += globalDecWall;

    if (nChanged > 0)
    {
        lastAutoRebalanceStep_ = currentStep;
        ++autoRebalanceCount_;
    }

    // Reset SAR state: w1_=w2_ keeps current idle rate as baseline
    tidl_ = 0.0;
    nEvalPeriods_ = 0;
    sar_ = 0.0;
    ndecps_ = 0;
    sarEvalSteps_ = 0;
    productiveTime_ = 0.0;
    if (autoDLBTriggerMode_ != "legacyWindow")
    {
        lastEvolveTime_ = evolveStepTime_;
        lastMigrationTime_ = migrationWallTime_;
    }

    // Reset per-cell move iteration counter after DLB uses it
    cloud_.moveItersPerCell() = 0;

    if (nChanged > 0)
    {
        Info<< "Phase C auto DLB complete: rebalance #"
            << autoRebalanceCount_ << nl << endl;
    }
    else
    {
        Info<< "Phase C auto DLB complete: ParMETIS kept the current cell "
            << "owners; actual rebalance count remains "
            << autoRebalanceCount_ << nl << endl;
    }

    // Inter-DLB load summary + arm post-DLB snapshot
    const auto tPostDiag0 = std::chrono::steady_clock::now();
    {
        const bool skipPostDiag = mesh_.time().controlDict()
            .lookupOrDefault<bool>("replicatedMeshDLBSkipPostDiag", false);
        const bool needPostDiag = dlbProfile || adaptiveAlpha_ || !skipPostDiag;
        scalarList allMoveT, allCollT;
        if (needPostDiag)
        {
            const scalar myMoveT = cloud_.evolveMoveWallTime();
            const scalar myCollT = cloud_.evolveCollisionWallTime();
            allMoveT.setSize(nProcs_, 0.0);
            allCollT.setSize(nProcs_, 0.0);
            const auto tPostDiagComm0 = std::chrono::steady_clock::now();
            MPI_Allgather(&myMoveT, 1, MPI_DOUBLE,
                          allMoveT.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgather(&myCollT, 1, MPI_DOUBLE,
                          allCollT.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
            addAutoRebalanceCommTime
            (
                std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now() - tPostDiagComm0
                ).count()
            );
        }
        if (dlbProfile && lastDLBMoveT_.size() == nProcs_)
        {
            Info<< "  Inter-DLB actual load (since last rebalance):" << nl;
            scalar maxA = 0.0, minA = GREAT;
            for (label i = 0; i < nProcs_; ++i)
            {
                const scalar dM = allMoveT[i] - lastDLBMoveT_[i];
                const scalar dC = allCollT[i] - lastDLBCollT_[i];
                Info<< "    rank" << i << ": move=" << dM
                    << " coll=" << dC << " total=" << (dM+dC) << nl;
                if (dM+dC > maxA) maxA = dM+dC;
                if (dM+dC < minA) minA = dM+dC;
            }
            Info<< "    max/min=" << maxA / max(minA, SMALL)
                << " (max=" << maxA << " min=" << minA << ")" << nl;
        }

        if (adaptiveAlpha_ && lastDLBMoveT_.size() == nProcs_)
        {
            scalar maxWork = 0.0;
            scalar minWork = GREAT;
            label bottleneckRank = 0;
            scalar sumMove = 0.0;
            scalar sumColl = 0.0;
            scalarList dMove(nProcs_, 0.0);
            scalarList dColl(nProcs_, 0.0);

            for (label i = 0; i < nProcs_; ++i)
            {
                dMove[i] = max
                (
                    scalar(0.0),
                    allMoveT[i] - lastDLBMoveT_[i]
                );
                dColl[i] = max
                (
                    scalar(0.0),
                    allCollT[i] - lastDLBCollT_[i]
                );
                const scalar work = dMove[i] + dColl[i];
                sumMove += dMove[i];
                sumColl += dColl[i];
                if (work > maxWork)
                {
                    maxWork = work;
                    bottleneckRank = i;
                }
                if (work < minWork)
                {
                    minWork = work;
                }
            }

            const scalar totalWork = sumMove + sumColl;
            const scalar workImbalance = maxWork / max(minWork, SMALL);
            const scalar oldAlpha = dlbAlpha_;
            scalar direction = 0.0;
            scalar step = 0.0;

            if (totalWork > SMALL && workImbalance > scalar(1.0) + SMALL)
            {
                const scalar avgMove = sumMove / scalar(nProcs_);
                const scalar avgColl = sumColl / scalar(nProcs_);
                const scalar bnMoveExcess = dMove[bottleneckRank] - avgMove;
                const scalar bnCollExcess = dColl[bottleneckRank] - avgColl;
                const scalar excessNorm =
                    mag(bnMoveExcess) + mag(bnCollExcess) + SMALL;

                // Move-heavy bottleneck: compress particle weights by reducing
                // alpha. Collision-heavy bottleneck: increase alpha.
                direction = (bnCollExcess - bnMoveExcess) / excessNorm;

                if
                (
                    workImbalance > adaptiveAlphaLastImbalance_
                  + adaptiveAlphaWorsenTol_
                 && mag(adaptiveAlphaLastStep_) > SMALL
                )
                {
                    direction = adaptiveAlphaLastStep_ > 0 ? -1.0 : 1.0;
                }

                const scalar severity =
                    min(max(workImbalance - scalar(1.0), scalar(0.0)), scalar(2.0));
                step = direction * min
                (
                    adaptiveAlphaMaxStep_,
                    adaptiveAlphaGain_ * severity
                );

                dlbAlpha_ += step;
                dlbAlpha_ = max(adaptiveAlphaMin_, min(adaptiveAlphaMax_, dlbAlpha_));
            }

            adaptiveAlphaLastStep_ = dlbAlpha_ - oldAlpha;
            adaptiveAlphaLastImbalance_ = workImbalance;

            Info<< "Phase C adaptive alpha: bottleneck=rank" << bottleneckRank
                << " workImbalance=" << workImbalance
                << " moveDelta=" << dMove[bottleneckRank]
                << " collDelta=" << dColl[bottleneckRank]
                << " direction=" << direction
                << " step=" << adaptiveAlphaLastStep_
                << " alpha: " << oldAlpha << " -> " << dlbAlpha_ << nl;
        }

        if (needPostDiag)
        {
            lastDLBMoveT_ = allMoveT;
            lastDLBCollT_ = allCollT;
            postDLBSnapshotCountdown_ = mesh_.time().controlDict().lookupOrDefault<label>
                ("replicatedMeshDLBProfileSteps", 5);
        }
        else
        {
            lastDLBMoveT_.clear();
            lastDLBCollT_.clear();
            postDLBSnapshotCountdown_ = 0;
        }
    }
    const auto tAutoEnd = std::chrono::steady_clock::now();
    autoRebalancePostDiagWallTime_ +=
        std::chrono::duration<scalar>(tAutoEnd - tPostDiag0).count();
    autoRebalanceWallTime_ +=
        std::chrono::duration<scalar>(tAutoEnd - tAuto0).count();
}


void dsmcReplicatedMesh::writeCellOwner() const
{
    if (!active_ || myRank_ != 0) return;

    volScalarField cellOwnerField
    (
        IOobject
        (
            "cellOwner",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );

    forAll(cellOwner_, cellI)
    {
        cellOwnerField[cellI] = scalar(cellOwner_[cellI]);
    }

    cellOwnerField.write();
    Info<< "Written cellOwner field at output time "
        << mesh_.time().timeName() << endl;
}


void dsmcReplicatedMesh::writeGatheredCloudOnRank0() const
{
    if (!active_ || myRank_ != 0) return;

    IOdictionary uniformPropsDict
    (
        IOobject
        (
            Cloud<dsmcParcel>::cloudPropertiesName,
            mesh_.time().timeName(),
            "uniform"/cloud::prefix/cloud_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    for (label procI = 0; procI < nProcs_; ++procI)
    {
        const word procName("processor" + Foam::name(procI));
        uniformPropsDict.add(procName, dictionary());
        uniformPropsDict.subDict(procName).add
        (
            "particleCount",
            dsmcParcel::particleCount_
        );
    }

    uniformPropsDict.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        mesh_.time().writeCompression()
    );

    if (!cloud_.size())
    {
        Info<< "Replicated mesh: rank0 gathered cloud is empty at output time "
            << mesh_.time().timeName() << endl;
        return;
    }

    particle::writeFields(cloud_);

    const label np = cloud_.size();

    IOField<vector> U(cloud_.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> RWF
    (
        cloud_.fieldIOobject("radialWeight", IOobject::NO_READ),
        np
    );
    IOField<scalar> ERot
    (
        cloud_.fieldIOobject("ERot", IOobject::NO_READ),
        np
    );
    IOField<labelField> vibLevel
    (
        cloud_.fieldIOobject("vibLevel", IOobject::NO_READ),
        np
    );
    IOField<label> ELevel
    (
        cloud_.fieldIOobject("ELevel", IOobject::NO_READ),
        np
    );
    IOField<label> typeId
    (
        cloud_.fieldIOobject("typeId", IOobject::NO_READ),
        np
    );
    IOField<label> newParcel
    (
        cloud_.fieldIOobject("newParcel", IOobject::NO_READ),
        np
    );
    IOField<label> classification
    (
        cloud_.fieldIOobject("classification", IOobject::NO_READ),
        np
    );

    IOField<label> stuckToWall
    (
        cloud_.fieldIOobject("stuckToWall", IOobject::NO_READ),
        np
    );
    IOField<scalarField> wallTemperature
    (
        cloud_.fieldIOobject("wallTemperature", IOobject::NO_READ),
        np
    );
    IOField<vectorField> wallVectors
    (
        cloud_.fieldIOobject("wallVectors", IOobject::NO_READ),
        np
    );

    IOField<label> isTracked
    (
        cloud_.fieldIOobject("isTracked", IOobject::NO_READ),
        np
    );
    IOField<label> inPatchId
    (
        cloud_.fieldIOobject("inPatchId", IOobject::NO_READ),
        np
    );
    IOField<scalar> tracerInitialTime
    (
        cloud_.fieldIOobject("tracerInitialTime", IOobject::NO_READ),
        np
    );
    IOField<vector> tracerInitialPosition
    (
        cloud_.fieldIOobject("tracerInitialPosition", IOobject::NO_READ),
        np
    );
    IOField<vector> tracerCurrentPosition
    (
        cloud_.fieldIOobject("tracerCurrentPosition", IOobject::NO_READ),
        np
    );
    IOField<vector> tracerDistanceTravelled
    (
        cloud_.fieldIOobject("tracerDistanceTravelled", IOobject::NO_READ),
        np
    );

    bool hasRWF = false;
    bool hasERot = false;
    bool hasELevel = false;
    bool hasStuck = false;
    bool hasTracked = false;

    label i = 0;
    forAllConstIter(Cloud<dsmcParcel>, cloud_, iter)
    {
        const dsmcParcel& p = iter();

        U[i] = p.U();
        RWF[i] = p.RWF();
        ERot[i] = p.ERot();
        vibLevel[i] = p.vibLevel();
        ELevel[i] = p.ELevel();
        typeId[i] = p.typeId();
        newParcel[i] = p.newParcel();
        classification[i] = p.classification();

        stuckToWall[i] = p.isStuck();
        if (stuckToWall[i])
        {
            wallTemperature[i] = p.stuck().wallTemperature();
            wallVectors[i] = p.stuck().wallVectors();
            hasStuck = true;
        }

        isTracked[i] = p.isTracked();
        if (isTracked[i])
        {
            inPatchId[i] = p.tracked().inPatchId();
            tracerInitialTime[i] = p.tracked().initialTime();
            tracerInitialPosition[i] = p.tracked().initialPosition();
            tracerCurrentPosition[i] = p.tracked().currentPosition();
            tracerDistanceTravelled[i] = p.tracked().distanceTravelledVector();
            hasTracked = true;
        }
        else
        {
            inPatchId[i] = -1;
            tracerInitialTime[i] = 0;
            tracerInitialPosition[i] = vector::zero;
            tracerCurrentPosition[i] = vector::zero;
            tracerDistanceTravelled[i] = vector::zero;
        }

        hasRWF = hasRWF || RWF[i] > 1.0;
        hasERot = hasERot || ERot[i] > 0.0;
        hasELevel = hasELevel || ELevel[i] > 0;

        ++i;
    }

    U.write();

    if (hasRWF)
    {
        RWF.write();
    }

    if (hasERot)
    {
        ERot.write();
    }

    if (hasELevel)
    {
        ELevel.write();
    }

    typeId.write();
    newParcel.write();
    classification.write();

    if (hasStuck)
    {
        stuckToWall.write();
        wallTemperature.write();
        wallVectors.write();
    }

    if (hasTracked)
    {
        isTracked.write();
        inPatchId.write();
        tracerInitialTime.write();
        tracerInitialPosition.write();
        tracerCurrentPosition.write();
        tracerDistanceTravelled.write();
    }

    vibLevel.write();

    Info<< "Replicated mesh: wrote gathered cloud on rank0 at output time "
        << mesh_.time().timeName() << " with " << np << " parcels" << endl;
}


bool dsmcReplicatedMesh::processorWriteEnabled() const
{
    return active_ && writeMode_ == "processor";
}


void dsmcReplicatedMesh::writeProcessorOutput() const
{
    if (!active_) return;

    const Time& runTime = mesh_.time();
    const word timeName = runTime.timeName();
    const bool writeTimeMesh = runTime.controlDict().lookupOrDefault<bool>
    (
        "replicatedMeshProcessorWriteTimeMesh",
        false
    );
    const fileName meshReadInstance = mesh_.facesInstance();
    const label ownerVersion = rebalanceCount_ + autoRebalanceCount_;
    word procMeshInstance =
        writeTimeMesh || ownerVersion > 0 ? timeName : runTime.constant();
    if
    (
        processorWriteMeshInstance_ != word::null
     && processorWriteDecompVersion_ == ownerVersion
    )
    {
        procMeshInstance = processorWriteMeshInstance_;
    }

    if (myRank_ == 0)
    {
        Info<< "Replicated mesh: processor output begin at time "
            << timeName << " (read mesh instance " << meshReadInstance
            << ", processor mesh instance " << procMeshInstance << ")"
            << endl;
    }

    Info<< "Replicated mesh: constructing processor mesh/addressing"
        << endl;

    const fileName processorDir =
        runTime.path()/fileName(word("processor") + Foam::name(myRank_));

    auto processorPolyMeshComplete =
        [&](const word& instance)
        {
            const fileName procPolyMeshDir =
                processorDir/instance/polyMesh::meshSubDir;

            return
                isFile(procPolyMeshDir/"cellProcAddressing")
             && isFile(procPolyMeshDir/"faceProcAddressing")
             && isFile(procPolyMeshDir/"boundaryProcAddressing")
             && isFile(procPolyMeshDir/"faces")
             && isFile(procPolyMeshDir/"owner")
             && isFile(procPolyMeshDir/"neighbour")
             && isFile(procPolyMeshDir/"boundary")
             && isFile(procPolyMeshDir/"points");
        };

    auto writeProcessorMesh =
        [&](const word& instance, const char* role)
        {
            domainDecomposition decomposition
            (
                IOobject
                (
                    mesh_.name(),
                    meshReadInstance,
                    runTime,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

            if (decomposition.nProcs() != nProcs_)
            {
                FatalErrorInFunction
                    << "decomposeParDict numberOfSubdomains "
                    << decomposition.nProcs()
                    << " does not match replicated mesh MPI ranks "
                    << nProcs_ << exit(FatalError);
            }

            decomposition.setCellToProc(cellOwner_);
            decomposition.setProcessorMeshInstance(instance);
            decomposition.setProcessorMeshWriteProc(myRank_);
            decomposition.setProcessorMeshAllProcPatches(true);
            decomposition.decomposeMesh();
            decomposition.writeDecomposition(false);

            Info<< "Replicated mesh: rank " << myRank_
                << " wrote processor" << myRank_ << " " << role
                << " mesh/addressing at instance " << instance << endl;
        };

    const word constantInstance = runTime.constant();
    if
    (
        procMeshInstance != constantInstance
     && !processorPolyMeshComplete(constantInstance)
    )
    {
        // reconstructPar constructs processorMeshes before stepping through
        // output times, so it requires a bootstrap processor constant mesh even
        // when the first valid owner partition is time-specific after DLB.
        writeProcessorMesh(constantInstance, "bootstrap");
    }

    if
    (
        processorWriteDecompVersion_ != ownerVersion
     || !processorPolyMeshComplete(procMeshInstance)
    )
    {
        writeProcessorMesh(procMeshInstance, "output");

        processorWriteDecompVersion_ = ownerVersion;
        processorWriteMeshInstance_ = procMeshInstance;

        Info<< "Replicated mesh: rank " << myRank_
            << " wrote processor" << myRank_ << " mesh/addressing"
            << endl;
    }
    else
    {
        Info<< "Replicated mesh: rank " << myRank_
            << " reusing processor" << myRank_ << " mesh/addressing"
            << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (myRank_ == 0)
    {
        Info<< "Replicated mesh: processor mesh/addressing barrier complete"
            << endl;
    }

    fileName processorCasePath
    (
        runTime.caseName()/fileName(word("processor") + Foam::name(myRank_))
    );

    Time processorDb
    (
        Time::controlDictName,
        runTime.rootPath(),
        processorCasePath,
        word("system"),
        word("constant")
    );
    processorDb.setTime(runTime);

    fvMesh procMesh
    (
        IOobject
        (
            mesh_.name(),
            procMeshInstance,
            processorDb,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    Info<< "Replicated mesh: rank " << myRank_
        << " loaded processor" << myRank_ << " mesh from instance "
        << procMeshInstance << endl;

    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellProcAddressing",
            procMesh.facesInstance(),
            procMesh.meshSubDir,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceProcAddressing",
            procMesh.facesInstance(),
            procMesh.meshSubDir,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    labelIOList boundaryProcAddressing
    (
        IOobject
        (
            "boundaryProcAddressing",
            procMesh.facesInstance(),
            procMesh.meshSubDir,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Replicated mesh: rank " << myRank_
        << " loaded processor addressing with "
        << cellProcAddressing.size() << " cells, "
        << faceProcAddressing.size() << " faces, "
        << boundaryProcAddressing.size() << " patches" << endl;

    fvFieldDecomposer fieldDecomposer
    (
        mesh_,
        procMesh,
        faceProcAddressing,
        cellProcAddressing,
        boundaryProcAddressing
    );

    const label nVolScalar =
        writeProcessorVolumeFields<volScalarField>(mesh_, fieldDecomposer);
    const label nVolVector =
        writeProcessorVolumeFields<volVectorField>(mesh_, fieldDecomposer);
    const label nVolSphericalTensor =
        writeProcessorVolumeFields<volSphericalTensorField>
        (
            mesh_,
            fieldDecomposer
        );
    const label nVolSymmTensor =
        writeProcessorVolumeFields<volSymmTensorField>(mesh_, fieldDecomposer);
    const label nVolTensor =
        writeProcessorVolumeFields<volTensorField>(mesh_, fieldDecomposer);

    Info<< "Replicated mesh: rank " << myRank_
        << " wrote processor" << myRank_ << " volume fields at output time "
        << timeName
        << " (scalar=" << nVolScalar
        << ", vector=" << nVolVector
        << ", sphericalTensor=" << nVolSphericalTensor
        << ", symmTensor=" << nVolSymmTensor
        << ", tensor=" << nVolTensor << ")" << endl;

    const bool writeCloud = runTime.controlDict().lookupOrDefault<bool>
    (
        "replicatedMeshProcessorWriteCloud",
        true
    );
    if (!writeCloud)
    {
        Info<< "Replicated mesh: rank " << myRank_
            << " skipped processor" << myRank_
            << " lagrangian cloud at output time " << timeName << endl;
        return;
    }

    labelList globalToLocal(mesh_.nCells(), -1);
    forAll(cellProcAddressing, procCellI)
    {
        const label globalCellI = cellProcAddressing[procCellI];
        if (globalCellI >= 0 && globalCellI < globalToLocal.size())
        {
            globalToLocal[globalCellI] = procCellI;
        }
    }

    label nParcels = 0;
    label nInvalid = 0;
    bool hasRWF = false;
    bool hasERot = false;
    bool hasELevel = false;
    bool hasStuck = false;
    bool hasTracked = false;

    forAllConstIter(Cloud<dsmcParcel>, cloud_, iter)
    {
        const dsmcParcel& p = iter();
        if
        (
            p.cell() >= 0
         && p.cell() < globalToLocal.size()
         && globalToLocal[p.cell()] >= 0
        )
        {
            ++nParcels;
            hasRWF = hasRWF || p.RWF() > 1.0;
            hasERot = hasERot || p.ERot() > 0.0;
            hasELevel = hasELevel || p.ELevel() > 0;
            hasStuck = hasStuck || p.isStuck();
            hasTracked = hasTracked || p.isTracked();
        }
        else
        {
            ++nInvalid;
        }
    }

    if (nInvalid)
    {
        FatalErrorInFunction
            << "Rank " << myRank_ << " has " << nInvalid
            << " parcels whose global cells are absent from processor"
            << myRank_ << " cellProcAddressing at time " << timeName
            << ". Aborting to avoid invalid reconstructPar output."
            << exit(FatalError);
    }

    int localFieldFlags[5] =
    {
        hasRWF ? 1 : 0,
        hasERot ? 1 : 0,
        hasELevel ? 1 : 0,
        hasStuck ? 1 : 0,
        hasTracked ? 1 : 0
    };
    int globalFieldFlags[5] = {0, 0, 0, 0, 0};
    MPI_Allreduce
    (
        localFieldFlags,
        globalFieldFlags,
        5,
        MPI_INT,
        MPI_MAX,
        MPI_COMM_WORLD
    );

    hasRWF = globalFieldFlags[0] != 0;
    hasERot = globalFieldFlags[1] != 0;
    hasELevel = globalFieldFlags[2] != 0;
    hasStuck = globalFieldFlags[3] != 0;
    hasTracked = globalFieldFlags[4] != 0;

    passiveParticleCloud positions
    (
        procMesh,
        cloud_.name(),
        IDLList<passiveParticle>()
    );

    IOField<vector> U(positions.fieldIOobject("U", IOobject::NO_READ), nParcels);
    IOField<labelField> vibLevel
    (
        positions.fieldIOobject("vibLevel", IOobject::NO_READ),
        nParcels
    );
    IOField<label> typeId
    (
        positions.fieldIOobject("typeId", IOobject::NO_READ),
        nParcels
    );
    IOField<label> newParcel
    (
        positions.fieldIOobject("newParcel", IOobject::NO_READ),
        nParcels
    );
    IOField<label> classification
    (
        positions.fieldIOobject("classification", IOobject::NO_READ),
        nParcels
    );
    IOField<label> origProcId
    (
        positions.fieldIOobject("origProcId", IOobject::NO_READ),
        nParcels
    );
    IOField<label> origId
    (
        positions.fieldIOobject("origId", IOobject::NO_READ),
        nParcels
    );

    IOField<scalar>* RWFPtr = hasRWF
      ? new IOField<scalar>
        (
            positions.fieldIOobject("radialWeight", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<scalar>* ERotPtr = hasERot
      ? new IOField<scalar>
        (
            positions.fieldIOobject("ERot", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<label>* ELevelPtr = hasELevel
      ? new IOField<label>
        (
            positions.fieldIOobject("ELevel", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<label>* stuckToWallPtr = hasStuck
      ? new IOField<label>
        (
            positions.fieldIOobject("stuckToWall", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<scalarField>* wallTemperaturePtr = hasStuck
      ? new IOField<scalarField>
        (
            positions.fieldIOobject("wallTemperature", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<vectorField>* wallVectorsPtr = hasStuck
      ? new IOField<vectorField>
        (
            positions.fieldIOobject("wallVectors", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<label>* isTrackedPtr = hasTracked
      ? new IOField<label>
        (
            positions.fieldIOobject("isTracked", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<label>* inPatchIdPtr = hasTracked
      ? new IOField<label>
        (
            positions.fieldIOobject("inPatchId", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<scalar>* tracerInitialTimePtr = hasTracked
      ? new IOField<scalar>
        (
            positions.fieldIOobject("tracerInitialTime", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<vector>* tracerInitialPositionPtr = hasTracked
      ? new IOField<vector>
        (
            positions.fieldIOobject("tracerInitialPosition", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<vector>* tracerCurrentPositionPtr = hasTracked
      ? new IOField<vector>
        (
            positions.fieldIOobject("tracerCurrentPosition", IOobject::NO_READ),
            nParcels
        )
      : nullptr;
    IOField<vector>* tracerDistanceTravelledPtr = hasTracked
      ? new IOField<vector>
        (
            positions.fieldIOobject("tracerDistanceTravelled", IOobject::NO_READ),
            nParcels
        )
      : nullptr;

    label i = 0;
    forAllConstIter(Cloud<dsmcParcel>, cloud_, iter)
    {
        const dsmcParcel& p = iter();
        const label localCellI = globalToLocal[p.cell()];

        positions.append
        (
            new passiveParticle
            (
                procMesh,
                p.position(),
                localCellI,
                true
            )
        );

        U[i] = p.U();
        vibLevel[i] = p.vibLevel();
        typeId[i] = p.typeId();
        newParcel[i] = p.newParcel();
        classification[i] = p.classification();
        origProcId[i] = p.origProc();
        origId[i] = p.origId();

        if (RWFPtr)
        {
            (*RWFPtr)[i] = p.RWF();
        }
        if (ERotPtr)
        {
            (*ERotPtr)[i] = p.ERot();
        }
        if (ELevelPtr)
        {
            (*ELevelPtr)[i] = p.ELevel();
        }

        if (stuckToWallPtr)
        {
            (*stuckToWallPtr)[i] = p.isStuck();
            if ((*stuckToWallPtr)[i])
            {
                (*wallTemperaturePtr)[i] = p.stuck().wallTemperature();
                (*wallVectorsPtr)[i] = p.stuck().wallVectors();
            }
            else
            {
                (*wallTemperaturePtr)[i] = scalarField(4, 0.0);
                (*wallVectorsPtr)[i] = vectorField(4, vector::zero);
            }
        }

        if (isTrackedPtr)
        {
            (*isTrackedPtr)[i] = p.isTracked();
            if ((*isTrackedPtr)[i])
            {
                (*inPatchIdPtr)[i] = p.tracked().inPatchId();
                (*tracerInitialTimePtr)[i] = p.tracked().initialTime();
                (*tracerInitialPositionPtr)[i] =
                    p.tracked().initialPosition();
                (*tracerCurrentPositionPtr)[i] =
                    p.tracked().currentPosition();
                (*tracerDistanceTravelledPtr)[i] =
                    p.tracked().distanceTravelledVector();
            }
            else
            {
                (*inPatchIdPtr)[i] = -1;
                (*tracerInitialTimePtr)[i] = 0;
                (*tracerInitialPositionPtr)[i] = vector::zero;
                (*tracerCurrentPositionPtr)[i] = vector::zero;
                (*tracerDistanceTravelledPtr)[i] = vector::zero;
            }
        }

        ++i;
    }

    IOdictionary uniformPropsDict
    (
        IOobject
        (
            Cloud<dsmcParcel>::cloudPropertiesName,
            processorDb.timeName(),
            "uniform"/cloud::prefix/cloud_.name(),
            procMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelList processorParticleCounts(nProcs_, 0);
    label myParticleCount = dsmcParcel::particleCount_;
    MPI_Allgather
    (
        &myParticleCount,
        1,
        MPI_INT,
        processorParticleCounts.data(),
        1,
        MPI_INT,
        MPI_COMM_WORLD
    );

    forAll(processorParticleCounts, procI)
    {
        const word procName("processor" + Foam::name(procI));
        uniformPropsDict.add(procName, dictionary());
        uniformPropsDict.subDict(procName).add
        (
            "particleCount",
            processorParticleCounts[procI]
        );
    }
    uniformPropsDict.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        processorDb.writeCompression()
    );

    Info<< "Replicated mesh: rank " << myRank_
        << " writing processor" << myRank_ << " lagrangian cloud"
        << endl;

    IOPosition<Cloud<passiveParticle>>(positions).write();
    origProcId.write();
    origId.write();
    U.write();
    if (RWFPtr)
    {
        RWFPtr->write();
    }
    if (ERotPtr)
    {
        ERotPtr->write();
    }
    if (ELevelPtr)
    {
        ELevelPtr->write();
    }
    typeId.write();
    newParcel.write();
    classification.write();
    if (stuckToWallPtr)
    {
        stuckToWallPtr->write();
        wallTemperaturePtr->write();
        wallVectorsPtr->write();
    }
    if (isTrackedPtr)
    {
        isTrackedPtr->write();
        inPatchIdPtr->write();
        tracerInitialTimePtr->write();
        tracerInitialPositionPtr->write();
        tracerCurrentPositionPtr->write();
        tracerDistanceTravelledPtr->write();
    }
    vibLevel.write();

    delete RWFPtr;
    delete ERotPtr;
    delete ELevelPtr;
    delete stuckToWallPtr;
    delete wallTemperaturePtr;
    delete wallVectorsPtr;
    delete isTrackedPtr;
    delete inPatchIdPtr;
    delete tracerInitialTimePtr;
    delete tracerInitialPositionPtr;
    delete tracerCurrentPositionPtr;
    delete tracerDistanceTravelledPtr;

    Info<< "Replicated mesh: rank " << myRank_
        << " wrote processor" << myRank_ << " cloud at output time "
        << timeName << " with " << nParcels << " parcels" << endl;

    MPI_Barrier(MPI_COMM_WORLD);
}


// ============================================================================
// distributeInitialParticles — delete non-owned particles on each rank
// Used when all ranks load the same initial particle data (replicated mesh
// with masterUncollated file handler). No MPI communication needed since
// each rank already has all particles.
// ============================================================================

void dsmcReplicatedMesh::distributeInitialParticles()
{
    if (!active_) return;

    DynamicList<dsmcParcel*> toDelete(cloud_.size() / 2);

    forAllIter(Cloud<dsmcParcel>, cloud_, iter)
    {
        dsmcParcel& p = iter();
        if (cellOwner_[p.cell()] != myRank_)
        {
            toDelete.append(&p);
        }
    }

    forAll(toDelete, i) { cloud_.deleteParcel(toDelete[i]); }

    Info<< "Replicated mesh: rank " << myRank_
        << " kept " << cloud_.size() << " parcels, deleted "
        << toDelete.size() << " non-owned parcels" << endl;

    ++migrationCalls_;
}


// ============================================================================
// migrateParticlesByCellOwner — Phase A
// ============================================================================

void dsmcReplicatedMesh::migrateParticlesByCellOwner()
{
    if (!active_) return;
    const auto tStart = std::chrono::steady_clock::now();

    if (useFlatTransfer_)
    {
        // ---- Flat POD transfer path ----
        List<DynamicList<dsmcParcel::TransferData>> sendTD(nProcs_);
        DynamicList<dsmcParcel*> toDelete(cloud_.size() / 4);
        DynamicList<dsmcParcel*> kept(cloud_.size());
        label nMigratedOut = 0;
        label nRecv = 0;
        const bool maintainMoveOrdered =
            cloud_.openmpEnabled()
         && cloud_.openmpMoveEnabled()
         && cloud_.ompNumThreads() > 1;
        const bool useOrderedTraversal =
            maintainMoveOrdered
         && cloud_.hasMoveOrderedParcels()
         && cloud_.moveOrderedParcels().size() == cloud_.size();

        if (useOrderedTraversal)
        {
#ifdef _OPENMP
            const label moveThreads = max(cloud_.ompNumThreads(), label(1));

            if (moveThreads > 1)
            {
                const auto& ordered = cloud_.moveOrderedParcels();
                const labelList& orderedOffsets =
                    cloud_.moveOrderedThreadOffsets();
                List<labelList> threadDestinationCounts(moveThreads);
                List<labelList> threadDestinationOffsets(moveThreads);
                labelList threadKeptCounts(moveThreads, 0);
                labelList threadDeleteCounts(moveThreads, 0);

                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    threadDestinationCounts[threadI].setSize(nProcs_, 0);
                    threadDestinationOffsets[threadI].setSize(nProcs_, 0);
                }

                #pragma omp parallel num_threads(moveThreads)
                {
                    const label threadI = omp_get_thread_num();
                    const label begin =
                        orderedOffsets.size() == moveThreads + 1
                      ? orderedOffsets[threadI]
                      : label
                        (
                            (
                                static_cast<long long>(threadI)
                              * static_cast<long long>(ordered.size())
                            )
                          / static_cast<long long>(moveThreads)
                        );
                    const label end =
                        orderedOffsets.size() == moveThreads + 1
                      ? orderedOffsets[threadI + 1]
                      : label
                        (
                            (
                                static_cast<long long>(threadI + 1)
                              * static_cast<long long>(ordered.size())
                            )
                          / static_cast<long long>(moveThreads)
                        );
                    labelList& destinationCounts =
                        threadDestinationCounts[threadI];

                    for (label i = begin; i < end; ++i)
                    {
                        const label dstRank =
                            cellOwner_[ordered[i]->cell()];
                        ++destinationCounts[dstRank];
                        if (dstRank == myRank_)
                        {
                            ++threadKeptCounts[threadI];
                        }
                        else
                        {
                            ++threadDeleteCounts[threadI];
                        }
                    }
                }

                labelList destinationCounts(nProcs_, 0);
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    for (label dst = 0; dst < nProcs_; ++dst)
                    {
                        destinationCounts[dst] +=
                            threadDestinationCounts[threadI][dst];
                    }
                    nMigratedOut += threadDeleteCounts[threadI];
                }

                label nKept = 0;
                labelList threadKeptOffsets(moveThreads + 1, 0);
                labelList threadDeleteOffsets(moveThreads + 1, 0);
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    threadKeptOffsets[threadI + 1] =
                        threadKeptOffsets[threadI] + threadKeptCounts[threadI];
                    threadDeleteOffsets[threadI + 1] =
                        threadDeleteOffsets[threadI]
                      + threadDeleteCounts[threadI];
                    nKept += threadKeptCounts[threadI];
                }

                kept.setSize(nKept);
                toDelete.setSize(nMigratedOut);
                for (label dst = 0; dst < nProcs_; ++dst)
                {
                    sendTD[dst].setSize(destinationCounts[dst]);
                }

                labelList destinationOffset(nProcs_, 0);
                for (label dst = 0; dst < nProcs_; ++dst)
                {
                    for (label threadI = 0; threadI < moveThreads; ++threadI)
                    {
                        threadDestinationOffsets[threadI][dst] =
                            destinationOffset[dst];
                        destinationOffset[dst] +=
                            threadDestinationCounts[threadI][dst];
                    }
                }

                #pragma omp parallel num_threads(moveThreads)
                {
                    const label threadI = omp_get_thread_num();
                    const label begin =
                        orderedOffsets.size() == moveThreads + 1
                      ? orderedOffsets[threadI]
                      : label
                        (
                            (
                                static_cast<long long>(threadI)
                              * static_cast<long long>(ordered.size())
                            )
                          / static_cast<long long>(moveThreads)
                        );
                    const label end =
                        orderedOffsets.size() == moveThreads + 1
                      ? orderedOffsets[threadI + 1]
                      : label
                        (
                            (
                                static_cast<long long>(threadI + 1)
                              * static_cast<long long>(ordered.size())
                            )
                          / static_cast<long long>(moveThreads)
                        );
                    labelList& destinationOffsets =
                        threadDestinationOffsets[threadI];
                    label keptI = threadKeptOffsets[threadI];
                    label deleteI = threadDeleteOffsets[threadI];

                    for (label i = begin; i < end; ++i)
                    {
                        dsmcParcel* pPtr = ordered[i];
                        const label dstRank = cellOwner_[pPtr->cell()];
                        if (dstRank != myRank_)
                        {
                            pPtr->packTransfer
                            (
                                sendTD[dstRank][destinationOffsets[dstRank]++]
                            );
                            toDelete[deleteI++] = pPtr;
                        }
                        else
                        {
                            kept[keptI++] = pPtr;
                        }
                    }
                }
            }
            else
#endif
            {
                const auto& ordered = cloud_.moveOrderedParcels();
                for (label i = 0; i < ordered.size(); ++i)
                {
                    dsmcParcel& p = *ordered[i];
                    const label dstRank = cellOwner_[p.cell()];
                    if (dstRank != myRank_)
                    {
                        sendTD[dstRank].append
                        (
                            dsmcParcel::TransferData()
                        );
                        p.packTransfer(sendTD[dstRank].last());
                        toDelete.append(ordered[i]);
                        ++nMigratedOut;
                    }
                    else
                    {
                        kept.append(ordered[i]);
                    }
                }
            }
        }
        else
        {
            forAllIter(Cloud<dsmcParcel>, cloud_, iter)
            {
                dsmcParcel& p = iter();
                const label dstRank = cellOwner_[p.cell()];
                if (dstRank != myRank_)
                {
                    sendTD[dstRank].append(dsmcParcel::TransferData());
                    p.packTransfer(sendTD[dstRank].last());
                    toDelete.append(&p);
                    ++nMigratedOut;
                }
                else
                {
                    kept.append(&p);
                }
            }
        }
        const auto tPackEnd = std::chrono::steady_clock::now();

        List<std::uint64_t> sendSizes
        (
            nProcs_,
            static_cast<std::uint64_t>(0)
        );
        for (label i = 0; i < nProcs_; ++i)
        {
            sendSizes[i] =
                static_cast<std::uint64_t>(sendTD[i].size())
              * static_cast<std::uint64_t>
                (sizeof(dsmcParcel::TransferData));
        }

        List<std::uint64_t> recvSizes
        (
            nProcs_,
            static_cast<std::uint64_t>(0)
        );
        MPI_Request sizeExchangeReq = MPI_REQUEST_NULL;
        if (!useNoAlltoall_ && overlapSizeExchange_)
        {
            MPI_Ialltoall
            (
                sendSizes.data(),
                1,
                MPI_UINT64_T,
                recvSizes.data(),
                1,
                MPI_UINT64_T,
                MPI_COMM_WORLD,
                &sizeExchangeReq
            );
        }

        forAll(toDelete, i) { cloud_.deleteParcel(toDelete[i]); }
        const auto tLocalPrepEnd = std::chrono::steady_clock::now();

        // Exchange 64-bit total byte counts.  The no-Alltoall mode retains
        // point-to-point headers; the regular mode uses the optimized
        // collective and can overlap it with local deletion.
        if (useNoAlltoall_)
        {
            DynamicList<MPI_Request> headerReqs(2*(nProcs_ - 1));
            for (label i = 0; i < nProcs_; ++i)
            {
                if (i != myRank_)
                {
                    MPI_Request req;
                    MPI_Irecv
                    (
                        &recvSizes[i],
                        1,
                        MPI_UINT64_T,
                        i,
                        100,
                        MPI_COMM_WORLD,
                        &req
                    );
                    headerReqs.append(req);
                }
            }
            for (label i = 0; i < nProcs_; ++i)
            {
                if (i != myRank_)
                {
                    MPI_Request req;
                    MPI_Isend
                    (
                        &sendSizes[i],
                        1,
                        MPI_UINT64_T,
                        i,
                        100,
                        MPI_COMM_WORLD,
                        &req
                    );
                    headerReqs.append(req);
                }
            }
            if (headerReqs.size() > 0)
            {
                MPI_Waitall
                (
                    headerReqs.size(),
                    headerReqs.data(),
                    MPI_STATUSES_IGNORE
                );
            }
        }
        else if (sizeExchangeReq != MPI_REQUEST_NULL)
        {
            MPI_Wait(&sizeExchangeReq, MPI_STATUS_IGNORE);
        }
        else
        {
            MPI_Alltoall
            (
                sendSizes.data(),
                1,
                MPI_UINT64_T,
                recvSizes.data(),
                1,
                MPI_UINT64_T,
                MPI_COMM_WORLD
            );
        }

        // MPI exchange: all data transfers are bounded chunks.  Two receive
        // buffers provide a small finite pipeline without MPI_THREAD_MULTIPLE.
        std::chrono::steady_clock::time_point tSizeExchangeEnd;
        std::chrono::steady_clock::time_point tRequestPostEnd;
        tSizeExchangeEnd = std::chrono::steady_clock::now();

        const label maxChunkRecords =
            max
            (
                label(1),
                transferChunkBytes_
              / label(sizeof(dsmcParcel::TransferData))
            );
        scalar chunkWaitWall = 0.0;
        scalar chunkDeserializeWall = 0.0;
        std::uint64_t maxIncomingRecords = 0;
        const std::uint64_t recordBytes =
            static_cast<std::uint64_t>
            (sizeof(dsmcParcel::TransferData));
        for (label src = 0; src < nProcs_; ++src)
        {
            if (src == myRank_ || recvSizes[src] == 0) continue;
            if (recvSizes[src] % recordBytes != 0)
            {
                FatalErrorInFunction
                    << "Received byte count "
                    << scalar(recvSizes[src])
                    << " from rank " << src
                    << " is not a multiple of TransferData size "
                    << scalar(recordBytes) << abort(FatalError);
            }
            const std::uint64_t nRecords =
                recvSizes[src] / recordBytes;
            if (nRecords > maxIncomingRecords)
            {
                maxIncomingRecords = nRecords;
            }
        }

        DynamicList<MPI_Request> sendReqs;
        for (label dst = 0; dst < nProcs_; ++dst)
        {
            if (dst == myRank_ || sendSizes[dst] == 0) continue;

            const std::uint64_t nRecords =
                sendSizes[dst]
              / static_cast<std::uint64_t>
                (sizeof(dsmcParcel::TransferData));
            std::uint64_t offset = 0;
            while (offset < nRecords)
            {
                const std::uint64_t chunkLimit =
                    static_cast<std::uint64_t>(maxChunkRecords);
                const std::uint64_t remaining = nRecords - offset;
                const std::uint64_t nChunk =
                    remaining < chunkLimit ? remaining : chunkLimit;
                MPI_Request req;
                MPI_Isend
                (
                    sendTD[dst].data() + label(offset),
                    int(nChunk*sizeof(dsmcParcel::TransferData)),
                    MPI_BYTE,
                    dst,
                    101,
                    MPI_COMM_WORLD,
                    &req
                );
                sendReqs.append(req);
                offset += nChunk;
            }
        }
        tRequestPostEnd = std::chrono::steady_clock::now();

        const std::uint64_t chunkLimit =
            static_cast<std::uint64_t>(maxChunkRecords);
        const label recvChunkRecords =
            label
            (
                maxIncomingRecords < chunkLimit
              ? maxIncomingRecords
              : chunkLimit
            );
        List<List<dsmcParcel::TransferData>> recvChunks(2);
        recvChunks[0].setSize(recvChunkRecords);
        recvChunks[1].setSize(recvChunkRecords);
        List<dsmcParcel*> received(recvChunkRecords);
        MPI_Request recvReqs[2] =
        {
            MPI_REQUEST_NULL,
            MPI_REQUEST_NULL
        };
        std::uint64_t postedRecords[2] = {0, 0};

        for (label src = 0; src < nProcs_; ++src)
        {
            if (src == myRank_ || recvSizes[src] == 0) continue;

            const std::uint64_t nRecords =
                recvSizes[src] / recordBytes;
            std::uint64_t nextRecord = 0;
            label activeReceives = 0;
            label expectedSlot = 0;

            for (label slot = 0; slot < 2 && nextRecord < nRecords; ++slot)
            {
                const std::uint64_t remaining = nRecords - nextRecord;
                const std::uint64_t nChunk =
                    remaining < chunkLimit ? remaining : chunkLimit;
                postedRecords[slot] = nChunk;
                MPI_Irecv
                (
                    recvChunks[slot].data(),
                    int(nChunk*sizeof(dsmcParcel::TransferData)),
                    MPI_BYTE,
                    src,
                    101,
                    MPI_COMM_WORLD,
                    &recvReqs[slot]
                );
                nextRecord += nChunk;
                ++activeReceives;
            }

            while (activeReceives > 0)
            {
                const label completedSlot = expectedSlot;
                const auto tChunkWaitStart =
                    std::chrono::steady_clock::now();
                MPI_Wait
                (
                    &recvReqs[completedSlot],
                    MPI_STATUS_IGNORE
                );
                chunkWaitWall += std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now() - tChunkWaitStart
                ).count();

                const label nChunk = label(postedRecords[completedSlot]);
                const auto tChunkDeserializeStart =
                    std::chrono::steady_clock::now();
#ifdef _OPENMP
                if (maintainMoveOrdered && cloud_.ompNumThreads() > 1)
                {
                    #pragma omp parallel for num_threads(cloud_.ompNumThreads()) schedule(static)
                    for (label j = 0; j < nChunk; ++j)
                    {
                        received[j] = dsmcParcel::unpackTransfer
                        (
                            mesh_,
                            recvChunks[completedSlot][j]
                        );
                    }
                }
                else
#endif
                {
                    for (label j = 0; j < nChunk; ++j)
                    {
                        received[j] = dsmcParcel::unpackTransfer
                        (
                            mesh_,
                            recvChunks[completedSlot][j]
                        );
                    }
                }
                for (label j = 0; j < nChunk; ++j)
                {
                    // Base-class raw append.  The per-parcel cache
                    // bookkeeping in dsmcCloud::addParticle is redundant on
                    // this path: the post-transfer code below re-establishes
                    // the ordered view (setMoveOrderedParcels) or clears all
                    // traversal caches (clearMoveOrderedParcels).
                    cloud_.Cloud<dsmcParcel>::addParticle(received[j]);
                    kept.append(received[j]);
                }
                nRecv += nChunk;
                chunkDeserializeWall += std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now()
                  - tChunkDeserializeStart
                ).count();

                if (nextRecord < nRecords)
                {
                    const std::uint64_t remaining =
                        nRecords - nextRecord;
                    const std::uint64_t nNextChunk =
                        remaining < chunkLimit ? remaining : chunkLimit;
                    postedRecords[completedSlot] = nNextChunk;
                    MPI_Irecv
                    (
                        recvChunks[completedSlot].data(),
                        int(nNextChunk*sizeof(dsmcParcel::TransferData)),
                        MPI_BYTE,
                        src,
                        101,
                        MPI_COMM_WORLD,
                        &recvReqs[completedSlot]
                    );
                    nextRecord += nNextChunk;
                }
                else
                {
                    recvReqs[completedSlot] = MPI_REQUEST_NULL;
                    postedRecords[completedSlot] = 0;
                    --activeReceives;
                }
                expectedSlot = 1 - expectedSlot;
            }
        }

        if (sendReqs.size() > 0)
        {
            const auto tSendWaitStart =
                std::chrono::steady_clock::now();
            MPI_Waitall
            (
                sendReqs.size(),
                sendReqs.data(),
                MPI_STATUSES_IGNORE
            );
            chunkWaitWall += std::chrono::duration<scalar>
            (
                std::chrono::steady_clock::now() - tSendWaitStart
            ).count();
        }
        const auto tTransferEnd = std::chrono::steady_clock::now();

        if (maintainMoveOrdered && kept.size() == cloud_.size())
        {
            cloud_.setMoveOrderedParcels(kept);
        }
        else
        {
            cloud_.clearMoveOrderedParcels();
        }
        const auto tPostEnd = std::chrono::steady_clock::now();

        // Exchange candidate counts for optional offload diagnostics/planners.
        // Current replicated-mesh DLB does not consume allProcCandidates_.
        {
            if (gatherCandidates_)
            {
                const labelList& nCandPerCell = cloud_.nCandidatesPerCell();
                label localCands = 0;
                for (label i = 0; i < nCandPerCell.size(); ++i)
                    localCands += nCandPerCell[i];
                allProcCandidates_.setSize(nProcs_, 0);
                MPI_Allgather
                (
                    &localCands,
                    1,
                    MPI_INT,
                    allProcCandidates_.data(),
                    1,
                    MPI_INT,
                    MPI_COMM_WORLD
                );
            }
        }

        const auto tEnd = std::chrono::steady_clock::now();
        migrationWallTime_ += std::chrono::duration<scalar>(tEnd - tStart).count();
        migrationPackWallTime_ +=
            std::chrono::duration<scalar>(tPackEnd - tStart).count();
        migrationLocalPrepWallTime_ +=
            std::chrono::duration<scalar>(tLocalPrepEnd - tPackEnd).count();
        migrationSizeExchangeWallTime_ +=
            std::chrono::duration<scalar>
            (
                tSizeExchangeEnd - tLocalPrepEnd
            ).count();
        migrationRequestPostWallTime_ +=
            std::chrono::duration<scalar>
            (
                tRequestPostEnd - tSizeExchangeEnd
            ).count();
        migrationWaitWallTime_ += chunkWaitWall;
        migrationDeserializeWallTime_ += chunkDeserializeWall;
        migrationPostWallTime_ +=
            std::chrono::duration<scalar>(tPostEnd - tTransferEnd).count();
        migrationCandidateGatherWallTime_ +=
            std::chrono::duration<scalar>(tEnd - tPostEnd).count();
        const scalar migrationWall =
            std::chrono::duration<scalar>(tEnd - tStart).count();
        const scalar migrationComm =
            std::chrono::duration<scalar>
            (
                tSizeExchangeEnd - tLocalPrepEnd
            ).count()
          + std::chrono::duration<scalar>
            (
                tRequestPostEnd - tSizeExchangeEnd
            ).count()
          + chunkWaitWall
          + std::chrono::duration<scalar>(tEnd - tPostEnd).count();
        addMigrationProfileSample(migrationWall, migrationComm);
        ++migrationCalls_;
        totalParcelsMigrated_ += nMigratedOut;

        if (cloud_.emitStepDiagnostics())
        {
            label globalMigrated = 0;
            MPI_Allreduce(&nMigratedOut, &globalMigrated, 1, MPI_INT,
                          MPI_SUM, MPI_COMM_WORLD);
            Info<< "Replicated mesh migration[" << migrationCalls_ << "]: "
                << "rank " << myRank_
                << " sent " << nMigratedOut << " parcels, received " << nRecv
                << " (global " << globalMigrated << ")"
                << " wall " << std::chrono::duration<scalar>(tEnd - tStart).count() << "s"
                << " [pack " << std::chrono::duration<scalar>
                    (tPackEnd - tStart).count() << "s"
                << " prep " << std::chrono::duration<scalar>
                    (tLocalPrepEnd - tPackEnd).count() << "s"
                << " sizeX " << std::chrono::duration<scalar>
                    (tSizeExchangeEnd - tLocalPrepEnd).count() << "s"
                << " postReq " << std::chrono::duration<scalar>
                    (tRequestPostEnd - tSizeExchangeEnd).count() << "s"
                << " wait " << chunkWaitWall << "s"
                << " deser " << chunkDeserializeWall << "s"
                << " post " << std::chrono::duration<scalar>
                    (tPostEnd - tTransferEnd).count() << "s"
                << " candX " << std::chrono::duration<scalar>
                    (tEnd - tPostEnd).count() << "s]"
                << endl;
        }
        return;
    }

    // ---- Stream-based transfer path (original) ----
    PtrList<OStringStream> sendStreams(nProcs_);
    for (label i = 0; i < nProcs_; ++i)
        if (i != myRank_) sendStreams.set(i, new OStringStream(IOstream::BINARY));

    DynamicList<dsmcParcel*> toDelete(cloud_.size() / 4);
    DynamicList<dsmcParcel*> kept(cloud_.size());
    label nMigratedOut = 0;
    const bool maintainMoveOrdered =
        cloud_.openmpEnabled()
     && cloud_.openmpMoveEnabled()
     && cloud_.ompNumThreads() > 1;
    const bool useOrderedTraversal =
        maintainMoveOrdered
     && cloud_.hasMoveOrderedParcels()
     && cloud_.moveOrderedParcels().size() == cloud_.size();

    if (useOrderedTraversal)
    {
        const auto& ordered = cloud_.moveOrderedParcels();
        for (label i = 0; i < ordered.size(); ++i)
        {
            dsmcParcel& p = *ordered[i];
            const label dstRank = cellOwner_[p.cell()];
            if (dstRank != myRank_)
            {
                p.writeBinaryFast(sendStreams[dstRank]);
                toDelete.append(ordered[i]);
                ++nMigratedOut;
            }
            else
            {
                kept.append(ordered[i]);
            }
        }
    }
    else
    {
        forAllIter(Cloud<dsmcParcel>, cloud_, iter)
        {
            dsmcParcel& p = iter();
            const label dstRank = cellOwner_[p.cell()];
            if (dstRank != myRank_)
            {
                p.writeBinaryFast(sendStreams[dstRank]);
                toDelete.append(&p);
                ++nMigratedOut;
            }
            else
            {
                kept.append(&p);
            }
        }
    }
    const auto t1 = std::chrono::steady_clock::now();

    // ---- Phase 2: delete migrated parcels and extract buffers ----------------
    forAll(toDelete, i) { cloud_.deleteParcel(toDelete[i]); }

    List<DynamicList<char>> sendBufs(nProcs_);
    labelList sendSizes(nProcs_, 0);
    for (label i = 0; i < nProcs_; ++i)
    {
        if (i != myRank_ && sendStreams.set(i))
        {
            const string s = sendStreams[i].str();
            sendBufs[i].setSize(s.size());
            forAll(s, si)
            {
                sendBufs[i][si] = s[si];
            }
            sendSizes[i] = sendBufs[i].size();
        }
    }
    sendStreams.clear();
    const auto t2 = std::chrono::steady_clock::now();

    // ---- Phase 3: point-to-point exchange (no global sync) -------------------
    // For 2 ranks: single MPI_Sendrecv for sizes, then for data.
    // For N ranks: Alltoall sizes then Isend/Irecv data.
    label nRecv = 0;
    List<char> recvBuf;

    if (nProcs_ == 2)
    {
        const label peer = 1 - myRank_;
        const label sendSize = sendSizes[peer];

        // Exchange sizes (non-blocking friendly — Sendrecv doesn't require prior sync)
        label recvSize = 0;
        MPI_Sendrecv(&sendSize, 1, MPI_INT, peer, 1,
                     &recvSize, 1, MPI_INT, peer, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Exchange data
        recvBuf.setSize(recvSize);
        MPI_Sendrecv(
            sendBufs[peer].data(), sendSize, MPI_BYTE, peer, 2,
            recvBuf.data(), recvSize, MPI_BYTE, peer, 2,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Deserialize
        if (recvSize > 0)
        {
            ISpanStream is(recvBuf.data(), recvSize, IOstream::BINARY);
            while (!is.eof())
            {
                auto* newp = new dsmcParcel(mesh_, is);
                cloud_.addParticle(newp);
                kept.append(newp);
                ++nRecv;
            }
        }
    }
    else if (useNoAlltoall_)
    {
        // No-Alltoall: post large Irecv per peer, then Isend, then Waitall+Get_count
        const label minBufSize = 65536;
        DynamicList<MPI_Request> allReqs(2 * (nProcs_ - 1));

        for (label i = 0; i < nProcs_; ++i)
        {
            if (i == myRank_) continue;
            label needed = prevRecvSizes_[i] * 2 + minBufSize;
            if (needed > perPeerRecvCapacity_[i])
            {
                perPeerRecvCapacity_[i] = needed;
                perPeerRecvBuf_[i].setSize(needed);
            }
            MPI_Request req;
            MPI_Irecv(perPeerRecvBuf_[i].data(), perPeerRecvCapacity_[i],
                      MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
            allReqs.append(req);
        }
        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && sendSizes[i] > 0)
            {
                MPI_Request req;
                MPI_Isend(sendBufs[i].data(), sendSizes[i],
                          MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
                allReqs.append(req);
            }
        }

        List<MPI_Status> statuses(allReqs.size());
        if (allReqs.size() > 0)
            MPI_Waitall(allReqs.size(), allReqs.data(), statuses.data());

        // Deserialize: first nProcs_-1 requests are Irecv
        label reqIdx = 0;
        for (label i = 0; i < nProcs_; ++i)
        {
            if (i == myRank_) continue;
            int actualBytes = 0;
            MPI_Get_count(&statuses[reqIdx], MPI_BYTE, &actualBytes);
            prevRecvSizes_[i] = actualBytes;
            if (actualBytes > 0)
            {
                ISpanStream is(perPeerRecvBuf_[i].data(), actualBytes,
                               IOstream::BINARY);
                while (!is.eof())
                {
                    auto* newp = new dsmcParcel(mesh_, is);
                    cloud_.addParticle(newp);
                    kept.append(newp);
                    ++nRecv;
                }
            }
            ++reqIdx;
        }
    }
    else
    {
        // Multi-rank: exchange sizes via Alltoall, then Isend/Irecv
        labelList recvSizes(nProcs_, 0);
        MPI_Alltoall(sendSizes.data(), 1, MPI_INT,
                     recvSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

        label totalRecvSize = 0;
        for (label i = 0; i < nProcs_; ++i)
            if (i != myRank_) totalRecvSize += recvSizes[i];
        recvBuf.setSize(totalRecvSize);

        labelList recvOffsets(nProcs_, 0);
        {
            label off = 0;
            for (label i = 0; i < nProcs_; ++i)
            {
                recvOffsets[i] = off;
                if (i != myRank_) off += recvSizes[i];
            }
        }

        DynamicList<MPI_Request> allReqs(2 * (nProcs_ - 1));
        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && recvSizes[i] > 0)
            {
                MPI_Request req;
                MPI_Irecv(recvBuf.data() + recvOffsets[i], recvSizes[i],
                          MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
                allReqs.append(req);
            }
        }
        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && sendSizes[i] > 0)
            {
                MPI_Request req;
                MPI_Isend(sendBufs[i].data(), sendSizes[i],
                          MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
                allReqs.append(req);
            }
        }
        if (allReqs.size() > 0)
            MPI_Waitall(allReqs.size(), allReqs.data(), MPI_STATUSES_IGNORE);

        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && recvSizes[i] > 0)
            {
                ISpanStream is(recvBuf.data() + recvOffsets[i], recvSizes[i],
                               IOstream::BINARY);
                while (!is.eof())
                {
                    auto* newp = new dsmcParcel(mesh_, is);
                    cloud_.addParticle(newp);
                    kept.append(newp);
                    ++nRecv;
                }
            }
        }
    }
    sendBufs.clear();
    const auto t3 = std::chrono::steady_clock::now();

    // ---- Phase 4: set moveOrderedParcels_ directly (no linked-list rebuild) ---
    if (maintainMoveOrdered && kept.size() == cloud_.size())
    {
        cloud_.setMoveOrderedParcels(kept);
    }
    else
    {
        cloud_.clearMoveOrderedParcels();
    }

    // ---- Phase 5: optional per-rank candidate counts for offload planner ----
    if (gatherCandidates_)
    {
        const labelList& nCandPerCell = cloud_.nCandidatesPerCell();
        label localCands = 0;
        for (label i = 0; i < nCandPerCell.size(); ++i)
        {
            localCands += nCandPerCell[i];
        }
        allProcCandidates_.setSize(nProcs_, 0);
        MPI_Allgather
        (
            &localCands,
            1,
            MPI_INT,
            allProcCandidates_.data(),
            1,
            MPI_INT,
            MPI_COMM_WORLD
        );
    }

    const auto tEnd = std::chrono::steady_clock::now();

    // ---- Timing -----------------------------------------------------------------
    const scalar migrationWall =
        std::chrono::duration<scalar>(tEnd - tStart).count();
    const scalar migrationComm =
        std::chrono::duration<scalar>(t3 - t2).count();
    migrationWallTime_ += migrationWall;
    addMigrationProfileSample(migrationWall, migrationComm);
    ++migrationCalls_;

    totalParcelsMigrated_ += nMigratedOut;

    if (cloud_.emitStepDiagnostics())
    {
        label globalMigrated = 0;
        MPI_Allreduce(&nMigratedOut, &globalMigrated, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        Info<< "Replicated mesh migration[" << migrationCalls_ << "]: "
            << "rank " << myRank_
            << " sent " << nMigratedOut << " parcels, received " << nRecv
            << " (global " << globalMigrated << ")"
            << " wall " << std::chrono::duration<scalar>(tEnd - tStart).count() << "s"
            << " [pack " << std::chrono::duration<scalar>(t2 - t1).count() << "s"
            << " MPI "  << std::chrono::duration<scalar>(t3 - t2).count() << "s"
            << " deser "<< std::chrono::duration<scalar>(tEnd - t3).count() << "s]"
            << endl;
    }
}


// ============================================================================
// migrateBegin — async Phase A: serialize + delete + Isend/Irecv (non-blocking)
// ============================================================================

void dsmcReplicatedMesh::migrateBegin()
{
    if (!active_) return;

    if (useFlatTransfer_)
    {
        migrateParticlesByCellOwner();
        return;
    }

    const auto tStart = std::chrono::steady_clock::now();

    // Phase 1: partition kept vs migrating
    PtrList<OStringStream> sendStreams(nProcs_);
    for (label i = 0; i < nProcs_; ++i)
        if (i != myRank_) sendStreams.set(i, new OStringStream(IOstream::BINARY));

    DynamicList<dsmcParcel*> toDelete(cloud_.size() / 4);
    DynamicList<dsmcParcel*> kept(cloud_.size());
    label nMigratedOut = 0;
    const bool maintainMoveOrdered =
        cloud_.openmpEnabled()
     && cloud_.openmpMoveEnabled()
     && cloud_.ompNumThreads() > 1;
    const bool useOrderedTraversal =
        maintainMoveOrdered
     && cloud_.hasMoveOrderedParcels()
     && cloud_.moveOrderedParcels().size() == cloud_.size();

    if (useOrderedTraversal)
    {
        const auto& ordered = cloud_.moveOrderedParcels();
        for (label i = 0; i < ordered.size(); ++i)
        {
            dsmcParcel& p = *ordered[i];
            const label dstRank = cellOwner_[p.cell()];
            if (dstRank != myRank_)
            {
                p.writeBinaryFast(sendStreams[dstRank]);
                toDelete.append(ordered[i]);
                ++nMigratedOut;
            }
            else
            {
                kept.append(ordered[i]);
            }
        }
    }
    else
    {
        forAllIter(Cloud<dsmcParcel>, cloud_, iter)
        {
            dsmcParcel& p = iter();
            const label dstRank = cellOwner_[p.cell()];
            if (dstRank != myRank_)
            {
                p.writeBinaryFast(sendStreams[dstRank]);
                toDelete.append(&p);
                ++nMigratedOut;
            }
            else
            {
                kept.append(&p);
            }
        }
    }

    // Phase 2: delete + extract buffers
    forAll(toDelete, i) { cloud_.deleteParcel(toDelete[i]); }

    List<DynamicList<char>> sendBufs(nProcs_);
    labelList sendSizes(nProcs_, 0);
    for (label i = 0; i < nProcs_; ++i)
    {
        if (i != myRank_ && sendStreams.set(i))
        {
            const string s = sendStreams[i].str();
            sendBufs[i].setSize(s.size());
            forAll(s, si)
            {
                sendBufs[i][si] = s[si];
            }
            sendSizes[i] = sendBufs[i].size();
        }
    }
    sendStreams.clear();

    // Phase 3: non-blocking MPI exchange
    asyncReqs_.clear();

    if (nProcs_ == 2)
    {
        const label peer = 1 - myRank_;
        const label sendSize = sendSizes[peer];

        // Exchange sizes (blocking — cheap)
        label recvSize = 0;
        MPI_Sendrecv(&sendSize, 1, MPI_INT, peer, 1,
                     &recvSize, 1, MPI_INT, peer, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        asyncRecvSize_ = recvSize;
        asyncRecvBuf_.setSize(recvSize);

        // Non-blocking data exchange
        if (recvSize > 0)
        {
            MPI_Request req;
            MPI_Irecv(asyncRecvBuf_.data(), recvSize, MPI_BYTE, peer, 2,
                      MPI_COMM_WORLD, &req);
            asyncReqs_.append(req);
        }
        if (sendSize > 0)
        {
            asyncSendBuf_.setSize(sendSize);
            std::memcpy(asyncSendBuf_.data(), sendBufs[peer].data(), sendSize);
            MPI_Request req;
            MPI_Isend(asyncSendBuf_.data(), sendSize, MPI_BYTE, peer, 2,
                      MPI_COMM_WORLD, &req);
            asyncReqs_.append(req);
        }
    }
    else if (useNoAlltoall_)
    {
        // No-Alltoall path: use per-peer pre-allocated buffers.
        const label minBufSize = 65536;

        for (label i = 0; i < nProcs_; ++i)
        {
            if (i == myRank_) continue;
            label needed = prevRecvSizes_[i] * 2 + minBufSize;
            if (needed > perPeerRecvCapacity_[i])
            {
                perPeerRecvCapacity_[i] = needed;
                perPeerRecvBuf_[i].setSize(needed);
            }
            MPI_Request req;
            MPI_Irecv(perPeerRecvBuf_[i].data(), perPeerRecvCapacity_[i],
                      MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
            asyncReqs_.append(req);
        }

        label totalSendSize = 0;
        for (label i = 0; i < nProcs_; ++i) totalSendSize += sendSizes[i];
        asyncSendBuf_.setSize(totalSendSize);
        labelList sendOffsets(nProcs_, 0);
        {
            label off = 0;
            for (label i = 0; i < nProcs_; ++i)
            {
                sendOffsets[i] = off;
                if (i != myRank_ && sendSizes[i] > 0)
                {
                    std::memcpy(asyncSendBuf_.data() + off,
                                sendBufs[i].data(), sendSizes[i]);
                }
                off += sendSizes[i];
            }
        }

        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && sendSizes[i] > 0)
            {
                MPI_Request req;
                MPI_Isend(asyncSendBuf_.data() + sendOffsets[i], sendSizes[i],
                          MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
                asyncReqs_.append(req);
            }
        }
    }
    else
    {
        // Multi-rank: exchange sizes via Alltoall, then Isend/Irecv
        labelList recvSizes(nProcs_, 0);
        MPI_Alltoall(sendSizes.data(), 1, MPI_INT,
                     recvSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

        label totalRecvSize = 0;
        for (label i = 0; i < nProcs_; ++i)
            if (i != myRank_) totalRecvSize += recvSizes[i];
        asyncRecvBuf_.setSize(totalRecvSize);
        asyncRecvSize_ = totalRecvSize;

        labelList recvOffsets(nProcs_, 0);
        {
            label off = 0;
            for (label i = 0; i < nProcs_; ++i)
            {
                recvOffsets[i] = off;
                if (i != myRank_) off += recvSizes[i];
            }
        }

        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && recvSizes[i] > 0)
            {
                MPI_Request req;
                MPI_Irecv(asyncRecvBuf_.data() + recvOffsets[i], recvSizes[i],
                          MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
                asyncReqs_.append(req);
            }
        }

        // Consolidate send buffers (must stay alive until Waitall)
        label totalSendSize = 0;
        for (label i = 0; i < nProcs_; ++i) totalSendSize += sendSizes[i];
        asyncSendBuf_.setSize(totalSendSize);
        labelList sendOffsets(nProcs_, 0);
        {
            label off = 0;
            for (label i = 0; i < nProcs_; ++i)
            {
                sendOffsets[i] = off;
                if (i != myRank_ && sendSizes[i] > 0)
                {
                    std::memcpy(asyncSendBuf_.data() + off,
                                sendBufs[i].data(), sendSizes[i]);
                }
                off += sendSizes[i];
            }
        }

        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && sendSizes[i] > 0)
            {
                MPI_Request req;
                MPI_Isend(asyncSendBuf_.data() + sendOffsets[i], sendSizes[i],
                          MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
                asyncReqs_.append(req);
            }
        }
    }

    // Set moveOrdered with kept parcels; received parcels are appended in
    // migrateFinish(). This keeps the next occupancy build off the linked list.
    if (maintainMoveOrdered && kept.size() == cloud_.size())
    {
        cloud_.setMoveOrderedParcels(kept);
    }
    else
    {
        cloud_.clearMoveOrderedParcels();
    }

    asyncMigrationPending_ = true;
    totalParcelsMigrated_ += nMigratedOut;
    ++migrationCalls_;

    const scalar migrationWall = std::chrono::duration<scalar>
    (
        std::chrono::steady_clock::now() - tStart
    ).count();
    migrationWallTime_ += migrationWall;
    addMigrationProfileSample(migrationWall, migrationWall);
}


// ============================================================================
// migrateFinish — async Phase A: Waitall + deserialize received parcels
// ============================================================================

void dsmcReplicatedMesh::migrateFinish()
{
    if (!asyncMigrationPending_) return;

    const bool maintainMoveOrdered =
        cloud_.openmpEnabled()
     && cloud_.openmpMoveEnabled()
     && cloud_.ompNumThreads() > 1;
    const auto tFinish0 = std::chrono::steady_clock::now();

    if (useNoAlltoall_ && asyncReqs_.size() > 0)
    {
        // No-Alltoall path: use MPI_Get_count per peer
        List<MPI_Status> statuses(asyncReqs_.size());
        MPI_Waitall(asyncReqs_.size(), asyncReqs_.data(), statuses.data());

        DynamicList<dsmcParcel*> received(1024);

        // First nProcs_-1 requests are Irecv (one per peer != myRank_)
        label reqIdx = 0;
        for (label i = 0; i < nProcs_; ++i)
        {
            if (i == myRank_) continue;
            int actualBytes = 0;
            MPI_Get_count(&statuses[reqIdx], MPI_BYTE, &actualBytes);
            prevRecvSizes_[i] = actualBytes;
            if (actualBytes > 0)
            {
                ISpanStream is(perPeerRecvBuf_[i].data(), actualBytes,
                               IOstream::BINARY);
                while (!is.eof())
                {
                    auto* newp = new dsmcParcel(mesh_, is);
                    cloud_.addParticle(newp);
                    received.append(newp);
                }
            }
            ++reqIdx;
        }
        if
        (
            received.size() > 0
         && maintainMoveOrdered
         && cloud_.hasMoveOrderedParcels()
        )
        {
            cloud_.appendBatchToMoveOrdered(received);
        }
        else if (received.size() > 0)
        {
            cloud_.clearMoveOrderedParcels();
        }
    }
    else if (asyncReqs_.size() > 0)
    {
        MPI_Waitall(asyncReqs_.size(), asyncReqs_.data(), MPI_STATUSES_IGNORE);

        if (asyncRecvSize_ > 0)
        {
            DynamicList<dsmcParcel*> received(asyncRecvSize_ / 100);
            ISpanStream is(asyncRecvBuf_.data(), asyncRecvSize_,
                           IOstream::BINARY);
            while (!is.eof())
            {
                auto* newp = new dsmcParcel(mesh_, is);
                cloud_.addParticle(newp);
                received.append(newp);
            }
            if
            (
                received.size() > 0
             && maintainMoveOrdered
             && cloud_.hasMoveOrderedParcels()
            )
            {
                cloud_.appendBatchToMoveOrdered(received);
            }
            else if (received.size() > 0)
            {
                cloud_.clearMoveOrderedParcels();
            }
        }
    }

    asyncReqs_.clear();
    asyncSendBuf_.clear();
    asyncRecvBuf_.clear();
    asyncMigrationPending_ = false;
    addMigrationProfileSample
    (
        0.0,
        std::chrono::duration<scalar>
        (
            std::chrono::steady_clock::now() - tFinish0
        ).count()
    );
}


// ============================================================================
// updateParticleCounts — Phase A
// ============================================================================

void dsmcReplicatedMesh::updateParticleCounts()
{
    const auto tStart = std::chrono::steady_clock::now();
    localParticleCount_ = cloud_.size();
    allParticleCounts_.setSize(nProcs_);
    MPI_Allgather(&localParticleCount_, 1, MPI_INT,
                  allParticleCounts_.data(), 1, MPI_INT, MPI_COMM_WORLD);
    updateParticleCountsWallTime_ += std::chrono::duration<scalar>
    (
        std::chrono::steady_clock::now() - tStart
    ).count();
    addMigrationProfileSample
    (
        0.0,
        std::chrono::duration<scalar>
        (
            std::chrono::steady_clock::now() - tStart
        ).count()
    );
}


// ============================================================================
// reassignCellOwner — Phase B: manual block/reverse-block
// ============================================================================

void dsmcReplicatedMesh::reassignCellOwner()
{
    if (!active_ || nProcs_ < 2) return;
    const label nCells = mesh_.nCells();
    labelList oldCellOwner(cellOwner_);

    if (rebalanceCount_ % 2 == 0)
    {
        for (label cellI = 0; cellI < nCells; ++cellI)
            cellOwner_[cellI] = nProcs_ - 1 - oldCellOwner[cellI];
    }
    else
    {
        for (label cellI = 0; cellI < nCells; ++cellI)
            cellOwner_[cellI] = label((scalar(cellI)/scalar(nCells))*scalar(nProcs_));
        if (cellOwner_[nCells-1] >= nProcs_) cellOwner_[nCells-1] = nProcs_ - 1;
    }

    rebuildMyCells();
    localMesh_.build(myCells_);
    label nChanged = 0;
    for (label cellI = 0; cellI < nCells; ++cellI)
        if (cellOwner_[cellI] != oldCellOwner[cellI]) ++nChanged;

    totalCellsChanged_ += nChanged;
    ++rebalanceCount_;

    Info<< "Phase B[" << rebalanceCount_ << "]: cellOwner_ reassigned. "
        << nChanged << " / " << nCells << " cells changed (R_changed="
        << scalar(nChanged)/scalar(nCells) << ")" << endl;
}


// ============================================================================
// gatherParcelsToRank0 — collect all particles on rank 0 for writing
// ============================================================================

void dsmcReplicatedMesh::gatherParcelsToRank0()
{
    if (!active_ || nProcs_ <= 1) return;

    cloud_.clearMoveOrderedParcels();

    if (useFlatTransfer_)
    {
        if (myRank_ != 0)
        {
            DynamicList<dsmcParcel::TransferData> sendTD(cloud_.size());
            forAllIter(Cloud<dsmcParcel>, cloud_, iter)
            {
                sendTD.append(dsmcParcel::TransferData());
                iter().packTransfer(sendTD.last());
            }

            const std::uint64_t sendBytes =
                static_cast<std::uint64_t>(sendTD.size())
              * static_cast<std::uint64_t>
                (sizeof(dsmcParcel::TransferData));
            MPI_Send
            (
                &sendBytes, 1, MPI_UINT64_T, 0, 10,
                MPI_COMM_WORLD
            );

            const label maxChunkRecords =
                max
                (
                    label(1),
                    transferChunkBytes_
                  / label(sizeof(dsmcParcel::TransferData))
                );
            DynamicList<MPI_Request> sendReqs;
            std::uint64_t offset = 0;
            const std::uint64_t nRecords = sendTD.size();
            while (offset < nRecords)
            {
                const std::uint64_t chunkLimit =
                    static_cast<std::uint64_t>(maxChunkRecords);
                const std::uint64_t remaining = nRecords - offset;
                const std::uint64_t nChunk =
                    remaining < chunkLimit ? remaining : chunkLimit;
                MPI_Request req;
                MPI_Isend
                (
                    sendTD.data() + label(offset),
                    int(nChunk*sizeof(dsmcParcel::TransferData)),
                    MPI_BYTE,
                    0,
                    11,
                    MPI_COMM_WORLD,
                    &req
                );
                sendReqs.append(req);
                offset += nChunk;
            }
            if (sendReqs.size() > 0)
            {
                MPI_Waitall
                (
                    sendReqs.size(),
                    sendReqs.data(),
                    MPI_STATUSES_IGNORE
                );
            }

            cloud_.clear();
        }
        else
        {
            for (label src = 1; src < nProcs_; ++src)
            {
                std::uint64_t recvBytes = 0;
                MPI_Recv
                (
                    &recvBytes,
                    1,
                    MPI_UINT64_T,
                    src,
                    10,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );

                const std::uint64_t recordBytes =
                    static_cast<std::uint64_t>
                    (sizeof(dsmcParcel::TransferData));
                if (recvBytes % recordBytes != 0)
                {
                    FatalErrorInFunction
                        << "Received output byte count "
                        << scalar(recvBytes)
                        << " from rank " << src
                        << " is not a multiple of TransferData size "
                        << scalar(recordBytes) << abort(FatalError);
                }

                std::uint64_t receivedRecords = 0;
                const std::uint64_t nRecords = recvBytes / recordBytes;
                if (nRecords == 0) continue;

                const label maxChunkRecords =
                    max
                    (
                        label(1),
                        transferChunkBytes_
                      / label(sizeof(dsmcParcel::TransferData))
                    );
                const std::uint64_t chunkLimit =
                    static_cast<std::uint64_t>(maxChunkRecords);
                const label recvChunkRecords =
                    label
                    (
                        nRecords < chunkLimit ? nRecords : chunkLimit
                    );
                List<dsmcParcel::TransferData> recvChunk(recvChunkRecords);
                while (receivedRecords < nRecords)
                {
                    const std::uint64_t remaining =
                        nRecords - receivedRecords;
                    const std::uint64_t nChunk =
                        remaining < chunkLimit ? remaining : chunkLimit;
                    MPI_Request req;
                    MPI_Irecv
                    (
                        recvChunk.data(),
                        int(nChunk*sizeof(dsmcParcel::TransferData)),
                        MPI_BYTE,
                        src,
                        11,
                        MPI_COMM_WORLD,
                        &req
                    );
                    MPI_Wait(&req, MPI_STATUS_IGNORE);

                    for (std::uint64_t i = 0; i < nChunk; ++i)
                    {
                        cloud_.addParticle
                        (
                            dsmcParcel::unpackTransfer
                            (mesh_, recvChunk[label(i)])
                        );
                    }
                    receivedRecords += nChunk;
                }
            }
        }

        return;
    }

    if (myRank_ != 0)
    {
        OStringStream os(IOstream::BINARY);
        forAllIter(Cloud<dsmcParcel>, cloud_, iter)
        {
            iter().writeBinaryFast(os);
        }
        cloud_.clear();

        const string s = os.str();
        DynamicList<char> buf(s.size());
        forAll(s, si)
        {
            buf[si] = s[si];
        }
        label sendSize = buf.size();
        MPI_Send(&sendSize, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
        if (sendSize > 0)
        {
            MPI_Send(buf.data(), sendSize, MPI_BYTE, 0, 11, MPI_COMM_WORLD);
        }
    }
    else
    {
        for (label src = 1; src < nProcs_; ++src)
        {
            label recvSize = 0;
            MPI_Recv(&recvSize, 1, MPI_INT, src, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (recvSize > 0)
            {
                List<char> recvBuf(recvSize);
                MPI_Recv(recvBuf.data(), recvSize, MPI_BYTE, src, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                ISpanStream is(recvBuf.data(), recvSize, IOstream::BINARY);
                while (!is.eof())
                {
                    auto* newp = new dsmcParcel(mesh_, is);
                    cloud_.addParticle(newp);
                }
            }
        }
    }
}

// ============================================================================
// printPerRankStepWall — per-rank iteration wall diagnostics at check time
// ============================================================================

void dsmcReplicatedMesh::printPerRankStepWall() const
{
    if (!triggerWallTime_) return;

    scalarList wallPerRank(nProcs_, 0.0);
    wallPerRank[myRank_] = lastStepWallTime_;
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        wallPerRank.data(),
        nProcs_,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );

    scalarList workPerRank(nProcs_, 0.0);
    workPerRank[myRank_] = max(lastStepWallTime_ - lastStepWait_, scalar(0));
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        workPerRank.data(),
        nProcs_,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );

    if (myRank_ == 0)
    {
        Info<< "Phase C per-rank step wall [s]:";
        forAll(wallPerRank, i)
        {
            Info<< (i ? " " : "") << wallPerRank[i];
        }
        Info<< nl << "Phase C per-rank step work [s]:";
        forAll(workPerRank, i)
        {
            Info<< (i ? " " : "") << workPerRank[i];
        }
        Info<< nl;
    }
}


// ============================================================================
// report — profiling summary
// ============================================================================

void dsmcReplicatedMesh::report() const
{
    if (!active_) return;

    const auto tReport0 = std::chrono::steady_clock::now();
    scalarList allTimes(nProcs_, 0.0);
    allTimes[myRank_] = evolveStepTime_;
    MPI_Allreduce(MPI_IN_PLACE, allTimes.data(), nProcs_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    const scalar commDLBLocal =
        dlbMigrationCommWallTime_
      + autoRebalanceExplicitCommWallTime_;
    const scalar commOutputLocal =
        outputMigrationCommWallTime_
      + outputGatherCommWallTime_;
    const scalar commTotalLocal =
        regularMigrationCommWallTime_
      + commDLBLocal
      + commOutputLocal
      + manualMigrationCommWallTime_
      + initialMigrationCommWallTime_
      + otherMigrationCommWallTime_;
    const scalar dlbNonCommLocal = max
    (
        autoRebalanceWallTime_
      - autoRebalanceExplicitCommWallTime_
      - dlbMigrationCommWallTime_,
        scalar(0.0)
    );

    const int nProfileValues = 35;
    scalar profileValues[nProfileValues] =
    {
        migrationWallTime_,
        migrationCommWallTime_,
        regularMigrationWallTime_,
        regularMigrationCommWallTime_,
        dlbMigrationWallTime_,
        dlbMigrationCommWallTime_,
        outputMigrationWallTime_,
        outputMigrationCommWallTime_,
        outputGatherCommWallTime_,
        manualMigrationWallTime_,
        manualMigrationCommWallTime_,
        initialMigrationWallTime_,
        initialMigrationCommWallTime_,
        otherMigrationWallTime_,
        otherMigrationCommWallTime_,
        updateParticleCountsWallTime_,
        autoRebalanceWallTime_,
        autoRebalanceCheckWallTime_,
        autoRebalanceRepartWallTime_,
        autoRebalanceMigrationWallTime_,
        autoRebalanceWriteWallTime_,
        autoRebalancePostDiagWallTime_,
        autoRebalanceExplicitCommWallTime_,
        commTotalLocal,
        commDLBLocal,
        commOutputLocal,
        dlbNonCommLocal,
        migrationPackWallTime_,
        migrationLocalPrepWallTime_,
        migrationSizeExchangeWallTime_,
        migrationRequestPostWallTime_,
        migrationWaitWallTime_,
        migrationDeserializeWallTime_,
        migrationPostWallTime_,
        migrationCandidateGatherWallTime_
    };

    MPI_Allreduce
    (
        MPI_IN_PLACE,
        profileValues,
        nProfileValues,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    const scalar profileReportComm =
        std::chrono::duration<scalar>
        (
            std::chrono::steady_clock::now() - tReport0
        ).count();

    // Per-rank full-iteration wall (§5.6): slot-style SUM reduce — must run
    // on every rank BEFORE the rank0 early return below.
    scalarList allStepWalls(nProcs_, 0.0);
    allStepWalls[myRank_] = stepFullWallCum_;
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        allStepWalls.data(),
        nProcs_,
        MPI_DOUBLE,
        MPI_SUM,
        MPI_COMM_WORLD
    );

    if (myRank_ != 0) return;

    Info<< nl << "Replicated mesh comm/DLB profile v2:" << nl
        << "    migration calls              = " << migrationCalls_ << nl
        << "    local particles (final)      = " << localParticleCount_ << nl
        << "    comm total max [s]           = " << profileValues[23] << nl
        << "    comm_regular_migration [s]   = " << profileValues[3] << nl
        << "    comm_dlb [s]                 = " << profileValues[24] << nl
        << "    comm_output [s]              = " << profileValues[25] << nl
        << "    comm_manual_rebalance [s]    = " << profileValues[10] << nl
        << "    comm_initial_distribution [s]= " << profileValues[12] << nl
        << "    comm_other [s]               = " << profileValues[14] << nl
        << "    comm_profile_report local [s]= " << profileReportComm << nl
        << "    migration wall max [s]       = " << profileValues[0] << nl
        << "    migration comm max [s]       = " << profileValues[1] << nl
        << "    updateParticleCounts max [s] = " << profileValues[15] << nl
        << "    regular migration wall [s]   = " << profileValues[2] << nl
        << "    migration pack max [s]       = " << profileValues[27] << nl
        << "    migration local prep max [s] = " << profileValues[28] << nl
        << "    migration size exchange [s]  = " << profileValues[29] << nl
        << "    migration request post [s]   = " << profileValues[30] << nl
        << "    migration wait max [s]       = " << profileValues[31] << nl
        << "    migration deserialize [s]    = " << profileValues[32] << nl
        << "    migration post max [s]       = " << profileValues[33] << nl
        << "    migration candidate gather [s]= " << profileValues[34] << nl
        << "    DLB migration wall [s]       = " << profileValues[4] << nl
        << "    output migration wall [s]    = " << profileValues[6] << nl;
    if (allParticleCounts_.size() > 0)
    {
        const label minP = min(allParticleCounts_);
        const label maxP = max(allParticleCounts_);
        Info<< "    particles per rank           = min " << minP
            << " max " << maxP
            << " max/min " << (minP > 0 ? scalar(maxP)/scalar(minP) : 0) << nl;
    }
    Info<< "    rank wall time (evolve, " << evolveTimeSteps_ << " steps):"
        << " min " << min(allTimes) << " max " << max(allTimes)
        << " max/min " << (min(allTimes) > 0 ? max(allTimes)/min(allTimes) : 0) << nl;

    if (triggerWallTime_)
    {
        Info<< "    rank step wall (full iter, " << evolveTimeSteps_
            << " steps): min " << min(allStepWalls)
            << " max " << max(allStepWalls)
            << " max/min "
            << (min(allStepWalls) > 0
              ? max(allStepWalls)/min(allStepWalls) : 0) << nl;
    }
    if (backOffEnabled_)
    {
        Info<< "    DLB back-off skips remaining = " << backOffRemaining_
            << ", futile count = " << futileCount_
            << ", last changed ratio = " << lastExecChangedRatio_ << nl;
    }
    if (rebalanceCount_ > 0)
    {
        const label nCells = mesh_.nCells();
        Info<< "    manual rebalances            = " << rebalanceCount_ << nl
            << "    R_changed (total)            = " << totalCellsChanged_
            << " / " << nCells << " (" << scalar(totalCellsChanged_)/scalar(nCells*nProcs_) << " per rank)" << nl
            << "    R_mig (total)                = " << totalParcelsMigrated_ << nl;
    }
    if (autoDLBEnabled_)
    {
        Info<< "    DLB checks                   = " << autoRebalanceChecks_ << nl
            << "    DLB rebalances               = " << autoRebalanceCount_ << nl
            << "    DLB triggered checks         = "
            << autoRebalanceTriggeredChecks_ << nl
            << "    DLB wall max [s]             = "
            << profileValues[16] << nl
            << "    DLB noncomm max [s]          = "
            << profileValues[26] << nl
            << "    DLB explicit comm max [s]    = "
            << profileValues[22] << nl
            << "    DLB check max [s]            = "
            << profileValues[17] << nl
            << "    DLB ParMETIS/repart max [s]  = "
            << profileValues[18] << nl
            << "    DLB migration total max [s]  = "
            << profileValues[19] << nl
            << "    DLB cellOwner write max [s]  = "
            << profileValues[20] << nl
            << "    DLB post-diagnostic max [s]  = "
            << profileValues[21] << nl;
        if (adaptiveAlpha_)
        {
            Info<< "    DLB adaptive alpha final     = " << dlbAlpha_ << nl
                << "    DLB adaptive alpha last imbalance = "
                << adaptiveAlphaLastImbalance_ << nl
                << "    DLB adaptive alpha last step = "
                << adaptiveAlphaLastStep_ << nl;
        }
    }
    Info<< endl;
}

} // End namespace Foam
