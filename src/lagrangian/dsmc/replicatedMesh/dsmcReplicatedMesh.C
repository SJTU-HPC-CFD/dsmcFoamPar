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
#include "passiveParticleCloud.H"
#include "volFields.H"
#include "scotchDecomp.H"
#include "domainDecomposition.H"
#include "parmetis.h"
#include <mpi.h>
#include <chrono>
#include <cstring>
#include <limits>

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

// ============================================================================
// Constructor / Destructor
// ============================================================================

dsmcReplicatedMesh::dsmcReplicatedMesh(dsmcCloud& cloud, const fvMesh& mesh)
:
    cloud_(cloud), mesh_(mesh),
    cellOwner_(mesh.nCells(), -1), localMesh_(mesh), myCells_(0),
    localParticleCount_(0), allParticleCounts_(0),
    migrationWallTime_(0.0), migrationCalls_(0),
    evolveStepTime_(0.0), evolveTimeSteps_(0),
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
    collDivisor_(128),
    alpha_(1.0),
    asyncMigrationPending_(false),
    asyncRecvSize_(0),
    useNoAlltoall_(false),
    useFlatTransfer_(false),
    writeMode_("gathered")
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
    for (label cellI = 0; cellI < nCells; ++cellI)
        if (cellOwner_[cellI] == myRank_) myCells_.append(cellI);
    Info<< "Replicated mesh: rank " << myRank_ << " owns "
        << myCells_.size() << " / " << nCells << " cells" << endl;
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
    const label defaultK = (nProcs_ >= 8) ? 128 : 64;
    collDivisor_ = mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBInitialK", defaultK);
    alpha_ = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBInitialAlpha", 1.0);

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
        const bool skipFixedKPostDiag = mesh_.time().controlDict()
            .lookupOrDefault<bool>("replicatedMeshDLBSkipFixedKPostDiag", false);
        const label minGapSteps = mesh_.time().controlDict()
            .lookupOrDefault<label>("replicatedMeshDLBMinGapSteps", 0);

        Info<< "Phase C auto DLB: enabled (ParMETIS AdaptiveRepart)"
            << ", imbalanceThreshold=" << imbalanceThreshold_
            << ", dlbSteps=" << dlbSteps_
            << ", sarSteps=" << sarSteps_
            << ", triggerMode=" << autoDLBTriggerMode_
            << ", minGapSteps=" << minGapSteps
            << ", checkCollective=" << checkCollective
            << ", skipFixedKPostDiag=" << skipFixedKPostDiag << endl;
        if (forcedDLBSteps_.size())
        {
            Info<< "Phase C auto DLB: forced trigger steps "
                << forcedDLBSteps_ << endl;
        }
    }

    active_ = true;
    Info<< "Replicated mesh: initialized with " << nProcs_ << " MPI ranks"
        << ", migrate interval " << migrateInterval_ << endl;
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
    forAll(cloud_.cellOccupancy(), cellI)
    {
        if (cellI < nCells)
        {
            localCellParticles[cellI] =
                idx_t(cloud_.cellOccupancy()[cellI].size());
        }
    }
    MPI_Allreduce
    (
        localCellParticles.data(),
        globalCellParticles.data(),
        nCells,
        MPI_INT,
        MPI_SUM,
        MPI_COMM_WORLD
    );

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
        forAll(globalCellParticles, cellI)
        {
            const idx_t nPart = globalCellParticles[cellI];
            if (nPart > 0)
            {
                ++activeCells;
                totalParticles += label(nPart);
                maxCellParticles = max(maxCellParticles, nPart);
            }
        }
        Info<< "Phase C ParMETIS weights: particles=" << totalParticles
            << " activeCells=" << activeCells
            << " maxCellParticles=" << maxCellParticles << endl;
    }

    // ---- Vertex weights: dual constraint (move + collision) ------------------
    // ncon=2: constraint 0 = N^alpha for move balance
    //         constraint 1 = N*(N-1) for collision balance

    const scalar moveTime = cloud_.evolveMoveWallTime();
    const scalar collTime = cloud_.evolveCollisionWallTime();

    // Dual-constraint DLB: ncon=2
    // Constraint 0: compressed particle count (balance move)
    // Constraint 1: N*(N-1) collision proxy
    const bool useDualConstraint = mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshDLBDualConstraint", false);
    idx_t ncon = useDualConstraint ? 2 : 1;
    idx_t wgtflag = 2;  // vertex weights only

    List<idx_t> vwgt(myN * ncon, 1);
    for (label i = 0; i < myN; ++i)
    {
        const label gi = myStart + i;
        const idx_t nPart = globalCellParticles[gi];

        if (useDualConstraint)
        {
            // Constraint 0: N^alpha (move balance)
            const scalar wMove = (nPart > 1)
                ? std::pow(scalar(nPart), alpha_) : scalar(nPart);
            vwgt[i*2 + 0] = positiveWeight(wMove);
            // Constraint 1: N*(N-1) (collision balance, compressed range)
            const scalar wColl =
                scalar(nPart)*max(scalar(nPart - 1), scalar(0)) + scalar(1);
            vwgt[i*2 + 1] = positiveWeight(wColl);
        }
        else
        {
            const scalar w = (nPart > 1)
                ? std::pow(scalar(nPart), alpha_) : scalar(nPart);
            vwgt[i] = positiveWeight(w);
        }
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
    for (label i = 0; i < myN; ++i)
        part[i] = cellOwner_[myStart + i];

    // ---- ParMETIS parameters -----------------------------------------------
    idx_t numflag = 0;  // C-style numbering
    idx_t nparts = nProcs_;
    List<real_t> tpwgts(ncon * nparts, real_t(1.0) / real_t(nparts));
    const real_t ubvecVal = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBUbvec", 1.05);
    List<real_t> ubvec(ncon, ubvecVal);
    if (useDualConstraint && ncon == 2)
    {
        ubvec[1] = mesh_.time().controlDict().lookupOrDefault<scalar>
            ("replicatedMeshDLBUbvec1", 1.5);
    }

    Info<< "Phase C ParMETIS: alpha=" << alpha_ << " ncon=" << ncon
        << " ubvec=[";
    for (label c = 0; c < ncon; ++c)
        Info<< (c > 0 ? "," : "") << ubvec[c];
    Info<< "] moveRatio=" << (moveTime / max(moveTime + collTime, SMALL))
        << " (moveT=" << moveTime << "s, collT=" << collTime << "s)" << endl;

    real_t itr = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBItr", 100.0);
    idx_t options[4] = {1, 0, 0, 42};  // options[0]=1: use custom, [3]=seed
    idx_t edgecut = 0;
    // vsize: 2-rank uses N (conservative), 4+ rank uses N^vsExp (vsExp=0 means free migration)
    const scalar vsExp = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBVsizeExp", 0.5);
    List<idx_t> vsize(myN, 1);
    if (nProcs_ <= 2)
    {
        for (label i = 0; i < myN; ++i)
        {
            const label gi = myStart + i;
            const idx_t nPart = globalCellParticles[gi];
            vsize[i] = positiveWeight(scalar(nPart));
        }
    }
    else if (vsExp > SMALL)
    {
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

    ParMETIS_V3_AdaptiveRepart
    (
        vtxdist.data(), xadj.data(), adjncy.data(),
        vwgt.data(), vsize.data(), nullptr,
        &wgtflag, &numflag, &ncon, &nparts,
        tpwgts.data(), ubvec.data(), &itr,
        options, &edgecut, part.data(), &comm
    );

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
    MPI_Allgatherv(part.data(), myN, MPI_INT,
                   fullPart.data(), recvCounts.data(), displs.data(),
                   MPI_INT, MPI_COMM_WORLD);

    // ---- Update cellOwner_ ------------------------------------------------
    label nChanged = 0;
    forAll(cellOwner_, i)
    {
        const label newOwner = label(fullPart[i]);
        if (newOwner != cellOwner_[i])
        {
            cellOwner_[i] = newOwner;
            ++nChanged;
        }
    }

    // ---- Rebuild local structures -----------------------------------------
    rebuildMyCells();
    localMesh_.build(myCells_);

    totalCellsChanged_ += nChanged;

    Info<< "Phase C ParMETIS AdaptiveRepart: " << nChanged
        << " / " << nCells << " cells changed ("
        << scalar(nChanged) / scalar(nCells) << ")" << nl;

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
        };

    // Post-DLB snapshot: per-rank load for 5 steps after DLB
    if (dlbProfile && postDLBSnapshotCountdown_ > 0)
    {
        const scalar myMoveT = cloud_.evolveMoveWallTime();
        const scalar myCollT = cloud_.evolveCollisionWallTime();
        scalarList snapMove(nProcs_, 0.0), snapColl(nProcs_, 0.0);
        MPI_Allgather(&myMoveT, 1, MPI_DOUBLE,
                      snapMove.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(&myCollT, 1, MPI_DOUBLE,
                      snapColl.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
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

    // ---- Accumulate productive time (move+collision, excluding migration) ----
    const scalar stepEvolve = evolveStepTime_ - lastEvolveTime_;
    const scalar stepMig = migrationWallTime_ - lastMigrationTime_;
    const scalar stepProductive = stepEvolve - stepMig;
    productiveTime_ += stepProductive;
    lastEvolveTime_ = evolveStepTime_;
    lastMigrationTime_ = migrationWallTime_;
    ++ndecps_;
    ++sarEvalSteps_;

    // Minimum step guard (ParMETIS needs stable particle distribution)
    if (currentStep < 30) return;

    bool triggered = false;
    bool forcedTriggered = false;

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

    if (forcedDLBSteps_.size())
    {
        if (sarEvalSteps_ >= sarSteps_)
        {
            MPI_Barrier(MPI_COMM_WORLD);
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
            const scalar localT = productiveTime_;
            scalar maxT = 0.0;
            scalar minT = GREAT;
            gatherLoadExtrema(localT, maxT, minT);

            if (maxT > SMALL)
            {
                const scalar loadImbalance = maxT / max(minT, SMALL);

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
                    (ndecps_ >= dlbSteps_) && (loadImbalance > imbalanceThreshold_);

                triggered = minGapSatisfied && (sarTriggered || thresholdTriggered);
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
                        << " threshold=" << imbalanceThreshold_ << nl;
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

        const scalar localT = productiveTime_;
        scalar maxT = 0.0;
        scalar minT = GREAT;
        gatherLoadExtrema(localT, maxT, minT);

        if (maxT > SMALL)
        {
            const scalar loadImbalance = maxT / max(minT, SMALL);

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

            triggered = minGapSatisfied && (sarTriggered || thresholdTriggered);
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
                    << " threshold=" << imbalanceThreshold_ << nl;
            }
        }

        // Current SAR sampling window has been consumed, but productiveTime_
        // remains cumulative until an actual DLB trigger completes.
        sarEvalSteps_ = 0;
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
    const auto tRepartEnd = std::chrono::steady_clock::now();
    autoRebalanceRepartWallTime_ +=
        std::chrono::duration<scalar>(tRepartEnd - tDecStart).count();

    // Skip migration if partition unchanged
    if (nChanged > 0)
    {
        const auto tMigrate0 = std::chrono::steady_clock::now();
        migrateParticlesByCellOwner();
        updateParticleCounts();
        const auto tMigrate1 = std::chrono::steady_clock::now();
        autoRebalanceMigrationWallTime_ +=
            std::chrono::duration<scalar>(tMigrate1 - tMigrate0).count();
    }

    const auto tDecEnd = std::chrono::steady_clock::now();
    tdecps_ += std::chrono::duration<scalar>(tDecEnd - tDecStart).count();

    lastAutoRebalanceStep_ = currentStep;
    ++autoRebalanceCount_;

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

    Info<< "Phase C auto DLB complete: rebalance #" << autoRebalanceCount_
        << nl << endl;

    // Inter-DLB load summary + arm post-DLB snapshot
    const auto tPostDiag0 = std::chrono::steady_clock::now();
    {
        const label fixedK = mesh_.time().controlDict().lookupOrDefault<label>
            ("replicatedMeshDLBFixedK", 0);
        const bool skipFixedKPostDiag = mesh_.time().controlDict()
            .lookupOrDefault<bool>("replicatedMeshDLBSkipFixedKPostDiag", false);
        const bool needPostDiag =
            dlbProfile || (fixedK <= 0) || !skipFixedKPostDiag;
        scalarList allMoveT, allCollT;
        if (needPostDiag)
        {
            const scalar myMoveT = cloud_.evolveMoveWallTime();
            const scalar myCollT = cloud_.evolveCollisionWallTime();
            allMoveT.setSize(nProcs_, 0.0);
            allCollT.setSize(nProcs_, 0.0);
            MPI_Allgather(&myMoveT, 1, MPI_DOUBLE,
                          allMoveT.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgather(&myCollT, 1, MPI_DOUBLE,
                          allCollT.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
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

        // ---- Adaptive K PID: use inter-DLB delta for error signal ----
        const label adaptiveKMode = mesh_.time().controlDict().lookupOrDefault<label>
            ("replicatedMeshDLBAdaptiveKMode", 1);
        if (fixedK <= 0 && lastDLBMoveT_.size() == nProcs_)
        {
            // Compute bottleneck + imbalance from delta (max/min)
            scalar maxWork = 0.0, minWork = GREAT;
            label bottleneckRank = 0;
            scalarList dMove(nProcs_), dColl(nProcs_);
            for (label i = 0; i < nProcs_; ++i)
            {
                dMove[i] = allMoveT[i] - lastDLBMoveT_[i];
                dColl[i] = allCollT[i] - lastDLBCollT_[i];
                const scalar w = dMove[i] + dColl[i];
                if (w > maxWork) { maxWork = w; bottleneckRank = i; }
                if (w < minWork) { minWork = w; }
            }
            const scalar workImbalance = maxWork / max(minWork, SMALL);

            // Direction from moveRatio
            const scalar bnMoveRatio = dMove[bottleneckRank]
                / max(dMove[bottleneckRank] + dColl[bottleneckRank], SMALL);
            scalar avgMoveRatio = 0.0;
            for (label i = 0; i < nProcs_; ++i)
                avgMoveRatio += dMove[i] / max(dMove[i] + dColl[i], SMALL);
            avgMoveRatio /= scalar(nProcs_);

            const scalar totalDelta = sum(dMove) + sum(dColl);
            const label oldK = collDivisor_;

            if (adaptiveKMode == 0)
            {
                // Adaptive alpha with hill-climbing:
                // weight = N^alpha, alpha in [0.5, 2.0]
                // Adjust alpha to minimize workImbalance.
                // Direction hint from move/coll excess, but reverse if
                // imbalance worsened after last adjustment.

                scalar sumMove = 0.0, sumColl = 0.0;
                for (label i = 0; i < nProcs_; ++i)
                {
                    sumMove += dMove[i];
                    sumColl += dColl[i];
                }
                const scalar avgMove = sumMove / nProcs_;
                const scalar avgColl = sumColl / nProcs_;
                const scalar bnMoveExcess = dMove[bottleneckRank] - avgMove;
                const scalar bnCollExcess = dColl[bottleneckRank] - avgColl;

                static scalar lastImbalance = GREAT;
                static scalar lastStep = 0.0;

                if (totalDelta > 2.0 * nProcs_)
                {
                    const scalar oldAlpha = alpha_;

                    // Determine direction from move/coll excess
                    // move excess > coll excess → alpha should decrease
                    // coll excess > move excess → alpha should increase
                    const scalar totalExcess = mag(bnMoveExcess) + mag(bnCollExcess) + SMALL;
                    scalar direction = (bnCollExcess - bnMoveExcess) / totalExcess;

                    // If imbalance worsened since last adjustment, reverse direction
                    if (workImbalance > lastImbalance + 0.05 && mag(lastStep) > SMALL)
                    {
                        direction = -Foam::sign(lastStep);
                    }

                    // Step size proportional to imbalance severity, but small
                    const scalar severity = min(workImbalance - 1.0, scalar(3.0));
                    const scalar step = direction * min(severity, scalar(2.0)) * 0.04;

                    alpha_ += step;
                    alpha_ = max(scalar(0.5), min(scalar(2.0), alpha_));

                    lastStep = alpha_ - oldAlpha;
                    lastImbalance = workImbalance;

                    Info<< "Phase C adaptive alpha: bottleneck=rank"
                        << bottleneckRank
                        << " imbal=" << workImbalance
                        << " bnMoveExc=" << bnMoveExcess
                        << " bnCollExc=" << bnCollExcess
                        << " dir=" << direction
                        << " alpha: " << oldAlpha << " -> " << alpha_ << nl;
                }
            }
            else
            {
                // Hill-Climbing controller (mode=1, default)
                static scalar lastImbalance = GREAT;
                static label lastDirection = -1;

                if (workImbalance > 1.3 && totalDelta > 2.0 * nProcs_)
                {
                    label direction = (bnMoveRatio > avgMoveRatio + 0.03) ? 1 : -1;

                    if (workImbalance > lastImbalance + 0.05)
                        direction = -lastDirection;

                    const scalar step = direction * min(scalar(0.4),
                        (workImbalance - 1.0) * scalar(0.4));

                    collDivisor_ = label(
                        scalar(collDivisor_) * (scalar(1.0) + step) + scalar(0.5)
                    );
                    collDivisor_ = max(label(1), min(label(2048), collDivisor_));
                    lastDirection = direction;
                }
                lastImbalance = workImbalance;

                Info<< "Phase C adaptive K (Hill): bottleneck=rank" << bottleneckRank
                    << " max/min=" << workImbalance
                    << " (moveRatio=" << bnMoveRatio << " vs avg=" << avgMoveRatio << ")"
                    << " K: " << oldK << " -> " << collDivisor_ << nl;
            }
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
    const word procMeshInstance = writeTimeMesh ? timeName : runTime.constant();

    if (myRank_ == 0)
    {
        Info<< "Replicated mesh: processor output begin at time "
            << timeName << " (read mesh instance " << meshReadInstance
            << ", processor mesh instance " << procMeshInstance << ")"
            << endl;
    }

    Info<< "Replicated mesh: constructing processor mesh/addressing"
        << endl;

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
    decomposition.setProcessorMeshInstance(procMeshInstance);
    decomposition.setProcessorMeshWriteProc(myRank_);
    decomposition.decomposeMesh();
    decomposition.writeDecomposition(false);

    Info<< "Replicated mesh: rank " << myRank_
        << " wrote processor" << myRank_ << " mesh/addressing"
        << endl;

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

    Info<< "Replicated mesh: rank " << myRank_
        << " loaded cellProcAddressing with "
        << cellProcAddressing.size() << " cells" << endl;

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

    passiveParticleCloud positions
    (
        procMesh,
        cloud_.name(),
        IDLList<passiveParticle>()
    );

    IOField<vector> U(positions.fieldIOobject("U", IOobject::NO_READ), nParcels);
    IOField<scalar> RWF
    (
        positions.fieldIOobject("radialWeight", IOobject::NO_READ),
        nParcels
    );
    IOField<scalar> ERot
    (
        positions.fieldIOobject("ERot", IOobject::NO_READ),
        nParcels
    );
    IOField<labelField> vibLevel
    (
        positions.fieldIOobject("vibLevel", IOobject::NO_READ),
        nParcels
    );
    IOField<label> ELevel
    (
        positions.fieldIOobject("ELevel", IOobject::NO_READ),
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
    IOField<label> stuckToWall
    (
        positions.fieldIOobject("stuckToWall", IOobject::NO_READ),
        nParcels
    );
    IOField<scalarField> wallTemperature
    (
        positions.fieldIOobject("wallTemperature", IOobject::NO_READ),
        nParcels
    );
    IOField<vectorField> wallVectors
    (
        positions.fieldIOobject("wallVectors", IOobject::NO_READ),
        nParcels
    );
    IOField<label> isTracked
    (
        positions.fieldIOobject("isTracked", IOobject::NO_READ),
        nParcels
    );
    IOField<label> inPatchId
    (
        positions.fieldIOobject("inPatchId", IOobject::NO_READ),
        nParcels
    );
    IOField<scalar> tracerInitialTime
    (
        positions.fieldIOobject("tracerInitialTime", IOobject::NO_READ),
        nParcels
    );
    IOField<vector> tracerInitialPosition
    (
        positions.fieldIOobject("tracerInitialPosition", IOobject::NO_READ),
        nParcels
    );
    IOField<vector> tracerCurrentPosition
    (
        positions.fieldIOobject("tracerCurrentPosition", IOobject::NO_READ),
        nParcels
    );
    IOField<vector> tracerDistanceTravelled
    (
        positions.fieldIOobject("tracerDistanceTravelled", IOobject::NO_READ),
        nParcels
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
        else
        {
            wallTemperature[i] = scalarField(4, 0.0);
            wallVectors[i] = vectorField(4, vector::zero);
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
    label myParticleCount = nParcels;
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

    if (nParcels)
    {
        Info<< "Replicated mesh: rank " << myRank_
            << " writing processor" << myRank_ << " lagrangian cloud"
            << endl;
        IOPosition<Cloud<passiveParticle>>(positions).write();
        U.write();
        if (hasRWF) RWF.write();
        if (hasERot) ERot.write();
        if (hasELevel) ELevel.write();
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
    }

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

        if (cloud_.hasMoveOrderedParcels())
        {
            const auto& ordered = cloud_.moveOrderedParcels();
            for (label i = 0; i < ordered.size(); ++i)
            {
                dsmcParcel& p = *ordered[i];
                const label dstRank = cellOwner_[p.cell()];
                if (dstRank != myRank_)
                {
                    sendTD[dstRank].append(dsmcParcel::TransferData());
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
        const auto t1 = std::chrono::steady_clock::now();

        forAll(toDelete, i) { cloud_.deleteParcel(toDelete[i]); }

        labelList sendSizes(nProcs_, 0);
        for (label i = 0; i < nProcs_; ++i)
            sendSizes[i] = sendTD[i].size() * sizeof(dsmcParcel::TransferData);
        const auto t2 = std::chrono::steady_clock::now();

        // MPI exchange (Alltoall sizes then Isend/Irecv data)
        labelList recvSizes(nProcs_, 0);
        MPI_Alltoall(sendSizes.data(), 1, MPI_INT,
                     recvSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

        label totalRecvBytes = 0;
        for (label i = 0; i < nProcs_; ++i)
            if (i != myRank_) totalRecvBytes += recvSizes[i];

        List<char> recvBuf(totalRecvBytes);
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
                MPI_Isend(sendTD[i].data(), sendSizes[i],
                          MPI_BYTE, i, 0, MPI_COMM_WORLD, &req);
                allReqs.append(req);
            }
        }
        if (allReqs.size() > 0)
            MPI_Waitall(allReqs.size(), allReqs.data(), MPI_STATUSES_IGNORE);

        // Deserialize from flat buffer
        label nRecv = 0;
        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && recvSizes[i] > 0)
            {
                const label nParcels =
                    recvSizes[i] / sizeof(dsmcParcel::TransferData);
                const auto* tdArr = reinterpret_cast<const dsmcParcel::TransferData*>
                    (recvBuf.data() + recvOffsets[i]);
                for (label j = 0; j < nParcels; ++j)
                {
                    auto* newp = dsmcParcel::unpackTransfer(mesh_, tdArr[j]);
                    cloud_.addParticle(newp);
                    kept.append(newp);
                    ++nRecv;
                }
            }
        }
        const auto t3 = std::chrono::steady_clock::now();

        cloud_.setMoveOrderedParcels(kept);

        // Exchange candidate counts
        {
            const labelList& nCandPerCell = cloud_.nCandidatesPerCell();
            label localCands = 0;
            for (label i = 0; i < nCandPerCell.size(); ++i)
                localCands += nCandPerCell[i];
            allProcCandidates_.setSize(nProcs_, 0);
            MPI_Allgather(&localCands, 1, MPI_INT,
                          allProcCandidates_.data(), 1, MPI_INT, MPI_COMM_WORLD);
        }

        const auto tEnd = std::chrono::steady_clock::now();
        migrationWallTime_ += std::chrono::duration<scalar>(tEnd - tStart).count();
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
                << " [pack " << std::chrono::duration<scalar>(t2 - t1).count() << "s"
                << " MPI "  << std::chrono::duration<scalar>(t3 - t2).count() << "s"
                << " deser "<< std::chrono::duration<scalar>(tEnd - t3).count() << "s]"
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

    if (cloud_.hasMoveOrderedParcels())
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
    cloud_.setMoveOrderedParcels(kept);

    // ---- Phase 5: exchange per-rank candidate counts for offload planner ----
    {
        const labelList& nCandPerCell = cloud_.nCandidatesPerCell();
        label localCands = 0;
        for (label i = 0; i < nCandPerCell.size(); ++i)
        {
            localCands += nCandPerCell[i];
        }
        allProcCandidates_.setSize(nProcs_, 0);
        MPI_Allgather(&localCands, 1, MPI_INT,
                      allProcCandidates_.data(), 1, MPI_INT, MPI_COMM_WORLD);
    }

    const auto tEnd = std::chrono::steady_clock::now();

    // ---- Timing -----------------------------------------------------------------
    migrationWallTime_ += std::chrono::duration<scalar>(tEnd - tStart).count();
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

    if (cloud_.hasMoveOrderedParcels())
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

    // Set moveOrdered with kept parcels (local only, no received yet)
    cloud_.setMoveOrderedParcels(kept);

    asyncMigrationPending_ = true;
    totalParcelsMigrated_ += nMigratedOut;
    ++migrationCalls_;

    migrationWallTime_ += std::chrono::duration<scalar>(
        std::chrono::steady_clock::now() - tStart).count();
}


// ============================================================================
// migrateFinish — async Phase A: Waitall + deserialize received parcels
// ============================================================================

void dsmcReplicatedMesh::migrateFinish()
{
    if (!asyncMigrationPending_) return;

    if (useNoAlltoall_ && asyncReqs_.size() > 0)
    {
        // No-Alltoall path: use MPI_Get_count per peer
        List<MPI_Status> statuses(asyncReqs_.size());
        MPI_Waitall(asyncReqs_.size(), asyncReqs_.data(), statuses.data());

        // First nProcs_-1 requests are Irecv (one per peer != myRank_)
        DynamicList<dsmcParcel*> received(1024);
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
        if (received.size() > 0)
        {
            cloud_.appendBatchToMoveOrdered(received);
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
            cloud_.appendBatchToMoveOrdered(received);
        }
    }

    asyncReqs_.clear();
    asyncSendBuf_.clear();
    asyncRecvBuf_.clear();
    asyncMigrationPending_ = false;
}


// ============================================================================
// updateParticleCounts — Phase A
// ============================================================================

void dsmcReplicatedMesh::updateParticleCounts()
{
    localParticleCount_ = cloud_.size();
    allParticleCounts_.setSize(nProcs_);
    MPI_Allgather(&localParticleCount_, 1, MPI_INT,
                  allParticleCounts_.data(), 1, MPI_INT, MPI_COMM_WORLD);
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

            const label sendCount = sendTD.size();
            MPI_Send(&sendCount, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
            if (sendCount > 0)
            {
                MPI_Send
                (
                    sendTD.data(),
                    sendCount*sizeof(dsmcParcel::TransferData),
                    MPI_BYTE,
                    0,
                    11,
                    MPI_COMM_WORLD
                );
            }

            cloud_.clear();
        }
        else
        {
            for (label src = 1; src < nProcs_; ++src)
            {
                label recvCount = 0;
                MPI_Recv
                (
                    &recvCount,
                    1,
                    MPI_INT,
                    src,
                    10,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );

                if (recvCount > 0)
                {
                    List<dsmcParcel::TransferData> recvTD(recvCount);
                    MPI_Recv
                    (
                        recvTD.data(),
                        recvCount*sizeof(dsmcParcel::TransferData),
                        MPI_BYTE,
                        src,
                        11,
                        MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE
                    );

                    forAll(recvTD, i)
                    {
                        cloud_.addParticle
                        (
                            dsmcParcel::unpackTransfer(mesh_, recvTD[i])
                        );
                    }
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
// report — profiling summary
// ============================================================================

void dsmcReplicatedMesh::report() const
{
    if (!active_) return;

    scalarList allTimes(nProcs_, 0.0);
    allTimes[myRank_] = evolveStepTime_;
    MPI_Allreduce(MPI_IN_PLACE, allTimes.data(), nProcs_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    scalar autoRebalanceWallMax = autoRebalanceWallTime_;
    scalar autoRebalanceCheckWallMax = autoRebalanceCheckWallTime_;
    scalar autoRebalanceRepartWallMax = autoRebalanceRepartWallTime_;
    scalar autoRebalanceMigrationWallMax = autoRebalanceMigrationWallTime_;
    scalar autoRebalanceWriteWallMax = autoRebalanceWriteWallTime_;
    scalar autoRebalancePostDiagWallMax = autoRebalancePostDiagWallTime_;

    MPI_Allreduce
    (
        MPI_IN_PLACE,
        &autoRebalanceWallMax,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        &autoRebalanceCheckWallMax,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        &autoRebalanceRepartWallMax,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        &autoRebalanceMigrationWallMax,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        &autoRebalanceWriteWallMax,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );
    MPI_Allreduce
    (
        MPI_IN_PLACE,
        &autoRebalancePostDiagWallMax,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        MPI_COMM_WORLD
    );

    if (myRank_ != 0) return;

    Info<< nl << "Replicated mesh profiling summary:" << nl
        << "    migration calls             = " << migrationCalls_ << nl
        << "    migration wall time [s]     = " << migrationWallTime_ << nl
        << "    local particles (final)     = " << localParticleCount_ << nl;
    if (allParticleCounts_.size() > 0)
    {
        const label minP = min(allParticleCounts_);
        const label maxP = max(allParticleCounts_);
        Info<< "    particles per rank          = min " << minP
            << " max " << maxP
            << " max/min " << (minP > 0 ? scalar(maxP)/scalar(minP) : 0) << nl;
    }
    Info<< "    rank wall time (evolve, " << evolveTimeSteps_ << " steps):"
        << " min " << min(allTimes) << " max " << max(allTimes)
        << " max/min " << (min(allTimes) > 0 ? max(allTimes)/min(allTimes) : 0) << nl;
    if (rebalanceCount_ > 0)
    {
        const label nCells = mesh_.nCells();
        Info<< "    Phase B rebalances          = " << rebalanceCount_ << nl
            << "    R_changed (total)           = " << totalCellsChanged_
            << " / " << nCells << " (" << scalar(totalCellsChanged_)/scalar(nCells*nProcs_) << " per rank)" << nl
            << "    R_mig (total)              = " << totalParcelsMigrated_ << nl;
    }
    if (autoDLBEnabled_)
    {
        const scalar autoRebalanceAccounted =
            autoRebalanceCheckWallMax
          + autoRebalanceRepartWallMax
          + autoRebalanceMigrationWallMax
          + autoRebalanceWriteWallMax
          + autoRebalancePostDiagWallMax;

        Info<< "    Phase C auto DLB checks      = " << autoRebalanceChecks_ << nl
            << "    Phase C auto DLB rebalances = " << autoRebalanceCount_ << nl
            << "    Phase C auto DLB triggered checks = "
            << autoRebalanceTriggeredChecks_ << nl
            << "    Phase C auto DLB wall max [s]= "
            << autoRebalanceWallMax << nl
            << "    Phase C auto DLB check max [s]= "
            << autoRebalanceCheckWallMax << nl
            << "    Phase C ParMETIS max [s]     = "
            << autoRebalanceRepartWallMax << nl
            << "    Phase C migration max [s]    = "
            << autoRebalanceMigrationWallMax << nl
            << "    Phase C cellOwner write max [s] = "
            << autoRebalanceWriteWallMax << nl
            << "    Phase C post-diagnostic max [s] = "
            << autoRebalancePostDiagWallMax << nl
            << "    Phase C auto DLB accounted max [s] = "
            << autoRebalanceAccounted << nl
            << "    Phase C auto DLB residual max [s] = "
            << autoRebalanceWallMax - autoRebalanceAccounted << nl;
    }
    Info<< endl;
}

} // End namespace Foam
