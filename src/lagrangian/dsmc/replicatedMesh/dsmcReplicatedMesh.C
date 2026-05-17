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
#include "OCharStream.H"
#include "ISpanStream.H"
#include "IStringStream.H"
#include "scotchDecomp.H"
#include "scotch.h"
#include "parmetis.h"
#include "PrecisionAdaptor.H"
#include "floatScalar.H"
#include <mpi.h>
#include <chrono>
#include <vector>
#include <cstring>
#include <algorithm>
#include <cmath>

namespace Foam
{

// ============================================================================
// 3D Hilbert curve — state tables
// ============================================================================

// Hilbert sub-cube index for (state, octant) → local position [0..7]
// 12 states × 8 octants (ordered by xyz bits)
static const uint8_t hilbertTable[12][8] = {
    {0,1,3,2,7,6,4,5},  // state  0
    {0,1,3,2,7,6,4,5},  // state  1
    {0,3,1,2,7,4,6,5},  // state  2
    {0,1,5,4,7,6,2,3},  // state  3
    {0,1,5,4,7,6,2,3},  // state  4
    {0,5,1,4,7,2,6,3},  // state  5
    {6,2,3,7,5,1,0,4},  // state  6
    {6,2,3,7,5,1,0,4},  // state  7
    {6,4,5,1,2,3,0,7},  // state  8
    {0,3,7,4,1,2,6,5},  // state  9
    {0,3,7,4,1,2,6,5},  // state 10
    {6,4,5,1,2,3,0,7}   // state 11
};

// Next state for (state, octant)
static const uint8_t hilbertNextState[12][8] = {
    { 2, 1, 3, 0, 0,11, 0, 9},  // state  0
    { 3, 1, 1, 0,10, 2, 9, 9},  // state  1
    { 2, 0, 3, 1, 0, 0,11, 0},  // state  2
    { 5, 4, 4, 1,11, 3, 0, 5},  // state  3
    { 5, 3, 1, 0, 1, 5,11, 0},  // state  4
    { 4, 0, 2, 1, 0, 2, 5, 4},  // state  5
    { 8, 7, 7,10, 9, 0,11, 0},  // state  6
    { 6, 4, 4, 9, 8,10, 7, 7},  // state  7
    { 7, 9, 5, 4, 5, 5, 7, 8},  // state  8
    { 1, 1, 2, 3, 6, 3, 7, 0},  // state  9
    { 1, 0, 6, 3, 9, 6, 3,10},  // state 10
    { 8,11,11, 9, 4, 9,10, 1}   // state 11
};


uint64_t dsmcReplicatedMesh::hilbert3DIndex
(
    uint32_t x, uint32_t y, uint32_t z,
    label bits
)
{
    uint64_t h = 0;
    uint8_t state = 0;

    for (label i = bits - 1; i >= 0; --i)
    {
        const uint8_t xi = (x >> i) & 1;
        const uint8_t yi = (y >> i) & 1;
        const uint8_t zi = (z >> i) & 1;
        const uint8_t octant = (xi << 2) | (yi << 1) | zi;

        h |= static_cast<uint64_t>(hilbertTable[state][octant]) << (3 * i);
        state = hilbertNextState[state][octant];
    }
    return h;
}


// ============================================================================
// Constructor / Destructor
// ============================================================================

dsmcReplicatedMesh::dsmcReplicatedMesh(dsmcCloud& cloud, const fvMesh& mesh)
:
    cloud_(cloud), mesh_(mesh),
    cellOwner_(mesh.nCells(), -1), localMesh_(mesh), myCells_(0),
    localParticleCount_(0), allParticleCounts_(0),
    migrationWallTime_(0.0), migrationCalls_(0),
    cellCost_(0), cellCostSteps_(0),
    cellCollisionCost_(0), cellCollisionCostSteps_(0),
    costMoveWeight_(1.0), costCollisionWeight_(1.0),
    evolveStepTime_(0.0), evolveTimeSteps_(0),
    active_(false), nProcs_(1), myRank_(0),
    migrateInterval_(1), stepCounter_(0),
    rebalanceCount_(0), totalCellsChanged_(0), totalParcelsMigrated_(0),
    // Phase C defaults
    autoDLBEnabled_(false),
    nSuperCells_(0), alpha_(8),
    imbalanceThreshold_(0.15),
    cooldownSteps_(100),
    lastAutoRebalanceStep_(-1000),
    autoRebalanceCount_(0),
    superCellCost_(0),
    superCellOfCell_(0),
    superCellOwner_(0),
    superCellCentres_(0),
    // ParDSMC3D trigger state
    productiveTime_(0.0),
    lastEvolveTime_(0.0),
    lastMigrationTime_(0.0),
    tidl_(0.0),
    tdecps_(0.0),
    w1_(1.0e6), w2_(0.0),  // w1 large → sar negative until idle accumulates
    ndecps_(0),
    nEvalPeriods_(0),
    sar_(0.0),
    asyncMigrationPending_(false),
    asyncRecvSize_(0)
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
        }
        MPI_Bcast(part.data(), nCells, MPI_INT, 0, MPI_COMM_WORLD);
        forAll(cellOwner_, i) cellOwner_[i] = label(part[i]);

        Info<< "Replicated mesh: METIS decomposition complete, "
            << nCells << " cells -> " << nProcs_ << " ranks" << endl;
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
    if (!mpiInit) MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs_);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
    allParticleCounts_.setSize(nProcs_, 0);

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
    migrateInterval_ = 1;
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

    if (autoDLBEnabled_)
    {
        imbalanceThreshold_ = mesh_.time().controlDict()
            .lookupOrDefault<scalar>("replicatedMeshDLBImbalanceThreshold", 0.15);
        cooldownSteps_ = mesh_.time().controlDict()
            .lookupOrDefault<label>("replicatedMeshDLBCooldownSteps", 100);
        if (cooldownSteps_ < 10) cooldownSteps_ = 10;

        // TACF cost weights: move (particle count) vs collision (candidates)
        costMoveWeight_ = mesh_.time().controlDict()
            .lookupOrDefault<scalar>("replicatedMeshCostMoveWeight", 3.7);
        costCollisionWeight_ = mesh_.time().controlDict()
            .lookupOrDefault<scalar>("replicatedMeshCostCollisionWeight", 1.0);

        Info<< "Phase C auto DLB: enabled (ParMETIS AdaptiveRepart)"
            << ", imbalanceThreshold=" << imbalanceThreshold_
            << ", cooldownSteps=" << cooldownSteps_ << endl;
    }

    active_ = true;
    Info<< "Replicated mesh: initialized with " << nProcs_ << " MPI ranks"
        << ", migrate interval " << migrateInterval_ << endl;
}


// ============================================================================
// buildSuperCells — Phase C: Scotch over-decomposition
// ============================================================================

void dsmcReplicatedMesh::buildSuperCells()
{
    const label nCells = mesh_.nCells();

    // Scotch decompose to M = alpha*P domains
    pointField cellCentres(nCells);
    forAll(mesh_.cells(), cellI) cellCentres[cellI] = mesh_.cellCentres()[cellI];

    dictionary decompDict;
    decompDict.set("method", "scotch");
    decompDict.set("numberOfSubdomains", nSuperCells_);
    dictionary scotchCoeffs;
    decompDict.set("scotchCoeffs", scotchCoeffs);
    scotchDecomp decomposer(decompDict);

    superCellOfCell_ = decomposer.decompose(mesh_, cellCentres);

    Info<< "Phase C: superCell decomposition complete, "
        << nCells << " cells -> " << nSuperCells_ << " superCells" << endl;
}


// ============================================================================
// updateSuperCellCosts — Phase C: aggregate cell costs to superCells
// ============================================================================

void dsmcReplicatedMesh::updateSuperCellCosts()
{
    const label nCells = mesh_.nCells();

    if (superCellCentres_.size() != nSuperCells_)
        superCellCentres_.setSize(nSuperCells_, vector::zero);
    if (superCellCost_.size() != nSuperCells_)
        superCellCost_.setSize(nSuperCells_, 0.0);

    superCellCost_ = 0.0;
    superCellCentres_ = vector::zero;
    labelList superCellCellCounts(nSuperCells_, 0);

    // TACF cell costs: weighted combination of move (particle count)
    // and collision (candidate pairs)
    const bool haveTACF = (cellCost_.size() == nCells && cellCostSteps_ > 0);
    const bool haveCollisionCost =
        (cellCollisionCost_.size() == nCells && cellCollisionCostSteps_ > 0);

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        const label sc = superCellOfCell_[cellI];
        superCellCentres_[sc] += mesh_.cellCentres()[cellI];
        ++superCellCellCounts[sc];

        if (haveTACF)
        {
            // move proxy (pre-weighted by costMoveWeight_)
            superCellCost_[sc] += cellCost_[cellI];

            // collision proxy (pre-weighted by costCollisionWeight_)
            if (haveCollisionCost)
                superCellCost_[sc] += cellCollisionCost_[cellI];
        }
    }

    // Compute centroids and fallback costs
    for (label sc = 0; sc < nSuperCells_; ++sc)
    {
        if (superCellCellCounts[sc] > 0)
        {
            superCellCentres_[sc] /= scalar(superCellCellCounts[sc]);
        }
        if (!haveTACF)
        {
            // No TACF yet — use cell count as proxy cost
            superCellCost_[sc] = scalar(superCellCellCounts[sc]);
        }
    }
}


// ============================================================================
// reassignByHilbertSuperCell — Phase C: Hilbert SFC load-balanced assignment
// ============================================================================

label dsmcReplicatedMesh::reassignByHilbertSuperCell()
{
    const label nCells = mesh_.nCells();

    // Update superCell data from current TACF
    updateSuperCellCosts();

    // Compute bounding box of superCell centroids for Hilbert normalization
    vector bbMin = vector(GREAT, GREAT, GREAT);
    vector bbMax = vector(-GREAT, -GREAT, -GREAT);
    forAll(superCellCentres_, sc)
    {
        const vector& c = superCellCentres_[sc];
        bbMin.x() = min(bbMin.x(), c.x());
        bbMin.y() = min(bbMin.y(), c.y());
        bbMin.z() = min(bbMin.z(), c.z());
        bbMax.x() = max(bbMax.x(), c.x());
        bbMax.y() = max(bbMax.y(), c.y());
        bbMax.z() = max(bbMax.z(), c.z());
    }

    const vector bbLen = bbMax - bbMin;
    const scalar eps = 1e-12;

    // Compute Hilbert keys for each superCell
    const label bits = 21; // ~2M resolution per axis, 63-bit key
    const uint32_t maxCoord = (uint32_t(1) << bits) - 1;

    // (hilbertKey, superCellIndex) pairs for sorting
    std::vector<std::pair<uint64_t, label>> sortedSCs;
    sortedSCs.reserve(nSuperCells_);

    for (label sc = 0; sc < nSuperCells_; ++sc)
    {
        const vector& c = superCellCentres_[sc];
        const uint32_t ix = uint32_t(
            (bbLen.x() > eps)
                ? ((c.x() - bbMin.x()) / bbLen.x()) * scalar(maxCoord)
                : 0
        );
        const uint32_t iy = uint32_t(
            (bbLen.y() > eps)
                ? ((c.y() - bbMin.y()) / bbLen.y()) * scalar(maxCoord)
                : 0
        );
        const uint32_t iz = uint32_t(
            (bbLen.z() > eps)
                ? ((c.z() - bbMin.z()) / bbLen.z()) * scalar(maxCoord)
                : 0
        );

        const uint64_t key = hilbert3DIndex(ix, iy, iz, bits);
        sortedSCs.emplace_back(key, sc);
    }

    // Sort superCells by Hilbert key
    std::sort(sortedSCs.begin(), sortedSCs.end());

    // Compute total cost and target per rank
    scalar totalCost = 0.0;
    forAll(superCellCost_, sc) totalCost += superCellCost_[sc];
    const scalar targetPerRank = totalCost / scalar(nProcs_);

    // Walk sorted superCells, assign to ranks by cumulative cost
    labelList newSuperCellOwner(nSuperCells_, -1);
    scalar accumCost = 0.0;
    label rankIdx = 0;

    for (const auto& [key, sc] : sortedSCs)
    {
        newSuperCellOwner[sc] = rankIdx;
        accumCost += superCellCost_[sc];

        // Move to next rank when we've exceeded this rank's target
        if (accumCost >= targetPerRank * scalar(rankIdx + 1) && rankIdx < nProcs_ - 1)
        {
            ++rankIdx;
        }
    }

    // Ensure last rank gets remaining superCells
    for (label sc = 0; sc < nSuperCells_; ++sc)
    {
        if (newSuperCellOwner[sc] < 0) newSuperCellOwner[sc] = nProcs_ - 1;
    }

    // Count changes
    label nSuperCellsChanged = 0;
    for (label sc = 0; sc < nSuperCells_; ++sc)
    {
        if (newSuperCellOwner[sc] != superCellOwner_[sc]) ++nSuperCellsChanged;
    }
    superCellOwner_ = newSuperCellOwner;

    // Derive cellOwner_ from superCellOwner_
    labelList oldCellOwner(cellOwner_);
    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        cellOwner_[cellI] = superCellOwner_[superCellOfCell_[cellI]];
    }

    label nCellsChanged = 0;
    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        if (cellOwner_[cellI] != oldCellOwner[cellI]) ++nCellsChanged;
    }

    rebuildMyCells();

    Info<< "Phase C Hilbert: " << nSuperCellsChanged << " / " << nSuperCells_
        << " superCells changed owner, "
        << nCellsChanged << " / " << nCells << " cells changed ("
        << scalar(nCellsChanged)/scalar(nCells) << ")"
        << " totalCost=" << totalCost << " targetPerRank=" << targetPerRank
        << endl;

    return nSuperCellsChanged;
}


// ============================================================================
// reassignByScotchRemap — Phase C: SCOTCH graph remap for load balancing
// ============================================================================

void dsmcReplicatedMesh::reassignByScotchRemap()
{
    const label nCells = mesh_.nCells();

    // ---- Build mesh adjacency graph (CSR format) --------------------------
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNei = mesh_.faceNeighbour();
    const label nIntFaces = mesh_.nInternalFaces();

    labelList cellDeg(nCells, 0);
    for (label fi = 0; fi < nIntFaces; ++fi)
    {
        ++cellDeg[faceOwner[fi]];
        ++cellDeg[faceNei[fi]];
    }

    labelList xadj(nCells + 1);
    xadj[0] = 0;
    for (label i = 0; i < nCells; ++i)
        xadj[i+1] = xadj[i] + cellDeg[i];

    labelList adjncy(xadj[nCells]);
    labelList off = xadj;
    for (label fi = 0; fi < nIntFaces; ++fi)
    {
        const label own = faceOwner[fi];
        const label nei = faceNei[fi];
        adjncy[off[own]++] = nei;
        adjncy[off[nei]++] = own;
    }

    // ---- Vertex weights from TACF costs ----------------------------------
    const bool haveTACF = (cellCost_.size() == nCells && cellCostSteps_ > 0);
    const bool haveColl =
        (cellCollisionCost_.size() == nCells && cellCollisionCostSteps_ > 0);

    List<scalar> cWeights(nCells, 1.0);
    if (haveTACF)
    {
        forAll(cWeights, i)
        {
            cWeights[i] = cellCost_[i];
            if (haveColl) cWeights[i] += cellCollisionCost_[i];
        }
    }

    // Normalise and convert to SCOTCH_Num integers
    const scalar minW = max(min(cWeights), SMALL);
    scalar rangeScale = 1.0;
    {
        const scalar wSum = sum(cWeights) / minW;
        const scalar upper = scalar(std::numeric_limits<SCOTCH_Num>::max() - 1);
        if (wSum > upper) rangeScale = 0.9 * upper / wSum;
    }

    List<SCOTCH_Num> velotab(nCells);
    forAll(velotab, i)
    {
        velotab[i] = static_cast<SCOTCH_Num>
        (
            ((cWeights[i] / minW - 1.0) * rangeScale) + 1.0
        );
    }

    // ---- Precision adaptors for SCOTCH API --------------------------------
    ConstPrecisionAdaptor<SCOTCH_Num, label, List> xadjParam(xadj);
    ConstPrecisionAdaptor<SCOTCH_Num, label, List> adjncyParam(adjncy);

    // ---- Build graph and remap on rank 0 only (memory-efficient) -----------
    List<SCOTCH_Num> outPart(nCells, 0);

    if (myRank_ == 0)
    {
        Info<< "Phase C SCOTCH: building graph (rank 0, "
            << nCells << " cells, " << adjncy.size() << " edges)" << endl;

        SCOTCH_Graph grafdat;
        SCOTCH_graphInit(&grafdat);

        {
            const int ret = SCOTCH_graphBuild
        (
            &grafdat, 0,
            SCOTCH_Num(nCells),
            xadjParam().cdata(), nullptr,
            velotab.cdata(), nullptr,
            SCOTCH_Num(adjncy.size()),
            adjncyParam().cdata(), nullptr
        );
        if (ret)
        {
            FatalErrorInFunction
                << "SCOTCH_graphBuild failed (" << ret << ")" << nl
                << exit(FatalError);
        }
    }

    // ---- Architecture (complete graph, fully connected) -------------------
    SCOTCH_Arch archdat;
    SCOTCH_archInit(&archdat);
    SCOTCH_archCmplt(&archdat, SCOTCH_Num(nProcs_));

    // ---- Current partition (input) ----------------------------------------
    List<SCOTCH_Num> inPart(nCells);
    forAll(cellOwner_, i) inPart[i] = SCOTCH_Num(cellOwner_[i]);

    // Output partition — separate array; SCOTCH_graphRemap takes
    // const input and non-const output
    List<SCOTCH_Num> outPart(nCells, 0);

    // ---- Strategy ---------------------------------------------------------
    SCOTCH_Strat stradat;
    SCOTCH_stratInit(&stradat);

    // ---- Compute new partition (SCOTCH_graphRemap) -------------------------
    // graphRemap minimises data movement while rebalancing.
    // 3rd param = vmlotab (vertex migration cost), NULL = uniform.
    // Vertex LOADS are already set in graphBuild via velotab.
    Info<< "Phase C SCOTCH: calling graphRemap..." << endl;

    #ifdef FE_NOMASK_ENV
    int oldExcepts = fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif

    {
        const int ret = SCOTCH_graphRemap
        (
            &grafdat,               // graph (velotab already set)
            &archdat,               // target architecture
            nullptr,                // vmlotab: vertex migration cost (NULL=uniform)
            1.05,                   // imbalance ratio (5%)
            inPart.cdata(),         // old partition (input, const)
            &stradat,               // strategy
            outPart.data()          // new partition (output)
        );
        if (ret)
        {
            FatalErrorInFunction
                << "SCOTCH_graphRemap failed (" << ret << ")" << nl
                << exit(FatalError);
        }
    }

    #ifdef FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif

    // ---- Cleanup SCOTCH ---------------------------------------------------
    SCOTCH_graphExit(&grafdat);
    SCOTCH_archExit(&archdat);
    SCOTCH_stratExit(&stradat);

    } // end rank 0 only

    // Broadcast new partition from rank 0 to all ranks
    MPI_Bcast(outPart.data(), nCells, MPI_INT, 0, MPI_COMM_WORLD);

    // ---- Update cellOwner_ ------------------------------------------------
    label nChanged = 0;
    forAll(cellOwner_, i)
    {
        const label newOwner = label(outPart[i]);
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

    Info<< "Phase C SCOTCH remap: " << nChanged
        << " / " << nCells << " cells changed ("
        << scalar(nChanged) / scalar(nCells) << ")" << nl;
}


// ============================================================================
// reassignByParMetisAdaptiveRepart — Phase C: ParMETIS incremental repartition
// ============================================================================

void dsmcReplicatedMesh::reassignByParMetisAdaptiveRepart()
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

    // ---- Vertex weights: dual constraint (move + collision) ------------------
    // ncon=2: constraint 1 = particle count (move proxy),
    //         constraint 2 = collision candidates (collision proxy)
    // Real-time scaling: normalize both to wall-time contribution.

    const scalar moveTime = cloud_.evolveMoveWallTime();
    const scalar collTime = cloud_.evolveCollisionWallTime();

    label totalParticles = 0;
    label totalCandidates = 0;
    const labelList& nCandPerCell = cloud_.nCandidatesPerCell();
    for (label i = 0; i < myN; ++i)
    {
        const label gi = myStart + i;
        if (cloud_.cellOccupancy().size() > gi)
            totalParticles += cloud_.cellOccupancy()[gi].size();
        if (nCandPerCell.size() > gi)
            totalCandidates += nCandPerCell[gi];
    }

    scalar collScale = 1.0;
    if (moveTime > SMALL && collTime > SMALL
        && totalParticles > 0 && totalCandidates > 0)
    {
        const scalar costPerParticle = moveTime / scalar(totalParticles);
        const scalar costPerCandidate = collTime / scalar(totalCandidates);
        collScale = costPerCandidate / costPerParticle;
    }

    idx_t ncon = 2;
    List<idx_t> vwgt(myN * ncon, 1);

    for (label i = 0; i < myN; ++i)
    {
        const label gi = myStart + i;
        label nPart = 0;
        if (cloud_.cellOccupancy().size() > gi)
            nPart = cloud_.cellOccupancy()[gi].size();
        vwgt[i * ncon] = max(idx_t(1), idx_t(nPart));

        idx_t nCand = 1;
        if (nCandPerCell.size() > gi)
            nCand = max(idx_t(1), idx_t(scalar(nCandPerCell[gi]) * collScale));
        vwgt[i * ncon + 1] = nCand;
    }

    // ---- Edge weights: particle density proxy for migration cost ----------
    // Higher weight = more expensive to cut = fewer particles will migrate.
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
            label nOwn = 0, nNei = 0;
            if (cloud_.cellOccupancy().size() > own)
                nOwn = cloud_.cellOccupancy()[own].size();
            if (cloud_.cellOccupancy().size() > nei)
                nNei = cloud_.cellOccupancy()[nei].size();
            const idx_t ew = max(idx_t(1), idx_t(nOwn + nNei));

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
    idx_t wgtflag = 2;  // weights on vertices only
    idx_t numflag = 0;  // C-style numbering
    idx_t nparts = nProcs_;
    List<real_t> tpwgts(ncon * nparts, real_t(1.0) / real_t(nparts));
    List<real_t> ubvec(ncon, real_t(1.05));
    real_t itr = 1000.0;
    idx_t options[4] = {1, 0, 0, 42};  // options[0]=1: use custom, [3]=seed
    idx_t edgecut = 0;
    List<idx_t> vsize(myN, 1);
    for (label i = 0; i < myN; ++i)
    {
        const label gi = myStart + i;
        label nPart = 0;
        if (cloud_.cellOccupancy().size() > gi)
            nPart = cloud_.cellOccupancy()[gi].size();
        vsize[i] = max(idx_t(1), idx_t(nPart));
    }

    MPI_Comm comm = MPI_COMM_WORLD;

    Info<< "Phase C ParMETIS: AdaptiveRepart (" << nCells << " cells, "
        << ncon << " constraints, collScale=" << collScale << ")" << endl;

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
}


// ============================================================================
// autoRebalance — Phase C: automatic DLB trigger + execution
// ============================================================================

void dsmcReplicatedMesh::autoRebalance()
{
    if (!autoDLBEnabled_ || !active_ || nProcs_ < 2) return;

    const label currentStep = stepCounter_;

    // ---- Accumulate productive time (move+collision, excluding migration) ----
    // ParDSMC3D: toper = move + collision (excludes communication cost)
    // evolveStepTime_ / migrationWallTime_ are cumulative → use deltas
    const scalar stepEvolve = evolveStepTime_ - lastEvolveTime_;
    const scalar stepMig = migrationWallTime_ - lastMigrationTime_;
    productiveTime_ += stepEvolve - stepMig;
    lastEvolveTime_ = evolveStepTime_;
    lastMigrationTime_ = migrationWallTime_;
    ++ndecps_;

    // Evaluate every cooldownSteps_ steps (ParDSMC3D evaluates every 10)
    if (ndecps_ < cooldownSteps_) return;

    // Need TACF data for SCOTCH vertex weights
    const label nCells = mesh_.nCells();
    if (cellCost_.size() != nCells || cellCostSteps_ < 2) return;

    // ---- ParDSMC3D-style trend-based trigger ----
    scalarList allProdTimes(nProcs_, 0.0);
    MPI_Allgather(&productiveTime_, 1, MPI_DOUBLE,
                  allProdTimes.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

    const scalar avgT = sum(allProdTimes) / scalar(nProcs_);
    const scalar maxT = max(allProdTimes);
    if (avgT < SMALL) return;

    // Accumulate idle across evaluation periods
    tidl_ += maxT - avgT;
    ++nEvalPeriods_;
    w2_ = (tidl_ + tdecps_) / scalar(nEvalPeriods_);
    sar_ = w2_ - w1_;
    w1_ = w2_;

    // Reset for next evaluation window
    productiveTime_ = 0.0;
    ndecps_ = 0;

    // ParDSMC3D: trigger when idle rate trend is increasing (sar > 0)
    if (sar_ <= 0.0) return;

    // ---- Triggered ----
    Info<< "\nPhase C auto DLB triggered at step " << currentStep
        << ": sar=" << sar_ << " tidl=" << tidl_
        << " nPeriods=" << nEvalPeriods_
        << " prodTimes " << allProdTimes << nl;

    // Time the decomposition for cost accounting
    const auto tDecStart = std::chrono::steady_clock::now();

    // Reassign by ParMETIS AdaptiveRepart
    reassignByParMetisAdaptiveRepart();

    // Migrate particles to new owners
    migrateParticlesByCellOwner();
    updateParticleCounts();

    const auto tDecEnd = std::chrono::steady_clock::now();
    tdecps_ += std::chrono::duration<scalar>(tDecEnd - tDecStart).count();

    lastAutoRebalanceStep_ = currentStep;
    ++autoRebalanceCount_;

    // ParDSMC3D: reset trigger state after repartition
    // w1_ set large so next sar stays negative until idle re-accumulates
    tidl_ = 0.0;
    nEvalPeriods_ = 0;
    w1_ = 1.0e6;
    w2_ = 0.0;
    sar_ = 0.0;

    Info<< "Phase C auto DLB complete: rebalance #" << autoRebalanceCount_
        << nl << endl;
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

    // ---- Phase 1: partition parcels into kept vs migrating ----
    // Use flat array (moveOrderedParcels_) if available, else fall back to linked list.
    PtrList<OCharStream> sendStreams(nProcs_);
    for (label i = 0; i < nProcs_; ++i)
        if (i != myRank_) sendStreams.set(i, new OCharStream(IOstreamOption::BINARY));

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
            sendBufs[i] = sendStreams[i].release();
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
            ISpanStream is(recvBuf.data(), recvSize, IOstreamOption::BINARY);
            while (!is.eof())
            {
                auto* newp = new dsmcParcel(mesh_, is);
                cloud_.addParticle(newp);
                kept.append(newp);
                ++nRecv;
            }
        }
    }
    else
    {
        // Multi-rank: exchange sizes via Alltoall, then Isend/Irecv
        labelList recvSizes(nProcs_, 0);
        MPI_Alltoall(sendSizes.data(), 1, MPI_INT,
                     recvSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

        // Post Irecv for each source
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

        // Deserialize from all sources
        for (label i = 0; i < nProcs_; ++i)
        {
            if (i != myRank_ && recvSizes[i] > 0)
            {
                ISpanStream is(recvBuf.data() + recvOffsets[i], recvSizes[i],
                               IOstreamOption::BINARY);
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

    const auto tEnd = std::chrono::steady_clock::now();

    // ---- Timing -----------------------------------------------------------------
    migrationWallTime_ += std::chrono::duration<scalar>(tEnd - tStart).count();
    ++migrationCalls_;

    totalParcelsMigrated_ += nMigratedOut;

    if (migrationCalls_ <= 3 || rebalanceCount_ > 0 || autoRebalanceCount_ > 0)
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
    const auto tStart = std::chrono::steady_clock::now();

    // Phase 1: partition kept vs migrating
    PtrList<OCharStream> sendStreams(nProcs_);
    for (label i = 0; i < nProcs_; ++i)
        if (i != myRank_) sendStreams.set(i, new OCharStream(IOstreamOption::BINARY));

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
            sendBufs[i] = sendStreams[i].release();
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
    else
    {
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

    if (asyncReqs_.size() > 0)
        MPI_Waitall(asyncReqs_.size(), asyncReqs_.data(), MPI_STATUSES_IGNORE);

    // Deserialize received parcels and batch-append to moveOrdered
    if (asyncRecvSize_ > 0)
    {
        DynamicList<dsmcParcel*> received(asyncRecvSize_ / 100);
        ISpanStream is(asyncRecvBuf_.data(), asyncRecvSize_, IOstreamOption::BINARY);
        while (!is.eof())
        {
            auto* newp = new dsmcParcel(mesh_, is);
            cloud_.addParticle(newp);
            received.append(newp);
        }
        cloud_.appendBatchToMoveOrdered(received);
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
// accumulateCellCosts — Phase A/B TACF
// ============================================================================

void dsmcReplicatedMesh::accumulateCellCosts()
{
    if (!active_) return;
    const label nCells = mesh_.nCells();
    if (cellCost_.size() != nCells)
    {
        cellCost_.setSize(nCells, 0.0);
        cellCollisionCost_.setSize(nCells, 0.0);
        cellCostSteps_ = 0;
        cellCollisionCostSteps_ = 0;
    }

    // Move cost proxy: particle count per cell
    forAllConstIter(Cloud<dsmcParcel>, cloud_, iter)
    {
        const label cellI = iter().cell();
        if (cellI >= 0 && cellI < nCells) cellCost_[cellI] += costMoveWeight_;
    }
    ++cellCostSteps_;

    // Collision cost proxy: collision candidates per cell
    const labelList& nCand = cloud_.nCandidatesPerCell();
    if (nCand.size() == nCells)
    {
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            if (nCand[cellI] > 0)
                cellCollisionCost_[cellI] += costCollisionWeight_ * scalar(nCand[cellI]);
        }
        ++cellCollisionCostSteps_;
    }
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
// report — profiling summary
// ============================================================================

void dsmcReplicatedMesh::report() const
{
    if (!active_) return;
    Info<< nl << "Replicated mesh profiling summary:" << nl
        << "    migration calls             = " << migrationCalls_ << nl
        << "    migration wall time [s]     = " << migrationWallTime_ << nl
        << "    local particles (final)     = " << localParticleCount_ << nl;
    if (allParticleCounts_.size() > 0)
    {
        const label minP = min(allParticleCounts_);
        const label maxP = max(allParticleCounts_);
        const scalar avgP = scalar(sum(allParticleCounts_))/allParticleCounts_.size();
        Info<< "    particles per rank          = min " << minP
            << " max " << maxP << " avg " << avgP
            << " imbalance " << (avgP > 0 ? maxP/avgP : 0) << nl;
    }
    {
        scalarList allTimes(nProcs_, 0.0);
        allTimes[myRank_] = evolveStepTime_;
        MPI_Allreduce(MPI_IN_PLACE, allTimes.data(), nProcs_, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        const scalar avgT = sum(allTimes)/nProcs_;
        Info<< "    rank wall time (evolve, " << evolveTimeSteps_ << " steps):"
            << " min " << min(allTimes) << " max " << max(allTimes)
            << " avg " << avgT
            << " imbalance " << (avgT > 0 ? max(allTimes)/avgT : 0) << nl;
    }
    if (cellCostSteps_ > 0)
    {
        Info<< "    TACF cell cost (" << cellCostSteps_ << " samples):"
            << " min " << min(cellCost_) << " max " << max(cellCost_)
            << " avg " << sum(cellCost_)/cellCost_.size()
            << " (move weight " << costMoveWeight_ << ")" << nl;
    }
    if (cellCollisionCostSteps_ > 0)
    {
        Info<< "    TACF collision cost (" << cellCollisionCostSteps_ << " samples):"
            << " min " << min(cellCollisionCost_) << " max " << max(cellCollisionCost_)
            << " avg " << sum(cellCollisionCost_)/cellCollisionCost_.size()
            << " (coll weight " << costCollisionWeight_ << ")" << nl;
    }
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
        Info<< "    Phase C auto DLB rebalances = " << autoRebalanceCount_ << nl;
    }
    Info<< endl;
}

} // End namespace Foam
