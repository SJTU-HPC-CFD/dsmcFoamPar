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
#include "IFstream.H"
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
#include <cstring>
#include <algorithm>
#include <fstream>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>

namespace Foam
{

// ---- Flat-path async migration overlap (worklog comm_further_opt.md §2) ----
// State carried across the begin -> (gap: buildCellOccupancy over kept) ->
// finish split.  File-static by design (§12 red line: zero header/layout
// change).  One replicated mesh per process.
struct dsmcFlatAsyncState
{
    bool active;
    bool gapBuildDone;   // caller ran buildCellOccupancy() between begin/finish
    List<DynamicList<dsmcParcel::TransferData>> sendTD;
    List<std::uint64_t> recvSizes;
    DynamicList<MPI_Request> sendReqs;
    DynamicList<dsmcParcel*> kept;
    label nMigratedOut;
    std::chrono::steady_clock::time_point tStart;
    std::chrono::steady_clock::time_point tPackEnd;
    std::chrono::steady_clock::time_point tLocalPrepEnd;
    std::chrono::steady_clock::time_point tSizeExEnd;
    std::chrono::steady_clock::time_point tReqPostEnd;

    dsmcFlatAsyncState()
    :
        active(false),
        gapBuildDone(false),
        nMigratedOut(0)
    {}
};

dsmcFlatAsyncState dsmcFlatAsync_;

// §15.5: EMA-smoothed per-cell particle counts for DLB weights.
// Replaces the instantaneous globalCellParticles with a time-averaged
// field to avoid chasing transient fluctuations (worklog comm_further_opt.md
// §8: the density field evolves faster than DLB can correct — smoothing
// gives ParMETIS a stable weight field).  File-static: zero header change.
List<scalar> dsmcSmoothCellParticles_;
bool dsmcSmoothWeightEnabled_ = false;
bool dsmcSmoothWeightRead_ = false;
scalar dsmcSmoothWeightTheta_ = 0.3;

// §15.4: particles + cumulative label counters at DLB snapshot time
// (baseline for the post-DLB window deltas).  File-static: zero header
// change.
List<scalar> dsmcLastDLBLabels_;   // double-backed: int32 overflows at 1200w scale

// §19.1: per-rank cumulative step wall / collective wait, for the
// Inter-DLB actual-load table redefined as (wall - comm wait) — the
// same pairing as "Phase C per-rank step work" (§17.5).  Accumulated
// per step in the autoRebalance accumulate block; baselines captured
// at each DLB event.
scalar dsmcStepWallCum_ = 0.0;
scalar dsmcStepWaitCum_ = 0.0;
List<scalar> dsmcLastDLBWallCum_;
List<scalar> dsmcLastDLBWaitCum_;

// §21: proposal-acceptance state (runtime candidate-rate calibration;
// same per-event re-derivation discipline as the reverted §18 — never
// cached across events).
scalar dsmcAccCandLastTotal_ = 0.0;
label dsmcAccCandLastTimeIdx_ = 0;

// §31: local rebuild wall of the last accepted rebalance (rebuildMyCells +
// localMesh_.build), for the acceptance payback estimate.  Per-rank; the
// estimate is reduced globally before the payback comparison so every
// rank takes the same accept/reject branch.
scalar dsmcRebuildLastWall_ = 0.0;

// §34.5: pending particle-count gather posted at the migration tail and
// completed by updateParticleCounts().  File-static: zero header change.
MPI_Request dsmcCountsReq_ = MPI_REQUEST_NULL;
bool dsmcCountsPosted_ = false;

// §35 P0: same-rank regular-migration stage accumulators.  Each is
// accumulated PER RANK and gathered into the report's per-rank table so
// the stages are addable on one rank (no cross-rank phase-max mixing).
double dsmcMigClassifyWall_ = 0.0;   // ordered view + classification + offsets
double dsmcMigPackWall_ = 0.0;       // packTransfer into sendTD (+ kept/toDelete fill)
double dsmcMigRebuildWall_ = 0.0;    // moveOrdered rebuild / clear
double dsmcMigUpcPostWall_ = 0.0;    // MPI_Iallgather post
double dsmcMigUpcWaitWall_ = 0.0;    // MPI_Wait for the gather
double dsmcMigSendBytes_ = 0.0;      // sum over steps
double dsmcMigRecvBytes_ = 0.0;
double dsmcMigPeerCount_ = 0.0;      // sum over steps of active peers
double dsmcMigChunkCount_ = 0.0;     // sum over steps of chunks processed
double dsmcMigTotalWall_ = 0.0;      // sum of (tEnd - tStart) for the flat path
// §35 P0 sub-items: unpack/append are inside deser; rebuild is inside post;
// upcWait happens outside the migration interval.  parcels out/in per rank.
double dsmcMigUnpackWall_ = 0.0;
double dsmcMigAppendWall_ = 0.0;
// H4/M2: async-path incremental occupancy append (and the debug validation
// inside it), timed separately from the synchronous chunk append above.
static double dsmcMigOccAppendWall_ = 0.0;
static double dsmcMigOccValidateWall_ = 0.0;

// §35 P2 (clean C2b): async-overlap path stages, kept separate from the
// mixed trigger accumulators so each path explains its own parent interval.
// classify+pack share one interval in the async begin; sizeX is inside prep
// (tSizeExEnd == tLocalPrepEnd), and finish covers wait+unpack+append+
// occupancy rebuild+post.
static double dsmcAsyncClassifyPackWall_ = 0.0;
static double dsmcAsyncPrepWall_ = 0.0;
static double dsmcAsyncPostReqWall_ = 0.0;
static double dsmcAsyncFinishWall_ = 0.0;
static double dsmcAsyncTotalWall_ = 0.0;
static double dsmcAsyncParcelsOut_ = 0.0;
static double dsmcAsyncParcelsIn_ = 0.0;
static double dsmcAsyncSendBytes_ = 0.0;
static double dsmcAsyncRecvBytes_ = 0.0;
static label dsmcAsyncCalls_ = 0;
// Sub-items of the async finish interval (needed to subtract the async
// contribution from the shared trigger accumulators for the sync table).
static double dsmcAsyncWaitWall_ = 0.0;
static double dsmcAsyncDeserWall_ = 0.0;
static double dsmcAsyncPostWall_ = 0.0;
double dsmcMigParcelsOut_ = 0.0;
double dsmcMigParcelsIn_ = 0.0;
// last printed values for the per-step upc deltas (one-step lag).
double dsmcMigUpcPostPrev_ = 0.0;
double dsmcMigUpcWaitPrev_ = 0.0;

// §36.2: bounded, persistent migration scratch buffers (flat path).
// DynamicList::setSize() keeps its allocation when the requested size is
// <= capacity (DynamicListI.H:186), so these survive across steps and
// remove the per-step allocate/free churn (the per-destination send lists
// plus the two receive chunk buffers).  On the 1200w production case the
// measured traffic is ~10-25 MB per call per rank (sendMB/recvMB over 502
// calls), so the retained scratch is ~50-100 MB per rank.  Capacity above
// replicatedMeshMigBufCapMB is released after each transfer as a guard
// against a DLB storm pinning memory.  `kept` is deliberately NOT
// persistent: its storage is handed to the cloud by
// takeMoveOrderedParcels().
List<DynamicList<dsmcParcel::TransferData>> dsmcMigSendBuf_;
List<DynamicList<dsmcParcel::TransferData>> dsmcMigRecvChunks_;
DynamicList<dsmcParcel*> dsmcMigDeleteBuf_;

// Lazy-read switch cache for the accessor
bool dsmcAsyncOverlapEnabledCache_ = false;
bool dsmcAsyncOverlapRead_ = false;

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


// * * * * * * * * clean N1: non-flat record framing  * * * * * * * * * * * //

// The legacy binary parcel stream was a bare concatenation of records, so
// the reader could not tell "end of buffer" from a truncated record: at the
// exact end of an IStringStream the next parcel read walked past the buffer
// and the token reader silently marked the stream bad (surfacing later as
// "Attempt to get back from bad stream").  Frame the stream explicitly:
//
//     <parcel count> <record length><record bytes> ...
//
// The count terminates the loop without relying on eof(), and the length
// isolates each record so a serializer mismatch is reported at its own
// record instead of corrupting the following ones.
static bool dsmcVerifyTransferRecords_ = false;

// Per-rank non-flat transfer accounting (reported by report()).
static double dsmcNFSentParcels_ = 0.0;
static double dsmcNFRecvParcels_ = 0.0;
static double dsmcNFSendBytes_ = 0.0;
static double dsmcNFRecvBytes_ = 0.0;
static double dsmcNFCalls_ = 0.0;
// Records whose write was immediately re-read and compared, and the
// synthetic stuck/tracked self-test counter.
static double dsmcNFRecordsVerified_ = 0.0;
static double dsmcNFStuckTrackedVerified_ = 0.0;
static Foam::label dsmcNFVerifyLimit_ = 0;
static Foam::label dsmcNFSyntheticDone_ = 0;
// H4 self-test mode (none|empty|singleCell|lastCell), run once at report.
static Foam::word dsmcOccSelfTestMode_("none");


//- Append one parcel as <length><record>; returns the record byte size.
static Foam::label dsmcWriteNonFlatRecord(Ostream& os, const dsmcParcel& p)
{
    OStringStream rec(IOstream::BINARY);
    p.writeBinaryFast(rec);
    const string s = rec.str();

    os << label(s.size()) << token::SPACE;
    if (s.size())
    {
        os.write(s.c_str(), std::streamsize(s.size()));
    }

    if
    (
        dsmcVerifyTransferRecords_
     && dsmcNFRecordsVerified_ < double(dsmcNFVerifyLimit_)
    )
    {
        ++dsmcNFRecordsVerified_;

        ISpanStream rs(s.c_str(), label(s.size()), IOstream::BINARY);
        dsmcParcel* back = dsmcParcel::readTransferRecord(p.mesh(), rs);

        const char* firstName = nullptr;
        scalar firstA = 0;
        scalar firstB = 0;
        const label nDiff =
            dsmcParcel::compareTransferFields
            (
                p, *back, firstName, firstA, firstB
            );
        delete back;
        if (nDiff != 0)
        {
            FatalErrorInFunction
                << "N1 non-flat record round-trip mismatch: " << nDiff
                << " field(s) differ, first " << (firstName ? firstName : "?")
                << " read " << firstA << " != written " << firstB
                << " (record " << s.size() << " bytes)"
                << abort(FatalError);
        }
    }

    // Synthetic state coverage: the flat record refuses stuck/tracked
    // parcels, so a real run never exercises those fields.  Once per rank
    // (and only with verification on) push a scratch copy carrying both
    // states through the same record path and compare every field.
    if
    (
        dsmcVerifyTransferRecords_
     && dsmcNFSyntheticDone_ < dsmcNFVerifyLimit_
    )
    {
        ++dsmcNFSyntheticDone_;

        dsmcParcel scratch(p);
        scratch.tracked().setMigrationState
        (
            true,                    // tracked
            17,                      // inPatchId
            true,                    // storePositions
            1.2345e-7,               // t0
            vector(1e-4, 2e-4, 3e-4),   // initial position
            vector(1.1e-4, 2.2e-4, 3.3e-4),
            vector(0.5e-5, -0.5e-5, 0.25e-5),
            vector(1e-9, 2e-9, 3e-9)
        );

        scalarField wt(2, 350.0);
        vectorField wv(2, vector(1, 0, 0));
        wt[1] = 400.0;
        wv[1] = vector(0, 1, 0);
        scratch.setStuck();
        scratch.stuck() = dsmcParcel::StuckParcel(wt, wv);

        OStringStream srec(IOstream::BINARY);
        scratch.writeBinaryFast(srec);
        const string ss = srec.str();

        ISpanStream srs(ss.c_str(), label(ss.size()), IOstream::BINARY);
        dsmcParcel* sback = dsmcParcel::readTransferRecord(p.mesh(), srs);

        const char* sFirst = nullptr;
        scalar sA = 0;
        scalar sB = 0;
        const label sDiff = dsmcParcel::compareTransferFields
        (
            scratch, *sback, sFirst, sA, sB
        );
        delete sback;

        if (sDiff != 0)
        {
            FatalErrorInFunction
                << "N1 stuck/tracked record round-trip mismatch: " << sDiff
                << " field(s) differ, first " << (sFirst ? sFirst : "?")
                << " read " << sA << " != written " << sB
                << abort(FatalError);
        }
        ++dsmcNFStuckTrackedVerified_;
    }

    return label(s.size());
}


//- Read the parcel count of one framed stream.
static Foam::label dsmcReadNonFlatCount(Istream& is, const label maxBytes)
{
    label nParcels = 0;
    is >> nParcels;
    if (nParcels < 0 || nParcels > maxBytes)
    {
        FatalErrorInFunction
            << "Invalid non-flat parcel count " << nParcels
            << " (buffer " << maxBytes << " bytes)"
            << abort(FatalError);
    }
    return nParcels;
}


//- Read one length-framed record and return the parcel.
static Foam::dsmcParcel* dsmcReadNonFlatRecord
(
    const polyMesh& mesh,
    Istream& is,
    const label maxBytes,
    const label recordi,
    const label nRecords
)
{
    static thread_local DynamicList<char> recScratch;

    label recLen = 0;
    is >> recLen;
    if (recLen <= 0 || recLen > maxBytes)
    {
        FatalErrorInFunction
            << "Invalid non-flat record length " << recLen
            << " at record " << recordi << " / " << nRecords
            << " (buffer " << maxBytes << " bytes)"
            << abort(FatalError);
    }

    recScratch.setSize(recLen);
    is.read(recScratch.data(), recLen);

    dsmcNFRecvParcels_ += 1.0;
    dsmcNFRecvBytes_ += double(recLen);

    ISpanStream rs(recScratch.data(), recLen, IOstream::BINARY);
    return dsmcParcel::readTransferRecord(mesh, rs);
}


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
    offPrep1_(-1.0), offSx1_(-1.0), offMig1_(-1.0),
    offSent_(-1), offRecv_(-1), offSentB_(-1.0), offRecvB_(-1.0),
    offloadEnabled_(false),
    offloadPlanMode_("off"),
    offloadMaxPairs_(2),
    offloadPairBytesCap_(256*1024*1024),
    offloadChunkParcels_(400000),
    offloadMinGain_(2.0),
    offloadSkipDLBSteps_(true),
    offloadStateDumpStep_(-1),
    offloadDumpDone_(false),
    offloadControlsRead_(false),
    offloadPlanId_(0),
    offloadLastRebalanceCount_(-1),
    offloadPairsFileOpen_(false),
    offloadIAmDonor_(false),
    offloadIAmHelper_(false),
    offloadPairDonor_(-1),
    offloadPairHelper_(-1),
    offloadPlanNPairs_(0),
    offloadPlanReject_(0),
    offloadNSelected_(0),
    offloadNRecv_(0),
    offloadNExecViolation_(0),
    offloadNInserted_(0),
    offloadPlanBytes_(0),
    offloadNTotal_(0),
    offloadNChunks_(0),
    offloadTDispatch_(0.0), offloadTSendDone_(0.0), offloadTRecvHeader_(0.0),
    offloadTRecvDone_(0.0), offloadTExecDone_(0.0), offloadTInsertDone_(0.0),
    offloadHdrSendReq_(MPI_REQUEST_NULL),
    offloadChunkSendReqs_(),
    offloadHdrRecvReq_(MPI_REQUEST_NULL),
    offloadChunkRecvReqs_(),
    offloadSendsPending_(false),
    offloadDispatchRate_(60e-9),
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
    migTakeOrderedParcels_(true),
    migReuseMigBuffers_(true),
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
    // §44 Route A: trend-aware rank-k
    trendMode_("off"),
    trendWindowSteps_(50),
    trendHorizonSteps_(50),
    trendTheta_(0.5),
    trendGainGd_(0.0),
    trendRelaxation_(0.5),
    trendMaxDeltaK_(1.0),
    trendResetChangedRatio_(0.0),
    rankKActive_(false),
    trendWorkSum_(0.0),
    trendSamples_(0),
    trendLastWindowStep_(0),
    trendSkipSample_(false),
    trendPrevBarEValid_(false),
    cellCostMode_("off"),
    cellCostDumpEvery_(5),
    cellCostDumpCount_(0),
    cellCostEventCount_(0),
    cellCostInterveneEvent_(2),
    cellCostInterveneDone_(false),
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

    // ---- §15.5: smoothed weight field ----
    dsmcSmoothWeightEnabled_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshDLBSmoothWeight", false);
    dsmcSmoothWeightTheta_ =
        mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBSmoothWeightTheta", 0.3);
    if (dsmcSmoothWeightEnabled_)
    {
        Info<< "Replicated mesh: DLB smoothed weight field enabled"
            << " (theta=" << dsmcSmoothWeightTheta_ << ")" << endl;
    }

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

    // ---- §44 Route A: trend-aware rank-k configuration ----
    // §44.5: rankK 模式复用 WorkFeedbackGain 作为 Gp、KMin/KMax 作为限幅；
    // legacy 比例反馈仅在 TrendMode=off 时生效（避免双重叠加以往反馈）。
    workFeedbackK_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshDLBWorkFeedbackK", false);
    cellWeightK_ = mesh_.time().controlDict().lookupOrDefault<label>
    (
        "replicatedMeshDLBCellWeight",
        0
    );
    trendMode_ = mesh_.time().controlDict().lookupOrDefault<word>
    (
        "replicatedMeshDLBTrendMode",
        "off"
    );
    if (trendMode_ != "off" && trendMode_ != "rankK")
    {
        FatalErrorInFunction
            << "Unknown replicatedMeshDLBTrendMode '" << trendMode_
            << "'. Valid values are off and rankK."
            << exit(FatalError);
    }
    trendWindowSteps_ = max
    (
        label(1),
        mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBTrendWindowSteps", 50)
    );
    trendHorizonSteps_ = max
    (
        label(1),
        mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBTrendHorizonSteps", 50)
    );
    trendTheta_ = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBTrendTheta",
        0.5
    );
    trendGainGd_ = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBTrendGain",
        0.0
    );
    trendRelaxation_ = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBTrendRelaxation",
        0.5
    );
    trendMaxDeltaK_ = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBTrendMaxDeltaK",
        1.0
    );
    trendResetChangedRatio_ =
        mesh_.time().controlDict().lookupOrDefault<scalar>
        (
            "replicatedMeshDLBTrendResetChangedRatio",
            0.0
        );
    if (trendTheta_ <= 0 || trendTheta_ > 1)
    {
        FatalErrorInFunction
            << "replicatedMeshDLBTrendTheta must be in (0,1], got "
            << trendTheta_ << exit(FatalError);
    }
    if (trendRelaxation_ <= 0 || trendRelaxation_ > 1)
    {
        FatalErrorInFunction
            << "replicatedMeshDLBTrendRelaxation must be in (0,1], got "
            << trendRelaxation_ << exit(FatalError);
    }
    if (trendMaxDeltaK_ < 0 || trendGainGd_ < 0)
    {
        FatalErrorInFunction
            << "replicatedMeshDLBTrendMaxDeltaK/Gain must be >= 0"
            << exit(FatalError);
    }
    if (trendMode_ == "rankK" && !triggerWallTime_)
    {
        Info<< "Replicated mesh: DLB trend rankK requires the wall trigger"
            << " metric; falling back to off" << endl;
        trendMode_ = "off";
    }
    // ---- §45 P0: per-cell cost 观测配置 ----
    cellCostMode_ = mesh_.time().controlDict().lookupOrDefault<word>
    (
        "replicatedMeshDLBCellCostMode",
        "off"
    );
    if
    (
        cellCostMode_ != "off"
     && cellCostMode_ != "observe"
     && cellCostMode_ != "weighted"
    )
    {
        FatalErrorInFunction
            << "Unknown replicatedMeshDLBCellCostMode '" << cellCostMode_
            << "'. Valid values are off, observe and weighted."
            << exit(FatalError);
    }
    // §45 P3: weighted 模式系数——vwgt = A*N + B*N^2 + C*Q + D（相对权重，
    // 与拟合窗口同长使用；B 建议置 0，P1 拟合的 SN2 与 Np 共线；
    // D 为每 cell 基重（build/post 空 cell 定价），**不可为 0**——首轮
    // A/B 无 D 时 4.3M 空 cells 权重 0，build/post 全堆少数 rank，
    // evolveMax 14.5×（§17 空 cell 教训在 per-cell 权重下的再现）。
    cellCostCoeffs_ = mesh_.time().controlDict().lookupOrDefault<List<scalar>>
    (
        "replicatedMeshDLBCellCostCoeffs",
        List<scalar>()
    );
    if (cellCostMode_ == "weighted" && nProcs_ > 1)
    {
        // clean H3: the Q term is accumulated per owning rank and read at
        // ParMETIS vertex-slice positions, so a multi-rank weighted run
        // would use Q=0 for every cell the rank does not own.  Refuse
        // instead of silently degrading the weight formula.
        FatalErrorInFunction
            << "replicatedMeshDLBCellCostMode weighted requires nProcs 1: "
            << "the per-cell Q term is not gathered across ranks yet "
            << "(see clean-20261101 H3).  Use the default k+N weights or a "
            << "single rank."
            << exit(FatalError);
    }
    if (cellCostMode_ == "weighted")
    {
        bool coeffsOk = (cellCostCoeffs_.size() == 4);
        forAll(cellCostCoeffs_, k)
        {
            const scalar v = cellCostCoeffs_[k];
            coeffsOk = coeffsOk && (v == v) && (v > -VGREAT) && (v < VGREAT);
        }
        if (!coeffsOk)
        {
            FatalErrorInFunction
                << "replicatedMeshDLBCellCostCoeffs must be 4 finite "
                << "scalars (A B C D) for weighted mode, got "
                << cellCostCoeffs_ << exit(FatalError);
        }
    }
    cellCostDumpEvery_ = max
    (
        label(1),
        mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBCellCostDumpEvery", 1)
    );
    // §45 P2: 干预入口配置（observe-only 实验用；默认空 = 关闭）。
    // 用 string 读取：路径含 "/" 时 word/fileName token 会拒绝。
    cellCostInterveneFile_ = mesh_.time().controlDict().lookupOrDefault<string>
    (
        "replicatedMeshDLBCellCostInterveneFile",
        ""
    );
    cellCostInterveneEvent_ = max
    (
        label(1),
        mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBCellCostInterveneEvent", 2)
    );
    cellCostInterveneDone_ = false;
    if (!cellCostInterveneFile_.empty() && cellCostMode_ != "observe")
    {
        FatalErrorInFunction
            << "replicatedMeshDLBCellCostInterveneFile requires "
            << "replicatedMeshDLBCellCostMode observe." << exit(FatalError);
    }
    if (cellCostMode_ != "off")
    {
        cellCostQ_.setSize(mesh_.nCells(), 0.0);
        cellCostTouchedMask_.setSize(mesh_.nCells(), false);
        Info<< "Replicated mesh: cellCost " << cellCostMode_ << " enabled"
            << " (dumpEvery=" << cellCostDumpEvery_
            << ", nCells=" << mesh_.nCells()
            << (cellCostMode_ == "observe"
              ? ", binary dump) "
              : ", per-cell vwgt) ")
            << endl;
        if (!cellCostInterveneFile_.empty())
        {
            Info<< "Replicated mesh: cellCost intervention armed"
                << " (file=" << cellCostInterveneFile_
                << ", event=" << cellCostInterveneEvent_ << ")" << endl;
        }
    }
    rankKActive_ =
        (trendMode_ == "rankK")
     || (workFeedbackK_ && trendMode_ == "off");
    if (trendMode_ == "rankK")
    {
        Info<< "Replicated mesh: DLB trend rankK enabled"
            << " (window=" << trendWindowSteps_
            << " horizon=" << trendHorizonSteps_
            << " theta=" << trendTheta_
            << " Gd=" << trendGainGd_
            << " relax=" << trendRelaxation_
            << " maxDeltaK=" << trendMaxDeltaK_
            << " resetChangedRatio=" << trendResetChangedRatio_ << ")" << endl;
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
    if (!useFlatTransfer_)
    {
        // clean N1: the legacy binary serializer round-trips member by
        // member now (dsmcParcel/TrackedParcel), so this path is supported;
        // it is slower than the flat POD path and only used when selected.
        Info<< "Replicated mesh: non-flat (binary stream) transfer path"
            << endl;
    }
    if (useFlatTransfer_)
    {
        Info<< "Replicated mesh: flat POD transfer enabled" << endl;
    }

    // H6/N1: verify every transfer record field by field (default off).
    // For the flat path the check packs and immediately unpacks a parcel;
    // for the non-flat path it re-reads each framed record.  Both compare
    // topology, velocity, weights, internal energy, species, vibrational
    // and stuck/tracked state, and abort on the first difference.
    {
        const bool verifyTransfer =
            mesh_.time().controlDict().lookupOrDefault<bool>
            ("replicatedMeshVerifyTransfer", false);
        const label verifyLimit =
            mesh_.time().controlDict().lookupOrDefault<label>
            ("replicatedMeshVerifyTransferLimit", 8);

        dsmcVerifyTransferRecords_ = verifyTransfer;
        dsmcNFVerifyLimit_ = max(label(0), verifyLimit);
        dsmcOccSelfTestMode_ =
            mesh_.time().controlDict().lookupOrDefault<word>
            (
                "replicatedMeshOccSelfTest",
                word("none")
            );
        dsmcNFRecordsVerified_ = 0.0;
        dsmcNFStuckTrackedVerified_ = 0.0;
        dsmcNFSyntheticDone_ = 0;
        dsmcParcel::verifyTransferRoundTrip(verifyTransfer, verifyLimit);

        if (verifyTransfer)
        {
            Info<< "Replicated mesh: transfer round-trip verification on"
                << " (first " << verifyLimit << " records per rank)" << endl;
        }
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

    // §36 A/B switches (default on).  A single binary can therefore
    // measure each §36 item against its pre-§36 control on the same
    // hardware, which node-placement variance makes necessary.
    migTakeOrderedParcels_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshMigTakeOrdered", true);
    migReuseMigBuffers_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        ("replicatedMeshMigReuseBuffers", true);
    Info<< "Replicated mesh: §36 migTakeOrdered=" << migTakeOrderedParcels_
        << " migReuseBuffers=" << migReuseMigBuffers_ << endl;

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

label dsmcReplicatedMesh::reassignByParMetisAdaptiveRepart(
    const scalar checkImbalance
)
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

    // Build LOCAL xadj/adjncy for my portion of the graph.
    // §21: the topology graph is static in replicated-mesh mode (the
    // local slice depends only on nCells/nProcs_) — cache it across DLB
    // events instead of re-scanning all internal faces every time.
    static List<idx_t> xadj, adjncy;
    static label graphBuiltCells = -1;
    if (graphBuiltCells != nCells)
    {
        labelList localDeg(myN, 0);
        for (label fi = 0; fi < nIntFaces; ++fi)
        {
            const label own = faceOwner[fi];
            const label nei = faceNei[fi];
            if (own >= myStart && own < myEnd) ++localDeg[own - myStart];
            if (nei >= myStart && nei < myEnd) ++localDeg[nei - myStart];
        }

        xadj.setSize(myN + 1);
        xadj[0] = 0;
        for (label i = 0; i < myN; ++i)
            xadj[i+1] = xadj[i] + localDeg[i];

        adjncy.setSize(xadj[myN]);
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
        graphBuiltCells = nCells;
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

    // §15.5: EMA smoothing of per-cell particle counts.  The smoothed
    // field replaces the instantaneous globalCellParticles for vwgt
    // construction, giving ParMETIS a stable weight field that captures
    // the underlying density trend without transient noise.
    if (dsmcSmoothWeightEnabled_)
    {
        if (dsmcSmoothCellParticles_.size() != nCells)
        {
            dsmcSmoothCellParticles_.setSize(nCells, 0.0);
            for (label c = 0; c < nCells; ++c)
            {
                dsmcSmoothCellParticles_[c] = scalar(globalCellParticles[c]);
            }
        }
        else
        {
            const scalar th = dsmcSmoothWeightTheta_;
            for (label c = 0; c < nCells; ++c)
            {
                dsmcSmoothCellParticles_[c] =
                    th*scalar(globalCellParticles[c])
                  + (1.0 - th)*dsmcSmoothCellParticles_[c];
            }
        }
        for (label c = 0; c < nCells; ++c)
        {
            globalCellParticles[c] =
                idx_t(max(dsmcSmoothCellParticles_[c], scalar(0)));
        }
    }

    if (dlbProfile && myRank_ == 0)
    {
        Info<< "Phase C ParMETIS profile: particle-count allreduce complete"
            << endl;
    }

    // §45 P0: 观测模式保存本帧 per-cell 粒子快照（scalar 副本仅用于
    // 诊断 dump；接受门/预算/守恒仍用 globalCellParticles 原值）
    if (cellCostMode_ == "observe")
    {
        if (cellCostP_.size() != nCells) cellCostP_.setSize(nCells, 0.0);
        for (label c = 0; c < nCells; ++c)
        {
            cellCostP_[c] = scalar(globalCellParticles[c]);
        }
    }

    // §45 P0: 观测模式——每 DLB 事件 dump 本窗口 per-cell (P, Q) 与
    // 每 rank 窗口净工作代理，然后清零窗口（§45.4；dumpEvery 控制 I/O 频率）。
    // §45 P3: weighted 模式不 dump，但每事件把已用于 vwgt 的 Q 窗口
    // 按触及表清零并推进事件计数。
    if (cellCostMode_ == "observe")
    {
        ++cellCostEventCount_;
        dumpCellCostWindow(label(mesh_.time().timeIndex()), globalCellParticles);
    }
    else if (cellCostMode_ == "weighted")
    {
        ++cellCostEventCount_;
        for (label i = 0; i < cellCostQtouched_.size(); ++i)
        {
            const label c = cellCostQtouched_[i];
            if (c >= 0 && c < cellCostQ_.size())
            {
                cellCostQ_[c] = 0.0;
            }
        }
        cellCostQtouched_.clear();
        for (label c = 0; c < cellCostTouchedMask_.size(); ++c)
        {
            cellCostTouchedMask_[c] = false;
        }
    }

    const idx_t maxParMetisWeight = std::numeric_limits<idx_t>::max()/4;
    auto positiveWeight = [&](const scalar value) -> idx_t
    {
        // §15.4 fix: empty cells (N=0) get weight 0, not 1.  Assigning
        // weight 1 to 4.3M empty cells creates 4.3M phantom weight that
        // distorts ParMETIS's balance: vwgt sums (with phantom) ≠ actual
        // particle counts (without phantom), so the post-migration
        // particle distribution is worse than the vwgt balance suggests.
        if (value <= scalar(0))
        {
            return idx_t(0);
        }
        if (value <= scalar(1))
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

    // §17: per-cell base-cost weight k (replicatedMeshDLBCellWeight).
    // Empty cells still cost build+post cell iteration (~354 ns/cell/step
    // vs ~26 ns/particle/step, 4826282 medians) — one empty cell is worth
    // ~14 particles of work.  w = k + N^alpha prices both channels in a
    // single constraint: the k*ownedCells term forbids empty-cell dumping,
    // the N term keeps particle balance.  k=0 restores weight-0 behaviour.
    // §44: k0 在构造时读取（cellWeightK_），趋势重置也需要它
    const label dlbCellWeightK = cellWeightK_;

    // §20: per-rank work feedback on the cell-weight base.
    // The uniform k=15 overprices each cell by (k - b/a) ~ 13.9 particle
    // equivalents; the measured work spread between the parked rank and
    // the gate equals a*(k - b/a)*dC (predicted 0.40 s vs 0.34 s
    // measured, worklog §19.2).  Feedback: measure each rank's inter-DLB
    // window work (wall - comm wait, §19.1 pairing, gathered here on the
    // partition this event is about to re-assign) and scale its
    // cell-weight base proportionally — heavy ranks' cells get pricier
    // so ParMETIS sheds them (edge-cut aware), light ranks' cells get
    // cheaper so neighbouring bulk cells flow in.  The global mean k
    // (=15, the anti-dumping behavioural term) is unchanged; only the
    // inter-rank differential moves.  Off by default: disabled == the
    // exact k15 production behaviour (§18 physical model reverted).
    // §44: legacy 比例反馈仅在 TrendMode=off 时生效；rankK 模式由
    // trendRankKUpdate() 维护 workFeedbackKR（含趋势项），避免双重叠加。
    const bool workFeedbackK =
        workFeedbackK_ && trendMode_ == "off";
    const scalar wfbGain = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBWorkFeedbackGain",
        2.0
    );
    const scalar wfbKMin = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBWorkFeedbackKMin",
        8.0
    );
    const scalar wfbKMax = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBWorkFeedbackKMax",
        25.0
    );

    List<scalar> workFeedbackKR(nProcs_, scalar(cellWeightK_));
    if (trendMode_ == "rankK")
    {
        workFeedbackKR = trendK_;   // §44: 趋势更新器维护的 per-rank k
    }
    else if (workFeedbackK && dsmcLastDLBWallCum_.size() == nProcs_)
    {
        double myWallCum = dsmcStepWallCum_;
        double myWaitCum = dsmcStepWaitCum_;
        List<scalar> wallNow(nProcs_, 0.0), waitNow(nProcs_, 0.0);
        MPI_Allgather
        (
            &myWallCum, 1, MPI_DOUBLE,
            wallNow.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD
        );
        MPI_Allgather
        (
            &myWaitCum, 1, MPI_DOUBLE,
            waitNow.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD
        );
        List<scalar> wR(nProcs_, 0.0);
        scalar wSum = 0.0;
        for (label r = 0; r < nProcs_; ++r)
        {
            wR[r] = max
            (
                (wallNow[r] - dsmcLastDLBWallCum_[r])
              - (waitNow[r] - dsmcLastDLBWaitCum_[r]),
                scalar(0)
            );
            wSum += wR[r];
        }
        const scalar wMean = wSum / nProcs_;
        for (label r = 0; r < nProcs_; ++r)
        {
            const scalar frac = (wMean > SMALL)
              ? (wR[r] - wMean) / wMean : 0.0;
            workFeedbackKR[r] = min
            (
                wfbKMax,
                max(wfbKMin, scalar(dlbCellWeightK)*(1.0 + wfbGain*frac))
            );
        }
        if (myRank_ == 0)
        {
            Info<< "Phase C ParMETIS work-feedback k:";
            for (label r = 0; r < nProcs_; ++r)
            {
                Info<< " r" << r << "=" << workFeedbackKR[r];
            }
            Info<< endl;
        }
    }

    // ---- §21: proposal acceptance machinery --------------------------------
    // Calibrated three-channel prediction (worklog §17.6.2):
    //   What_r = (a*N_r + b*(beta*SumN2_r) + c*C_r) * nStepsInWindow
    // with a=31.6 ns/particle, b=124 ns/candidate, c=34 ns/cell and beta
    // the runtime candidate rate (global candidates per step / sum_c
    // N_c^2, re-derived per event).  The measured per-rank window work
    // (wall - comm wait, §19.1) minus the current-partition prediction
    // gives a per-rank residual that is transferred with the cells, so
    // the proposal prediction carries the un-modelled spatial price.
    // The §20 lesson is encoded in the acceptance test: the predicted
    // MAX work must drop (a tighter spread with a larger total is
    // rejected), not the spread itself.
    // §26: production default ON.  3x3 A/B on the magazine+scratch
    // binary: accON 598.1 +-13 s vs accOFF 639.2 +-27 s (total_loop
    // max) — -6.4%, non-overlapping distributions, mechanism verified
    // (2 skipped migrations x ~20 s storm cost).  See worklog §26.4.
    const bool acceptCheck = mesh_.time().controlDict().lookupOrDefault<bool>
    (
        "replicatedMeshDLBAcceptCheck",
        true
    );
    const scalar accCoeffA = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBAcceptCoeffA",
        31.6
    );
    const scalar accCoeffB = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBAcceptCoeffB",
        124.0
    );
    const scalar accCoeffC = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBAcceptCoeffC",
        34.0
    );
    const scalar accEpsRel = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBAcceptEpsRel",
        0.01
    );
    const scalar accMaxMigFrac = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBAcceptMaxMigFrac",
        0.3
    );
    // §21.7: migration budget in particles (DSMC migration cost scales
    // with parcels, not cells) and the payback factor (predicted gate
    // gain over the window vs the measured average DLB execution cost).
    const scalar accMaxMigPartFrac = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBAcceptMaxMigPartFrac",
        0.15
    );
    const scalar accPaybackFactor = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBAcceptPaybackFactor",
        1.0
    );

    scalar accBeta = 0.0;
    List<scalar> accWorkMeas(nProcs_, 0.0);    // measured window work (s)
    List<scalar> accResidPerCell(nProcs_, 0.0);
    List<label> accNCur(nProcs_, 0);
    List<label> accCCur(nProcs_, 0);
    List<scalar> accSN2Cur(nProcs_, 0.0);
    label accWindowSteps = 0;
    double accCandDelta = 0.0;
    bool accHaveMeas = false;
    if (acceptCheck)
    {
        // §21.7 fix (user audit point 3): the candidate/time baselines
        // advance at EVERY DLB event — including the first, where no
        // work baseline exists yet — so candDelta, accWindowSteps and
        // the measured-work window always share the same start/end.
        // The previous version left them at zero until the first
        // acceptance, so event 2 scaled the counts model by the
        // startup-to-event step count against a post-event measured
        // window (150-step model vs ~100-step work at events 50/150).
        const double candNow = double(cloud_.cumulativeCollisionCandidates());
        double candTotal = 0.0;
        MPI_Allreduce
        (
            const_cast<double*>(&candNow), &candTotal,
            1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD
        );
        accCandDelta = candTotal - dsmcAccCandLastTotal_;
        accWindowSteps =
            label(mesh_.time().timeIndex()) - dsmcAccCandLastTimeIdx_;
        dsmcAccCandLastTotal_ = candTotal;
        dsmcAccCandLastTimeIdx_ = label(mesh_.time().timeIndex());

        if (dsmcLastDLBWallCum_.size() == nProcs_ && accWindowSteps > 0)
        {
            // One serial pass: window-sum_c N_c^2 (for beta) + current
            // per-partition counts.  Serial by choice — 13M simple
            // iterations once per event, no reduction complexity.
            double sumN2Acc = 0.0;
            for (label c = 0; c < nCells; ++c)
            {
                const double n = double(globalCellParticles[c]);
                const label own = cellOwner_[c];
                sumN2Acc += n * n;
                accNCur[own] += label(n);
                ++accCCur[own];
                accSN2Cur[own] += n * n;
            }
            if (sumN2Acc > 0.0 && accCandDelta > 0.0)
            {
                accBeta = accCandDelta / (double(accWindowSteps) * sumN2Acc);
            }

            // measured per-rank window work (wall - comm wait)
            double myWallCum = dsmcStepWallCum_;
            double myWaitCum = dsmcStepWaitCum_;
            List<scalar> wallNow(nProcs_, 0.0), waitNow(nProcs_, 0.0);
            MPI_Allgather
            (
                &myWallCum, 1, MPI_DOUBLE,
                wallNow.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD
            );
            MPI_Allgather
            (
                &myWaitCum, 1, MPI_DOUBLE,
                waitNow.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD
            );
            for (label r = 0; r < nProcs_; ++r)
            {
                accWorkMeas[r] = max
                (
                    (wallNow[r] - dsmcLastDLBWallCum_[r])
                  - (waitNow[r] - dsmcLastDLBWaitCum_[r]),
                    scalar(0)
                );
            }
            // residual per owned cell (measured - counts-model, per window)
            for (label r = 0; r < nProcs_; ++r)
            {
                const scalar predR =
                    (
                        accCoeffA * scalar(accNCur[r])
                      + accCoeffB * accBeta * accSN2Cur[r]
                      + accCoeffC * scalar(accCCur[r])
                    ) * scalar(accWindowSteps) * 1.0e-9;
                accResidPerCell[r] =
                    (accCCur[r] > 0)
                  ? (accWorkMeas[r] - predR) / scalar(accCCur[r])
                  : scalar(0);
            }
            accHaveMeas = true;
        }
    }

    idx_t ncon = 1;
    idx_t wgtflag = 2;  // vertex weights only

    List<idx_t> vwgt(myN * ncon, 1);

    // §45 P3: weighted 模式——vwgt = A*N + B*N^2 + C*Q（cellCostCoeffs_）。
    // Q 为上一窗口的候选累计（平稳性假设，与 P2 干预同口径）；
    // 空且无候选的 cell 权重 0（positiveWeight 语义）。
    if (cellCostMode_ == "weighted")
    {
        const scalar cA = cellCostCoeffs_[0];
        const scalar cB = cellCostCoeffs_[1];
        const scalar cC = cellCostCoeffs_[2];
        const scalar cD = cellCostCoeffs_[3];
        #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static)
        for (label i = 0; i < myN; ++i)
        {
            const label gi = myStart + i;
            const double n = double(globalCellParticles[gi]);
            const double q = (gi < cellCostQ_.size())
              ? double(cellCostQ_[gi]) : 0.0;
            vwgt[i] = positiveWeight(cA*n + cB*n*n + cC*q + cD);
        }
    }
    else
    #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static)
    for (label i = 0; i < myN; ++i)
    {
        const label gi = myStart + i;
        const idx_t nPart = globalCellParticles[gi];

        const scalar w = (nPart > 1)
            ? std::pow(scalar(nPart), dlbAlpha_) : scalar(nPart);
        const label owner = cellOwner_[gi];
        const scalar kR =
            (rankKActive_ && owner >= 0 && owner < nProcs_)
          ? workFeedbackKR[owner] : scalar(dlbCellWeightK);
        vwgt[i] = idx_t(kR) + positiveWeight(w);
    }

    // clean M9: the old adjwgt build (per-internal-face sum of the two
    // cells' particle counts) was never passed to ParMETIS (wgtflag=2, the
    // edge-weight argument is nullptr) and cost O(nInternalFaces) per
    // event.  Removed; edge weights are not part of this weighting model.

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
        ("replicatedMeshDLBUbvec", 1.01);
    List<real_t> ubvec(ncon, ubvecVal);

    Info<< "Phase C ParMETIS: alpha=" << dlbAlpha_ << " ncon=" << ncon
        << " cellK=" << dlbCellWeightK
        << " rankK=" << (rankKActive_ ? "on" : "off")
        << " trendMode=" << trendMode_
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
    // comm_further_opt.md §7: hybrid ParMETIS mode.  AdaptiveRepart
    // refines the current partition but does NOT reliably enforce ubvec
    // when the current partition is severely imbalanced (zero-lag
    // measurement: max/min 1.48 at ubvec 1.05).  PartKway finds the
    // globally optimal balanced partition from scratch and strictly
    // enforces ubvec.  Switch: checkImbalance > hybridThreshold →
    // PartKway (severe imbalance, needs global rebalance); otherwise
    // AdaptiveRepart (mild, incremental is cheaper).
    const scalar hybridThreshold = mesh_.time().controlDict().lookupOrDefault<scalar>
        ("replicatedMeshDLBHybridThreshold", 1.6);
    const bool usePartKway =
        checkImbalance > hybridThreshold
     && nProcs_ > 1;

    // §45 P2: 干预事件——owner 覆盖文件直接作为本事件 proposal（跳过
    // ParMETIS），后续事件恢复自动路径。覆盖表每 rank 解析同一文件，
    // 只改自己 vtxdist 段的 part[]；迁移执行路径与 ParMETIS proposal 相同。
    bool interveneNow =
    (
        cellCostMode_ == "observe"
     && !cellCostInterveneFile_.empty()
     && !cellCostInterveneDone_
     && cellCostEventCount_ == cellCostInterveneEvent_
     && nProcs_ > 1
    );
    if (interveneNow)
    {
        cellCostInterveneDone_ = true;
        cellCostInterveneOwner_.setSize(nCells, -1);
        IFstream is(cellCostInterveneFile_);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot open cellCost intervention file: "
                << cellCostInterveneFile_ << exit(FatalError);
        }
        label ci, oi, nOver = 0;
        while (is.read(ci).read(oi).good())
        {
            if (ci < 0 || ci >= nCells || oi < 0 || oi >= nProcs_)
            {
                FatalErrorInFunction
                    << "Intervention entry out of range: cell=" << ci
                    << " owner=" << oi << exit(FatalError);
            }
            cellCostInterveneOwner_[ci] = oi;
            ++nOver;
        }
        for (label i = 0; i < myN; ++i)
        {
            const label gi = myStart + i;
            if (cellCostInterveneOwner_[gi] >= 0)
            {
                part[i] = idx_t(cellCostInterveneOwner_[gi]);
            }
        }
        Info<< "Phase C cellCost intervention: applied at event "
            << cellCostEventCount_ << " (" << nOver << " cell overrides)"
            << endl;
    }
    else if (usePartKway)
    {
        // PartKway: edgecut minimization with HARD balance constraint.
        // vsize is not used (no migration-cost penalty — we WANT the
        // freedom to move cells anywhere to achieve balance).
        if (dlbProfile && myRank_ == 0)
        {
            Info<< "Phase C ParMETIS: PartKway (severe imbalance "
                << checkImbalance << " > " << hybridThreshold << ")" << endl;
        }
        ParMETIS_V3_PartKway
        (
            vtxdist.data(), xadj.data(), adjncy.data(),
            vwgt.data(), nullptr, &wgtflag, &numflag, &ncon, &nparts,
            tpwgts.data(), ubvec.data(), options, &edgecut,
            part.data(), &comm
        );
    }
    else
    {
        ParMETIS_V3_AdaptiveRepart
        (
            vtxdist.data(), xadj.data(), adjncy.data(),
            vwgt.data(), vsize.data(), nullptr,
            &wgtflag, &numflag, &ncon, &nparts,
            tpwgts.data(), ubvec.data(), &itr,
            options, &edgecut, part.data(), &comm
        );
    }
    if (dlbProfile && myRank_ == 0)
    {
        Info<< "Phase C ParMETIS profile: "
            << (usePartKway ? "PartKway" : "AdaptiveRepart")
            << " complete" << endl;
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

    // §8/clean M10: global vwgt for the diagnostic chain.  Reconstruction
    // of the k+N path only: multi-rank weighted mode is refused (H3), and
    // the integer truncation below mirrors the vwgt actually passed to
    // ParMETIS (idx_t(kR) + positiveWeight(w)).
    // §20: mirrors the work-feedback per-rank cell-weight base.
    List<scalar> globalVWgt(nCells, 0.0);
    for (label c = 0; c < nCells; ++c)
    {
        const idx_t nPart = globalCellParticles[c];
        scalar w = (nPart > 1)
            ? std::pow(scalar(nPart), dlbAlpha_) : scalar(nPart);
        const label owner = cellOwner_[c];
        const scalar kR =
            (rankKActive_ && owner >= 0 && owner < nProcs_)
          ? workFeedbackKR[owner] : scalar(dlbCellWeightK);
        globalVWgt[c] = scalar(idx_t(kR) + positiveWeight(w));
    }

    // ---- comm_further_opt.md §8: ParMETIS chain diagnostics ----
    // Four-point measurement to locate where balance is lost:
    // (A) ParMETIS input: per-part vwgt sums from the CURRENT partition
    // (B) ParMETIS output: per-part vwgt sums from the PROPOSED partition
    // (C) After remap: per-part vwgt sums after remap modification
    // Combined with the zero-lag post-migration table (actual particle
    // counts), this isolates: input problem / ParMETIS problem / remap
    // problem / migration problem.
    if (dlbProfile && myRank_ == 0)
    {
        // Per-cell vwgt is only on the local rank; assemble global vwgt
        // from the current partition (cellOwner_) and the proposed
        // partition (fullPart).  For each, sum vwgt per part.
        List<scalar> curSum(nProcs_, 0.0);
        List<scalar> propSum(nProcs_, 0.0);

        for (label c = 0; c < nCells; ++c)
        {
            const label curOwner = cellOwner_[c];
            const label propOwner = label(fullPart[c]);
            if (curOwner >= 0 && curOwner < nProcs_)
                curSum[curOwner] += globalVWgt[c];
            if (propOwner >= 0 && propOwner < nProcs_)
                propSum[propOwner] += globalVWgt[c];
        }

        scalar curMax = 0, curMin = GREAT;
        scalar propMax = 0, propMin = GREAT;
        for (label r = 0; r < nProcs_; ++r)
        {
            if (curSum[r] > curMax) curMax = curSum[r];
            if (curSum[r] < curMin) curMin = curSum[r];
            if (propSum[r] > propMax) propMax = propSum[r];
            if (propSum[r] < propMin) propMin = propSum[r];
        }

        Info<< "  [ParMETIS chain] input(current) vwgt max/min = "
            << curMax / max(curMin, SMALL) << nl
            << "  [ParMETIS chain] output(proposed) vwgt max/min = "
            << propMax / max(propMin, SMALL) << nl;

        // Also print per-part vwgt sums (proposed) for detailed inspection
        Info<< "    per-part vwgt (proposed):";
        for (label r = 0; r < nProcs_; ++r)
        {
            Info<< " r" << r << "=" << propSum[r];
        }
        Info<< nl;

        // §8 addendum: DIRECT particle-count balance (not vwgt).
        // This is the ground truth: does ParMETIS's proposed partition
        // balance the ACTUAL per-cell particle counts?
        List<label> propPartSum(nProcs_, 0);
        List<label> curPartSum(nProcs_, 0);
        for (label c = 0; c < nCells; ++c)
        {
            const label propOwner = label(fullPart[c]);
            const label curOwner = cellOwner_[c];
            const label actualN = label(globalCellParticles[c]);
            if (propOwner >= 0 && propOwner < nProcs_)
                propPartSum[propOwner] += actualN;
            if (curOwner >= 0 && curOwner < nProcs_)
                curPartSum[curOwner] += actualN;
        }
        label propPMax = 0, propPMin = labelMax;
        label curPMax = 0, curPMin = labelMax;
        for (label r = 0; r < nProcs_; ++r)
        {
            if (propPartSum[r] > propPMax) propPMax = propPartSum[r];
            if (propPartSum[r] < propPMin) propPMin = propPartSum[r];
            if (curPartSum[r] > curPMax) curPMax = curPartSum[r];
            if (curPartSum[r] < curPMin) curPMin = curPartSum[r];
        }
        Info<< "  [ParMETIS chain] INPUT particles max/min = "
            << scalar(curPMax)/scalar(max(curPMin, 1)) << nl
            << "  [ParMETIS chain] OUTPUT particles max/min = "
            << scalar(propPMax)/scalar(max(propPMin, 1)) << nl;
        Info<< "    per-part particles (proposed):";
        for (label r = 0; r < nProcs_; ++r)
        {
            Info<< " r" << r << "=" << propPartSum[r];
        }
        Info<< nl;

        // §15.24: total owned cells per rank (NOT ownedCollCells — this
        // is the total cells assigned to each rank by the partition).
        // Needed to separate "more cells" from "more particles per cell"
        // as causes of move imbalance.
        List<label> propCellCount(nProcs_, 0);
        List<label> curCellCount(nProcs_, 0);
        for (label c = 0; c < nCells; ++c)
        {
            const label propOwner = label(fullPart[c]);
            const label curOwner = cellOwner_[c];
            if (propOwner >= 0 && propOwner < nProcs_)
                propCellCount[propOwner]++;
            if (curOwner >= 0 && curOwner < nProcs_)
                curCellCount[curOwner]++;
        }
        Info<< "    per-part total owned cells (proposed):";
        for (label r = 0; r < nProcs_; ++r)
        {
            Info<< " r" << r << "=" << propCellCount[r];
        }
        Info<< nl;
        Info<< "    per-part total owned cells (current):";
        for (label r = 0; r < nProcs_; ++r)
        {
            Info<< " r" << r << "=" << curCellCount[r];
        }
        Info<< nl;

        // §17.6.1: per-part ACTIVE cells (N>0).  The cell-driven cost
        // (build/post) prices empty cells at roughly half, so the
        // three-channel work model needs the active/empty split per
        // part; total-only cells are collinear with N across parked
        // ranks and cannot be identified in the regression.
        List<label> propActCellCount(nProcs_, 0);
        List<label> curActCellCount(nProcs_, 0);
        for (label c = 0; c < nCells; ++c)
        {
            if (globalCellParticles[c] <= 0) continue;
            const label propOwner = label(fullPart[c]);
            const label curOwner = cellOwner_[c];
            if (propOwner >= 0 && propOwner < nProcs_)
                propActCellCount[propOwner]++;
            if (curOwner >= 0 && curOwner < nProcs_)
                curActCellCount[curOwner]++;
        }
        Info<< "    per-part active owned cells (proposed):";
        for (label r = 0; r < nProcs_; ++r)
        {
            Info<< " r" << r << "=" << propActCellCount[r];
        }
        Info<< nl;
        Info<< "    per-part active owned cells (current):";
        for (label r = 0; r < nProcs_; ++r)
        {
            Info<< " r" << r << "=" << curActCellCount[r];
        }
        Info<< nl;
    }


    // §45 P2: 干预事件的 proposal 是显式 per-cell 指定，remap 的标签
    // 置换会把它整体扭曲（首验：103093 指定被 remap 扩成 573759 变更、
    // 干预目标 rank 被置换抵消）→ 干预事件禁用 remap。
    const bool remapEnabled =
        !interveneNow
     && mesh_.time().controlDict().lookupOrDefault<bool>
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
        // §23: particle-aware greedyOverlap (roadmap direction #2/#3).
        // The historical metric counted CELLS; the actual event-migration
        // cost is carried by the PARTICLES, so each cell now contributes
        // (1 + alpha*N_c) to the weighted overlap.  The label bijection
        // does not change the partition shape nor the per-step boundary
        // faces (steady traffic is labeling-invariant), so this is a
        // pure event-migration saving with no balance downside.  The
        // unit-weight labelling is still computed as an IN-RUN
        // counterfactual (same partition, same run — immune to the +-3%
        // cross-run noise), and the acceptance check remains the safety
        // net for the resulting proposal.
        const scalar remapPartW = mesh_.time().controlDict().lookupOrDefault<scalar>
        (
            "replicatedMeshDLBRemapParticleWeight",
            1.0
        );

        // Per-thread flattened nProcs x nProcs overlap matrices (unit and
        // particle-weighted), merged serially (M1 two-pass pattern) to
        // avoid write conflicts on the shared matrix entries.
        List<labelList> threadOverlap(nT);
        List<scalarList> threadOverlapP(nT);
        for (label tI = 0; tI < nT; ++tI)
        {
            threadOverlap[tI].setSize(nProcs_*nProcs_, 0);
            threadOverlapP[tI].setSize(nProcs_*nProcs_, 0);
        }

        #pragma omp parallel num_threads(nT) if (nT > 1)
        {
            #ifdef _OPENMP
            const label tid = omp_get_thread_num();
            #else
            const label tid = 0;
            #endif
            labelList& myOverlap = threadOverlap[tid];
            scalarList& myOverlapP = threadOverlapP[tid];

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
                    myOverlapP[newPart*nProcs_ + oldRank] +=
                        scalar(1) + remapPartW*scalar(globalCellParticles[i]);
                }
            }
        }

        List<scalarList> overlap(nProcs_);
        List<scalarList> overlapP(nProcs_);
        for (label newPart = 0; newPart < nProcs_; ++newPart)
        {
            overlap[newPart].setSize(nProcs_, 0.0);
            overlapP[newPart].setSize(nProcs_, 0.0);
        }
        for (label tI = 0; tI < nT; ++tI)
        {
            const labelList& myOverlap = threadOverlap[tI];
            const scalarList& myOverlapP = threadOverlapP[tI];
            for (label newPart = 0; newPart < nProcs_; ++newPart)
            {
                for (label oldRank = 0; oldRank < nProcs_; ++oldRank)
                {
                    overlap[newPart][oldRank] +=
                        scalar(myOverlap[newPart*nProcs_ + oldRank]);
                    overlapP[newPart][oldRank] +=
                        myOverlapP[newPart*nProcs_ + oldRank];
                }
            }
        }

        // §23: one greedy matcher, applied to both matrices — the
        // weighted one chooses the applied labelling, the unit one is
        // the counterfactual baseline.
        auto greedyMatch =
            [&](const List<scalarList>& w, labelList& out)
        {
            out.setSize(nProcs_, -1);
            labelList newPartUsed(nProcs_, 0);
            labelList oldRankUsed(nProcs_, 0);

            for (label matchI = 0; matchI < nProcs_; ++matchI)
            {
                label bestNewPart = -1;
                label bestOldRank = -1;
                scalar bestVal = -1.0;

                for (label newPart = 0; newPart < nProcs_; ++newPart)
                {
                    if (newPartUsed[newPart]) continue;

                    for (label oldRank = 0; oldRank < nProcs_; ++oldRank)
                    {
                        if (oldRankUsed[oldRank]) continue;

                        const scalar value = w[newPart][oldRank];
                        if
                        (
                            value > bestVal
                         || (
                                value == bestVal
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
                            bestVal = value;
                            bestNewPart = newPart;
                            bestOldRank = oldRank;
                        }
                    }
                }

                if (bestNewPart >= 0 && bestOldRank >= 0)
                {
                    out[bestNewPart] = bestOldRank;
                    newPartUsed[bestNewPart] = 1;
                    oldRankUsed[bestOldRank] = 1;
                }
            }

            for (label newPart = 0; newPart < nProcs_; ++newPart)
            {
                if (out[newPart] >= 0) continue;

                for (label oldRank = 0; oldRank < nProcs_; ++oldRank)
                {
                    if (!oldRankUsed[oldRank])
                    {
                        out[newPart] = oldRank;
                        oldRankUsed[oldRank] = 1;
                        break;
                    }
                }
            }
        };

        labelList remapCells, remapPart;
        greedyMatch(overlap, remapCells);
        greedyMatch(overlapP, remapPart);

        // §23 counterfactual: migrated particles under each labelling
        // (raw ParMETIS labels still intact — the apply pass is below).
        // §47: diagnostic-only scan, parallelised (per-cell accumulation
        // is independent; reduction over the four counters).
        double migPartCells = 0.0;
        double migPartWeighted = 0.0;
        label changedCellsCF = 0;
        label changedCellsW = 0;
        #pragma omp parallel for num_threads(nT) if (nT > 1) schedule(static) \
            reduction(+:migPartCells, migPartWeighted, changedCellsCF, changedCellsW)
        for (label c = 0; c < nCells; ++c)
        {
            const label rawPart = label(fullPart[c]);
            if (rawPart < 0 || rawPart >= nProcs_) continue;

            const double n = double(globalCellParticles[c]);
            const label cur = cellOwner_[c];
            if (remapCells[rawPart] != cur)
            {
                ++changedCellsCF;
                migPartCells += n;
            }
            if (remapPart[rawPart] != cur)
            {
                ++changedCellsW;
                migPartWeighted += n;
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
             && remapPart[newPart] >= 0
            )
            {
                const label newOwner = remapPart[newPart];
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
            << nl
            << "  [§23 remap counterfactual] cell-count label: cellsChanged="
            << changedCellsCF << " migratedParticles=" << migPartCells
            << " | particle-weighted label: cellsChanged="
            << changedCellsW << " migratedParticles=" << migPartWeighted
            << " | particle reduction="
            << (migPartCells > 0.0
                ? (1.0 - migPartWeighted/migPartCells) : 0.0)
            << nl
            << "  remap=(";
        forAll(remapPart, i)
        {
            Info<< (i ? " " : "") << i << "->" << remapPart[i];
        }
        Info<< ")" << endl;
    }
    else if (remapEnabled)
    {
        WarningInFunction
            << "Unsupported replicatedMeshDLBRemapMode '" << remapMode
            << "'. Falling back to ParMETIS part labels." << endl;
    }

    // ---- §21: proposal acceptance decision ---------------------------------
    // Predict per-rank work for the proposed partition (counts + beta +
    // residual transferred with the cells) and compare against the
    // MEASURED max of the current window.  §21.7 (user audit): the
    // predicted total is conserved by construction — the acceptance can
    // only test max-redistribution under fixed per-cell costs, so the
    // decision also requires an explicit payback: the predicted gate
    // gain must exceed the estimated migration+rebuild cost, and the
    // migration is budgeted in PARTICLES (not cells — DSMC migration
    // cost scales with parcels).
    label nChanged = 0;
    bool proposalAccepted = true;
    if (acceptCheck && accHaveMeas)
    {
        List<label> nProp(nProcs_, 0);
        List<label> cProp(nProcs_, 0);
        List<scalar> sn2Prop(nProcs_, 0.0);
        List<scalar> residProp(nProcs_, 0.0);
        label nDiff = 0;
        double nMigParticles = 0.0;
        for (label c = 0; c < nCells; ++c)
        {
            const label prop = label(fullPart[c]);
            const double n = double(globalCellParticles[c]);
            nProp[prop] += label(n);
            ++cProp[prop];
            sn2Prop[prop] += n * n;
            residProp[prop] += accResidPerCell[cellOwner_[c]];
            if (prop != cellOwner_[c])
            {
                ++nDiff;
                nMigParticles += n;
            }
        }
        scalar wMaxOld = 0.0;
        scalar wSumOld = 0.0;
        for (label r = 0; r < nProcs_; ++r)
        {
            wMaxOld = max(wMaxOld, accWorkMeas[r]);
            wSumOld += accWorkMeas[r];
        }
        const scalar eps = accEpsRel * wSumOld / nProcs_;
        scalar wMaxNew = 0.0;
        for (label r = 0; r < nProcs_; ++r)
        {
            const scalar wNew =
                (
                    accCoeffA * scalar(nProp[r])
                  + accCoeffB * accBeta * sn2Prop[r]
                  + accCoeffC * scalar(cProp[r])
                ) * scalar(accWindowSteps) * 1.0e-9
              + residProp[r];
            wMaxNew = max(wMaxNew, wNew);
        }
        const bool migOkCells =
            scalar(nDiff) <= accMaxMigFrac * scalar(nCells);
        // §21.7: particle-based migration budget
        double nPartTotal = 0.0;
        for (label r = 0; r < nProcs_; ++r) nPartTotal += nProp[r];
        const bool migOkPart =
            nMigParticles
         <= accMaxMigPartFrac * max(nPartTotal, 1.0);
        // §21.7: payback — the predicted gate gain over the window must
        // exceed the measured average migration+rebuild cost of one DLB
        // execution (margins folded into accEpsRel/paybackFactor).
        const scalar wGain = wMaxOld - wMaxNew;
        // §31: the payback gate must be GLOBALLY consistent.  The
        // per-rank migration+rebuild estimate uses local wall clocks and
        // differs across ranks, so a borderline proposal could pass on
        // one rank and fail on another — the owner map would split and
        // validateCellOwnerMap() would abort the run.  Reduce the
        // threshold (MAX = worst-case cost) so every rank takes the same
        // branch.  The estimate also includes the local rebuild
        // (rebuildMyCells + localMesh_.build) measured at the previous
        // accepted event, which the migration timer does not cover.
        const scalar tMigLocal =
            autoRebalanceMigrationWallTime_
          / max(autoRebalanceCount_, label(1))
          + dsmcRebuildLastWall_;
        scalar tMigPerEvent = tMigLocal;
        scalar tMigLocalMin = tMigLocal;
        MPI_Allreduce
        (
            MPI_IN_PLACE, &tMigPerEvent, 1, MPI_DOUBLE,
            MPI_MAX, MPI_COMM_WORLD
        );
        MPI_Allreduce
        (
            MPI_IN_PLACE, &tMigLocalMin, 1, MPI_DOUBLE,
            MPI_MIN, MPI_COMM_WORLD
        );
        const bool paybackOk = wGain > accPaybackFactor * tMigPerEvent;
        // §33: consistency/exposure evidence.  tMigLocalMax (= used) and
        // tMigLocalMin span the per-rank estimate spread; the band
        // [Factor*min, Factor*max] is where a per-rank (unreduced)
        // threshold would have split the decision across ranks.  Logging
        // the spread plus the position of wGain documents both the
        // reduction actually in force and this proposal's exposure.
        const scalar thrUsed = accPaybackFactor * tMigPerEvent;
        const scalar thrPermissive = accPaybackFactor * tMigLocalMin;
        const bool inDivergenceBand =
            (wGain > thrPermissive) && (wGain <= thrUsed);
        proposalAccepted = (nDiff == 0)
            || ((wMaxNew < wMaxOld - eps) && migOkCells && migOkPart
                && paybackOk);
        if (myRank_ == 0)
        {
            Info<< "Phase C ParMETIS acceptance:"
                << " windowSteps=" << accWindowSteps
                << " candDelta=" << accCandDelta
                << " beta=" << accBeta
                << " Wmax_old=" << wMaxOld
                << " Wmax_pred_new=" << wMaxNew
                << " gain=" << wGain
                << " tMigPerEvent=" << tMigPerEvent
                << " (perRankMin=" << tMigLocalMin
                << " spread=" << (tMigPerEvent - tMigLocalMin) << ")"
                << " thrUsed=" << thrUsed
                << " thrMin=" << thrPermissive
                << " inDivergenceBand=" << (inDivergenceBand ? 1 : 0)
                << " eps=" << eps
                << " nChanged=" << nDiff
                << " migFrac=" << scalar(nDiff)/scalar(nCells)
                << " migParticles=" << nMigParticles
                << " -> " << (proposalAccepted ? "ACCEPT" : "REJECT")
                << nl;
        }
        if (!proposalAccepted)
        {
            Info<< "Phase C ParMETIS proposal rejected — keeping the"
                << " current partition (predicted max gain does not"
                << " pay back the migration cost)" << endl;
        }
    }
    if (proposalAccepted)
    {
        // ---- Update cellOwner_ --------------------------------------------
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
    }

    // ---- §8: after remap, measure per-part vwgt sums ----
    if (dlbProfile && myRank_ == 0 && nChanged > 0)
    {
        List<scalar> postRemapSum(nProcs_, 0.0);
        for (label c = 0; c < nCells; ++c)
        {
            const label own = cellOwner_[c];
            if (own >= 0 && own < nProcs_)
                postRemapSum[own] += globalVWgt[c];
        }
        scalar pMax = 0, pMin = GREAT;
        for (label r = 0; r < nProcs_; ++r)
        {
            if (postRemapSum[r] > pMax) pMax = postRemapSum[r];
            if (postRemapSum[r] < pMin) pMin = postRemapSum[r];
        }
        Info<< "  [ParMETIS chain] after-remap vwgt max/min = "
            << pMax / max(pMin, SMALL) << nl;
    }

    // ---- Rebuild local structures only when ownership changed --------------
    // ParMETIS can return the current partition (especially after a remap).
    // Avoid invalidating local state and rebuilding the local mesh in that case.
    if (nChanged > 0)
    {
        // §31: measured for the acceptance payback estimate (the DLB
        // migration timer covers migrate+updateParticleCounts only).
        const auto dsmcRebuild0 = std::chrono::steady_clock::now();
        rebuildMyCells();
        localMesh_.build(myCells_);
        dsmcRebuildLastWall_ =
            std::chrono::duration<scalar>
            (
                std::chrono::steady_clock::now() - dsmcRebuild0
            ).count();
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

// §45 P0: 碰撞相每 cell 候选数 → 窗口和累加（传入的 touched 可含重复；
// 持久触及表 cellCostQtouched_ 用掩码去重，dump 后按表归零清窗口）。
// 仅在 observe 模式且复制网格活动时生效。
void dsmcReplicatedMesh::accumulateCellCostCandidates
(
    const labelList& counts,
    const labelList& touched
)
{
    if (cellCostMode_ == "off") return;
    if (cellCostQ_.size() != mesh_.nCells()) return;
    if (cellCostTouchedMask_.size() != mesh_.nCells())
    {
        cellCostTouchedMask_.setSize(mesh_.nCells(), false);
    }
    for (label i = 0; i < touched.size(); ++i)
    {
        const label c = touched[i];
        if (c >= 0 && c < cellCostQ_.size())
        {
            cellCostQ_[c] += scalar(counts[c]);
            if (!cellCostTouchedMask_[c])
            {
                cellCostTouchedMask_[c] = true;
                cellCostQtouched_.append(c);
            }
        }
    }
}


// §45 P0: 窗口 dump——只在 dump 事件（每 dumpEvery 个 DLB 事件一次）写：
//   r<R>.P.f64   连续 vtxdist 切片的 P（globalCellParticles 快照，全场
//                复制 → 切片拼起来即完整全局数组）；
//   r<R>.Qs.bin  稀疏 Q 对 (int32 全局 cell 号, double 窗口候选和)——
//                cellCostQ_ 只在本 rank 拥有的 cell 上非零（所有者散布
//                全网格），写连续切片会得到近似全零的文件（§45 修复四）；
//   rankwork_r<R>.dat  全 rank 表：窗口净工作、等待、及每 rank 拥有
//                cell 的特征和（Np、ΣN²、owned 数、窗口候选和）——拟合
//                所需特征自包含，不再依赖连续切片与 owner 图对照。
// 写后按触及表清零窗口。中间事件只推进 §19.1 基线，不清窗口——否则
// dumpEvery>1 时累计被逐事件丢弃（§45 修复一）。
// 文件：<case>/cellCost/w<dump 序号>/。
void dsmcReplicatedMesh::dumpCellCostWindow
(
    const label currentStep,
    const labelList& globalCellParticles
)
{
    if (cellCostMode_ != "observe") return;
    const bool dumpThis = (cellCostDumpCount_ % cellCostDumpEvery_ == 0);

    if (!dumpThis) return;

    const fileName dir =
        mesh_.time().path()/"cellCost"
       /("w" + name(cellCostDumpCount_));
    mkDir(dir);

    // 二进制 P：本 rank vtxdist 段的 P 各 myN 个 double（cell 顺序 =
    // 全局编号 myStart..myEnd-1，升序）。先写后清。
    const label baseN = mesh_.nCells()/nProcs_;
    const label myStart = myRank_*baseN;
    // 尾 rank 吸收整除余数，保证 0..nCells-1 全部落在某 rank 段内
    const label myN =
        (myRank_ == nProcs_ - 1) ? (mesh_.nCells() - myStart) : baseN;

    std::ofstream posFile
    (
        (dir/("r" + name(myRank_) + ".P.f64")).c_str(),
        std::ios::out | std::ios::binary
    );
    posFile.write
    (
        reinterpret_cast<const char*>(&cellCostP_[myStart]),
        static_cast<std::streamsize>(myN)*sizeof(scalar)
    );
    posFile.close();

    // 稀疏 Q 对：触及表（已按掩码去重，即本 rank 拥有且有累加的 cell）
    {
        std::ofstream qFile
        (
            (dir/("r" + name(myRank_) + ".Qs.bin")).c_str(),
            std::ios::out | std::ios::binary
        );
        for (label i = 0; i < cellCostQtouched_.size(); ++i)
        {
            const label c = cellCostQtouched_[i];
            if (c >= 0 && c < cellCostQ_.size())
            {
                const int32_t ci = static_cast<int32_t>(c);
                const double qv = static_cast<double>(cellCostQ_[c]);
                qFile.write
                (
                    reinterpret_cast<const char*>(&ci),
                    sizeof(int32_t)
                );
                qFile.write
                (
                    reinterpret_cast<const char*>(&qv),
                    sizeof(double)
                );
            }
        }
        qFile.close();
    }

    // 每 rank 窗口净工作代理（§19.1 配对，与 feedback 同口径：
    // Allgather 累计 wall/wait 后按 rank 差分；首个 DLB 事件时
    // dsmcLastDLB* 尚未武装 → 记 -1 表示不可用）。
    // 特征列：每 rank 拥有 cell 的 Np=ΣN、SN2=ΣN²、CC=owned 数，以及
    // 本 rank 窗口候选和 Qw（触及表 cellCostQ_ 求和，含 0 值触及）。
    List<scalar> wallNow(nProcs_, 0.0), waitNow(nProcs_, 0.0);
    MPI_Allgather
    (
        &dsmcStepWallCum_, 1, MPI_DOUBLE,
        wallNow.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD
    );
    MPI_Allgather
    (
        &dsmcStepWaitCum_, 1, MPI_DOUBLE,
        waitNow.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD
    );

    double feat[4] = {0.0, 0.0, 0.0, 0.0};
    const label nCells = mesh_.nCells();
    for (label c = 0; c < nCells; ++c)
    {
        const label own = cellOwner_[c];
        if (own == myRank_)
        {
            const double n = double(globalCellParticles[c]);
            feat[0] += n;
            feat[1] += n*n;
            feat[2] += 1.0;
        }
    }
    for (label i = 0; i < cellCostQtouched_.size(); ++i)
    {
        const label c = cellCostQtouched_[i];
        if (c >= 0 && c < cellCostQ_.size())
        {
            feat[3] += static_cast<double>(cellCostQ_[c]);
        }
    }

    List<double> featAll(nProcs_*4, 0.0);
    MPI_Allgather
    (
        feat, 4, MPI_DOUBLE,
        featAll.data(), 4, MPI_DOUBLE, MPI_COMM_WORLD
    );

    OFstream workFile(dir/("rankwork_r" + name(myRank_) + ".dat"));
    workFile << "# step " << currentStep
             << " cols: r work wait Np SN2 CC Qw" << nl;
    for (label r = 0; r < nProcs_; ++r)
    {
        scalar work = -1.0, wait = -1.0;
        if (dsmcLastDLBWallCum_.size() == nProcs_)
        {
            work = max
            (
                (wallNow[r] - dsmcLastDLBWallCum_[r])
              - (waitNow[r] - dsmcLastDLBWaitCum_[r]),
                scalar(0)
            );
            wait = waitNow[r] - dsmcLastDLBWaitCum_[r];
        }
        workFile << "r" << r
                 << " " << work
                 << " " << wait
                 << " " << featAll[4*r + 0]
                 << " " << featAll[4*r + 1]
                 << " " << featAll[4*r + 2]
                 << " " << featAll[4*r + 3]
                 << nl;
    }

    // 诊断统计在清零前算，否则永远打印 touched=0 qSum=0（§45 修复三）
    if (myRank_ == 0)
    {
        scalar qSum = 0.0; label qNz = 0;
        for (label i = 0; i < cellCostQtouched_.size(); ++i)
        {
            const label c = cellCostQtouched_[i];
            if (c >= 0 && c < cellCostQ_.size())
            {
                qSum += cellCostQ_[c];
                if (cellCostQ_[c] > 0) ++qNz;
            }
        }
        Info<< "Phase C cellCost dump w" << cellCostDumpCount_
            << ": touched=" << cellCostQtouched_.size()
            << " qSum=" << qSum << " qNz=" << qNz << endl;
    }

    cellCostQtouched_.clear();
    for (label c = 0; c < cellCostTouchedMask_.size(); ++c)
    {
        cellCostTouchedMask_[c] = false;
    }
    ++cellCostDumpCount_;
}


// §44 Route A: 趋势感知 rank-k 更新器（worklog §44.3）。
// 每 TrendWindowSteps 步闭合一次窗口：收集每 rank 窗口有效工作
// （wall − 集合等待，事件步已在采样端剔除）→ 相对均值偏差 e_r →
// EMA(bar_e) → 趋势 d_r（限幅）→ q_r = Gp*bar_e + Gd*H*d_r（中心化）→
// kTarget → 阻尼更新 trendK_。重分区（owner 变更比例超限）后重置
// 历史与 trendK_，等待两个完整窗口再启用趋势项。
void dsmcReplicatedMesh::trendRankKUpdate(const label currentStep)
{
    if (trendMode_ != "rankK") return;
    if (trendK_.size() != nProcs_)
    {
        trendK_.setSize(nProcs_, cellWeightK_);
    }

    // 窗口未满：继续积累（样本在 setLastStepWallTime 逐步累加）
    if (trendSamples_ < trendWindowSteps_) return;

    // ---- 窗口闭合：收集每 rank（工作和，有效样本数） ----
    List<scalar> workNow(nProcs_, 0.0);
    List<label> sampNow(nProcs_, 0);
    MPI_Allgather
    (
        &trendWorkSum_, 1, MPI_DOUBLE,
        workNow.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD
    );
    MPI_Allgather
    (
        &trendSamples_, 1, MPI_INT,
        sampNow.data(), 1, MPI_INT, MPI_COMM_WORLD
    );

    // 每 rank 每步均值；任一 rank 整窗无有效样本 → 回退本窗（§44.2）
    List<scalar> u(nProcs_, 0.0);
    label minSamp = labelMax;
    for (label r = 0; r < nProcs_; ++r)
    {
        u[r] = sampNow[r] > 0
          ? workNow[r]/scalar(sampNow[r]) : scalar(0);
        minSamp = min(minSamp, sampNow[r]);
    }
    scalar uMean = 0.0;
    for (label r = 0; r < nProcs_; ++r) uMean += u[r];
    uMean /= scalar(nProcs_);

    label nClampK = 0;
    if (uMean <= SMALL || minSamp <= 0)
    {
        // 窗口无效：重开窗口，历史保持
        trendWorkSum_ = 0.0;
        trendSamples_ = 0;
        if (myRank_ == 0)
        {
            Info<< "Phase C DLB trend (rankK): window invalid"
                << " (minSamples=" << minSamp << ") - restart window"
                << endl;
        }
        return;
    }

    // ---- e_r / bar_e / d_r / q_r ----
    List<scalar> e(nProcs_), barE(nProcs_), q(nProcs_);
    scalar eMean = 0.0;
    for (label r = 0; r < nProcs_; ++r)
    {
        e[r] = u[r]/uMean - 1.0;
        eMean += e[r];
    }
    eMean /= scalar(nProcs_);
    for (label r = 0; r < nProcs_; ++r) e[r] -= eMean;   // 中心化

    bool trendActive = false;
    scalar dSteps = scalar(max(currentStep - trendLastWindowStep_, label(1)));
    if (trendPrevBarEValid_ && trendPrevBarE_.size() == nProcs_)
    {
        for (label r = 0; r < nProcs_; ++r)
        {
            barE[r] =
                trendTheta_*e[r] + (1.0 - trendTheta_)*trendPrevBarE_[r];
        }
        trendActive = true;   // 第二个完整窗口起启用趋势项
    }
    else
    {
        barE = e;             // 窗口 1：播种
    }
    const label kMin = max
    (
        label(1),
        mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBWorkFeedbackKMin", 8)
    );
    const label kMax = max
    (
        kMin,
        mesh_.time().controlDict().lookupOrDefault<label>
        ("replicatedMeshDLBWorkFeedbackKMax", 25)
    );
    const scalar gp = mesh_.time().controlDict().lookupOrDefault<scalar>
    (
        "replicatedMeshDLBWorkFeedbackGain",
        2.0
    );

    for (label r = 0; r < nProcs_; ++r)
    {
        q[r] = gp*barE[r];
        if (trendActive)
        {
            // 趋势项：d 限幅 ±0.5/步（实现选择，见 worklog §44）
            const scalar dRaw =
                (barE[r] - trendPrevBarE_[r])/max(dSteps, scalar(1));
            const scalar d = max(-0.5, min(0.5, dRaw));
            q[r] += trendGainGd_*trendHorizonSteps_*d;
        }
    }
    // §44.3/clean H1: store bar_e(j) as the next window's history, AFTER
    // the trend term above has consumed the previous window's value.
    // Assigning it before the dRaw loop made the trend contribution
    // identically zero (rankK silently degraded to proportional feedback).
    trendPrevBarE_ = barE;
    trendPrevBarEValid_ = true;

    scalar qMean = 0.0;
    for (label r = 0; r < nProcs_; ++r) qMean += q[r];
    qMean /= scalar(nProcs_);
    for (label r = 0; r < nProcs_; ++r) q[r] -= qMean;   // 中心化（§44.3）

    // ---- kTarget → 阻尼更新 trendK_ ----
    for (label r = 0; r < nProcs_; ++r)
    {
        const scalar kTarget =
            max(scalar(kMin), min(scalar(kMax),
                cellWeightK_*(1.0 + q[r])));
        const scalar dKraw = trendRelaxation_*(kTarget - trendK_[r]);
        const scalar dK = max(-trendMaxDeltaK_, min(trendMaxDeltaK_, dKraw));
        const scalar kNew = max(scalar(kMin), min(scalar(kMax),
            trendK_[r] + dK));
        if (kNew != trendK_[r]) ++nClampK;
        trendK_[r] = kNew;
    }

    if (myRank_ == 0)
    {
        Info<< "Phase C DLB trend (rankK): window=[" 
            << currentStep - trendSamples_ << "," << currentStep
            << "] samples=" << trendSamples_
            << " trend=" << (trendActive ? "on" : "seed")
            << " kClamped=" << nClampK << " k:";
        for (label r = 0; r < nProcs_; ++r)
        {
            Info<< " r" << r << "=" << trendK_[r];
        }
        Info<< endl;
    }

    // ---- 窗口重开 ----
    trendWorkSum_ = 0.0;
    trendSamples_ = 0;
    trendLastWindowStep_ = currentStep;
}


void dsmcReplicatedMesh::autoRebalance()
{
    if (!autoDLBEnabled_ || !active_ || nProcs_ < 2) return;

    const auto tAuto0 = std::chrono::steady_clock::now();
    ++autoRebalanceChecks_;

    const label currentStep = stepCounter_;

    // §44 Route A: 趋势窗口评估（与触发解耦，每 TrendWindowSteps 步一次；
    // 内含两次小 Allgather，所有 rank 均匀到达，集合安全）
    trendRankKUpdate(currentStep);

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
        // §17.5: per-rank load = per-step work (step wall - collective
        // wait), the same definition as "Phase C per-rank step work".
        // The old metric (move+coll wall delta since the DLB snapshot)
        // excluded build/post/comm/wait, so it could not see either the
        // dump-cell tax on the gate rank or the waiting that equalises
        // per-rank walls.  Value lags one step, like productiveTime_.
        const scalar myStepWork =
            max(lastStepWallTime_ - lastStepWait_, scalar(0));
        scalarList snapWork(nProcs_, 0.0);
        const auto tPostSnapshotComm0 = std::chrono::steady_clock::now();
        MPI_Allgather(&myStepWork, 1, MPI_DOUBLE,
                      snapWork.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);

        // comm_further_opt.md §15.4: per-rank instantaneous load table
        // (particles absolute; moveParcels/trackCalls/collisions/
        // candidates as window deltas since the DLB snapshot).  Zero
        // extra collectives: particles come from allParticleCounts_
        // (already gathered by updateParticleCounts), the counters are
        // cumulative label values gathered once here.
        double myLabels[5] =
        {
            double(allParticleCounts_[myRank_]),
            cloud_.moveDetailCumulative(0),
            cloud_.moveDetailCumulative(1),
            cloud_.cumulativeCollisions(),
            cloud_.cumulativeCollisionCandidates()
        };
        List<scalar> snapLabels(5*nProcs_, scalar(0));
        MPI_Allgather
        (
            myLabels,
            5,
            MPI_DOUBLE,
            snapLabels.data(),
            5,
            MPI_DOUBLE,
            MPI_COMM_WORLD
        );
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
            Info<< "  [post-DLB step " << step
                << "] per-rank step work (wall-wait):";
            for (label i = 0; i < nProcs_; ++i)
            {
                const scalar load = snapWork[i];
                Info<< " r" << i << "=" << load;
                if (load > maxLoad) maxLoad = load;
                if (load < minLoad) minLoad = load;
            }
            Info<< " | max/min=" << maxLoad / max(minLoad, SMALL) << nl;

            // Per-rank load detail: particles (absolute) + window deltas
            // of moveParcels/trackCalls/collisions/candidates vs the DLB
            // snapshot — the direct read on whether the rebalance
            // actually equalised the particle distribution.
            scalar maxPart = 0, minPart = GREAT;
            Info<< "  [post-DLB step " << step
                << "] per-rank particles | moveParcels trackCalls"
                << " collisions candidates (delta):" << nl;
            for (label i = 0; i < nProcs_; ++i)
            {
                const scalar part = snapLabels[5*i + 0];
                maxPart = max(maxPart, part);
                minPart = min(minPart, part);
                Info<< "    r" << i << " " << part;
                if (dsmcLastDLBLabels_.size() == 5*nProcs_)
                {
                    for (label k = 1; k < 5; ++k)
                    {
                        Info<< " "
                            << snapLabels[5*i + k] - dsmcLastDLBLabels_[5*i + k];
                    }
                }
                Info<< nl;
            }
            Info<< "  [post-DLB step " << step
                << "] particles max/min = "
                << scalar(maxPart)/scalar(max(minPart, 1)) << nl;
        }
        --postDLBSnapshotCountdown_;
    }

    // ---- Accumulate productive time ----
    // Legacy metric: CPU time of move+collision minus migration wall.
    // Wall metric (§5.6): per-iteration full wall reported by the solver
    // main loop minus that iteration's collective waits.  §17.5 fix: the
    // wait delta is now captured in setLastStepWallTime (end of each
    // iteration, after every collective) instead of here, so wall(N) and
    // wait(N) pair exactly — a DLB's migration storm no longer leaks into
    // the next step's wait and shows up as phantom work / clamped-to-0.
    // Both inputs still lag one step here: the pair refers to iteration
    // N-1 while this code runs during iteration N.
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
    productiveTime_ += stepProductive;
    // §19.1: cumulative step wall / collective wait for the Inter-DLB
    // actual-load table (wall - comm wait).  The pair is the same
    // (N-1)-lagged one as stepProductive, so each step is counted once.
    dsmcStepWallCum_ += lastStepWallTime_;
    dsmcStepWaitCum_ += lastStepWait_;
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
    // §43: particleMaxMin diagnostic.  The measured per-rank particle
    // max/min is computed unconditionally (no extra collective: the
    // counts are already gathered by updateParticleCounts), so the
    // logged "particleMaxMin=" is a real measurement even when the
    // particle gate is disabled — previously it stayed at the 1.0
    // placeholder and misled tuning (worklog §37/§40/§41, §1.4).
    // The gate condition itself is byte-for-byte unchanged: when the
    // gate is enabled it compares the same ratio against the same
    // threshold, and a missing counts table still allows the check.
    bool particleCountsAvailable = (allParticleCounts_.size() == nProcs_);
    scalar particleMaxMinMeas = -1.0;    // -1 = counts not available yet
    if (particleCountsAvailable)
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
            particleMaxMinMeas =
                scalar(maxParticles) / max(scalar(minParticles), SMALL);
        }
    }
    if
    (
        particleGateEnabled
     && !forcedTriggered
     && particleCountsAvailable
     && particleMaxMinMeas > 0
    )
    {
        particleGateImbalance = particleMaxMinMeas;
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
                        << ": particle max/min=" << particleMaxMinMeas
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
                            << " particleMaxMin=" << particleMaxMinMeas
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
                    << ": particle max/min=" << particleMaxMinMeas
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
                        << " particleMaxMin=" << particleMaxMinMeas
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
    const label nChanged = reassignByParMetisAdaptiveRepart(checkImbalance);

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

        // comm_further_opt.md §7: zero-lag post-migration per-rank
        // particle table.  Verifies whether ParMETIS's partition actually
        // balances the particle distribution AT THE INSTANT of
        // rebalancing (before any drift).  allParticleCounts_ was just
        // gathered by updateParticleCounts() — no extra collective.
        if (dlbProfile && myRank_ == 0)
        {
            label immMax = 0, immMin = labelMax;
            for (label r = 0; r < nProcs_; ++r)
            {
                const label cnt = allParticleCounts_[r];
                if (cnt > immMax) immMax = cnt;
                if (cnt < immMin) immMin = cnt;
            }
            Info<< "  [post-rebalance zero-lag] per-rank particles"
                << " max/min = "
                << scalar(immMax)/scalar(max(immMin, 1))
                << " (max=" << immMax << " min=" << immMin << ")" << nl;
            Info<< "    ";
            for (label r = 0; r < nProcs_; ++r)
            {
                Info<< " r" << r << "=" << allParticleCounts_[r];
            }
            Info<< nl;
        }
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

        // §44 Route A: 重分区执行步的样本含迁移风暴 → 剔除；
        // owner 变更比例超限时重置趋势历史与累积 k（§44.3），
        // 小修正（≤ 阈值）保留历史继续积累。
        trendSkipSample_ = true;
        if (trendMode_ == "rankK")
        {
            const scalar cr =
                scalar(nChanged)/max(scalar(mesh_.nCells()), scalar(1));
            if (cr > trendResetChangedRatio_)
            {
                if (trendK_.size() != nProcs_) trendK_.setSize(nProcs_, cellWeightK_);
                for (label r = 0; r < nProcs_; ++r) trendK_[r] = cellWeightK_;
                trendPrevBarEValid_ = false;
            }
        }
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
        scalarList allWallCum, allWaitCum;
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
            // §19.1: per-rank cumulative step wall / comm wait for the
            // Inter-DLB actual-load table (wall - comm wait).
            double myWallCum = dsmcStepWallCum_;
            double myWaitCum = dsmcStepWaitCum_;
            allWallCum.setSize(nProcs_, 0.0);
            allWaitCum.setSize(nProcs_, 0.0);
            MPI_Allgather(&myWallCum, 1, MPI_DOUBLE,
                          allWallCum.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgather(&myWaitCum, 1, MPI_DOUBLE,
                          allWaitCum.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
            addAutoRebalanceCommTime
            (
                std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now() - tPostDiagComm0
                ).count()
            );
        }
        if (dlbProfile && dsmcLastDLBWallCum_.size() == nProcs_)
        {
            // §19.1: redefined as wall - collective wait over the window
            // since the last rebalance — the same definition as "Phase C
            // per-rank step work", replacing the move+coll wall deltas
            // (which could not see build/post or the comm-wait that
            // equalises per-rank walls).
            Info<< "  Inter-DLB actual load (since last rebalance,"
                << " wall - comm wait):" << nl;
            scalar maxA = 0.0, minA = GREAT;
            for (label i = 0; i < nProcs_; ++i)
            {
                const scalar dW = allWallCum[i] - dsmcLastDLBWallCum_[i];
                const scalar dC = allWaitCum[i] - dsmcLastDLBWaitCum_[i];
                const scalar work = dW - dC;
                Info<< "    rank" << i << ": wall=" << dW
                    << " commWait=" << dC << " work=" << work << nl;
                if (work > maxA) maxA = work;
                if (work < minA) minA = work;
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
            // §19.1: wall/wait cumulative baselines for the Inter-DLB
            // actual-load table (wall - comm wait).
            dsmcLastDLBWallCum_ = allWallCum;
            dsmcLastDLBWaitCum_ = allWaitCum;
            postDLBSnapshotCountdown_ = mesh_.time().controlDict().lookupOrDefault<label>
                ("replicatedMeshDLBProfileSteps", 5);

            // §15.4: label snapshot (particles + cumulative counters) at
            // DLB time — the baseline for the post-DLB window deltas.
            double myLabels[5] =
            {
                double(allParticleCounts_[myRank_]),
                cloud_.moveDetailCumulative(0),
                cloud_.moveDetailCumulative(1),
                cloud_.cumulativeCollisions(),
                cloud_.cumulativeCollisionCandidates()
            };
            dsmcLastDLBLabels_.setSize(5*nProcs_, 0.0);
            MPI_Allgather
            (
                myLabels,
                5,
                MPI_DOUBLE,
                dsmcLastDLBLabels_.data(),
                5,
                MPI_DOUBLE,
                MPI_COMM_WORLD
            );
        }
        else
        {
            lastDLBMoveT_.clear();
            lastDLBCollT_.clear();
            dsmcLastDLBLabels_.clear();
            dsmcLastDLBWallCum_.clear();
            dsmcLastDLBWaitCum_.clear();
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
        // §36.2: bounded, persistent migration buffers.  DynamicList
        // keeps its storage when setSize() is called with a size <= its
        // §36.2: bounded, persistent migration scratch buffers.  DynamicList
        // keeps its storage when setSize() is called with a size <= its
        // capacity (DynamicListI.H:186 — unlike List, which always
        // reallocates), so persistent buffers remove the per-step
        // allocate/free churn: the per-destination send lists and the two
        // receive chunk buffers.  Measured on the 1200w case the traffic is
        // ~10-25 MB per call per rank (sendMB/recvMB over 502 calls), so the
        // retained scratch is ~50-100 MB per rank.  Capacity exceeding
        // replicatedMeshMigBufCapMB is released after the transfer, as a
        // guard against a DLB storm pinning memory.  `kept` stays per-step
        // because its storage is handed to the cloud by
        // takeMoveOrderedParcels().
        // §36.2 A/B: with buffers disabled, release all retained capacity
        // so the setSize() calls below reallocate exactly as the pre-§36
        // per-step locals did.
        if (!migReuseMigBuffers_)
        {
            dsmcMigSendBuf_.clear();      // List::clear() frees the array
            dsmcMigRecvChunks_.clear();
            dsmcMigDeleteBuf_.clearStorage();
        }
        if (dsmcMigSendBuf_.size() != nProcs_)
        {
            dsmcMigSendBuf_.setSize(nProcs_);
        }
        List<DynamicList<dsmcParcel::TransferData>>& sendTD = dsmcMigSendBuf_;

        // §36.2: reset the addressed size (capacity is kept).  The serial
        // OMP=1 classification path only appends to sendTD[dst], so without
        // this the persistent buffers would carry last step's records
        // forward and re-send them (duplicating parcels).  The pre-§36 code
        // got this for free from the per-step local construction.
        forAll(sendTD, dst)
        {
            sendTD[dst].clear();
        }

        if (dsmcMigRecvChunks_.size() != 2)
        {
            dsmcMigRecvChunks_.setSize(2);
        }
        List<DynamicList<dsmcParcel::TransferData>>& recvChunks =
            dsmcMigRecvChunks_;

        dsmcMigDeleteBuf_.clear();
        DynamicList<dsmcParcel*>& toDelete = dsmcMigDeleteBuf_;
        if (!migReuseMigBuffers_)
        {
            // pre-§36 initial sizing of the delete scratch.  The original
            // was DynamicList(n), which RESERVES n slots and leaves the
            // addressed size at 0 (DynamicListI.H:36-46); setSize() here
            // would expose n uninitialised pointers to the append path.
            toDelete.reserve(cloud_.size() / 4);
        }
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

        // §35 P0: classify/count vs pack split (set in each branch).
        std::chrono::steady_clock::time_point tClassifyEnd = tStart;

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
                // §34.1: skip the self entry.  destinationCounts[myRank_]
                // is non-zero (the classification counts every parcel) but
                // the pack loop skips self-destined parcels, so
                // sendTD[myRank_] would be allocated and never written —
                // in a normal step most parcels stay local, so this array
                // can far exceed the real send volume.  The async path
                // already excludes self.
                for (label dst = 0; dst < nProcs_; ++dst)
                {
                    if (dst != myRank_)
                    {
                        sendTD[dst].setSize(destinationCounts[dst]);
                    }
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

                // §35 P0: classify/count ends here, the pack region begins.
                tClassifyEnd = std::chrono::steady_clock::now();
                dsmcMigClassifyWall_ +=
                    std::chrono::duration<scalar>
                    (tClassifyEnd - tStart).count();

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
            // §35 P0: the fallback path classifies and packs in one loop.
            tClassifyEnd = std::chrono::steady_clock::now();
            dsmcMigClassifyWall_ +=
                std::chrono::duration<scalar>(tClassifyEnd - tStart).count();
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
        // §35 P0: pack stage (packTransfer + kept/toDelete fill).
        dsmcMigPackWall_ +=
            std::chrono::duration<scalar>(tPackEnd - tClassifyEnd).count();

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

        // §35 P0: per-step transfer volume / active peers (sizes are final
        // after the exchange; every rank knows its own row and column).
        double stepSendBytes = 0.0;
        double stepRecvBytes = 0.0;
        label stepPeers = 0;
        {
            for (label i = 0; i < nProcs_; ++i)
            {
                stepSendBytes += double(sendSizes[i]);
                stepRecvBytes += double(recvSizes[i]);
                if
                (
                    i != myRank_
                 && (sendSizes[i] != 0 || recvSizes[i] != 0)
                )
                {
                    ++stepPeers;
                }
            }
            dsmcMigSendBytes_ += stepSendBytes;
            dsmcMigRecvBytes_ += stepRecvBytes;
            dsmcMigPeerCount_ += double(stepPeers);
        }

        const label maxChunkRecords =
            max
            (
                label(1),
                transferChunkBytes_
              / label(sizeof(dsmcParcel::TransferData))
            );
        scalar chunkWaitWall = 0.0;
        scalar chunkDeserializeWall = 0.0;
        scalar chunkUnpackWall = 0.0;
        scalar chunkAppendWall = 0.0;
        label stepChunks = 0;   // §35 P0
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
                ++stepChunks;
                ++dsmcMigChunkCount_;   // §35 P0
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
                // §24.6 probe: split the deser interval — the parallel
                // unpack allocates through the M6 global arena spinlock
                // (design assumed ~1% churn spread over a step; a
                // migration is a 48-thread allocation burst), the append
                // loop is serial (addParticle + kept growth).
                const auto tChunkUnpackEnd = std::chrono::steady_clock::now();
                chunkUnpackWall += std::chrono::duration<scalar>
                (
                    tChunkUnpackEnd - tChunkDeserializeStart
                ).count();
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
                chunkAppendWall += std::chrono::duration<scalar>
                (
                    std::chrono::steady_clock::now() - tChunkUnpackEnd
                ).count();
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

        // §35 P0: per-step sub-items (unpack/append ⊂ deser) and parcels.
        dsmcMigUnpackWall_ += double(chunkUnpackWall);
        dsmcMigAppendWall_ += double(chunkAppendWall);
        dsmcMigParcelsIn_ += double(nRecv);

        // §34.5: post the particle-count gather here (all received parcels
        // are appended, so cloud_.size() is final) so its convergence
        // overlaps the ordered-view rebuild below.  The cloud's
        // updateParticleCounts() completes it.  (The async-finish path
        // leaves the flag unset and updateParticleCounts falls back to
        // the blocking gather — no overlap, no hazard.)
        postParticleCountsGather();

        // §35 P0: moveOrdered rebuild / clear stage.
        const auto tRebuild0 = std::chrono::steady_clock::now();
        if (maintainMoveOrdered && kept.size() == cloud_.size())
        {
            if (migTakeOrderedParcels_)
            {
                // §36.1: zero-copy hand-off — the cloud steals this list's
                // storage instead of re-walking it into a fresh array.
                cloud_.takeMoveOrderedParcels(kept);
            }
            else
            {
                cloud_.setMoveOrderedParcels(kept);
            }
        }
        else
        {
            cloud_.clearMoveOrderedParcels();
        }
        dsmcMigRebuildWall_ +=
            std::chrono::duration<scalar>
            (std::chrono::steady_clock::now() - tRebuild0).count();
        const auto tPostEnd = std::chrono::steady_clock::now();

        // §36.2: release persistent buffer capacity that grew past the
        // configured cap (a one-off DLB storm must not pin hundreds of MB
        // per rank for the rest of the run).  capacity() reads only — no
        // allocation, so this is off the allocation budget.
        {
            const scalar capMB =
                mesh_.time().controlDict().lookupOrDefault<scalar>
                (
                    "replicatedMeshMigBufCapMB", 1024.0
                );
            const double capBytes = double(capMB)*1024.0*1024.0;
            const double tdBytes = double(sizeof(dsmcParcel::TransferData));
            const double ptrBytes = double(sizeof(dsmcParcel*));

            forAll(sendTD, i)
            {
                if (double(sendTD[i].capacity())*tdBytes > capBytes)
                {
                    sendTD[i].clearStorage();
                }
            }
            forAll(recvChunks, j)
            {
                if (double(recvChunks[j].capacity())*tdBytes > capBytes)
                {
                    recvChunks[j].clearStorage();
                }
            }
            if (double(toDelete.capacity())*ptrBytes > capBytes)
            {
                toDelete.clearStorage();
            }
        }

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
        dsmcMigTotalWall_ +=
            std::chrono::duration<scalar>(tEnd - tStart).count();   // §35 P0
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
            // §24 probe: split the "sizeX" interval into the exchange
            // itself and the arrival skew.  If the exchange is cheap and
            // the skew carries most of sizeX, the comm cost is BSP
            // synchronisation of step-work imbalance, not data volume.
            const scalar myPrep =
                std::chrono::duration<scalar>(tLocalPrepEnd - tStart).count();
            double prepMaxD = double(myPrep);
            MPI_Allreduce
            (
                MPI_IN_PLACE, &prepMaxD, 1, MPI_DOUBLE,
                MPI_MAX, MPI_COMM_WORLD
            );
            const scalar prepMax = scalar(prepMaxD);
            Info<< "Replicated mesh migration[" << migrationCalls_ << "]: "
                << "rank " << myRank_
                << " sent " << nMigratedOut << " parcels, received " << nRecv
                << " (global " << globalMigrated << ")"
                << " wall " << std::chrono::duration<scalar>(tEnd - tStart).count() << "s"
                << " [classify " << std::chrono::duration<scalar>
                    (tClassifyEnd - tStart).count() << "s"
                << " pack " << std::chrono::duration<scalar>
                    (tPackEnd - tClassifyEnd).count() << "s"
                << " rebuild " << std::chrono::duration<scalar>
                    (std::chrono::steady_clock::now() - tRebuild0).count() << "s"
                << " upcPost+ " << (dsmcMigUpcPostWall_ - dsmcMigUpcPostPrev_)
                << "s upcWait+ " << (dsmcMigUpcWaitWall_ - dsmcMigUpcWaitPrev_)
                << "s peers " << stepPeers
                << " chunks " << stepChunks
                << " sendB " << stepSendBytes
                << " recvB " << stepRecvBytes
                << " prep " << std::chrono::duration<scalar>
                    (tLocalPrepEnd - tPackEnd).count() << "s"
                << " prepMax " << prepMax << "s"
                << " sizeX " << std::chrono::duration<scalar>
                    (tSizeExchangeEnd - tLocalPrepEnd).count() << "s"
                << " skew "
                << (prepMax - std::chrono::duration<scalar>
                    (tLocalPrepEnd - tStart).count()) << "s"
                << " postReq " << std::chrono::duration<scalar>
                    (tRequestPostEnd - tSizeExchangeEnd).count() << "s"
                << " wait " << chunkWaitWall << "s"
                << " deser " << chunkDeserializeWall << "s"
                << " unpack " << chunkUnpackWall << "s"
                << " append " << chunkAppendWall << "s"
                << " post " << std::chrono::duration<scalar>
                    (tPostEnd - tTransferEnd).count() << "s"
                << " candX " << std::chrono::duration<scalar>
                    (tEnd - tPostEnd).count() << "s]"
                << endl;
            // §35 P0: remember the printed upc totals for the next delta.
            dsmcMigUpcPostPrev_ = dsmcMigUpcPostWall_;
            dsmcMigUpcWaitPrev_ = dsmcMigUpcWaitWall_;
        }
        // dlb_offload.md P0: per-step migration stamps on the common
        // cloud probe epoch (sync flat path).  Stored for the cloud's
        // per-rank probe file — per-rank Info does not reach the job
        // stdout in replicated-mesh mode, only the cloud master's does.
        if (cloud_.offloadProbeEnabled())
        {
            const double probeEpoch = cloud_.offloadProbeEpoch();
            offPrep1_ = std::chrono::duration<double>
                (tLocalPrepEnd.time_since_epoch()).count() - probeEpoch;
            offSx1_ = std::chrono::duration<double>
                (tSizeExchangeEnd.time_since_epoch()).count() - probeEpoch;
            offMig1_ = std::chrono::duration<double>
                (tEnd.time_since_epoch()).count() - probeEpoch;
            offSent_ = nMigratedOut;
            offRecv_ = nRecv;
            offSentB_ = stepSendBytes;
            offRecvB_ = stepRecvBytes;
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

    // clean N1: count records per destination first, then write the count
    // header ahead of the framed records.
    labelList nSendParcels(nProcs_, 0);
    if (useOrderedTraversal)
    {
        const auto& ordered = cloud_.moveOrderedParcels();
        forAll(ordered, i)
        {
            const label dstRank = cellOwner_[ordered[i]->cell()];
            if (dstRank != myRank_) ++nSendParcels[dstRank];
        }
    }
    else
    {
        forAllConstIter(Cloud<dsmcParcel>, cloud_, iter)
        {
            const label dstRank = cellOwner_[iter().cell()];
            if (dstRank != myRank_) ++nSendParcels[dstRank];
        }
    }
    for (label i = 0; i < nProcs_; ++i)
    {
        if (i != myRank_)
        {
            sendStreams[i] << nSendParcels[i] << token::SPACE;
        }
    }

    if (useOrderedTraversal)
    {
        const auto& ordered = cloud_.moveOrderedParcels();
        for (label i = 0; i < ordered.size(); ++i)
        {
            dsmcParcel& p = *ordered[i];
            const label dstRank = cellOwner_[p.cell()];
            if (dstRank != myRank_)
            {
                dsmcNFSentParcels_ += 1.0;
                dsmcNFSendBytes_ += double
                (
                    dsmcWriteNonFlatRecord(sendStreams[dstRank], p)
                );
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
                dsmcNFSentParcels_ += 1.0;
                dsmcNFSendBytes_ += double
                (
                    dsmcWriteNonFlatRecord(sendStreams[dstRank], p)
                );
                toDelete.append(&p);
                ++nMigratedOut;
            }
            else
            {
                kept.append(&p);
            }
        }
    }
    ++dsmcNFCalls_;
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
            const label nRecvParcels = dsmcReadNonFlatCount(is, recvSize);
            for (label k = 0; k < nRecvParcels; ++k)
            {
                auto* newp = dsmcReadNonFlatRecord
                (
                    mesh_, is, recvSize, k, nRecvParcels
                );
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
                const label nRecvParcels =
                    dsmcReadNonFlatCount(is, actualBytes);
                for (label k = 0; k < nRecvParcels; ++k)
                {
                    auto* newp = dsmcReadNonFlatRecord
                    (
                        mesh_, is, actualBytes, k, nRecvParcels
                    );
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
                const label nRecvParcels =
                    dsmcReadNonFlatCount(is, recvSizes[i]);
                for (label k = 0; k < nRecvParcels; ++k)
                {
                    auto* newp = dsmcReadNonFlatRecord
                    (
                        mesh_, is, recvSizes[i], k, nRecvParcels
                    );
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

// ============================================================================
// comm_further_opt.md §2.2 — flat-path async migration overlap
// ============================================================================

bool dsmcReplicatedMesh::asyncOverlapEnabled() const
{
    if (!dsmcAsyncOverlapRead_)
    {
        dsmcAsyncOverlapEnabledCache_ =
            mesh_.time().controlDict().lookupOrDefault<bool>
            (
                "replicatedMeshAsyncMigration",
                false
            );
        dsmcAsyncOverlapRead_ = true;
        if (dsmcAsyncOverlapEnabledCache_ && myRank_ == 0)
        {
            Info<< "Replicated mesh: async migration overlap enabled"
                << " (begin / gap-build / finish)" << endl;
        }
    }
    return dsmcAsyncOverlapEnabledCache_;
}

void dsmcReplicatedMesh::markAsyncOverlapGapBuild()
{
    dsmcFlatAsync_.gapBuildDone = true;
}

// flatBegin — pack + delete migrated-out + size exchange + post chunked
// Isends.  Everything after this point is overlappable with the caller's
// kept-parcel occupancy build; state lives in dsmcFlatAsync_.
// Mirrors the sync flat path in migrateParticlesByCellOwner(); that path
// is intentionally left untouched as the fallback.
void dsmcReplicatedMesh::flatMigrateBegin()
{
    dsmcFlatAsyncState& st = dsmcFlatAsync_;
    const auto tStart = std::chrono::steady_clock::now();
    st.tStart = tStart;
    st.active = true;
    st.gapBuildDone = false;
    st.nMigratedOut = 0;
    st.sendReqs.clear();
    st.sendTD.setSize(nProcs_);
    st.kept.clear();

    DynamicList<dsmcParcel*> toDelete(cloud_.size() / 4);
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
                }
            }
            #pragma omp barrier

            // Serial prepass: per-destination write offsets, per-thread
            // kept/delete base positions, exact sendTD sizes.
            {
                // Total per-destination counts across threads
                labelList destTotal(nProcs_, 0);
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    for (label dst = 0; dst < nProcs_; ++dst)
                    {
                        destTotal[dst] +=
                            threadDestinationCounts[threadI][dst];
                    }
                }
                label totalDelete = 0;
                for (label dst = 0; dst < nProcs_; ++dst)
                {
                    if (dst != myRank_)
                    {
                        totalDelete += destTotal[dst];
                        st.sendTD[dst].setSize(destTotal[dst]);
                    }
                }
                toDelete.setSize(totalDelete);
                st.kept.setSize(ordered.size() - totalDelete);

                // Per-destination write cursors per thread (concatenate
                // each thread's counts per destination)
                labelList destCursor(nProcs_, 0);
                for (label threadI = 0; threadI < moveThreads; ++threadI)
                {
                    for (label dst = 0; dst < nProcs_; ++dst)
                    {
                        if (dst == myRank_) continue;
                        threadDestinationOffsets[threadI][dst] =
                            destCursor[dst];
                        destCursor[dst] +=
                            threadDestinationCounts[threadI][dst];
                    }
                }

                // Per-thread kept/delete base positions
                labelList threadKeptOffsets(moveThreads, 0);
                labelList threadDeleteOffsets(moveThreads, 0);
                {
                    label keptBase = 0;
                    label deleteBase = 0;
                    for (label threadI = 0; threadI < moveThreads; ++threadI)
                    {
                        label threadKept = 0;
                        label threadDelete = 0;
                        for (label dst = 0; dst < nProcs_; ++dst)
                        {
                            const label cnt =
                                threadDestinationCounts[threadI][dst];
                            if (dst == myRank_)
                            {
                                threadKept += cnt;
                            }
                            else
                            {
                                threadDelete += cnt;
                            }
                        }
                        threadKeptOffsets[threadI] = keptBase;
                        threadDeleteOffsets[threadI] = deleteBase;
                        keptBase += threadKept;
                        deleteBase += threadDelete;
                    }
                }

                #pragma omp parallel num_threads(moveThreads)
                {
                    const label threadI = omp_get_thread_num();
                    const label begin2 =
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
                    const label end2 =
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

                    for (label i = begin2; i < end2; ++i)
                    {
                        dsmcParcel* pPtr = ordered[i];
                        const label dstRank = cellOwner_[pPtr->cell()];
                        if (dstRank != myRank_)
                        {
                            pPtr->packTransfer
                            (
                                st.sendTD[dstRank]
                                        [destinationOffsets[dstRank]++]
                            );
                            toDelete[deleteI++] = pPtr;
                        }
                        else
                        {
                            st.kept[keptI++] = pPtr;
                        }
                    }
                }
            }
            nMigratedOut = toDelete.size();
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
                    st.sendTD[dstRank].append
                    (
                        dsmcParcel::TransferData()
                    );
                    p.packTransfer(st.sendTD[dstRank].last());
                    toDelete.append(ordered[i]);
                    ++nMigratedOut;
                }
                else
                {
                    st.kept.append(ordered[i]);
                }
            }
        }
    }
    else
    {
        forAllIter(dsmcCloud, cloud_, iter)
        {
            dsmcParcel& p = iter();
            const label dstRank = cellOwner_[p.cell()];
            if (dstRank != myRank_)
            {
                st.sendTD[dstRank].append(dsmcParcel::TransferData());
                p.packTransfer(st.sendTD[dstRank].last());
                toDelete.append(&p);
                ++nMigratedOut;
            }
            else
            {
                st.kept.append(&p);
            }
        }
    }
    st.nMigratedOut = nMigratedOut;
    const auto tPackEnd = std::chrono::steady_clock::now();
    st.tPackEnd = tPackEnd;

    List<std::uint64_t> sendSizes
    (
        nProcs_,
        static_cast<std::uint64_t>(0)
    );
    for (label i = 0; i < nProcs_; ++i)
    {
        sendSizes[i] =
            static_cast<std::uint64_t>(st.sendTD[i].size())
          * static_cast<std::uint64_t>
            (sizeof(dsmcParcel::TransferData));
    }

    st.recvSizes.setSize
    (
        nProcs_,
        static_cast<std::uint64_t>(0)
    );

    // Size exchange (synchronous collective — all ranks must arrive)
    MPI_Alltoall
    (
        sendSizes.data(),
        1,
        MPI_UINT64_T,
        st.recvSizes.data(),
        1,
        MPI_UINT64_T,
        MPI_COMM_WORLD
    );

    // Delete migrated-out parcels AFTER packing (transfer data copied)
    forAll(toDelete, i) { cloud_.deleteParcel(toDelete[i]); }
    st.tLocalPrepEnd = std::chrono::steady_clock::now();
    st.tSizeExEnd = st.tLocalPrepEnd;

    // Post chunked Isends; sendTD must stay alive until finish's Waitall
    const label maxChunkRecords =
        max
        (
            label(1),
            transferChunkBytes_
          / label(sizeof(dsmcParcel::TransferData))
        );
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
                st.sendTD[dst].data() + label(offset),
                int(nChunk*sizeof(dsmcParcel::TransferData)),
                MPI_BYTE,
                dst,
                101,
                MPI_COMM_WORLD,
                &req
            );
            st.sendReqs.append(req);
            offset += nChunk;
        }
    }
    st.tReqPostEnd = std::chrono::steady_clock::now();

    // §35 P2: async-path stage accounting (separate from the mixed
    // trigger accumulators above).
    {
        double sendBytes = 0.0;
        forAll(st.sendTD, peer)
        {
            sendBytes +=
                double(st.sendTD[peer].size())
              * double(sizeof(dsmcParcel::TransferData));
        }
        dsmcAsyncClassifyPackWall_ +=
            std::chrono::duration<double>(st.tPackEnd - st.tStart).count();
        dsmcAsyncPrepWall_ +=
            std::chrono::duration<double>
            (st.tLocalPrepEnd - st.tPackEnd).count();
        dsmcAsyncPostReqWall_ +=
            std::chrono::duration<double>
            (st.tReqPostEnd - st.tLocalPrepEnd).count();
        dsmcAsyncParcelsOut_ += double(nMigratedOut);
        dsmcAsyncSendBytes_ += sendBytes;
    }

    // moveOrdered = kept (received appended in finish).  cloud_.size()
    // == kept now, so the gap buildCellOccupancy() ordered path is valid.
    if (maintainMoveOrdered && st.kept.size() == cloud_.size())
    {
        cloud_.setMoveOrderedParcels(st.kept);
    }
    else
    {
        cloud_.clearMoveOrderedParcels();
    }

    asyncMigrationPending_ = true;
    totalParcelsMigrated_ += nMigratedOut;
    ++migrationCalls_;
}

// flatMigrateFinish — receive pipeline + deserialize + hand received to
// the cloud for incremental occupancy + release send buffers.
void dsmcReplicatedMesh::flatMigrateFinish()
{
    dsmcFlatAsyncState& st = dsmcFlatAsync_;
    if (!st.active) return;

    const bool maintainMoveOrdered =
        cloud_.openmpEnabled()
     && cloud_.openmpMoveEnabled()
     && cloud_.ompNumThreads() > 1;

    const std::uint64_t recordBytes =
        static_cast<std::uint64_t>
        (sizeof(dsmcParcel::TransferData));
    const label maxChunkRecords =
        max
        (
            label(1),
            transferChunkBytes_
          / label(sizeof(dsmcParcel::TransferData))
        );
    std::uint64_t maxIncomingRecords = 0;
    for (label src = 0; src < nProcs_; ++src)
    {
        if (src == myRank_ || st.recvSizes[src] == 0) continue;
        if (st.recvSizes[src] % recordBytes != 0)
        {
            FatalErrorInFunction
                << "Received byte count " << scalar(st.recvSizes[src])
                << " from rank " << src
                << " is not a multiple of TransferData size "
                << scalar(recordBytes) << abort(FatalError);
        }
        const std::uint64_t nRecords = st.recvSizes[src] / recordBytes;
        if (nRecords > maxIncomingRecords)
        {
            maxIncomingRecords = nRecords;
        }
    }

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
    MPI_Request recvReqs[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    std::uint64_t postedRecords[2] = {0, 0};

    scalar chunkWaitWall = 0.0;
    scalar chunkDeserializeWall = 0.0;
    label nRecv = 0;
    DynamicList<dsmcParcel*> receivedAll;

    for (label src = 0; src < nProcs_; ++src)
    {
        if (src == myRank_ || st.recvSizes[src] == 0) continue;

        const std::uint64_t nRecords = st.recvSizes[src] / recordBytes;
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
                // Base-class raw append: the post-transfer code below
                // re-establishes the ordered view / incremental occupancy.
                cloud_.Cloud<dsmcParcel>::addParticle(received[j]);
                receivedAll.append(received[j]);
                st.kept.append(received[j]);
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

    // Sends: wait then release the packed buffers (M6-fix semantics)
    if (st.sendReqs.size() > 0)
    {
        const auto tSendWaitStart = std::chrono::steady_clock::now();
        MPI_Waitall
        (
            st.sendReqs.size(),
            st.sendReqs.data(),
            MPI_STATUSES_IGNORE
        );
        chunkWaitWall += std::chrono::duration<scalar>
        (
            std::chrono::steady_clock::now() - tSendWaitStart
        ).count();
    }
    forAll(st.sendTD, dst)
    {
        st.sendTD[dst].clear();
        st.sendTD[dst].shrink();
    }

    const auto tTransferEnd = std::chrono::steady_clock::now();

    // moveOrdered = kept + received
    if (maintainMoveOrdered && st.kept.size() == cloud_.size())
    {
        cloud_.setMoveOrderedParcels(st.kept);
    }
    else
    {
        cloud_.clearMoveOrderedParcels();
    }

    // Incremental occupancy update: the caller built occupancy over kept
    // in the gap; append the received parcels into it.  The entry check
    // uses gapBuildDone (NOT occupancyOrderedParcelsValid_, which the
    // moveOrdered set/clear above just invalidated - the ordered DATA is
    // still the gap build's product).  Falls back to a full rebuild when
    // the received fraction is large.
    const auto tOccAppend0 = std::chrono::steady_clock::now();
    const double occValidate0 = cloud_.occupancyValidationWall();
    const bool incremental =
        st.gapBuildDone
     && cloud_.appendReceivedToOccupancy(receivedAll);
    dsmcMigOccAppendWall_ += std::chrono::duration<double>
    (
        std::chrono::steady_clock::now() - tOccAppend0
    ).count();
    dsmcMigOccValidateWall_ +=
        cloud_.occupancyValidationWall() - occValidate0;

    // Exchange candidate counts for optional offload diagnostics.
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
    {
        // §35 P2: async-path stages (finish covers wait+unpack+append+
        // occupancy rebuild+post as one parent).
        double recvBytes = 0.0;
        forAll(st.recvSizes, peer)
        {
            recvBytes += double(st.recvSizes[peer]);
        }
        dsmcAsyncFinishWall_ +=
            std::chrono::duration<double>(tEnd - st.tReqPostEnd).count();
        dsmcAsyncTotalWall_ +=
            std::chrono::duration<double>(tEnd - st.tStart).count();
        dsmcAsyncParcelsIn_ += double(nRecv);
        dsmcAsyncRecvBytes_ += recvBytes;
        ++dsmcAsyncCalls_;
    }
    const scalar migrationWall =
        std::chrono::duration<scalar>(tEnd - st.tStart).count();
    migrationWallTime_ += migrationWall;
    migrationPackWallTime_ +=
        std::chrono::duration<scalar>(st.tPackEnd - st.tStart).count();
    migrationLocalPrepWallTime_ +=
        std::chrono::duration<scalar>
        (
            st.tLocalPrepEnd - st.tPackEnd
        ).count();
    migrationSizeExchangeWallTime_ +=
        std::chrono::duration<scalar>
        (
            st.tSizeExEnd - st.tLocalPrepEnd
        ).count();
    migrationRequestPostWallTime_ +=
        std::chrono::duration<scalar>
        (
            st.tReqPostEnd - st.tSizeExEnd
        ).count();
    dsmcAsyncWaitWall_ += chunkWaitWall;
    dsmcAsyncDeserWall_ += chunkDeserializeWall;
    dsmcAsyncPostWall_ +=
        std::chrono::duration<double>(tTransferEnd - st.tReqPostEnd).count();
    migrationWaitWallTime_ += chunkWaitWall;
    migrationDeserializeWallTime_ += chunkDeserializeWall;
    migrationPostWallTime_ +=
        std::chrono::duration<scalar>(tTransferEnd - st.tReqPostEnd).count()
      + std::chrono::duration<scalar>(tEnd - tTransferEnd).count();
    // Note: the gap build runs between begin and finish and is accounted
    // to buildCellOccupancy, not to migration — by design (overlap).
    const scalar migrationComm =
        std::chrono::duration<scalar>
        (
            st.tSizeExEnd - st.tLocalPrepEnd
        ).count()
      + std::chrono::duration<scalar>
        (
            st.tReqPostEnd - st.tSizeExEnd
        ).count()
      + chunkWaitWall;
    addMigrationProfileSample(migrationWall, migrationComm);

    if (cloud_.emitStepDiagnostics())
    {
        label globalMigrated = 0;
        MPI_Allreduce
        (
            &st.nMigratedOut,
            &globalMigrated,
            1,
            MPI_INT,
            MPI_SUM,
            MPI_COMM_WORLD
        );
        Info<< "Replicated mesh migration[" << migrationCalls_ << "]: "
            << "rank " << myRank_
            << " sent " << st.nMigratedOut << " parcels, received " << nRecv
            << " (global " << globalMigrated << ")"
            << " wall " << migrationWall << "s"
            << " [async overlap, incremental=" << incremental << "]"
            << endl;
    }

    st.active = false;
    asyncMigrationPending_ = false;
}

void dsmcReplicatedMesh::migrateBegin()
{
    if (!active_) return;

    if (useFlatTransfer_)
    {
        if
        (
            asyncOverlapEnabled()
         && migrationProfileContext_ == migrationProfileRegular
         && stepCounter_ != 0
        )
        {
            flatMigrateBegin();
            return;
        }
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

    // clean N1: count records per destination first, then write the count
    // header ahead of the framed records.
    labelList nSendParcels(nProcs_, 0);
    if (useOrderedTraversal)
    {
        const auto& ordered = cloud_.moveOrderedParcels();
        forAll(ordered, i)
        {
            const label dstRank = cellOwner_[ordered[i]->cell()];
            if (dstRank != myRank_) ++nSendParcels[dstRank];
        }
    }
    else
    {
        forAllConstIter(Cloud<dsmcParcel>, cloud_, iter)
        {
            const label dstRank = cellOwner_[iter().cell()];
            if (dstRank != myRank_) ++nSendParcels[dstRank];
        }
    }
    for (label i = 0; i < nProcs_; ++i)
    {
        if (i != myRank_)
        {
            sendStreams[i] << nSendParcels[i] << token::SPACE;
        }
    }

    if (useOrderedTraversal)
    {
        const auto& ordered = cloud_.moveOrderedParcels();
        for (label i = 0; i < ordered.size(); ++i)
        {
            dsmcParcel& p = *ordered[i];
            const label dstRank = cellOwner_[p.cell()];
            if (dstRank != myRank_)
            {
                dsmcNFSentParcels_ += 1.0;
                dsmcNFSendBytes_ += double
                (
                    dsmcWriteNonFlatRecord(sendStreams[dstRank], p)
                );
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
                dsmcNFSentParcels_ += 1.0;
                dsmcNFSendBytes_ += double
                (
                    dsmcWriteNonFlatRecord(sendStreams[dstRank], p)
                );
                toDelete.append(&p);
                ++nMigratedOut;
            }
            else
            {
                kept.append(&p);
            }
        }
    }
    ++dsmcNFCalls_;

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

        if (dsmcVerifyTransferRecords_ && dsmcNFCalls_ <= 1)
        {
            Info<< "N1 debug: rank " << myRank_ << " send " << sendSize
                << " bytes, head:";
            for (label b = 0; b < min(label(48), sendSize); ++b)
            {
                const unsigned char ch =
                    static_cast<unsigned char>(sendBufs[peer][b]);
                if (ch >= 32 && ch < 127) Info<< char(ch);
                else Info<< "\\x" << std::hex << int(ch) << std::dec;
            }
            Info<< endl;
        }

        // Exchange sizes (blocking — cheap)
        label recvSize = 0;
        MPI_Sendrecv(&sendSize, 1, MPI_INT, peer, 1,
                     &recvSize, 1, MPI_INT, peer, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        asyncRecvSize_ = recvSize;
        asyncRecvBuf_.setSize(recvSize);
        asyncRecvSizes_.setSize(nProcs_, 0);
        asyncRecvSizes_[peer] = recvSize;

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

        // This branch obtains the per-peer receive sizes with MPI_Get_count
        // in migrateFinish, so no stored sizes are needed here.
        asyncRecvSizes_.clear();

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
        asyncRecvSizes_ = recvSizes;

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

    // comm_further_opt.md §2: flat-path async overlap
    if (dsmcFlatAsync_.active)
    {
        flatMigrateFinish();
        return;
    }

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
                const label nRecvParcels =
                    dsmcReadNonFlatCount(is, actualBytes);
                for (label k = 0; k < nRecvParcels; ++k)
                {
                    auto* newp = dsmcReadNonFlatRecord
                    (
                        mesh_, is, actualBytes, k, nRecvParcels
                    );
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

            // clean N1: asyncRecvBuf_ holds one framed stream per peer (the
            // multi-rank branch concatenates them, the 2-rank branch has a
            // single peer).  Each stream starts with its own parcel count,
            // so the buffer must be split by asyncRecvSizes_ — reading one
            // count over the whole buffer would drop every other peer.
            if (asyncRecvSizes_.size() == nProcs_)
            {
                label off = 0;
                for (label i = 0; i < nProcs_; ++i)
                {
                    if (i == myRank_ || asyncRecvSizes_[i] <= 0) continue;

                    ISpanStream is
                    (
                        asyncRecvBuf_.data() + off,
                        asyncRecvSizes_[i],
                        IOstream::BINARY
                    );
                    const label nRecvParcels =
                        dsmcReadNonFlatCount(is, asyncRecvSizes_[i]);
                    for (label k = 0; k < nRecvParcels; ++k)
                    {
                        auto* newp = dsmcReadNonFlatRecord
                        (
                            mesh_, is, asyncRecvSizes_[i], k, nRecvParcels
                        );
                        cloud_.addParticle(newp);
                        received.append(newp);
                    }
                    off += asyncRecvSizes_[i];
                }
            }
            else
            {
                ISpanStream is(asyncRecvBuf_.data(), asyncRecvSize_,
                               IOstream::BINARY);
                const label nRecvParcels =
                    dsmcReadNonFlatCount(is, asyncRecvSize_);
                for (label k = 0; k < nRecvParcels; ++k)
                {
                    auto* newp = dsmcReadNonFlatRecord
                    (
                        mesh_, is, asyncRecvSize_, k, nRecvParcels
                    );
                    cloud_.addParticle(newp);
                    received.append(newp);
                }
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
    asyncRecvSizes_.clear();
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
    if (dsmcCountsPosted_)
    {
        // §34.5: complete the gather posted at the migration tail; the
        // ordered-view rebuild (setMoveOrderedParcels) ran in between and
        // overlapped the collective's convergence.  Note the wait may
        // have moved rather than vanished — judge by total_loop.
        const auto tUpcWait0 = std::chrono::steady_clock::now();
        MPI_Wait(&dsmcCountsReq_, MPI_STATUS_IGNORE);
        dsmcMigUpcWaitWall_ +=
            std::chrono::duration<scalar>
            (std::chrono::steady_clock::now() - tUpcWait0).count();
        dsmcCountsPosted_ = false;
    }
    else
    {
        localParticleCount_ = cloud_.size();
        allParticleCounts_.setSize(nProcs_);
        MPI_Allgather(&localParticleCount_, 1, MPI_INT,
                      allParticleCounts_.data(), 1, MPI_INT, MPI_COMM_WORLD);
    }
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


void dsmcReplicatedMesh::postParticleCountsGather()
{
    // §34.5: nonblocking gather posted after all received parcels have
    // been appended (cloud_.size() is final for this step) and before the
    // ordered-view rebuild.
    const auto tUpcPost0 = std::chrono::steady_clock::now();
    localParticleCount_ = cloud_.size();
    allParticleCounts_.setSize(nProcs_);
    MPI_Iallgather
    (
        &localParticleCount_, 1, MPI_INT,
        allParticleCounts_.data(), 1, MPI_INT,
        MPI_COMM_WORLD, &dsmcCountsReq_
    );
    dsmcCountsPosted_ = true;
    dsmcMigUpcPostWall_ +=
        std::chrono::duration<scalar>
        (std::chrono::steady_clock::now() - tUpcPost0).count();
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
        // clean N1: count header before the framed records.
        os << label(cloud_.size()) << token::SPACE;
        forAllIter(Cloud<dsmcParcel>, cloud_, iter)
        {
            dsmcNFSentParcels_ += 1.0;
            dsmcNFSendBytes_ += double(dsmcWriteNonFlatRecord(os, iter()));
        }
        ++dsmcNFCalls_;
        cloud_.clear();

        const string s = os.str();
        // H7: DynamicList(size) only reserves capacity (addressed size 0);
        // the previous code then sent a zero-length message and dropped the
        // particles.  Use a sized buffer so sendSize equals the payload.
        List<char> buf(s.size());
        forAll(s, si)
        {
            buf[si] = s[si];
        }
        const label sendSize = buf.size();
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
                const label nRecvParcels =
                    dsmcReadNonFlatCount(is, recvSize);
                for (label k = 0; k < nRecvParcels; ++k)
                {
                    auto* newp = dsmcReadNonFlatRecord
                    (
                        mesh_, is, recvSize, k, nRecvParcels
                    );
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
        0.0,   // clean M4: write phase not instrumented (reported n/a)
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

    // §35 P0: per-rank stage rows (must also run on every rank before the
    // early return).  Layout: classify pack unpack* append* prep sizeX
    // postReq wait deser post rebuild* upcPost upcWait+ parcelsOut
    // parcelsIn sendMB recvMB totalWall chunks peerAvg   (* = sub-item of
    // deser / post, + = outside the migration interval; excluded from the
    // sum).
    const int nStage = 22;
    scalarList stageRows(nStage*nProcs_, 0.0);
    {
        double myStage[nStage] =
        {
            dsmcMigClassifyWall_,
            dsmcMigPackWall_,
            dsmcMigUnpackWall_,
            dsmcMigAppendWall_,
            double(migrationLocalPrepWallTime_),
            double(migrationSizeExchangeWallTime_),
            double(migrationRequestPostWallTime_),
            double(migrationWaitWallTime_),
            double(migrationDeserializeWallTime_),
            double(migrationPostWallTime_),
            dsmcMigRebuildWall_,
            dsmcMigUpcPostWall_,
            dsmcMigUpcWaitWall_,
            double(totalParcelsMigrated_),
            dsmcMigParcelsIn_,
            dsmcMigSendBytes_/1.0e6,
            dsmcMigRecvBytes_/1.0e6,
            dsmcMigTotalWall_,
            dsmcMigChunkCount_,
            dsmcMigPeerCount_/double(max(migrationCalls_, label(1))),
            dsmcMigOccAppendWall_,
            dsmcMigOccValidateWall_
        };
        MPI_Allgather
        (
            myStage, nStage, MPI_DOUBLE,
            stageRows.data(), nStage, MPI_DOUBLE,
            MPI_COMM_WORLD
        );
    }

    // §36.3: per-rank unpack allocation counters.  The magazine and the
    // per-thread counters live in the parcel process, so these are process
    // totals (cumulative like the P0 stage rows); gathered here because
    // the report body only runs on rank 0.
    const int nAlloc = 6;
    scalarList allocRows(nAlloc*nProcs_, 0.0);
    {
        double myAlloc[nAlloc] =
        {
            double(dsmcParcel::unpackAllocCount()),
            double(dsmcParcel::unpackMagazineRefillCount()),
            double(dsmcParcel::unpackArenaBlockCount()),
            double(dsmcParcel::unpackArenaInUse()),
            double(dsmcParcel::unpackVibReallocCount()),
            double(dsmcParcel::unpackVibMaxModes())
        };
        MPI_Allgather
        (
            myAlloc, nAlloc, MPI_DOUBLE,
            allocRows.data(), nAlloc, MPI_DOUBLE,
            MPI_COMM_WORLD
        );
    }

    // §35 P2: async-path per-rank stages (collective before the return).
    const int nAsyncStage = 12;
    scalarList asyncRows(nAsyncStage*nProcs_, 0.0);
    {
        double myAsync[nAsyncStage] =
        {
            dsmcAsyncClassifyPackWall_,
            dsmcAsyncPrepWall_,
            dsmcAsyncPostReqWall_,
            dsmcAsyncFinishWall_,
            dsmcAsyncTotalWall_,
            dsmcAsyncParcelsOut_,
            dsmcAsyncParcelsIn_,
            dsmcAsyncSendBytes_/1.0e6,
            dsmcAsyncRecvBytes_/1.0e6,
            dsmcAsyncWaitWall_,
            dsmcAsyncDeserWall_,
            dsmcAsyncPostWall_
        };
        MPI_Allgather
        (
            myAsync, nAsyncStage, MPI_DOUBLE,
            asyncRows.data(), nAsyncStage, MPI_DOUBLE,
            MPI_COMM_WORLD
        );
    }

    // H6/N1: transfer verification and non-flat transfer totals.  The
    // counters are per rank, so they are reduced here, before the
    // rank-0-only report body.
    double transferStats[10] =
    {
        0.0, 0.0,
        dsmcNFSentParcels_, dsmcNFRecvParcels_,
        dsmcNFSendBytes_, dsmcNFRecvBytes_,
        dsmcNFCalls_, 0.0,
        dsmcNFRecordsVerified_, dsmcNFStuckTrackedVerified_
    };
    {
        label nVerified = 0;
        label nMismatch = 0;
        dsmcParcel::transferRoundTripStats(nVerified, nMismatch);
        transferStats[0] = double(nVerified);
        transferStats[1] = double(nMismatch);
    }
    MPI_Allreduce
    (
        MPI_IN_PLACE, transferStats, 10, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD
    );

    // H4: incremental occupancy append coverage (collective).
    double occCoverage[4] = {0.0, 0.0, 0.0, 0.0};
    cloud_.occupancyAppendCoverage
    (
        occCoverage[0], occCoverage[1], occCoverage[2], occCoverage[3]
    );
    MPI_Allreduce
    (
        MPI_IN_PLACE, occCoverage, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD
    );

    // H4 self-test: run the real append with the input shapes that a
    // production schedule rarely produces (empty, single cell, last cell).
    if (dsmcOccSelfTestMode_ != "none")
    {
        label nScratch = 0;
        label nErrors = 0;
        const bool ok =
            cloud_.occupancySelfTest(dsmcOccSelfTestMode_, nScratch, nErrors);

        double selfStats[3] =
        {
            ok ? 1.0 : 0.0, double(nScratch), double(nErrors)
        };
        MPI_Allreduce
        (
            MPI_IN_PLACE, selfStats, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD
        );
        if (myRank_ == 0)
        {
            Info<< "    H4 self-test (" << dsmcOccSelfTestMode_ << "): ranks ok "
                << selfStats[0] << "/" << nProcs_
                << ", scratch parcels " << selfStats[1]
                << ", layout errors " << selfStats[2] << nl;
        }
    }

    if (myRank_ != 0) return;

    // §35 P0: per-rank regular-migration stage table.  The "max" phase
    // values above come from different ranks and cannot be added; this
    // table gives each rank's own stage totals so the decomposition is
    // addable on one rank and the total-vs-sum gap is explicit.
    {
        const int nStage = 22;
        // NOTE (clean C2b): when the async overlap path is active, the
        // columns prep/sizeX/postReq/wait/deser/post are shared with it and
        // cannot be tiled against the sync totalWall; use the P2 table for
        // the async interval and a sync-only run (async off, the production
        // configuration) for this decomposition.
        Info<< nl << "    §35 P0 per-rank migration stage totals (s / "
            << migrationCalls_ << " migration calls; * sub-item, + outside"
            << " migration; async path active on a rank -> shared columns"
            << " cannot be tiled, see P2):" << nl
            << "      rank classify pack unpack* append* prep sizeX postReq"
            << " wait deser post rebuild* upcPost upcWait+ parcelsOut"
            << " parcelsIn sendMB recvMB totalWall chunks peerAvg"
            << " occAppend* occValid*" << nl;
        for (label r = 0; r < nProcs_; ++r)
        {
            const scalar* row = &stageRows[nStage*r];
            Info<< "      r" << r;
            for (int k = 0; k < nStage; ++k)
            {
                Info<< " " << row[k];
            }
            Info<< nl;
        }
        // Tiling stages on one rank: classify+pack+prep+sizeX+postReq+wait
        // +deser+post+upcPost == totalWall (rebuild ⊂ post, unpack/append
        // ⊂ deser, upcWait happens after the migration returns).
        Info<< "      reconciliation: tiling sum vs totalWall per rank:"
            << nl;
        const int idx[9] = {0, 1, 4, 5, 6, 7, 8, 9, 11};
        for (label r = 0; r < nProcs_; ++r)
        {
            const scalar* row = &stageRows[nStage*r];
            scalar sum = 0.0;
            for (int k = 0; k < 9; ++k)
            {
                sum += row[idx[k]];
            }
            Info<< "      r" << r
                << " sum=" << sum
                << " totalWall=" << row[17]
                << " diff=" << (sum - row[17]) << nl;
        }
    }

    // §35 P2 (clean C2b): async-overlap path, its own parent interval.
    // classify+pack share the begin interval; sizeX is inside prep
    // (tSizeExEnd == tLocalPrepEnd); finish is one parent covering
    // wait+unpack+append+occupancy rebuild+post.  Sync-only totals are the
    // P0 columns minus these values.
    {
        Info<< nl << "    §35 P2 async-path migration stages (s / "
            << dsmcAsyncCalls_ << " async calls; classify+pack combined,"
            << " sizeX inside prep, finish = wait+unpack+append+rebuild+post):"
            << nl
            << "      rank classifyPack prep postReq finish totalWall"
            << " parcelsOut parcelsIn sendMB recvMB wait* deser* post*" << nl;
        for (label r = 0; r < nProcs_; ++r)
        {
            const scalar* row = &asyncRows[nAsyncStage*r];
            Info<< "      r" << r;
            for (int k = 0; k < nAsyncStage; ++k)
            {
                Info<< " " << row[k];
            }
            Info<< nl;
        }
        Info<< "      reconciliation (async): tiling sum vs totalWall per rank:"
            << nl;
        for (label r = 0; r < nProcs_; ++r)
        {
            const scalar* row = &asyncRows[nAsyncStage*r];
            // wait/deser/post are sub-items of finish (not added twice).
            const scalar sum = row[0] + row[1] + row[2] + row[3];
            Info<< "      r" << r
                << " sum=" << sum
                << " totalWall=" << row[4]
                << " diff=" << (sum - row[4]) << nl;
        }
    }

    // §36.3/§36.4: what the unpack path actually allocated.  arenaSlots is
    // the real allocation count (magazine-batched at 64/refill); refills is
    // the number of arena lock acquisitions; arenaBlocks is the ::malloc
    // count (4096 parcels per block); vibReallocs counts the vibScratch
    // growths.  The counters are updated outside the arena lock (§36.4) —
    // inside it, queueing on the lock amplified them into a 5x unpack
    // slowdown.
    {
        const int nAlloc = 6;
        Info<< "      §36.4 unpack allocations (per-rank process totals;"
            << " arenaSlots = real count, magazine-batched):" << nl
            << "      rank arenaSlots refills arenaBlocks liveSlots"
            << " vibReallocs vibMaxModes" << nl;
        double totSlots = 0.0;
        double totBlocks = 0.0;
        double totRefills = 0.0;
        double totVib = 0.0;
        for (label r = 0; r < nProcs_; ++r)
        {
            const scalar* row = &allocRows[nAlloc*r];
            Info<< "      r" << r;
            for (int k = 0; k < nAlloc; ++k)
            {
                Info<< " " << row[k];
            }
            Info<< nl;
            totSlots += row[0];
            totRefills += row[1];
            totBlocks += row[2];
            totVib += row[4];
        }
        Info<< "      sum arenaSlots=" << totSlots
            << " refills=" << totRefills
            << " arenaBlocks=" << totBlocks
            << " vibReallocs=" << totVib
            << " blockBytes=" << scalar(sizeof(dsmcParcel)*4096)
            << " arenaBytes=" << scalar(totBlocks*sizeof(dsmcParcel)*4096)
            << nl;
    }

    // H6/N1: transfer round-trip verification and non-flat transfer totals
    // (reduced above, before the rank-0-only report body).  The non-flat
    // path is the legacy serializer (replicatedMeshFlatTransfer false), so
    // its global send/receive equality is a direct conservation check.
    {
        if
        (
            transferStats[0] > 0.0
         || transferStats[2] > 0.0
         || transferStats[3] > 0.0
        )
        {
            Info<< nl << "    transfer verification:"
                << nl << "      flat round-trip checked = " << transferStats[0]
                << " records (mismatches " << transferStats[1] << ")"
                << nl << "      non-flat calls          = " << transferStats[6]
                << nl << "      non-flat sent           = " << transferStats[2]
                << " parcels (" << transferStats[4] << " bytes)"
                << nl << "      non-flat received       = " << transferStats[3]
                << " parcels (" << transferStats[5] << " bytes)"
                << nl << "      non-flat conservation   = "
                << ((transferStats[2] == transferStats[3])
                    ? "sent == received" : "MISMATCH")
                << nl << "      non-flat round-trip     = "
                << transferStats[8] << " records checked"
                << " (stuck/tracked synthetic "
                << transferStats[9] << ")"
                << nl;
        }
    }

    // H4: measured shapes of the incremental occupancy append inputs.
    {
        Info<< "    H4 occupancy append coverage (rank sum):"
            << " empty " << occCoverage[0]
            << " singleCell " << occCoverage[1]
            << " multiCell " << occCoverage[2]
            << " lastCell " << occCoverage[3]
            << nl;
    }

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
            << "    DLB cellOwner write max [s]  = (not instrumented)"
            << nl
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

// * * * * * * dlb_offload.md P2-P5: per-step move offload * * * * * * * //

namespace
{

constexpr label offloadTagBase = 4100;
constexpr label offloadStepSlots = 64;
constexpr label offloadMaxChunksPerPair = 7;

inline label offloadTagHdr(const label step)
{
    return offloadTagBase + (step % offloadStepSlots)*(offloadMaxChunksPerPair + 1);
}

inline label offloadTagChunk(const label step, const label chunkI)
{
    return offloadTagHdr(step) + 1 + chunkI;
}

inline double offloadNodeHash()
{
    char name[MPI_MAX_PROCESSOR_NAME];
    int len = 0;
    MPI_Get_processor_name(name, &len);
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < len; ++i)
    {
        h ^= static_cast<unsigned char>(name[i]);
        h *= 1099511628211ULL;
    }
    return double(h);
}

// Order-independent per-particle hash over the transfer-relevant state
// (field-wise, so struct padding cannot leak into the checksum).
inline uint64_t offloadParcelHash(const dsmcParcel& p)
{
    uint64_t h = 1469598103934665603ULL;
    auto mix = [&](const unsigned char* b, size_t n)
    {
        for (size_t i = 0; i < n; ++i)
        {
            h ^= b[i];
            h *= 1099511628211ULL;
        }
    };
    const scalar px[3] = {p.position().x(), p.position().y(), p.position().z()};
    const scalar uu[3] = {p.U().x(), p.U().y(), p.U().z()};
    const label cl = p.cell();
    const label tf = p.tetFace();
    const label tp = p.tetPt();
    const label ty = p.typeId();
    const label np = p.newParcel();
    const label cf = p.classification();
    const label el = p.ELevel();
    const scalar rw = p.RWF();
    const scalar er = p.ERot();
    mix(reinterpret_cast<const unsigned char*>(px), sizeof(px));
    mix(reinterpret_cast<const unsigned char*>(uu), sizeof(uu));
    mix(reinterpret_cast<const unsigned char*>(&cl), sizeof(cl));
    mix(reinterpret_cast<const unsigned char*>(&tf), sizeof(tf));
    mix(reinterpret_cast<const unsigned char*>(&tp), sizeof(tp));
    mix(reinterpret_cast<const unsigned char*>(&ty), sizeof(ty));
    mix(reinterpret_cast<const unsigned char*>(&np), sizeof(np));
    mix(reinterpret_cast<const unsigned char*>(&cf), sizeof(cf));
    mix(reinterpret_cast<const unsigned char*>(&el), sizeof(el));
    mix(reinterpret_cast<const unsigned char*>(&rw), sizeof(rw));
    mix(reinterpret_cast<const unsigned char*>(&er), sizeof(er));
    return h;
}

struct offloadRankInfo
{
    double work;
    double nP;
    double node;
    double dispatchRate;   // donor EMA of dispatch seconds per parcel
};

} // anonymous namespace


void Foam::dsmcReplicatedMesh::offloadReadControls()
{
    const dictionary& cd = mesh_.time().controlDict();
    offloadEnabled_ = cd.lookupOrDefault<bool>("offloadEnable", false);
    offloadPlanMode_ = cd.lookupOrDefault<word>("offloadPlanMode", "singlePairSync");
    if (offloadPlanMode_ == "dumpOnly")
    {
        // clean M11: dumpOnly keeps raw parcel pointers across the move and
        // dereferences them afterwards; boundary models can delete parcels
        // during the move, so the snapshot is not lifetime-safe.
        FatalErrorInFunction
            << "offloadPlanMode dumpOnly is disabled: it holds raw parcel "
            << "pointers across the move (use-after-free when a boundary "
            << "deletes a parcel).  See clean-20261101 M11."
            << exit(FatalError);
    }
    if
    (
        offloadPlanMode_ != "off"
     && offloadPlanMode_ != "dumpOnly"
     && offloadPlanMode_ != "singlePairSync"
     && offloadPlanMode_ != "fixedPairs"
     && offloadPlanMode_ != "greedy"
    )
    {
        FatalErrorInFunction
            << "Unknown offloadPlanMode '" << offloadPlanMode_
            << "' (off|dumpOnly|singlePairSync|fixedPairs|greedy)"
            << exit(FatalError);
    }
    offloadMaxPairs_ =
        max(1, cd.lookupOrDefault<label>("offloadMaxPairs", 2));
    offloadPairBytesCap_ =
        max(1, cd.lookupOrDefault<label>("offloadPairBytesCapMB", 256))
       *1024*1024;
    offloadChunkParcels_ =
        max(1000, cd.lookupOrDefault<label>("offloadChunkParcels", 400000));
    offloadMinGain_ =
        cd.lookupOrDefault<scalar>("offloadMinGainMs", 2.0)*1e-3;
    offloadSkipDLBSteps_ =
        cd.lookupOrDefault<bool>("offloadSkipDLBSteps", true);
    offloadStateDumpStep_ =
        cd.lookupOrDefault<label>("offloadStateDumpStep", -1);
    offloadControlsRead_ = true;

    if (offloadEnabled_ && myRank_ == 0)
    {
        Info<< "Offload: enabled (mode " << offloadPlanMode_
            << ", maxPairs " << offloadMaxPairs_
            << ", pairBytesCapMB "
            << offloadPairBytesCap_/(1024*1024)
            << ", chunkParcels " << offloadChunkParcels_
            << ", dumpStep " << offloadStateDumpStep_ << ")" << endl;
    }
}


double Foam::dsmcReplicatedMesh::offloadRelNow() const
{
    return std::chrono::duration<double>
    (
        std::chrono::steady_clock::now().time_since_epoch()
    ).count() - cloud_.offloadProbeEpoch();
}


void Foam::dsmcReplicatedMesh::offloadPairsCsvLine
(
    const label step,
    const char* role,
    const double tDispatch,
    const double tSendDone
)
{
    if (!offloadPairsFileOpen_)
    {
        offloadPairsFile_.open
        (
            "offloadPairs-r" + std::to_string(myRank_) + ".csv"
        );
        offloadPairsFile_
            << "step,mode,planId,reject,role,pairDonor,pairHelper,"
               "planPairs,nSelected,nRecv,nExecViol,nInserted,planBytes,"
               "tDispatch,tSendDone,tRecvHeader,tRecvDone,tExecDone,"
               "tInsertDone\n";
        offloadPairsFileOpen_ = true;
    }
    offloadPairsFile_
        << step << ',' << offloadPlanMode_ << ',' << offloadPlanId_
        << ',' << offloadPlanReject_ << ',' << role
        << ',' << offloadPairDonor_ << ',' << offloadPairHelper_
        << ',' << offloadPlanNPairs_
        << ',' << offloadNSelected_ << ',' << offloadNRecv_
        << ',' << offloadNExecViolation_ << ',' << offloadNInserted_
        << ',' << offloadPlanBytes_
        << ',' << tDispatch << ',' << tSendDone
        << ',' << offloadTRecvHeader_ << ',' << offloadTRecvDone_
        << ',' << offloadTExecDone_ << ',' << offloadTInsertDone_
        << '\n';
    // Lustre 上逐行 endl 刷盘代价过高：缓冲写，每 50 行 flush 一次
    if ((step % 50) == 0)
    {
        offloadPairsFile_.flush();
    }
}


void Foam::dsmcReplicatedMesh::offloadDispatch()
{
    if (!offloadControlsRead_)
    {
        offloadReadControls();
    }
    if (!offloadEnabled_ || offloadPlanMode_ == "off")
    {
        return;
    }

    // timeIndex is restart-safe (stepCounter restarts at 0 on fork runs)
    const label step = mesh_.time().timeIndex();

    // reset per-step view
    offloadIAmDonor_ = false;
    offloadIAmHelper_ = false;
    offloadPairDonor_ = -1;
    offloadPairHelper_ = -1;
    offloadPlanNPairs_ = 0;
    offloadPlanReject_ = 0;
    offloadNSelected_ = 0;
    offloadNRecv_ = 0;
    offloadNExecViolation_ = 0;
    offloadNInserted_ = 0;
    offloadPlanBytes_ = 0;
    offloadNTotal_ = 0;
    offloadNChunks_ = 0;
    offloadTDispatch_ = offloadRelNow();
    offloadTSendDone_ = -1.0;
    offloadTRecvHeader_ = -1.0;
    offloadTRecvDone_ = -1.0;
    offloadTExecDone_ = -1.0;
    offloadTInsertDone_ = -1.0;

    // prior sends must be complete before tag reuse
    if (offloadSendsPending_)
    {
        MPI_Waitall
        (
            offloadChunkSendReqs_.size(),
            offloadChunkSendReqs_.begin(),
            MPI_STATUSES_IGNORE
        );
        MPI_Wait(&offloadHdrSendReq_, MPI_STATUS_IGNORE);
        offloadSendsPending_ = false;
    }

    // ---- plan (rank 0 decides, broadcast) ----
    offloadRankInfo info;
    info.work = double(offloadPredictWork());
    info.nP = double(cloud_.size());
    info.node = offloadNodeHash();
    info.dispatchRate = offloadDispatchRate_;
    List<offloadRankInfo> all(nProcs_);
    MPI_Allgather
    (
        &info, sizeof(offloadRankInfo), MPI_BYTE,
        all.begin(), sizeof(offloadRankInfo), MPI_BYTE,
        MPI_COMM_WORLD
    );

    // fresh-run fill step: no useful predictions and no owner-aligned
    // particle distribution yet (fork runs re-align on their first
    // step, which also has zeroed predictions and is skipped the same
    // way by the planner's viability check)
    bool skip = (stepCounter() == 0 || cloud_.size() == 0);
    if (offloadPlanMode_ == "dumpOnly")
    {
        skip = (step != offloadStateDumpStep_);
    }
    else if (offloadSkipDLBSteps_)
    {
        if (postDLBSnapshotCountdown_ > 0)
        {
            skip = true;
        }
        if (lastAutoRebalanceStep_ >= 0 && lastAutoRebalanceStep_ == stepCounter())
        {
            skip = true;
        }
        if
        (
            offloadLastRebalanceCount_ >= 0
         && autoRebalanceCount_ != offloadLastRebalanceCount_
        )
        {
            skip = true;
        }
    }
    offloadLastRebalanceCount_ = autoRebalanceCount_;
    if (asyncMigrationPending_)
    {
        skip = true;
    }

    // Fixed-size plan so the Bcast count is identical on every rank
    // (root used to send 3+3*nPairs while receivers expected 3 — a
    // collective mismatch that deadlocked the planner rank).
    const label planLen = 3 + 3*max(offloadMaxPairs_, label(1));
    labelList plan(planLen, -1);
    plan[0] = 0;
    plan[1] = 0;
    plan[2] = 0;
    if (myRank_ == 0)
    {
        DynamicList<label> pf;
        if (!skip)
        {
            if
            (
                offloadPlanMode_ == "singlePairSync"
             || offloadPlanMode_ == "dumpOnly"
            )
            {
                label d = -1;
                label h = -1;
                double wmax = -GREAT;
                double wmin = GREAT;
                for (label r = 0; r < nProcs_; ++r)
                {
                    if (all[r].work > wmax) { wmax = all[r].work; d = r; }
                    if (all[r].work < wmin) { wmin = all[r].work; h = r; }
                }
                if (d >= 0 && h >= 0 && d != h && (wmax - wmin) > offloadMinGain_)
                {
                    const double rate =
                        all[d].work/max(all[d].nP, 1.0);
                    // clean M3: bound the selection by the byte cap and
                    // reject non-positive plans before any allocation or
                    // communication (the greedy branch already did).
                    const label capN =
                        offloadPairBytesCap_
                      / label(sizeof(dsmcParcel::TransferData));
                    const label nSel = min
                    (
                        label((wmax - wmin)/2.0/max(rate, 1e-12)),
                        capN
                    );
                    if (nSel <= 0)
                    {
                        pf.clear();
                        plan[1] = 4;
                    }
                    else
                    {
                    pf.append(d);
                    pf.append(h);
                    pf.append(nSel);
                    }
                }
                else
                {
                    pf.clear();
                    plan[1] = 3;
                }
            }
            else // fixedPairs / greedy: the same greedy loop
            {
                List<double> a(nProcs_);
                forAll(all, r)
                {
                    a[r] = all[r].work;
                }
                boolList used(nProcs_, false);
                for (label p = 0; p < offloadMaxPairs_; ++p)
                {
                    label d = -1;
                    label h = -1;
                    double wmax = -GREAT;
                    double wmin = GREAT;
                    for (label r = 0; r < nProcs_; ++r)
                    {
                        if (used[r]) continue;
                        if (a[r] > wmax) { wmax = a[r]; d = r; }
                        if (a[r] < wmin) { wmin = a[r]; h = r; }
                    }
                    if (d < 0 || h < 0 || d == h) break;
                    if (wmax - wmin <= offloadMinGain_) break;
                    const double rateD = all[d].work/max(all[d].nP, 1.0);
                    const double rateH = all[h].work/max(all[h].nP, 1.0);
                    const label capN =
                        offloadPairBytesCap_
                      / label(sizeof(dsmcParcel::TransferData));
                    // donor 净收益 = nSel*(rateD - dispatchRate)：派发
                    // （扫描/打包/删除/张贴）按实测 EMA 计入，模型自限
                    const double netRate = rateD - all[d].dispatchRate;
                    if (netRate <= 0) break;
                    const label nSel = min
                    (
                        label(((wmax - wmin)/2.0)/netRate),
                        capN
                    );
                    if (nSel <= 0) break;
                    pf.append(d);
                    pf.append(h);
                    pf.append(nSel);
                    used[d] = true;
                    used[h] = true;
                    a[d] = max(wmin, wmax - nSel*netRate);
                    a[h] = wmin + nSel*rateH;
                }
                double oldMax = 0.0;
                double newMax = 0.0;
                forAll(all, r)
                {
                    oldMax = max(oldMax, all[r].work);
                }
                forAll(a, r)
                {
                    newMax = max(newMax, a[r]);
                }
                if (pf.empty() || (oldMax - newMax) < offloadMinGain_)
                {
                    pf.clear();
                    plan[1] = 2;
                }
            }
        }
        else
        {
            plan[1] = 1;
        }
        plan[0] = offloadPlanId_ + 1;
        // broadcast: [planId, reject, nPairs, (d,h,nSel) x maxPairs]
        plan[1] = skip ? 1 : plan[1];
        plan[2] = pf.size()/3;
        forAll(pf, i)
        {
            plan[3 + i] = pf[i];
        }
    }
    MPI_Bcast(plan.begin(), plan.size(), MPI_INT, 0, MPI_COMM_WORLD);
    offloadPlanId_ = plan[0];
    offloadPlanReject_ = plan[1];
    offloadPlanNPairs_ = plan[2];
    for (label p = 0; p < offloadPlanNPairs_; ++p)
    {
        const label d = plan[3 + 3*p];
        const label h = plan[4 + 3*p];
        const label nSel = plan[5 + 3*p];
        if (d == myRank_)
        {
            offloadIAmDonor_ = true;
            offloadPairDonor_ = d;
            offloadPairHelper_ = h;
            offloadPlanBytes_ = nSel*label(sizeof(dsmcParcel::TransferData));
        }
        if (h == myRank_)
        {
            offloadIAmHelper_ = true;
            offloadPairDonor_ = d;
            offloadPairHelper_ = h;
            offloadPlanBytes_ = nSel*label(sizeof(dsmcParcel::TransferData));
        }
    }

    // ---- donor: select + (dump | pack+send+delete) ----
    if (offloadIAmDonor_ && offloadPlanReject_ == 0)
    {
        const label nSelMax = offloadPlanBytes_
            / label(sizeof(dsmcParcel::TransferData));
        DynamicList<dsmcParcel*> sel(min(nSelMax, cloud_.size()) + 1);
        // §4.2: scan the ordered array when valid (cache-friendly); the
        // linked list is the fallback and its pointer chasing dominated
        // the first P3 dispatch measurements
        const bool orderedScan = cloud_.hasMoveOrderedParcels();
        label scanned = 0;
        // clean M8: the parcel-side safety test uses the per-cell time
        // step, so the selector uses the same rule.
        if (orderedScan)
        {
            const List<dsmcParcel*>& arr = cloud_.moveOrderedParcels();
            for (label i = 0; i < arr.size() && sel.size() < nSelMax; ++i)
            {
                dsmcParcel& p = *arr[i];
                ++scanned;
                if (!p.isFree()) continue;
                if (p.newParcel() != -1) continue;
                const label c = p.cell();
                const scalar clearance = cloud_.offloadCellClearance(c);
                if (clearance <= 0) continue;
                const scalar dt = cloud_.deltaTValue(c);
                const scalar m2 = magSqr(p.U());
                if (clearance*clearance <= 4*m2*dt*dt) continue;
                sel.append(&p);
            }
        }
        if (sel.size() < nSelMax)
        {
            forAllIter(Cloud<dsmcParcel>, cloud_, iter)
            {
                if (sel.size() >= nSelMax) break;
                dsmcParcel& p = iter();
                if (!p.isFree()) continue;
                if (p.newParcel() != -1) continue;
                const label c = p.cell();
                const scalar clearance = cloud_.offloadCellClearance(c);
                if (clearance <= 0) continue;
                const scalar dt = cloud_.deltaTValue(c);
                const scalar m2 = magSqr(p.U());
                if (clearance*clearance <= 4*m2*dt*dt) continue;
                sel.append(&p);
            }
        }
        offloadNSelected_ = sel.size();

        // P2 dump: task input records (identical selection rule in the
        // dumpOnly reference and the offload run — same checkpoint state)
        if
        (
            step == offloadStateDumpStep_
         && (offloadPlanMode_ == "dumpOnly"
          || offloadPlanMode_ == "singlePairSync")
        )
        {
            std::ofstream os
            (
                "offloadTaskIn-r" + std::to_string(myRank_) + ".bin",
                std::ios::binary
            );
            const label n = sel.size();
            os.write
            (
                reinterpret_cast<const char*>(&n), sizeof(label)
            );
            for (label i = 0; i < n; ++i)
            {
                dsmcParcel::TransferData tdr;
                sel[i]->packTransfer(tdr);
                os.write
                (
                    reinterpret_cast<const char*>(&tdr),
                    sizeof(tdr)
                );
            }
        }

        if (offloadPlanMode_ == "dumpOnly")
        {
            // keep pointers: after the normal local move the same
            // particles are dumped as the bitwise reference
            offloadDumpTaskPtrs_ = List<dsmcParcel*>(nSelMax);
            label nKeep = 0;
            for (label i = 0; i < sel.size(); ++i)
            {
                offloadDumpTaskPtrs_[nKeep++] = sel[i];
            }
            offloadDumpTaskPtrs_.setSize(nKeep);
        }
        else
        {
            // pack + delete + front send
            offloadNTotal_ = sel.size();
            offloadNChunks_ =
                max
                (
                    1,
                    (offloadNTotal_ + offloadChunkParcels_ - 1)
                   /offloadChunkParcels_
                );
            offloadSendBuf_ =
                List<double>
                (
                    (offloadNTotal_*sizeof(dsmcParcel::TransferData))
                   /sizeof(double) + 8
                );
            dsmcParcel::TransferData* tds =
                reinterpret_cast<dsmcParcel::TransferData*>
                (
                    offloadSendBuf_.begin()
                );
            for (label i = 0; i < offloadNTotal_; ++i)
            {
                sel[i]->packTransfer(tds[i]);
            }
            for (label i = 0; i < offloadNTotal_; ++i)
            {
                cloud_.deleteParcel(sel[i]);
            }
            label hdr[6] =
            {
                step, offloadPlanId_, 0,
                offloadNTotal_, offloadNChunks_, myRank_
            };
            MPI_Isend
            (
                hdr, 6, MPI_INT, offloadPairHelper_,
                offloadTagHdr(step), MPI_COMM_WORLD,
                &offloadHdrSendReq_
            );
            offloadChunkSendReqs_ = List<MPI_Request>(offloadNChunks_);
            for (label ci = 0; ci < offloadNChunks_; ++ci)
            {
                const label begin = ci*offloadChunkParcels_;
                const label cnt =
                    min(offloadChunkParcels_, offloadNTotal_ - begin);
                MPI_Isend
                (
                    tds + begin,
                    cnt*label(sizeof(dsmcParcel::TransferData)),
                    MPI_BYTE, offloadPairHelper_,
                    offloadTagChunk(step, ci), MPI_COMM_WORLD,
                    offloadChunkSendReqs_.begin() + ci
                );
            }
            offloadSendsPending_ = true;
            if (offloadPlanMode_ == "singlePairSync")
            {
                // P2 prototype: complete the whole transfer before any
                // local move on either side
                MPI_Waitall
                (
                    offloadNChunks_,
                    offloadChunkSendReqs_.begin(),
                    MPI_STATUSES_IGNORE
                );
                MPI_Wait(&offloadHdrSendReq_, MPI_STATUS_IGNORE);
                offloadSendsPending_ = false;
            }
            offloadTSendDone_ = offloadRelNow();
            // self-calibrating dispatch cost (scan+pack+delete+post)
            const double dtDis = offloadTSendDone_ - offloadTDispatch_;
            if (offloadNSelected_ > 0 && dtDis > 0)
            {
                const double r = dtDis/offloadNSelected_;
                offloadDispatchRate_ =
                    0.7*offloadDispatchRate_ + 0.3*r;
            }
        }
    }

    // ---- helper: pre-post header receive ----
    if
    (
        offloadIAmHelper_
     && offloadPlanReject_ == 0
     && offloadPlanMode_ != "dumpOnly"
    )
    {
        offloadHdrBuf_ = List<double>(8, 0.0);
        offloadRecvBuf_ =
            List<double>
            (
                offloadPlanBytes_/sizeof(double) + 8
            );
        MPI_Irecv
        (
            offloadHdrBuf_.begin(), 6, MPI_INT, offloadPairDonor_,
            offloadTagHdr(step), MPI_COMM_WORLD, &offloadHdrRecvReq_
        );
        if (offloadPlanMode_ == "singlePairSync")
        {
            // complete the whole transfer before any local move
            MPI_Wait(&offloadHdrRecvReq_, MPI_STATUS_IGNORE);
            label* hdr = reinterpret_cast<label*>(offloadHdrBuf_.begin());
            offloadNTotal_ = hdr[3];
            offloadNChunks_ = hdr[4];
            offloadRecvBuf_ =
                List<double>
                (
                    (offloadNTotal_*sizeof(dsmcParcel::TransferData))
                   /sizeof(double) + 8
                );
            offloadChunkRecvReqs_ = List<MPI_Request>(offloadNChunks_);
            dsmcParcel::TransferData* tds =
                reinterpret_cast<dsmcParcel::TransferData*>
                (
                    offloadRecvBuf_.begin()
                );
            for (label ci = 0; ci < offloadNChunks_; ++ci)
            {
                const label begin = ci*offloadChunkParcels_;
                const label cnt =
                    min(offloadChunkParcels_, offloadNTotal_ - begin);
                MPI_Irecv
                (
                    tds + begin,
                    cnt*label(sizeof(dsmcParcel::TransferData)),
                    MPI_BYTE, offloadPairDonor_,
                    offloadTagChunk(step, ci), MPI_COMM_WORLD,
                    offloadChunkRecvReqs_.begin() + ci
                );
            }
            MPI_Waitall
            (
                offloadNChunks_,
                offloadChunkRecvReqs_.begin(),
                MPI_STATUSES_IGNORE
            );
            offloadTRecvHeader_ = offloadRelNow();
            offloadTRecvDone_ = offloadRelNow();
            offloadNRecv_ = offloadNTotal_;
        }
    }
}


void Foam::dsmcReplicatedMesh::offloadExecuteAndInsert()
{
    if (!offloadControlsRead_)
    {
        offloadReadControls();
    }
    if (!offloadEnabled_ || offloadPlanMode_ == "off")
    {
        return;
    }
    const label step = mesh_.time().timeIndex();
    const double t0 = offloadRelNow();

    // ---- dumpOnly reference: dump the same particles after the local
    // move; they were never dispatched ----
    if (offloadIAmDonor_ && offloadPlanMode_ == "dumpOnly")
    {
        if (offloadDumpTaskPtrs_.size())
        {
            std::ofstream os
            (
                "offloadTaskOut-r" + std::to_string(myRank_) + ".bin",
                std::ios::binary
            );
            const label n = offloadDumpTaskPtrs_.size();
            os.write(reinterpret_cast<const char*>(&n), sizeof(label));
            for (label i = 0; i < n; ++i)
            {
                dsmcParcel::TransferData tdr;
                offloadDumpTaskPtrs_[i]->packTransfer(tdr);
                os.write
                (
                    reinterpret_cast<const char*>(&tdr), sizeof(tdr)
                );
            }
            offloadTExecDone_ = offloadRelNow();
        }
        offloadDumpTaskPtrs_.clear();
        offloadPairsCsvLine(step, "dumpRef", offloadTDispatch_, -1.0);
        return;
    }

    if
    (
        !offloadIAmHelper_
     || offloadPlanReject_ != 0
     || offloadPlanMode_ == "dumpOnly"
    )
    {
        offloadPairsCsvLine
        (
            step,
            offloadPlanReject_ ? "reject" : "idle",
            offloadTDispatch_,
            offloadTSendDone_
        );
        return;
    }

    // ---- helper: ensure task data ----
    if (offloadPlanMode_ != "singlePairSync")
    {
        // front-posted receive + sender-side push: data normally arrived
        // during the local move; wait here for header then chunks
        MPI_Wait(&offloadHdrRecvReq_, MPI_STATUS_IGNORE);
        offloadTRecvHeader_ = offloadRelNow();
        label* hdr = reinterpret_cast<label*>(offloadHdrBuf_.begin());
        offloadNTotal_ = hdr[3];
        offloadNChunks_ = hdr[4];
        offloadRecvBuf_ =
            List<double>
            (
                (offloadNTotal_*sizeof(dsmcParcel::TransferData))
               /sizeof(double) + 8
            );
        offloadChunkRecvReqs_ = List<MPI_Request>(offloadNChunks_);
        dsmcParcel::TransferData* tds =
            reinterpret_cast<dsmcParcel::TransferData*>
            (
                offloadRecvBuf_.begin()
            );
        for (label ci = 0; ci < offloadNChunks_; ++ci)
        {
            const label begin = ci*offloadChunkParcels_;
            const label cnt =
                min(offloadChunkParcels_, offloadNTotal_ - begin);
            MPI_Irecv
            (
                tds + begin,
                cnt*label(sizeof(dsmcParcel::TransferData)),
                MPI_BYTE, offloadPairDonor_,
                offloadTagChunk(step, ci), MPI_COMM_WORLD,
                offloadChunkRecvReqs_.begin() + ci
            );
        }
        MPI_Waitall
        (
            offloadNChunks_,
            offloadChunkRecvReqs_.begin(),
            MPI_STATUSES_IGNORE
        );
        offloadTRecvDone_ = offloadRelNow();
    }
    offloadNRecv_ = offloadNTotal_;

    // ---- unpack (parallel) ----
    const label nT = max(cloud_.ompNumThreads(), label(1));
    List<dsmcParcel*> received(offloadNTotal_, nullptr);
    dsmcParcel::TransferData* tds =
        reinterpret_cast<dsmcParcel::TransferData*>
        (
            offloadRecvBuf_.begin()
        );
    #pragma omp parallel for num_threads(nT) schedule(static)
    for (label j = 0; j < offloadNTotal_; ++j)
    {
        received[j] = dsmcParcel::unpackTransfer(mesh_, tds[j]);
    }

    // ---- task move with the helper's OMP team ----
    // Safe-task path is RNG-free (no newParcel sampling, no boundary
    // models), so per-particle results are bitwise identical to local
    // execution regardless of thread scheduling.
    label violations = 0;
    #pragma omp parallel num_threads(nT) reduction(+:violations)
    {
        dsmcParcel::trackingData td(cloud_);
        td.fastRng = false;
        td.moveDetailProfile = false;
        #pragma omp for schedule(runtime)
        for (label j = 0; j < offloadNTotal_; ++j)
        {
            dsmcParcel* p = received[j];
            p->stepFraction() = 0;
            td.switchProcessor = false;
            td.keepParticle = true;
            const bool keep = p->move(td, cloud_.deltaTValue());
            if (!keep || td.switchProcessor)
            {
                ++violations;
            }
        }
    }
    offloadNExecViolation_ = violations;
    offloadTExecDone_ = offloadRelNow();

    // ---- insert results into the cloud; the regular migration routes
    // them to their post-move cell owners ----
    label inserted = 0;
    for (label j = 0; j < offloadNTotal_; ++j)
    {
        if (received[j])
        {
            cloud_.Cloud<dsmcParcel>::addParticle(received[j]);
            ++inserted;
        }
    }
    if (offloadNTotal_ > 0)
    {
        // moveOrdered/occupancy views no longer match the cloud
        cloud_.clearMoveOrderedParcels();
    }
    offloadNInserted_ = inserted;
    offloadTInsertDone_ = offloadRelNow();

    offloadPairsCsvLine(step, "helper", offloadTDispatch_, offloadTSendDone_);

    // P2 dump: post-execution records of the executed task (helper side)
    if
    (
        step == offloadStateDumpStep_
     && offloadPlanMode_ == "singlePairSync"
    )
    {
        std::ofstream os
        (
            "offloadTaskOut-r" + std::to_string(myRank_) + ".bin",
            std::ios::binary
        );
        const label n = offloadNTotal_;
        os.write(reinterpret_cast<const char*>(&n), sizeof(label));
        for (label j = 0; j < n; ++j)
        {
            if (received[j])
            {
                dsmcParcel::TransferData tdr;
                received[j]->packTransfer(tdr);
                os.write
                (
                    reinterpret_cast<const char*>(&tdr), sizeof(tdr)
                );
            }
        }
    }
}


void Foam::dsmcReplicatedMesh::offloadDumpChecksum()
{
    if
    (
        !offloadEnabled_
     || offloadDumpDone_
     || offloadStateDumpStep_ < 0
     || mesh_.time().timeIndex() != offloadStateDumpStep_
    )
    {
        return;
    }
    offloadDumpDone_ = true;
    uint64_t h = 1469598103934665603ULL;
    label n = 0;
    forAllIter(Cloud<dsmcParcel>, cloud_, iter)
    {
        h += offloadParcelHash(iter());
        ++n;
    }
    std::ofstream os
    (
        "offloadChecksum-r" + std::to_string(myRank_) + ".txt"
    );
    os << myRank_ << " " << n << " " << h << "\n";
}


} // End namespace Foam
