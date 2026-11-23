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

#ifndef dsmcReplicatedMesh_H
#define dsmcReplicatedMesh_H

#include "fvMesh.H"
#include "labelList.H"
#include "Pstream.H"
#include "dsmcLocalMesh.H"
#include <mpi.h>
#include <fstream>

namespace Foam
{
class dsmcCloud;
class dsmcParcel;

class dsmcReplicatedMesh
{
public:
    enum MigrationProfileContext
    {
        migrationProfileOther = 0,
        migrationProfileRegular = 1,
        migrationProfileDLB = 2,
        migrationProfileOutput = 3,
        migrationProfileManual = 4,
        migrationProfileInitial = 5
    };

private:
    dsmcCloud& cloud_;
    const fvMesh& mesh_;
    dsmcLocalMesh localMesh_;
    labelList cellOwner_;
    DynamicList<label> myCells_;
    label localParticleCount_;
    labelList allParticleCounts_;

    scalar migrationWallTime_;
    scalar migrationPackWallTime_;
    scalar migrationLocalPrepWallTime_;
    scalar migrationSizeExchangeWallTime_;
    scalar migrationRequestPostWallTime_;
    scalar migrationWaitWallTime_;
    scalar migrationDeserializeWallTime_;
    scalar migrationCandidateGatherWallTime_;
    scalar migrationPostWallTime_;
    scalar updateParticleCountsWallTime_;
    label migrationCalls_;
    scalar evolveStepTime_;
    label evolveTimeSteps_;

    label migrationProfileContext_;
    scalar migrationCommWallTime_;
    scalar regularMigrationWallTime_;
    scalar regularMigrationCommWallTime_;
    scalar dlbMigrationWallTime_;
    scalar dlbMigrationCommWallTime_;
    scalar outputMigrationWallTime_;
    scalar outputMigrationCommWallTime_;
    scalar outputGatherCommWallTime_;
    scalar manualMigrationWallTime_;
    scalar manualMigrationCommWallTime_;
    scalar initialMigrationWallTime_;
    scalar initialMigrationCommWallTime_;
    scalar otherMigrationWallTime_;
    scalar otherMigrationCommWallTime_;

    bool active_;
    label nProcs_;
    label myRank_;
    label migrateInterval_;
    label stepCounter_;
    labelList rebalanceSteps_;
    label rebalanceCount_;
    label totalCellsChanged_;
    label totalParcelsMigrated_;

    // ---- Phase C: auto DLB ----
    bool autoDLBEnabled_;
    scalar imbalanceThreshold_;         // threshold trigger: max/min > this
    label dlbSteps_;                    // evaluation interval (steps)
    label sarSteps_;                    // minimum accumulated steps before SAR checks
    word autoDLBTriggerMode_;           // cumulative or legacyWindow
    labelList forcedDLBSteps_;          // optional deterministic DLB trigger steps
    label lastAutoRebalanceStep_;      // step of last auto rebalance
    label autoRebalanceCount_;         // number of auto rebalances performed

    // ParDSMC3D-style trend-based trigger (SAR)
    scalar productiveTime_;            // cumulative move+coll time in eval window
    scalar lastEvolveTime_;            // evolveStepTime_ snapshot for delta
    scalar lastMigrationTime_;         // migrationWallTime_ snapshot for delta
    scalar tidl_;                      // cumulative idle time across eval periods
    scalar tdecps_;                    // cumulative decomposition wall time
    scalar w1_, w2_;                   // previous / current avg idle rate
    label ndecps_;                     // steps in current eval window
    label sarEvalSteps_;               // steps in current SAR sampling window
    label nEvalPeriods_;              // number of eval periods since last trigger
    scalar sar_;                       // trend metric: sar > 0 triggers DLB

    // DLB diagnostics
    label postDLBSnapshotCountdown_;
    scalarList lastDLBMoveT_;
    scalarList lastDLBCollT_;
    label autoRebalanceChecks_;
    label autoRebalanceTriggeredChecks_;
    scalar autoRebalanceWallTime_;
    scalar autoRebalanceCheckWallTime_;
    scalar autoRebalanceRepartWallTime_;
    scalar autoRebalanceMigrationWallTime_;
    scalar autoRebalanceWriteWallTime_;
    scalar autoRebalancePostDiagWallTime_;
    scalar autoRebalanceExplicitCommWallTime_;

    // Particle-count exponent for the DLB proxy: w = N^alpha
    scalar dlbAlpha_;
    bool adaptiveAlpha_;
    scalar adaptiveAlphaMin_;
    scalar adaptiveAlphaMax_;
    scalar adaptiveAlphaGain_;
    scalar adaptiveAlphaMaxStep_;
    scalar adaptiveAlphaWorsenTol_;
    scalar adaptiveAlphaLastImbalance_;
    scalar adaptiveAlphaLastStep_;

    // ---- Cached per-rank candidate counts (updated during migration) ----
    labelList allProcCandidates_;

    // ---- dlb_offload.md P2-P5: per-step move offload (default off) ----
    // One rank plays at most one role (donor, helper or idle) per step;
    // results are inserted into the helper's cloud before the regular
    // migration, which routes them to their post-move cell owners — no
    // dedicated result protocol.  Roles: offloadPlanMode = off | dumpOnly
    // | singlePairSync | fixedPairs | greedy.
    bool offloadEnabled_;
    word offloadPlanMode_;
    label offloadMaxPairs_;
    label offloadPairBytesCap_;
    label offloadChunkParcels_;
    scalar offloadMinGain_;
    bool offloadSkipDLBSteps_;
    label offloadStateDumpStep_;
    bool offloadDumpDone_;
    bool offloadControlsRead_;
    label offloadPlanId_;
    label offloadLastRebalanceCount_;
    bool offloadPairsFileOpen_;
    std::ofstream offloadPairsFile_;

    // Owned buffers: one send stream (donor role) and one receive
    // stream (helper role); reused across steps, sized to the plan cap
    List<double> offloadSendBuf_;
    List<double> offloadRecvBuf_;
    List<double> offloadHdrBuf_;
    List<dsmcParcel*> offloadDumpTaskPtrs_;
    MPI_Request offloadHdrSendReq_;
    List<MPI_Request> offloadChunkSendReqs_;
    MPI_Request offloadHdrRecvReq_;
    List<MPI_Request> offloadChunkRecvReqs_;
    bool offloadSendsPending_;

    //- Donor-side EMA of dispatch seconds per selected parcel, fed into
    //  the planner so tasks are only priced when net-positive
    double offloadDispatchRate_;

    // Per-step pair view of THIS rank (0 or 1 donor-role + 0 or 1
    // helper-role pair; roles are exclusive per the revised §4).
    bool offloadIAmDonor_;
    bool offloadIAmHelper_;
    label offloadPairDonor_;
    label offloadPairHelper_;
    label offloadPlanNPairs_;
    label offloadPlanReject_;
    label offloadNSelected_;
    label offloadNRecv_;
    label offloadNExecViolation_;
    label offloadNInserted_;
    label offloadPlanBytes_;
    label offloadNTotal_;      // header: actual parcels in this task
    label offloadNChunks_;     // header: actual chunk count

    scalar offloadTDispatch_, offloadTSendDone_, offloadTRecvHeader_;
    scalar offloadTRecvDone_, offloadTExecDone_, offloadTInsertDone_;

    // ---- Async migration state ----
    bool asyncMigrationPending_;
    List<char> asyncSendBuf_;          // kept alive until Waitall
    List<char> asyncRecvBuf_;
    label asyncRecvSize_;
    // clean N1: per-peer received byte counts.  In the multi-rank non-flat
    // path asyncRecvBuf_ is a concatenation of per-peer streams, each with
    // its own <parcel count> header, so the reader must split it by these
    // sizes (reading a single count would drop the other peers' records).
    labelList asyncRecvSizes_;
    DynamicList<MPI_Request> asyncReqs_;

    // ---- Phase 2: per-peer predicted receive buffers (no Alltoall) ----
    labelList prevRecvSizes_;             // previous step's recv bytes per peer
    labelList perPeerRecvCapacity_;       // allocated capacity per peer
    List<List<char>> perPeerRecvBuf_;     // per-peer receive buffers
    bool useNoAlltoall_;                  // controlDict flag
    bool useFlatTransfer_;                // flat POD buffer serialization
    label transferChunkBytes_;            // maximum bytes per MPI transfer chunk
    bool gatherCandidates_;               // gather per-rank collision candidates
    bool overlapSizeExchange_;            // overlap Alltoall sizes with local prep
    // §36 A/B switches (default true = the §36 behaviour).  Both default
    // on; set false in controlDict to emulate the pre-§36 code path so the
    // same binary can measure each item against its own control.
    bool migTakeOrderedParcels_;          // §36.1 zero-copy ordered table
    bool migReuseMigBuffers_;             // §36.2 persistent migrate buffers
    word writeMode_;                      // gathered or processor
    mutable label processorWriteDecompVersion_;
    mutable word processorWriteMeshInstance_;

    // ---- Wall-time trigger metric (worklog dlb_further_opt.md §5.6) ----
    bool triggerWallTime_;                // replicatedMeshDLBTriggerWallTime
    scalar lastStepWallTime_;             // full wall of the previous iteration
                                          // (set by the solver main loop)
    scalar stepFullWallCum_;              // cumulative per-rank full iteration wall
    scalar lastStepWait_;                 // collective-wait delta of the previous
                                          // accounted step (paired with
                                          // lastStepWallTime_)
    scalar lastSizeExchangeWall_;         // wait-accumulator snapshots
    scalar lastWaitWall_;
    scalar lastUpcWall_;

    // ---- Futile-rebalance back-off (worklog dlb_further_opt.md §5.3) ----
    bool backOffEnabled_;                 // replicatedMeshDLBBackOff
    scalar backOffChangedRatio_;          // eps: ownership-change ratio floor
    scalar backOffTol_;                   // imbalance improvement tolerance
    label backOffFutileRuns_;             // consecutive futile rebalances
    label backOffSkipChecks_;             // checks to skip once in back-off
    scalar backOffEscalate_;              // safety valve (x threshold)
    label backOffRemaining_;              // checks remaining in back-off
    label futileCount_;
    scalar lastExecChangedRatio_;         // nChanged/nCells of last executed DLB
    scalar lastExecImbalance_;            // imbalance that triggered last DLB

    // ---- §45 P0: per-cell cost 观测（worklog §45） ----
    word cellCostMode_;                   // off | observe | weighted (§45 P3)
    label cellCostDumpEvery_;             // 每 N 个 DLB 事件 dump 一次（1=每次）
    label cellCostDumpCount_;             // 已 dump 窗口数
    label cellCostEventCount_;            // 已发生 DLB 事件数（含未 dump 事件）
    scalarList cellCostCoeffs_;           // §45 P3: weighted vwgt 系数 (A,B,C)
    // §45 P2: 受控迁移干预入口（observe-only 实验用，默认空=关闭）
    fileName cellCostInterveneFile_;      // owner 覆盖文件（"cellIdx owner" 两列）
    label cellCostInterveneEvent_;        // 第几个 DLB 事件应用（默认 2）
    bool cellCostInterveneDone_;          // 已应用标志
    List<label> cellCostInterveneOwner_;  // 覆盖表（全局 cell 号 → owner）
    List<scalar> cellCostQ_;              // 窗口内每 cell 候选和
    List<scalar> cellCostP_;              // 事件帧 per-cell 粒子快照（诊断）
    List<label> cellCostQtouched_;        // 触及过的 cell（掩码去重后追加）
    List<bool> cellCostTouchedMask_;      // 本窗口触及掩码（防列表无界增长）
    // ---- §44 Route A: trend-aware rank-k (worklog §44) ----
    word trendMode_;                      // off | rankK（cellParticles 路线已回退，§44.14）
    label trendWindowSteps_;              // 观测窗口步数
    label trendHorizonSteps_;             // 外推步数 H
    scalar trendTheta_;                   // rank 相对工作 EMA 系数
    scalar trendGainGd_;                  // 趋势增益 Gd
    scalar trendRelaxation_;              // k 更新阻尼 lambda
    scalar trendMaxDeltaK_;               // 每窗口最大 k 改变量
    scalar trendResetChangedRatio_;       // owner 变更超过该比例才重置趋势历史
    bool rankKActive_;                    // 权重构造是否使用 per-rank k
    scalar trendWorkSum_;                 // 本窗口有效工作（事件步除外）
    label trendSamples_;                  // 本窗口有效样本数
    label trendLastWindowStep_;           // 上次窗口端点步号
    bool trendSkipSample_;                // 事件步标记：跳过本步样本
    bool trendPrevBarEValid_;             // bar_e 历史有效（趋势项就绪条件之一）
        // 重置后完成的窗口数
    List<scalar> trendPrevBarE_;          // 上一窗口 bar_e
    List<scalar> trendK_;                 // 当前 per-rank k（跨窗口持续）
    bool workFeedbackK_;                  // §20 比例反馈（TrendMode=off 时生效）
    label cellWeightK_;                   // 基础 cell 权重 k（§17/§42）

    // (Phase B) construct myCells_ from cellOwner_
    void rebuildMyCells();

    // Per-rank step-wall diagnostic at check time (§5.6)
    void printPerRankStepWall() const;

    // Verify that the replicated owner map and this rank's local view agree.
    // This is intentionally called only after an executed DLB repartition.
    void validateCellOwnerMap(const char* context) const;

    void addMigrationProfileSample(const scalar wallTime, const scalar commTime);
    void addAutoRebalanceCommTime(const scalar commTime);

    // Phase C: ParMETIS AdaptiveRepart — parallel incremental repartitioning
    // that minimises data movement while rebalancing vertex weights.
    //- comm_further_opt.md §7: checkImbalance selects the ParMETIS mode
    //- (PartKway when > hybridThreshold, AdaptiveRepart otherwise)
    label reassignByParMetisAdaptiveRepart(const scalar checkImbalance = -1.0);

    // Restore a processor-write checkpoint selected by dsmcFoam+ at startup.
    bool restoreProcessorCheckpoint();

public:
    dsmcReplicatedMesh(dsmcCloud& cloud, const fvMesh& mesh);
    ~dsmcReplicatedMesh();

    inline const labelList& cellOwner() const { return cellOwner_; }
    inline bool isMyCell(const label cellI) const
        { return cellOwner_[cellI] == myRank_; }
    inline const DynamicList<label>& myCells() const { return myCells_; }
    inline bool active() const { return active_; }
    inline label myRank() const { return myRank_; }
    inline label nProcs() const { return nProcs_; }
    inline const labelList& allProcCandidates() const { return allProcCandidates_; }
    inline const dsmcLocalMesh& localMesh() const { return localMesh_; }
    inline dsmcLocalMesh& localMesh() { return localMesh_; }

    void initialize();

    // comm_further_opt.md §2: flat-path async migration overlap
    //- Whether the flat POD transfer path is active (production default)
    inline bool asyncOverlapCapable() const { return useFlatTransfer_; }
    //- Whether the async begin/gap-build/finish split is enabled
    //  (controlDict "replicatedMeshAsyncMigration", default false).
    //  Non-inline, state cached in a file-static: no data member added.
    bool asyncOverlapEnabled() const;
    //- Mark that the caller completed the gap-time occupancy build between
    //  migrateBegin() and migrateFinish() (enables the incremental append)
    void markAsyncOverlapGapBuild();

    // Phase A: static cellOwner decomposition
    void computeCellOwnerScotch();
    // Phase A: initial particle distribution (delete non-owned, no MPI)
    void distributeInitialParticles();
    // Phase A: particle migration (synchronous)
    void migrateParticlesByCellOwner();
    // Phase A: async migration — begin (serialize + Isend/Irecv, non-blocking)
    void migrateBegin();
    // Phase A: async migration — finish (Waitall + deserialize + append)
    void migrateFinish();
    // comm_further_opt.md §2: flat-path async overlap halves (private)
    void flatMigrateBegin();
    void flatMigrateFinish();
    // Phase A: particle count sync
    void updateParticleCounts();

    inline bool asyncMigrationPending() const { return asyncMigrationPending_; }

    // Phase B: manual cellOwner reassignment (block/reverse-block)
    void reassignCellOwner();

    // Phase C: automatic DLB — called once per step; checks trigger & executes
    void autoRebalance();
    void trendRankKUpdate(const label currentStep);   // §44 Route A
    void accumulateCellCostCandidates(const labelList& counts, const labelList& touched);   // §45 P0
    void dumpCellCostWindow
    (
        const label currentStep,
        const labelList& globalCellParticles
    );                                                        // §45 P0
    inline bool cellCostObserve() const { return cellCostMode_ == "observe"; }
    void writeCellOwner() const;
    void writeGatheredCloudOnRank0() const;
    bool processorWriteEnabled() const;
    void writeProcessorOutput() const;

    // ---- inline accessors ----
    inline scalar migrationWallTime() const { return migrationWallTime_; }
    inline label migrationCalls() const { return migrationCalls_; }
    inline label migrateInterval() const { return migrateInterval_; }
    inline label stepCounter() const { return stepCounter_; }
    inline void advanceStepCounter() { ++stepCounter_; }

    // dlb_offload.md P0: per-step migration stamps for the cloud probe
    // file — seconds relative to the cloud probe epoch (cloud epoch is
    // shared so cloud- and migration-side stamps join), set in
    // migrateParticlesByCellOwner when the probe is enabled.  Written by
    // the cloud into its per-rank CSV; public data to keep the probe out
    // of the mesh's hot-path call graph.

    double offPrep1_;
    double offSx1_;
    double offMig1_;
    label offSent_;
    label offRecv_;
    double offSentB_;
    double offRecvB_;

    // ---- dlb_offload.md P2-P5 methods ----

    //- Predicted own work of this rank for planning (last step wall
    //  minus collective wait), seconds
    inline scalar offloadPredictWork() const
    {
        return max(lastStepWallTime_ - lastStepWait_, scalar(0));
    }

    //- Read offload* controlDict keys (lazily on first access; the
    //  cloud-side hooks gate on this, so it must self-initialise)
    inline bool offloadEnabled()
    {
        if (!offloadControlsRead_)
        {
            offloadReadControls();
        }
        return offloadEnabled_;
    }

    void offloadReadControls();

    //- Planning + front dispatch (donor packs/sends/deletes; helper
    //  posts receives; singlePairSync completes transfer here)
    void offloadDispatch();

    //- Helper side: wait for task data, execute the task move with the
    //  OMP team, insert results into the cloud before the regular
    //  migration; donor/reference side dumps as configured
    void offloadExecuteAndInsert();

    //- Per-rank order-independent particle-multiset checksum at the
    //  dump step (P2 off/on comparison), written to a per-rank file
    void offloadDumpChecksum();

    //- Seconds relative to the cloud probe epoch (shared clock)
    double offloadRelNow() const;

    //- One per-rank per-step line of the §6 offload CSV
    void offloadPairsCsvLine
    (
        const label step,
        const char* role,
        const double tDispatch,
        const double tSendDone
    );

    inline const labelList& rebalanceSteps() const { return rebalanceSteps_; }
    inline label rebalanceCount() const { return rebalanceCount_; }
    inline void addEvolveTime(scalar dt) { evolveStepTime_ += dt; ++evolveTimeSteps_; }
    inline bool autoDLBEnabled() const { return autoDLBEnabled_; }
    inline label autoRebalanceCount() const { return autoRebalanceCount_; }
    inline scalar commDLBWallTime() const
    {
        return dlbMigrationCommWallTime_ + autoRebalanceExplicitCommWallTime_;
    }
    inline scalar commOutputWallTime() const
    {
        return outputMigrationCommWallTime_ + outputGatherCommWallTime_;
    }
    inline scalar commTotalWallTime() const
    {
        return regularMigrationCommWallTime_
             + commDLBWallTime()
             + commOutputWallTime()
             + manualMigrationCommWallTime_
             + initialMigrationCommWallTime_
             + otherMigrationCommWallTime_;
    }
    inline scalar dlbNonCommWallTime() const
    {
        const scalar nonComm =
            autoRebalanceWallTime_
          - autoRebalanceExplicitCommWallTime_
          - dlbMigrationCommWallTime_;
        return nonComm > scalar(0.0) ? nonComm : scalar(0.0);
    }
    inline label migrationProfileContext() const { return migrationProfileContext_; }
    // §28: evolve-residual attribution (read by the cloud profile).
    inline scalar updateParticleCountsWallTime() const
    {
        return updateParticleCountsWallTime_;
    }
    // §33: per-rank residual table (DLB migration wall, total).
    inline scalar autoRebalanceMigrationWallTime() const
    {
        return autoRebalanceMigrationWallTime_;
    }
    inline void setMigrationProfileContext(const label context)
    {
        migrationProfileContext_ = context;
    }
    inline void setLastStepWallTime(const scalar t)
    {
        lastStepWallTime_ = t;
        // §17.5 fix: pair the collective wait with the SAME iteration as
        // the wall.  The wait delta was previously computed inside
        // autoRebalance (mid-iteration), so a DLB executed later in that
        // iteration deposited its migration-storm wait into the NEXT
        // step's delta: post-DLB step 1 showed the storm wall as "work"
        // (~30s) and step 2 clamped to 0 (storm wait > normal wall).
        // This setter is called at the end of each iteration, after every
        // collective, so wall(N) and wait(N) pair exactly.
        lastStepWait_ =
            (migrationSizeExchangeWallTime_ - lastSizeExchangeWall_)
          + (migrationWaitWallTime_ - lastWaitWall_)
          + (updateParticleCountsWallTime_ - lastUpcWall_);
        lastSizeExchangeWall_ = migrationSizeExchangeWallTime_;
        lastWaitWall_ = migrationWaitWallTime_;
        lastUpcWall_ = updateParticleCountsWallTime_;
        // §44 Route A: 趋势窗口采样（事件步——本步内执行过重分区——整步剔除，
        // 样本数一并扣除，§44.2）。wall/wait 配对与 §17.5 相同（一步滞后）。
        if (trendMode_ == "rankK")
        {
            if (trendSkipSample_)
            {
                trendSkipSample_ = false;
            }
            else
            {
                trendWorkSum_ += max(t - lastStepWait_, scalar(0));
                ++trendSamples_;
            }
        }
    }
    inline void addOutputGatherCommTime(const scalar dt)
    {
        outputGatherCommWallTime_ += dt;
    }

    void gatherParcelsToRank0();

    // §34.5: post the particle-count gather at the migration tail so its
    // collective convergence overlaps the ordered-view rebuild; the
    // cloud's updateParticleCounts() completes it.
    void postParticleCountsGather();

    void report() const;
};

} // End namespace Foam
#endif
