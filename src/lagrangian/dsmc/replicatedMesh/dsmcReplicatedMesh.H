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

namespace Foam
{
class dsmcCloud;

class dsmcReplicatedMesh
{
    dsmcCloud& cloud_;
    const fvMesh& mesh_;
    dsmcLocalMesh localMesh_;
    labelList cellOwner_;
    DynamicList<label> myCells_;
    label localParticleCount_;
    labelList allParticleCounts_;

    scalar migrationWallTime_;
    label migrationCalls_;
    scalar evolveStepTime_;
    label evolveTimeSteps_;

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

    // Adaptive K (persistent across DLB calls)
    label collDivisor_;

    // Adaptive alpha for weight formula w = N^alpha
    scalar alpha_;

    // ---- Cached per-rank candidate counts (updated during migration) ----
    labelList allProcCandidates_;

    // ---- Async migration state ----
    bool asyncMigrationPending_;
    List<char> asyncSendBuf_;          // kept alive until Waitall
    List<char> asyncRecvBuf_;
    label asyncRecvSize_;
    DynamicList<MPI_Request> asyncReqs_;

    // ---- Phase 2: per-peer predicted receive buffers (no Alltoall) ----
    labelList prevRecvSizes_;             // previous step's recv bytes per peer
    labelList perPeerRecvCapacity_;       // allocated capacity per peer
    List<List<char>> perPeerRecvBuf_;     // per-peer receive buffers
    bool useNoAlltoall_;                  // controlDict flag
    bool useFlatTransfer_;                // flat POD buffer serialization
    word writeMode_;                      // gathered or processor

    // (Phase B) construct myCells_ from cellOwner_
    void rebuildMyCells();

    // Phase C: ParMETIS AdaptiveRepart — parallel incremental repartitioning
    // that minimises data movement while rebalancing vertex weights.
    label reassignByParMetisAdaptiveRepart();

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
    // Phase A: particle count sync
    void updateParticleCounts();

    inline bool asyncMigrationPending() const { return asyncMigrationPending_; }

    // Phase B: manual cellOwner reassignment (block/reverse-block)
    void reassignCellOwner();

    // Phase C: automatic DLB — called once per step; checks trigger & executes
    void autoRebalance();
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
    inline const labelList& rebalanceSteps() const { return rebalanceSteps_; }
    inline label rebalanceCount() const { return rebalanceCount_; }
    inline void addEvolveTime(scalar dt) { evolveStepTime_ += dt; ++evolveTimeSteps_; }
    inline bool autoDLBEnabled() const { return autoDLBEnabled_; }
    inline label autoRebalanceCount() const { return autoRebalanceCount_; }

    void gatherParcelsToRank0();

    void report() const;
};

} // End namespace Foam
#endif
