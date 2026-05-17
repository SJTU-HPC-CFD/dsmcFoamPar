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
    scalarField cellCost_;
    label cellCostSteps_;
    scalarField cellCollisionCost_;
    label cellCollisionCostSteps_;
    scalar costMoveWeight_;
    scalar costCollisionWeight_;
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
    label nSuperCells_;
    label alpha_;                       // over-decomposition factor M/P
    scalar imbalanceThreshold_;         // trigger: accumulated idle > threshold
    label cooldownSteps_;              // min steps between rebalances
    label lastAutoRebalanceStep_;      // step of last auto rebalance
    label autoRebalanceCount_;         // number of auto rebalances performed
    scalarField superCellCost_;        // cost per superCell
    labelList superCellOfCell_;        // cell -> superCell mapping
    labelList superCellOwner_;         // superCell -> rank mapping
    vectorField superCellCentres_;     // centroid per superCell

    // ParDSMC3D-style trend-based trigger
    scalar productiveTime_;            // cumulative move+coll time in eval window
    scalar lastEvolveTime_;            // evolveStepTime_ snapshot for delta
    scalar lastMigrationTime_;         // migrationWallTime_ snapshot for delta
    scalar tidl_;                      // cumulative idle time across eval periods
    scalar tdecps_;                    // cumulative decomposition wall time
    scalar w1_, w2_;                   // previous / current avg idle rate
    label ndecps_;                     // steps in current eval window
    label nEvalPeriods_;              // number of eval periods since last trigger
    scalar sar_;                       // trend metric: sar > 0 triggers DLB

    // ---- Async migration state ----
    bool asyncMigrationPending_;
    List<char> asyncSendBuf_;          // kept alive until Waitall
    List<char> asyncRecvBuf_;
    label asyncRecvSize_;
    DynamicList<MPI_Request> asyncReqs_;

    // (Phase B) construct myCells_ from cellOwner_
    void rebuildMyCells();

    // Hilbert curve helper
    static uint64_t hilbert3DIndex
    (
        uint32_t x, uint32_t y, uint32_t z,
        label bits
    );

    // Phase C: build superCells by Scotch over-decomposition
    void buildSuperCells();

    // Phase C: compute superCell centroids and costs from TACF data
    void updateSuperCellCosts();

    // Phase C: sort superCells by Hilbert key, assign to ranks by cost balance
    // returns number of superCells that changed owner
    label reassignByHilbertSuperCell();

    // Phase C: SCOTCH graphRemap — directly remaps cellOwner_ using
    // per-cell TACF costs.  Minimises data movement while rebalancing.
    void reassignByScotchRemap();

    // Phase C: ParMETIS AdaptiveRepart — parallel incremental repartitioning
    // that minimises data movement while rebalancing vertex weights.
    void reassignByParMetisAdaptiveRepart();

public:
    dsmcReplicatedMesh(dsmcCloud& cloud, const fvMesh& mesh);
    ~dsmcReplicatedMesh();

    inline const labelList& cellOwner() const { return cellOwner_; }
    inline bool isMyCell(const label cellI) const
        { return cellOwner_[cellI] == myRank_; }
    inline const DynamicList<label>& myCells() const { return myCells_; }
    inline bool active() const { return active_; }
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

    // Phase A/B: TACF cell cost accumulation
    void accumulateCellCosts();
    // Phase B: manual cellOwner reassignment (block/reverse-block)
    void reassignCellOwner();

    // Phase C: automatic DLB — called once per step; checks trigger & executes
    void autoRebalance();

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

    void report() const;
};

} // End namespace Foam
#endif
