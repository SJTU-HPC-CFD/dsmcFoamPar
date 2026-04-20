/*---------------------------------------------------------------------------*\\
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

Class
    Foam::dsmcCloud

Description
    Minimal v2506 DSMC cloud extended with the subset of stage-2 features
    required by the orion107kmNR case.

SourceFiles
    dsmcCloudI.H
    dsmcCloud.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCloud_H
#define dsmcCloud_H

#include "Cloud.H"
#include "IOdictionary.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "dsmcParcel.H"
#include "BinaryCollisionModel.H"
#include "dsmcReactions.H"
#include "dsmcCoordinateSystem.H"
#include "cellMeasurements.H"
#include "dsmcControllers.H"
#include "boundaryMeasurements.H"
#include "dsmcFieldProperties.H"
#include "dsmcBoundaries.H"
#include "dsmcFaceTracker.H"
#include "collisionPartnerSelection.H"
#include "porousMeasurements.H"

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace Foam
{

class dsmcCloud
:
    public Cloud<dsmcParcel>
{
    struct diffuseWallSpec
    {
        word patchName;
        label patchId = -1;
        vector velocity = vector::zero;
        scalar temperature = 0.0;
    };

    struct inflowSpec
    {
        word patchName;
        label patchId = -1;
        labelList speciesIds;
        vector velocity = vector::zero;
        scalar translationalTemperature = 0.0;
        scalar rotationalTemperature = 0.0;
        scalar vibrationalTemperature = 0.0;
        scalar electronicTemperature = 0.0;
        scalarList numberDensities;
        labelList faces;
        labelList cells;
        List<scalarField> accumulatedParcelsToInsert;
    };

    struct fieldSpec
    {
        word fieldName;
        labelList speciesIds;
        bool measureMeanFreePath = false;
        bool resetAtOutput = false;
        scalar resetAtOutputUntilTime = -GREAT;
        bool pendingReset = false;
    };

    word cloudName_;
    const fvMesh& mesh_;
    IOdictionary particleProperties_;
    List<word> typeIdList_;
    autoPtr<dsmcCoordinateSystem> dsmcCoordinateSystem_;
    scalar nParticle_;
    label nTerminalOutputs_;
    word collisionPartnerSelectionModel_;
    List<DynamicList<dsmcParcel*>> cellOccupancy_;
    volScalarField sigmaTcRMax_;
    scalarField collisionSelectionRemainder_;
    List<dsmcParcel::constantProperties> constProps_;
    Random rndGen_;
    bool openmpEnabled_;
    bool openmpMoveEnabled_;
    bool trackerActive_;
    label ompNumThreads_;
    word openmpCollisionStrategy_;
    label openmpCollisionChunk_;
    word openmpMoveSchedule_;
    label openmpMoveChunk_;
    label openmpMoveGuardLayers_;
    boolList openmpMoveGuardCells_;
    label openmpMoveGuardCellCount_;
    bool openmpAdaptivePartition_;
    label openmpPartitionMinUpdateInterval_;
    scalar openmpPartitionImbalanceThreshold_;
    scalar openmpPartitionDriftThreshold_;
    scalar openmpCollisionCostCandidateWeight_;
    scalar openmpCollisionCostActiveCellWeight_;
    bool collisionProfileEnabled_;
    bool evolveProfileEnabled_;
    bool profilingDetailEnabled_;
    bool emitStepDiagnostics_;
    List<Random> ompRndGens_;
    labelList particleLoadStart_;
    labelList particleLoadEnd_;
    label particlePartitionStep_;
    label particlePartitionLastRebuildStep_;
    labelList particlePartitionReferenceLoad_;
    List<dsmcParcel*> moveOrderedParcels_;
    labelList moveOrderedThreadOffsets_;
    bool moveOrderedParcelsValid_;
    bool moveAppendCaptureActive_;
    DynamicList<dsmcParcel*> moveAppendedParcels_;
    DynamicList<dsmcParcel*> pendingMoveParcels_;
    DynamicList<dsmcParcel*> occupancyOrderedParcels_;
    labelList occupancyCellOffsets_;
    labelList occupancyActiveCells_;
    labelList occupancyCollisionCells_;
    bool occupancyOrderedParcelsValid_;
    mutable bool cellOccupancyMaterialized_;
    List<labelList> occupancyThreadCellCounts_;
    List<DynamicList<label>> occupancyThreadActiveCells_;
    scalarField selectedPairsPerCell_;
    labelList nCandidatesPerCell_;
    labelList collisionCandidateCells_;
    labelList collisionLoadStart_;
    labelList collisionLoadEnd_;
    label collisionPartitionStep_;
    label collisionPartitionLastRebuildStep_;
    scalarField collisionPartitionReferenceLoad_;
    scalar collisionPrecomputeWallTime_;
    scalar collisionPartitionWallTime_;
    scalar collisionSelectionWallTime_;
    label collisionProfileCalls_;
    scalar buildOccupancyExtractWallTime_;
    scalar buildOccupancyCountWallTime_;
    scalar buildOccupancyAssembleWallTime_;
    label buildOccupancyProfileCalls_;
    label buildOccupancyMoveOrderedHits_;
    label buildOccupancyFallbackHits_;
    label buildOccupancyMoveOrderedInvalidHits_;
    label buildOccupancyMoveOrderedOffsetMismatchHits_;
    label buildOccupancyMoveOrderedSizeMismatchHits_;
    scalar buildOccupancyMoveOrderedSizeDeltaSum_;
    label buildOccupancyMoveOrderedSizeDeltaMin_;
    label buildOccupancyMoveOrderedSizeDeltaMax_;
    scalar buildOccupancyMoveOrderedParcelsSum_;
    scalar buildOccupancyMoveAppendedParcelsSum_;
    scalar buildOccupancyCloudSizeSum_;
    scalar buildOccupancyCloudValidCellSum_;
    scalar buildOccupancyCloudInvalidCellSum_;
    scalar buildOccupancyMoveOrderedValidCellSum_;
    scalar buildOccupancyMoveOrderedInvalidCellSum_;
    scalar buildOccupancyPartitionAssignedParcelsSum_;
    scalar buildOccupancyPartitionGapCellsSum_;
    scalar buildOccupancyPartitionOverlapCellsSum_;
    scalar movePreControlWallTime_;
    scalar moveResetSetupWallTime_;
    scalar moveExtractWallTime_;
    scalar moveExtractDeferredWallTime_;
    scalar moveExtractOrderedReuseWallTime_;
    scalar moveExtractFullScanWallTime_;
    scalar moveKernelWallTime_;
    scalar moveCommitWallTime_;
    scalar moveTransferFinalizeWallTime_;
    label moveProfileCalls_;
    label moveCommitCountCalls_;
    label moveExtractDeferredPasses_;
    label moveExtractOrderedReusePasses_;
    label moveExtractFullScanPasses_;
    scalar moveExtractedParcelsSum_;
    scalar moveExtractDeferredParcelsSum_;
    scalar moveExtractOrderedReuseParcelsSum_;
    scalar moveExtractFullScanParcelsSum_;
    label moveFirstPassReuseChecks_;
    label moveFirstPassHasOrderedCount_;
    label moveFirstPassMissingOrderedCount_;
    label moveFirstPassOffsetMismatchCount_;
    label moveFirstPassSizeMatchCount_;
    label moveFirstPassPriorPlusPendingMatchCount_;
    scalar moveFirstPassCurrentSizeSum_;
    scalar moveFirstPassPriorSizeSum_;
    scalar moveFirstPassPendingSizeSum_;
    scalar moveFirstPassSizeDeltaSum_;
    scalar moveFirstPassDeltaMinusPendingSum_;
    label previousMoveEndCloudSize_;
    scalar moveInterStepCloudDeltaSum_;
    label moveInterStepCloudDeltaMin_;
    label moveInterStepCloudDeltaMax_;
    scalar movePostStepCloudDeltaSum_;
    label movePostStepCloudDeltaMin_;
    label movePostStepCloudDeltaMax_;
    scalar moveSurvivorParcelsSum_;
    scalar moveTransferredParcelsSum_;
    scalar moveDeletedParcelsSum_;
    scalar moveLoopPassesSum_;
    label moveLoopPassesMax_;
    scalar moveDeferredParcelsSum_;
    label moveDeferredParcelsMax_;
    scalar moveReceivedParcelsSum_;
    label moveReceivedParcelsMax_;
    scalar moveTrackWallTime_;
    scalar moveTrackerCallbackWallTime_;
    scalar moveBoundaryControlWallTime_;
    scalar moveFaceHitsSum_;
    scalar moveCyclicHitsSum_;
    scalar moveStuckHitsSum_;
    scalar movePatchHitsSum_;
    scalar moveProcessorHitsSum_;
    scalar evolveMoveWallTime_;
    scalar evolveBuildWallTime_;
    scalar evolveCoordWallTime_;
    scalar evolveCollisionWallTime_;
    scalar evolveReactionWallTime_;
    scalar evolvePostWallTime_;
    label evolveProfileCalls_;
    labelList moveThreadParticleCounts_;
    scalarField moveThreadWallTimes_;
    scalarField moveThreadTrackWallTimes_;
    scalarField moveThreadTrackerWallTimes_;
    scalarField moveThreadBoundaryWallTimes_;
    labelList moveThreadFaceHitCounts_;
    labelList moveThreadCyclicHitCounts_;
    labelList moveThreadStuckHitCounts_;
    labelList moveThreadPatchHitCounts_;
    labelList moveThreadProcessorHitCounts_;
    labelList moveLastThreadParticleCounts_;
    scalarField moveLastThreadWallTimes_;
    labelList collisionThreadCandidateCounts_;
    labelList collisionThreadAcceptedCounts_;
    labelList collisionThreadActiveCellCounts_;
    labelList collisionThreadReactionHitCounts_;
    scalarField collisionThreadWallTimes_;
    autoPtr<porousMeasurements> porousMeasurements_;
    dsmcControllers controllers_;
    boundaryMeasurements boundaryMeas_;
    dsmcFieldProperties fields_;
    dsmcBoundaries boundaries_;
    dsmcFaceTracker trackingInfo_;
    autoPtr<BinaryCollisionModel> binaryCollisionModel_;
    autoPtr<collisionPartnerSelection> collisionPartnerSelectionPtr_;
    autoPtr<dsmcReactions> reactions_;
    cellMeasurements cellMeas_;

    List<diffuseWallSpec> diffuseWalls_;
    DynamicList<label> specularWallPatchIds_;
    labelList wallPatchToModelId_;
    List<inflowSpec> inflows_;

    List<fieldSpec> fieldSpecs_;
    PtrList<volScalarField> dsmcNFields_;
    PtrList<volScalarField> rhoNFields_;
    PtrList<volScalarField> rhoMFields_;
    PtrList<volScalarField> pFields_;
    PtrList<volScalarField> TtraFields_;
    PtrList<volVectorField> UMeanFields_;
    PtrList<volVectorField> fDFields_;
    List<List<vectorField>> speciesBoundaryForceDensity_;
    List<label> fieldAverageSteps_;
    List<scalarField> rhoNCumulative_;
    List<scalarField> rhoMCumulative_;
    List<vectorField> momentumCumulative_;
    List<scalarField> kineticCumulative_;

    void buildConstProps();
    void buildCellOccupancy(const bool rebuildParticlePartition = true);
    void materializeCellOccupancy();
    void readBoundaryModels();
    void readFieldSpecs();
    void createFields();
    void resetBoundaryForceDensity();
    void initOpenMPMoveGuardCells();
    void collisions();
    void noTimeCounterCollisions();
    void insertInflowParcels();
    void sampleFields();
    void initialiseSpeciesBoundaryForceDensity();
    void initOpenMP();
    void refreshTrackerUsage();
    void precomputeCollisionCandidates();
    void rebuildCollisionLoadPartition();
    void rebuildParticleLoadPartition();
    void rebuildParticleLoadPartition
    (
        const labelList& currentCellLoads,
        const label totalParticles
    );
    scalar initialSigmaTcRMax() const;
    bool speciesSelected(const label fieldI, const label typeId) const;

    dsmcCloud(const dsmcCloud&);
    void operator=(const dsmcCloud&);

public:

    dsmcCloud
    (
        Time& t,
        const word& cloudName,
        const fvMesh& mesh,
        bool readFields = true
    );

    virtual ~dsmcCloud();

    inline const word& cloudName() const;
    inline const fvMesh& mesh() const;
    inline const IOdictionary& particleProperties() const;
    inline const List<word>& typeIdList() const;
    inline const dsmcCoordinateSystem& coordSystem() const;
    inline dsmcCoordinateSystem& coordSystem();
    inline const porousMeasurements& porousMeas() const;
    inline porousMeasurements& porousMeas();
    inline scalar nParticle() const;
    inline scalar nParticles(const label celli = 0) const;
    inline scalar nParticles(const label patchi, const label facei) const;
    inline scalar deltaTValue(const label celli = 0) const;
    inline label nTerminalOutputs() const;
    inline bool openmpEnabled() const;
    inline bool openmpMoveEnabled() const;
    inline bool trackerActive() const;
    inline bool evolveProfileEnabled() const;
    inline bool profilingDetailEnabled() const;
    inline void setStepDiagnosticOutput(const bool enabled);
    inline bool emitStepDiagnostics() const;
    inline label ompNumThreads() const;
    inline const word& openmpCollisionStrategy() const;
    inline label openmpCollisionChunk() const;
    inline const word& openmpMoveSchedule() const;
    inline label openmpMoveChunk() const;
    inline bool openmpMoveGuardCell(const label celli) const;
    inline label currentThreadId() const;
    inline const labelList& particleLoadStart() const;
    inline const labelList& particleLoadEnd() const;
    void beginMoveAppendCapture();
    void endMoveAppendCapture();
    void recordMoveAppendedParcel(dsmcParcel* pPtr);
    void clearMoveAppendedParcels();
    void recordPendingMoveParcel(dsmcParcel* pPtr);
    void clearPendingMoveParcels();
    void recordMoveCommitCounts
    (
        const label extracted,
        const label survivors,
        const label transferred,
        const label deleted
    );
    void recordMoveLoopPasses(const label nPasses);
    void recordMoveDeferredParcels(const label nDeferred);
    void recordMoveReceivedParcels(const label nReceived);
    void recordMoveExtractDetail
    (
        const scalar deferredWallTime,
        const scalar orderedReuseWallTime,
        const scalar fullScanWallTime,
        const label deferredParcels,
        const label orderedReuseParcels,
        const label fullScanParcels
    );
    void recordMoveFirstPassReuseCheck
    (
        const bool hasOrdered,
        const bool offsetsOk,
        const label currentSize,
        const label priorSize,
        const label appendedSize
    );
    void storeMoveOrderedParcels(const List<dsmcParcel*>& parcels, const labelList& threadOffsets);
    void clearMoveOrderedParcels();
    inline bool hasMoveOrderedParcels() const;
    inline const List<dsmcParcel*>& moveOrderedParcels() const;
    inline const labelList& moveOrderedThreadOffsets() const;
    inline const DynamicList<dsmcParcel*>& moveAppendedParcels() const;
    inline const DynamicList<dsmcParcel*>& pendingMoveParcels() const;
    inline bool hasOccupancyOrderedParcels() const;
    inline const List<dsmcParcel*>& occupancyOrderedParcels() const;
    inline const labelList& occupancyCellOffsets() const;
    inline const labelList& occupancyActiveCells() const;
    inline const labelList& occupancyCollisionCells() const;
    inline label occupancyStart(const label celli) const;
    inline label occupancyEnd(const label celli) const;
    inline label occupancyCount(const label celli) const;
    inline dsmcParcel* occupancyParcel(const label index) const;
    inline const scalarField& selectedPairsPerCell() const;
    inline const labelList& nCandidatesPerCell() const;
    inline const labelList& collisionCandidateCells() const;
    inline const labelList& collisionLoadStart() const;
    inline const labelList& collisionLoadEnd() const;
    inline label randomLabel(const label min, const label max);
    inline const List<DynamicList<dsmcParcel*>>& cellOccupancy() const;
    inline volScalarField& sigmaTcRMax();
    inline scalarField& collisionSelectionRemainder();
    inline const List<dsmcParcel::constantProperties>& constProps() const;
    inline const dsmcParcel::constantProperties& constProps(label typeId) const;
    inline Random& rndGen();
    inline const BinaryCollisionModel& binaryCollision() const;
    inline BinaryCollisionModel& binaryCollision();
    inline const dsmcControllers& controllers() const;
    inline dsmcControllers& controllers();
    inline const boundaryMeasurements& boundaryFluxMeasurements() const;
    inline boundaryMeasurements& boundaryFluxMeasurements();
    inline const dsmcFieldProperties& fields() const;
    inline dsmcFieldProperties& fields();
    inline const dsmcBoundaries& boundaries() const;
    inline dsmcBoundaries& boundaries();
    inline const dsmcFaceTracker& tracker() const;
    inline dsmcFaceTracker& tracker();
    inline const cellMeasurements& cellPropMeasurements() const;
    inline cellMeasurements& cellPropMeasurements();

    inline bool reactionsActive() const
    {
        return reactions_.valid() && reactions_->nReactions() > 0;
    }

    inline const dsmcReactions& reactions() const
    {
        return reactions_();
    }

    inline dsmcReactions& reactions()
    {
        return reactions_();
    }

    inline scalar maxwellianMostProbableSpeed
    (
        const scalar temperature,
        const scalar mass
    ) const;

    vector equipartitionLinearVelocity
    (
        const scalar temperature,
        const scalar mass
    );

    scalar equipartitionRotationalEnergy
    (
        const scalar temperature,
        const scalar rotationalDof
    );

    labelList equipartitionVibrationalEnergyLevel
    (
        const scalar temperature,
        const label nVibrationalModes,
        const label typeId
    );

    label equipartitionElectronicLevel
    (
        const scalar temperature,
        const labelList& electronicDegeneracyList,
        const scalarList& electronicEnergyList
    );

    scalar postCollisionRotationalEnergy
    (
        const scalar rotationalDof,
        const scalar chiB
    );

    label postCollisionVibrationalEnergyLevel
    (
        bool postReaction,
        const label vibLevel,
        const label iMax,
        const scalar thetaV,
        const scalar thetaD,
        const scalar refTempZv,
        const scalar omega,
        const scalar Zref,
        const scalar Ec,
        const scalar fixedZv,
        const label invZvFormulation,
        const label cellI
    );

    label postCollisionElectronicEnergyLevel
    (
        const scalar Ec,
        const label jMax,
        const scalar omega,
        const scalarList& EElist,
        const labelList& gList
    );

    void handleWallInteraction(dsmcParcel& p, dsmcParcel::trackingData& td);
    void calculateAndWriteTimeStep(const Time& runTime, const dictionary& dsmcInitialiseDict);
    void initialiseFromDict(const dictionary& dsmcInitialiseDict);

    void info() const;
    void reportProfiling() const;
    void recordMoveThreadCounts(const labelList& counts);
    void recordMoveThreadProfile
    (
        const labelList& counts,
        const scalarField& wallTimes
    );
    void recordMoveInnerProfile
    (
        const scalarField& trackWallTimes,
        const scalarField& trackerWallTimes,
        const scalarField& boundaryWallTimes,
        const labelList& faceHitCounts,
        const labelList& cyclicHitCounts,
        const labelList& stuckHitCounts,
        const labelList& patchHitCounts,
        const labelList& processorHitCounts
    );
    void recordMovePhaseProfile
    (
        const scalar preControlWallTime,
        const scalar resetSetupWallTime,
        const scalar extractWallTime,
        const scalar kernelWallTime,
        const scalar commitWallTime,
        const scalar transferFinalizeWallTime
    );
    void recordCollisionThreadCounts(const labelList& candidateCounts, const labelList& acceptedCounts);
    void recordCollisionThreadProfile
    (
        const labelList& candidateCounts,
        const labelList& acceptedCounts,
        const labelList& activeCellCounts,
        const labelList& reactionHitCounts,
        const scalarField& wallTimes
    );
    void resetLoadStats();
    void loadBalanceCheck();
    void evolve();

    void addNewParcel
    (
        const barycentric& coordinates,
        const vector& U,
        const scalar RWF,
        const scalar ERot,
        const label ELevel,
        const label cellI,
        const label tetFaceI,
        const label tetPtI,
        const label typeId,
        const label newParcel,
        const label classification,
        const labelList& vibLevel
    );

    void addNewParcel
    (
        const vector& position,
        const vector& U,
        const scalar RWF,
        const scalar ERot,
        const label ELevel,
        const label cellI,
        const label tetFaceI,
        const label tetPtI,
        const label typeId,
        const label newParcel,
        const label classification,
        const labelList& vibLevel
    );

    void addNewStuckParcel
    (
        const vector& position,
        const vector& U,
        const scalar RWF,
        const scalar ERot,
        const label ELevel,
        const label cellI,
        const label tetFaceI,
        const label tetPtI,
        const label typeId,
        const label newParcel,
        const label classification,
        const labelList& vibLevel,
        const scalarField& wallTemperature,
        const vectorField& wallVectors
    );

    inline void reBuildCellOccupancy()
    {
        buildCellOccupancy();
    }

    inline void clear();
};

} // End namespace Foam

#include "dsmcCloudI.H"

#endif

// ************************************************************************* //










