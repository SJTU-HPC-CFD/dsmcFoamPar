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

\*---------------------------------------------------------------------------*/

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

inline const Foam::word& Foam::dsmcCloud::cloudName() const
{
    return cloudName_;
}

inline const Foam::fvMesh& Foam::dsmcCloud::mesh() const
{
    return mesh_;
}

inline const Foam::IOdictionary& Foam::dsmcCloud::particleProperties() const
{
    return particleProperties_;
}

inline const Foam::List<Foam::word>& Foam::dsmcCloud::typeIdList() const
{
    return typeIdList_;
}

inline const Foam::dsmcCoordinateSystem& Foam::dsmcCloud::coordSystem() const
{
    return dsmcCoordinateSystem_();
}

inline Foam::dsmcCoordinateSystem& Foam::dsmcCloud::coordSystem()
{
    return dsmcCoordinateSystem_();
}

inline const Foam::porousMeasurements& Foam::dsmcCloud::porousMeas() const
{
    return porousMeasurements_();
}

inline Foam::porousMeasurements& Foam::dsmcCloud::porousMeas()
{
    return porousMeasurements_();
}

inline Foam::scalar Foam::dsmcCloud::nParticle() const
{
    return dsmcCoordinateSystem_.valid() ? coordSystem().nParticlesOrg() : nParticle_;
}

inline Foam::scalar Foam::dsmcCloud::nParticles(const label celli) const
{
    return dsmcCoordinateSystem_.valid() ? coordSystem().nParticles(celli) : nParticle_;
}

inline Foam::scalar Foam::dsmcCloud::nParticles(const label patchi, const label facei) const
{
    return dsmcCoordinateSystem_.valid() ? coordSystem().nParticles(patchi, facei) : nParticle_;
}

inline Foam::scalar Foam::dsmcCloud::deltaTValue(const label celli) const
{
    return dsmcCoordinateSystem_.valid() ? coordSystem().deltaTValue(celli) : mesh_.time().deltaTValue();
}

inline Foam::label Foam::dsmcCloud::nTerminalOutputs() const
{
    return nTerminalOutputs_;
}

inline bool Foam::dsmcCloud::replicatedMeshActive() const
{
    return replicatedMesh_.valid();
}

inline const Foam::dsmcReplicatedMesh& Foam::dsmcCloud::replicatedMesh() const
{
    return replicatedMesh_();
}

inline bool Foam::dsmcCloud::openmpEnabled() const
{
    return openmpEnabled_;
}

inline bool Foam::dsmcCloud::openmpMoveEnabled() const
{
    return openmpMoveEnabled_;
}

inline bool Foam::dsmcCloud::trackerActive() const
{
    return trackerActive_;
}

inline bool Foam::dsmcCloud::evolveProfileEnabled() const
{
    return evolveProfileEnabled_;
}

inline bool Foam::dsmcCloud::profilingDetailEnabled() const
{
    return profilingDetailEnabled_;
}

inline Foam::scalar Foam::dsmcCloud::evolveMoveWallTime() const
{
    return evolveMoveWallTime_;
}

inline Foam::scalar Foam::dsmcCloud::evolveCollisionWallTime() const
{
    return evolveCollisionWallTime_;
}

inline Foam::label Foam::dsmcCloud::evolveProfileCalls() const
{
    return evolveProfileCalls_;
}

inline bool Foam::dsmcCloud::moveStageProbeEnabled() const
{
    return moveStageProbe_;
}

inline void Foam::dsmcCloud::setStepDiagnosticOutput(const bool enabled)
{
    emitStepDiagnostics_ = enabled;
}

inline bool Foam::dsmcCloud::emitStepDiagnostics() const
{
    return emitStepDiagnostics_;
}

inline Foam::label Foam::dsmcCloud::ompNumThreads() const
{
    return ompNumThreads_;
}

inline const Foam::word& Foam::dsmcCloud::openmpCollisionStrategy() const
{
    return openmpCollisionStrategy_;
}

inline Foam::label Foam::dsmcCloud::openmpCollisionChunk() const
{
    return openmpCollisionChunk_;
}

inline const Foam::word& Foam::dsmcCloud::openmpMoveSchedule() const
{
    return openmpMoveSchedule_;
}

inline Foam::label Foam::dsmcCloud::openmpMoveChunk() const
{
    return openmpMoveChunk_;
}

inline bool Foam::dsmcCloud::openmpMoveGuardCell(const label celli) const
{
    return
        celli >= 0
     && celli < openmpMoveGuardCells_.size()
     && openmpMoveGuardCells_[celli];
}

inline Foam::label Foam::dsmcCloud::currentThreadId() const
{
    #ifdef _OPENMP
    if (openmpEnabled_ && omp_in_parallel())
    {
        return omp_get_thread_num();
    }
    #endif

    return 0;
}

inline const Foam::labelList& Foam::dsmcCloud::particleLoadStart() const
{
    return particleLoadStart_;
}

inline const Foam::labelList& Foam::dsmcCloud::particleLoadEnd() const
{
    return particleLoadEnd_;
}

inline bool Foam::dsmcCloud::hasMoveOrderedParcels() const
{
    return moveOrderedParcelsValid_;
}

inline const Foam::List<Foam::dsmcParcel*>& Foam::dsmcCloud::moveOrderedParcels() const
{
    return moveOrderedParcels_;
}

inline const Foam::labelList& Foam::dsmcCloud::moveOrderedThreadOffsets() const
{
    return moveOrderedThreadOffsets_;
}

inline const Foam::DynamicList<Foam::dsmcParcel*>&
Foam::dsmcCloud::moveAppendedParcels() const
{
    return moveAppendedParcels_;
}

inline const Foam::DynamicList<Foam::dsmcParcel*>&
Foam::dsmcCloud::pendingMoveParcels() const
{
    return pendingMoveParcels_;
}

inline void Foam::dsmcCloud::recordPendingMigration(dsmcParcel* p)
{
    pendingMigration_.append(p);
}

inline const Foam::DynamicList<Foam::dsmcParcel*>&
Foam::dsmcCloud::pendingMigration() const
{
    return pendingMigration_;
}

inline void Foam::dsmcCloud::clearPendingMigration()
{
    pendingMigration_.clear();
}

inline bool Foam::dsmcCloud::hasOccupancyOrderedParcels() const
{
    return occupancyOrderedParcelsValid_;
}

inline const Foam::List<Foam::dsmcParcel*>& Foam::dsmcCloud::occupancyOrderedParcels() const
{
    return occupancyOrderedParcels_;
}

inline const Foam::labelList& Foam::dsmcCloud::occupancyCellOffsets() const
{
    return occupancyCellOffsets_;
}

inline const Foam::labelList& Foam::dsmcCloud::occupancyActiveCells() const
{
    return occupancyActiveCells_;
}

inline const Foam::labelList& Foam::dsmcCloud::occupancyCollisionCells() const
{
    return occupancyCollisionCells_;
}

inline Foam::label Foam::dsmcCloud::occupancyStart(const label celli) const
{
    if (occupancyOrderedParcelsValid_ || !cellOccupancyMaterialized_)
    {
        return occupancyCellOffsets_[celli];
    }

    return 0;
}

inline Foam::label Foam::dsmcCloud::occupancyEnd(const label celli) const
{
    if (occupancyOrderedParcelsValid_ || !cellOccupancyMaterialized_)
    {
        return occupancyCellOffsets_[celli + 1];
    }

    return 0;
}

inline Foam::label Foam::dsmcCloud::occupancyCount(const label celli) const
{
    if (occupancyOrderedParcelsValid_ || !cellOccupancyMaterialized_)
    {
        return occupancyCellOffsets_[celli + 1] - occupancyCellOffsets_[celli];
    }

    return cellOccupancy_[celli].size();
}

inline Foam::dsmcParcel* Foam::dsmcCloud::occupancyParcel(const label index) const
{
    return occupancyOrderedParcels_[index];
}

inline const Foam::scalarField& Foam::dsmcCloud::selectedPairsPerCell() const
{
    return selectedPairsPerCell_;
}

inline const Foam::labelList& Foam::dsmcCloud::nCandidatesPerCell() const
{
    return nCandidatesPerCell_;
}

inline const Foam::labelList& Foam::dsmcCloud::collisionCandidateCells() const
{
    return collisionCandidateCells_;
}

inline const Foam::labelList& Foam::dsmcCloud::collisionLoadStart() const
{
    return collisionLoadStart_;
}

inline const Foam::labelList& Foam::dsmcCloud::collisionLoadEnd() const
{
    return collisionLoadEnd_;
}

inline Foam::label Foam::dsmcCloud::randomLabel(const label min, const label max)
{
    return rndGen().position<label>(min, max);
}

inline const Foam::List<Foam::DynamicList<Foam::dsmcParcel*>>&
Foam::dsmcCloud::cellOccupancy() const
{
    if (!cellOccupancyMaterialized_)
    {
        const_cast<dsmcCloud&>(*this).materializeCellOccupancy();
    }

    return cellOccupancy_;
}

inline bool Foam::dsmcCloud::isCellOccupancyValid() const
{
    return cellOccupancyMaterialized_;
}

inline Foam::volScalarField& Foam::dsmcCloud::sigmaTcRMax()
{
    return sigmaTcRMax_;
}

inline Foam::scalarField& Foam::dsmcCloud::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}

inline const Foam::List<Foam::dsmcParcel::constantProperties>&
Foam::dsmcCloud::constProps() const
{
    return constProps_;
}

inline const Foam::dsmcParcel::constantProperties& Foam::dsmcCloud::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorInFunction
            << "constantProperties for requested typeId index " << typeId
            << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}

inline Foam::Random& Foam::dsmcCloud::rndGen()
{
    if (openmpEnabled_ && ompRndGens_.size())
    {
        return ompRndGens_[currentThreadId()];
    }

    return rndGen_;
}

inline const Foam::BinaryCollisionModel& Foam::dsmcCloud::binaryCollision() const
{
    return binaryCollisionModel_();
}

inline Foam::BinaryCollisionModel& Foam::dsmcCloud::binaryCollision()
{
    return binaryCollisionModel_();
}

inline const Foam::dsmcControllers& Foam::dsmcCloud::controllers() const
{
    return controllers_;
}

inline Foam::dsmcControllers& Foam::dsmcCloud::controllers()
{
    return controllers_;
}

inline const Foam::boundaryMeasurements& Foam::dsmcCloud::boundaryFluxMeasurements() const
{
    return boundaryMeas_;
}

inline Foam::boundaryMeasurements& Foam::dsmcCloud::boundaryFluxMeasurements()
{
    return boundaryMeas_;
}

inline const Foam::dsmcFieldProperties& Foam::dsmcCloud::fields() const
{
    return fields_;
}

inline Foam::dsmcFieldProperties& Foam::dsmcCloud::fields()
{
    return fields_;
}

inline const Foam::dsmcBoundaries& Foam::dsmcCloud::boundaries() const
{
    return boundaries_;
}

inline Foam::dsmcBoundaries& Foam::dsmcCloud::boundaries()
{
    return boundaries_;
}

inline const Foam::dsmcFaceTracker& Foam::dsmcCloud::tracker() const
{
    return trackingInfo_;
}

inline Foam::dsmcFaceTracker& Foam::dsmcCloud::tracker()
{
    return trackingInfo_;
}

inline const Foam::cellMeasurements& Foam::dsmcCloud::cellPropMeasurements() const
{
    return cellMeas_;
}

inline Foam::cellMeasurements& Foam::dsmcCloud::cellPropMeasurements()
{
    return cellMeas_;
}

inline Foam::scalar Foam::dsmcCloud::maxwellianMostProbableSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}

inline void Foam::dsmcCloud::clear()
{
    Cloud<dsmcParcel>::clear();
}

// ************************************************************************* //




