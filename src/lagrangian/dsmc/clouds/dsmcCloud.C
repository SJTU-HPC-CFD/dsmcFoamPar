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

#include "dsmcCloud.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"
#include <chrono>
#include <mpi.h>

using namespace Foam::constant;

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<dsmcParcel>, 0);
};


namespace
{
    typedef std::chrono::steady_clock steadyWallClock;

    Foam::scalar elapsedWallSeconds
    (
        const steadyWallClock::time_point& start
    )
    {
        return std::chrono::duration_cast
        <
            std::chrono::duration<Foam::scalar>
        >(steadyWallClock::now() - start).count();
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dsmcCloud::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< tab << id << endl;

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] = dsmcParcel::constantProperties(molDict);
    }
}


void Foam::dsmcCloud::buildCellOccupancy()
{
    const steadyWallClock::time_point wallStart = steadyWallClock::now();
    const scalar cpu0 = profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;
    const label nCells = mesh_.nCells();

    if (cellOccupancy_.size() != nCells)
    {
        cellOccupancy_.setSize(nCells);
    }

    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;

    const bool moveOrderedReady =
        openmpEnabled_
     && openmpMoveEnabled_
     && moveOrderedParcelsValid_
     && moveOrderedThreadOffsets_.size() == ompNumThreads_ + 1;

    const label appendedParcels = moveAppendedParcels_.size();
    const bool useMoveOrderedOnly =
        moveOrderedReady
     && moveOrderedParcels_.size() == this->size();
    const bool useMoveOrderedWithAppended =
        moveOrderedReady
     && !useMoveOrderedOnly
     && moveOrderedParcels_.size() + appendedParcels == this->size();
    const bool useMoveOrderedParcels =
        useMoveOrderedOnly || useMoveOrderedWithAppended;

    const label nParcels =
        useMoveOrderedParcels
      ? moveOrderedParcels_.size()
      + (useMoveOrderedWithAppended ? appendedParcels : 0)
      : this->size();

    #ifdef _OPENMP
    if (openmpEnabled_ && ompNumThreads_ > 1 && nParcels > 0)
    {
        const List<dsmcParcel*>* parcelsPtr = nullptr;
        const labelList* parcelThreadOffsetsPtr = nullptr;
        const DynamicList<dsmcParcel*>* appendedParcelsPtr = nullptr;
        labelList appendedThreadOffsets;
        List<dsmcParcel*> gatheredParcels;
        labelList generatedThreadOffsets;

        if (useMoveOrderedParcels)
        {
            parcelsPtr = &moveOrderedParcels_;
            parcelThreadOffsetsPtr = &moveOrderedThreadOffsets_;

            if (useMoveOrderedWithAppended)
            {
                appendedParcelsPtr = &moveAppendedParcels_;
                appendedThreadOffsets.setSize(ompNumThreads_ + 1, 0);

                for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
                {
                    appendedThreadOffsets[threadI] =
                        threadI*appendedParcels/ompNumThreads_;
                }

                appendedThreadOffsets[ompNumThreads_] = appendedParcels;
            }
        }
        else
        {
            gatheredParcels.setSize(nParcels);
            label parcelI = 0;

            forAllIter(dsmcCloud, *this, iter)
            {
                gatheredParcels[parcelI++] = &iter();
            }

            gatheredParcels.setSize(parcelI);
            generatedThreadOffsets.setSize(ompNumThreads_ + 1, 0);

            for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
            {
                generatedThreadOffsets[threadI] =
                    threadI*gatheredParcels.size()/ompNumThreads_;
            }

            generatedThreadOffsets[ompNumThreads_] = gatheredParcels.size();
            parcelsPtr = &gatheredParcels;
            parcelThreadOffsetsPtr = &generatedThreadOffsets;
        }

        const List<dsmcParcel*>& parcels = *parcelsPtr;
        const labelList& parcelThreadOffsets = *parcelThreadOffsetsPtr;
        const bool useAppendedParcels = appendedParcelsPtr != nullptr;

        if (occupancyThreadCellCounts_.size() != ompNumThreads_)
        {
            occupancyThreadCellCounts_.setSize(ompNumThreads_);
        }

        if (occupancyThreadActiveCells_.size() != ompNumThreads_)
        {
            occupancyThreadActiveCells_.setSize(ompNumThreads_);
        }

        forAll(occupancyThreadCellCounts_, threadI)
        {
            labelList& localCounts = occupancyThreadCellCounts_[threadI];
            DynamicList<label>& localActiveCells =
                occupancyThreadActiveCells_[threadI];

            if (localCounts.size() != nCells)
            {
                localCounts.setSize(nCells, 0);
            }
            else if (!useMoveOrderedParcels)
            {
                forAll(localCounts, cellI)
                {
                    localCounts[cellI] = 0;
                }
            }
            else
            {
                forAll(localActiveCells, activeI)
                {
                    localCounts[localActiveCells[activeI]] = 0;
                }
            }

            localActiveCells.clear();
            localActiveCells.setCapacity
            (
                min
                (
                    nCells,
                    (parcelThreadOffsets[threadI + 1]
                   - parcelThreadOffsets[threadI])
                  + (useAppendedParcels
                    ? appendedThreadOffsets[threadI + 1]
                    - appendedThreadOffsets[threadI]
                    : 0)
                )
            );
        }

        #pragma omp parallel num_threads(ompNumThreads_)
        {
            const label threadI = omp_get_thread_num();
            labelList& localCounts = occupancyThreadCellCounts_[threadI];
            DynamicList<label>& localActiveCells =
                occupancyThreadActiveCells_[threadI];

            for
            (
                label i = parcelThreadOffsets[threadI];
                i < parcelThreadOffsets[threadI + 1];
                ++i
            )
            {
                const label cellI = parcels[i]->cell();

                if (cellI >= 0 && cellI < nCells)
                {
                    label& count = localCounts[cellI];

                    if (count == 0)
                    {
                        localActiveCells.append(cellI);
                    }

                    ++count;
                }
            }

            if (useAppendedParcels)
            {
                const DynamicList<dsmcParcel*>& appended = *appendedParcelsPtr;

                for
                (
                    label i = appendedThreadOffsets[threadI];
                    i < appendedThreadOffsets[threadI + 1];
                    ++i
                )
                {
                    const label cellI = appended[i]->cell();

                    if (cellI >= 0 && cellI < nCells)
                    {
                        label& count = localCounts[cellI];

                        if (count == 0)
                        {
                            localActiveCells.append(cellI);
                        }

                        ++count;
                    }
                }
            }
        }

        labelList totalCounts(nCells, 0);

        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            const labelList& localCounts = occupancyThreadCellCounts_[threadI];
            const DynamicList<label>& activeCells =
                occupancyThreadActiveCells_[threadI];

            forAll(activeCells, activeI)
            {
                const label cellI = activeCells[activeI];
                totalCounts[cellI] += localCounts[cellI];
            }
        }

        label activeCellCount = 0;
        label collisionCellCount = 0;

        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            const label count = totalCounts[cellI];

            if (count > 0)
            {
                ++activeCellCount;

                if (count > 1)
                {
                    ++collisionCellCount;
                }
            }
        }

        occupancyActiveCells_.setSize(activeCellCount);
        occupancyCollisionCells_.setSize(collisionCellCount);

        activeCellCount = 0;
        collisionCellCount = 0;

        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            const label count = totalCounts[cellI];

            if (count > 0)
            {
                occupancyActiveCells_[activeCellCount++] = cellI;

                if (count > 1)
                {
                    occupancyCollisionCells_[collisionCellCount++] = cellI;
                }
            }
        }

        occupancyCellOffsets_.setSize(nCells + 1, 0);

        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            occupancyCellOffsets_[cellI + 1] =
                occupancyCellOffsets_[cellI] + totalCounts[cellI];
        }

        occupancyOrderedParcels_.resize(occupancyCellOffsets_.last());
        labelList nextCellOffsets(occupancyCellOffsets_);

        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            labelList& localCounts = occupancyThreadCellCounts_[threadI];
            const DynamicList<label>& activeCells =
                occupancyThreadActiveCells_[threadI];

            forAll(activeCells, activeI)
            {
                const label cellI = activeCells[activeI];
                const label count = localCounts[cellI];
                localCounts[cellI] = nextCellOffsets[cellI];
                nextCellOffsets[cellI] += count;
            }
        }

        #pragma omp parallel num_threads(ompNumThreads_)
        {
            const label threadI = omp_get_thread_num();
            labelList& localOffsets = occupancyThreadCellCounts_[threadI];

            for
            (
                label i = parcelThreadOffsets[threadI];
                i < parcelThreadOffsets[threadI + 1];
                ++i
            )
            {
                dsmcParcel* pPtr = parcels[i];
                const label cellI = pPtr->cell();

                if (cellI >= 0 && cellI < nCells)
                {
                    const label slot = localOffsets[cellI]++;
                    occupancyOrderedParcels_[slot] = pPtr;
                }
            }

            if (useAppendedParcels)
            {
                const DynamicList<dsmcParcel*>& appended = *appendedParcelsPtr;

                for
                (
                    label i = appendedThreadOffsets[threadI];
                    i < appendedThreadOffsets[threadI + 1];
                    ++i
                )
                {
                    dsmcParcel* pPtr = appended[i];
                    const label cellI = pPtr->cell();

                    if (cellI >= 0 && cellI < nCells)
                    {
                        const label slot = localOffsets[cellI]++;
                        occupancyOrderedParcels_[slot] = pPtr;
                    }
                }
            }
        }

        occupancyOrderedParcelsValid_ = true;
        cellOccupancyMaterialized_ = false;

        if (profileSummary_ && profileTimingActive_)
        {
            profileBuildCellOccupancyWall_ += elapsedWallSeconds(wallStart);
            profileBuildCellOccupancyCpu_ +=
                mesh_.time().elapsedCpuTime() - cpu0;
        }

        return;
    }
    #endif

    forAll(cellOccupancy_, celli)
    {
        cellOccupancy_[celli].clear();
    }

    forAllIter(dsmcCloud, *this, iter)
    {
        const label cellI = iter().cell();

        if (cellI >= 0 && cellI < nCells)
        {
            cellOccupancy_[cellI].append(&iter());
        }
    }

    labelList totalCounts(nCells, 0);
    occupancyCellOffsets_.setSize(nCells + 1, 0);

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        totalCounts[cellI] = cellOccupancy_[cellI].size();
        occupancyCellOffsets_[cellI + 1] =
            occupancyCellOffsets_[cellI] + totalCounts[cellI];
    }

    label activeCellCount = 0;
    label collisionCellCount = 0;

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        const label count = totalCounts[cellI];

        if (count > 0)
        {
            ++activeCellCount;

            if (count > 1)
            {
                ++collisionCellCount;
            }
        }
    }

    occupancyActiveCells_.setSize(activeCellCount);
    occupancyCollisionCells_.setSize(collisionCellCount);
    activeCellCount = 0;
    collisionCellCount = 0;

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        const label count = totalCounts[cellI];

        if (count > 0)
        {
            occupancyActiveCells_[activeCellCount++] = cellI;

            if (count > 1)
            {
                occupancyCollisionCells_[collisionCellCount++] = cellI;
            }
        }
    }

    occupancyOrderedParcels_.resize(occupancyCellOffsets_.last());

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        label offset = occupancyCellOffsets_[cellI];
        const DynamicList<dsmcParcel*>& cellParcels = cellOccupancy_[cellI];

        forAll(cellParcels, i)
        {
            occupancyOrderedParcels_[offset++] = cellParcels[i];
        }
    }

    occupancyOrderedParcelsValid_ = true;
    cellOccupancyMaterialized_ = true;

    if (profileSummary_ && profileTimingActive_)
    {
        profileBuildCellOccupancyWall_ += elapsedWallSeconds(wallStart);
        profileBuildCellOccupancyCpu_ += mesh_.time().elapsedCpuTime() - cpu0;
    }
}


void Foam::dsmcCloud::materializeCellOccupancy()
{
    if (cellOccupancyMaterialized_)
    {
        return;
    }

    const label nCells = mesh_.nCells();

    if (cellOccupancy_.size() != nCells)
    {
        cellOccupancy_.setSize(nCells);
    }

    forAll(cellOccupancy_, cellI)
    {
        cellOccupancy_[cellI].clear();
    }

    if (occupancyOrderedParcelsValid_)
    {
        forAll(occupancyActiveCells_, activeI)
        {
            const label cellI = occupancyActiveCells_[activeI];
            DynamicList<dsmcParcel*>& cellParcels = cellOccupancy_[cellI];

            for
            (
                label occI = occupancyCellOffsets_[cellI];
                occI < occupancyCellOffsets_[cellI + 1];
                ++occI
            )
            {
                cellParcels.append(occupancyOrderedParcels_[occI]);
            }
        }
    }

    cellOccupancyMaterialized_ = true;
}


void Foam::dsmcCloud::relocateStuckParcels()
{
    forAllIter(dsmcCloud, *this, iter)
    {
        if (iter().isStuck())
        {
            bool success = iter().relocateStuckParcel(mesh_);
            if (!success)
            {
                FatalErrorInFunction
                    << "Could not relocate stuck parcel!"
                    << exit(FatalError);
            }
        }
    }
}


void Foam::dsmcCloud::buildCellOccupancyFromScratch()
{
    cellOccupancy_.clear();
    cellOccupancy_.setSize(mesh_.nCells());

    nCandidatesPerCell_.setSize(mesh_.nCells(), 0);
    moveItersPerCell_.setSize(mesh_.nCells(), 0);
    moveItersPerCellCumulative_.setSize(mesh_.nCells(), 0);

    buildCellOccupancy();
    relocateStuckParcels();
}


void Foam::dsmcCloud::buildCollisionSelectionRemainderFromScratch()
{
    collisionSelectionRemainder_.clear();
    collisionSelectionRemainder_.setSize(mesh_.nCells());

    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, cO)
    {
        collisionSelectionRemainder_[cO] = rndGen_.sample01<scalar>();
    }
}


void Foam::dsmcCloud::resetBoundaries()
{
    boundaryMeas_.reset();
    boundaries_.setNewConfig();
}


void Foam::dsmcCloud::resetMeasurementTools()
{
    refreshTrackerUsage();
    trackingInfo_.reset();
    fields_.resetFields();

    cellMeas_.reset();
}


void Foam::dsmcCloud::refreshTrackerUsage()
{
    trackerActive_ = false;

    const List< autoPtr<dsmcField> >& configuredFields = fields_.fields();

    forAll(configuredFields, fieldI)
    {
        if (configuredFields[fieldI].valid())
        {
            const word fieldType(configuredFields[fieldI]->type());

            if (fieldType == "dsmcFluxSurface")
            {
                trackerActive_ = true;
                break;
            }
        }
    }
}


void Foam::dsmcCloud::removeElectrons()
{
    forAll(cellOccupancy_, c)
    {
        const DynamicList<dsmcParcel*>& molsInCell = cellOccupancy_[c];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];

            const dsmcParcel::constantProperties& constProp =
                constProps(p->typeId());

            const label& charge = constProp.charge();

            const scalar& RWF = coordSystem().RWF(c);

            const scalar mass = constProps(p->typeId()).mass();

            momentumMean_[c] += mass*RWF*p->U();
            rhoMMean_[c] += mass*RWF;

            if (charge == -1)
            {
                rhoNMeanElectron_[c] += 1.0*RWF;
                rhoMMeanElectron_[c] += mass*RWF;
                momentumMeanElectron_[c] += mass*RWF*p->U();
                linearKEMeanElectron_[c] += mass*RWF*(p->U() & p->U());

                //- found an electron
                deleteParticle(*p);
            }
        }
    }
}


void Foam::dsmcCloud::addElectrons()
{
    label electronTypeId = -1;

    //- find electron typeId
    forAll(constProps_, cP)
    {
        const label& particleCharge = constProps_[cP].charge();

        if (particleCharge == -1)
        {
            electronTypeId = cP;
            break;
        }
    }

    forAll(cellOccupancy_, c)
    {
        if (rhoMMeanElectron_[c] > VSMALL)
        {
            scalar V = mesh_.cellVolumes()[c];

            scalar rhoMMeanElectron = rhoMMeanElectron_[c]*nParticles(c)/V;
            scalar rhoNMeanElectron = rhoNMeanElectron_[c]*nParticles(c)/V;
            vector UElectron = momentumMeanElectron_[c] /(rhoMMeanElectron*V);
            scalar linearKEMeanElectron =
		            (0.5*linearKEMeanElectron_[c]*nParticles(c))/V;

            electronTemperature_[c] = 2.0/(3.0*physicoChemical::k.value()
                * rhoNMeanElectron)*(linearKEMeanElectron
                - 0.5*rhoMMeanElectron
                * (UElectron & UElectron));
        }

        const DynamicList<dsmcParcel*>& molsInCell = cellOccupancy_[c];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];

            const dsmcParcel::constantProperties& constProp =
                constProps(p->typeId());

            label charge = constProp.charge();

            if (charge == 1)
            {
                const label& cellI = p->cell();

                //- found an ion, add an electron here

                //- electron temperature will be zero if there have been no
                //  electrons in the cell during the simulation

                if (electronTemperature_[cellI] < SMALL)
                {
                    electronTemperature_[cellI] = 6000.0;
                }
                if (electronTemperature_[cellI] > 8.0e4)
                {
                    electronTemperature_[cellI] = 30000.0;
                }

                vector electronVelocity = equipartitionLinearVelocity
                    (
                        electronTemperature_[cellI],
                        constProps_[electronTypeId].mass()
                    );

                if (rhoMMean_[cellI] > VSMALL)
                {
                    cellVelocity_[cellI] = momentumMean_[cellI]
                        /rhoMMean_[cellI];
                }

                labelList vibLevel;

                electronVelocity += cellVelocity_[cellI];

                addNewParcel
                (
                    p->position(),
                    electronVelocity,
                    p->RWF(),
                    0.0,
                    0,
                    p->cell(),
                    p->tetFace(),
                    p->tetPt(),
                    electronTypeId,
                    -1,
                    0,
                    vibLevel
                );
            }
        }
    }
}


Foam::label Foam::dsmcCloud::pickFromCandidateList
(
    DynamicList<label>& candidatesInCell
)
{
    label entry = -1;
    label size = candidatesInCell.size();

    if (size > 0)
    {
        // choose a random number between 0 and the size of the candidateList
        //label randomIndex = rndGen_.position<label>(0, size - 1); OLD
        label randomIndex = randomLabel(0, size-1);
        entry = candidatesInCell[randomIndex];

        // build a new list without the chosen entry
        DynamicList<label> newCandidates(0);

        forAll(candidatesInCell, i)
        {
            if (i != randomIndex)
            {
                newCandidates.append(candidatesInCell[i]);
            }
        }

        // transfer the new list
        candidatesInCell.transfer(newCandidates);
        candidatesInCell.shrink();
    }

    return entry;
}

void Foam::dsmcCloud::updateCandidateSubList
(
    const label& candidate,
    DynamicList<label>& candidatesInSubCell
)
{
    label newIndex = findIndex(candidatesInSubCell, candidate);

    DynamicList<label> newCandidates(0);

    forAll(candidatesInSubCell, i)
    {
        if (i != newIndex)
        {
            newCandidates.append(candidatesInSubCell[i]);
        }
    }

    // transfer the new list
    candidatesInSubCell.transfer(newCandidates);
    candidatesInSubCell.shrink();

//     Info <<  " list (after) " << candidatesInSubCell << endl;
}


Foam::label Foam::dsmcCloud::pickFromCandidateSubList
(
    DynamicList<label>& candidatesInCell,
    DynamicList<label>& candidatesInSubCell
)
{
//     Info << " list (before) " << candidatesInCell << endl;
//     Info << " sub list (before) " << candidatesInSubCell << endl;


    label entry = -1;
    label subCellSize = candidatesInSubCell.size();

    if (subCellSize > 0)
    {
        //label randomIndex = rndGen_.position<label>(0, subCellSize - 1); OLD
        label randomIndex = randomLabel(0, subCellSize-1);
        entry = candidatesInSubCell[randomIndex];

//         Info<< "random index: " << randomIndex <<" entry "
//             << entry << endl;

        DynamicList<label> newSubCellList(0);

        forAll(candidatesInSubCell, i)
        {
            if (i != randomIndex)
            {
                newSubCellList.append(candidatesInSubCell[i]);
            }
        }

        candidatesInSubCell.transfer(newSubCellList);
        candidatesInSubCell.shrink();

//         Info <<  " sub list (after) " << candidatesInSubCell << endl;

        label newIndex = findIndex(candidatesInCell, entry);

        DynamicList<label> newList(0);

        forAll(candidatesInCell, i)
        {
            if (i != newIndex)
            {
                newList.append(candidatesInCell[i]);
            }
        }

        candidatesInCell.transfer(newList);
        candidatesInCell.shrink();

//         Info <<  " list (after) " << candidatesInCell << endl;
    }

    return entry;
}

void Foam::dsmcCloud::collisions()
{
    collisionPartnerSelectionModel_->collide();
}


void Foam::dsmcCloud::addNewParcel
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
)
{
    if
    (
        replicatedMeshActive()
     && cellI >= 0
     && cellI < replicatedMesh().cellOwner().size()
     && replicatedMesh().cellOwner()[cellI] != replicatedMesh().myRank()
    )
    {
        return;
    }

    dsmcParcel* pPtr = new dsmcParcel
    (
        mesh_,
        position,
        U,
        RWF,
        ERot,
        ELevel,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        newParcel,
        classification,
        vibLevel
    );

    // parcels that have been freshly injected on boundary patch faces should
    // be tracked as having crossed that boundary patch face in the time step
    // in which they have been inserted. This is justified as the trackFraction
    // is initialized to a random value in the interval [0, 1] (cf.
    // dsmcParcel::move newParcel handling).
    if (newParcel != -1 && trackerActive())
    {
        tracker().trackFaceTransition(typeId, U, RWF, tetFaceI);
    }

    porousMeas().additionInteraction(*pPtr, newParcel);

    addParticle(pPtr);
    recordMoveAppendedParcel(pPtr);
}


void Foam::dsmcCloud::addNewStuckParcel
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
)
{
    if
    (
        replicatedMeshActive()
     && cellI >= 0
     && cellI < replicatedMesh().cellOwner().size()
     && replicatedMesh().cellOwner()[cellI] != replicatedMesh().myRank()
    )
    {
        return;
    }

    dsmcParcel* pPtr = new dsmcParcel
    (
        mesh_,
        position,
        U,
        RWF,
        ERot,
        ELevel,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        newParcel,
        classification,
        vibLevel
    );

    // set parcel stuck
    dsmcParcel& p = *pPtr;
    p.setStuck();
    p.stuck().wallTemperature() = wallTemperature;
    p.stuck().wallVectors() = wallVectors;

    porousMeas().additionInteraction(p, newParcel);

    addParticle(pPtr);
    recordMoveAppendedParcel(pPtr);
}


Foam::scalar Foam::dsmcCloud::energyRatio
(
    scalar ChiA,
    scalar ChiB
)
{
    scalar ChiAMinusOne = ChiA - 1;

    scalar ChiBMinusOne = ChiB - 1;

    if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
    {
        return rndGen_.sample01<scalar>();
    }

    scalar energyRatio;

    scalar P;

    do
    {
        P = 0;

        energyRatio = rndGen_.sample01<scalar>();

        if (ChiAMinusOne < SMALL)
        {
            P = pow((1.0 - energyRatio),ChiBMinusOne);
        }
        else if (ChiBMinusOne < SMALL)
        {
            P = pow((1.0 - energyRatio),ChiAMinusOne);
        }
        else
        {
            P =
                pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne,
                    ChiAMinusOne
                )
               *pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*(1 - energyRatio)
                    /ChiBMinusOne,
                    ChiBMinusOne
                );
        }
    } while (P < rndGen_.sample01<scalar>());

    return energyRatio;
}

Foam::scalar Foam::dsmcCloud::PSIm
(
    scalar DOFm,
    scalar DOFtot
)
{
    if (DOFm == DOFtot)
    {
        return 1.0;
    }

    if (DOFm == 2.0 && DOFtot == 4.0)
    {
        return rndGen_.sample01<scalar>();
    }

    if (DOFtot < 4.0)
    {
        return (DOFm/DOFtot);
    }

    scalar rPSIm = 0.0;
    scalar prob = 0.0;

    scalar h1 = 0.5*DOFtot - 2.0;
    scalar h2 = 0.5*DOFm - 1.0 + 1.0e-5;
    scalar h3 = 0.5*(DOFtot-DOFm)-1.0 + 1.0e-5;

    do
    {
        rPSIm = rndGen_.sample01<scalar>();
        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);
    } while (prob < rndGen_.sample01<scalar>());

    return rPSIm;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// for running dsmcFoam+
Foam::dsmcCloud::dsmcCloud
(
    Time& t,
    const word& cloudName,
    const dynamicFvMesh& mesh,
    bool readFields
)
:
    Cloud<dsmcParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    controlDict_(t.controlDict()),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    dsmcCoordinateSystem_(dsmcCoordinateSystem::New(t, mesh, *this)),
    uniformDeltaT_(dsmcCoordinateSystem_->dtModel().uniformDeltaT()),
    porousMeasurements_(porousMeasurements::New(t, mesh, *this)),
    nTerminalOutputs_
    (
        controlDict_.lookupOrDefault<label>("nTerminalOutputs", 1)
    ),
    profileSummary_
    (
        controlDict_.lookupOrDefault<bool>("profileSummary", true)
    ),
    profileDetail_
    (
        controlDict_.lookupOrDefault<bool>("profileDetail", false)
    ),
    profileTimingActive_(false),
    profileSteps_(0),
    moveDetailParcels_(0),
    moveDetailTrackCalls_(0),
    moveDetailSameTetNoFaceHits_(0),
    moveDetailInternalTetNoFaceHits_(0),
    moveDetailFaceHits_(0),
    moveDetailProcessorHits_(0),
    moveDetailCyclicHits_(0),
    moveDetailPatchHits_(0),
    moveDetailStuckHits_(0),
    moveDetailTrackWallTime_(0.0),
    moveDetailTrackerWallTime_(0.0),
    moveDetailBoundaryWallTime_(0.0),
    openmpEnabled_
    (
        controlDict_.lookupOrDefault<bool>("useOpenMP", false)
    ),
    openmpMoveEnabled_
    (
        controlDict_.lookupOrDefault<bool>("openmpMove", openmpEnabled_)
    ),
    trackerActive_(true),
    ompNumThreads_(1),
    openmpCollisionSchedule_
    (
        controlDict_.lookupOrDefault<word>("openmpCollisionSchedule", "dynamic")
    ),
    openmpCollisionChunk_
    (
        max(label(1), controlDict_.lookupOrDefault<label>("openmpCollisionChunk", 1))
    ),
    openmpMoveSchedule_
    (
        controlDict_.lookupOrDefault<word>("openmpMoveSchedule", "static")
    ),
    openmpMoveChunk_
    (
        max(label(1), controlDict_.lookupOrDefault<label>("openmpMoveChunk", 1))
    ),
    nCandidatesPerCell_(mesh_.nCells(), 0),
    moveItersPerCell_(mesh_.nCells(), 0),
    moveItersPerCellCumulative_(mesh_.nCells(), 0),
    moveOrderedParcels_(),
    moveOrderedThreadOffsets_(),
    moveOrderedParcelsValid_(false),
    moveAppendCaptureActive_(false),
    moveAppendedParcels_(),
    profileFullEvolveWall_(0.0),
    profileMoveAndCollideWall_(0.0),
    profileMoveWall_(0.0),
    profileBuildCellOccupancyWall_(0.0),
    profileCollisionWall_(0.0),
    profilePostFieldsWall_(0.0),
    profileFullEvolveCpu_(0.0),
    profileMoveAndCollideCpu_(0.0),
    profileMoveCpu_(0.0),
    profileBuildCellOccupancyCpu_(0.0),
    profileCollisionCpu_(0.0),
    profilePostFieldsCpu_(0.0),
    occupancyOrderedParcels_(),
    occupancyCellOffsets_(),
    occupancyActiveCells_(),
    occupancyCollisionCells_(),
    occupancyOrderedParcelsValid_(false),
    cellOccupancyMaterialized_(true),
    occupancyThreadCellCounts_(),
    occupancyThreadActiveCells_(),
    cellOccupancy_(),
    rhoNMeanElectron_(mesh_.nCells(), 0.0),
    rhoMMeanElectron_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    momentumMeanElectron_(mesh_.nCells(), vector::zero),
    momentumMean_(mesh_.nCells(), vector::zero),
    linearKEMeanElectron_(mesh_.nCells(), 0.0),
    electronTemperature_(mesh_.nCells(), 0.0),
    cellVelocity_(mesh_.nCells(), vector::zero),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    collisionSelectionRemainder_(),
    constProps_(),
    rndGen_
    (
        particleProperties_.lookupOrDefault<label>
        (
            "seedNumber",
            label(clock::getTime()) + 7183*Pstream::myProcNo()
        )
    ), 
    controllers_(t, mesh, *this),
    dynamicLoadBalancing_(t, mesh, *this),
    replicatedMesh_(),
    boundaryMeas_(mesh, *this, true),
    fields_(t, mesh, *this),
    boundaries_(t, mesh, *this),
    trackingInfo_(mesh, *this, true),
    binaryCollisionModel_
    (
        BinaryCollisionModel::New
        (
            particleProperties_,
            *this
        )
    ),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh, *this),
    cellMeas_(mesh, *this, true)
{
    if
    (
        openmpCollisionSchedule_ != "static"
     && openmpCollisionSchedule_ != "dynamic"
     && openmpCollisionSchedule_ != "guided"
     && openmpCollisionSchedule_ != "partition"
    )
    {
        WarningInFunction
            << "Unknown openmpCollisionSchedule '"
            << openmpCollisionSchedule_
            << "', falling back to dynamic" << endl;
        openmpCollisionSchedule_ = "dynamic";
    }

    if
    (
        openmpMoveSchedule_ != "static"
     && openmpMoveSchedule_ != "dynamic"
     && openmpMoveSchedule_ != "guided"
    )
    {
        WarningInFunction
            << "Unknown openmpMoveSchedule '"
            << openmpMoveSchedule_
            << "', falling back to static" << endl;
        openmpMoveSchedule_ = "static";
    }

    #ifdef _OPENMP
    ompNumThreads_ =
        max
        (
            label(1),
            controlDict_.lookupOrDefault<label>
            (
                "openmpThreads",
                label(omp_get_max_threads())
            )
        );
    omp_set_num_threads(int(ompNumThreads_));
    #else
    openmpEnabled_ = false;
    openmpMoveEnabled_ = false;
    ompNumThreads_ = 1;
    #endif

    if (readFields)
    {
        dsmcParcel::readFields(*this);
    }

    buildConstProps();

    coordSystem().checkCoordinateSystemInputs();
    porousMeas().checkPorousMeasurementsInputs();

    reactions_.initialConfiguration();

    buildCellOccupancyFromScratch();
    buildCollisionSelectionRemainderFromScratch();

    collisionPartnerSelectionModel_ = autoPtr<collisionPartnerSelection>
    (
        collisionPartnerSelection::New(mesh, *this, particleProperties_)
    );

    collisionPartnerSelectionModel_->initialConfiguration();

    fields_.createFields();
    refreshTrackerUsage();
    boundaryMeas_.setInitialConfig();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();

    if (controlDict_.lookupOrDefault<bool>("replicatedMesh", false))
    {
        replicatedMesh_.reset(new dsmcReplicatedMesh(*this, mesh_));
        replicatedMesh_->initialize();
    }
}


// running dsmcInitialise+
Foam::dsmcCloud::dsmcCloud
(
    Time& t,
    const word& cloudName,
    const dynamicFvMesh& mesh,
    const IOdictionary& dsmcInitialiseDict,
    const bool& clearFields
)
:
    Cloud<dsmcParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    controlDict_(t.controlDict()),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    dsmcCoordinateSystem_(dsmcCoordinateSystem::New(t, mesh, *this)),
    uniformDeltaT_(dsmcCoordinateSystem_->dtModel().uniformDeltaT()),
    porousMeasurements_(porousMeasurements::New(t, mesh, *this)),
    nTerminalOutputs_
    (
        controlDict_.lookupOrDefault<label>("nTerminalOutputs", 1)
    ),
    profileSummary_
    (
        controlDict_.lookupOrDefault<bool>("profileSummary", true)
    ),
    profileDetail_
    (
        controlDict_.lookupOrDefault<bool>("profileDetail", false)
    ),
    profileTimingActive_(false),
    profileSteps_(0),
    moveDetailParcels_(0),
    moveDetailTrackCalls_(0),
    moveDetailSameTetNoFaceHits_(0),
    moveDetailInternalTetNoFaceHits_(0),
    moveDetailFaceHits_(0),
    moveDetailProcessorHits_(0),
    moveDetailCyclicHits_(0),
    moveDetailPatchHits_(0),
    moveDetailStuckHits_(0),
    moveDetailTrackWallTime_(0.0),
    moveDetailTrackerWallTime_(0.0),
    moveDetailBoundaryWallTime_(0.0),
    openmpEnabled_(false),
    openmpMoveEnabled_(false),
    trackerActive_(true),
    ompNumThreads_(1),
    openmpCollisionSchedule_("dynamic"),
    openmpCollisionChunk_(1),
    openmpMoveSchedule_("static"),
    openmpMoveChunk_(1),
    nCandidatesPerCell_(),
    moveItersPerCell_(),
    moveItersPerCellCumulative_(),
    moveOrderedParcels_(),
    moveOrderedThreadOffsets_(),
    moveOrderedParcelsValid_(false),
    moveAppendCaptureActive_(false),
    moveAppendedParcels_(),
    profileFullEvolveWall_(0.0),
    profileMoveAndCollideWall_(0.0),
    profileMoveWall_(0.0),
    profileBuildCellOccupancyWall_(0.0),
    profileCollisionWall_(0.0),
    profilePostFieldsWall_(0.0),
    profileFullEvolveCpu_(0.0),
    profileMoveAndCollideCpu_(0.0),
    profileMoveCpu_(0.0),
    profileBuildCellOccupancyCpu_(0.0),
    profileCollisionCpu_(0.0),
    profilePostFieldsCpu_(0.0),
    occupancyOrderedParcels_(),
    occupancyCellOffsets_(),
    occupancyActiveCells_(),
    occupancyCollisionCells_(),
    occupancyOrderedParcelsValid_(false),
    cellOccupancyMaterialized_(true),
    occupancyThreadCellCounts_(),
    occupancyThreadActiveCells_(),
    cellOccupancy_(),
    rhoNMeanElectron_(),
    rhoMMeanElectron_(),
    rhoMMean_(),
    momentumMeanElectron_(),
    momentumMean_(),
    linearKEMeanElectron_(),
    electronTemperature_(),
    cellVelocity_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 3, -1, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    collisionSelectionRemainder_(),
    constProps_(),
    rndGen_
    (
        particleProperties_.lookupOrDefault<label>
        (
            "seedNumber",
            label(clock::getTime()) + 1526*Pstream::myProcNo()
        )
    ),
    controllers_(t, mesh),
    dynamicLoadBalancing_(t, mesh, *this),
    replicatedMesh_(),
    boundaryMeas_(mesh, *this),
    fields_(t, mesh),
    boundaries_(t, mesh),
    trackingInfo_(mesh, *this),
    binaryCollisionModel_(),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh),
    cellMeas_(mesh, *this)
{
    if (!clearFields)
    {
        dsmcParcel::readFields(*this);
    }

    label initialParcels = this->size();

    if (Pstream::parRun())
    {
        reduce(initialParcels, sumOp<label>());
    }

    if (clearFields)
    {
        Info << "clearing existing field of parcels " << endl;

        clear();

        initialParcels = 0;
    }

    buildConstProps();

    coordSystem().checkCoordinateSystemInputs(true);

    dsmcAllConfigurations conf(dsmcInitialiseDict, *this);
    conf.setInitialConfig();

    label finalParcels = this->size();

    if (Pstream::parRun())
    {
        reduce(finalParcels, sumOp<label>());
    }

    Info << nl << "Initial no. of parcels: " << initialParcels
         << " added parcels: " << finalParcels - initialParcels
         << ", total no. of parcels: " << finalParcels
         << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dsmcCloud::~dsmcCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcCloud::evolve()
{
    const steadyWallClock::time_point wallStart = steadyWallClock::now();
    const scalar cpu0 = profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;

    profileTimingActive_ = profileSummary_;

    evolve_moveAndCollide();
    evolve_fields();

    profileTimingActive_ = false;

    if (profileSummary_)
    {
        profileFullEvolveWall_ += elapsedWallSeconds(wallStart);
        profileFullEvolveCpu_ += mesh_.time().elapsedCpuTime() - cpu0;
        ++profileSteps_;
    }
}


void Foam::dsmcCloud::evolve_moveAndCollide()
{
    const steadyWallClock::time_point moveAndCollideWallStart =
        steadyWallClock::now();
    const scalar moveAndCollideCpuStart =
        profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;

    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    dsmcParcel::trackingData td(*this);
    td.moveDetailProfile =
        profileDetail_
     && controlDict_.lookupOrDefault<bool>("moveDetailProfile", false);

    if (debug)
    {
        this->dumpParticlePositions();
    }

    if (openmpEnabled_ && openmpMoveEnabled_)
    {
        beginMoveAppendCapture();
    }

    controllers_.controlBeforeMove();
    boundaries_.controlBeforeMove();

    if (replicatedMeshActive() && replicatedMesh_->migrationCalls() == 0)
    {
        Info<< "Replicated mesh: initial particle distribution" << nl << endl;
        if (Pstream::parRun())
        {
            replicatedMesh_->distributeInitialParticles();
        }
        else
        {
            replicatedMesh_->migrateParticlesByCellOwner();
        }
        replicatedMesh_->updateParticleCounts();
        clearMoveOrderedParcels();
    }

    const bool replicatedMeshDelayedReceive =
        replicatedMeshActive()
     && controlDict_.lookupOrDefault<bool>
        (
            "replicatedMeshDelayedReceive",
            false
        );

    if
    (
        replicatedMeshDelayedReceive
     && replicatedMesh_->asyncMigrationPending()
    )
    {
        replicatedMesh_->migrateFinish();
        replicatedMesh_->updateParticleCounts();
        clearMoveOrderedParcels();
    }

    //- Remove electrons
    if (findIndex(typeIdList_, "e-") != -1)
    {
        // TODO VINCENT: there is a clever way than rebuilding entire cell occ.
        removeElectrons();
        buildCellOccupancy();
    }

    //- Move the particles ballistically with their current velocities
    // Note: The parcels radial weighting factor (RWF) will stay constant over
    // the _entire_ move step. It will be updated by coordSystem().evolve (see
    // below). Any function that operates on the parcel in between has to
    // consider this. This is especially relevant for boundary measurements
    // that are performed when a parcel hits a wall during the move step.
    // Consider the following situation:
    //  1. parcel starts in cell = x
    //  2. parcel is moved to cell = y
    //  3. parcel hits a wall face y1 that belongs to cell = y. This hit now has
    //     to be counted with the RWF that the parcel had at the beginning of
    //     its move step, i.e. RWF(cell = x), _neither_ RWF(cell = y) _nor_
    //     RWF(face = y1)

    //scalar timer = mesh_.time().elapsedCpuTime();
    const steadyWallClock::time_point moveWallStart = steadyWallClock::now();
    const scalar moveCpuStart =
        profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;
    Cloud<dsmcParcel>::move(td, deltaTValue());
    if (openmpEnabled_ && openmpMoveEnabled_)
    {
        endMoveAppendCapture();
    }
    if (td.moveDetailProfile)
    {
        moveDetailParcels_ += td.moveParcels;
        moveDetailTrackCalls_ += td.moveTrackCalls;
        moveDetailSameTetNoFaceHits_ += td.moveSameTetNoFaceHits;
        moveDetailInternalTetNoFaceHits_ += td.moveInternalTetNoFaceHits;
        moveDetailFaceHits_ += td.moveFaceHits;
        moveDetailProcessorHits_ += td.moveProcessorHits;
        moveDetailCyclicHits_ += td.moveCyclicHits;
        moveDetailPatchHits_ += td.movePatchHits;
        moveDetailStuckHits_ += td.moveStuckHits;
        moveDetailTrackWallTime_ += td.moveTrackWallTime;
        moveDetailTrackerWallTime_ += td.moveTrackerWallTime;
        moveDetailBoundaryWallTime_ += td.moveBoundaryWallTime;
    }
    if (profileSummary_)
    {
        profileMoveWall_ += elapsedWallSeconds(moveWallStart);
        profileMoveCpu_ += mesh_.time().elapsedCpuTime() - moveCpuStart;
    }
    //Info<< "move" << tab << mesh_.time().elapsedCpuTime() - timer << " s" << endl;

    if (replicatedMeshActive())
    {
        if
        (
            replicatedMesh_->stepCounter() == 0
         || (
                replicatedMesh_->migrateInterval() > 0
             && replicatedMesh_->stepCounter() % replicatedMesh_->migrateInterval() == 0
            )
        )
        {
            if (replicatedMeshDelayedReceive)
            {
                replicatedMesh_->migrateBegin();
                replicatedMesh_->migrateFinish();
            }
            else
            {
                replicatedMesh_->migrateParticlesByCellOwner();
            }
            replicatedMesh_->updateParticleCounts();
            clearMoveOrderedParcels();
        }

        replicatedMesh_->advanceStepCounter();
    }

    if (replicatedMeshActive() && replicatedMesh_->rebalanceSteps().size())
    {
        const label currentStep = replicatedMesh_->stepCounter();
        const labelList& steps = replicatedMesh_->rebalanceSteps();

        forAll(steps, i)
        {
            if (steps[i] == currentStep)
            {
                if (replicatedMesh_->asyncMigrationPending())
                {
                    replicatedMesh_->migrateFinish();
                    replicatedMesh_->updateParticleCounts();
                    clearMoveOrderedParcels();
                }

                Info<< nl
                    << "Replicated mesh: manual cell owner reassignment at step "
                    << currentStep << nl << endl;
                replicatedMesh_->reassignCellOwner();
                replicatedMesh_->migrateParticlesByCellOwner();
                replicatedMesh_->updateParticleCounts();
                clearMoveOrderedParcels();
                break;
            }
        }
    }

    //- Update cell occupancy
    //timer = mesh_.time().elapsedCpuTime();
    buildCellOccupancy();
    //Info<< "buildCellOccupancy" << tab << mesh_.time().elapsedCpuTime() - timer << " s " << endl;

    //- Add electrons back after the move function
    if (findIndex(typeIdList_, "e-") != -1)
    {
        // TODO VINCENT: there is a clever way than rebuilding entire cell occ.
        addElectrons();
        buildCellOccupancy();
    }

    //- Radial weighting for non-Cartesian flows (e.g., axisymmetric). This is
    // where parcels will receive their new RWF and will possibly be cloned or
    // deleted.
    coordSystem().evolve();

    controllers_.controlBeforeCollisions();
    boundaries_.controlBeforeCollisions();

    //- Calculate new velocities via stochastic collisions
    //timer = mesh_.time().elapsedCpuTime();
    const steadyWallClock::time_point collisionWallStart =
        steadyWallClock::now();
    const scalar collisionCpuStart =
        profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;
    collisions();
    if (profileSummary_)
    {
        profileCollisionWall_ += elapsedWallSeconds(collisionWallStart);
        profileCollisionCpu_ +=
            mesh_.time().elapsedCpuTime() - collisionCpuStart;
    }
    //Info<< "collisions" << tab << mesh_.time().elapsedCpuTime() - timer << " s" << endl;

    //- Reactions may have changed cell occupancy, update if any reaction
    if (reactions_.nReactions() != 0)
    {
        buildCellOccupancy();
    }

    controllers_.controlAfterCollisions();
    boundaries_.controlAfterCollisions();

    if (replicatedMeshActive())
    {
        if (replicatedMesh_->asyncMigrationPending())
        {
            replicatedMesh_->migrateFinish();
            replicatedMesh_->updateParticleCounts();
            clearMoveOrderedParcels();
            buildCellOccupancy();
        }

        const label prevRebalances = replicatedMesh_->autoRebalanceCount();
        replicatedMesh_->addEvolveTime
        (
            mesh_.time().elapsedCpuTime() - moveAndCollideCpuStart
        );
        replicatedMesh_->autoRebalance();

        if (replicatedMesh_->autoRebalanceCount() > prevRebalances)
        {
            clearMoveOrderedParcels();
            buildCellOccupancy();
        }
    }

    if (profileSummary_)
    {
        profileMoveAndCollideWall_ +=
            elapsedWallSeconds(moveAndCollideWallStart);
        profileMoveAndCollideCpu_ +=
            mesh_.time().elapsedCpuTime() - moveAndCollideCpuStart;
    }
}


void Foam::dsmcCloud::evolve_fields()
{
    const steadyWallClock::time_point wallStart = steadyWallClock::now();
    const scalar cpu0 = profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;

    reactions_.outputData();

    fields_.calculateFields();

    //timer = mesh_.time().elapsedCpuTime();
    fields_.writeFields();
    //Info<< "fields W" << tab << mesh_.time().elapsedCpuTime() - timer << " s" << endl;

    controllers_.calculateProps();
    controllers_.outputResults();

    boundaries_.calculateProps();
    boundaries_.outputResults();

    boundaryMeas_.outputResults();

    trackingInfo_.clean();
    boundaryMeas_.clean();
    cellMeas_.clean();

    if (profileSummary_)
    {
        profilePostFieldsWall_ += elapsedWallSeconds(wallStart);
        profilePostFieldsCpu_ += mesh_.time().elapsedCpuTime() - cpu0;
    }
}


void Foam::dsmcCloud::rebuildMoveOrderedParcels()
{
    moveOrderedParcels_.setSize(this->size());

    label i = 0;
    forAllIter(dsmcCloud, *this, iter)
    {
        moveOrderedParcels_[i++] = &iter();
    }

    moveOrderedParcels_.setSize(i);
    const label nThreads = max(ompNumThreads_, label(1));
    moveOrderedThreadOffsets_.setSize(nThreads + 1);

    for (label threadI = 0; threadI <= nThreads; ++threadI)
    {
        moveOrderedThreadOffsets_[threadI] = threadI*i/nThreads;
    }
    moveOrderedParcelsValid_ = true;
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::storeMoveOrderedParcels
(
    const List<dsmcParcel*>& parcels,
    const labelList& threadOffsets
)
{
    moveOrderedParcels_ = parcels;
    moveOrderedThreadOffsets_ = threadOffsets;
    moveOrderedParcelsValid_ = true;
}


void Foam::dsmcCloud::transferMoveOrderedParcels
(
    List<dsmcParcel*>& parcels,
    const labelList& threadOffsets
)
{
    moveOrderedParcels_.transfer(parcels);
    moveOrderedThreadOffsets_ = threadOffsets;
    moveOrderedParcelsValid_ = true;
}


void Foam::dsmcCloud::beginMoveAppendCapture()
{
    if (!moveAppendCaptureActive_)
    {
        moveAppendedParcels_.clear();
    }

    moveAppendCaptureActive_ = true;
}


void Foam::dsmcCloud::endMoveAppendCapture()
{
    moveAppendCaptureActive_ = false;
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::recordMoveAppendedParcel(dsmcParcel* pPtr)
{
    if (!pPtr || !moveAppendCaptureActive_)
    {
        return;
    }

    #ifdef _OPENMP
    if (openmpEnabled_ && openmpMoveEnabled_ && omp_in_parallel())
    {
        #pragma omp critical(moveAppendedParcelsAppend)
        {
            moveAppendedParcels_.append(pPtr);
        }
        return;
    }
    #endif

    moveAppendedParcels_.append(pPtr);
}


void Foam::dsmcCloud::clearMoveAppendedParcels()
{
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::setMoveOrderedParcels
(
    const DynamicList<dsmcParcel*>& parcels
)
{
    moveOrderedParcels_.setSize(parcels.size());

    forAll(parcels, i)
    {
        moveOrderedParcels_[i] = parcels[i];
    }

    const label nThreads = max(ompNumThreads_, label(1));
    moveOrderedThreadOffsets_.setSize(nThreads + 1);

    for (label threadI = 0; threadI <= nThreads; ++threadI)
    {
        moveOrderedThreadOffsets_[threadI] =
            threadI*moveOrderedParcels_.size()/nThreads;
    }
    moveOrderedParcelsValid_ = true;
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::appendBatchToMoveOrdered
(
    const DynamicList<dsmcParcel*>& parcels
)
{
    const label oldSize = moveOrderedParcels_.size();
    moveOrderedParcels_.setSize(oldSize + parcels.size());

    forAll(parcels, i)
    {
        moveOrderedParcels_[oldSize + i] = parcels[i];
    }

    const label nThreads = max(ompNumThreads_, label(1));
    moveOrderedThreadOffsets_.setSize(nThreads + 1);

    for (label threadI = 0; threadI <= nThreads; ++threadI)
    {
        moveOrderedThreadOffsets_[threadI] =
            threadI*moveOrderedParcels_.size()/nThreads;
    }

    moveOrderedParcelsValid_ = true;
}


void Foam::dsmcCloud::deleteParcel(dsmcParcel* p)
{
    if (p)
    {
        deleteParticle(*p);
        moveOrderedParcelsValid_ = false;
        moveOrderedThreadOffsets_.clear();
        occupancyOrderedParcelsValid_ = false;
        cellOccupancyMaterialized_ = false;
    }
}


void Foam::dsmcCloud::printProfileSummary() const
{
    if (!profileSummary_)
    {
        return;
    }

    label steps = profileSteps_;
    scalar fullEvolveWall = profileFullEvolveWall_;
    scalar moveAndCollideWall = profileMoveAndCollideWall_;
    scalar moveWall = profileMoveWall_;
    scalar buildCellOccupancyWall = profileBuildCellOccupancyWall_;
    scalar collisionWall = profileCollisionWall_;
    scalar postFieldsWall = profilePostFieldsWall_;
    scalar fullEvolveCpu = profileFullEvolveCpu_;
    scalar moveAndCollideCpu = profileMoveAndCollideCpu_;
    scalar moveCpu = profileMoveCpu_;
    scalar buildCellOccupancyCpu = profileBuildCellOccupancyCpu_;
    scalar collisionCpu = profileCollisionCpu_;
    scalar postFieldsCpu = profilePostFieldsCpu_;
    label moveDetailParcels = moveDetailParcels_;
    label moveDetailTrackCalls = moveDetailTrackCalls_;
    label moveDetailSameTetNoFaceHits = moveDetailSameTetNoFaceHits_;
    label moveDetailInternalTetNoFaceHits = moveDetailInternalTetNoFaceHits_;
    label moveDetailFaceHits = moveDetailFaceHits_;
    label moveDetailProcessorHits = moveDetailProcessorHits_;
    label moveDetailCyclicHits = moveDetailCyclicHits_;
    label moveDetailPatchHits = moveDetailPatchHits_;
    label moveDetailStuckHits = moveDetailStuckHits_;
    scalar moveDetailTrackWallTime = moveDetailTrackWallTime_;
    scalar moveDetailTrackerWallTime = moveDetailTrackerWallTime_;
    scalar moveDetailBoundaryWallTime = moveDetailBoundaryWallTime_;

    const bool replicatedRawMpi = replicatedMeshActive() && !Pstream::parRun();

    if (replicatedRawMpi)
    {
        int mpiInit = 0;
        MPI_Initialized(&mpiInit);

        if (mpiInit)
        {
            MPI_Allreduce
            (
                MPI_IN_PLACE,
                &steps,
                1,
                MPI_INT,
                MPI_MAX,
                MPI_COMM_WORLD
            );

            scalar profileValues[] =
            {
                fullEvolveWall,
                moveAndCollideWall,
                moveWall,
                buildCellOccupancyWall,
                collisionWall,
                postFieldsWall,
                fullEvolveCpu,
                moveAndCollideCpu,
                moveCpu,
                buildCellOccupancyCpu,
                collisionCpu,
                postFieldsCpu
            };

            MPI_Allreduce
            (
                MPI_IN_PLACE,
                profileValues,
                12,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD
            );

            fullEvolveWall = profileValues[0];
            moveAndCollideWall = profileValues[1];
            moveWall = profileValues[2];
            buildCellOccupancyWall = profileValues[3];
            collisionWall = profileValues[4];
            postFieldsWall = profileValues[5];
            fullEvolveCpu = profileValues[6];
            moveAndCollideCpu = profileValues[7];
            moveCpu = profileValues[8];
            buildCellOccupancyCpu = profileValues[9];
            collisionCpu = profileValues[10];
            postFieldsCpu = profileValues[11];

            label moveDetailCounts[] =
            {
                moveDetailParcels,
                moveDetailTrackCalls,
                moveDetailSameTetNoFaceHits,
                moveDetailInternalTetNoFaceHits,
                moveDetailFaceHits,
                moveDetailProcessorHits,
                moveDetailCyclicHits,
                moveDetailPatchHits,
                moveDetailStuckHits
            };

            MPI_Allreduce
            (
                MPI_IN_PLACE,
                moveDetailCounts,
                9,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD
            );

            moveDetailParcels = moveDetailCounts[0];
            moveDetailTrackCalls = moveDetailCounts[1];
            moveDetailSameTetNoFaceHits = moveDetailCounts[2];
            moveDetailInternalTetNoFaceHits = moveDetailCounts[3];
            moveDetailFaceHits = moveDetailCounts[4];
            moveDetailProcessorHits = moveDetailCounts[5];
            moveDetailCyclicHits = moveDetailCounts[6];
            moveDetailPatchHits = moveDetailCounts[7];
            moveDetailStuckHits = moveDetailCounts[8];

            scalar moveDetailTimes[] =
            {
                moveDetailTrackWallTime,
                moveDetailTrackerWallTime,
                moveDetailBoundaryWallTime
            };

            MPI_Allreduce
            (
                MPI_IN_PLACE,
                moveDetailTimes,
                3,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD
            );

            moveDetailTrackWallTime = moveDetailTimes[0];
            moveDetailTrackerWallTime = moveDetailTimes[1];
            moveDetailBoundaryWallTime = moveDetailTimes[2];
        }
    }
    else if (Pstream::parRun())
    {
        reduce(steps, maxOp<label>());
        reduce(fullEvolveWall, maxOp<scalar>());
        reduce(moveAndCollideWall, maxOp<scalar>());
        reduce(moveWall, maxOp<scalar>());
        reduce(buildCellOccupancyWall, maxOp<scalar>());
        reduce(collisionWall, maxOp<scalar>());
        reduce(postFieldsWall, maxOp<scalar>());
        reduce(fullEvolveCpu, maxOp<scalar>());
        reduce(moveAndCollideCpu, maxOp<scalar>());
        reduce(moveCpu, maxOp<scalar>());
        reduce(buildCellOccupancyCpu, maxOp<scalar>());
        reduce(collisionCpu, maxOp<scalar>());
        reduce(postFieldsCpu, maxOp<scalar>());
        reduce(moveDetailParcels, sumOp<label>());
        reduce(moveDetailTrackCalls, sumOp<label>());
        reduce(moveDetailSameTetNoFaceHits, sumOp<label>());
        reduce(moveDetailInternalTetNoFaceHits, sumOp<label>());
        reduce(moveDetailFaceHits, sumOp<label>());
        reduce(moveDetailProcessorHits, sumOp<label>());
        reduce(moveDetailCyclicHits, sumOp<label>());
        reduce(moveDetailPatchHits, sumOp<label>());
        reduce(moveDetailStuckHits, sumOp<label>());
        reduce(moveDetailTrackWallTime, maxOp<scalar>());
        reduce(moveDetailTrackerWallTime, maxOp<scalar>());
        reduce(moveDetailBoundaryWallTime, maxOp<scalar>());
    }

    if ((replicatedRawMpi && isOutputRank()) || (!replicatedRawMpi && Pstream::master()))
    {
        const scalar accountedWall =
            moveWall
          + buildCellOccupancyWall
          + collisionWall
          + postFieldsWall;
        const scalar accountedCpu =
            moveCpu
          + buildCellOccupancyCpu
          + collisionCpu
          + postFieldsCpu;

        Info<< nl
            << "DSMC solver profile summary" << nl
            << "    solver profile steps          = " << steps << nl
            << "    move+collide wall [s]         = " << moveAndCollideWall << nl
            << "    move only [s]                 = " << moveWall << nl
            << "    buildCellOccupancy [s]        = " << buildCellOccupancyWall << nl
            << "    collision phase [s]           = " << collisionWall << nl
            << "    post fields/output [s]        = " << postFieldsWall << nl
            << "    total profiled [s]            = " << accountedWall << nl
            << "    full evolve wall [s]          = " << fullEvolveWall << nl
            << "    move+collide cpu [s]          = " << moveAndCollideCpu << nl
            << "    move only cpu [s]             = " << moveCpu << nl
            << "    buildCellOccupancy cpu [s]    = " << buildCellOccupancyCpu << nl
            << "    collision phase cpu [s]       = " << collisionCpu << nl
            << "    post fields/output cpu [s]    = " << postFieldsCpu << nl
            << "    total profiled cpu [s]        = " << accountedCpu << nl
            << "    full evolve cpu [s]           = " << fullEvolveCpu;

        if (profileDetail_)
        {
            Info<< nl
                << "    profile detail                = basic-stage timers only";
        }

        if (moveDetailParcels > 0 || moveDetailTrackCalls > 0)
        {
            Info<< nl
                << "    move detail parcels          = " << moveDetailParcels << nl
                << "    move detail track calls      = " << moveDetailTrackCalls << nl
                << "    move detail same-tet no-face = "
                << moveDetailSameTetNoFaceHits << nl
                << "    move detail internal tet only= "
                << moveDetailInternalTetNoFaceHits << nl
                << "    move detail face hits        = " << moveDetailFaceHits << nl
                << "    move detail processor hits   = " << moveDetailProcessorHits << nl
                << "    move detail cyclic hits      = " << moveDetailCyclicHits << nl
                << "    move detail patch hits       = " << moveDetailPatchHits << nl
                << "    move detail stuck hits       = " << moveDetailStuckHits << nl
                << "    move detail track max [s]    = "
                << moveDetailTrackWallTime << nl
                << "    move detail tracker max [s]  = "
                << moveDetailTrackerWallTime << nl
                << "    move detail boundary max [s] = "
                << moveDetailBoundaryWallTime;
        }

        Info<< nl
            << "    OpenMP enabled                = " << openmpEnabled_ << nl
            << "    OpenMP max threads            = " << ompNumThreads_ << nl
            << "    OpenMP move                   = " << openmpMoveEnabled_
            << " (" << openmpMoveSchedule_ << ", chunk "
            << openmpMoveChunk_ << ")" << nl
            << "    OpenMP collision              = "
            << openmpCollisionSchedule_ << ", chunk "
            << openmpCollisionChunk_;

        Info<< nl << endl;
    }

    if (replicatedMeshActive())
    {
        replicatedMesh_->report();
    }
}


Foam::label Foam::dsmcCloud::nTerminalOutputs()
{
    return nTerminalOutputs_;
}


void Foam::dsmcCloud::info()
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    const scalarList& iM = infoMeasurements();

    scalar nMol = iM[6];
    reduce(nMol, sumOp<scalar>());

    scalar linearKineticEnergy = iM[1];
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar rotationalEnergy = iM[2];
    reduce(rotationalEnergy, sumOp<scalar>());

    scalar vibrationalEnergy = iM[3];
    reduce(vibrationalEnergy, sumOp<scalar>());

    scalar electronicEnergy = iM[4];
    reduce(electronicEnergy, sumOp<scalar>());

    scalar stuckMolecules = iM[5];
    reduce(stuckMolecules, sumOp<scalar>());

    Info<< "    Number of DSMC particles        = "
        << nDsmcParticles
        << endl;

    if (nDsmcParticles > VSMALL)
    {
       Info << "    Number of stuck particles       = "
            << stuckMolecules/nParticle() << nl
            << "    Number of free particles        = "
            << nMol/nParticle() << nl
            << "    Average linear kinetic energy   = "
            << linearKineticEnergy/nMol << nl
            << "    Average rotational energy       = "
            << rotationalEnergy/nMol << nl
            << "    Average vibrational energy      = "
            << vibrationalEnergy/nMol << nl
            << "    Average electronic energy       = "
            << electronicEnergy/nMol << nl
            << "    Total energy                    = "
            << (linearKineticEnergy + rotationalEnergy
                + vibrationalEnergy + electronicEnergy)
            << endl;

        porousMeas().writePorousMeasurementsInfo();
    }
}


void Foam::dsmcCloud::loadBalanceCheck()
{
    dynamicLoadBalancing_.update();
}


void Foam::dsmcCloud::loadBalance(const label noRefinement)
{
    dynamicLoadBalancing_.perform(noRefinement);
}


void Foam::dsmcCloud::autoMap(const mapPolyMesh& mapper)
{
    dsmcParcel::trackingData td(*this);

    Cloud<dsmcParcel>::autoMap(td, mapper);

    coordSystem().dtModel().update();

    buildCellOccupancyFromScratch();
    buildCollisionSelectionRemainderFromScratch();
    resetBoundaries();
    resetMeasurementTools();
}


Foam::label Foam::dsmcCloud::randomLabel
(
    const label valOne,
    const label valTwo
)
{
    if (valOne == valTwo)
    {
        return valOne;
    }
    else
    {
        const label start = Foam::min(valOne, valTwo);
        const label end = Foam::max(valOne, valTwo);

        label val = start + label(rndGen_.sample01<scalar>()*(end - start + 1));

        // Rare case when scalar01() returns exactly 1.000 and the truncated
        // value would be out of range.
        if(val == end + 1)
        {
            val = randomLabel(start, end);
        }
        return val;
    }
}


Foam::vector Foam::dsmcCloud::equipartitionLinearVelocity
(
    const scalar temperature,
    const scalar mass
)
{
    return sqrt(physicoChemical::k.value()*temperature/mass)
        *rndGen_.GaussNormal<vector>();
}


Foam::vector Foam::dsmcCloud::chapmanEnskogVelocity
(
    const scalar temperature,
    const scalar mass,
    const vector& q,
    const tensor& tau
)
{
    const scalar B = max(mag(q), mag(tau));
    const scalar A = 1.0 + 30.0*B;

    bool repeatTry = true;

    vector CTry = vector::zero;

    while (repeatTry)
    {
        CTry = rndGen_.GaussNormal<vector>()/sqrt(2.0);

        const scalar gammaTry = 1.0 + (q & CTry)*(0.4*(CTry & CTry) - 1.0)
            - (CTry & (tau & CTry));

        if (gammaTry >= A*rndGen_.sample01<scalar>())
        {
            repeatTry = false;
        }
    }

    return CTry*sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}


void Foam::dsmcCloud::generalisedChapmanEnskog
(
    const label& typeID,
    const scalar& translationalTemperature,
    const scalar& rotationalTemperature,
    const scalar& vibrationalTemperature,
    const scalar& mass,
    const vector& D,
    const vector& qTra,
    const vector& qRot,
    const vector& qVib,
    const tensor& tau,
    scalar& ERot,
    labelList& vibLevel,
    vector& U
)
{
    const scalar kB = physicoChemical::k.value();
    const label nVibModes = constProps(typeID).nVibrationalModes();

    scalar B = max(mag(D), mag(tau));
    B = max(B, mag(qTra));
    B = max(B, mag(qRot));
    B = max(B, mag(qVib));
    const scalar A = 1.0 + 30.0*B;

    const scalarList& thetaV = constProps(typeID).thetaV();
    const scalar epsRotAv = rotationalTemperature/translationalTemperature;

    scalar epsVibAv = 0.0;
    if (nVibModes > 0 && vibrationalTemperature > 5.)
    {
        forAll(vibLevel, mode)
        {
            epsVibAv += thetaV[mode]/vibrationalTemperature
                /(exp(thetaV[mode]/vibrationalTemperature) - 1.0);
        }

    }

    vector CTry = vector::zero;

    scalar epsRot = 0.0;
    scalar epsVib = 0.0;

    bool repeatTry = true;

    while (repeatTry)
    {
        ERot =
            equipartitionRotationalEnergy
            (
                rotationalTemperature,
                constProps(typeID).rotationalDegreesOfFreedom()
            );

        epsRot = ERot/(kB*translationalTemperature);

        if (nVibModes > 0 && vibrationalTemperature > 5.)
        {
            vibLevel =
                equipartitionVibrationalEnergyLevel
                (
                    vibrationalTemperature,
                    nVibModes,
                    typeID
                );

            scalar epsVib = 0.0;
            forAll(vibLevel, mode)
            {
                epsVib += vibLevel[mode]*thetaV[mode];
            }

            epsVib /= vibrationalTemperature;
        }

        CTry = rndGen_.GaussNormal<vector>()/sqrt(2.0);

        const scalar gammaTry = 1.0 + 2.0*(D & CTry)
            + (qTra & CTry)*(0.4*(CTry & CTry) - 1.0)
            + (qRot & CTry)*(epsRot - epsRotAv)
            + (qVib & CTry)*(epsVib - epsVibAv)
            - (CTry & (tau & CTry));

        if (gammaTry >= A*rndGen_.sample01<scalar>())
        {
            repeatTry = false;
        }
    }

    U = CTry*sqrt(2.0*kB*translationalTemperature/mass);
}


Foam::scalar Foam::dsmcCloud::equipartitionRotationalEnergy
(
    const scalar temperature,
    const scalar rotationalDof
)
{
    scalar ERot = 0.0;

    if (rotationalDof < SMALL)
    {
        return ERot;
    }
    else if (rotationalDof < 2.0 + SMALL && rotationalDof > 2.0 - SMALL)
    {
        // Special case for rDof = 2, i.e. diatomics;
        ERot = -log(rndGen_.sample01<scalar>())*physicoChemical::k.value()*temperature;
    }
    else
    {
        scalar a = 0.5*rotationalDof - 1;

        scalar energyRatio;

        scalar P = -1;

        do
        {
            energyRatio = 10*rndGen_.sample01<scalar>();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < rndGen_.sample01<scalar>());

        ERot = energyRatio*physicoChemical::k.value()*temperature;
    }

    return ERot;
}


Foam::labelList Foam::dsmcCloud::equipartitionVibrationalEnergyLevel
(
    const scalar temperature,
    const label nVibrationalModes,
    const label typeId
)
{
    labelList vibLevel(nVibrationalModes, 0);

    if (nVibrationalModes == 0)
    {
        return vibLevel;
    }
    else
    {
        forAll(vibLevel, mode)
        {
            vibLevel[mode] = -log(rndGen_.sample01<scalar>())*temperature
                /constProps(typeId).thetaV_m(mode);
        }
    }

    return vibLevel;
}


Foam::label Foam::dsmcCloud::equipartitionElectronicLevel
(
    const scalar temperature,
    const labelList& electronicDegeneracyList,
    const scalarList& electronicEnergyList
)
{
    const scalar EMax = physicoChemical::k.value()*temperature;
    const label jMax = electronicDegeneracyList.size() - 1;

    //- Random integer between 0 and jMax
    label jDash = 0;
    //- Maximum possible electronic energy level within list based on k*TElec
    scalar EJ = 0.0;
    //- Maximum possible degeneracy level within list
    label gJ = 0;
    //- Selected intermediate integer electronic level (0 to jMax)
    label jSelect = 0;
    //- Maximum denominator value in Liechty pdf (see below)
    scalar expMax = 0.0;
    //- Summation term based on random electronic level
    scalar expSum = 0.0;
    //- Boltzmann distribution of Eq. 3.1.1 of Liechty thesis
    scalar boltz = 0.0;
    //- Distribution function Eq. 3.1.2 of Liechty thesis
    scalar func = 0.0;

    if (jMax > 0 and temperature > SMALL)
    {
        //- Calculate summation term in denominator of Eq. 3.1.1 in Liechty
        //  thesis
        forAll(electronicDegeneracyList, i)
        {
            expSum += electronicDegeneracyList[i]
                *exp(-electronicEnergyList[i]/EMax);
        }

        //- Select maximum integer energy level based on boltz value.
        //  Note that this depends on the temperature.
        scalar boltzMax = 0.0;

        forAll(electronicDegeneracyList, i)
        {
            //- Eq. 3.1.1 of Liechty thesis.
            boltz =
                electronicDegeneracyList[i]
               *exp(-electronicEnergyList[i]/EMax)
               /expSum;

            if (boltzMax < boltz)
            {
                boltzMax = boltz;
                jSelect = i;
            }
        }

        //- Max. poss energy in list: list goes from 0 to jMax
        EJ = electronicEnergyList[jSelect];
        //- Max. poss degeneracy in list: list goes from 0 to jMax
        gJ = electronicDegeneracyList[jSelect];
        //- Max. in denominator of Liechty pdf for initialisation/wall
        //  bcs/freestream EEle etc..
        expMax = gJ*exp(-EJ/EMax);

        //- Acceptance - rejection based on Eq. 3.1.2 of Liechty thesis
        do
        {
          //jDash = rndGen_.position<label>(0,jMax); OLD
            jDash = randomLabel(0, jMax);
            func =
                electronicDegeneracyList[jDash]
               *exp(-electronicEnergyList[jDash]/EMax)
               /expMax;
        } while(func < rndGen_.sample01<scalar>());
    }

    return jDash;
}


Foam::scalar Foam::dsmcCloud::postCollisionRotationalEnergy
(
    const scalar rotationalDof,
    const scalar ChiB
)
{
    scalar energyRatio = 0.0;

    if (rotationalDof == 2.0)
    {
        energyRatio = 1.0 - pow(rndGen_.sample01<scalar>(), 1.0/ChiB);
    }
    else
    {
        const scalar ChiA = 0.5*rotationalDof;

        scalar ChiAMinusOne = ChiA - 1.;

        scalar ChiBMinusOne = ChiB - 1.;

        if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
        {
            return rndGen_.sample01<scalar>();
        }

        scalar P = 0.0;

        do
        {
            P = 0;

            energyRatio = rndGen_.sample01<scalar>();

            if (ChiAMinusOne < SMALL)
            {
                P = pow(1.0 - energyRatio, ChiBMinusOne);
            }
            else if (ChiBMinusOne < SMALL)
            {
                P = pow(1.0 - energyRatio, ChiAMinusOne);
            }
            else
            {
                P =
                    pow
                    (
                        (ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne,
                        ChiAMinusOne
                    )
                *pow
                    (
                        (ChiAMinusOne + ChiBMinusOne)*(1 - energyRatio)
                        /ChiBMinusOne,
                        ChiBMinusOne
                    );
            }
        } while (P < rndGen_.sample01<scalar>());
    }

    return energyRatio;
}


Foam::label Foam::dsmcCloud::postCollisionVibrationalEnergyLevel
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
)
{
    label iDash = vibLevel;

    if (postReaction)
    {
        // post-collision quantum number
        scalar func = 0.0;
        scalar EVib = 0.0;

        do // acceptance - rejection
        {
            //iDash = rndGen_.position<label>(0, iMax); OLD
            iDash = randomLabel(0, iMax);
            EVib = iDash*physicoChemical::k.value()*thetaV;

            // - equation 5.61, Bird
            func = pow(1.0 - EVib/Ec, 1.5 - omega);

        } while(func < rndGen_.sample01<scalar>());
    }
    else
    {
        scalar inverseVibrationalCollisionNumber = 1.0;

        if (fixedZv == 0)
        {
            //- Temperature used to calculate Zv
            scalar T = 0;

            if (invZvFormulation == 0)
            {
                //- "Quantised collision temperature" (equation 3, Bird 2010)
                //  denominator from Bird 5.42
                T = iMax*thetaV/(3.5 - omega);
            }
            else if (invZvFormulation == 1)
            {
                //- Macroscopic (overall) temperature
                const scalar TMacro = fields().overallT(cellI);

                if (TMacro > SMALL)
                {
                    T = TMacro;
                }
                else
                {
                    //- Collision temperature used instead
                    //  the pre-2008 formulation is recovered
                    T = iMax*thetaV/(3.5 - omega);
                }

            }
            else
            {
                //- Macroscopic (translational) temperature
                /*const scalar TMacro = fields().translationalT(cellI);

                if (TMacro > SMALL)
                {
                    T = TMacro;
                }
                else
                {
                    //- Collision temperature used instead
                    //  the pre-2008 formulation is recovered
                    T = iMax*thetaV/(3.5 - omega);
                }*/ //TODO
                // Collision temperature for the time being
                // it gives a better agreement
                T = iMax*thetaV/(3.5 - omega);
            }

            const scalar pow1 = pow(thetaD/T, 1./3.) - 1.0;

            const scalar pow2 = pow(thetaD/refTempZv, 1./3.) - 1.0;

            //- vibrational collision number (equation 2, Bird 2010)
            const scalar ZvP1 = pow(thetaD/T, omega);

            const scalar ZvP2 =
                pow
                (
                    Zref*pow(thetaD/refTempZv, -omega),
                    pow1/pow2
                );

            const scalar Zv = ZvP1*ZvP2;

            //- In order to obtain the relaxation rate corresponding to Zv with the collision
            //  energy-based procedure, the inelastic fraction should be set to about 1/(5Zv)
            //  Bird 2008 RGD "A Comparison of Collision Energy-Based and Temperature-Based..."
            if (invZvFormulation == 2)
            {
                inverseVibrationalCollisionNumber = 1.0/(5.0*Zv);
            }
            else
            {
                inverseVibrationalCollisionNumber = 1.0/Zv;
            }
        }
        else
        {
            inverseVibrationalCollisionNumber = 1.0/fixedZv;
        }

        if (inverseVibrationalCollisionNumber > rndGen_.sample01<scalar>())
        {
            // post-collision quantum number
            scalar func = 0.0;
            scalar EVib = 0.0;

            do // acceptance - rejection
            {
                //iDash = rndGen_.position<label>(0, iMax); OLD
                iDash = randomLabel(0, iMax);

                EVib = iDash*physicoChemical::k.value()*thetaV;

                // - equation 5.61, Bird
                func = pow(1.0 - EVib/Ec, 1.5 - omega);

            } while(func < rndGen_.sample01<scalar>());
        }
    }

    return iDash;
}


Foam::label Foam::dsmcCloud::postCollisionElectronicEnergyLevel
(
    const scalar Ec,
    const label jMax,
    const scalar omega,
    const scalarList& EElist,
    const labelList& gList
)
{
    /*label nPossibleStates = 0;

    //- Post collision electronic level uniformly selected, taking the
    //  degeneracies of the different energy levels into account.

    //- Summation for all levels with energy below the collision energy
    if (jMax == 1)
    {
        nPossibleStates = gList[0];
    }
    else
    {
        forAll(EElist, i)
        {
            if (Ec > EElist[i])
            {
                nPossibleStates += gList[i];
            }
        }
    }

    label II = 0;

    //- Post-collision electronic energy
    label jDash = -1;

    do
    {
        const label nState = randomLabel(1, nPossibleStates);
        label nAvailableStates = 0;
        label nLevel = -1;

        forAll(EElist, n)
        {
            nAvailableStates += gList[n];

            if (nState <= nAvailableStates && nLevel < 0)
            {
                nLevel = n;
                break;
            }
        }

        //- Acceptance - rejection based on Eq. 3.1.8 of Liechty thesis
        if (Ec > EElist[nLevel])
        {
            scalar prob = pow(1.0 - EElist[nLevel]/Ec, 1.5 - omega);

            if (prob > rndGen_.sample01<scalar>())
            {
                II = 1;
                jDash = nLevel;
            }
        }

    } while (II == 0);

    return jDash;*/

    //- Maximum allowable electronic level obtainable from Ecoll
    label jSelectA = 0;
    //- Energy level maximazing expression gList[j]*pow(Ec - EElist[j], 1.5 - omega)
    label jSelectB = 0;

    scalar g = 0.0;
    scalar gMax = 0.0;

    //- Determine the maximum possible integer energy level immediately below Ec
    forAll(gList, i)
    {
       if (EElist[i] > Ec)
       {
           break;
       }

       jSelectA = i;

       //- Eq. 3.1.6 of Liechty thesis.
       g = gList[i]*pow(Ec - EElist[i], 1.5 - omega);

       if (gMax < g)
       {
           gMax = g;
           jSelectB = i;
       }
    }

    //- Liechty's procedure - step 3
    //  Minimum of jSelectA and jSelectB
    const label jSelect = min(jSelectA, jSelectB);

    //- Max. poss energy in list: list goes from 0 to jSelect
    const scalar EJ = EElist[jSelect];
    //- Max. poss degeneracy in list: list goes from 0 to jSelect
    const label gJ = gList[jSelect];
    //- Max. denominator of Liechty pdf for post-collision pdf
    const scalar denomMax = gJ*pow(Ec - EJ, 1.5 - omega);

    //- Acceptance - rejection based on Eq. 3.1.8 of Liechty thesis
    //- Post-collision electronic energy
    label jDash = 0;
    scalar prob = 0.0;

    do
    {
     //jDash = rndGen_.position<label>(0,jSelectA); OLD
       jDash = randomLabel(0, jSelectA);
       prob = gList[jDash]*pow(Ec - EElist[jDash], 1.5 - omega)/denomMax;

    } while(prob < rndGen_.sample01<scalar>());

    return jDash;

}


void Foam::dsmcCloud::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        pObj<< "v " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << nl;
    }

    pObj.flush();
}


void Foam::dsmcCloud::reBuildCellOccupancy()
{
    buildCellOccupancy();
}


void Foam::dsmcCloud::insertParcelInCellOccupancy(dsmcParcel* p)
{
    cellOccupancy_[p->cell()].append(p);
    cellOccupancy_[p->cell()].shrink();
}


void Foam::dsmcCloud::removeParcelFromCellOccupancy
(
    const label& cellMolId,
    const label& cell
)
{
    DynamicList<dsmcParcel*> molsInCell(0);

    forAll(cellOccupancy_[cell], c)
    {
        if (c != cellMolId)
        {
            molsInCell.append(cellOccupancy_[cell][c]);
        }
    }

    molsInCell.shrink();
    cellOccupancy_[cell].clear();
    cellOccupancy_[cell].transfer(molsInCell);
}


// NEW DANIEL *****************************************************************
/*void Foam::dsmcCloud::resetHybrid
(
    volScalarField& TtrInitial,
    volVectorField& UInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& qInitial,
    PtrList<volTensorField>& tauInitial,
    dimensionedScalar& B,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if (zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];

                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );

                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtrInitial[cellI];
                        scalar rotationalTemperature = TtrInitial[cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.sample01<scalar>()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = this->chapmanEnskogVelocityMiu
                            (
                                translationalTemperature,
                                cP.mass(),
                                B.value(),
                                qInitial[i][cellI],
                                tauInitial[i][cellI]
                            );

                            scalar ERot = this->equipartitionRotationalEnergy
                            (
                                rotationalTemperature,
                                cP.rotationalDegreesOfFreedom()
                            );

                            scalar EVib = this->equipartitionVibrationalEnergy
                            (
                                vibrationalTemperature,
                                cP.nVibrationalModes(),
                                i
                            );

                            U += velocity;

                            label newParcel = -1;

                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}

void Foam::dsmcCloud::resetHybrid2
(
    volScalarField& TtrInitial,
    volVectorField& UInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& qInitial,
    PtrList<volTensorField>& tauInitial,
    dimensionedScalar& B,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if (zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];

                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );

                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtrInitial[cellI];
                        scalar rotationalTemperature = TtrInitial[cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.sample01<scalar>()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = this->chapmanEnskogVelocity
                            (
                                translationalTemperature,
                                cP.mass(),
                                qInitial[i][cellI],
                                tauInitial[i][cellI]
                            );

                            scalar ERot = this->equipartitionRotationalEnergy
                            (
                                rotationalTemperature,
                                cP.rotationalDegreesOfFreedom()
                            );

                            scalar EVib = this->equipartitionVibrationalEnergy
                            (
                                vibrationalTemperature,
                                cP.nVibrationalModes(),
                                i
                            );

                            U += velocity;

                            label newParcel = -1;

                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}

// Hybrid Reset----------------------------------------------------------------
void Foam::dsmcCloud::resetHybridMax
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    for(label zoneToReset = 0; zoneToReset < zonesToReset.size(); zoneToReset++)
    {
        word regionName(zonesToReset[zoneToReset]);
        label zoneId = cellZones.findZoneID(regionName);

        if (zoneId == -1)
        {
            FatalErrorIn("resetHybridChapEnsk")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];

                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );

                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = TrInitial[i][cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.sample01<scalar>()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = this->equipartitionLinearVelocity
                            (
                                translationalTemperature,
                                cP.mass()
                            );

                            scalar ERot = this->equipartitionRotationalEnergy
                            (
                                rotationalTemperature,
                                cP.rotationalDegreesOfFreedom()
                            );

                            scalar EVib = this->equipartitionVibrationalEnergy
                            (
                                vibrationalTemperature,
                                cP.nVibrationalModes(),
                                i
                            );

                            U += velocity;

                            label newParcel = -1;

                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}


void Foam::dsmcCloud::resetHybridTra
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& qtInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    for(label zoneToReset = 0; zoneToReset < zonesToReset.size(); zoneToReset++)
    {
        word regionName(zonesToReset[zoneToReset]);
        label zoneId = cellZones.findZoneID(regionName);

        if (zoneId == -1)
        {
            FatalErrorIn("resetHybridChapEnsk")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];

                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );

                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = TrInitial[i][cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.sample01<scalar>()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = this->chapmanEnskogVelocity
                            (
                                translationalTemperature,
                                cP.mass(),
                                qtInitial[i][cellI],
                                tauInitial[i][cellI]
                            );

                            scalar ERot = this->equipartitionRotationalEnergy
                            (
                                rotationalTemperature,
                                cP.rotationalDegreesOfFreedom()
                            );

                            scalar EVib = this->equipartitionVibrationalEnergy
                            (
                                vibrationalTemperature,
                                cP.nVibrationalModes(),
                                i
                            );

                            U += velocity;

                            label newParcel = -1;

                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}*/

void Foam::dsmcCloud::resetHybridTraRotVib
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& DInitial,
    PtrList<volVectorField>& qtInitial,
    PtrList<volVectorField>& qrInitial,
    PtrList<volVectorField>& qvInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    //Info << "Deleting (" << typeOfReset << ")" << endl;
    //scalar time1_ = mesh_.time().elapsedCpuTime();
    //molsToDelete(mesh_, *this, typeOfReset);
    //scalar time2_ = mesh_.time().elapsedCpuTime();
    //Pout << "Proc " << UPstream::myProcNo() << " deletion time: " << time2_ - time1_ << "s (" << time1_ << ", " << time2_ << ")" << nl << endl;

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        List<scalar> massToIntroduce(TvInitial.size(), 0.0);//////////////
        List<scalar> massIntroduced(TvInitial.size(), 0.0);///////////////
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if (zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "\nInserting particles in" << endl;

            forAll(zone, c)
            {
                const label cellI = zone[c];
                forAll(typeIdList_, i)  ////////////////////////////////
                {                       ////////////////////////////////
                    massToIntroduce[i] += numberDensitiesField[i][cellI]*mesh_.V()[cellI];
                }                       ////////////////////////////////
                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );

                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    const scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = TrInitial[i][cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.sample01<scalar>()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        massIntroduced[i] += nParticlesToInsert;///////////////

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = vector::zero;

                            scalar ERot = 0.0;

                            labelList vibLevel
                            (
                                constProps(i).thetaV().size(),
                                0
                            );

                            label ELevel = 0; // TODO by generalisedChapmanEnskog

                            generalisedChapmanEnskog
                            (
                                i,
                                translationalTemperature,
                                rotationalTemperature,
                                vibrationalTemperature,
                                constProps(i).mass(),
                                DInitial[i][cellI],
                                qtInitial[i][cellI],
                                qrInitial[i][cellI],
                                qvInitial[i][cellI],
                                tauInitial[i][cellI],
                                ERot,
                                vibLevel,
                                U
                            );

                            U += velocity;

                            label newParcel = -1;

                            label classification = 0;

                            const scalar RWF = coordSystem().RWF(cellI);

                            addNewParcel
                            (
                                p,
                                U,
                                RWF,
                                ERot,
                                ELevel,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification,
                                vibLevel
                            );
                        }
                    }
                }
            }
        }

        Info<< "      Zone: " << regionName << endl;

        forAll(typeIdList_, i)
        {
            const scalar mass = this->constProps(i).mass();
            massToIntroduce[i] *= mass * nParticle();
            massIntroduced[i] *= mass * nParticle();
            Info<< "        Specie " << typeIdList_[i]
                << ", mTI: "
                << massToIntroduce[i] << "; mI: " << massIntroduced[i]
                << " (" << 100.0 * (massIntroduced[i]
                / (massToIntroduce[i] + VSMALL) - 1.0) << "% off)" << endl;
        }

        Info<< endl;
    }

    buildCellOccupancy();
}

/*void Foam::dsmcCloud::resetHybridTraRotVib2
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& DInitial,
    PtrList<volVectorField>& qtInitial,
    PtrList<volVectorField>& qrInitial,
    PtrList<volVectorField>& qvInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    Info << "Deleting (" << typeOfReset << ")" << endl;
    molsToDeleteHybrid(mesh_, *this, typeOfReset);

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if (zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "Lattice in zone: " << regionName << endl;

            forAll(zone, c)
            {
                const label& cellI = zone[c];

                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );

                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        const dsmcParcel::constantProperties& cP = this->constProps(i);

                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = TrInitial[i][cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

                        // Only integer numbers of particles can be inserted
                        label nParticlesToInsert = label(particlesRequired);

                        // Add another particle with a probability proportional to the
                        // remainder of taking the integer part of particlesRequired
                        if
                        (
                            (particlesRequired - nParticlesToInsert)
                                > rndGen_.sample01<scalar>()
                        )
                        {
                            nParticlesToInsert++;
                        }

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U;

                            scalar ERot;
                            scalar EVib;

                            this->generalisedChapmanEnskog2
                            (
                                i,
                                translationalTemperature,
                                rotationalTemperature,
                                vibrationalTemperature,
                                cP.mass(),
                                DInitial[i][cellI],
                                qtInitial[i][cellI],
                                qrInitial[i][cellI],
                                qvInitial[i][cellI],
                                tauInitial[i][cellI],
                                ERot,
                                EVib,
                                U
                            );

                            U += velocity;

                            label newParcel = -1;

                            label classification = 0;

                            this->addNewParcel
                            (
                                p,
                                U,
                                ERot,
                                EVib,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification
                            );
                        }
                    }
                }
            }
        }
    }
    buildCellOccupancy();
}*/


void Foam::dsmcCloud::resetHybridWhenUpdated
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& DInitial,
    PtrList<volVectorField>& qtInitial,
    PtrList<volVectorField>& qrInitial,
    PtrList<volVectorField>& qvInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    word& zoneToReset
)
{
    /*Info << "Deleting (" << typeOfReset << ")" << endl;
    scalar time1_ = mesh_.time().elapsedCpuTime();
    molsToDeleteHybrid(mesh_, *this, typeOfReset);
    scalar time2_ = mesh_.time().elapsedCpuTime();

    const cellZoneMesh& cellZones = mesh_.cellZones();
    List<scalar> massToIntroduce(TvInitial.size(), 0.0);//////////////
    List<scalar> massIntroduced(TvInitial.size(), 0.0);///////////////

    label zoneId = cellZones.findZoneID(zoneToReset);

    if(zoneId == -1)
    {
        FatalErrorIn("resetHybrid")
            << "Cannot find region: " << zoneToReset << nl << "in: "
            << mesh_.time().constant()/"cellZones"
            << exit(FatalError);
    }

    const cellZone& zone = cellZones[zoneId];

    if (zone.size())
    {
        Info << "Lattice in zone: " << zoneToReset << endl;

        forAll(zone, c)
        {
            const label& cellI = zone[c];
            forAll(typeIdList_, i)  ////////////////////////////////
            {                       ////////////////////////////////
                massToIntroduce[i] += numberDensitiesField[i][cellI]
                    * mesh_.V()[cellI];
            }                       ////////////////////////////////
            List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
            (
                mesh_,
                cellI
            );

            forAll(cellTets, tetI)
            {
                const tetIndices& cellTetIs = cellTets[tetI];

                tetPointRef tet = cellTetIs.tet(mesh_);

                scalar tetVolume = tet.mag();

                forAll(typeIdList_, i)
                {
                    const dsmcParcel::constantProperties& cP = this->constProps(i);

                    scalar numberDensity = numberDensitiesField[i][cellI];
                    scalar translationalTemperature = TtInitial[i][cellI];
                    scalar rotationalTemperature = TrInitial[i][cellI];
                    scalar vibrationalTemperature = TvInitial[i][cellI];
                    vector velocity = UInitial[cellI];

                    // Calculate the number of particles required
                    scalar particlesRequired = numberDensity*tetVolume;

                    // Only integer numbers of particles can be inserted
                    label nParticlesToInsert = label(particlesRequired);

                    // Add another particle with a probability proportional to the
                    // remainder of taking the integer part of particlesRequired
                    if
                    (
                        (particlesRequired - nParticlesToInsert)
                            > rndGen_.scalar01()
                    )
                    {
                        nParticlesToInsert++;
                    }
                    massIntroduced[i] += nParticlesToInsert;///////////////
                    for (label pI = 0; pI < nParticlesToInsert; pI++)
                    {
                        point p = tet.randomPoint(rndGen_);

                        vector U;

                        scalar ERot;
                        scalar EVib;

                        this->generalisedChapmanEnskog
                        (
                            i,
                            translationalTemperature,
                            rotationalTemperature,
                            vibrationalTemperature,
                            cP.mass(),
                            DInitial[i][cellI],
                            qtInitial[i][cellI],
                            qrInitial[i][cellI],
                            qvInitial[i][cellI],
                            tauInitial[i][cellI],
                            ERot,
                            EVib,
                            U
                        );

                        U += velocity;

                        label newParcel = 0;

                        label classification = 0;

                        this->addNewParcel
                        (
                            p,
                            U,
                            ERot,
                            EVib,
                            cellI,
                            cellTetIs.face(),
                            cellTetIs.tetPt(),
                            i,
                            newParcel,
                            classification
                        );
                    }
                }
            }
        }
    }
    Info << "For zone " + zoneToReset + ":" << endl;
    forAll(typeIdList_, i)
    {
        scalar mass = this->constProps(i).mass();
        massToIntroduce[i] *= mass*nParticle();
        massIntroduced[i] *= mass*nParticle();
        Info << "  Specie " << i << ", mTI: "
            << massToIntroduce[i] << "; mI: " << massIntroduced[i]
            << nl << "\t(" << 100.0 * (massIntroduced[i]
            / (massToIntroduce[i] + VSMALL) - 1.0) << "% off)" << endl;
    }

    buildCellOccupancy();*/
}


void Foam::dsmcCloud::shockReset()
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    const IOdictionary& shockDict
    (
        IOobject
        (
            "shockDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar initialParticles = readScalar(shockDict.lookup("initialParticles"));
    scalar maxPercentageParticles = readScalar(shockDict.lookup("maxPercentageParticles"));
    scalar rhoLeft = readScalar(shockDict.lookup("rhoLeft"));
    scalar rhoRight = readScalar(shockDict.lookup("rhoRight"));
    scalar xLeft = readScalar(shockDict.lookup("xLeft"));
    scalar xRight = readScalar(shockDict.lookup("xRight"));
    scalar areaTube = readScalar(shockDict.lookup("areaTube"));

    scalar maxDeltaParticles = maxPercentageParticles * initialParticles / 100.0;

    scalar deltaParticles = nDsmcParticles - initialParticles;
    Info << "deltaParticles [%]: " << 100.0 * deltaParticles / initialParticles
        << nl << "updating:" << (mag(deltaParticles) >= maxDeltaParticles)
        << endl;

    scalar deltaXTonParticles = deltaParticles * this->constProps(0).mass()
        / areaTube / (rhoRight - rhoLeft);
    Info << "deltaX/nEqParticles: " << deltaXTonParticles << endl;

    if (mag(deltaParticles) >= maxDeltaParticles)
    {
        forAll(mesh_.cells(), cellI)
        {
            const List<dsmcParcel*>& parcelsInCell
                = cellOccupancy_[cellI];

            scalar deltaX = deltaXTonParticles*nParticles(cellI); // NEW VINCENT

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                // deltaX > 0: if position > xRight + deltaX => out of domain
                // deltaX < 0: if position < xLeft + deltaX => out of domain
                if
                (
                    (p->classification() < 1000)
                    &&
                    (
                        (p->position().x() < xLeft + deltaX)
                            || (p->position().x() > xRight + deltaX)
                    )
                )
                {
                    p->classification() += 2000;
                    dsmcParcel* pPtr = new dsmcParcel(*p);
                    addParticle(pPtr);
                    p->classification() -= 2000;
                }
            }
        }

        forAllIter(dsmcCloud, *this, iter)
        {
            dsmcParcel& p = iter();

            if (p.classification() > 1000)
            {
                p.classification() -= 2000;
            }
            else
            {
                p.position().x() += deltaXTonParticles*nParticles(p.cell()); // NEW VINCENT

                if ((p.position().x() > xRight) || (p.position().x() < xLeft))
                {
                    deleteParticle(p);
                }
            }
        }

        buildCellOccupancy();
    }
}
// END NEW DANIEL *************************************************************


bool Foam::dsmcCloud::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
