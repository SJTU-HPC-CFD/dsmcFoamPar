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

#include "dsmcCloud.H"
#include "constants.H"
#include "polyMeshTetDecomposition.H"
#include "tetPointRef.H"
#include "entry.H"
#include "wallPolyPatch.H"
#include "processorPolyPatch.H"
#include "zeroGradientFvPatchFields.H"
#include <mpi.h>
#include <chrono>
#include <unordered_set>

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

namespace
{
Foam::vector makeTangential(const Foam::fvMesh& mesh, const Foam::label faceI, const Foam::vector& n)
{
    const Foam::vector fC = mesh.faceCentres()[faceI];
    const Foam::face& f = mesh.faces()[faceI];

    Foam::vector t1 = fC - mesh.points()[f[0]];

    if (Foam::mag(t1) < Foam::VSMALL)
    {
        Foam::vector ref = Foam::vector(1, 0, 0);

        if (Foam::mag(n & ref) > 0.9)
        {
            ref = Foam::vector(0, 1, 0);
        }

        t1 = ref - (ref & n)*n;
    }

    t1 /= Foam::max(Foam::mag(t1), Foam::VSMALL);
    return t1;
}

Foam::barycentric randomTetCoordinates(Foam::Random& rnd)
{
    const Foam::scalar e0 = -Foam::log(Foam::max(rnd.sample01<Foam::scalar>(), Foam::VSMALL));
    const Foam::scalar e1 = -Foam::log(Foam::max(rnd.sample01<Foam::scalar>(), Foam::VSMALL));
    const Foam::scalar e2 = -Foam::log(Foam::max(rnd.sample01<Foam::scalar>(), Foam::VSMALL));
    const Foam::scalar e3 = -Foam::log(Foam::max(rnd.sample01<Foam::scalar>(), Foam::VSMALL));
    const Foam::scalar sum = e0 + e1 + e2 + e3;
    return Foam::barycentric(e0/sum, e1/sum, e2/sum, e3/sum);
}

Foam::barycentric randomFaceTetCoordinates(Foam::Random& rnd, const Foam::scalar eps)
{
    const Foam::scalar e0 = -Foam::log(Foam::max(rnd.sample01<Foam::scalar>(), Foam::VSMALL));
    const Foam::scalar e1 = -Foam::log(Foam::max(rnd.sample01<Foam::scalar>(), Foam::VSMALL));
    const Foam::scalar e2 = -Foam::log(Foam::max(rnd.sample01<Foam::scalar>(), Foam::VSMALL));
    const Foam::scalar sum = e0 + e1 + e2;
    const Foam::scalar scale = 1.0 - eps;
    return Foam::barycentric(eps, scale*e0/sum, scale*e1/sum, scale*e2/sum);
}
}

void Foam::dsmcCloud::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    const dictionary moleculeProperties(particleProperties_.subDict("moleculeProperties"));

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);
        Info<< "    " << id << endl;
        constProps_[i] = dsmcParcel::constantProperties(moleculeProperties.subDict(id));
    }
}


void Foam::dsmcCloud::buildCellOccupancy(const bool rebuildParticlePartition)
{
    using clock_type = std::chrono::steady_clock;

    const auto t0 = clock_type::now();
    const bool emitOccupancyDiagnostics =
        profilingDetailEnabled_ && emitStepDiagnostics_;
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;

    const bool moveOrderedReady =
        openmpEnabled_
     && openmpMoveEnabled_
     && !moveOrderedReuseDisabled_
     && moveOrderedParcelsValid_
     && moveOrderedThreadOffsets_.size() == ompNumThreads_ + 1;
    DynamicList<dsmcParcel*> filteredAppendedParcels;
    const DynamicList<dsmcParcel*>* appendedParcelsSource = &moveAppendedParcels_;
    label appendedParcels = moveAppendedParcels_.size();

    if (Pstream::parRun() && moveOrderedReady && appendedParcels > 0)
    {
        std::unordered_set<dsmcParcel*> remainingAppended;
        remainingAppended.reserve(appendedParcels*2);

        forAll(moveAppendedParcels_, i)
        {
            remainingAppended.insert(moveAppendedParcels_[i]);
        }

        forAll(moveOrderedParcels_, i)
        {
            remainingAppended.erase(moveOrderedParcels_[i]);
        }

        if (remainingAppended.size() != moveAppendedParcels_.size())
        {
            filteredAppendedParcels.setCapacity(remainingAppended.size());

            forAll(moveAppendedParcels_, i)
            {
                dsmcParcel* pPtr = moveAppendedParcels_[i];

                if (remainingAppended.erase(pPtr))
                {
                    filteredAppendedParcels.append(pPtr);
                }
            }

            appendedParcelsSource = &filteredAppendedParcels;
            appendedParcels = filteredAppendedParcels.size();
        }
    }

    const bool useMoveOrderedParcels =
        moveOrderedReady
     && moveOrderedParcels_.size() + appendedParcels == this->size();
    const label moveOrderedSizeDelta =
        moveOrderedParcels_.size() + appendedParcels - this->size();

    label cloudValidCellParcels = 0;
    label cloudInvalidCellParcels = 0;
    label moveOrderedValidCellParcels = 0;
    label moveOrderedInvalidCellParcels = 0;

    if (emitOccupancyDiagnostics)
    {
        forAllIter(dsmcCloud, *this, iter)
        {
            const label celli = iter().cell();

            if (celli >= 0 && celli < cellOccupancy_.size())
            {
                ++cloudValidCellParcels;
            }
            else
            {
                ++cloudInvalidCellParcels;
            }
        }

        forAll(moveOrderedParcels_, i)
        {
            const label celli = moveOrderedParcels_[i]->cell();

            if (celli >= 0 && celli < cellOccupancy_.size())
            {
                ++moveOrderedValidCellParcels;
            }
            else
            {
                ++moveOrderedInvalidCellParcels;
            }
        }

        forAll(moveAppendedParcels_, i)
        {
            const label celli = moveAppendedParcels_[i]->cell();

            if (celli >= 0 && celli < cellOccupancy_.size())
            {
                ++moveOrderedValidCellParcels;
            }
            else
            {
                ++moveOrderedInvalidCellParcels;
            }
        }
    }

    if (profilingDetailEnabled_)
    {
        buildOccupancyMoveOrderedParcelsSum_ += moveOrderedParcels_.size();
        buildOccupancyMoveAppendedParcelsSum_ += appendedParcels;
        buildOccupancyCloudSizeSum_ += this->size();
        buildOccupancyMoveOrderedSizeDeltaSum_ += moveOrderedSizeDelta;
        buildOccupancyMoveOrderedSizeDeltaMin_ =
            min(buildOccupancyMoveOrderedSizeDeltaMin_, moveOrderedSizeDelta);
        buildOccupancyMoveOrderedSizeDeltaMax_ =
            max(buildOccupancyMoveOrderedSizeDeltaMax_, moveOrderedSizeDelta);

        if (emitOccupancyDiagnostics)
        {
            buildOccupancyCloudValidCellSum_ += cloudValidCellParcels;
            buildOccupancyCloudInvalidCellSum_ += cloudInvalidCellParcels;
            buildOccupancyMoveOrderedValidCellSum_ += moveOrderedValidCellParcels;
            buildOccupancyMoveOrderedInvalidCellSum_ += moveOrderedInvalidCellParcels;
        }

        if (useMoveOrderedParcels)
        {
            ++buildOccupancyMoveOrderedHits_;
        }
        else
        {
            ++buildOccupancyFallbackHits_;

            if (!moveOrderedParcelsValid_)
            {
                ++buildOccupancyMoveOrderedInvalidHits_;
            }
            else if (moveOrderedThreadOffsets_.size() != ompNumThreads_ + 1)
            {
                ++buildOccupancyMoveOrderedOffsetMismatchHits_;
            }
            else if (moveOrderedParcels_.size() + appendedParcels != this->size())
            {
                ++buildOccupancyMoveOrderedSizeMismatchHits_;
            }
        }
    }

    const label nParcels =
        useMoveOrderedParcels
      ? moveOrderedParcels_.size() + appendedParcels
      : this->size();
    const label nCells = cellOccupancy_.size();

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
            appendedParcelsPtr = appendedParcelsSource;
            appendedThreadOffsets.setSize(ompNumThreads_ + 1, 0);
            for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
            {
                appendedThreadOffsets[threadI] =
                    threadI*appendedParcels/ompNumThreads_;
            }
            appendedThreadOffsets[ompNumThreads_] = appendedParcels;
        }
        else
        {
            gatheredParcels.setSize(nParcels);
            label parcelI = 0;

            forAllIter(dsmcCloud, *this, iter)
            {
                gatheredParcels[parcelI++] = &iter();
            }

            generatedThreadOffsets.setSize(ompNumThreads_ + 1, 0);
            for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
            {
                generatedThreadOffsets[threadI] = threadI*nParcels/ompNumThreads_;
            }
            generatedThreadOffsets[ompNumThreads_] = nParcels;

            parcelsPtr = &gatheredParcels;
            parcelThreadOffsetsPtr = &generatedThreadOffsets;
        }

        const List<dsmcParcel*>& parcels = *parcelsPtr;
        const labelList& parcelThreadOffsets = *parcelThreadOffsetsPtr;
        const bool useAppendedParcels = useMoveOrderedParcels && appendedParcelsPtr;

        const auto t1 = clock_type::now();

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
            DynamicList<label>& localActiveCells = occupancyThreadActiveCells_[threadI];

            if (localCounts.size() != nCells)
            {
                localCounts.setSize(nCells, 0);
            }
            else if (!useMoveOrderedParcels)
            {
                forAll(localCounts, celli)
                {
                    localCounts[celli] = 0;
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
                    (parcelThreadOffsets[threadI + 1] - parcelThreadOffsets[threadI])
                  + (useAppendedParcels
                    ? appendedThreadOffsets[threadI + 1] - appendedThreadOffsets[threadI]
                    : 0)
                )
            );
        }

        #pragma omp parallel num_threads(ompNumThreads_)
        {
            const label threadI = currentThreadId();
            labelList& localCounts = occupancyThreadCellCounts_[threadI];
            DynamicList<label>& localActiveCells = occupancyThreadActiveCells_[threadI];

            for (label i = parcelThreadOffsets[threadI]; i < parcelThreadOffsets[threadI + 1]; ++i)
            {
                const label celli = parcels[i]->cell();

                if (celli >= 0 && celli < nCells)
                {
                    label& count = localCounts[celli];

                    if (count == 0)
                    {
                        localActiveCells.append(celli);
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
                    const label celli = appended[i]->cell();

                    if (celli >= 0 && celli < nCells)
                    {
                        label& count = localCounts[celli];

                        if (count == 0)
                        {
                            localActiveCells.append(celli);
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
            const DynamicList<label>& activeCells = occupancyThreadActiveCells_[threadI];

            forAll(activeCells, activeI)
            {
                const label celli = activeCells[activeI];
                totalCounts[celli] += localCounts[celli];
            }
        }

        const auto t2 = clock_type::now();

        label activeCellCount = 0;
        label collisionCellCount = 0;

        for (label celli = 0; celli < nCells; ++celli)
        {
            const label count = totalCounts[celli];

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

        for (label celli = 0; celli < nCells; ++celli)
        {
            const label count = totalCounts[celli];

            if (count > 0)
            {
                occupancyActiveCells_[activeCellCount++] = celli;

                if (count > 1)
                {
                    occupancyCollisionCells_[collisionCellCount++] = celli;
                }
            }
        }

        occupancyCellOffsets_.setSize(nCells + 1, 0);
        for (label celli = 0; celli < nCells; ++celli)
        {
            occupancyCellOffsets_[celli + 1] =
                occupancyCellOffsets_[celli] + totalCounts[celli];
        }
        occupancyOrderedParcels_.resize(occupancyCellOffsets_.last());

        labelList nextCellOffsets(occupancyCellOffsets_);

        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            labelList& localCounts = occupancyThreadCellCounts_[threadI];
            const DynamicList<label>& activeCells = occupancyThreadActiveCells_[threadI];

            forAll(activeCells, activeI)
            {
                const label celli = activeCells[activeI];
                const label count = localCounts[celli];
                localCounts[celli] = nextCellOffsets[celli];
                nextCellOffsets[celli] += count;
            }
        }

        if (rebuildParticlePartition)
        {
            rebuildParticleLoadPartition
            (
                totalCounts,
                occupancyCellOffsets_.last()
            );

            if (emitOccupancyDiagnostics)
            {
                label assignedParcels = 0;
                label gapCells = 0;
                label overlapCells = 0;
                label prevEnd = 0;

                for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
                {
                    const label start =
                        max(min(particleLoadStart_[threadI], nCells), label(0));
                    const label end =
                        max(min(particleLoadEnd_[threadI], nCells), start);

                    if (start > prevEnd)
                    {
                        gapCells += start - prevEnd;
                    }
                    else if (start < prevEnd)
                    {
                        overlapCells += prevEnd - start;
                    }

                    for (label celli = start; celli < end; ++celli)
                    {
                        assignedParcels += totalCounts[celli];
                    }

                    prevEnd = max(prevEnd, end);
                }

                if (prevEnd < nCells)
                {
                    gapCells += nCells - prevEnd;
                }

                buildOccupancyPartitionAssignedParcelsSum_ += assignedParcels;
                buildOccupancyPartitionGapCellsSum_ += gapCells;
                buildOccupancyPartitionOverlapCellsSum_ += overlapCells;
            }
        }

        #pragma omp parallel num_threads(ompNumThreads_)
        {
            const label threadI = currentThreadId();
            labelList& localOffsets = occupancyThreadCellCounts_[threadI];

            for (label i = parcelThreadOffsets[threadI]; i < parcelThreadOffsets[threadI + 1]; ++i)
            {
                dsmcParcel* pPtr = parcels[i];
                const label celli = pPtr->cell();

                if (celli >= 0 && celli < nCells)
                {
                    const label slot = localOffsets[celli]++;
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
                    const label celli = pPtr->cell();

                    if (celli >= 0 && celli < nCells)
                    {
                        const label slot = localOffsets[celli]++;
                        occupancyOrderedParcels_[slot] = pPtr;
                    }
                }
            }
        }

        const auto t3 = clock_type::now();
        occupancyOrderedParcelsValid_ = true;

        if (emitOccupancyDiagnostics)
        {
            scalar localPairPotential = 0.0;

            for (label celli = 0; celli < nCells; ++celli)
            {
                const scalar nC = scalar(totalCounts[celli]);

                if (nC > 1.0)
                {
                    localPairPotential += 0.5*nC*(nC - 1.0);
                }
            }

            label globalOccupancyParcels = occupancyCellOffsets_.last();
            label globalActiveCells = occupancyActiveCells_.size();
            label globalCollisionCells = occupancyCollisionCells_.size();
            scalar globalPairPotential = localPairPotential;

            if (Pstream::parRun())
            {
                reduce(globalOccupancyParcels, sumOp<label>());
                reduce(globalActiveCells, sumOp<label>());
                reduce(globalCollisionCells, sumOp<label>());
                reduce(globalPairPotential, sumOp<scalar>());
            }

            if (Pstream::master())
            {
                Info<< "    BuildCellOccupancy step summary:" << nl
                    << "        source                    = "
                    << (useMoveOrderedParcels ? "moveOrdered" : "fallbackGather") << nl
                    << "        cloud valid parcels       = "
                    << cloudValidCellParcels << nl
                    << "        cloud invalid parcels     = "
                    << cloudInvalidCellParcels << nl
                    << "        move-order valid parcels  = "
                    << moveOrderedValidCellParcels << nl
                    << "        move-order invalid parcels= "
                    << moveOrderedInvalidCellParcels << nl
                    << "        occupancy parcels global  = "
                    << globalOccupancyParcels << nl
                    << "        active cells global       = "
                    << globalActiveCells << nl
                    << "        collision cells global    = "
                    << globalCollisionCells << nl
                    << "        pair potential global     = "
                    << globalPairPotential << nl
                    << endl;
            }
        }

        if (evolveProfileEnabled_)
        {
            buildOccupancyExtractWallTime_ += std::chrono::duration<scalar>(t1 - t0).count();
            buildOccupancyCountWallTime_ += std::chrono::duration<scalar>(t2 - t1).count();
            buildOccupancyAssembleWallTime_ += std::chrono::duration<scalar>(t3 - t2).count();
            ++buildOccupancyProfileCalls_;

        }

        pendingMoveParcels_.clear();

        return;
    }
    #endif

    // Sparse clear: only clear cells that had particles last step.
    // For replicated mesh with 60K cells but only ~7.5K active, this
    // avoids ~52K unnecessary DynamicList::clear() calls.
    if (occupancyActiveCells_.size() > 0
        && occupancyActiveCells_.size() < cellOccupancy_.size() / 2)
    {
        forAll(occupancyActiveCells_, i)
        {
            cellOccupancy_[occupancyActiveCells_[i]].clear();
        }
    }
    else
    {
        forAll(cellOccupancy_, celli)
        {
            cellOccupancy_[celli].clear();
        }
    }

    // Use linked list iteration (parcelArray disabled due to memory)
    forAllIter(dsmcCloud, *this, iter)
    {
        if (iter().cell() >= 0 && iter().cell() < cellOccupancy_.size())
        {
            cellOccupancy_[iter().cell()].append(&iter());
        }
    }

    labelList totalCounts(nCells, 0);
    occupancyCellOffsets_.setSize(nCells + 1, 0);
    for (label celli = 0; celli < nCells; ++celli)
    {
        totalCounts[celli] = cellOccupancy_[celli].size();
        occupancyCellOffsets_[celli + 1] =
            occupancyCellOffsets_[celli] + totalCounts[celli];
    }

    label activeCellCount = 0;
    label collisionCellCount = 0;

    for (label celli = 0; celli < nCells; ++celli)
    {
        const label count = totalCounts[celli];

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

    for (label celli = 0; celli < nCells; ++celli)
    {
        const label count = totalCounts[celli];

        if (count > 0)
        {
            occupancyActiveCells_[activeCellCount++] = celli;

            if (count > 1)
            {
                occupancyCollisionCells_[collisionCellCount++] = celli;
            }
        }
    }

    occupancyOrderedParcels_.resize(occupancyCellOffsets_.last());

    for (label celli = 0; celli < nCells; ++celli)
    {
        label offset = occupancyCellOffsets_[celli];
        const auto& cellParcels = cellOccupancy_[celli];

        forAll(cellParcels, i)
        {
            occupancyOrderedParcels_[offset++] = cellParcels[i];
        }
    }

    if (rebuildParticlePartition)
    {
        rebuildParticleLoadPartition
        (
            totalCounts,
            occupancyCellOffsets_.last()
        );
    }

    occupancyOrderedParcelsValid_ = true;
    cellOccupancyMaterialized_ = true;
    pendingMoveParcels_.clear();

    if (emitOccupancyDiagnostics)
    {
        scalar localPairPotential = 0.0;

        for (label celli = 0; celli < nCells; ++celli)
        {
            const scalar nC = scalar(totalCounts[celli]);

            if (nC > 1.0)
            {
                localPairPotential += 0.5*nC*(nC - 1.0);
            }
        }

        label globalOccupancyParcels = occupancyCellOffsets_.last();
        label globalActiveCells = occupancyActiveCells_.size();
        label globalCollisionCells = occupancyCollisionCells_.size();
        scalar globalPairPotential = localPairPotential;

        if (Pstream::parRun())
        {
            reduce(globalOccupancyParcels, sumOp<label>());
            reduce(globalActiveCells, sumOp<label>());
            reduce(globalCollisionCells, sumOp<label>());
            reduce(globalPairPotential, sumOp<scalar>());
        }

        if (Pstream::master())
        {
            Info<< "    BuildCellOccupancy step summary:" << nl
                << "        source                    = fallbackSerial" << nl
                << "        cloud valid parcels       = "
                << cloudValidCellParcels << nl
                << "        cloud invalid parcels     = "
                << cloudInvalidCellParcels << nl
                << "        move-order valid parcels  = "
                << moveOrderedValidCellParcels << nl
                << "        move-order invalid parcels= "
                << moveOrderedInvalidCellParcels << nl
                << "        occupancy parcels global  = "
                << globalOccupancyParcels << nl
                << "        active cells global       = "
                << globalActiveCells << nl
                << "        collision cells global    = "
                << globalCollisionCells << nl
                << "        pair potential global     = "
                << globalPairPotential << nl
                << endl;
        }
    }

    const auto t1 = clock_type::now();

    if (evolveProfileEnabled_)
    {
        buildOccupancyExtractWallTime_ += std::chrono::duration<scalar>(t1 - t0).count();
        ++buildOccupancyProfileCalls_;

    }
}


void Foam::dsmcCloud::beginMoveAppendCapture()
{
    if (!moveAppendCaptureActive_)
    {
        moveAppendedParcels_.clear();
        moveAppendToPending_ = false;
    }

    moveAppendCaptureActive_ = true;
}


void Foam::dsmcCloud::beginMoveDeferredAppendStage()
{
    if (!moveAppendCaptureActive_)
    {
        return;
    }

    moveAppendToPending_ = true;
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::endMoveAppendCapture()
{
    moveAppendCaptureActive_ = false;
    moveAppendToPending_ = false;
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::beginCollisionPhase()
{
    if (!openmpEnabled_) return;

    collisionPhaseActive_ = true;
    const label nThreads = max(ompNumThreads_, label(1));
    collisionNewParcels_.setSize(nThreads);
    forAll(collisionNewParcels_, i)
    {
        collisionNewParcels_[i].clear();
    }
}


void Foam::dsmcCloud::endCollisionPhase()
{
    if (!collisionPhaseActive_) return;

    collisionPhaseActive_ = false;

    label totalNew = 0;
    forAll(collisionNewParcels_, threadI)
    {
        totalNew += collisionNewParcels_[threadI].size();
    }

    if (totalNew == 0) return;

    forAll(collisionNewParcels_, threadI)
    {
        forAll(collisionNewParcels_[threadI], i)
        {
            dsmcParcel* pPtr = collisionNewParcels_[threadI][i];
            Cloud<dsmcParcel>::addParticle(pPtr);
        }
        collisionNewParcels_[threadI].clear();
    }

    buildCellOccupancy();
}


void Foam::dsmcCloud::recordMoveAppendedParcel(dsmcParcel* pPtr)
{
    if (!pPtr)
    {
        return;
    }

    if (moveAppendCaptureActive_)
    {
        if (moveAppendToPending_ && openmpEnabled_ && openmpMoveEnabled_ && Pstream::parRun())
        {
            #ifdef _OPENMP
            #pragma omp critical(pendingMoveParcelsAppend)
            #endif
            {
                pendingMoveParcels_.append(pPtr);
            }
        }
        else
        {
            #ifdef _OPENMP
            #pragma omp critical(moveAppendedParcelsAppend)
            #endif
            {
                moveAppendedParcels_.append(pPtr);
            }
        }
    }
    else if (openmpEnabled_ && openmpMoveEnabled_ && Pstream::parRun())
    {
        #ifdef _OPENMP
        #pragma omp critical(pendingMoveParcelsAppend)
        #endif
        {
            pendingMoveParcels_.append(pPtr);
        }
    }
}


void Foam::dsmcCloud::clearMoveAppendedParcels()
{
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::recordPendingMoveParcel(dsmcParcel* pPtr)
{
    if (!pPtr)
    {
        return;
    }

    #ifdef _OPENMP
    #pragma omp critical(pendingMoveParcelsAppend)
    #endif
    {
        pendingMoveParcels_.append(pPtr);
    }
}


void Foam::dsmcCloud::clearPendingMoveParcels()
{
    pendingMoveParcels_.clear();
}


void Foam::dsmcCloud::recordMoveCommitCounts
(
    const label extracted,
    const label survivors,
    const label transferred,
    const label deleted
)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    ++moveCommitCountCalls_;
    moveExtractedParcelsSum_ += extracted;
    moveSurvivorParcelsSum_ += survivors;
    moveTransferredParcelsSum_ += transferred;
    moveDeletedParcelsSum_ += deleted;
}


void Foam::dsmcCloud::recordMoveLoopPasses(const label nPasses)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    moveLoopPassesSum_ += nPasses;
    moveLoopPassesMax_ = max(moveLoopPassesMax_, nPasses);
}


void Foam::dsmcCloud::recordMoveDeferredParcels(const label nDeferred)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    moveDeferredParcelsSum_ += nDeferred;
    moveDeferredParcelsMax_ = max(moveDeferredParcelsMax_, nDeferred);
}


void Foam::dsmcCloud::recordMoveReceivedParcels(const label nReceived)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    moveReceivedParcelsSum_ += nReceived;
    moveReceivedParcelsMax_ = max(moveReceivedParcelsMax_, nReceived);
}


void Foam::dsmcCloud::recordMoveExtractDetail
(
    const scalar deferredWallTime,
    const scalar orderedReuseWallTime,
    const scalar fullScanWallTime,
    const label deferredParcels,
    const label orderedReuseParcels,
    const label fullScanParcels
)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    if (deferredWallTime > 0)
    {
        moveExtractDeferredWallTime_ += deferredWallTime;
        ++moveExtractDeferredPasses_;
        moveExtractDeferredParcelsSum_ += deferredParcels;
    }

    if (orderedReuseWallTime > 0)
    {
        moveExtractOrderedReuseWallTime_ += orderedReuseWallTime;
        ++moveExtractOrderedReusePasses_;
        moveExtractOrderedReuseParcelsSum_ += orderedReuseParcels;
    }

    if (fullScanWallTime > 0)
    {
        moveExtractFullScanWallTime_ += fullScanWallTime;
        ++moveExtractFullScanPasses_;
        moveExtractFullScanParcelsSum_ += fullScanParcels;
    }
}


void Foam::dsmcCloud::recordMoveFirstPassReuseCheck
(
    const bool hasOrdered,
    const bool offsetsOk,
    const label currentSize,
    const label priorSize,
    const label appendedSize
)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    ++moveFirstPassReuseChecks_;

    if (hasOrdered)
    {
        ++moveFirstPassHasOrderedCount_;
    }
    else
    {
        ++moveFirstPassMissingOrderedCount_;
    }

    if (!offsetsOk)
    {
        ++moveFirstPassOffsetMismatchCount_;
    }

    moveFirstPassCurrentSizeSum_ += currentSize;
    moveFirstPassPriorSizeSum_ += priorSize;
    moveFirstPassPendingSizeSum_ += appendedSize;

    const label sizeDelta = currentSize - priorSize;
    moveFirstPassSizeDeltaSum_ += sizeDelta;
    moveFirstPassDeltaMinusPendingSum_ += sizeDelta - appendedSize;

    if (sizeDelta == 0)
    {
        ++moveFirstPassSizeMatchCount_;
    }

    if (sizeDelta == appendedSize)
    {
        ++moveFirstPassPriorPlusPendingMatchCount_;
    }
}


void Foam::dsmcCloud::materializeCellOccupancy()
{
    if (cellOccupancyMaterialized_ || !occupancyOrderedParcelsValid_)
    {
        return;
    }

    forAll(cellOccupancy_, celli)
    {
        const label count =
            occupancyCellOffsets_[celli + 1] - occupancyCellOffsets_[celli];

        cellOccupancy_[celli].clear();

        if (count > 0)
        {
            if (cellOccupancy_[celli].capacity() < count)
            {
                cellOccupancy_[celli].setCapacity(count);
            }

            cellOccupancy_[celli].setSize(count);

            for (label i = 0; i < count; ++i)
            {
                cellOccupancy_[celli][i] =
                    occupancyOrderedParcels_[occupancyCellOffsets_[celli] + i];
            }
        }
    }

    cellOccupancyMaterialized_ = true;
}


void Foam::dsmcCloud::resetAfterMeshDistribution()
{
    const label nCells = mesh_.nCells();

    cellOccupancy_.setSize(nCells);
    forAll(cellOccupancy_, celli)
    {
        cellOccupancy_[celli].clear();
    }

    collisionSelectionRemainder_.setSize(nCells);
    forAll(collisionSelectionRemainder_, celli)
    {
        collisionSelectionRemainder_[celli] = rndGen_.sample01<scalar>();
    }

    selectedPairsPerCell_.setSize(nCells, 0.0);
    nCandidatesPerCell_.setSize(nCells, 0);
    collisionCandidateCells_.clear();

    clearMoveOrderedParcels();
    clearMoveAppendedParcels();
    clearPendingMoveParcels();
    moveOrderedReuseDisabled_ = false;

    occupancyOrderedParcels_.clear();
    occupancyCellOffsets_.setSize(nCells + 1, 0);
    occupancyActiveCells_.clear();
    occupancyCollisionCells_.clear();
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = true;
    occupancyThreadCellCounts_.clear();
    occupancyThreadActiveCells_.clear();

    particleLoadStart_.setSize(max(ompNumThreads_, label(1)), 0);
    particleLoadEnd_.setSize(max(ompNumThreads_, label(1)), nCells);
    collisionLoadStart_.setSize(max(ompNumThreads_, label(1)), 0);
    collisionLoadEnd_.setSize(max(ompNumThreads_, label(1)), nCells);
    particlePartitionStep_ = 0;
    particlePartitionLastRebuildStep_ = -1;
    particlePartitionReferenceLoad_.clear();
    collisionPartitionStep_ = 0;
    collisionPartitionLastRebuildStep_ = -1;
    collisionPartitionReferenceLoad_.clear();

    initOpenMPMoveGuardCells();
}


void Foam::dsmcCloud::rebuildCellOccupancyAfterMeshDistribution()
{
    buildCellOccupancy();
}


void Foam::dsmcCloud::refreshAfterMeshDistribution()
{
    refreshTrackerUsage();
    trackingInfo_.reset();
    boundaryMeas_.setInitialConfig();
    cellMeas_.reset();
    fields_.refreshAfterMeshDistribution();
    boundaries_.refreshAfterMeshDistribution();
    controllers_.refreshAfterMeshDistribution();

    if (reactions_)
    {
        reactions_->initialConfiguration();
    }

    if (collisionPartnerSelectionPtr_)
    {
        collisionPartnerSelectionPtr_->initialConfiguration();
    }
}


Foam::scalar Foam::dsmcCloud::initialSigmaTcRMax() const
{
    scalar sigmaTcRMax = SMALL;

    forAll(constProps_, i)
    {
        const scalar estimate =
            constProps_[i].sigmaT()
           *maxwellianMostProbableSpeed(300.0, constProps_[i].mass());

        sigmaTcRMax = max(sigmaTcRMax, estimate);
    }

    return max(sigmaTcRMax, SMALL);
}


bool Foam::dsmcCloud::speciesSelected(const label fieldI, const label typeId) const
{
    return fieldSpecs_[fieldI].speciesIds.found(typeId);
}


void Foam::dsmcCloud::readBoundaryModels()
{
    wallPatchToModelId_.setSize(mesh_.boundaryMesh().size(), -1);
    diffuseWalls_.clear();
    specularWallPatchIds_.clear();
    inflows_.clear();

    IOobject io
    (
        "boundariesDict",
        mesh_.time().system(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return;
    }

    IOdictionary boundariesDict(io);

    if (boundariesDict.found("dsmcPatchBoundaries"))
    {
        const PtrList<entry> patchEntries(Foam::hyCompat::lookup(Foam::hyCompat::lookup(boundariesDict, "dsmcPatchBoundaries")));

        forAll(patchEntries, i)
        {
            if (!patchEntries[i].isDict())
            {
                continue;
            }

            const dictionary& entryDict = patchEntries[i].dict();
            const word model(entryDict.get<word>("boundaryModel"));
            const dictionary& patchProps = entryDict.subDict("patchBoundaryProperties");
            const word patchName = patchProps.get<word>("patchName");
            const label patchId = mesh_.boundaryMesh().findPatchID(patchName);

            if (patchId < 0)
            {
                FatalIOErrorInFunction(boundariesDict)
                    << "Unknown wall patch " << patchName << nl
                    << exit(FatalIOError);
            }

            if (model == "dsmcSpecularWallPatch")
            {
                specularWallPatchIds_.append(patchId);
                continue;
            }

            if (model != "dsmcDiffuseWallPatch")
            {
                continue;
            }

            const dictionary& wallProps = entryDict.subDict("dsmcDiffuseWallPatchProperties");

            diffuseWallSpec spec;
            spec.patchName = patchName;
            spec.patchId = patchId;
            spec.velocity = wallProps.lookupOrDefault<vector>("velocity", vector::zero);
            spec.temperature = wallProps.lookupOrDefault<scalar>("temperature", 300.0);

            wallPatchToModelId_[spec.patchId] = diffuseWalls_.size();
            diffuseWalls_.append(spec);
        }
    }

    if (boundariesDict.found("dsmcGeneralBoundaries"))
    {
        const PtrList<entry> patchEntries(Foam::hyCompat::lookup(Foam::hyCompat::lookup(boundariesDict, "dsmcGeneralBoundaries")));

        forAll(patchEntries, i)
        {
            if (!patchEntries[i].isDict())
            {
                continue;
            }

            const dictionary& entryDict = patchEntries[i].dict();
            const word model(entryDict.get<word>("boundaryModel"));

            if (model != "dsmcFreeStreamInflowPatch")
            {
                continue;
            }

            const dictionary& generalProps = entryDict.subDict("generalBoundaryProperties");
            const dictionary& inflowProps = entryDict.subDict("dsmcFreeStreamInflowPatchProperties");

            inflowSpec spec;
            spec.patchName = generalProps.get<word>("patchName");
            spec.patchId = mesh_.boundaryMesh().findPatchID(spec.patchName);

            if (spec.patchId < 0)
            {
                FatalIOErrorInFunction(boundariesDict)
                    << "Unknown inflow patch " << spec.patchName << nl
                    << exit(FatalIOError);
            }

            const wordList molecules(inflowProps.get<wordList>("typeIds"));
            spec.speciesIds.setSize(molecules.size(), -1);
            spec.numberDensities.setSize(molecules.size(), 0.0);

            const dictionary& nd = inflowProps.subDict("numberDensities");

            forAll(molecules, j)
            {
                const label typeId = typeIdList_.find(molecules[j]);

                if (typeId < 0)
                {
                    FatalIOErrorInFunction(boundariesDict)
                        << "Unknown inflow species " << molecules[j] << nl
                        << exit(FatalIOError);
                }

                spec.speciesIds[j] = typeId;
                spec.numberDensities[j] = nd.get<scalar>(molecules[j]);
            }

            spec.velocity = inflowProps.lookupOrDefault<vector>("velocity", vector::zero);
            spec.translationalTemperature = inflowProps.lookupOrDefault<scalar>("translationalTemperature", 300.0);
            spec.rotationalTemperature = inflowProps.lookupOrDefault<scalar>("rotationalTemperature", spec.translationalTemperature);
            spec.vibrationalTemperature = inflowProps.lookupOrDefault<scalar>("vibrationalTemperature", spec.translationalTemperature);
            spec.electronicTemperature = inflowProps.lookupOrDefault<scalar>("electronicTemperature", 0.0);

            const polyPatch& pp = mesh_.boundaryMesh()[spec.patchId];
            spec.faces.setSize(pp.size(), -1);
            spec.cells.setSize(pp.size(), -1);
            spec.accumulatedParcelsToInsert.setSize(spec.speciesIds.size());

            forAll(spec.faces, f)
            {
                spec.faces[f] = pp.start() + f;
                spec.cells[f] = pp.faceCells()[f];
            }

            forAll(spec.accumulatedParcelsToInsert, m)
            {
                spec.accumulatedParcelsToInsert[m].setSize(pp.size(), 0.0);
            }

            inflows_.append(spec);
        }
    }
}


void Foam::dsmcCloud::readFieldSpecs()
{
    fieldSpecs_.clear();

    IOobject io
    (
        "fieldPropertiesDict",
        mesh_.time().system(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return;
    }

    IOdictionary fieldProperties(io);

    if (!fieldProperties.found("dsmcFields"))
    {
        return;
    }

    const PtrList<entry> fieldEntries(Foam::hyCompat::lookup(fieldProperties, "dsmcFields"));
    DynamicList<fieldSpec> specs;

    forAll(fieldEntries, i)
    {
        if (!fieldEntries[i].isDict())
        {
            continue;
        }

        const dictionary& entryDict = fieldEntries[i].dict();

        if (entryDict.get<word>("fieldModel") != "dsmcVolFields")
        {
            continue;
        }

        const dictionary& props = entryDict.subDict("dsmcVolFieldsProperties");
        fieldSpec spec;
        spec.fieldName = props.get<word>("fieldName");
        spec.measureMeanFreePath = props.lookupOrDefault<bool>("measureMeanFreePath", false);

        if (entryDict.found("timeProperties"))
        {
            const dictionary& timeProps = entryDict.subDict("timeProperties");
            spec.resetAtOutput = timeProps.lookupOrDefault<Switch>("resetAtOutput", false);
            spec.resetAtOutputUntilTime = timeProps.lookupOrDefault<scalar>("resetAtOutputUntilTime", -GREAT);
        }

        const wordList molecules(props.get<wordList>("typeIds"));
        spec.speciesIds.setSize(molecules.size(), -1);

        forAll(molecules, j)
        {
            const label typeId = typeIdList_.find(molecules[j]);

            if (typeId < 0)
            {
                FatalIOErrorInFunction(fieldProperties)
                    << "Unknown field species " << molecules[j] << nl
                    << exit(FatalIOError);
            }

            spec.speciesIds[j] = typeId;
        }

        specs.append(spec);
    }

    fieldSpecs_.transfer(specs);
}

void Foam::dsmcCloud::initialiseSpeciesBoundaryForceDensity()
{
    speciesBoundaryForceDensity_.setSize(fieldSpecs_.size());

    forAll(speciesBoundaryForceDensity_, fieldI)
    {
        speciesBoundaryForceDensity_[fieldI].setSize(mesh_.boundaryMesh().size());

        forAll(speciesBoundaryForceDensity_[fieldI], patchI)
        {
            speciesBoundaryForceDensity_[fieldI][patchI].setSize(mesh_.boundaryMesh()[patchI].size(), vector::zero);
        }
    }
}


void Foam::dsmcCloud::resetBoundaryForceDensity()
{
    forAll(speciesBoundaryForceDensity_, fieldI)
    {
        forAll(speciesBoundaryForceDensity_[fieldI], patchI)
        {
            speciesBoundaryForceDensity_[fieldI][patchI] = vector::zero;
        }
    }
}


void Foam::dsmcCloud::createFields()
{
    dsmcNFields_.setSize(fieldSpecs_.size());
    rhoNFields_.setSize(fieldSpecs_.size());
    rhoMFields_.setSize(fieldSpecs_.size());
    pFields_.setSize(fieldSpecs_.size());
    TtraFields_.setSize(fieldSpecs_.size());
    UMeanFields_.setSize(fieldSpecs_.size());
    fDFields_.setSize(fieldSpecs_.size());
    fieldAverageSteps_.setSize(fieldSpecs_.size(), 0);
    rhoNCumulative_.setSize(fieldSpecs_.size());
    rhoMCumulative_.setSize(fieldSpecs_.size());
    momentumCumulative_.setSize(fieldSpecs_.size());
    kineticCumulative_.setSize(fieldSpecs_.size());

    forAll(fieldSpecs_, fieldI)
    {
        const word& name = fieldSpecs_[fieldI].fieldName;

        dsmcNFields_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject("dsmcN_" + name, mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        rhoNFields_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject("rhoN_" + name, mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
                mesh_,
                dimensionedScalar("zero", dimDensity/dimMass, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        rhoMFields_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject("rhoM_" + name, mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
                mesh_,
                dimensionedScalar("zero", dimDensity, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        pFields_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject("p_" + name, mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
                mesh_,
                dimensionedScalar("zero", dimPressure, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        TtraFields_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject("Ttra_" + name, mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
                mesh_,
                dimensionedScalar("zero", dimTemperature, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        UMeanFields_.set
        (
            fieldI,
            new volVectorField
            (
                IOobject("UMean_" + name, mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
                mesh_,
                dimensionedVector("zero", dimVelocity, vector::zero),
                zeroGradientFvPatchVectorField::typeName
            )
        );

        fDFields_.set
        (
            fieldI,
            new volVectorField
            (
                IOobject("fD_" + name, mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
                mesh_,
                dimensionedVector("zero", dimPressure, vector::zero),
                zeroGradientFvPatchVectorField::typeName
            )
        );

        rhoNCumulative_[fieldI].setSize(mesh_.nCells(), 0.0);
        rhoMCumulative_[fieldI].setSize(mesh_.nCells(), 0.0);
        momentumCumulative_[fieldI].setSize(mesh_.nCells(), vector::zero);
        kineticCumulative_[fieldI].setSize(mesh_.nCells(), 0.0);
    }

    initialiseSpeciesBoundaryForceDensity();
}


Foam::vector Foam::dsmcCloud::equipartitionLinearVelocity
(
    const scalar temperature,
    const scalar mass
)
{
    return sqrt(physicoChemical::k.value()*temperature/mass)*rndGen_.GaussNormal<vector>();
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
        ERot = -log(rndGen_.sample01<scalar>())*physicoChemical::k.value()*temperature;
    }
    else
    {
        const scalar a = 0.5*rotationalDof - 1.0;
        scalar energyRatio = 0.0;
        scalar P = -1.0;

        do
        {
            energyRatio = 10.0*rndGen_.sample01<scalar>();
            P = pow(energyRatio/a, a)*exp(a - energyRatio);
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

    forAll(vibLevel, mode)
    {
        vibLevel[mode] = -log(max(rndGen_.sample01<scalar>(), VSMALL))*temperature/constProps(typeId).thetaV_m(mode);
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
    if (temperature <= SMALL || electronicDegeneracyList.size() <= 1)
    {
        return 0;
    }

    const scalar EMax = physicoChemical::k.value()*temperature;
    scalar expSum = 0.0;
    scalar expMax = 0.0;
    label jSelect = 0;

    forAll(electronicDegeneracyList, i)
    {
        expSum += electronicDegeneracyList[i]*exp(-electronicEnergyList[i]/EMax);
    }

    forAll(electronicDegeneracyList, i)
    {
        const scalar value = electronicDegeneracyList[i]*exp(-electronicEnergyList[i]/EMax)/expSum;

        if (value > expMax)
        {
            expMax = value;
            jSelect = i;
        }
    }

    const scalar denom = electronicDegeneracyList[jSelect]*exp(-electronicEnergyList[jSelect]/EMax);
    label jDash = 0;
    scalar func = 0.0;

    do
    {
        jDash = randomLabel(0, electronicDegeneracyList.size() - 1);
        func = electronicDegeneracyList[jDash]*exp(-electronicEnergyList[jDash]/EMax)/denom;
    } while (func < rndGen_.sample01<scalar>());

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
        const scalar ChiAMinusOne = ChiA - 1.0;
        const scalar ChiBMinusOne = ChiB - 1.0;

        if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
        {
            return rndGen_.sample01<scalar>();
        }

        scalar P = 0.0;

        do
        {
            energyRatio = rndGen_.sample01<scalar>();

            if (ChiAMinusOne < SMALL)
            {
                P = pow(1.0 - energyRatio, ChiBMinusOne);
            }
            else if (ChiBMinusOne < SMALL)
            {
                P = pow(energyRatio, ChiAMinusOne);
            }
            else
            {
                P =
                    pow((ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne, ChiAMinusOne)
                   *pow((ChiAMinusOne + ChiBMinusOne)*(1.0 - energyRatio)/ChiBMinusOne, ChiBMinusOne);
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
    const label
)
{
    label iDash = vibLevel;

    auto sampleLevel = [&](const label maxLevel)
    {
        scalar func = 0.0;
        scalar EVib = 0.0;

        do
        {
            iDash = randomLabel(0, maxLevel);
            EVib = iDash*physicoChemical::k.value()*thetaV;
            func = pow(max(1.0 - EVib/max(Ec, VSMALL), 0.0), 1.5 - omega);
        } while (func < rndGen_.sample01<scalar>());
    };

    if (postReaction)
    {
        sampleLevel(iMax);
        return iDash;
    }

    scalar inverseVibrationalCollisionNumber = 1.0;

    if (fixedZv == 0.0)
    {
        const scalar T = iMax*thetaV/(3.5 - omega);
        const scalar pow1 = pow(thetaD/max(T, VSMALL), 1.0/3.0) - 1.0;
        const scalar pow2 = pow(thetaD/max(refTempZv, VSMALL), 1.0/3.0) - 1.0;
        const scalar ZvP1 = pow(thetaD/max(T, VSMALL), omega);
        const scalar ZvP2 = pow(Zref*pow(thetaD/max(refTempZv, VSMALL), -omega), pow1/max(pow2, VSMALL));
        const scalar Zv = ZvP1*ZvP2;

        inverseVibrationalCollisionNumber = invZvFormulation == 2 ? 1.0/(5.0*Zv) : 1.0/Zv;
    }
    else
    {
        inverseVibrationalCollisionNumber = 1.0/fixedZv;
    }

    if (inverseVibrationalCollisionNumber > rndGen_.sample01<scalar>())
    {
        sampleLevel(iMax);
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
    if (jMax <= 0 || Ec <= SMALL)
    {
        return 0;
    }

    label jSelectA = 0;
    label jSelectB = 0;
    scalar gMax = 0.0;

    forAll(gList, i)
    {
        if (EElist[i] > Ec)
        {
            break;
        }

        jSelectA = i;

        const scalar g = gList[i]*pow(Ec - EElist[i], 1.5 - omega);

        if (g > gMax)
        {
            gMax = g;
            jSelectB = i;
        }
    }

    const label jSelect = min(jSelectA, jSelectB);
    const scalar denomMax = gList[jSelect]*pow(max(Ec - EElist[jSelect], VSMALL), 1.5 - omega);

    label jDash = 0;
    scalar prob = 0.0;

    do
    {
        jDash = randomLabel(0, jSelectA);
        prob = gList[jDash]*pow(max(Ec - EElist[jDash], 0.0), 1.5 - omega)/denomMax;
    } while (prob < rndGen_.sample01<scalar>());

    return jDash;
}


void Foam::dsmcCloud::noTimeCounterCollisions()
{
    if (!binaryCollision().active())
    {
        return;
    }

    const scalar deltaT = mesh_.time().deltaTValue();
    label collisionCandidates = 0;
    label collisions = 0;

    forAll(cellOccupancy_, celli)
    {
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy_[celli]);
        const label nC = cellParcels.size();

        if (nC < 2)
        {
            continue;
        }

        const scalar sigmaTcRMaxCell = max(sigmaTcRMax_[celli], SMALL);
        const scalar selectedPairs =
            collisionSelectionRemainder_[celli]
          + 0.5*nC*(nC - 1)*nParticle_*sigmaTcRMaxCell*deltaT/mesh_.cellVolumes()[celli];

        const label nCandidates = label(selectedPairs);
        collisionSelectionRemainder_[celli] = selectedPairs - nCandidates;
        collisionCandidates += nCandidates;

        for (label c = 0; c < nCandidates; ++c)
        {
            const label candidateP = rndGen_.position<label>(0, nC - 1);
            label candidateQ = candidateP;

            while (candidateQ == candidateP)
            {
                candidateQ = rndGen_.position<label>(0, nC - 1);
            }

            dsmcParcel& parcelP = *cellParcels[candidateP];
            dsmcParcel& parcelQ = *cellParcels[candidateQ];

            const scalar sigmaTcR = binaryCollision().sigmaTcR(parcelP, parcelQ);

            if (sigmaTcR > sigmaTcRMax_[celli])
            {
                sigmaTcRMax_[celli] = sigmaTcR;
            }

            if (sigmaTcR > 0 && (sigmaTcR/sigmaTcRMaxCell) > rndGen_.sample01<scalar>())
            {
                if (reactionsActive())
                {
                    const label reactionModelId = reactions().returnModelId(parcelP, parcelQ);

                    if (reactionModelId != -1)
                    {
                        reactions().reactions()[reactionModelId]->reaction(parcelP, parcelQ);

                        if (reactions().reactions()[reactionModelId]->relax())
                        {
                            binaryCollision().collide(parcelP, parcelQ, celli);
                        }
                    }
                    else
                    {
                        binaryCollision().collide(parcelP, parcelQ, celli);
                    }
                }
                else
                {
                    binaryCollision().collide(parcelP, parcelQ, celli);
                }

                ++collisions;
            }
        }
    }

    reduce(collisions, sumOp<label>());
    reduce(collisionCandidates, sumOp<label>());
    sigmaTcRMax_.correctBoundaryConditions();

    if (collisionCandidates)
    {
        Info<< "    Collisions                      = " << collisions << nl
            << "    Acceptance rate                 = "
            << scalar(collisions)/scalar(collisionCandidates) << nl << endl;
    }
}


void Foam::dsmcCloud::collisions()
{
    if (!collisionPartnerSelectionPtr_.valid())
    {
        FatalErrorInFunction
            << "collisionPartnerSelection model was not initialised" << nl
            << exit(FatalError);
    }

    using clock_type = std::chrono::steady_clock;

    const auto t0 = clock_type::now();
    precomputeCollisionCandidates();
    const auto t1 = clock_type::now();

    if (openmpEnabled_ && openmpCollisionSchedule_ == "partition")
    {
        rebuildCollisionLoadPartition();
    }
    const auto t2 = clock_type::now();
    collisionPartnerSelectionPtr_->collide();
    const auto t3 = clock_type::now();

    if (collisionProfileEnabled_)
    {
        collisionPrecomputeWallTime_ += std::chrono::duration<scalar>(t1 - t0).count();
        collisionPartitionWallTime_ += std::chrono::duration<scalar>(t2 - t1).count();
        collisionSelectionWallTime_ += std::chrono::duration<scalar>(t3 - t2).count();
        ++collisionProfileCalls_;

    }
}



void Foam::dsmcCloud::insertInflowParcels()
{
    const scalar sqrtPi = sqrt(pi);

    forAll(inflows_, inflowI)
    {
        inflowSpec& inflow = inflows_[inflowI];

        forAll(inflow.accumulatedParcelsToInsert, specieI)
        {
            const label typeId = inflow.speciesIds[specieI];
            const scalar mass = constProps(typeId).mass();
            const scalar mostProbableSpeed = maxwellianMostProbableSpeed(inflow.translationalTemperature, mass);

            forAll(inflow.faces, f)
            {
                const label faceI = inflow.faces[f];
                const vector sF = mesh_.faceAreas()[faceI];
                const scalar fA = mag(sF);
                const scalar sCosTheta = (inflow.velocity & (-sF/fA))/mostProbableSpeed;

                inflow.accumulatedParcelsToInsert[specieI][f] +=
                    fA*inflow.numberDensities[specieI]*mesh_.time().deltaTValue()*mostProbableSpeed
                   *(exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1.0 + erf(sCosTheta)))
                   /(2.0*sqrtPi*nParticle_);
            }
        }

        forAll(inflow.faces, f)
        {
            const label faceI = inflow.faces[f];
            const label cellI = inflow.cells[f];
            const vector sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);
            const vector n = -sF/fA;
            const vector t1 = makeTangential(mesh_, faceI, n);
            const vector t2 = (n ^ t1)/max(mag(n ^ t1), VSMALL);
            const List<tetIndices> faceTets(polyMeshTetDecomposition::faceTetIndices(mesh_, faceI, cellI));
            List<scalar> cTriAFracs(faceTets.size(), 0.0);
            scalar cumulative = 0.0;

            forAll(faceTets, triI)
            {
                cumulative += faceTets[triI].faceTri(mesh_).mag()/fA;
                cTriAFracs[triI] = cumulative;
            }

            if (cTriAFracs.size())
            {
                cTriAFracs.last() = 1.0;
            }

            forAll(inflow.speciesIds, specieI)
            {
                scalar& faceAccumulator = inflow.accumulatedParcelsToInsert[specieI][f];
                label nInsert = max(label(faceAccumulator), 0);

                if ((faceAccumulator - nInsert) > rndGen_.sample01<scalar>())
                {
                    ++nInsert;
                }

                faceAccumulator -= nInsert;

                const label typeId = inflow.speciesIds[specieI];
                const scalar mass = constProps(typeId).mass();
                const scalar mostProbableSpeed = maxwellianMostProbableSpeed(inflow.translationalTemperature, mass);
                const scalar sCosTheta = (inflow.velocity & n)/mostProbableSpeed;
                const scalar coeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);
                const scalar coeffB = 0.5*(1.0 + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0)));
                scalar randomScaling = 3.0;

                if (sCosTheta < -3.0)
                {
                    randomScaling = mag(sCosTheta) + 1.0;
                }

                for (label i = 0; i < nInsert; ++i)
                {
                    const scalar triSelection = rndGen_.sample01<scalar>();
                    label selectedTriI = 0;

                    forAll(cTriAFracs, triI)
                    {
                        selectedTriI = triI;

                        if (cTriAFracs[triI] >= triSelection)
                        {
                            break;
                        }
                    }

                    const tetIndices& faceTetIs = faceTets[selectedTriI];
                    const scalar eps = 1.0e-3;
                    const barycentric coordinates = randomFaceTetCoordinates(rndGen_, eps);
                    scalar uNormal = 0.0;
                    scalar P = -1.0;

                    if (mag(inflow.velocity & n) > VSMALL)
                    {
                        do
                        {
                            const scalar uNormalThermal = randomScaling*(2.0*rndGen_.sample01<scalar>() - 1.0);
                            uNormal = uNormalThermal + sCosTheta;

                            if (uNormal < 0.0)
                            {
                                P = -1.0;
                            }
                            else
                            {
                                P = 2.0*uNormal/coeffA*exp(coeffB - sqr(uNormalThermal));
                            }
                        } while (P < rndGen_.sample01<scalar>());
                    }
                    else
                    {
                        uNormal = sqrt(-log(max(rndGen_.sample01<scalar>(), VSMALL)));
                    }

                    vector U =
                        sqrt(physicoChemical::k.value()*inflow.translationalTemperature/mass)
                       *(rndGen_.GaussNormal<scalar>()*t1 + rndGen_.GaussNormal<scalar>()*t2)
                      + (t1 & inflow.velocity)*t1
                      + (t2 & inflow.velocity)*t2
                      + mostProbableSpeed*uNormal*n;

                    const scalar ERot = equipartitionRotationalEnergy(inflow.rotationalTemperature, constProps(typeId).rotationalDegreesOfFreedom());
                    const labelList vibLevel = equipartitionVibrationalEnergyLevel(inflow.vibrationalTemperature, constProps(typeId).nVibrationalModes(), typeId);
                    const label ELevel = equipartitionElectronicLevel(inflow.electronicTemperature, constProps(typeId).electronicDegeneracyList(), constProps(typeId).electronicEnergyList());

                    addNewParcel(coordinates, U, 1.0, ERot, ELevel, cellI, faceI, faceTetIs.tetPt(), typeId, inflow.patchId, 0, vibLevel);
                }
            }
        }
    }
}


void Foam::dsmcCloud::sampleFields()
{
    resetBoundaryForceDensity();

    forAll(fieldSpecs_, fieldI)
    {
        if (fieldSpecs_[fieldI].pendingReset)
        {
            fieldAverageSteps_[fieldI] = 0;
            rhoNCumulative_[fieldI] = 0.0;
            rhoMCumulative_[fieldI] = 0.0;
            momentumCumulative_[fieldI] = vector::zero;
            kineticCumulative_[fieldI] = 0.0;
            fieldSpecs_[fieldI].pendingReset = false;
        }

        scalarField parcelCount(mesh_.nCells(), 0.0);
        scalarField realCount(mesh_.nCells(), 0.0);
        scalarField mass(mesh_.nCells(), 0.0);
        vectorField momentum(mesh_.nCells(), vector::zero);
        scalarField kinetic(mesh_.nCells(), 0.0);

        forAllIter(dsmcCloud, *this, iter)
        {
            const dsmcParcel& p = iter();

            if (!speciesSelected(fieldI, p.typeId()) || p.cell() < 0)
            {
                continue;
            }

            const scalar m = constProps(p.typeId()).mass();
            const scalar w = nParticle_;
            const label cellI = p.cell();

            parcelCount[cellI] += 1.0;
            realCount[cellI] += w;
            mass[cellI] += m*w;
            momentum[cellI] += m*w*p.U();
            kinetic[cellI] += 0.5*m*w*magSqr(p.U());
        }

        ++fieldAverageSteps_[fieldI];
        rhoNCumulative_[fieldI] += realCount;
        rhoMCumulative_[fieldI] += mass;
        momentumCumulative_[fieldI] += momentum;
        kineticCumulative_[fieldI] += kinetic;

        const scalar nAverageSteps = max(fieldAverageSteps_[fieldI], 1);

        volScalarField& dsmcN = dsmcNFields_[fieldI];
        volScalarField& rhoN = rhoNFields_[fieldI];
        volScalarField& rhoM = rhoMFields_[fieldI];
        volScalarField& p = pFields_[fieldI];
        volScalarField& Ttra = TtraFields_[fieldI];
        volVectorField& UMean = UMeanFields_[fieldI];
        volVectorField& fD = fDFields_[fieldI];

        dsmcN = dimensionedScalar("zero", dimless, 0.0);
        rhoN = dimensionedScalar("zero", dimDensity/dimMass, 0.0);
        rhoM = dimensionedScalar("zero", dimDensity, 0.0);
        p = dimensionedScalar("zero", dimPressure, 0.0);
        Ttra = dimensionedScalar("zero", dimTemperature, 0.0);
        UMean = dimensionedVector("zero", dimVelocity, vector::zero);
        fD = dimensionedVector("zero", dimPressure, vector::zero);

        forAll(parcelCount, cellI)
        {
            dsmcN[cellI] = parcelCount[cellI];

            const scalar averagedRealCount = rhoNCumulative_[fieldI][cellI]/nAverageSteps;
            const scalar averagedMass = rhoMCumulative_[fieldI][cellI]/nAverageSteps;

            if (averagedRealCount > SMALL)
            {
                const scalar vol = mesh_.cellVolumes()[cellI];
                const vector Uc = momentumCumulative_[fieldI][cellI]/max(rhoMCumulative_[fieldI][cellI], VSMALL);
                const scalar kineticMean = kineticCumulative_[fieldI][cellI]/nAverageSteps;
                const scalar rhoNValue = averagedRealCount/vol;
                const scalar rhoMValue = averagedMass/vol;
                const scalar eth = max(kineticMean/vol - 0.5*rhoMValue*magSqr(Uc), 0.0);
                const scalar T = (2.0/3.0)*eth/(physicoChemical::k.value()*rhoNValue);

                rhoN[cellI] = rhoNValue;
                rhoM[cellI] = rhoMValue;
                UMean[cellI] = Uc;
                Ttra[cellI] = T;
                p[cellI] = rhoNValue*physicoChemical::k.value()*T;
            }
        }

        dsmcN.correctBoundaryConditions();
        rhoN.correctBoundaryConditions();
        rhoM.correctBoundaryConditions();
        p.correctBoundaryConditions();
        Ttra.correctBoundaryConditions();
        UMean.correctBoundaryConditions();
        fD.correctBoundaryConditions();

        if
        (
            fieldSpecs_[fieldI].resetAtOutput
         && mesh_.time().writeTime()
         && mesh_.time().value() <= fieldSpecs_[fieldI].resetAtOutputUntilTime + SMALL
        )
        {
            fieldSpecs_[fieldI].pendingReset = true;
        }
    }
}


void Foam::dsmcCloud::initOpenMP()
{
    const dictionary& controlDict = mesh_.time().controlDict();

    openmpEnabled_ = controlDict.lookupOrDefault<bool>("useOpenMP", false);
    openmpMoveEnabled_ = controlDict.lookupOrDefault<bool>("openmpMove", openmpEnabled_);
    ompNumThreads_ = controlDict.lookupOrDefault<label>("openmpThreads", 0);
    openmpCollisionSchedule_ =
        controlDict.lookupOrDefault<word>("openmpCollisionSchedule", "dynamic");
    openmpCollisionChunk_ =
        controlDict.lookupOrDefault<label>("openmpCollisionChunk", 1);
    openmpMoveSchedule_ =
        controlDict.lookupOrDefault<word>("openmpMoveSchedule", "static");
    openmpMoveChunk_ =
        controlDict.lookupOrDefault<label>("openmpMoveChunk", 64);
    openmpMoveGuardLayers_ =
        controlDict.lookupOrDefault<label>("openmpMoveGuardLayers", 2);
    openmpAdaptivePartition_ =
        controlDict.lookupOrDefault<bool>("openmpAdaptivePartition", true);
    openmpPartitionMinUpdateInterval_ =
        controlDict.lookupOrDefault<label>("openmpPartitionMinUpdateInterval", 5);
    openmpPartitionImbalanceThreshold_ =
        controlDict.lookupOrDefault<scalar>("openmpPartitionImbalanceThreshold", 1.10);
    openmpPartitionDriftThreshold_ =
        controlDict.lookupOrDefault<scalar>("openmpPartitionDriftThreshold", 0.10);
    openmpCollisionCostCandidateWeight_ =
        controlDict.lookupOrDefault<scalar>("openmpCollisionCostCandidateWeight", 1.0);
    openmpCollisionCostActiveCellWeight_ =
        controlDict.lookupOrDefault<scalar>("openmpCollisionCostActiveCellWeight", 0.0);
    const bool legacyCollisionProfile =
        controlDict.lookupOrDefault<bool>("profileCollisionPhases", false);
    const bool legacyEvolveProfile =
        controlDict.lookupOrDefault<bool>("profileEvolvePhases", false);
    const bool summaryProfileEnabled =
        controlDict.lookupOrDefault<bool>
        (
            "profileSummary",
            legacyCollisionProfile || legacyEvolveProfile
        );
    const bool detailProfileEnabled =
        controlDict.lookupOrDefault<bool>("profileDetail", false);
    evolveStageProbe_ =
        controlDict.lookupOrDefault<bool>("evolveStageProbe", false);
    moveStageProbe_ =
        controlDict.lookupOrDefault<bool>("moveStageProbe", false);
    collisionProfileEnabled_ = summaryProfileEnabled || detailProfileEnabled;
    evolveProfileEnabled_ = summaryProfileEnabled || detailProfileEnabled;
    profilingDetailEnabled_ = detailProfileEnabled;

    #ifdef _OPENMP
    if (openmpEnabled_)
    {
        if
        (
            openmpCollisionSchedule_ != "dynamic"
         && openmpCollisionSchedule_ != "partition"
        )
        {
            WarningInFunction
                << "Unknown openmpCollisionSchedule '"
                << openmpCollisionSchedule_
                << "'. Falling back to 'dynamic'." << endl;

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
                << "'. Falling back to 'static'." << endl;

            openmpMoveSchedule_ = "static";
        }

        if (openmpMoveChunk_ < 1)
        {
            openmpMoveChunk_ = 1;
        }

        if (openmpMoveGuardLayers_ < 0)
        {
            openmpMoveGuardLayers_ = 0;
        }

        if (openmpCollisionChunk_ < 1)
        {
            openmpCollisionChunk_ = 1;
        }

        if (openmpPartitionMinUpdateInterval_ < 1)
        {
            openmpPartitionMinUpdateInterval_ = 1;
        }

        if (openmpPartitionImbalanceThreshold_ < 1.0)
        {
            openmpPartitionImbalanceThreshold_ = 1.0;
        }

        if (openmpPartitionDriftThreshold_ < 0.0)
        {
            openmpPartitionDriftThreshold_ = 0.0;
        }

        if (openmpCollisionCostCandidateWeight_ < 0.0)
        {
            openmpCollisionCostCandidateWeight_ = 0.0;
        }

        if (openmpCollisionCostActiveCellWeight_ < 0.0)
        {
            openmpCollisionCostActiveCellWeight_ = 0.0;
        }

        if
        (
            openmpCollisionCostCandidateWeight_ < SMALL
         && openmpCollisionCostActiveCellWeight_ < SMALL
        )
        {
            WarningInFunction
                << "Both collision partition cost weights are near zero. "
                << "Falling back to candidate-only weighting." << endl;
            openmpCollisionCostCandidateWeight_ = 1.0;
            openmpCollisionCostActiveCellWeight_ = 0.0;
        }

        if (ompNumThreads_ <= 0)
        {
            ompNumThreads_ = omp_get_max_threads();
        }

        if (ompNumThreads_ < 1)
        {
            ompNumThreads_ = 1;
        }

        ompRndGens_.setSize(ompNumThreads_);
        particleLoadStart_.setSize(ompNumThreads_, 0);
        particleLoadEnd_.setSize(ompNumThreads_, mesh_.nCells());
        collisionLoadStart_.setSize(ompNumThreads_, 0);
        collisionLoadEnd_.setSize(ompNumThreads_, mesh_.nCells());
        particlePartitionStep_ = 0;
        particlePartitionLastRebuildStep_ = -1;
        particlePartitionReferenceLoad_.clear();
        collisionPartitionStep_ = 0;
        collisionPartitionLastRebuildStep_ = -1;
        collisionPartitionReferenceLoad_.clear();

        forAll(ompRndGens_, threadI)
        {
            const label seed = 104729*Pstream::myProcNo() + threadI + 1;
            ompRndGens_[threadI].reset(seed);
        }

        Info<< "OpenMP enabled for dsmcCloud with "
            << ompNumThreads_ << " thread-local RNG streams"
            << " using collision schedule '" << openmpCollisionSchedule_
            << "' (chunk " << openmpCollisionChunk_ << ")"
            << "', move kernel "
            << (openmpMoveEnabled_ ? "enabled" : "disabled")
            << " (" << openmpMoveSchedule_ << ", chunk "
            << openmpMoveChunk_ << ")"
            << ", adaptive partition "
            << (openmpAdaptivePartition_ ? "enabled" : "disabled")
            << " (interval " << openmpPartitionMinUpdateInterval_
            << ", imbalance " << openmpPartitionImbalanceThreshold_
            << ", drift " << openmpPartitionDriftThreshold_ << ")"
            << ", collision partition cost (candidate x"
            << openmpCollisionCostCandidateWeight_
            << ", active-cell x" << openmpCollisionCostActiveCellWeight_
            << ")"
            << endl;
    }
    else
    {
        openmpMoveEnabled_ = false;
        ompNumThreads_ = 1;
        ompRndGens_.clear();
        particleLoadStart_.setSize(1, 0);
        particleLoadEnd_.setSize(1, mesh_.nCells());
        collisionLoadStart_.setSize(1, 0);
        collisionLoadEnd_.setSize(1, mesh_.nCells());
        particlePartitionStep_ = 0;
        particlePartitionLastRebuildStep_ = -1;
        particlePartitionReferenceLoad_.clear();
        collisionPartitionStep_ = 0;
        collisionPartitionLastRebuildStep_ = -1;
        collisionPartitionReferenceLoad_.clear();
    }
    #else
    if (openmpEnabled_)
    {
        WarningInFunction
            << "OpenMP requested via controlDict entry 'useOpenMP', "
            << "but the code was built without OpenMP support. "
            << "Falling back to serial execution." << endl;
    }

    openmpEnabled_ = false;
    openmpMoveEnabled_ = false;
    ompNumThreads_ = 1;
    ompRndGens_.clear();
    particleLoadStart_.setSize(1, 0);
    particleLoadEnd_.setSize(1, mesh_.nCells());
    collisionLoadStart_.setSize(1, 0);
    collisionLoadEnd_.setSize(1, mesh_.nCells());
    particlePartitionStep_ = 0;
    particlePartitionLastRebuildStep_ = -1;
    particlePartitionReferenceLoad_.clear();
    collisionPartitionStep_ = 0;
    collisionPartitionLastRebuildStep_ = -1;
    collisionPartitionReferenceLoad_.clear();
    #endif
}


void Foam::dsmcCloud::initOpenMPMoveGuardCells()
{
    openmpMoveGuardCells_.setSize(mesh_.nCells(), false);
    openmpMoveGuardCellCount_ = 0;

    const bool replicated = replicatedMesh_.valid() && replicatedMesh_->active();
    if (!(openmpEnabled_ && openmpMoveEnabled_ && Pstream::parRun()))
    {
        return;
    }

    DynamicList<label> frontierCells;
    frontierCells.setCapacity(min(mesh_.nCells(), mesh_.boundaryMesh().size()*8));

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        // Standard parallel: guard processor-patch cells for transfer.
        // Replicated mesh or pure OpenMP: guard physical-boundary cells.
        if (Pstream::parRun() && !replicated)
        {
            if (!isA<processorPolyPatch>(pp)) continue;
        }
        else
        {
            if (isA<processorPolyPatch>(pp)) continue;
        }

        const labelUList& faceCells = pp.faceCells();

        forAll(faceCells, i)
        {
            const label celli = faceCells[i];

            if (!openmpMoveGuardCells_[celli])
            {
                openmpMoveGuardCells_[celli] = true;
                frontierCells.append(celli);
                ++openmpMoveGuardCellCount_;
            }
        }
    }

    for (label layer = 0; layer < openmpMoveGuardLayers_; ++layer)
    {
        DynamicList<label> nextFrontier;
        nextFrontier.setCapacity(frontierCells.size()*2 + 1);

        forAll(frontierCells, frontierI)
        {
            const label celli = frontierCells[frontierI];
            const cell& cFaces = mesh_.cells()[celli];

            forAll(cFaces, faceI)
            {
                const label meshFaceI = cFaces[faceI];

                if (!mesh_.isInternalFace(meshFaceI))
                {
                    continue;
                }

                const label owner = mesh_.faceOwner()[meshFaceI];
                const label neighbour = mesh_.faceNeighbour()[meshFaceI];
                const label otherCelli = owner == celli ? neighbour : owner;

                if (otherCelli >= 0 && !openmpMoveGuardCells_[otherCelli])
                {
                    openmpMoveGuardCells_[otherCelli] = true;
                    nextFrontier.append(otherCelli);
                    ++openmpMoveGuardCellCount_;
                }
            }
        }

        frontierCells.transfer(nextFrontier);

        if (!frontierCells.size())
        {
            break;
        }
    }

    Info<< "OpenMP mixed-move guard cells enabled on "
        << openmpMoveGuardCellCount_ << " / " << mesh_.nCells()
        << " cells (processor halo layers " << openmpMoveGuardLayers_
        << ')' << endl;
}


void Foam::dsmcCloud::precomputeCollisionCandidates()
{
    if (selectedPairsPerCell_.size() != mesh_.nCells())
    {
        selectedPairsPerCell_.setSize(mesh_.nCells(), 0.0);
    }

    if (nCandidatesPerCell_.size() != mesh_.nCells())
    {
        nCandidatesPerCell_.setSize(mesh_.nCells(), 0);
    }

    forAll(collisionCandidateCells_, candidateI)
    {
        const label celli = collisionCandidateCells_[candidateI];
        selectedPairsPerCell_[celli] = 0.0;
        nCandidatesPerCell_[celli] = 0;
    }

    collisionCandidateCells_.clear();

    const labelList& collisionCells = occupancyCollisionCells_;

    #ifdef _OPENMP
    if (openmpEnabled_)
    {
        labelList candidateFlags(collisionCells.size(), 0);

        #pragma omp parallel for schedule(static)
        forAll(collisionCells, collisionCellI)
        {
            const label celli = collisionCells[collisionCellI];
            const label nC = occupancyCount(celli);
            const scalar sigmaTcRMaxCell = max(sigmaTcRMax_[celli], SMALL);
            const scalar selectedPairs =
                collisionSelectionRemainder_[celli]
              + 0.5*nC*(nC - 1)*nParticles(celli)*sigmaTcRMaxCell*deltaTValue(celli)
               /mesh_.cellVolumes()[celli];

            selectedPairsPerCell_[celli] = selectedPairs;
            nCandidatesPerCell_[celli] = label(selectedPairs);
            collisionSelectionRemainder_[celli] =
                selectedPairs - scalar(nCandidatesPerCell_[celli]);

            if (nCandidatesPerCell_[celli] > 0)
            {
                candidateFlags[collisionCellI] = 1;
            }
        }

        label candidateCount = 0;

        forAll(candidateFlags, collisionCellI)
        {
            candidateCount += candidateFlags[collisionCellI];
        }

        collisionCandidateCells_.setSize(candidateCount);
        candidateCount = 0;

        forAll(candidateFlags, collisionCellI)
        {
            if (candidateFlags[collisionCellI])
            {
                collisionCandidateCells_[candidateCount++] =
                    collisionCells[collisionCellI];
            }
        }
    }
    else
    #endif
    {
        collisionCandidateCells_.setSize(collisionCells.size());
        label candidateCount = 0;

        forAll(collisionCells, collisionCellI)
        {
            const label celli = collisionCells[collisionCellI];
            const label nC = occupancyCount(celli);
            const scalar sigmaTcRMaxCell = max(sigmaTcRMax_[celli], SMALL);
            const scalar selectedPairs =
                collisionSelectionRemainder_[celli]
              + 0.5*nC*(nC - 1)*nParticles(celli)*sigmaTcRMaxCell*deltaTValue(celli)
               /mesh_.cellVolumes()[celli];

            selectedPairsPerCell_[celli] = selectedPairs;
            nCandidatesPerCell_[celli] = label(selectedPairs);
            collisionSelectionRemainder_[celli] =
                selectedPairs - scalar(nCandidatesPerCell_[celli]);

            if (nCandidatesPerCell_[celli] > 0)
            {
                collisionCandidateCells_[candidateCount++] = celli;
            }
        }

        collisionCandidateCells_.setSize(candidateCount);
    }
}


void Foam::dsmcCloud::rebuildCollisionLoadPartition()
{
    if (!openmpEnabled_ || ompNumThreads_ <= 1)
    {
        collisionLoadStart_.setSize(1, 0);
        collisionLoadEnd_.setSize(1, mesh_.nCells());
        return;
    }

    ++collisionPartitionStep_;

    const label nCells = mesh_.nCells();
    scalarField currentCellLoads(nCells, 0.0);
    scalar totalLoad = 0.0;

    forAll(collisionCandidateCells_, candidateI)
    {
        const label celli = collisionCandidateCells_[candidateI];
        const label candidateCount = nCandidatesPerCell_[celli];
        const scalar activeCell = candidateCount > 0 ? 1.0 : 0.0;
        const scalar load =
            openmpCollisionCostCandidateWeight_*scalar(candidateCount)
          + openmpCollisionCostActiveCellWeight_*activeCell;

        currentCellLoads[celli] = load;
        totalLoad += load;
    }

    if
    (
        openmpAdaptivePartition_
     && collisionPartitionLastRebuildStep_ >= 0
     && collisionLoadStart_.size() == ompNumThreads_
     && collisionLoadEnd_.size() == ompNumThreads_
     && collisionPartitionReferenceLoad_.size() == nCells
    )
    {
        scalar maxThreadLoad = 0.0;

        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            const label start =
                max(min(collisionLoadStart_[threadI], nCells), label(0));
            const label end =
                max(min(collisionLoadEnd_[threadI], nCells), start);

            scalar threadLoad = 0.0;

            for (label celli = start; celli < end; ++celli)
            {
                threadLoad += currentCellLoads[celli];
            }

            maxThreadLoad = max(maxThreadLoad, threadLoad);
        }

        const scalar avgThreadLoad =
            totalLoad/max(scalar(ompNumThreads_), scalar(1));
        const scalar imbalance =
            avgThreadLoad > SMALL ? maxThreadLoad/avgThreadLoad : 1.0;

        scalar driftAccum = 0.0;
        forAll(currentCellLoads, celli)
        {
            driftAccum += mag
            (
                scalar(currentCellLoads[celli] - collisionPartitionReferenceLoad_[celli])
            );
        }
        const scalar drift =
            totalLoad > SMALL ? driftAccum/totalLoad : 0.0;

        const label sinceLastRebuild =
            collisionPartitionStep_ - collisionPartitionLastRebuildStep_;

        if
        (
            sinceLastRebuild < openmpPartitionMinUpdateInterval_
         && imbalance <= openmpPartitionImbalanceThreshold_
         && drift <= openmpPartitionDriftThreshold_
        )
        {
            return;
        }

        if
        (
            imbalance <= openmpPartitionImbalanceThreshold_
         && drift <= openmpPartitionDriftThreshold_
        )
        {
            return;
        }
    }

    collisionLoadStart_.setSize(ompNumThreads_, mesh_.nCells());
    collisionLoadEnd_.setSize(ompNumThreads_, mesh_.nCells());

    if (totalLoad <= SMALL)
    {
        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            collisionLoadStart_[threadI] = threadI*mesh_.nCells()/ompNumThreads_;
            collisionLoadEnd_[threadI] = (threadI + 1)*mesh_.nCells()/ompNumThreads_;
        }

        collisionPartitionReferenceLoad_.transfer(currentCellLoads);
        collisionPartitionLastRebuildStep_ = collisionPartitionStep_;
        return;
    }

    const scalar avgLoad = totalLoad/scalar(ompNumThreads_);
    scalar accumulatedLoad = 0.0;
    label start = 0;

    for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
    {
        collisionLoadStart_[threadI] = start;

        if (threadI == ompNumThreads_ - 1)
        {
            collisionLoadEnd_[threadI] = mesh_.nCells();
            break;
        }

        const scalar targetLoad = avgLoad*scalar(threadI + 1);
        label end = start;

        while (end < mesh_.nCells() && accumulatedLoad < targetLoad)
        {
            accumulatedLoad += currentCellLoads[end];
            ++end;
        }

        collisionLoadEnd_[threadI] = end;
        start = end;
    }

    collisionPartitionReferenceLoad_.transfer(currentCellLoads);
    collisionPartitionLastRebuildStep_ = collisionPartitionStep_;
}


void Foam::dsmcCloud::rebuildParticleLoadPartition()
{
    if (!openmpEnabled_ || ompNumThreads_ <= 1)
    {
        particleLoadStart_.setSize(1, 0);
        particleLoadEnd_.setSize(1, mesh_.nCells());
        return;
    }

    // Replicated mesh: only iterate owned cells to avoid scanning
    // 60k cells when only ~7.5k have particles.
    if (replicatedMesh_.valid() && replicatedMesh_->active())
    {
        const auto& myCells = replicatedMesh_->myCells();
        const label nOwned = myCells.size();

        if (nOwned == 0)
        {
            particleLoadStart_.setSize(1, 0);
            particleLoadEnd_.setSize(1, 0);
            return;
        }

        // myCells_ is sorted by construction; divide across OMP threads
        particleLoadStart_.setSize(ompNumThreads_, 0);
        particleLoadEnd_.setSize(ompNumThreads_, nOwned);

        const label perThread = nOwned / ompNumThreads_;
        const label remainder = nOwned % ompNumThreads_;
        label start = 0;

        for (label t = 0; t < ompNumThreads_; ++t)
        {
            const label count = perThread + (t < remainder ? 1 : 0);
            label end = min(start + count, nOwned);

            // Convert indices-in-myCells_ to actual cell indices.
            // particleLoadEnd_ is exclusive (one past last cell).
            particleLoadStart_[t] =
                (start < nOwned) ? myCells[start] : myCells[nOwned-1];
            particleLoadEnd_[t] =
                (end < nOwned) ? myCells[end] : (myCells[nOwned-1] + 1);

            start = end;
        }
        return;
    }

    const label nCells = mesh_.nCells();
    labelList currentCellLoads(nCells, 0);
    label totalParticles = 0;

    for (label celli = 0; celli < nCells; ++celli)
    {
        const label load = occupancyCount(celli);
        currentCellLoads[celli] = load;
        totalParticles += load;
    }

    rebuildParticleLoadPartition(currentCellLoads, totalParticles);
}


void Foam::dsmcCloud::rebuildParticleLoadPartition
(
    const labelList& currentCellLoads,
    const label totalParticles
)
{
    if (!openmpEnabled_ || ompNumThreads_ <= 1)
    {
        particleLoadStart_.setSize(1, 0);
        particleLoadEnd_.setSize(1, mesh_.nCells());
        return;
    }

    ++particlePartitionStep_;

    const label nCells = mesh_.nCells();

    if
    (
        openmpAdaptivePartition_
     && particlePartitionLastRebuildStep_ >= 0
     && particleLoadStart_.size() == ompNumThreads_
     && particleLoadEnd_.size() == ompNumThreads_
     && particlePartitionReferenceLoad_.size() == nCells
    )
    {
        label maxThreadLoad = 0;

        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            const label start = max(min(particleLoadStart_[threadI], nCells), label(0));
            const label end = max(min(particleLoadEnd_[threadI], nCells), start);

            label threadLoad = 0;
            for (label celli = start; celli < end; ++celli)
            {
                threadLoad += currentCellLoads[celli];
            }

            maxThreadLoad = max(maxThreadLoad, threadLoad);
        }

        const scalar avgThreadLoad =
            scalar(totalParticles)/max(scalar(ompNumThreads_), scalar(1));
        const scalar imbalance =
            avgThreadLoad > SMALL ? scalar(maxThreadLoad)/avgThreadLoad : 1.0;

        scalar driftAccum = 0.0;
        forAll(currentCellLoads, celli)
        {
            driftAccum += mag
            (
                scalar(currentCellLoads[celli] - particlePartitionReferenceLoad_[celli])
            );
        }
        const scalar drift =
            totalParticles > 0 ? driftAccum/scalar(totalParticles) : 0.0;

        const label sinceLastRebuild =
            particlePartitionStep_ - particlePartitionLastRebuildStep_;

        if
        (
            sinceLastRebuild < openmpPartitionMinUpdateInterval_
         && imbalance <= openmpPartitionImbalanceThreshold_
         && drift <= openmpPartitionDriftThreshold_
        )
        {
            return;
        }

        if
        (
            imbalance <= openmpPartitionImbalanceThreshold_
         && drift <= openmpPartitionDriftThreshold_
        )
        {
            return;
        }
    }

    particleLoadStart_.setSize(ompNumThreads_, mesh_.nCells());
    particleLoadEnd_.setSize(ompNumThreads_, mesh_.nCells());

    if (totalParticles <= 0)
    {
        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            particleLoadStart_[threadI] = threadI*mesh_.nCells()/ompNumThreads_;
            particleLoadEnd_[threadI] = (threadI + 1)*mesh_.nCells()/ompNumThreads_;
        }

        particlePartitionReferenceLoad_ = currentCellLoads;
        particlePartitionLastRebuildStep_ = particlePartitionStep_;
        return;
    }

    scalarField threadWeights(ompNumThreads_, 1.0);
    scalar totalWeight = scalar(ompNumThreads_);

    if
    (
        moveLastThreadParticleCounts_.size() == ompNumThreads_
     && moveLastThreadWallTimes_.size() == ompNumThreads_
    )
    {
        totalWeight = 0.0;

        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            const label processed = moveLastThreadParticleCounts_[threadI];
            const scalar wallTime = moveLastThreadWallTimes_[threadI];

            if (processed > 0 && wallTime > SMALL)
            {
                // Use last-step throughput so faster threads receive more cells.
                threadWeights[threadI] = scalar(processed)/wallTime;
            }
            else
            {
                threadWeights[threadI] = 1.0;
            }

            totalWeight += threadWeights[threadI];
        }

        if (totalWeight <= SMALL)
        {
            threadWeights = scalarField(ompNumThreads_, 1.0);
            totalWeight = scalar(ompNumThreads_);
        }
    }

    scalar accumulatedParticles = 0.0;
    scalar accumulatedWeight = 0.0;
    label start = 0;

    for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
    {
        particleLoadStart_[threadI] = start;

        if (threadI == ompNumThreads_ - 1)
        {
            particleLoadEnd_[threadI] = mesh_.nCells();
            break;
        }

        accumulatedWeight += threadWeights[threadI];
        const scalar targetParticles =
            scalar(totalParticles)*(accumulatedWeight/totalWeight);
        label end = start;

        while (end < mesh_.nCells() && accumulatedParticles < targetParticles)
        {
            accumulatedParticles += currentCellLoads[end];
            ++end;
        }

        particleLoadEnd_[threadI] = end;
        start = end;
    }

    particlePartitionReferenceLoad_ = currentCellLoads;
    particlePartitionLastRebuildStep_ = particlePartitionStep_;
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


void Foam::dsmcCloud::clearMoveOrderedParcels()
{
    moveOrderedParcels_.clear();
    moveOrderedThreadOffsets_.clear();
    moveOrderedParcelsValid_ = false;
}


void Foam::dsmcCloud::rebuildMoveOrderedParcels()
{
    moveOrderedParcels_.setSize(this->size());
    label i = 0;
    forAllIter(Cloud<dsmcParcel>, *this, iter)
    {
        moveOrderedParcels_[i++] = &iter();
    }
    moveOrderedParcels_.setSize(i);

    const label nThreads = ompNumThreads_;
    moveOrderedThreadOffsets_.setSize(nThreads + 1);
    for (label t = 0; t <= nThreads; ++t)
    {
        moveOrderedThreadOffsets_[t] = t * i / nThreads;
    }
    moveOrderedParcelsValid_ = true;
    moveOrderedReuseDisabled_ = false;
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::setMoveOrderedParcels(const DynamicList<dsmcParcel*>& parcels)
{
    moveOrderedParcels_.setSize(parcels.size());
    forAll(parcels, i) moveOrderedParcels_[i] = parcels[i];

    const label nThreads = ompNumThreads_;
    moveOrderedThreadOffsets_.setSize(nThreads + 1);
    for (label t = 0; t <= nThreads; ++t)
    {
        moveOrderedThreadOffsets_[t] = t * parcels.size() / nThreads;
    }
    moveOrderedParcelsValid_ = true;
    moveOrderedReuseDisabled_ = false;
    moveAppendedParcels_.clear();
}


void Foam::dsmcCloud::appendToMoveOrdered(dsmcParcel* p)
{
    const label oldSize = moveOrderedParcels_.size();
    moveOrderedParcels_.setSize(oldSize + 1);
    moveOrderedParcels_[oldSize] = p;
}


void Foam::dsmcCloud::appendBatchToMoveOrdered(const DynamicList<dsmcParcel*>& parcels)
{
    const label oldSize = moveOrderedParcels_.size();
    moveOrderedParcels_.setSize(oldSize + parcels.size());
    forAll(parcels, i)
    {
        moveOrderedParcels_[oldSize + i] = parcels[i];
    }
}


void Foam::dsmcCloud::refreshTrackerUsage()
{
    trackerActive_ = false;

    const auto& configuredFields = fields_.fields();

    forAll(configuredFields, i)
    {
        if (configuredFields[i].valid())
        {
            const word fieldType(configuredFields[i]->type());

            if (fieldType == "dsmcFluxSurface")
            {
                trackerActive_ = true;
                break;
            }
        }
    }
}


Foam::dsmcCloud::dsmcCloud
(
    Time&,
    const word& cloudName,
    const fvMesh& mesh,
    bool readFields
)
:
    Cloud<dsmcParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_(IOobject(cloudName + "Properties", mesh_.time().constant(), mesh_, IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE)),
    typeIdList_(particleProperties_.get<wordList>("typeIdList")),
    dsmcCoordinateSystem_(dsmcCoordinateSystem::New(const_cast<Time&>(mesh_.time()), mesh_, *this)),
    nParticle_(particleProperties_.lookupOrDefault<scalar>("nEquivalentParticles", 1.0)),
    nTerminalOutputs_(mesh_.time().controlDict().lookupOrDefault<label>("nTerminalOutputs", 1)),
    collisionPartnerSelectionModel_(particleProperties_.lookupOrDefault<word>("collisionPartnerSelectionModel", "noTimeCounter")),
    cellOccupancy_(mesh_.nCells()),
    sigmaTcRMax_(IOobject(this->name() + "SigmaTcRMax", mesh_.time().timeName(), mesh_, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh_, dimensionedScalar("sigmaTcRMax", dimensionSet(0, 3, -1, 0, 0), SMALL)),
    collisionSelectionRemainder_(mesh_.nCells(), 0.0),
    constProps_(),
    rndGen_(Pstream::myProcNo()),
    openmpEnabled_(false),
    openmpMoveEnabled_(false),
    trackerActive_(true),
    ompNumThreads_(1),
    openmpCollisionSchedule_("dynamic"),
    openmpCollisionChunk_(1),
    openmpMoveSchedule_("static"),
    openmpMoveChunk_(64),
    openmpMoveGuardLayers_(2),
    openmpMoveGuardCells_(),
    openmpMoveGuardCellCount_(0),
    openmpAdaptivePartition_(true),
    openmpPartitionMinUpdateInterval_(5),
    openmpPartitionImbalanceThreshold_(1.10),
    openmpPartitionDriftThreshold_(0.10),
    openmpCollisionCostCandidateWeight_(1.0),
    openmpCollisionCostActiveCellWeight_(0.0),
    collisionProfileEnabled_(false),
    evolveProfileEnabled_(false),
    profilingDetailEnabled_(false),
    emitStepDiagnostics_(false),
    evolveStageProbe_(false),
    moveStageProbe_(false),
    ompRndGens_(),
    particleLoadStart_(),
    particleLoadEnd_(),
    particlePartitionStep_(0),
    particlePartitionLastRebuildStep_(-1),
    particlePartitionReferenceLoad_(),
    moveOrderedParcels_(),
    moveOrderedThreadOffsets_(),
    moveOrderedParcelsValid_(false),
    moveOrderedReuseDisabled_(false),
    moveAppendCaptureActive_(false),
    moveAppendToPending_(false),
    moveAppendedParcels_(),
    pendingMoveParcels_(),
    occupancyOrderedParcels_(),
    occupancyCellOffsets_(),
    occupancyActiveCells_(),
    occupancyCollisionCells_(),
    occupancyOrderedParcelsValid_(false),
    cellOccupancyMaterialized_(true),
    occupancyThreadCellCounts_(),
    occupancyThreadActiveCells_(),
    selectedPairsPerCell_(mesh_.nCells(), 0.0),
    nCandidatesPerCell_(mesh_.nCells(), 0),
    moveItersPerCell_(mesh_.nCells(), 0),
    moveItersPerCellCumulative_(mesh_.nCells(), 0),
    collisionCandidateCells_(),
    collisionPartitionStep_(0),
    collisionPartitionLastRebuildStep_(-1),
    collisionPartitionReferenceLoad_(),
    collisionPrecomputeWallTime_(0.0),
    collisionPartitionWallTime_(0.0),
    collisionSelectionWallTime_(0.0),
    collisionProfileCalls_(0),
    cumulativeCollisions_(0),
    cumulativeCollisionCandidates_(0),
    collisionPhaseActive_(false),
    collisionNewParcels_(),
    buildOccupancyExtractWallTime_(0.0),
    buildOccupancyCountWallTime_(0.0),
    buildOccupancyAssembleWallTime_(0.0),
    buildOccupancyProfileCalls_(0),
    buildOccupancyMoveOrderedHits_(0),
    buildOccupancyFallbackHits_(0),
    buildOccupancyMoveOrderedInvalidHits_(0),
    buildOccupancyMoveOrderedOffsetMismatchHits_(0),
    buildOccupancyMoveOrderedSizeMismatchHits_(0),
    buildOccupancyMoveOrderedSizeDeltaSum_(0.0),
    buildOccupancyMoveOrderedSizeDeltaMin_(labelMax),
    buildOccupancyMoveOrderedSizeDeltaMax_(labelMin),
    buildOccupancyMoveOrderedParcelsSum_(0.0),
    buildOccupancyMoveAppendedParcelsSum_(0.0),
    buildOccupancyCloudSizeSum_(0.0),
    buildOccupancyCloudValidCellSum_(0.0),
    buildOccupancyCloudInvalidCellSum_(0.0),
    buildOccupancyMoveOrderedValidCellSum_(0.0),
    buildOccupancyMoveOrderedInvalidCellSum_(0.0),
    buildOccupancyPartitionAssignedParcelsSum_(0.0),
    buildOccupancyPartitionGapCellsSum_(0.0),
    buildOccupancyPartitionOverlapCellsSum_(0.0),
    movePreControlWallTime_(0.0),
    moveResetSetupWallTime_(0.0),
    moveExtractWallTime_(0.0),
    moveExtractDeferredWallTime_(0.0),
    moveExtractOrderedReuseWallTime_(0.0),
    moveExtractFullScanWallTime_(0.0),
    moveKernelWallTime_(0.0),
    moveCommitWallTime_(0.0),
    moveTransferFinalizeWallTime_(0.0),
    moveProfileCalls_(0),
    moveCommitCountCalls_(0),
    moveExtractDeferredPasses_(0),
    moveExtractOrderedReusePasses_(0),
    moveExtractFullScanPasses_(0),
    moveExtractedParcelsSum_(0.0),
    moveExtractDeferredParcelsSum_(0.0),
    moveExtractOrderedReuseParcelsSum_(0.0),
    moveExtractFullScanParcelsSum_(0.0),
    moveFirstPassReuseChecks_(0),
    moveFirstPassHasOrderedCount_(0),
    moveFirstPassMissingOrderedCount_(0),
    moveFirstPassOffsetMismatchCount_(0),
    moveFirstPassSizeMatchCount_(0),
    moveFirstPassPriorPlusPendingMatchCount_(0),
    moveFirstPassCurrentSizeSum_(0.0),
    moveFirstPassPriorSizeSum_(0.0),
    moveFirstPassPendingSizeSum_(0.0),
    moveFirstPassSizeDeltaSum_(0.0),
    moveFirstPassDeltaMinusPendingSum_(0.0),
    previousMoveEndCloudSize_(-1),
    moveInterStepCloudDeltaSum_(0.0),
    moveInterStepCloudDeltaMin_(labelMax),
    moveInterStepCloudDeltaMax_(labelMin),
    movePostStepCloudDeltaSum_(0.0),
    movePostStepCloudDeltaMin_(labelMax),
    movePostStepCloudDeltaMax_(labelMin),
    moveSurvivorParcelsSum_(0.0),
    moveTransferredParcelsSum_(0.0),
    moveDeletedParcelsSum_(0.0),
    moveLoopPassesSum_(0.0),
    moveLoopPassesMax_(0),
    moveDeferredParcelsSum_(0.0),
    moveDeferredParcelsMax_(0),
    moveReceivedParcelsSum_(0.0),
    moveReceivedParcelsMax_(0),
    moveTrackWallTime_(0.0),
    moveTrackerCallbackWallTime_(0.0),
    moveBoundaryControlWallTime_(0.0),
    moveFaceHitsSum_(0.0),
    moveCyclicHitsSum_(0.0),
    moveStuckHitsSum_(0.0),
    movePatchHitsSum_(0.0),
    moveProcessorHitsSum_(0.0),
    evolveMoveWallTime_(0.0),
    evolveBuildWallTime_(0.0),
    evolveCoordWallTime_(0.0),
    evolveCollisionWallTime_(0.0),
    evolveReactionWallTime_(0.0),
    evolvePostWallTime_(0.0),
    evolvePreWallTime_(0.0),
    evolveRankTimeWallTime_(0.0),
    evolveAutoRebalanceWallTime_(0.0),
    evolvePostProfileResidualWallTime_(0.0),
    evolveFullWallTime_(0.0),
    evolveProfileCalls_(0),
    moveThreadParticleCounts_(),
    moveThreadWallTimes_(),
    moveThreadTrackWallTimes_(),
    moveThreadTrackerWallTimes_(),
    moveThreadBoundaryWallTimes_(),
    moveThreadFaceHitCounts_(),
    moveThreadCyclicHitCounts_(),
    moveThreadStuckHitCounts_(),
    moveThreadPatchHitCounts_(),
    moveThreadProcessorHitCounts_(),
    moveLastThreadParticleCounts_(),
    moveLastThreadWallTimes_(),
    collisionThreadCandidateCounts_(),
    collisionThreadAcceptedCounts_(),
    collisionThreadActiveCellCounts_(),
    collisionThreadReactionHitCounts_(),
    collisionThreadWallTimes_(),
    porousMeasurements_(porousMeasurements::New(const_cast<Time&>(mesh_.time()), mesh_, *this)),
    controllers_(const_cast<Time&>(mesh_.time()), mesh_, *this),
    boundaryMeas_(mesh, *this, true),
    fields_(const_cast<Time&>(mesh_.time()), mesh_, *this),
    boundaries_(const_cast<Time&>(mesh_.time()), mesh_, *this),
    trackingInfo_(mesh, *this, true),
    binaryCollisionModel_(BinaryCollisionModel::New(particleProperties_, *this)),
    collisionPartnerSelectionPtr_(nullptr),
    reactions_(nullptr),
    cellMeas_(mesh, *this, true),
    diffuseWalls_(),
    wallPatchToModelId_(),
    inflows_(),
    fieldSpecs_(),
    dsmcNFields_(),
    rhoNFields_(),
    rhoMFields_(),
    pFields_(),
    TtraFields_(),
    UMeanFields_(),
    fDFields_(),
    speciesBoundaryForceDensity_(),
    fieldAverageSteps_(),
    rhoNCumulative_(),
    rhoMCumulative_(),
    momentumCumulative_(),
    kineticCumulative_()
{
    if (!readFields)
    {
        // Avoid reusing stale lagrangian positions when the caller requests
        // a fresh initialisation-only cloud.
        this->clear();
    }

    initOpenMP();
    initOpenMPMoveGuardCells();

    coordSystem().checkCoordinateSystemInputs();
    buildConstProps();

    // Phase A: replicated mesh DLB (all ranks hold full mesh, cellOwner_ routing)
    const bool replicatedMesh =
        mesh_.time().controlDict().lookupOrDefault<bool>("replicatedMesh", false);
    if (replicatedMesh && readFields)
    {
        replicatedMesh_.reset(new dsmcReplicatedMesh(*this, mesh_));
        replicatedMesh_->initialize();

        // In replicated-mesh mode (no -parallel), all ranks read the full
        // particle set. Only rank 0 keeps them; others clear and receive
        // their share via the first migrateParticlesByCellOwner() call.
        // In decomposed mode (with -parallel), each rank already has its
        // own subset — no clearing needed.
        if (!Pstream::parRun())
        {
            int mpiInit = 0;
            MPI_Initialized(&mpiInit);
            if (!mpiInit) UPstream::initNull();
            int myRank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
            if (myRank != 0)
            {
                this->clear();
            }
            Info<< "Replicated mesh: rank " << myRank
                << " starts with " << this->size() << " parcels" << endl;
        }
        else
        {
            Info<< "Replicated mesh: rank " << Pstream::myProcNo()
                << " starts with " << this->size() << " parcels"
                << " (decomposed mesh)" << endl;
        }
    }

    IOobject chemReactIO
    (
        "chemReactDict",
        mesh_.time().system(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (chemReactIO.typeHeaderOk<IOdictionary>(true))
    {
        reactions_.reset(new dsmcReactions(mesh_.time(), mesh_, *this));
        reactions_->initialConfiguration();
    }

    collisionPartnerSelectionPtr_ = collisionPartnerSelection::New
    (
        mesh_,
        *this,
        particleProperties_
    );
    collisionPartnerSelectionPtr_->initialConfiguration();

    fields_.createFields();
    refreshTrackerUsage();
    boundaryMeas_.setInitialConfig();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();

    if (max(sigmaTcRMax_.primitiveField()) <= SMALL)
    {
        sigmaTcRMax_.primitiveFieldRef() = initialSigmaTcRMax();
    }

    sigmaTcRMax_.correctBoundaryConditions();

    forAll(collisionSelectionRemainder_, i)
    {
        collisionSelectionRemainder_[i] = rndGen_.sample01<scalar>();
    }

    if (readFields)
    {
        dsmcParcel::readFields(*this);
    }

    buildCellOccupancy(readFields);
}


Foam::dsmcCloud::~dsmcCloud()
{}


void Foam::dsmcCloud::addNewParcel
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
)
{
    dsmcParcel* pPtr = new dsmcParcel(mesh_, coordinates, cellI, tetFaceI, tetPtI, U, RWF, ERot, ELevel, typeId, newParcel, classification, vibLevel);

    if (collisionPhaseActive_)
    {
        const label threadI = currentThreadId();
        collisionNewParcels_[threadI].append(pPtr);
        return;
    }

    Cloud<dsmcParcel>::addParticle(pPtr);
    recordMoveAppendedParcel(pPtr);
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
    dsmcParcel* pPtr = new dsmcParcel(mesh_, position, cellI, U, RWF, ERot, ELevel, typeId, newParcel, classification, vibLevel);

    if (collisionPhaseActive_)
    {
        const label threadI = currentThreadId();
        collisionNewParcels_[threadI].append(pPtr);
        return;
    }

    Cloud<dsmcParcel>::addParticle(pPtr);
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
    dsmcParcel* pPtr = new dsmcParcel(mesh_, position, cellI, U, RWF, ERot, ELevel, typeId, newParcel, classification, vibLevel);
    pPtr->setStuck(wallTemperature, wallVectors);
    Cloud<dsmcParcel>::addParticle(pPtr);
    recordMoveAppendedParcel(pPtr);
}


void Foam::dsmcCloud::handleWallInteraction(dsmcParcel& p, dsmcParcel::trackingData& td)
{
    const label patchId = p.patch();

    if (specularWallPatchIds_.found(patchId))
    {
        const label faceI = p.face();
        const vector sF = mesh_.faceAreas()[faceI];
        const vector n = -sF/max(mag(sF), VSMALL);
        p.U() -= 2.0*(p.U() & n)*n;
        td.keepParticle = true;
        return;
    }

    if (patchId < 0 || patchId >= wallPatchToModelId_.size() || wallPatchToModelId_[patchId] < 0)
    {
        td.keepParticle = false;
        return;
    }

    const diffuseWallSpec& wall = diffuseWalls_[wallPatchToModelId_[patchId]];
    const label faceI = p.face();
    const vector sF = mesh_.faceAreas()[faceI];
    const scalar fA = mag(sF);
    const vector n = -sF/fA;
    const vector t1 = makeTangential(mesh_, faceI, n);
    const vector t2 = (n ^ t1)/max(mag(n ^ t1), VSMALL);
    const label typeId = p.typeId();
    const scalar mass = constProps(typeId).mass();

    p.U() =
        sqrt(physicoChemical::k.value()*wall.temperature/mass)
       *(rndGen_.GaussNormal<scalar>()*t1 + rndGen_.GaussNormal<scalar>()*t2
         + sqrt(-2.0*log(max(1.0 - rndGen_.sample01<scalar>(), VSMALL)))*n)
      + wall.velocity;

    p.ERot() = equipartitionRotationalEnergy(wall.temperature, constProps(typeId).rotationalDegreesOfFreedom());
    p.vibLevel() = equipartitionVibrationalEnergyLevel(wall.temperature, constProps(typeId).nVibrationalModes(), typeId);
    p.ELevel() = equipartitionElectronicLevel(wall.temperature, constProps(typeId).electronicDegeneracyList(), constProps(typeId).electronicEnergyList());
    td.keepParticle = true;
}

void Foam::dsmcCloud::calculateAndWriteTimeStep
(
    const Time& runTime,
    const dictionary& dsmcInitialiseDict
)
{
    if (!dsmcInitialiseDict.found("configurations"))
    {
        return;
    }

    const PtrList<entry> configurations(Foam::hyCompat::lookup(Foam::hyCompat::lookup(dsmcInitialiseDict, "configurations")));

    if (configurations.empty() || !configurations[0].isDict())
    {
        return;
    }

    const dictionary& configurationDict = configurations[0].dict();

    if (!configurationDict.found("numberDensities"))
    {
        return;
    }

    const scalar temperature =
        configurationDict.lookupOrDefault<scalar>("translationalTemperature", 0.0);
    const vector flowVelocity =
        configurationDict.lookupOrDefault<vector>("velocity", vector::zero);
    const scalar flowSpeed = mag(flowVelocity);
    const dictionary& numberDensitiesDict = configurationDict.subDict("numberDensities");
    const dictionary& moleculePropertiesDict = particleProperties_.subDict("moleculeProperties");

    scalar totalNumberDensity = 0.0;
    scalar totalMassDensity = 0.0;

    forAllConstIters(numberDensitiesDict, iter)
    {
        const word& speciesName = iter().keyword();
        const scalar numberDensity = numberDensitiesDict.get<scalar>(speciesName);

        totalNumberDensity += numberDensity;

        if (!moleculePropertiesDict.found(speciesName))
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Species '" << speciesName << "' not found in moleculeProperties. "
                    << "Skipping it in automatic time-step estimation." << endl;
            }

            continue;
        }

        const dictionary& speciesDict = moleculePropertiesDict.subDict(speciesName);
        const scalar molecularMass = speciesDict.get<scalar>("mass");
        totalMassDensity += numberDensity*molecularMass;

        if (Pstream::master())
        {
            Info<< "  Species " << speciesName
                << ": mass = " << molecularMass
                << " kg, numberDensity = " << numberDensity << " m^-3" << endl;
        }
    }

    if (Pstream::parRun())
    {
        reduce(totalNumberDensity, sumOp<scalar>());
        reduce(totalMassDensity, sumOp<scalar>());
    }

    if (totalNumberDensity <= SMALL || totalMassDensity <= SMALL)
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Invalid freestream density data. Skipping automatic time-step estimation."
                << endl;
        }

        return;
    }

    const scalar averageMolecularMass = totalMassDensity/totalNumberDensity;
    const scalar vThermal =
        max(maxwellianMostProbableSpeed(max(temperature, SMALL), averageMolecularMass), SMALL);
    const scalar characteristicSpeed = max(flowSpeed + vThermal, SMALL);

    const boundBox meshBounds(mesh_.points(), false);
    const vector domainSize = meshBounds.max() - meshBounds.min();
    const scalar xlen = domainSize.x();
    const scalar ylen = domainSize.y();
    const scalar zlen = domainSize.z();
    const scalar maxDomainLength = max(xlen, max(ylen, zlen));

    scalar totalVolume = 0.0;
    scalar minVolume = GREAT;
    scalar maxVolume = -GREAT;

    forAll(mesh_.cells(), cellI)
    {
        const scalar volume = mesh_.V()[cellI];
        totalVolume += volume;
        minVolume = min(minVolume, volume);
        maxVolume = max(maxVolume, volume);
    }

    if (Pstream::parRun())
    {
        reduce(totalVolume, sumOp<scalar>());
        reduce(minVolume, minOp<scalar>());
        reduce(maxVolume, maxOp<scalar>());
    }

    const label globalTotalCells = returnReduce(mesh_.nCells(), sumOp<label>());

    if (globalTotalCells <= 0 || totalVolume <= SMALL)
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Invalid mesh statistics. Skipping automatic time-step estimation." << endl;
        }

        return;
    }

    const scalar meanCellVolume = totalVolume/scalar(globalTotalCells);
    const scalar characteristicLength = cbrt(meanCellVolume);
    const dictionary& controlDict = runTime.controlDict();
    const scalar CAC = controlDict.lookupOrDefault<scalar>("CAC", 1.0);
    const scalar CTC = controlDict.lookupOrDefault<scalar>("CTC", 1.0);
    const label sampleSteps = controlDict.lookupOrDefault<label>("sampleSteps", 500);

    const scalar baseDeltaT = characteristicLength/characteristicSpeed;
    const scalar newDeltaT = CTC*baseDeltaT;

    if (newDeltaT <= SMALL)
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Calculated invalid time step: " << newDeltaT
                << ". Skipping automatic time-step estimation." << endl;
        }

        return;
    }

    const label steadyStateSteps =
        max(label(1), label(CAC*maxDomainLength/vThermal/newDeltaT + 1.0));
    const scalar steadyStateTime = steadyStateSteps*newDeltaT;
    const scalar totalTime = steadyStateTime + sampleSteps*newDeltaT;

    if (Pstream::master())
    {
        Info<< nl << "========== DSMC Time Step Calculation ==========" << endl;
        Info<< "Input parameters:" << endl;
        Info<< "  Translational temperature: " << temperature << " K" << endl;
        Info<< "  Flow velocity: " << flowSpeed << " m/s" << endl;
        Info<< "  Total number density: " << totalNumberDensity << " m^-3" << endl;
        Info<< "Domain dimensions:" << endl;
        Info<< "  x-length: " << xlen << " m" << endl;
        Info<< "  y-length: " << ylen << " m" << endl;
        Info<< "  z-length: " << zlen << " m" << endl;
        Info<< "  Max domain length: " << maxDomainLength << " m" << endl;
        Info<< "  Total number of cells: " << globalTotalCells << endl;
        Info<< "Calculated values:" << endl;
        Info<< "  Average molecular mass: " << averageMolecularMass << " kg" << endl;
        Info<< "  Thermal velocity: " << vThermal << " m/s" << endl;
        Info<< "  Characteristic speed: " << characteristicSpeed << " m/s" << endl;
        Info<< "  Mean cell volume: " << meanCellVolume << " m^3" << endl;
        Info<< "  Max cell volume: " << maxVolume << " m^3" << endl;
        Info<< "  Min cell volume: " << minVolume << " m^3" << endl;
        if (minVolume > SMALL)
        {
            Info<< "  Max/Min: " << maxVolume/minVolume << endl;
        }
        Info<< "  Characteristic length: " << characteristicLength << " m" << endl;
        Info<< "Time step adjustment:" << endl;
        Info<< "  Original time step: " << runTime.deltaTValue() << " s" << endl;
        Info<< "  New time step: " << newDeltaT
            << " s before CTC-" << CTC << " modification is " << baseDeltaT << endl;
        Info<< "  Steady state steps (naver): " << steadyStateSteps << endl;
        Info<< "  Sample steps: " << sampleSteps << endl;
        Info<< "  Total steps: " << steadyStateSteps + sampleSteps << endl;
        Info<< "------------------------------------------------" << endl;
        Info<< "For control file parameter replace:" << endl;
        Info<< "  deltaT               " << newDeltaT << endl;
        Info<< "  steadyStateTime      " << steadyStateTime << endl;
        Info<< "  totalTime            " << totalTime << endl;
        Info<< "================================================" << nl << endl;
    }
}

void Foam::dsmcCloud::initialiseFromDict(const dictionary& dsmcInitialiseDict)
{
    if (!dsmcInitialiseDict.found("configurations"))
    {
        FatalIOErrorInFunction(dsmcInitialiseDict)
            << "dsmcInitialiseDict is missing configurations" << nl
            << exit(FatalIOError);
    }

    calculateAndWriteTimeStep(mesh_.time(), dsmcInitialiseDict);

    const PtrList<entry> configs(Foam::hyCompat::lookup(Foam::hyCompat::lookup(dsmcInitialiseDict, "configurations")));

    forAll(configs, configI)
    {
        if (!configs[configI].isDict())
        {
            continue;
        }

        const dictionary& config = configs[configI].dict();

        if (config.get<word>("type") != "dsmcMeshFill")
        {
            FatalIOErrorInFunction(dsmcInitialiseDict)
                << "Only dsmcMeshFill is supported in stage 2" << nl
                << exit(FatalIOError);
        }

        const scalar translationalTemperature = config.get<scalar>("translationalTemperature");
        const scalar rotationalTemperature = config.lookupOrDefault<scalar>("rotationalTemperature", translationalTemperature);
        const scalar vibrationalTemperature = config.lookupOrDefault<scalar>("vibrationalTemperature", translationalTemperature);
        const scalar electronicTemperature = config.lookupOrDefault<scalar>("electronicTemperature", 0.0);
        const vector velocity(config.get<vector>("velocity"));
        const dictionary& numberDensitiesDict(config.subDict("numberDensities"));
        const wordList molecules(numberDensitiesDict.toc());
        scalarList numberDensities(molecules.size(), 0.0);

        forAll(molecules, i)
        {
            numberDensities[i] = numberDensitiesDict.get<scalar>(molecules[i]);
        }

        forAll(mesh_.cells(), cellI)
        {
            const List<tetIndices> cellTets(polyMeshTetDecomposition::cellTetIndices(mesh_, cellI));

            forAll(cellTets, tetI)
            {
                const tetIndices& cellTetIs = cellTets[tetI];
                const tetPointRef tet = cellTetIs.tet(mesh_);
                const scalar tetVolume = tet.mag();

                forAll(molecules, i)
                {
                    const label typeId = typeIdList_.find(molecules[i]);

                    if (typeId < 0)
                    {
                        FatalIOErrorInFunction(dsmcInitialiseDict)
                            << "Unknown typeId " << molecules[i] << nl
                            << exit(FatalIOError);
                    }

                    const dsmcParcel::constantProperties& cP = constProps(typeId);
                    const scalar particlesRequired = numberDensities[i]*tetVolume/nParticle_;
                    label nInsert = label(particlesRequired);

                    if ((particlesRequired - nInsert) > rndGen_.sample01<scalar>())
                    {
                        ++nInsert;
                    }

                    for (label pI = 0; pI < nInsert; ++pI)
                    {
                        const barycentric coordinates = randomTetCoordinates(rndGen_);
                        vector U = equipartitionLinearVelocity(translationalTemperature, cP.mass()) + velocity;
                        const scalar ERot = equipartitionRotationalEnergy(rotationalTemperature, cP.rotationalDegreesOfFreedom());
                        const labelList vibLevel = equipartitionVibrationalEnergyLevel(vibrationalTemperature, cP.nVibrationalModes(), typeId);
                        const label ELevel = equipartitionElectronicLevel(electronicTemperature, cP.electronicDegeneracyList(), cP.electronicEnergyList());

                        addNewParcel(coordinates, U, 1.0, ERot, ELevel, cellI, cellTetIs.face(), cellTetIs.tetPt(), typeId, -1, 0, vibLevel);
                    }
                }
            }
        }

        const label mostAbundantType = findMax(numberDensities);
        const dsmcParcel::constantProperties& cP = constProps(mostAbundantType);
        sigmaTcRMax_.primitiveFieldRef() = cP.sigmaT()*maxwellianMostProbableSpeed(translationalTemperature, cP.mass());
        sigmaTcRMax_.correctBoundaryConditions();
    }

    buildCellOccupancy();
    fields_.updateTimeInfo();
    fields_.calculateFields();
    fields_.writeFields();
}


void Foam::dsmcCloud::evolve()
{
    using clock_type = std::chrono::steady_clock;
    const auto tEvolveStart = clock_type::now();

    const auto logEvolveStage = [&](const char* stage)
    {
        if (evolveStageProbe_)
        {
            Pout<< "Evolve stage rank " << Pstream::myProcNo()
                << " timeIndex " << mesh_.time().timeIndex()
                << " time " << mesh_.time().value()
                << ": " << stage << nl << endl;
        }
    };

    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    dsmcParcel::trackingData td(*this);

    // Phase A: all ranks process boundaries.  Each boundary face is
    // handled only by the rank that owns its adjacent cell (checked
    // per-face in the boundary model / hitPatch).  This ensures each
    // injection / deletion happens exactly once.
    const bool processBoundaries = true;

    const auto t0 = clock_type::now();
    if (openmpEnabled_ && openmpMoveEnabled_)
    {
        beginMoveAppendCapture();
    }
    logEvolveStage("before controlBeforeMove");
    controllers_.controlBeforeMove();
    logEvolveStage("after controllers controlBeforeMove");
    if (processBoundaries)
    {
        boundaries_.controlBeforeMove();
    }
    logEvolveStage("after boundaries controlBeforeMove");
    if
    (
        openmpEnabled_
     && openmpMoveEnabled_
     && (
            particleLoadStart_.size() != ompNumThreads_
         || particleLoadEnd_.size() != ompNumThreads_
        )
    )
    {
        rebuildParticleLoadPartition();
    }
    const label preMoveCloudSize = this->size();
    if (profilingDetailEnabled_ && previousMoveEndCloudSize_ >= 0)
    {
        const label interStepDelta = preMoveCloudSize - previousMoveEndCloudSize_;
        moveInterStepCloudDeltaSum_ += interStepDelta;
        moveInterStepCloudDeltaMin_ = min(moveInterStepCloudDeltaMin_, interStepDelta);
        moveInterStepCloudDeltaMax_ = max(moveInterStepCloudDeltaMax_, interStepDelta);
    }
    const auto tMovePreEnd = clock_type::now();
    recordMovePhaseProfile
    (
        std::chrono::duration<scalar>(tMovePreEnd - t0).count(),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0
    );

    // Phase A: pre-move migration only needed on the first step (initial
    // distribution).  After that, the post-move migration from the previous
    // step already placed particles on their owning ranks.
    if (replicatedMeshActive() && replicatedMesh_->migrationCalls() == 0)
    {
        logEvolveStage("before pre-move migration (initial)");
        if (Pstream::parRun())
        {
            replicatedMesh_->distributeInitialParticles();
        }
        else
        {
            replicatedMesh_->migrateParticlesByCellOwner();
        }
        replicatedMesh_->updateParticleCounts();
        logEvolveStage("after pre-move migration (initial)");
    }

    // Delayed-receive: finish previous step's async migration
    bool postDoneInMigration = false;
    const bool delayedReceive = replicatedMeshActive()
        && mesh_.time().controlDict().lookupOrDefault<bool>
           ("replicatedMeshDelayedReceive", false);
    if (delayedReceive && replicatedMesh_->asyncMigrationPending())
    {
        replicatedMesh_->migrateFinish();
    }

    logEvolveStage("before move");
    Cloud<dsmcParcel>::move(*this, td, mesh_.time().deltaTValue());
    logEvolveStage("after move");
    if (openmpEnabled_ && openmpMoveEnabled_)
    {
        endMoveAppendCapture();
    }
    else if (replicatedMeshActive() && moveOrderedParcelsValid_)
    {
        moveOrderedParcelsValid_ = false;
    }
    const label moveEndCloudSize = this->size();
    const auto t1 = clock_type::now();

    // Phase A: migrate particles to their owning rank before building occupancy.
    if (replicatedMeshActive() && replicatedMesh_->stepCounter() > 0)
    {
        if (replicatedMesh_->stepCounter() % replicatedMesh_->migrateInterval() == 0)
        {
            if (delayedReceive)
            {
                // Post overlap: do post BEFORE migration to hide post time
                // behind slow ranks' move time.
                if (processBoundaries)
                {
                    fields_.calculateFields();
                    postDoneInMigration = true;
                }

                logEvolveStage("before migrateBegin");
                replicatedMesh_->migrateBegin();
                replicatedMesh_->updateParticleCounts();
                cellOccupancyMaterialized_ = false;
                logEvolveStage("after migrateBegin");

                logEvolveStage("before migrateFinish");
                replicatedMesh_->migrateFinish();
                logEvolveStage("after migrateFinish");
            }
            else
            {
                logEvolveStage("before migrateParticlesByCellOwner");
                replicatedMesh_->migrateParticlesByCellOwner();
                replicatedMesh_->updateParticleCounts();
                cellOccupancyMaterialized_ = false;
                logEvolveStage("after migrateParticlesByCellOwner");
            }
        }
        replicatedMesh_->advanceStepCounter();
    }
    else if (replicatedMeshActive() && replicatedMesh_->stepCounter() == 0)
    {
        logEvolveStage("before migrateParticlesByCellOwner (first step)");
        replicatedMesh_->migrateParticlesByCellOwner();
        replicatedMesh_->updateParticleCounts();
        replicatedMesh_->advanceStepCounter();
        logEvolveStage("after migrateParticlesByCellOwner (first step)");
    }

    // Phase B: reassign cellOwner_ at configured steps and redistribute
    if (replicatedMeshActive() && replicatedMesh_->rebalanceSteps().size())
    {
        const label currentStep = replicatedMesh_->stepCounter();
        const labelList& steps = replicatedMesh_->rebalanceSteps();
        forAll(steps, i)
        {
            if (steps[i] == currentStep)
            {
                Info<< "\nPhase B: reassigning cellOwner_ at step "
                    << currentStep << nl << endl;
                replicatedMesh_->reassignCellOwner();
                replicatedMesh_->migrateParticlesByCellOwner();
                replicatedMesh_->updateParticleCounts();
                Info<< "Phase B: redistribution complete\n" << endl;
                break;
            }
        }
    }

    logEvolveStage("before buildCellOccupancy");
    buildCellOccupancy();
    logEvolveStage("after buildCellOccupancy");
    const auto t2 = clock_type::now();

    logEvolveStage("before coordSystem evolve");
    coordSystem().evolve();
    logEvolveStage("after coordSystem evolve");
    const auto t3 = clock_type::now();

    logEvolveStage("before controlBeforeCollisions");
    controllers_.controlBeforeCollisions();
    logEvolveStage("after controllers controlBeforeCollisions");
    if (processBoundaries)
    {
        boundaries_.controlBeforeCollisions();
    }
    logEvolveStage("after boundaries controlBeforeCollisions");

    logEvolveStage("before collisions");
    beginCollisionPhase();
    collisions();
    endCollisionPhase();
    logEvolveStage("after collisions");

    const auto t4 = clock_type::now();

    if (reactionsActive() && emitStepDiagnostics_)
    {
        logEvolveStage("before reaction output");
        reactions().outputData();
        logEvolveStage("after reaction output");
    }
    const auto t5 = clock_type::now();

    logEvolveStage("before controlAfterCollisions");
    controllers_.controlAfterCollisions();
    logEvolveStage("after controllers controlAfterCollisions");
    if (processBoundaries)
    {
        boundaries_.controlAfterCollisions();
    }
    logEvolveStage("after boundaries controlAfterCollisions");

    logEvolveStage("before fields/output");
    if (processBoundaries && !postDoneInMigration)
    {
        fields_.calculateFields();
        fields_.writeFields();
    }
    else if (processBoundaries)
    {
        fields_.writeFields();
    }

    controllers_.calculateProps();
    controllers_.outputResults();

    if (processBoundaries)
    {
        boundaries_.calculateProps();
        boundaries_.outputResults();
    }

    boundaryMeas_.outputResults();
    logEvolveStage("after fields/output");

    trackingInfo_.clean();
    boundaryMeas_.clean();
    cellMeas_.clean();
    logEvolveStage("after cleanup");
    const label endStepCloudSize = this->size();

    if (profilingDetailEnabled_)
    {
        const label postStepDelta = endStepCloudSize - moveEndCloudSize;
        movePostStepCloudDeltaSum_ += postStepDelta;
        movePostStepCloudDeltaMin_ = min(movePostStepCloudDeltaMin_, postStepDelta);
        movePostStepCloudDeltaMax_ = max(movePostStepCloudDeltaMax_, postStepDelta);
    }

    previousMoveEndCloudSize_ = endStepCloudSize;

    const auto t6 = clock_type::now();
    scalar rankTimeWall = 0.0;
    scalar autoRebalanceWall = 0.0;

    if (evolveProfileEnabled_)
    {
        evolveMoveWallTime_ += std::chrono::duration<scalar>(t1 - t0).count();
        evolveBuildWallTime_ += std::chrono::duration<scalar>(t2 - t1).count();
        evolveCoordWallTime_ += std::chrono::duration<scalar>(t3 - t2).count();
        evolveCollisionWallTime_ += std::chrono::duration<scalar>(t4 - t3).count();
        evolveReactionWallTime_ += std::chrono::duration<scalar>(t5 - t4).count();
        evolvePostWallTime_ += std::chrono::duration<scalar>(t6 - t5).count();
        ++evolveProfileCalls_;

        const label localMoveParcels = sum(moveThreadParticleCounts_);
        const label localCollisionCandidates = sum(collisionThreadCandidateCounts_);
        const label localAcceptedCollisions = sum(collisionThreadAcceptedCounts_);

        if (profilingDetailEnabled_ && emitStepDiagnostics_ && Pstream::parRun())
        {
            Pout<< "Load stats rank " << Pstream::myProcNo() << ":" << nl
                << "    move particles processed      = " << localMoveParcels << nl
                << "    collision candidates          = " << localCollisionCandidates << nl
                << "    accepted collisions           = " << localAcceptedCollisions << nl
                << endl;
        }

        if (profilingDetailEnabled_ && emitStepDiagnostics_)
        {
            scalar moveMax = localMoveParcels;
            scalar moveMin = localMoveParcels;
            scalar candMax = localCollisionCandidates;
            scalar candMin = localCollisionCandidates;
            scalar collMax = localAcceptedCollisions;
            scalar collMin = localAcceptedCollisions;

            if (Pstream::parRun())
            {
                reduce(moveMax, maxOp<scalar>());
                reduce(moveMin, minOp<scalar>());
                reduce(candMax, maxOp<scalar>());
                reduce(candMin, minOp<scalar>());
                reduce(collMax, maxOp<scalar>());
                reduce(collMin, minOp<scalar>());
            }

            if (Pstream::master())
            {
                Info<< "Load balance summary:" << nl
                    << "    move particles max/min        = "
                    << moveMax << " / " << moveMin << nl
                    << "    move imbalance max/min        = "
                    << (moveMin > SMALL ? moveMax/moveMin : 0.0) << nl
                    << "    collision cand max/min        = "
                    << candMax << " / " << candMin << nl
                    << "    collision cand imbalance      = "
                    << (candMin > SMALL ? candMax/candMin : 0.0) << nl
                    << "    accepted coll max/min         = "
                    << collMax << " / " << collMin << nl
                    << "    accepted coll imbalance       = "
                    << (collMin > SMALL ? collMax/collMin : 0.0) << nl;

                if (!Pstream::parRun() && moveThreadParticleCounts_.size())
                {
                    Info<< "    thread move particles         = " << moveThreadParticleCounts_ << nl
                        << "    thread collision candidates   = " << collisionThreadCandidateCounts_ << nl
                        << "    thread accepted collisions    = " << collisionThreadAcceptedCounts_ << nl;
                }

                Info<< endl;
            }

        }
    }

    // Per-rank evolve wall time (for load balance diagnostics)
    if (replicatedMeshActive())
    {
        const auto tRankTime0 = clock_type::now();
        replicatedMesh_->addEvolveTime
        (
            std::chrono::duration<scalar>(clock_type::now() - tEvolveStart).count()
        );
        const auto tRankTime1 = clock_type::now();
        rankTimeWall = std::chrono::duration<scalar>(tRankTime1 - tRankTime0).count();

        // Phase C: automatic DLB — checks per-rank wall time imbalance
        // and triggers Hilbert SFC rebalancing if threshold exceeded.
        // Must be called after addEvolveTime so the current step's time
        // is included in the imbalance calculation.
        const label prevRebalances = replicatedMesh_->autoRebalanceCount();
        const auto tAutoRebalance0 = clock_type::now();
        replicatedMesh_->autoRebalance();
        const auto tAutoRebalance1 = clock_type::now();
        autoRebalanceWall =
            std::chrono::duration<scalar>(tAutoRebalance1 - tAutoRebalance0).count();
        if (replicatedMesh_->autoRebalanceCount() > prevRebalances)
        {
            // DLB triggered migration — cellOccupancy has dangling pointers
            cellOccupancyMaterialized_ = false;
            clearMoveOrderedParcels();
        }
    }

    const auto tEvolveEnd = clock_type::now();
    if (evolveProfileEnabled_)
    {
        const scalar postProfileWall =
            std::chrono::duration<scalar>(tEvolveEnd - t6).count();

        evolvePreWallTime_ += std::chrono::duration<scalar>(t0 - tEvolveStart).count();
        evolveRankTimeWallTime_ += rankTimeWall;
        evolveAutoRebalanceWallTime_ += autoRebalanceWall;
        evolvePostProfileResidualWallTime_ +=
            postProfileWall - rankTimeWall - autoRebalanceWall;
        evolveFullWallTime_ +=
            std::chrono::duration<scalar>(tEvolveEnd - tEvolveStart).count();
    }
}


void Foam::dsmcCloud::loadBalanceCheck()
{}


void Foam::dsmcCloud::info() const
{
    const label localParcels = this->size();
    label nParcels = localParcels;
    reduce(nParcels, sumOp<label>());

    scalar mass = 0.0;
    scalar linearKineticEnergy = 0.0;
    scalar rotationalEnergy = 0.0;
    scalar vibrationalEnergy = 0.0;
    scalar electronicEnergy = 0.0;

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();
        const dsmcParcel::constantProperties& cP = constProps(p.typeId());

        mass += cP.mass()*nParticle_;
        linearKineticEnergy += 0.5*cP.mass()*(p.U() & p.U())*nParticle_;
        rotationalEnergy += p.ERot()*nParticle_;
        vibrationalEnergy += cP.eVib_tot(p.vibLevel())*nParticle_;
        electronicEnergy += cP.electronicEnergyList()[p.ELevel()]*nParticle_;
    }

    reduce(mass, sumOp<scalar>());
    reduce(linearKineticEnergy, sumOp<scalar>());
    reduce(rotationalEnergy, sumOp<scalar>());
    reduce(vibrationalEnergy, sumOp<scalar>());
    reduce(electronicEnergy, sumOp<scalar>());

    if (replicatedMeshActive() && !Pstream::parRun())
    {
        const label myRank = replicatedMesh_->myRank();
        Info<< "Cloud name: " << this->name()
            << " [rank " << myRank << "]" << nl
            << "    Number of dsmc particles        = " << localParcels << nl;
    }
    else
    {
        if (!isOutputRank()) return;

        Info<< "Cloud name: " << this->name() << nl
            << "    Number of dsmc particles        = " << nParcels << nl;
    }

    if (nParcels)
    {
        const scalar nMol = nParcels*nParticle_;

        Info<< "    Number of molecules             = " << nMol << nl
            << "    Mass in system                  = " << mass << nl
            << "    Average linear kinetic energy   = " << linearKineticEnergy/nMol << nl
            << "    Average rotational energy       = " << rotationalEnergy/nMol << nl
            << "    Average vibrational energy      = " << vibrationalEnergy/nMol << nl
            << "    Average electronic energy       = " << electronicEnergy/nMol << nl
            << "    Average total energy            = "
            << (linearKineticEnergy + rotationalEnergy + vibrationalEnergy + electronicEnergy)/nMol
            << endl;
    }
}


void Foam::dsmcCloud::recordMoveThreadCounts(const labelList& counts)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    if (moveThreadParticleCounts_.size() != counts.size())
    {
        moveThreadParticleCounts_.setSize(counts.size(), 0);
    }

    forAll(counts, i)
    {
        moveThreadParticleCounts_[i] += counts[i];
    }
}


void Foam::dsmcCloud::recordMoveThreadProfile
(
    const labelList& counts,
    const scalarField& wallTimes
)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    recordMoveThreadCounts(counts);

    if (moveThreadWallTimes_.size() != wallTimes.size())
    {
        moveThreadWallTimes_.setSize(wallTimes.size(), 0.0);
    }

    forAll(wallTimes, i)
    {
        moveThreadWallTimes_[i] += wallTimes[i];
    }

    moveLastThreadParticleCounts_ = counts;
    moveLastThreadWallTimes_ = wallTimes;
}


void Foam::dsmcCloud::recordMoveInnerProfile
(
    const scalarField& trackWallTimes,
    const scalarField& trackerWallTimes,
    const scalarField& boundaryWallTimes,
    const labelList& faceHitCounts,
    const labelList& cyclicHitCounts,
    const labelList& stuckHitCounts,
    const labelList& patchHitCounts,
    const labelList& processorHitCounts
)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    auto accumulateScalarField =
    [](scalarField& total, const scalarField& values)
    {
        if (total.size() != values.size())
        {
            total.setSize(values.size(), 0.0);
        }

        forAll(values, i)
        {
            total[i] += values[i];
        }
    };

    auto accumulateLabelList =
    [](labelList& total, const labelList& values)
    {
        if (total.size() != values.size())
        {
            total.setSize(values.size(), 0);
        }

        forAll(values, i)
        {
            total[i] += values[i];
        }
    };

    accumulateScalarField(moveThreadTrackWallTimes_, trackWallTimes);
    accumulateScalarField(moveThreadTrackerWallTimes_, trackerWallTimes);
    accumulateScalarField(moveThreadBoundaryWallTimes_, boundaryWallTimes);
    accumulateLabelList(moveThreadFaceHitCounts_, faceHitCounts);
    accumulateLabelList(moveThreadCyclicHitCounts_, cyclicHitCounts);
    accumulateLabelList(moveThreadStuckHitCounts_, stuckHitCounts);
    accumulateLabelList(moveThreadPatchHitCounts_, patchHitCounts);
    accumulateLabelList(moveThreadProcessorHitCounts_, processorHitCounts);

    forAll(trackWallTimes, i)
    {
        moveTrackWallTime_ += trackWallTimes[i];
    }

    forAll(trackerWallTimes, i)
    {
        moveTrackerCallbackWallTime_ += trackerWallTimes[i];
    }

    forAll(boundaryWallTimes, i)
    {
        moveBoundaryControlWallTime_ += boundaryWallTimes[i];
    }

    forAll(faceHitCounts, i)
    {
        moveFaceHitsSum_ += faceHitCounts[i];
    }

    forAll(cyclicHitCounts, i)
    {
        moveCyclicHitsSum_ += cyclicHitCounts[i];
    }

    forAll(stuckHitCounts, i)
    {
        moveStuckHitsSum_ += stuckHitCounts[i];
    }

    forAll(patchHitCounts, i)
    {
        movePatchHitsSum_ += patchHitCounts[i];
    }

    forAll(processorHitCounts, i)
    {
        moveProcessorHitsSum_ += processorHitCounts[i];
    }
}


void Foam::dsmcCloud::recordMovePhaseProfile
(
    const scalar preControlWallTime,
    const scalar resetSetupWallTime,
    const scalar extractWallTime,
    const scalar kernelWallTime,
    const scalar commitWallTime,
    const scalar transferFinalizeWallTime
)
{
    if (!evolveProfileEnabled_)
    {
        return;
    }

    movePreControlWallTime_ += preControlWallTime;
    moveResetSetupWallTime_ += resetSetupWallTime;
    moveExtractWallTime_ += extractWallTime;
    moveKernelWallTime_ += kernelWallTime;
    moveCommitWallTime_ += commitWallTime;
    moveTransferFinalizeWallTime_ += transferFinalizeWallTime;
    ++moveProfileCalls_;
}


void Foam::dsmcCloud::recordCollisionThreadCounts
(
    const labelList& candidateCounts,
    const labelList& acceptedCounts
)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    if (collisionThreadCandidateCounts_.size() != candidateCounts.size())
    {
        collisionThreadCandidateCounts_.setSize(candidateCounts.size(), 0);
    }

    if (collisionThreadAcceptedCounts_.size() != acceptedCounts.size())
    {
        collisionThreadAcceptedCounts_.setSize(acceptedCounts.size(), 0);
    }

    forAll(candidateCounts, i)
    {
        collisionThreadCandidateCounts_[i] += candidateCounts[i];
    }

    forAll(acceptedCounts, i)
    {
        collisionThreadAcceptedCounts_[i] += acceptedCounts[i];
    }
}


void Foam::dsmcCloud::recordCollisionThreadProfile
(
    const labelList& candidateCounts,
    const labelList& acceptedCounts,
    const labelList& activeCellCounts,
    const labelList& reactionHitCounts,
    const scalarField& wallTimes
)
{
    if (!profilingDetailEnabled_)
    {
        return;
    }

    recordCollisionThreadCounts(candidateCounts, acceptedCounts);

    if (collisionThreadActiveCellCounts_.size() != activeCellCounts.size())
    {
        collisionThreadActiveCellCounts_.setSize(activeCellCounts.size(), 0);
    }

    if (collisionThreadReactionHitCounts_.size() != reactionHitCounts.size())
    {
        collisionThreadReactionHitCounts_.setSize(reactionHitCounts.size(), 0);
    }

    if (collisionThreadWallTimes_.size() != wallTimes.size())
    {
        collisionThreadWallTimes_.setSize(wallTimes.size(), 0.0);
    }

    forAll(activeCellCounts, i)
    {
        collisionThreadActiveCellCounts_[i] += activeCellCounts[i];
    }

    forAll(reactionHitCounts, i)
    {
        collisionThreadReactionHitCounts_[i] += reactionHitCounts[i];
    }

    forAll(wallTimes, i)
    {
        collisionThreadWallTimes_[i] += wallTimes[i];
    }
}


void Foam::dsmcCloud::accumulateCollisionCounts
(
    label candidates,
    label collisions
)
{
    cumulativeCollisionCandidates_ = candidates;
    cumulativeCollisions_ = collisions;
}


void Foam::dsmcCloud::resetLoadStats()
{
    moveThreadParticleCounts_.clear();
    moveThreadWallTimes_.clear();
    moveThreadTrackWallTimes_.clear();
    moveThreadTrackerWallTimes_.clear();
    moveThreadBoundaryWallTimes_.clear();
    moveThreadFaceHitCounts_.clear();
    moveThreadCyclicHitCounts_.clear();
    moveThreadStuckHitCounts_.clear();
    moveThreadPatchHitCounts_.clear();
    moveThreadProcessorHitCounts_.clear();
    collisionThreadCandidateCounts_.clear();
    collisionThreadAcceptedCounts_.clear();
    collisionThreadActiveCellCounts_.clear();
    collisionThreadReactionHitCounts_.clear();
    collisionThreadWallTimes_.clear();
}


void Foam::dsmcCloud::reportMoveCellHotspots() const
{
    if (!profilingDetailEnabled_ || moveItersPerCellCumulative_.empty())
    {
        return;
    }

    if (Pstream::parRun() && !replicatedMeshActive())
    {
        return;
    }

    labelList cumulativeMoveIters(moveItersPerCellCumulative_);
    labelList windowMoveIters(moveItersPerCell_);

    if (Pstream::parRun())
    {
        reduce(cumulativeMoveIters, sumOp<labelList>());
        reduce(windowMoveIters, sumOp<labelList>());
    }

    if (!isOutputRank())
    {
        return;
    }

    scalar totalMoveIters = 0.0;
    scalar windowTotalMoveIters = 0.0;
    label activeCells = 0;

    forAll(cumulativeMoveIters, cellI)
    {
        totalMoveIters += scalar(cumulativeMoveIters[cellI]);
        windowTotalMoveIters += scalar(windowMoveIters[cellI]);

        if (cumulativeMoveIters[cellI] > 0)
        {
            ++activeCells;
        }
    }

    if (totalMoveIters <= 0)
    {
        return;
    }

    const dictionary& controlDict = mesh_.time().controlDict();
    const label nTop =
        min
        (
            max(controlDict.lookupOrDefault<label>("moveHotspotTopCells", 20), label(0)),
            mesh_.nCells()
        );
    const label maxBoundaryLayer =
        max(controlDict.lookupOrDefault<label>("moveHotspotBoundaryLayers", 3), label(0));

    labelList hotspotBoundaryPatch(mesh_.boundaryMesh().size(), 0);
    label nHotspotBoundaryPatches = 0;

    const labelList& patchBoundaryIds = boundaries_.patchBoundaryIds();
    forAll(patchBoundaryIds, i)
    {
        const label patchI = patchBoundaryIds[i];

        if
        (
            patchI >= 0
         && patchI < hotspotBoundaryPatch.size()
         && !hotspotBoundaryPatch[patchI]
        )
        {
            hotspotBoundaryPatch[patchI] = 1;
            ++nHotspotBoundaryPatches;
        }
    }

    const labelList& generalBoundaryIds = boundaries_.generalBoundaryIds();
    forAll(generalBoundaryIds, i)
    {
        const label patchI = generalBoundaryIds[i];

        if
        (
            patchI >= 0
         && patchI < hotspotBoundaryPatch.size()
         && !hotspotBoundaryPatch[patchI]
        )
        {
            hotspotBoundaryPatch[patchI] = 1;
            ++nHotspotBoundaryPatches;
        }
    }

    const labelList& cyclicBoundaryIds = boundaries_.cyclicBoundaryIds();
    forAll(cyclicBoundaryIds, i)
    {
        const label patchI = cyclicBoundaryIds[i];

        if
        (
            patchI >= 0
         && patchI < hotspotBoundaryPatch.size()
         && !hotspotBoundaryPatch[patchI]
        )
        {
            hotspotBoundaryPatch[patchI] = 1;
            ++nHotspotBoundaryPatches;
        }
    }

    labelList boundaryLayer(mesh_.nCells(), -1);
    DynamicList<label> frontierCells;
    frontierCells.setCapacity(min(mesh_.nCells(), mesh_.boundaryMesh().size()*8));

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if (!hotspotBoundaryPatch[patchI] || isA<processorPolyPatch>(pp))
        {
            continue;
        }

        const labelUList& faceCells = pp.faceCells();

        forAll(faceCells, i)
        {
            const label cellI = faceCells[i];

            if (cellI >= 0 && cellI < boundaryLayer.size() && boundaryLayer[cellI] == -1)
            {
                boundaryLayer[cellI] = 0;
                frontierCells.append(cellI);
            }
        }
    }

    for (label layer = 1; layer <= maxBoundaryLayer && frontierCells.size(); ++layer)
    {
        DynamicList<label> nextFrontier;
        nextFrontier.setCapacity(frontierCells.size()*2 + 1);

        forAll(frontierCells, frontierI)
        {
            const label cellI = frontierCells[frontierI];
            const cell& cFaces = mesh_.cells()[cellI];

            forAll(cFaces, faceI)
            {
                const label meshFaceI = cFaces[faceI];

                if (!mesh_.isInternalFace(meshFaceI))
                {
                    continue;
                }

                const label owner = mesh_.faceOwner()[meshFaceI];
                const label neighbour = mesh_.faceNeighbour()[meshFaceI];
                const label otherCellI = owner == cellI ? neighbour : owner;

                if
                (
                    otherCellI >= 0
                 && otherCellI < boundaryLayer.size()
                 && boundaryLayer[otherCellI] == -1
                )
                {
                    boundaryLayer[otherCellI] = layer;
                    nextFrontier.append(otherCellI);
                }
            }
        }

        frontierCells.transfer(nextFrontier);
    }

    const label nLayerBuckets = maxBoundaryLayer + 2;
    scalarField layerMoveIters(nLayerBuckets, 0.0);
    labelList layerActiveCells(nLayerBuckets, 0);

    const bool haveCellOwners =
        replicatedMeshActive()
     && replicatedMesh().cellOwner().size() == cumulativeMoveIters.size();
    const labelList* cellOwnersPtr =
        haveCellOwners ? &replicatedMesh().cellOwner() : nullptr;
    scalarField ownerMoveIters
    (
        haveCellOwners ? replicatedMesh().nProcs() : 0,
        0.0
    );
    labelList ownerActiveCells
    (
        haveCellOwners ? replicatedMesh().nProcs() : 0,
        0
    );

    labelList topCells(nTop, -1);
    labelList topCounts(nTop, 0);

    forAll(cumulativeMoveIters, cellI)
    {
        const label count = cumulativeMoveIters[cellI];

        if (count <= 0)
        {
            continue;
        }

        const label layer = boundaryLayer[cellI];
        const label layerBucket =
            (layer >= 0 && layer <= maxBoundaryLayer)
          ? layer
          : maxBoundaryLayer + 1;

        layerMoveIters[layerBucket] += scalar(count);
        ++layerActiveCells[layerBucket];

        if (haveCellOwners)
        {
            const label owner = (*cellOwnersPtr)[cellI];

            if (owner >= 0 && owner < ownerMoveIters.size())
            {
                ownerMoveIters[owner] += scalar(count);
                ++ownerActiveCells[owner];
            }
        }

        forAll(topCounts, topI)
        {
            if (count > topCounts[topI])
            {
                for (label shiftI = topCounts.size() - 1; shiftI > topI; --shiftI)
                {
                    topCounts[shiftI] = topCounts[shiftI - 1];
                    topCells[shiftI] = topCells[shiftI - 1];
                }

                topCounts[topI] = count;
                topCells[topI] = cellI;
                break;
            }
        }
    }

    Info<< "Move cell hotspot summary:" << nl
        << "    cumulative move iters         = " << totalMoveIters << nl
        << "    current DLB-window move iters = " << windowTotalMoveIters << nl
        << "    active move cells             = "
        << activeCells << " / " << mesh_.nCells() << nl
        << "    hotspot boundary seed patches = "
        << nHotspotBoundaryPatches << nl
        << "    boundary layer limit          = " << maxBoundaryLayer << nl;

    if (nTop > 0)
    {
        Info<< "    top cells by cumulative move iters:" << nl;

        forAll(topCells, topI)
        {
            const label cellI = topCells[topI];

            if (cellI < 0)
            {
                continue;
            }

            const label count = topCounts[topI];
            const scalar share = 100.0*scalar(count)/max(totalMoveIters, scalar(1));
            const label layer = boundaryLayer[cellI];

            Info<< "        rank=" << topI + 1
                << " cell=" << cellI
                << " count=" << count
                << " share[%]=" << share
                << " boundaryLayer=";

            if (layer >= 0 && layer <= maxBoundaryLayer)
            {
                Info<< layer;
            }
            else
            {
                Info<< ">" << maxBoundaryLayer;
            }

            Info<< " owner=";

            if (haveCellOwners)
            {
                Info<< (*cellOwnersPtr)[cellI];
            }
            else
            {
                Info<< -1;
            }

            Info<< " centre=" << mesh_.cellCentres()[cellI]
                << " volume=" << mesh_.cellVolumes()[cellI]
                << nl;
        }
    }

    Info<< "    boundary-layer cumulative move iters:" << nl;

    forAll(layerMoveIters, bucketI)
    {
        const scalar share =
            100.0*layerMoveIters[bucketI]/max(totalMoveIters, scalar(1));

        Info<< "        layer=";

        if (bucketI <= maxBoundaryLayer)
        {
            Info<< bucketI;
        }
        else
        {
            Info<< ">" << maxBoundaryLayer;
        }

        Info<< " count=" << layerMoveIters[bucketI]
            << " share[%]=" << share
            << " activeCells=" << layerActiveCells[bucketI]
            << nl;
    }

    if (haveCellOwners)
    {
        scalarField ownerShares(ownerMoveIters.size(), 0.0);

        forAll(ownerMoveIters, ownerI)
        {
            ownerShares[ownerI] =
                100.0*ownerMoveIters[ownerI]/max(totalMoveIters, scalar(1));
        }

        Info<< "    current cellOwner cumulative move iters = "
            << ownerMoveIters << nl
            << "    current cellOwner share[%]              = "
            << ownerShares << nl
            << "    current cellOwner active move cells     = "
            << ownerActiveCells << nl;
    }

    Info<< endl;
}


void Foam::dsmcCloud::deleteParcel(dsmcParcel* p)
{
    // Cast to IDLList to resolve multiple-inheritance ambiguity for erase()
    IDLList<dsmcParcel>& list = static_cast<IDLList<dsmcParcel>&>(*this);
    list.remove(p);
    delete p;
}


void Foam::dsmcCloud::reportProfiling() const
{
    if (moveProfileCalls_ > 0 && isOutputRank())
    {
        const scalar totalProfiled =
            movePreControlWallTime_
          + moveResetSetupWallTime_
          + moveExtractWallTime_
          + moveKernelWallTime_
          + moveCommitWallTime_
          + moveTransferFinalizeWallTime_;

        Info<< "Move profiling summary:" << nl
            << "    move calls                    = " << moveProfileCalls_ << nl
            << "    move pre-control/partition [s]= " << movePreControlWallTime_ << nl
            << "    move reset/setup [s]          = " << moveResetSetupWallTime_ << nl
            << "    move extract parcels [s]      = " << moveExtractWallTime_ << nl
            << "    move parallel kernel wall [s] = " << moveKernelWallTime_ << nl
            << "    move commit/survivor rebuild [s] = " << moveCommitWallTime_ << nl
            << "    move transfer/delete finalize [s] = "
            << moveTransferFinalizeWallTime_ << nl
            << "    total profiled [s]            = " << totalProfiled << nl
            << endl;

        if (profilingDetailEnabled_)
        {
            Info<< "Move profiling detail:" << nl
                << "    move commit passes total      = "
                << moveCommitCountCalls_ << nl
                << "    avg move passes/call          = "
                << moveLoopPassesSum_/max(moveProfileCalls_, label(1)) << nl
                << "    max move passes/call          = "
                << moveLoopPassesMax_ << nl
                << "    avg move extracted parcels    = "
                << moveExtractedParcelsSum_/max(moveCommitCountCalls_, label(1)) << nl
                << "    extract deferred-only [s]     = "
                << moveExtractDeferredWallTime_ << nl
                << "    extract ordered-reuse [s]     = "
                << moveExtractOrderedReuseWallTime_ << nl
                << "    extract full-scan [s]         = "
                << moveExtractFullScanWallTime_ << nl
                << "    deferred extract passes       = "
                << moveExtractDeferredPasses_ << nl
                << "    ordered-reuse extract passes  = "
                << moveExtractOrderedReusePasses_ << nl
                << "    full-scan extract passes      = "
                << moveExtractFullScanPasses_ << nl
                << "    avg deferred-only parcels     = "
                << moveExtractDeferredParcelsSum_/max(moveExtractDeferredPasses_, label(1)) << nl
                << "    avg ordered-reuse parcels     = "
                << moveExtractOrderedReuseParcelsSum_/max(moveExtractOrderedReusePasses_, label(1)) << nl
                << "    avg full-scan parcels         = "
                << moveExtractFullScanParcelsSum_/max(moveExtractFullScanPasses_, label(1)) << nl
                << "    first-pass reuse checks       = "
                << moveFirstPassReuseChecks_ << nl
                << "    first-pass has ordered        = "
                << moveFirstPassHasOrderedCount_ << nl
                << "    first-pass missing ordered    = "
                << moveFirstPassMissingOrderedCount_ << nl
                << "    first-pass offset mismatch    = "
                << moveFirstPassOffsetMismatchCount_ << nl
                << "    first-pass exact size match   = "
                << moveFirstPassSizeMatchCount_ << nl
                << "    first-pass prior+appended match= "
                << moveFirstPassPriorPlusPendingMatchCount_ << nl
                << "    avg first-pass current size   = "
                << moveFirstPassCurrentSizeSum_/max(moveFirstPassReuseChecks_, label(1)) << nl
                << "    avg first-pass prior size     = "
                << moveFirstPassPriorSizeSum_/max(moveFirstPassReuseChecks_, label(1)) << nl
                << "    avg first-pass appended size  = "
                << moveFirstPassPendingSizeSum_/max(moveFirstPassReuseChecks_, label(1)) << nl
                << "    avg first-pass size delta     = "
                << moveFirstPassSizeDeltaSum_/max(moveFirstPassReuseChecks_, label(1)) << nl
                << "    avg first-pass (delta-appended)= "
                << moveFirstPassDeltaMinusPendingSum_/max(moveFirstPassReuseChecks_, label(1)) << nl
                << "    avg inter-step cloud delta    = "
                << moveInterStepCloudDeltaSum_/max(moveProfileCalls_ - 1, label(1)) << nl
                << "    inter-step cloud delta min/max= "
                << (moveInterStepCloudDeltaMin_ == labelMax ? 0 : moveInterStepCloudDeltaMin_)
                << " / "
                << (moveInterStepCloudDeltaMax_ == labelMin ? 0 : moveInterStepCloudDeltaMax_) << nl
                << "    avg post-step cloud delta     = "
                << movePostStepCloudDeltaSum_/max(moveProfileCalls_, label(1)) << nl
                << "    post-step cloud delta min/max = "
                << (movePostStepCloudDeltaMin_ == labelMax ? 0 : movePostStepCloudDeltaMin_)
                << " / "
                << (movePostStepCloudDeltaMax_ == labelMin ? 0 : movePostStepCloudDeltaMax_) << nl
                << "    avg deferred parcels/pass     = "
                << moveDeferredParcelsSum_/max(moveCommitCountCalls_, label(1)) << nl
                << "    max deferred parcels/pass     = "
                << moveDeferredParcelsMax_ << nl
                << "    avg received parcels/pass     = "
                << moveReceivedParcelsSum_/max(moveCommitCountCalls_, label(1)) << nl
                << "    max received parcels/pass     = "
                << moveReceivedParcelsMax_ << nl
                << "    avg move surviving parcels    = "
                << moveSurvivorParcelsSum_/max(moveCommitCountCalls_, label(1)) << nl
                << "    avg move transferred parcels  = "
                << moveTransferredParcelsSum_/max(moveCommitCountCalls_, label(1)) << nl
                << "    avg move deleted parcels      = "
                << moveDeletedParcelsSum_/max(moveCommitCountCalls_, label(1)) << nl
                << "    trackToAndHitFace total [s]   = "
                << moveTrackWallTime_ << nl
                << "    tracker callback total [s]    = "
                << moveTrackerCallbackWallTime_ << nl
                << "    boundary control total [s]    = "
                << moveBoundaryControlWallTime_ << nl
                << "    avg face hits/call            = "
                << moveFaceHitsSum_/max(moveProfileCalls_, label(1)) << nl
                << "    avg processor hits/call       = "
                << moveProcessorHitsSum_/max(moveProfileCalls_, label(1)) << nl
                << "    avg patch hits/call           = "
                << movePatchHitsSum_/max(moveProfileCalls_, label(1)) << nl
                << "    avg cyclic hits/call          = "
                << moveCyclicHitsSum_/max(moveProfileCalls_, label(1)) << nl
                << "    avg stuck hits/call           = "
                << moveStuckHitsSum_/max(moveProfileCalls_, label(1)) << nl
                << endl;
        }
    }

    if (moveProfileCalls_ > 0 && profilingDetailEnabled_)
    {
        reportMoveCellHotspots();
    }

    if (buildOccupancyProfileCalls_ > 0 && isOutputRank())
    {
        const scalar totalProfiled =
            buildOccupancyExtractWallTime_
          + buildOccupancyCountWallTime_
          + buildOccupancyAssembleWallTime_;

        Info<< "BuildCellOccupancy profiling summary:" << nl
            << "    buildCellOccupancy calls      = " << buildOccupancyProfileCalls_ << nl
            << "    extract parcels [s]           = " << buildOccupancyExtractWallTime_ << nl
            << "    count/reduce [s]              = " << buildOccupancyCountWallTime_ << nl
            << "    allocate/fill [s]             = " << buildOccupancyAssembleWallTime_ << nl
            << "    total profiled [s]            = " << totalProfiled << nl
            << endl;

        if (profilingDetailEnabled_)
        {
            Info<< "BuildCellOccupancy profiling detail:" << nl
                << "    move-ordered hits             = " << buildOccupancyMoveOrderedHits_ << nl
                << "    fallback gathers             = " << buildOccupancyFallbackHits_ << nl
                << "    fallback invalid-order       = "
                << buildOccupancyMoveOrderedInvalidHits_ << nl
                << "    fallback offset-mismatch     = "
                << buildOccupancyMoveOrderedOffsetMismatchHits_ << nl
                << "    fallback size-mismatch       = "
                << buildOccupancyMoveOrderedSizeMismatchHits_ << nl
                << "    avg move-ordered parcels      = "
                << buildOccupancyMoveOrderedParcelsSum_/max(buildOccupancyProfileCalls_, label(1))
                << nl
                << "    avg appended parcels          = "
                << buildOccupancyMoveAppendedParcelsSum_/max(buildOccupancyProfileCalls_, label(1))
                << nl
                << "    avg cloud size                = "
                << buildOccupancyCloudSizeSum_/max(buildOccupancyProfileCalls_, label(1))
                << nl
                << "    size delta avg/min/max        = "
                << buildOccupancyMoveOrderedSizeDeltaSum_/max(buildOccupancyProfileCalls_, label(1))
                << " / " << buildOccupancyMoveOrderedSizeDeltaMin_
                << " / " << buildOccupancyMoveOrderedSizeDeltaMax_ << nl
                << nl
                << endl;

            if (emitStepDiagnostics_)
            {
                Info<< "    avg cloud valid-cell parcels   = "
                    << buildOccupancyCloudValidCellSum_/max(buildOccupancyProfileCalls_, label(1))
                    << nl
                    << "    avg cloud invalid-cell parcels = "
                    << buildOccupancyCloudInvalidCellSum_/max(buildOccupancyProfileCalls_, label(1))
                    << nl
                    << "    avg move-order valid parcels   = "
                    << buildOccupancyMoveOrderedValidCellSum_/max(buildOccupancyProfileCalls_, label(1))
                    << nl
                    << "    avg move-order invalid parcels = "
                    << buildOccupancyMoveOrderedInvalidCellSum_/max(buildOccupancyProfileCalls_, label(1))
                    << nl
                    << "    avg partition assigned parcels = "
                    << buildOccupancyPartitionAssignedParcelsSum_/max(buildOccupancyProfileCalls_, label(1))
                    << nl
                    << "    avg partition gap cells        = "
                    << buildOccupancyPartitionGapCellsSum_/max(buildOccupancyProfileCalls_, label(1))
                    << nl
                    << "    avg partition overlap cells    = "
                    << buildOccupancyPartitionOverlapCellsSum_/max(buildOccupancyProfileCalls_, label(1))
                    << nl
                    << endl;
            }
        }
    }

    if (collisionProfileCalls_ > 0 && isOutputRank())
    {
        const scalar totalProfiled =
            collisionPrecomputeWallTime_
          + collisionPartitionWallTime_
          + collisionSelectionWallTime_;

        Info<< "Collision profiling summary:" << nl
            << "    collision calls               = " << collisionProfileCalls_ << nl
            << "    precompute candidates [s]     = " << collisionPrecomputeWallTime_ << nl
            << "    rebuild partition [s]         = " << collisionPartitionWallTime_ << nl
            << "    selection/collide [s]         = " << collisionSelectionWallTime_ << nl
            << "    total profiled [s]            = " << totalProfiled << nl
            << endl;
    }

    if (evolveProfileCalls_ > 0 && isOutputRank())
    {
        const scalar totalProfiled =
            evolveMoveWallTime_
          + evolveBuildWallTime_
          + evolveCoordWallTime_
          + evolveCollisionWallTime_
          + evolveReactionWallTime_
          + evolvePostWallTime_;
        const scalar fullAccounted =
            totalProfiled
          + evolvePreWallTime_
          + evolveRankTimeWallTime_
          + evolveAutoRebalanceWallTime_
          + evolvePostProfileResidualWallTime_;

        Info<< "Evolve profiling summary:" << nl
            << "    evolve calls                  = " << evolveProfileCalls_ << nl
            << "    pre/update setup [s]          = " << evolvePreWallTime_ << nl
            << "    move only [s]                 = " << evolveMoveWallTime_ << nl
            << "    buildCellOccupancy [s]        = " << evolveBuildWallTime_ << nl
            << "    coordSystem [s]               = " << evolveCoordWallTime_ << nl
            << "    collision phase [s]           = " << evolveCollisionWallTime_ << nl
            << "    reaction/output [s]           = " << evolveReactionWallTime_ << nl
            << "    post fields/output [s]        = " << evolvePostWallTime_ << nl
            << "    rank-time bookkeeping [s]     = " << evolveRankTimeWallTime_ << nl
            << "    auto DLB/rebalance [s]        = " << evolveAutoRebalanceWallTime_ << nl
            << "    post-profile residual [s]     = " << evolvePostProfileResidualWallTime_ << nl
            << "    total profiled [s]            = " << totalProfiled << nl
            << "    full evolve accounted [s]     = " << fullAccounted << nl
            << "    full evolve wall [s]          = " << evolveFullWallTime_ << nl
            << "    evolve residual [s]           = " << evolveFullWallTime_ - fullAccounted << nl
            << endl;
    }

    const label localMoveParcels = sum(moveThreadParticleCounts_);
    const label localCollisionCandidates = sum(collisionThreadCandidateCounts_);
    const label localAcceptedCollisions = sum(collisionThreadAcceptedCounts_);

    if
    (
        profilingDetailEnabled_
     && Pstream::parRun()
     && (localMoveParcels || localCollisionCandidates || localAcceptedCollisions)
    )
    {
        Pout<< "Load stats rank " << Pstream::myProcNo() << ":" << nl
            << "    move particles processed      = " << localMoveParcels << nl
            << "    collision candidates          = " << localCollisionCandidates << nl
            << "    accepted collisions           = " << localAcceptedCollisions << nl
            << endl;
    }

    if
    (
        profilingDetailEnabled_
     && (localMoveParcels || localCollisionCandidates || localAcceptedCollisions)
    )
    {
        scalar moveMax = localMoveParcels;
        scalar moveMin = localMoveParcels;
        scalar candMax = localCollisionCandidates;
        scalar candMin = localCollisionCandidates;
        scalar collMax = localAcceptedCollisions;
        scalar collMin = localAcceptedCollisions;

        if (Pstream::parRun())
        {
            reduce(moveMax, maxOp<scalar>());
            reduce(moveMin, minOp<scalar>());
            reduce(candMax, maxOp<scalar>());
            reduce(candMin, minOp<scalar>());
            reduce(collMax, maxOp<scalar>());
            reduce(collMin, minOp<scalar>());
        }
        else if (moveThreadParticleCounts_.size())
        {
            moveMax = scalar(max(moveThreadParticleCounts_));
            moveMin = scalar(min(moveThreadParticleCounts_));
            candMax = scalar(max(collisionThreadCandidateCounts_));
            candMin = scalar(min(collisionThreadCandidateCounts_));
            collMax = scalar(max(collisionThreadAcceptedCounts_));
            collMin = scalar(min(collisionThreadAcceptedCounts_));
        }

        if (isOutputRank())
        {
            Info<< "Profiling detail load-balance summary:" << nl
                << "    move particles max/min        = "
                << moveMax << " / " << moveMin << nl
                << "    move imbalance max/min        = "
                << (moveMin > SMALL ? moveMax/moveMin : 0.0) << nl
                << "    collision cand max/min        = "
                << candMax << " / " << candMin << nl
                << "    collision cand imbalance      = "
                << (candMin > SMALL ? candMax/candMin : 0.0) << nl
                << "    accepted coll max/min         = "
                << collMax << " / " << collMin << nl
                << "    accepted coll imbalance       = "
                << (collMin > SMALL ? collMax/collMin : 0.0) << nl;

            if (!Pstream::parRun() && moveThreadParticleCounts_.size())
            {
                Info<< "    thread move particles         = " << moveThreadParticleCounts_ << nl;

                if (moveThreadWallTimes_.size() == moveThreadParticleCounts_.size())
                {
                    scalarField moveNsPerParticle(moveThreadParticleCounts_.size(), 0.0);

                    forAll(moveThreadParticleCounts_, i)
                    {
                        if (moveThreadParticleCounts_[i] > 0)
                        {
                            moveNsPerParticle[i] =
                                1.0e9*moveThreadWallTimes_[i]
                               /scalar(moveThreadParticleCounts_[i]);
                        }
                    }

                    Info<< "    thread move wall time [s]     = " << moveThreadWallTimes_ << nl
                        << "    thread move ns/particle       = " << moveNsPerParticle << nl;
                }

                if
                (
                    moveThreadTrackWallTimes_.size() == moveThreadParticleCounts_.size()
                )
                {
                    scalarField trackNsPerParticle(moveThreadParticleCounts_.size(), 0.0);
                    scalarField boundaryNsPerParticle(moveThreadParticleCounts_.size(), 0.0);

                    forAll(moveThreadParticleCounts_, i)
                    {
                        if (moveThreadParticleCounts_[i] > 0)
                        {
                            trackNsPerParticle[i] =
                                1.0e9*moveThreadTrackWallTimes_[i]
                               /scalar(moveThreadParticleCounts_[i]);
                            boundaryNsPerParticle[i] =
                                1.0e9*moveThreadBoundaryWallTimes_[i]
                               /scalar(moveThreadParticleCounts_[i]);
                        }
                    }

                    Info<< "    thread track wall time [s]    = " << moveThreadTrackWallTimes_ << nl
                        << "    thread tracker wall time [s]  = " << moveThreadTrackerWallTimes_ << nl
                        << "    thread boundary wall time [s] = " << moveThreadBoundaryWallTimes_ << nl
                        << "    thread face hits              = " << moveThreadFaceHitCounts_ << nl
                        << "    thread processor hits         = " << moveThreadProcessorHitCounts_ << nl
                        << "    thread patch hits             = " << moveThreadPatchHitCounts_ << nl
                        << "    thread cyclic hits            = " << moveThreadCyclicHitCounts_ << nl
                        << "    thread stuck hits             = " << moveThreadStuckHitCounts_ << nl
                        << "    thread track ns/particle      = " << trackNsPerParticle << nl
                        << "    thread boundary ns/particle   = " << boundaryNsPerParticle << nl;
                }

                Info<< "    thread collision candidates   = " << collisionThreadCandidateCounts_ << nl
                    << "    thread accepted collisions    = " << collisionThreadAcceptedCounts_ << nl;

                if
                (
                    collisionThreadActiveCellCounts_.size() == collisionThreadCandidateCounts_.size()
                 && collisionThreadReactionHitCounts_.size() == collisionThreadCandidateCounts_.size()
                 && collisionThreadWallTimes_.size() == collisionThreadCandidateCounts_.size()
                )
                {
                    scalarField nsPerCandidate(collisionThreadCandidateCounts_.size(), 0.0);
                    scalarField nsPerAccepted(collisionThreadAcceptedCounts_.size(), 0.0);

                    forAll(collisionThreadCandidateCounts_, i)
                    {
                        if (collisionThreadCandidateCounts_[i] > 0)
                        {
                            nsPerCandidate[i] =
                                1.0e9*collisionThreadWallTimes_[i]
                               /scalar(collisionThreadCandidateCounts_[i]);
                        }

                        if (collisionThreadAcceptedCounts_[i] > 0)
                        {
                            nsPerAccepted[i] =
                                1.0e9*collisionThreadWallTimes_[i]
                               /scalar(collisionThreadAcceptedCounts_[i]);
                        }
                    }

                    Info<< "    thread collision active cells = " << collisionThreadActiveCellCounts_ << nl
                        << "    thread collision reaction hits= " << collisionThreadReactionHitCounts_ << nl
                        << "    thread collision wall time [s]= " << collisionThreadWallTimes_ << nl
                        << "    thread collision ns/candidate = " << nsPerCandidate << nl
                        << "    thread collision ns/accepted  = " << nsPerAccepted << nl;
                }
            }

            Info<< endl;
        }

        const_cast<dsmcCloud&>(*this).resetLoadStats();
    }

    if (replicatedMeshActive() && !Pstream::parRun())
    {
        const label myRank = replicatedMesh_->myRank();
        const label nProcs = replicatedMesh_->nProcs();
        const label localParcels = this->size();

        labelList allParcels(nProcs, 0);
        labelList allCollisions(nProcs, 0);
        labelList allCandidates(nProcs, 0);
        allParcels[myRank] = localParcels;
        allCollisions[myRank] = cumulativeCollisions_;
        allCandidates[myRank] = cumulativeCollisionCandidates_;

        MPI_Allreduce(MPI_IN_PLACE, allParcels.data(), nProcs, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, allCollisions.data(), nProcs, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, allCandidates.data(), nProcs, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (myRank == 0)
        {
            const label totalParcels = sum(allParcels);
            const label totalCollisions = sum(allCollisions);
            const label totalCandidates = sum(allCandidates);

            Info<< "Replicated mesh final summary:" << nl
                << "    total dsmc particles            = " << totalParcels << nl
                << "    last-step collisions            = " << totalCollisions << nl
                << "    last-step candidates            = " << totalCandidates << nl
                << "    last-step acceptance rate       = "
                << (totalCandidates > 0 ? scalar(totalCollisions)/scalar(totalCandidates) : 0)
                << nl;
            for (label r = 0; r < nProcs; ++r)
            {
                Info<< "    rank " << r << ": particles=" << allParcels[r]
                    << " collisions=" << allCollisions[r]
                    << " candidates=" << allCandidates[r] << nl;
            }
            Info<< endl;
        }
    }

    if (replicatedMeshActive() && evolveProfileEnabled_)
    {
        replicatedMesh_->report();
    }
}

// ************************************************************************* //
