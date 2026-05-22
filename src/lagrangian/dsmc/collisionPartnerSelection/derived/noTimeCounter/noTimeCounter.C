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

Class
    noTimeCounter

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounter.H"
#include "addToRunTimeSelectionTable.H"
#include <chrono>
#include <cstdint>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounter, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void noTimeCounter::readControlDictParams()
{
    const dictionary& cd = cloud_.mesh().time().controlDict();

    collisionFastRng_ =
        cd.lookupOrDefault<bool>("collisionFastRng", false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noTimeCounter::noTimeCounter
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    infoCounter_(0),
    threadWhichSubCell_(),
    threadSubCells_(),
    threadParcelPtrs_(),
    threadVelocities_(),
    threadTypeIds_(),
    threadCharges_(),
    collisionFastRng_(false)
{
    readControlDictParams();
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounter::~noTimeCounter()
{}

// Fast XorShiro128+ RNG for collision hot path (2-5x faster than Rand48)
struct FastRng
{
    uint64_t s0, s1;
    FastRng(uint64_t seed = 1)
    {
        // SplitMix64 seeding
        uint64_t z = seed + 0x9e3779b97f4a7c15ULL;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        s0 = z ^ (z >> 31);
        z = seed + 0x9e3779b97f4a7c15ULL + 1;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        s1 = z ^ (z >> 31);
    }
    inline uint64_t next()
    {
        const uint64_t s0_ = s0;
        uint64_t s1_ = s1;
        s1_ ^= s1_ << 23;
        s1 = s1_ ^ s0_ ^ (s1_ >> 18) ^ (s0_ >> 5);
        s0 = s1_;
        return s1 + s0_;
    }
    inline double sample01() { return (next() >> 11) * 0x1.0p-53; }
    inline label position(label n) { return label(next() % uint64_t(n)); }
};

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounter::initialConfiguration()
{

}

void noTimeCounter::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    const label statsThreads =
        cloud_.openmpEnabled() ? max(cloud_.ompNumThreads(), label(1)) : label(1);

    // Per-thread fast RNGs seeded from cloud master RNG
    List<FastRng> threadFastRng(statsThreads);
    for (label threadI = 0; threadI < statsThreads; ++threadI)
    {
        const label myProc = Pstream::myProcNo();
        threadFastRng[threadI] = FastRng
        (
            uint64_t(myProc)*1000000ULL + uint64_t(threadI)*10000ULL
          + uint64_t(cloud_.mesh().time().timeIndex())
        );
    }

    labelList threadCandidateCounts(statsThreads, 0);
    labelList threadAcceptedCounts(statsThreads, 0);
    labelList threadActiveCellCounts(statsThreads, 0);
    labelList threadReactionHitCounts(statsThreads, 0);
    scalarField threadWallTimes(statsThreads, 0.0);

    const bool useFlatOccupancy = cloud_.hasOccupancyOrderedParcels();
    const List<DynamicList<dsmcParcel*>>* cellOccupancyPtr =
        useFlatOccupancy ? nullptr : &cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();
    const label nCells = mesh.nCells();
    const label collisionChunk = max(cloud_.openmpCollisionChunk(), label(1));
    const bool collisionFastRng = collisionFastRng_;

    if (threadWhichSubCell_.size() != statsThreads)
    {
        threadWhichSubCell_.setSize(statsThreads);
        threadSubCells_.setSize(statsThreads);
        threadParcelPtrs_.setSize(statsThreads);
        threadVelocities_.setSize(statsThreads);
        threadTypeIds_.setSize(statsThreads);
        threadCharges_.setSize(statsThreads);
    }

    for (label threadI = 0; threadI < statsThreads; ++threadI)
    {
        if (threadSubCells_[threadI].size() != 8)
        {
            threadSubCells_[threadI].setSize(8);
        }
    }

    auto processCell =
    [&]
    (
        const label cellI,
        const label threadI,
        label& activeCellCount,
        label& reactionHitCount
    )
    -> label
    {
        FastRng& frng = threadFastRng[threadI];

        auto rngPos = [&](label n) -> label
        {
            return collisionFastRng
                 ? frng.position(n)
                 : cloud_.randomLabel(0, n - 1);
        };
        auto rng01 = [&]() -> scalar
        {
            return collisionFastRng
                 ? frng.sample01()
                 : cloud_.rndGen().sample01<scalar>();
        };

        const DynamicList<dsmcParcel*>* cellParcelsPtr =
            useFlatOccupancy ? nullptr : &(*cellOccupancyPtr)[cellI];
        const label occStart = useFlatOccupancy ? cloud_.occupancyStart(cellI) : 0;
        const label nC =
            useFlatOccupancy ? cloud_.occupancyCount(cellI) : cellParcelsPtr->size();
        const label nCandidates = cloud_.nCandidatesPerCell()[cellI];
        label acceptedCollisions = 0;

        if (nC > 1 && nCandidates > 0)
        {
            ++activeCellCount;
            DynamicList<label>& whichSubCell = threadWhichSubCell_[threadI];
            List<DynamicList<label>>& subCells = threadSubCells_[threadI];
            DynamicList<dsmcParcel*>& parcelPtrs = threadParcelPtrs_[threadI];
            DynamicList<vector>& velocities = threadVelocities_[threadI];
            DynamicList<label>& typeIds = threadTypeIds_[threadI];
            DynamicList<label>& charges = threadCharges_[threadI];

            const point& cC = mesh.cellCentres()[cellI];
            label subCellCounts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            label subCellOffsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};

            whichSubCell.setSize(nC);
            parcelPtrs.setSize(nC);
            velocities.setSize(nC);
            typeIds.setSize(nC);
            charges.setSize(nC);

            for (label i = 0; i < nC; ++i)
            {
                dsmcParcel* pPtr =
                    useFlatOccupancy
                  ? cloud_.occupancyParcel(occStart + i)
                  : (*cellParcelsPtr)[i];
                const label typeId = pPtr->typeId();
                const vector relPos = pPtr->position() - cC;
                const label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                parcelPtrs[i] = pPtr;
                velocities[i] = pPtr->U();
                typeIds[i] = typeId;
                charges[i] = cloud_.constProps(typeId).charge();
                whichSubCell[i] = subCell;
                ++subCellCounts[subCell];
            }

            for (label subCellI = 0; subCellI < 8; ++subCellI)
            {
                subCells[subCellI].setSize(subCellCounts[subCellI]);
            }

            for (label i = 0; i < nC; ++i)
            {
                const label subCell = whichSubCell[i];
                subCells[subCell][subCellOffsets[subCell]++] = i;
            }

            scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];

            for (label c = 0; c < nCandidates; c++)
            {
                label candidateP = rngPos(nC);
                label candidateQ = -1;

                const List<label>& subCellPs = subCells[whichSubCell[candidateP]];
                const label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    do
                    {
                        candidateQ = subCellPs[rngPos(nSC)];
                    } while (candidateP == candidateQ);
                }
                else
                {
                    do
                    {
                        candidateQ = rngPos(nC);
                    } while (candidateP == candidateQ);
                }

                const label typeIdP = typeIds[candidateP];
                const label typeIdQ = typeIds[candidateQ];
                const label chargeP = charges[candidateP];
                const label chargeQ = charges[candidateQ];

                if(!(chargeP == -1 && chargeQ == -1))
                {

                    const scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        *parcelPtrs[candidateP],
                        *parcelPtrs[candidateQ]
                    );

                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    if ((sigmaTcR/sigmaTcRMax) > rng01())
                    {
                        const label rMId =
                            cloud_.reactionsActive()
                          ? cloud_.reactions().pairModelAddressing()[typeIdP][typeIdQ]
                          : -1;

                        dsmcParcel& parcelP = *parcelPtrs[candidateP];
                        dsmcParcel& parcelQ = *parcelPtrs[candidateQ];

                        if(rMId != -1)
                        {
                            ++reactionHitCount;
                            cloud_.reactions().reactions()[rMId]->reaction
                            (
                                parcelP,
                                parcelQ
                            );
                            if(cloud_.reactions().reactions()[rMId]->relax())
                            {
                                cloud_.binaryCollision().collide
                                (
                                    parcelP,
                                    parcelQ,
                                    cellI
                                );
                            }
                        }
                        else
                        {
                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }

                        acceptedCollisions++;

                        velocities[candidateP] = parcelP.U();
                        velocities[candidateQ] = parcelQ.U();
                        typeIds[candidateP] = parcelP.typeId();
                        typeIds[candidateQ] = parcelQ.typeId();
                        charges[candidateP] =
                            cloud_.constProps(typeIds[candidateP]).charge();
                        charges[candidateQ] =
                            cloud_.constProps(typeIds[candidateQ]).charge();
                    }
                }
            }
        }
        return acceptedCollisions;
    };

    #ifdef _OPENMP
    if (cloud_.openmpEnabled())
    {
        if (cloud_.openmpCollisionSchedule() == "partition")
        {
            #pragma omp parallel num_threads(statsThreads)
            {
                using clock_type = std::chrono::steady_clock;
                const auto tBegin = clock_type::now();
                const label threadI = cloud_.currentThreadId();
                const label startCell = cloud_.collisionLoadStart()[threadI];
                const label endCell = cloud_.collisionLoadEnd()[threadI];
                label localCandidates = 0;
                label localCollisions = 0;
                label localActiveCells = 0;
                label localReactionHits = 0;

                for (label cellI = startCell; cellI < endCell; ++cellI)
                {
                    localCandidates += cloud_.nCandidatesPerCell()[cellI];
                    localCollisions += processCell
                    (
                        cellI,
                        threadI,
                        localActiveCells,
                        localReactionHits
                    );
                }

                threadCandidateCounts[threadI] = localCandidates;
                threadAcceptedCounts[threadI] = localCollisions;
                threadActiveCellCounts[threadI] = localActiveCells;
                threadReactionHitCounts[threadI] = localReactionHits;
                threadWallTimes[threadI] =
                    std::chrono::duration<scalar>(clock_type::now() - tBegin).count();
            }
        }
        else
        {
            #pragma omp parallel num_threads(statsThreads)
            {
                using clock_type = std::chrono::steady_clock;
                const auto tBegin = clock_type::now();
                const label threadI = cloud_.currentThreadId();
                label localCandidates = 0;
                label localCollisions = 0;
                label localActiveCells = 0;
                label localReactionHits = 0;

                #pragma omp for schedule(dynamic, collisionChunk)
                for (label cellI = 0; cellI < nCells; ++cellI)
                {
                    if (cloud_.nCandidatesPerCell()[cellI] == 0) continue;
                    localCandidates += cloud_.nCandidatesPerCell()[cellI];
                    localCollisions += processCell
                    (
                        cellI,
                        threadI,
                        localActiveCells,
                        localReactionHits
                    );
                }

                threadCandidateCounts[threadI] = localCandidates;
                threadAcceptedCounts[threadI] = localCollisions;
                threadActiveCellCounts[threadI] = localActiveCells;
                threadReactionHitCounts[threadI] = localReactionHits;
                threadWallTimes[threadI] =
                    std::chrono::duration<scalar>(clock_type::now() - tBegin).count();
            }
        }
    }
    else
    #endif
    {
        using clock_type = std::chrono::steady_clock;
        const auto tBegin = clock_type::now();
        label localActiveCells = 0;
        label localReactionHits = 0;

        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            if (cloud_.nCandidatesPerCell()[cellI] == 0) continue;
            threadCandidateCounts[0] += cloud_.nCandidatesPerCell()[cellI];
            threadAcceptedCounts[0] += processCell
            (
                cellI,
                0,
                localActiveCells,
                localReactionHits
            );
        }

        threadActiveCellCounts[0] = localActiveCells;
        threadReactionHitCounts[0] = localReactionHits;
        threadWallTimes[0] =
            std::chrono::duration<scalar>(clock_type::now() - tBegin).count();
    }

    cloud_.recordCollisionThreadProfile
    (
        threadCandidateCounts,
        threadAcceptedCounts,
        threadActiveCellCounts,
        threadReactionHitCounts,
        threadWallTimes
    );

    label collisionCandidates = sum(threadCandidateCounts);
    label collisions = sum(threadAcceptedCounts);
    label localCandidateCells = 0;
    scalar localSigmaSum = 0.0;
    scalar localSigmaMax = -GREAT;
    scalar localSigmaMin = GREAT;

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        const label nc = cloud_.nCandidatesPerCell()[cellI];
        if (nc > 0)
        {
            const scalar sigma = cloud_.sigmaTcRMax()[cellI];
            ++localCandidateCells;
            localSigmaSum += sigma;
            localSigmaMax = max(localSigmaMax, sigma);
            localSigmaMin = min(localSigmaMin, sigma);
        }
    }

    reduce(collisions, sumOp<label>());
    reduce(collisionCandidates, sumOp<label>());

    cloud_.accumulateCollisionCounts(collisionCandidates, collisions);

    label globalCandidateCells = localCandidateCells;
    scalar globalSigmaSum = localSigmaSum;
    scalar globalSigmaMax = localCandidateCells > 0 ? localSigmaMax : scalar(-GREAT);
    scalar globalSigmaMin = localCandidateCells > 0 ? localSigmaMin : scalar(GREAT);

    reduce(globalCandidateCells, sumOp<label>());
    reduce(globalSigmaSum, sumOp<scalar>());
    reduce(globalSigmaMax, maxOp<scalar>());
    reduce(globalSigmaMin, minOp<scalar>());

    const scalar globalAcceptanceRate =
        collisionCandidates > 0
      ? scalar(collisions)/scalar(collisionCandidates)
      : scalar(0);
    const scalar globalSigmaAvg =
        globalCandidateCells > 0
      ? globalSigmaSum/scalar(globalCandidateCells)
      : scalar(0);

    cloud_.sigmaTcRMax().correctBoundaryConditions();

    infoCounter_++;

    if(infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (cloud_.replicatedMeshActive() && !Pstream::parRun())
        {
            const label myRank = cloud_.replicatedMesh().myRank();
            if (collisionCandidates)
            {
                Info<< "    Collisions [rank " << myRank << "]"
                    << "              = " << collisions << nl
                    << "    Collision candidates           = "
                    << collisionCandidates << nl
                    << "    Collision acceptance rate      = "
                    << globalAcceptanceRate << nl
                    << endl;
            }
            else
            {
                Info<< "    No collisions [rank " << myRank << "]" << endl;
            }
        }
        else if (cloud_.isOutputRank())
        {
            if (collisionCandidates)
            {
                Info<< "    Collisions                      = "
                    << collisions << nl
                    << "    Collision candidates           = "
                    << collisionCandidates << nl
                    << "    Collision acceptance rate      = "
                    << globalAcceptanceRate << nl
                    << "    Candidate-cell sigmaTcRMax avg/max/min = "
                    << globalSigmaAvg << " / "
                    << (globalCandidateCells > 0 ? globalSigmaMax : scalar(0)) << " / "
                    << (globalCandidateCells > 0 ? globalSigmaMin : scalar(0)) << nl
                    << endl;
            }
            else
            {
                Info<< "    No collisions" << endl;
            }
        }

        infoCounter_ = 0;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
