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
#include <cstdio>
#include "addToRunTimeSelectionTable.H"
#include "../../../dsmcFastRng.H"
#include <chrono>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounter, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace
{

typedef std::chrono::steady_clock CollisionClock;

inline scalar elapsedCollisionSeconds
(
    const CollisionClock::time_point& start
)
{
    return scalar
    (
        std::chrono::duration<double>(CollisionClock::now() - start).count()
    );
}

}


void noTimeCounter::readControlDictParams()
{
    const dictionary& controlDict = cloud_.mesh().time().controlDict();

    fastRng_ = controlDict.lookupOrDefault<bool>("fastRng", false);

    const bool replicatedMpiDefault =
        controlDict.lookupOrDefault<bool>("replicatedMesh", false);
    collisionReduceOnlyOnOutput_ =
        controlDict.lookupOrDefault<bool>
        (
            "collisionReduceOnlyOnOutput",
            replicatedMpiDefault
        );
    collisionOutputGlobalReduce_ =
        controlDict.lookupOrDefault<bool>
        (
            "collisionOutputGlobalReduce",
            !replicatedMpiDefault
        );

    if (cloud_.isOutputRank() && controlDict.found("collisionReduceOnlyOnOutput"))
    {
        Info<< "Collision profile: collisionReduceOnlyOnOutput explicitly set to "
            << collisionReduceOnlyOnOutput_ << endl;
    }
    else if (cloud_.isOutputRank() && replicatedMpiDefault)
    {
        Info<< "Collision profile: collisionReduceOnlyOnOutput defaulted to true"
            << " for replicatedMesh run" << endl;
    }

    if (cloud_.isOutputRank() && controlDict.found("collisionOutputGlobalReduce"))
    {
        Info<< "Collision profile: collisionOutputGlobalReduce explicitly set to "
            << collisionOutputGlobalReduce_ << endl;
    }
    else if (cloud_.isOutputRank() && replicatedMpiDefault)
    {
        Info<< "Collision profile: collisionOutputGlobalReduce defaulted to false"
            << " for replicatedMesh run" << endl;
    }
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
    threadTypeIds_(),
    threadCharges_(),
    threadCandidateCells_(),
    collisionThreadFastRng_(),
    candidateCellsToClear_(),
    profileLocalCollisions_(0),
    profileLocalCollisionCandidates_(0),
    profileCollisionLocalLoopWall_(0.0),
    profileCollisionReduceWall_(0.0),
    profileCollisionSigmaWall_(0.0),
    profileCollisionPreOutputWall_(0.0),
    fastRng_(false),
    collisionReduceOnlyOnOutput_(false),
    collisionOutputGlobalReduce_(true)
//     propsDict_(dict.subDict(typeName + "Properties"))
{
    readControlDictParams();
}


void noTimeCounter::clearCandidateCounts
(
    labelList& nCandidatesPerCell,
    const label nCells
)
{
    if (nCandidatesPerCell.size() != nCells)
    {
        nCandidatesPerCell.setSize(nCells, 0);
        candidateCellsToClear_.clear();
        return;
    }

    forAll(candidateCellsToClear_, i)
    {
        const label cellI = candidateCellsToClear_[i];
        if (cellI >= 0 && cellI < nCells)
        {
            nCandidatesPerCell[cellI] = 0;
        }
    }

    candidateCellsToClear_.clear();
}


void noTimeCounter::ensureCollisionThreadFastRngs(const label nThreads)
{
    if (collisionThreadFastRng_.size() >= nThreads)
    {
        return;
    }

    const label oldSize = collisionThreadFastRng_.size();
    collisionThreadFastRng_.setSize(nThreads, nullptr);

    for (label threadI = oldSize; threadI < nThreads; ++threadI)
    {
        collisionThreadFastRng_[threadI] = new dsmcPaddedFastRng();
    }
}


void noTimeCounter::printReplicatedRankCollisionDetail
(
    const label localCollisions,
    const label localCollisionCandidates,
    const label globalCollisions,
    const label globalCollisionCandidates
) const
{
    if
    (
        !cloud_.emitStepDiagnostics()
     || !cloud_.replicatedMeshActive()
     || cloud_.replicatedMesh().nProcs() <= 1
    )
    {
        return;
    }

    int mpiInit = 0;
    MPI_Initialized(&mpiInit);
    if (!mpiInit)
    {
        return;
    }

    const label myRank = cloud_.replicatedMesh().myRank();
    const label nRanks = cloud_.replicatedMesh().nProcs();
    const label nLocalValues = 5;
    const label nScalarValues = 8;

    // clean M13: cumulative counters can exceed int32 at scale (this build
    // has WM_LABEL_SIZE=32); gather them as doubles (exact for integers
    // < 2^53) instead of MPI_INT.
    //
    // Verification hook (default 0): add a synthetic offset to the four
    // local counts to exercise this aggregation path near/across 2^31
    // without running billions of collisions.  The offset is applied to
    // the GATHERED copy only — the real accumulators are untouched.
    // Read as a double: a 32-bit label cannot even express 2^31, which is
    // exactly the boundary this verification targets.
    const double syntheticOffset =
        cloud_.mesh().time().controlDict().lookupOrDefault<double>
        (
            "collisionDetailSyntheticOffset",
            0.0
        );
    const double synthetic = syntheticOffset + double(myRank);

    double localValues[nLocalValues] =
    {
        double(localCollisions) + synthetic,
        double(localCollisionCandidates) + synthetic,
        profileLocalCollisions_ + synthetic,
        profileLocalCollisionCandidates_ + synthetic,
        double(cloud_.size())
    };

    scalar localScalarValues[nScalarValues] =
    {
        profileCollisionLocalLoopWall_,
        profileCollisionReduceWall_,
        profileCollisionSigmaWall_,
        profileCollisionPreOutputWall_,
        cloud_.moveDetailCumulative(0),
        cloud_.moveDetailCumulative(1)
    };

    List<scalar> allValues;
    List<scalar> allScalarValues;
    if (myRank == 0)
    {
        allValues.setSize(nRanks*nLocalValues, 0);
        allScalarValues.setSize(nRanks*nScalarValues, 0.0);
    }

    MPI_Gather
    (
        localValues,
        nLocalValues,
        MPI_DOUBLE,
        myRank == 0 ? allValues.data() : nullptr,
        nLocalValues,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );

    MPI_Gather
    (
        localScalarValues,
        nScalarValues,
        MPI_DOUBLE,
        myRank == 0 ? allScalarValues.data() : nullptr,
        nScalarValues,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD
    );

    if (myRank != 0)
    {
        return;
    }

    // clean M13: the sentinel must not be labelMax (and not any 32-bit
    // limit) — a synthetic or real count above 2^31 would never become the
    // minimum and the max/min ratio would be wrong (observed 2.0 instead
    // of ~1.0 with every count above 2^31).
    double minLocalCandidates = 1.0e300;
    double maxLocalCandidates = 0.0;
    double minCumulativeCandidates = 1.0e300;
    double maxCumulativeCandidates = 0.0;

    Info<< "    Replicated mesh collision detail by rank:" << nl
        << "        rank localColl localCand cumColl cumCand particles moveParcels trackCalls"
        << " localAcc cumAcc" << nl;

    for (label rankI = 0; rankI < nRanks; ++rankI)
    {
        const label base = rankI*nLocalValues;
        // clean M13: the gathered counters are double; convert to int64 for
        // printing (int32 label would truncate at 2^31 and hide overflow).
        const double rankLocalCollisions = allValues[base + 0];
        const double rankLocalCandidates = allValues[base + 1];
        const double rankCumulativeCollisions = allValues[base + 2];
        const double rankCumulativeCandidates = allValues[base + 3];
        const double rankParticles = allValues[base + 4];
        // moveDetail counters ride the scalar gather (int32 overflows at
        // 1200w scale within ~10 steps)
        const label sbase = rankI*nScalarValues;
        const double rankMoveParcels = allScalarValues[sbase + 4];
        const double rankTrackCalls = allScalarValues[sbase + 5];

        minLocalCandidates = min(minLocalCandidates, rankLocalCandidates);
        maxLocalCandidates = max(maxLocalCandidates, rankLocalCandidates);
        minCumulativeCandidates =
            min(minCumulativeCandidates, rankCumulativeCandidates);
        maxCumulativeCandidates =
            max(maxCumulativeCandidates, rankCumulativeCandidates);

        // Print the wide counters without going through label: in this
        // build label is int32, and OF's int64 stream operator writes
        // label(value), so the value would wrap at 2^31 (observed as a
        // negative count).  "%.0f" keeps every integer up to 2^53 exact.
        char collBuf[32];
        char candBuf[32];
        char cumCollBuf[32];
        char cumCandBuf[32];
        char partBuf[32];
        std::snprintf(collBuf, sizeof(collBuf), "%.0f", rankLocalCollisions);
        std::snprintf(candBuf, sizeof(candBuf), "%.0f", rankLocalCandidates);
        std::snprintf(cumCollBuf, sizeof(cumCollBuf), "%.0f",
                      rankCumulativeCollisions);
        std::snprintf(cumCandBuf, sizeof(cumCandBuf), "%.0f",
                      rankCumulativeCandidates);
        std::snprintf(partBuf, sizeof(partBuf), "%.0f", rankParticles);

        Info<< "        rank" << rankI
            << " " << collBuf
            << " " << candBuf
            << " " << cumCollBuf
            << " " << cumCandBuf
            << " " << partBuf
            << " " << rankMoveParcels
            << " " << rankTrackCalls
            << " " << rankLocalCollisions
                /max(rankLocalCandidates, SMALL)
            << " " << rankCumulativeCollisions
                /max(rankCumulativeCandidates, SMALL)
            << nl;
    }

    Info<< "        local candidates max/min      = "
        << scalar(maxLocalCandidates)
            /max(scalar(minLocalCandidates), SMALL) << nl
        << "        cumulative candidates max/min = "
        << scalar(maxCumulativeCandidates)
            /max(scalar(minCumulativeCandidates), SMALL) << nl
        << "        global collisions/candidates  = "
        << globalCollisions << " / " << globalCollisionCandidates << nl;

    scalar maxLocalLoopWall = 0.0;
    scalar maxReduceWall = 0.0;
    scalar maxSigmaWall = 0.0;
    scalar maxPreOutputWall = 0.0;

    Info<< "    Replicated mesh collision subphase wall by rank [s]:" << nl
        << "        rank localLoop reduce sigmaBC preOutputTotal accountedFrac" << nl;

    for (label rankI = 0; rankI < nRanks; ++rankI)
    {
        const label base = rankI*nScalarValues;
        const scalar localLoopWall = allScalarValues[base + 0];
        const scalar reduceWall = allScalarValues[base + 1];
        const scalar sigmaWall = allScalarValues[base + 2];
        const scalar preOutputWall = allScalarValues[base + 3];
        const scalar accountedWall =
            localLoopWall + reduceWall + sigmaWall;

        maxLocalLoopWall = max(maxLocalLoopWall, localLoopWall);
        maxReduceWall = max(maxReduceWall, reduceWall);
        maxSigmaWall = max(maxSigmaWall, sigmaWall);
        maxPreOutputWall = max(maxPreOutputWall, preOutputWall);

        Info<< "        rank" << rankI
            << " " << localLoopWall
            << " " << reduceWall
            << " " << sigmaWall
            << " " << preOutputWall
            << " " << accountedWall/max(preOutputWall, SMALL)
            << nl;
    }

    Info<< "        max localLoop [s]       = " << maxLocalLoopWall << nl
        << "        max reduce [s]          = " << maxReduceWall << nl
        << "        max sigmaBC [s]         = " << maxSigmaWall << nl
        << "        max preOutputTotal [s]  = " << maxPreOutputWall << nl;
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounter::~noTimeCounter()
{
    forAll(collisionThreadFastRng_, threadI)
    {
        delete collisionThreadFastRng_[threadI];
        collisionThreadFastRng_[threadI] = nullptr;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounter::initialConfiguration()
{

}

void noTimeCounter::collide()
{
    const bool subphaseTimers = cloud_.profileSummaryEnabled();
    const CollisionClock::time_point collideWallStart =
        subphaseTimers ? CollisionClock::now() : CollisionClock::time_point();

    if (!cloud_.binaryCollision().active())
    {
        labelList& nCandidatesPerCell = cloud_.nCandidatesPerCell();
        cloud_.accumulateCellCostCandidates(nCandidatesPerCell, candidateCellsToClear_);
        clearCandidateCounts(nCandidatesPerCell, cloud_.mesh().nCells());
        if (subphaseTimers)
        {
            profileCollisionPreOutputWall_ +=
                elapsedCollisionSeconds(collideWallStart);
        }
        return;
    }

    #ifdef _OPENMP
    if (cloud_.openmpEnabled())
    {
        const CollisionClock::time_point localLoopStart =
            subphaseTimers ? CollisionClock::now() : CollisionClock::time_point();
        const label statsThreads = max(cloud_.ompNumThreads(), label(1));
        const label collisionChunk = max(cloud_.openmpCollisionChunk(), label(1));
        const bool useFlatOccupancy = cloud_.hasOccupancyOrderedParcels();
        const DynamicList<DynamicList<dsmcParcel*>>* cellOccupancyPtr =
            useFlatOccupancy ? nullptr : &cloud_.cellOccupancy();
        const polyMesh& mesh = cloud_.mesh();
        const scalarField& cellVolumes = mesh.cellVolumes();
        const List<dsmcParcel::constantProperties>& constProps =
            cloud_.constProps();
        BinaryCollisionModel& binaryCollision = cloud_.binaryCollision();
        dsmcReactions& reactions = cloud_.reactions();
        const bool hasReactions = reactions.nReactions() > 0;
        List<autoPtr<dsmcReaction>>& reactionModels = reactions.reactions();
        const List<List<label>>& pairModelAddressing =
            reactions.pairModelAddressing();
        const volScalarField& nParticles = cloud_.nParticles();
        scalarField& collisionSelectionRemainder =
            cloud_.collisionSelectionRemainder();
        volScalarField& sigmaTcRMaxField = cloud_.sigmaTcRMax();
        const label nCells = mesh.nCells();
        const bool replicatedCollisionUseActiveCells =
            cloud_.replicatedMeshActive()
         && mesh.time().controlDict().lookupOrDefault<bool>
            (
                "replicatedMeshCollisionUseActiveCells",
                false
            );
        const labelList& collisionCells =
            replicatedCollisionUseActiveCells
          ? cloud_.occupancyActiveCells()
          : cloud_.occupancyOwnedCollisionCells();
        const label nCollisionCells = collisionCells.size();

        labelList& nCandidatesPerCell = cloud_.nCandidatesPerCell();
        cloud_.accumulateCellCostCandidates(nCandidatesPerCell, candidateCellsToClear_);
        clearCandidateCounts(nCandidatesPerCell, nCells);

        if
        (
            threadWhichSubCell_.size() != statsThreads
         || threadSubCells_.size() != statsThreads
         || threadParcelPtrs_.size() != statsThreads
         || threadTypeIds_.size() != statsThreads
         || threadCharges_.size() != statsThreads
         || threadCandidateCells_.size() != statsThreads
        )
        {
            threadWhichSubCell_.setSize(statsThreads);
            threadSubCells_.setSize(statsThreads);
            threadParcelPtrs_.setSize(statsThreads);
            threadTypeIds_.setSize(statsThreads);
            threadCharges_.setSize(statsThreads);
            threadCandidateCells_.setSize(statsThreads);
        }

        for (label threadI = 0; threadI < statsThreads; ++threadI)
        {
            if (threadSubCells_[threadI].size() != 8)
            {
                threadSubCells_[threadI].setSize(8);
            }
            threadCandidateCells_[threadI].clear();
        }

        ensureCollisionThreadFastRngs(statsThreads);

        for (label threadI = 0; threadI < statsThreads; ++threadI)
        {
            collisionThreadFastRng_[threadI]->reset
            (
                dsmcFastRng::collisionSeed
                (
                    Pstream::myProcNo(),
                    threadI,
                    mesh.time().timeIndex()
                )
            );
        }

        labelList threadCandidateCounts(statsThreads, 0);
        labelList threadAcceptedCounts(statsThreads, 0);

        auto processCell =
        [&]
        (
            const label cellI,
            const label threadI
        )
        -> label
        {
            dsmcFastRng& fastRng = collisionThreadFastRng_[threadI]->rng();

            auto randomIndex = [&](const label n) -> label
            {
                return fastRng_
                  ? fastRng.position(n)
                  : cloud_.collisionRandomLabel(0, n - 1);
            };

            auto random01 = [&]() -> scalar
            {
                return fastRng_
                  ? fastRng.sample01()
                  : cloud_.collisionSample01();
            };

            const DynamicList<dsmcParcel*>* cellParcelsPtr =
                useFlatOccupancy ? nullptr : &(*cellOccupancyPtr)[cellI];
            const label occStart =
                useFlatOccupancy ? cloud_.occupancyStart(cellI) : 0;
            const label nC =
                useFlatOccupancy
              ? cloud_.occupancyCount(cellI)
              : cellParcelsPtr->size();

            if (nC <= 1)
            {
                return 0;
            }

            const scalar selectedPairs =
                collisionSelectionRemainder[cellI]
              + 0.5*nC*(nC - 1)
               *nParticles[cellI]
               *sigmaTcRMaxField[cellI]
               *cloud_.deltaTValue(cellI)
               /cellVolumes[cellI];

            const label nCandidates(selectedPairs);
            collisionSelectionRemainder[cellI] =
                selectedPairs - nCandidates;
            nCandidatesPerCell[cellI] = nCandidates;

            if (nCandidates <= 0)
            {
                return 0;
            }

            DynamicList<label>& whichSubCell = threadWhichSubCell_[threadI];
            List<DynamicList<label>>& subCells = threadSubCells_[threadI];
            DynamicList<dsmcParcel*>& parcelPtrs = threadParcelPtrs_[threadI];
            DynamicList<label>& typeIds = threadTypeIds_[threadI];
            DynamicList<label>& charges = threadCharges_[threadI];
            DynamicList<label>& candidateCells = threadCandidateCells_[threadI];

            whichSubCell.setSize(nC);
            parcelPtrs.setSize(nC);
            typeIds.setSize(nC);
            charges.setSize(nC);
            candidateCells.append(cellI);

            label subCellCounts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            label subCellOffsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            const point& cC = mesh.cellCentres()[cellI];

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
                typeIds[i] = typeId;
                charges[i] = constProps[typeId].charge();
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
            label collisions = 0;

            for (label c = 0; c < nCandidates; ++c)
            {
                const label candidateP = randomIndex(nC);
                label candidateQ = -1;

                const List<label>& subCellPs =
                    subCells[whichSubCell[candidateP]];
                const label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    do
                    {
                        candidateQ = subCellPs[randomIndex(nSC)];
                    } while (candidateP == candidateQ);
                }
                else
                {
                    do
                    {
                        candidateQ = randomIndex(nC);
                    } while (candidateP == candidateQ);
                }

                const label typeIdP = typeIds[candidateP];
                const label typeIdQ = typeIds[candidateQ];
                const label chargeP = charges[candidateP];
                const label chargeQ = charges[candidateQ];

                if (chargeP == -1 && chargeQ == -1)
                {
                    continue;
                }

                dsmcParcel& parcelP = *parcelPtrs[candidateP];
                dsmcParcel& parcelQ = *parcelPtrs[candidateQ];

                const scalar sigmaTcR =
                    binaryCollision.sigmaTcR(parcelP, parcelQ);

                if (sigmaTcR > sigmaTcRMaxField[cellI])
                {
                    sigmaTcRMaxField[cellI] = sigmaTcR;
                }

                if ((sigmaTcR/sigmaTcRMax) > random01())
                {
                    const label rMId =
                        hasReactions
                      ? pairModelAddressing[typeIdP][typeIdQ]
                      : -1;

                    if (rMId != -1)
                    {
                        reactionModels[rMId]->reaction
                        (
                            parcelP,
                            parcelQ
                        );

                        if (reactionModels[rMId]->relax())
                        {
                            binaryCollision.collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }
                    }
                    else
                    {
                        binaryCollision.collide
                        (
                            parcelP,
                            parcelQ,
                            cellI
                        );
                    }

                    ++collisions;

                    typeIds[candidateP] = parcelP.typeId();
                    typeIds[candidateQ] = parcelQ.typeId();
                    charges[candidateP] =
                        constProps[typeIds[candidateP]].charge();
                    charges[candidateQ] =
                        constProps[typeIds[candidateQ]].charge();
                }
            }

            return collisions;
        };

        #pragma omp parallel num_threads(statsThreads)
        {
            const label threadI = omp_get_thread_num();
            label localCandidates = 0;
            label localCollisions = 0;
            dsmcFastRng& threadRng = collisionThreadFastRng_[threadI]->rng();

            if (fastRng_)
            {
                cloud_.setCollisionRngContext(&threadRng);
            }
            else
            {
                cloud_.clearCollisionRngContext();
            }

            if (cloud_.openmpCollisionSchedule() == "static")
            {
                #pragma omp for schedule(static, collisionChunk)
                for (label idx = 0; idx < nCollisionCells; ++idx)
                {
                    const label cellI = collisionCells[idx];
                    localCollisions += processCell(cellI, threadI);
                    localCandidates += nCandidatesPerCell[cellI];
                }
            }
            else if (cloud_.openmpCollisionSchedule() == "guided")
            {
                #pragma omp for schedule(guided, collisionChunk)
                for (label idx = 0; idx < nCollisionCells; ++idx)
                {
                    const label cellI = collisionCells[idx];
                    localCollisions += processCell(cellI, threadI);
                    localCandidates += nCandidatesPerCell[cellI];
                }
            }
            else
            {
                #pragma omp for schedule(dynamic, collisionChunk)
                for (label idx = 0; idx < nCollisionCells; ++idx)
                {
                    const label cellI = collisionCells[idx];
                    localCollisions += processCell(cellI, threadI);
                    localCandidates += nCandidatesPerCell[cellI];
                }
            }

            threadCandidateCounts[threadI] = localCandidates;
            threadAcceptedCounts[threadI] = localCollisions;
            cloud_.clearCollisionRngContext();
        }

        for (label threadI = 0; threadI < statsThreads; ++threadI)
        {
            const DynamicList<label>& candidateCells =
                threadCandidateCells_[threadI];

            forAll(candidateCells, i)
            {
                candidateCellsToClear_.append(candidateCells[i]);
            }
        }

        label collisionCandidates = sum(threadCandidateCounts);
        label collisions = sum(threadAcceptedCounts);
        const label localCollisionCandidates = collisionCandidates;
        const label localCollisions = collisions;

        profileLocalCollisionCandidates_ += localCollisionCandidates;
        profileLocalCollisions_ += localCollisions;
        cloud_.addToCumulativeCollisions(localCollisions, localCollisionCandidates);

        if (subphaseTimers)
        {
            profileCollisionLocalLoopWall_ +=
                elapsedCollisionSeconds(localLoopStart);
        }

        label globalCollisions = localCollisions;
        label globalCollisionCandidates = localCollisionCandidates;

        if (!collisionReduceOnlyOnOutput_)
        {
            const CollisionClock::time_point reduceStart =
                subphaseTimers
              ? CollisionClock::now()
              : CollisionClock::time_point();
            reduce(globalCollisions, sumOp<label>());
            reduce(globalCollisionCandidates, sumOp<label>());
            if (subphaseTimers)
            {
                profileCollisionReduceWall_ +=
                    elapsedCollisionSeconds(reduceStart);
            }
        }

        const CollisionClock::time_point sigmaStart =
            subphaseTimers ? CollisionClock::now() : CollisionClock::time_point();
        cloud_.sigmaTcRMax().correctBoundaryConditions();
        if (subphaseTimers)
        {
            profileCollisionSigmaWall_ += elapsedCollisionSeconds(sigmaStart);
        }

        infoCounter_++;

        if(infoCounter_ >= cloud_.nTerminalOutputs())
        {
            if (collisionReduceOnlyOnOutput_ && collisionOutputGlobalReduce_)
            {
                const CollisionClock::time_point reduceStart =
                    subphaseTimers
                  ? CollisionClock::now()
                  : CollisionClock::time_point();
                reduce(globalCollisions, sumOp<label>());
                reduce(globalCollisionCandidates, sumOp<label>());
                if (subphaseTimers)
                {
                    profileCollisionReduceWall_ +=
                        elapsedCollisionSeconds(reduceStart);
                }
            }

            if (subphaseTimers)
            {
                profileCollisionPreOutputWall_ +=
                    elapsedCollisionSeconds(collideWallStart);
            }

            const bool replicatedRawMpi =
                cloud_.replicatedMeshActive()
             && cloud_.replicatedMesh().nProcs() > 1;

            if (!replicatedRawMpi || collisionOutputGlobalReduce_)
            {
                printReplicatedRankCollisionDetail
                (
                    localCollisions,
                    localCollisionCandidates,
                    globalCollisions,
                    globalCollisionCandidates
                );
            }

            if (replicatedRawMpi)
            {
                if
                (
                    cloud_.replicatedMesh().myRank() == 0
                 && globalCollisionCandidates
                )
                {
                    if (!collisionOutputGlobalReduce_)
                    {
                        Info<< "    Collision global reduction skipped"
                            << " (rank 0 local diagnostics)" << nl;
                    }
                    Info<< "    Collisions                      = "
                        << globalCollisions << nl
                        << "    Collision candidates           = "
                        << globalCollisionCandidates << nl
                        << "    Collision acceptance rate      = "
                        << scalar(globalCollisions)
                           /scalar(globalCollisionCandidates) << nl
                        << endl;
                }
                else if (cloud_.replicatedMesh().myRank() == 0)
                {
                    Info<< "    No collisions" << endl;
                }
            }
            else if (cloud_.isOutputRank() && globalCollisionCandidates)
            {
                Info<< "    Collisions                      = "
                    << globalCollisions << nl
                    << "    Collision candidates           = "
                    << globalCollisionCandidates << nl
                    << "    Collision acceptance rate      = "
                    << scalar(globalCollisions)
                       /scalar(globalCollisionCandidates) << nl
                    << endl;
            }
            else if (cloud_.isOutputRank())
            {
                Info<< "    No collisions" << endl;
            }

            infoCounter_ = 0;
        }
        else if (subphaseTimers)
        {
            profileCollisionPreOutputWall_ +=
                elapsedCollisionSeconds(collideWallStart);
        }

        return;
    }
    #endif

    const CollisionClock::time_point localLoopStart =
        subphaseTimers ? CollisionClock::now() : CollisionClock::time_point();

    label collisionCandidates = 0;

    label collisions = 0;

    const polyMesh& mesh = cloud_.mesh();
    const label nCells = mesh.nCells();
    const bool useFlatOccupancy = cloud_.hasOccupancyOrderedParcels();
    const DynamicList<DynamicList<dsmcParcel*>>* cellOccupancyPtr =
        useFlatOccupancy ? nullptr : &cloud_.cellOccupancy();
    const scalarField& cellVolumes = mesh.cellVolumes();
    const List<dsmcParcel::constantProperties>& constProps =
        cloud_.constProps();
    BinaryCollisionModel& binaryCollision = cloud_.binaryCollision();
    dsmcReactions& reactions = cloud_.reactions();
    const bool hasReactions = reactions.nReactions() > 0;
    List<autoPtr<dsmcReaction>>& reactionModels = reactions.reactions();
    const List<List<label>>& pairModelAddressing =
        reactions.pairModelAddressing();
    const volScalarField& nParticles = cloud_.nParticles();
    scalarField& collisionSelectionRemainder =
        cloud_.collisionSelectionRemainder();
    volScalarField& sigmaTcRMaxField = cloud_.sigmaTcRMax();

    const bool replicatedCollisionUseActiveCells =
        cloud_.replicatedMeshActive()
     && mesh.time().controlDict().lookupOrDefault<bool>
        (
            "replicatedMeshCollisionUseActiveCells",
            false
        );
    const labelList& collisionCells =
        replicatedCollisionUseActiveCells
      ? cloud_.occupancyActiveCells()
      : cloud_.occupancyOwnedCollisionCells();

    labelList& nCandidatesPerCell = cloud_.nCandidatesPerCell();
    cloud_.accumulateCellCostCandidates(nCandidatesPerCell, candidateCellsToClear_);
    clearCandidateCounts(nCandidatesPerCell, nCells);

    if (threadWhichSubCell_.size() < 1)
    {
        threadWhichSubCell_.setSize(1);
        threadSubCells_.setSize(1);
        threadParcelPtrs_.setSize(1);
        threadTypeIds_.setSize(1);
        threadCharges_.setSize(1);
    }

    if (threadSubCells_[0].size() != 8)
    {
        threadSubCells_[0].setSize(8);
    }

    ensureCollisionThreadFastRngs(1);

    collisionThreadFastRng_[0]->reset
    (
        dsmcFastRng::collisionSeed
        (
            Pstream::myProcNo(),
            0,
            mesh.time().timeIndex()
        )
    );
    dsmcFastRng& fastRng = collisionThreadFastRng_[0]->rng();

    if (fastRng_)
    {
        cloud_.setCollisionRngContext(&fastRng);
    }
    else
    {
        cloud_.clearCollisionRngContext();
    }

    auto randomIndex = [&](const label n) -> label
    {
        if (!fastRng_)
        {
            return cloud_.collisionRandomLabel(0, n - 1);
        }

        return fastRng.position(n);
    };

    auto random01 = [&]() -> scalar
    {
        return fastRng_
             ? fastRng.sample01()
             : cloud_.collisionSample01();
    };

    auto processCell =
    [&]
    (
        const label cellI
    )
    {
        const scalar deltaT = cloud_.deltaTValue(cellI);

        const DynamicList<dsmcParcel*>* cellParcelsPtr =
            useFlatOccupancy ? nullptr : &(*cellOccupancyPtr)[cellI];
        const label occStart =
            useFlatOccupancy ? cloud_.occupancyStart(cellI) : 0;

        const scalar& cellVolume = cellVolumes[cellI];

        const label nC =
            useFlatOccupancy
          ? cloud_.occupancyCount(cellI)
          : cellParcelsPtr->size();

        if (nC <= 1)
        {
            return;
        }

        const scalar sigmaTcRMax = sigmaTcRMaxField[cellI];

        const scalar selectedPairs =
            collisionSelectionRemainder[cellI]
            + 0.5*nC*(nC - 1)*nParticles[cellI]*sigmaTcRMax*deltaT
            /cellVolume;

        const label nCandidates(selectedPairs);

        collisionSelectionRemainder[cellI] =
            selectedPairs - nCandidates;

        nCandidatesPerCell[cellI] = nCandidates;

        if (nCandidates <= 0)
        {
            return;
        }

        candidateCellsToClear_.append(cellI);
        collisionCandidates += nCandidates;

        DynamicList<label>& whichSubCell = threadWhichSubCell_[0];
        List<DynamicList<label>>& subCells = threadSubCells_[0];
        DynamicList<dsmcParcel*>& parcelPtrs = threadParcelPtrs_[0];
        DynamicList<label>& typeIds = threadTypeIds_[0];
        DynamicList<label>& charges = threadCharges_[0];

        whichSubCell.setSize(nC);
        parcelPtrs.setSize(nC);
        typeIds.setSize(nC);
        charges.setSize(nC);

        label subCellCounts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        label subCellOffsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        const point& cC = mesh.cellCentres()[cellI];

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
            typeIds[i] = typeId;
            charges[i] = constProps[typeId].charge();
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

        for (label c = 0; c < nCandidates; c++)
        {
            const label candidateP = randomIndex(nC);
            label candidateQ = -1;

            const List<label>& subCellPs = subCells[whichSubCell[candidateP]];

            const label nSC = subCellPs.size();

            if (nSC > 1)
            {
                do
                {
                    candidateQ = subCellPs[randomIndex(nSC)];

                } while (candidateP == candidateQ);
            }
            else
            {
                do
                {
                    candidateQ = randomIndex(nC);

                } while (candidateP == candidateQ);
            }

            const label typeIdP = typeIds[candidateP];
            const label typeIdQ = typeIds[candidateQ];
            const label chargeP = charges[candidateP];
            const label chargeQ = charges[candidateQ];

            // Do not allow electron-electron collisions.
            if (chargeP == -1 && chargeQ == -1)
            {
                continue;
            }

            dsmcParcel& parcelP = *parcelPtrs[candidateP];
            dsmcParcel& parcelQ = *parcelPtrs[candidateQ];

            const scalar sigmaTcR = binaryCollision.sigmaTcR
            (
                parcelP,
                parcelQ
            );

            if (sigmaTcR > sigmaTcRMaxField[cellI])
            {
                sigmaTcRMaxField[cellI] = sigmaTcR;
            }

            if ((sigmaTcR/sigmaTcRMax) > random01())
            {
                const label rMId =
                    hasReactions
                  ? pairModelAddressing[typeIdP][typeIdQ]
                  : -1;

                if(rMId != -1)
                {
                    reactionModels[rMId]->reaction
                    (
                        parcelP,
                        parcelQ
                    );

                    if(reactionModels[rMId]->relax())
                    {
                        binaryCollision.collide
                        (
                            parcelP,
                            parcelQ,
                            cellI
                        );
                    }
                }
                else
                {
                    binaryCollision.collide
                    (
                        parcelP,
                        parcelQ,
                        cellI
                    );
                }

                collisions++;

                typeIds[candidateP] = parcelP.typeId();
                typeIds[candidateQ] = parcelQ.typeId();
                charges[candidateP] =
                    constProps[typeIds[candidateP]].charge();
                charges[candidateQ] =
                    constProps[typeIds[candidateQ]].charge();
            }
        }
    };

    forAll(collisionCells, i)
    {
        processCell(collisionCells[i]);
    }

    cloud_.clearCollisionRngContext();

    const label localCollisionCandidates = collisionCandidates;
    const label localCollisions = collisions;

    profileLocalCollisionCandidates_ += localCollisionCandidates;
    profileLocalCollisions_ += localCollisions;
    cloud_.addToCumulativeCollisions(localCollisions, localCollisionCandidates);

    if (subphaseTimers)
    {
        profileCollisionLocalLoopWall_ +=
            elapsedCollisionSeconds(localLoopStart);
    }

    label globalCollisions = localCollisions;
    label globalCollisionCandidates = localCollisionCandidates;

    if (!collisionReduceOnlyOnOutput_)
    {
        const CollisionClock::time_point reduceStart =
            subphaseTimers ? CollisionClock::now() : CollisionClock::time_point();
        reduce(globalCollisions, sumOp<label>());
        reduce(globalCollisionCandidates, sumOp<label>());
        if (subphaseTimers)
        {
            profileCollisionReduceWall_ += elapsedCollisionSeconds(reduceStart);
        }
    }

    const CollisionClock::time_point sigmaStart =
        subphaseTimers ? CollisionClock::now() : CollisionClock::time_point();
    cloud_.sigmaTcRMax().correctBoundaryConditions();
    if (subphaseTimers)
    {
        profileCollisionSigmaWall_ += elapsedCollisionSeconds(sigmaStart);
    }

    infoCounter_++;

    if(infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisionReduceOnlyOnOutput_ && collisionOutputGlobalReduce_)
        {
            const CollisionClock::time_point reduceStart =
                subphaseTimers ? CollisionClock::now() : CollisionClock::time_point();
            reduce(globalCollisions, sumOp<label>());
            reduce(globalCollisionCandidates, sumOp<label>());
            if (subphaseTimers)
            {
                profileCollisionReduceWall_ +=
                    elapsedCollisionSeconds(reduceStart);
            }
        }

        if (subphaseTimers)
        {
            profileCollisionPreOutputWall_ +=
            elapsedCollisionSeconds(collideWallStart);
    }

        const bool replicatedRawMpi =
            cloud_.replicatedMeshActive()
         && cloud_.replicatedMesh().nProcs() > 1;

        if (!replicatedRawMpi || collisionOutputGlobalReduce_)
        {
            printReplicatedRankCollisionDetail
            (
                localCollisions,
                localCollisionCandidates,
                globalCollisions,
                globalCollisionCandidates
            );
        }

        if (replicatedRawMpi)
        {
            if
            (
                cloud_.replicatedMesh().myRank() == 0
             && globalCollisionCandidates
            )
            {
                if (!collisionOutputGlobalReduce_)
                {
                    Info<< "    Collision global reduction skipped"
                        << " (rank 0 local diagnostics)" << nl;
                }
                Info<< "    Collisions                      = "
                    << globalCollisions << nl
                    << "    Collision candidates           = "
                    << globalCollisionCandidates << nl
                    << "    Collision acceptance rate      = "
                    << scalar(globalCollisions)
                       /scalar(globalCollisionCandidates) << nl
                    << endl;
            }
            else if (cloud_.replicatedMesh().myRank() == 0)
            {
                Info<< "    No collisions" << endl;
            }

            infoCounter_ = 0;
        }
        else if (cloud_.isOutputRank() && globalCollisionCandidates)
        {
            Info<< "    Collisions                      = "
                << globalCollisions << nl
                << "    Collision candidates           = "
                << globalCollisionCandidates << nl
                << "    Collision acceptance rate      = "
                << scalar(globalCollisions)
                   /scalar(globalCollisionCandidates) << nl
                << endl;

            infoCounter_ = 0;
        }
        else if (cloud_.isOutputRank())
        {
            Info<< "    No collisions" << endl;

            infoCounter_ = 0;
        }
        else
        {
            infoCounter_ = 0;
        }
    }
    else if (subphaseTimers)
    {
        profileCollisionPreOutputWall_ +=
            elapsedCollisionSeconds(collideWallStart);
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
