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
#include "PstreamBuffers.H"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <vector>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounter, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


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
    threadCharges_()
//     propsDict_(dict.subDict(typeName + "Properties"))
{
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
    const dictionary& controlDict = mesh.time().controlDict();
    const bool collisionFastRng =
        controlDict.lookupOrDefault<bool>("collisionFastRng", false);
    const bool dlbOffloadPlanner =
        controlDict.lookupOrDefault<bool>("dlbOffloadPlanner", false);
    const bool dlbOffloadExecute =
        controlDict.lookupOrDefault<bool>("dlbOffloadExecute", false)
     || controlDict.lookupOrDefault<bool>("dlbOffloadExperimentalExecute", false);
    const bool dlbOffloadReport =
        controlDict.lookupOrDefault<bool>("dlbOffloadReport", false);
    const scalar dlbOffloadBalanceTarget =
        controlDict.lookupOrDefault<scalar>("dlbOffloadBalanceTarget", 1.30);
    const scalar dlbOffloadMaxFraction =
        controlDict.lookupOrDefault<scalar>("dlbOffloadMaxFraction", 0.30);
    const scalar dlbOffloadMaxParcelFraction =
        max(controlDict.lookupOrDefault<scalar>("dlbOffloadMaxParcelFraction", 1.0), scalar(0));
    const label dlbOffloadTaskCells =
        max(controlDict.lookupOrDefault<label>("dlbOffloadTaskCells", 256), label(1));
    const label dlbOffloadMinCandidates =
        max(controlDict.lookupOrDefault<label>("dlbOffloadMinCandidates", 1), label(1));
    const scalar dlbOffloadMinCandidatesPerParcel =
        max
        (
            controlDict.lookupOrDefault<scalar>("dlbOffloadMinCandidatesPerParcel", 0.0),
            scalar(0)
        );
    const bool dlbOffloadEfficiencySort =
        controlDict.lookupOrDefault<bool>("dlbOffloadEfficiencySort", false);
    const scalar dlbOffloadImbalance =
        max(controlDict.lookupOrDefault<scalar>("dlbOffloadImbalance", 1.0), scalar(1));
    const scalar dlbOffloadRemoteCostFactor =
        max(controlDict.lookupOrDefault<scalar>("dlbOffloadRemoteCostFactor", 1.0), scalar(1));
    const scalar dlbOffloadBudgetSafety =
        max(controlDict.lookupOrDefault<scalar>("dlbOffloadBudgetSafety", 1.0), scalar(1));
    const bool dlbOffloadUseTimerCost =
        controlDict.lookupOrDefault<bool>("dlbOffloadUseTimerCost", true);
    const bool dlbOffloadRequireTimerCost =
        controlDict.lookupOrDefault<bool>("dlbOffloadRequireTimerCost", false);
    const label dlbOffloadMinSelectedCandidates =
        max
        (
            controlDict.lookupOrDefault<label>("dlbOffloadMinSelectedCandidates", 1),
            label(1)
        );
    const label dlbOffloadIntervalInput =
        controlDict.lookupOrDefault<label>("dlbOffloadInterval", cloud_.nTerminalOutputs());
    const label dlbOffloadInterval =
        dlbOffloadIntervalInput > 0
      ? dlbOffloadIntervalInput
      : max(cloud_.nTerminalOutputs(), label(1));
    const bool dlbOffloadCompactSingleSpeciesPayload =
        controlDict.lookupOrDefault<bool>("dlbOffloadCompactSingleSpeciesPayload", false);
    const bool dlbOffloadPersistent =
        controlDict.lookupOrDefault<bool>("dlbOffloadPersistent", false);
    const label dlbOffloadPersistentLeaseSteps =
        max
        (
            controlDict.lookupOrDefault<label>("dlbOffloadPersistentLeaseSteps", 5),
            label(1)
        );
    const bool dlbOffloadIntervalHit =
        dlbOffloadInterval <= 1
     || (mesh.time().timeIndex() % dlbOffloadInterval) == 0;
    static labelList dlbPersistentDonorCells;
    static label dlbPersistentDonorPeer = -1;
    static label dlbPersistentDonorExpire = -1;
    static label dlbPersistentHelperPeer = -1;
    static label dlbPersistentHelperExpire = -1;
    static scalar dlbPrevRankCollisionWall = -1.0;
    static scalar dlbPrevSecondsPerCandidate = -1.0;

    const bool dlbBaseActive =
        dlbOffloadPlanner
     && dlbOffloadExecute
     && Pstream::parRun()
     && Pstream::nProcs() >= 2
     && !cloud_.reactionsActive();
    const bool dlbCompactStaticFields =
        dlbOffloadCompactSingleSpeciesPayload
     && cloud_.typeIdList().size() == 1;
    const label dlbCompactTypeId = 0;
    // Determine peer for offload communication.
    // For 2-rank: trivial peer = 1 - myProcNo.
    // For multi-rank: heaviest rank offloads to ALL below-average helpers.
    label peerProc = -1;
    labelList allProcCands;     // saved AllGather data for planner reuse
    labelList donorHelperList;  // multi-rank donor: all helpers
    label donorProc = -1;       // multi-rank helper: the donor
    if (dlbBaseActive)
    {
        if (Pstream::nProcs() == 2)
        {
            peerProc = 1 - Pstream::myProcNo();
        }
        else
        {
            label localCandQuick = 0;
            for (label cellI = 0; cellI < nCells; ++cellI)
            {
                localCandQuick += cloud_.nCandidatesPerCell()[cellI];
            }

            allProcCands.setSize(Pstream::nProcs(), 0);
            allProcCands[Pstream::myProcNo()] = localCandQuick;
            Pstream::allGatherList(allProcCands);

            scalar totalCand = 0;
            forAll(allProcCands, pi) totalCand += scalar(allProcCands[pi]);
            const scalar avgCand = totalCand / scalar(Pstream::nProcs());

            // Find heaviest rank as donor
            label heaviest = 0;
            for (label pi = 1; pi < Pstream::nProcs(); ++pi)
            {
                if (allProcCands[pi] > allProcCands[heaviest]) heaviest = pi;
            }

            if (scalar(allProcCands[heaviest]) > avgCand)
            {
                // Collect all below-average ranks as helpers
                DynamicList<label> helpers;
                for (label pi = 0; pi < Pstream::nProcs(); ++pi)
                {
                    if (pi != heaviest && scalar(allProcCands[pi]) < avgCand)
                    {
                        helpers.append(pi);
                    }
                }

                // Sort helpers by spare capacity (descending)
                std::sort
                (
                    helpers.begin(),
                    helpers.end(),
                    [&](const label a, const label b)
                    {
                        return (avgCand - scalar(allProcCands[a]))
                             > (avgCand - scalar(allProcCands[b]));
                    }
                );

                donorHelperList.transfer(helpers);

                // Multi-helper: all below-average ranks help the heaviest
                if (Pstream::myProcNo() == heaviest)
                {
                    peerProc = donorHelperList.size() > 0
                             ? donorHelperList[0] : -1;
                }
                else if (donorHelperList.found(Pstream::myProcNo()))
                {
                    donorProc = heaviest;
                    peerProc = heaviest;
                }
            }
        }
    }
    const label dlbTimeIndex = mesh.time().timeIndex();
    const bool dlbPersistentDonorActive =
        dlbBaseActive
     && dlbOffloadPersistent
     && dlbPersistentDonorPeer == peerProc
     && dlbPersistentDonorExpire >= dlbTimeIndex
     && dlbPersistentDonorCells.size() > 0;
    const bool dlbPersistentHelperActive =
        dlbBaseActive
     && dlbOffloadPersistent
     && dlbPersistentHelperPeer == peerProc
     && dlbPersistentHelperExpire >= dlbTimeIndex;
    const bool dlbActive =
        dlbBaseActive
     && (dlbOffloadIntervalHit || dlbPersistentDonorActive || dlbPersistentHelperActive);
    boolList offloadCells(nCells, false);
    DynamicList<label> offloadCellLabels;
    label offloadCandidateCount = 0;
    label offloadParcelCount = 0;
    label remoteAcceptedCount = 0;
    label remoteActiveCellCount = 0;
    scalar dlbOffloadWallTime = 0.0;
    label dlbLocalCandidateTotal = 0;
    label dlbPeerCandidateTotal = 0;
    scalar dlbCandidateImbalance = 0.0;
    label dlbDesiredOffloadCandidates = 0;
    label dlbSelectedCells = 0;
    label dlbSelectedCandidatesBeforeGate = 0;
    scalar dlbCurrentMaxCost = 0.0;
    scalar dlbPredictedDonorCost = 0.0;
    scalar dlbPredictedHelperCost = 0.0;
    scalar dlbPredictedMaxCost = 0.0;
    bool dlbDonorCandidate = false;
    bool dlbGateCancelled = false;
    label dlbSelectedParcels = 0;
    scalar dlbPreviousLocalCollisionWall = -1.0;
    scalar dlbPreviousPeerCollisionWall = -1.0;
    scalar dlbPreviousLocalSecondsPerCandidate = -1.0;
    scalar dlbPreviousPeerSecondsPerCandidate = -1.0;
    bool dlbCompactTaskPayload = false;
    bool dlbIncomingCompactPayload = false;
    bool dlbPersistentTaskPayload = false;
    label dlbPersistentTaskExpire = -1;
    bool dlbUsingPersistentLease = false;
    bool dlbResultsPending = false;
    PstreamBuffers resultBufs;

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

    auto collectCellParcels =
    [&]
    (
        const label cellI,
        DynamicList<dsmcParcel*>& parcelPtrs
    )
    {
        parcelPtrs.clear();

        if (useFlatOccupancy)
        {
            const label occStart = cloud_.occupancyStart(cellI);
            const label nC = cloud_.occupancyCount(cellI);
            parcelPtrs.setCapacity(nC);

            for (label i = 0; i < nC; ++i)
            {
                parcelPtrs.append(cloud_.occupancyParcel(occStart + i));
            }
        }
        else
        {
            const DynamicList<dsmcParcel*>& cellParcels =
                (*cellOccupancyPtr)[cellI];
            parcelPtrs.setCapacity(cellParcels.size());

            forAll(cellParcels, i)
            {
                parcelPtrs.append(cellParcels[i]);
            }
        }
    };

    auto writeParcelState =
    [&]
    (
        Ostream& os,
        const dsmcParcel& p,
        const bool compactStaticFields
    )
    {
        os  << p.U()
            << p.RWF()
            << p.ERot()
            << p.ELevel();

        if (compactStaticFields)
        {
            os << p.classification();
        }
        else
        {
            os  << p.typeId()
                << p.newParcel()
                << p.classification();
        }

        os  << p.vibLevel()
            << token::SPACE;
    };

    auto readParcelState =
    [&]
    (
        Istream& is,
        vector& U,
        scalar& RWF,
        scalar& ERot,
        label& ELevel,
        label& typeId,
        label& newParcel,
        label& classification,
        labelList& vibLevel,
        const bool compactStaticFields
    )
    {
        is  >> U
            >> RWF
            >> ERot
            >> ELevel;

        if (compactStaticFields)
        {
            typeId = dlbCompactTypeId;
            newParcel = -1;
            is >> classification;
        }
        else
        {
            is  >> typeId
                >> newParcel
                >> classification;
        }

        is >> vibLevel;
    };

    auto applyParcelState =
    [&]
    (
        Istream& is,
        dsmcParcel& p,
        const bool compactStaticFields
    )
    {
        vector U;
        scalar RWF = 1.0;
        scalar ERot = 0.0;
        label ELevel = 0;
        label typeId = -1;
        label newParcel = -1;
        label classification = 0;
        labelList vibLevel;

        readParcelState
        (
            is,
            U,
            RWF,
            ERot,
            ELevel,
            typeId,
            newParcel,
            classification,
            vibLevel,
            compactStaticFields
        );

        p.U() = U;
        p.RWF() = RWF;
        p.ERot() = ERot;
        p.ELevel() = ELevel;

        if (!compactStaticFields)
        {
            p.typeId() = typeId;
            p.newParcel() = newParcel;
        }

        p.classification() = classification;
        p.vibLevel() = vibLevel;
    };

    auto executeRemoteCell =
    [&]
    (
        const label donorCellI,
        DynamicList<dsmcParcel*>& parcelPtrs,
        const labelList& subCellIds,
        const label nCandidates,
        const scalar sigmaTcRMaxInitial,
        scalar& sigmaTcRMaxUpdated,
        boolList& dirtyParcels
    )
    {
        auto rngPosRemote = [&](label n) -> label
        {
            return collisionFastRng
                 ? threadFastRng[0].position(n)
                 : cloud_.randomLabel(0, n - 1);
        };
        auto rng01Remote = [&]() -> scalar
        {
            return collisionFastRng
                 ? threadFastRng[0].sample01()
                 : cloud_.rndGen().sample01<scalar>();
        };

        const label nC = parcelPtrs.size();
        label acceptedCollisions = 0;

        sigmaTcRMaxUpdated = sigmaTcRMaxInitial;
        dirtyParcels.setSize(nC);
        dirtyParcels = false;

        if (nC <= 1 || nCandidates <= 0)
        {
            return acceptedCollisions;
        }

        DynamicList<label> whichSubCell;
        List<DynamicList<label>> subCells(8);
        DynamicList<vector> velocities;
        DynamicList<label> typeIds;
        DynamicList<label> charges;
        label subCellCounts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        label subCellOffsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};

        whichSubCell.setSize(nC);
        velocities.setSize(nC);
        typeIds.setSize(nC);
        charges.setSize(nC);

        for (label i = 0; i < nC; ++i)
        {
            dsmcParcel* pPtr = parcelPtrs[i];
            const label typeId = pPtr->typeId();
            const label subCell = subCellIds[i];

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

        for (label c = 0; c < nCandidates; ++c)
        {
            const label candidateP = rngPosRemote(nC);
            label candidateQ = -1;

            const List<label>& subCellPs = subCells[whichSubCell[candidateP]];
            const label nSC = subCellPs.size();

            if (nSC > 1)
            {
                do
                {
                    candidateQ = subCellPs[rngPosRemote(nSC)];
                } while (candidateP == candidateQ);
            }
            else
            {
                do
                {
                    candidateQ = rngPosRemote(nC);
                } while (candidateP == candidateQ);
            }

            const label typeIdP = typeIds[candidateP];
            const label typeIdQ = typeIds[candidateQ];
            label chargeP = charges[candidateP];
            label chargeQ = charges[candidateQ];

            if (chargeP == -2)
            {
                chargeP = cloud_.constProps(typeIdP).charge();
            }

            if (chargeQ == -2)
            {
                chargeQ = cloud_.constProps(typeIdQ).charge();
            }

            if (chargeP == -1 && chargeQ == -1)
            {
                continue;
            }

            const scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
            (
                *parcelPtrs[candidateP],
                *parcelPtrs[candidateQ]
            );

            if (sigmaTcR > sigmaTcRMaxUpdated)
            {
                sigmaTcRMaxUpdated = sigmaTcR;
            }

            if ((sigmaTcR/sigmaTcRMaxInitial) > rng01Remote())
            {
                dsmcParcel& parcelP = *parcelPtrs[candidateP];
                dsmcParcel& parcelQ = *parcelPtrs[candidateQ];

                cloud_.binaryCollision().collide
                (
                    parcelP,
                    parcelQ,
                    donorCellI
                );

                ++acceptedCollisions;
                dirtyParcels[candidateP] = true;
                dirtyParcels[candidateQ] = true;

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

        return acceptedCollisions;
    };

    if (dlbActive)
    {
        using dlb_clock_type = std::chrono::steady_clock;
        const auto dlbBegin = dlb_clock_type::now();
        const bool dlbLeaseStep =
            dlbPersistentDonorActive || dlbPersistentHelperActive;
        const bool dlbPlannerThisStep =
            dlbOffloadIntervalHit && !dlbLeaseStep;

        label localCandidateTotal = 0;

        if (dlbPlannerThisStep)
        {
            for (label cellI = 0; cellI < nCells; ++cellI)
            {
                localCandidateTotal += cloud_.nCandidatesPerCell()[cellI];
            }

            label peerCandidateTotal = 0;

            if (Pstream::nProcs() == 2)
            {
                // Original 2-rank path: direct peer exchange
                PstreamBuffers totalsBufs;

                {
                    UOPstream os(peerProc, totalsBufs);
                    os  << localCandidateTotal
                        << dlbPrevRankCollisionWall
                        << dlbPrevSecondsPerCandidate;
                }

                totalsBufs.finishedSends();

                {
                    UIPstream is(peerProc, totalsBufs);
                    is  >> peerCandidateTotal
                        >> dlbPreviousPeerCollisionWall
                        >> dlbPreviousPeerSecondsPerCandidate;
                }
            }
            else if (peerProc >= 0)
            {
                // Multi-rank: use saved AllGather data from pairing phase.
                // No blocking exchange needed — avoids deadlock with
                // multi-helper where donor communicates with N helpers.
                peerCandidateTotal = allProcCands[peerProc];
                dlbPreviousPeerCollisionWall = -1.0;
                dlbPreviousPeerSecondsPerCandidate = -1.0;
            }

            dlbPreviousLocalCollisionWall = dlbPrevRankCollisionWall;
            dlbPreviousLocalSecondsPerCandidate = dlbPrevSecondsPerCandidate;

            const scalar averageCandidates =
                0.5*scalar(localCandidateTotal + peerCandidateTotal);
            const scalar candidateImbalance =
                averageCandidates > SMALL
              ? max(scalar(localCandidateTotal), scalar(peerCandidateTotal))/averageCandidates
              : scalar(0);
            const label targetLocalCandidates =
                label(dlbOffloadBalanceTarget*averageCandidates);
            const label maxOffloadCandidates =
                label(dlbOffloadMaxFraction*scalar(localCandidateTotal));
            const label balanceDesiredOffloadCandidates =
                max
                (
                    label(0),
                    min
                    (
                        localCandidateTotal - targetLocalCandidates,
                        maxOffloadCandidates
                    )
                );
            label desiredOffloadCandidates = balanceDesiredOffloadCandidates;
            const bool timerBudgetActive =
                dlbOffloadUseTimerCost
             && dlbPreviousLocalCollisionWall > SMALL
             && dlbPreviousPeerCollisionWall >= 0.0
             && dlbPreviousLocalSecondsPerCandidate > SMALL
             && dlbPreviousPeerSecondsPerCandidate > SMALL
             && dlbPreviousLocalCollisionWall
              > dlbPreviousPeerCollisionWall + SMALL;

            if (timerBudgetActive)
            {
                const scalar timerDesiredCandidatesScalar =
                    (
                        dlbPreviousLocalCollisionWall
                      - dlbPreviousPeerCollisionWall
                    )
                   /
                    max
                    (
                        dlbPreviousLocalSecondsPerCandidate
                      + dlbPreviousPeerSecondsPerCandidate
                       *dlbOffloadRemoteCostFactor,
                        SMALL
                    );
                const label timerDesiredOffloadCandidates =
                    max
                    (
                        label(0),
                        min(label(timerDesiredCandidatesScalar), maxOffloadCandidates)
                    );

                if (desiredOffloadCandidates > 0)
                {
                    desiredOffloadCandidates =
                        min(desiredOffloadCandidates, timerDesiredOffloadCandidates);
                }
                else
                {
                    desiredOffloadCandidates = timerDesiredOffloadCandidates;
                }
            }

            const bool donorByCandidates = localCandidateTotal > peerCandidateTotal;
            const bool donorByTimer =
                timerBudgetActive
             && dlbPreviousLocalCollisionWall
              > dlbPreviousPeerCollisionWall + SMALL;
            const bool donor =
                (!dlbOffloadRequireTimerCost || timerBudgetActive)
             &&
                (donorByCandidates || donorByTimer)
             && (candidateImbalance >= dlbOffloadImbalance || donorByTimer)
             && desiredOffloadCandidates >= dlbOffloadMinCandidates
             && desiredOffloadCandidates >= dlbOffloadMinSelectedCandidates;

            dlbLocalCandidateTotal = localCandidateTotal;
            dlbPeerCandidateTotal = peerCandidateTotal;
            dlbCandidateImbalance = candidateImbalance;
            dlbDesiredOffloadCandidates = desiredOffloadCandidates;
            dlbDonorCandidate = localCandidateTotal > peerCandidateTotal;

            if (donor)
            {
                struct offloadCellInfo
                {
                    label cellI;
                    label nCandidates;
                    label nParcels;
                    scalar efficiency;
                };

                std::vector<offloadCellInfo> candidateCellInfo;
                candidateCellInfo.reserve(nCells);
                const label localParcelBudget =
                    dlbOffloadMaxParcelFraction >= scalar(1)
                  ? max(cloud_.size(), label(1))
                  : max
                    (
                        label(1),
                        label
                        (
                            dlbOffloadMaxParcelFraction
                          * scalar(max(cloud_.size(), label(1)))
                        )
                    );

                for (label cellI = 0; cellI < nCells; ++cellI)
                {
                    const label nCandidates = cloud_.nCandidatesPerCell()[cellI];

                    if (nCandidates >= dlbOffloadMinCandidates)
                    {
                        const label nC =
                            useFlatOccupancy
                          ? cloud_.occupancyCount(cellI)
                          : (*cellOccupancyPtr)[cellI].size();

                        if (nC > 1)
                        {
                            const scalar efficiency =
                                scalar(nCandidates)/scalar(max(nC, label(1)));

                            if (efficiency + SMALL >= dlbOffloadMinCandidatesPerParcel)
                            {
                                candidateCellInfo.push_back
                                (
                                    {cellI, nCandidates, nC, efficiency}
                                );
                            }
                        }
                    }
                }

                std::sort
                (
                    candidateCellInfo.begin(),
                    candidateCellInfo.end(),
                    [&](const offloadCellInfo& a, const offloadCellInfo& b)
                    {
                        if (dlbOffloadEfficiencySort)
                        {
                            if (mag(a.efficiency - b.efficiency) > SMALL)
                            {
                                return a.efficiency > b.efficiency;
                            }
                        }

                        if (a.nCandidates != b.nCandidates)
                        {
                            return a.nCandidates > b.nCandidates;
                        }

                        return a.nParcels < b.nParcels;
                    }
                );

                for (const offloadCellInfo& info : candidateCellInfo)
                {
                    if
                    (
                        offloadCandidateCount >= desiredOffloadCandidates
                     || offloadCellLabels.size() >= dlbOffloadTaskCells
                    )
                    {
                        break;
                    }

                    if (offloadParcelCount + info.nParcels > localParcelBudget)
                    {
                        continue;
                    }

                    const label cellI = info.cellI;
                    offloadCells[cellI] = true;
                    offloadCellLabels.append(cellI);
                    offloadCandidateCount += info.nCandidates;
                    offloadParcelCount += info.nParcels;
                }

                dlbSelectedCells = offloadCellLabels.size();
                dlbSelectedCandidatesBeforeGate = offloadCandidateCount;
                dlbSelectedParcels = offloadParcelCount;
                if (timerBudgetActive)
                {
                    dlbCurrentMaxCost =
                        max
                        (
                            dlbPreviousLocalCollisionWall,
                            dlbPreviousPeerCollisionWall
                        );
                    dlbPredictedDonorCost =
                        max
                        (
                            scalar(0),
                            dlbPreviousLocalCollisionWall
                          - dlbPreviousLocalSecondsPerCandidate
                           *scalar(offloadCandidateCount)
                        );
                    dlbPredictedHelperCost =
                        dlbPreviousPeerCollisionWall
                      + dlbPreviousPeerSecondsPerCandidate
                       *dlbOffloadRemoteCostFactor
                       *scalar(offloadCandidateCount);
                    dlbPredictedMaxCost =
                        max(dlbPredictedDonorCost, dlbPredictedHelperCost);
                }
                else
                {
                    dlbCurrentMaxCost =
                        max(scalar(localCandidateTotal), scalar(peerCandidateTotal));
                    dlbPredictedDonorCost =
                        scalar(localCandidateTotal - offloadCandidateCount);
                    dlbPredictedHelperCost =
                        scalar(peerCandidateTotal)
                      + dlbOffloadRemoteCostFactor*scalar(offloadCandidateCount);
                    dlbPredictedMaxCost =
                        max(dlbPredictedDonorCost, dlbPredictedHelperCost);
                }

                if (dlbPredictedMaxCost*dlbOffloadBudgetSafety >= dlbCurrentMaxCost)
                {
                    dlbGateCancelled = true;
                    forAll(offloadCellLabels, i)
                    {
                        offloadCells[offloadCellLabels[i]] = false;
                    }

                    offloadCellLabels.clear();
                    offloadCandidateCount = 0;
                    offloadParcelCount = 0;
                }
                else if (dlbOffloadPersistent && offloadCellLabels.size() > 0)
                {
                    dlbPersistentDonorCells.setSize(offloadCellLabels.size());
                    forAll(offloadCellLabels, i)
                    {
                        dlbPersistentDonorCells[i] = offloadCellLabels[i];
                    }

                    dlbPersistentDonorPeer = peerProc;
                    dlbPersistentDonorExpire =
                        dlbTimeIndex + dlbOffloadPersistentLeaseSteps - 1;
                }
            }
        }

        if (dlbPersistentDonorActive)
        {
            dlbUsingPersistentLease = true;
            DynamicList<dsmcParcel*> cellParcels;

            forAll(dlbPersistentDonorCells, i)
            {
                const label cellI = dlbPersistentDonorCells[i];

                if (cellI < 0 || cellI >= nCells)
                {
                    continue;
                }

                collectCellParcels(cellI, cellParcels);

                if (cellParcels.size() <= 1 || cloud_.nCandidatesPerCell()[cellI] <= 0)
                {
                    continue;
                }

                offloadCells[cellI] = true;
                offloadCellLabels.append(cellI);
                offloadCandidateCount += cloud_.nCandidatesPerCell()[cellI];
                offloadParcelCount += cellParcels.size();
            }

            dlbSelectedCells = offloadCellLabels.size();
            dlbSelectedCandidatesBeforeGate = offloadCandidateCount;
            dlbSelectedParcels = offloadParcelCount;
            dlbCurrentMaxCost =
                scalar(localCandidateTotal);
            dlbPredictedDonorCost =
                scalar(localCandidateTotal - offloadCandidateCount);
            dlbPredictedHelperCost =
                dlbOffloadRemoteCostFactor*scalar(offloadCandidateCount);
            dlbPredictedMaxCost =
                max(dlbPredictedDonorCost, dlbPredictedHelperCost);
        }

        if (dlbOffloadIntervalHit && !dlbUsingPersistentLease && offloadCellLabels.size() == 0)
        {
            dlbPersistentDonorCells.setSize(0);
            dlbPersistentDonorPeer = -1;
            dlbPersistentDonorExpire = -1;
        }

        if (dlbCompactStaticFields && offloadCandidateCount > 0)
        {
            dlbCompactTaskPayload = true;
        }

        if (dlbOffloadPersistent && offloadCandidateCount > 0)
        {
            dlbPersistentTaskPayload = true;
            dlbPersistentTaskExpire =
                dlbUsingPersistentLease
              ? dlbPersistentDonorExpire
              : dlbTimeIndex + dlbOffloadPersistentLeaseSteps - 1;
        }

        PstreamBuffers taskBufs;

        if (peerProc >= 0)
        {
            const bool multiSend =
                donorHelperList.size() > 1 && offloadCellLabels.size() > 0;

            if (multiSend)
            {
                // Multi-helper: split offload cells among all helpers
                const label cellsPer =
                    max(label(1), offloadCellLabels.size()/donorHelperList.size());

                forAll(donorHelperList, hi)
                {
                    const label helper = donorHelperList[hi];
                    const label start = hi*cellsPer;
                    const label end =
                        (hi == donorHelperList.size()-1)
                      ? offloadCellLabels.size()
                      : min(start + cellsPer, offloadCellLabels.size());

                    if (start >= end) continue;

                    UOPstream os(helper, taskBufs);
                    os << dlbCompactTaskPayload << token::SPACE
                       << dlbPersistentTaskPayload << token::SPACE
                       << dlbPersistentTaskExpire << token::SPACE
                       << (end - start);

                    DynamicList<dsmcParcel*> cellParcels;
                    for (label ci = start; ci < end; ++ci)
                    {
                        const label cellI = offloadCellLabels[ci];
                        const point& cC = mesh.cellCentres()[cellI];
                        collectCellParcels(cellI, cellParcels);
                        os << cellI << cloud_.nCandidatesPerCell()[cellI]
                           << cloud_.sigmaTcRMax()[cellI] << cellParcels.size();
                        labelList subCellIds(cellParcels.size(), 0);
                        forAll(cellParcels, parcelI)
                        {
                            const vector relPos =
                                cellParcels[parcelI]->position() - cC;
                            subCellIds[parcelI] =
                                pos(relPos.x())+2*pos(relPos.y())+4*pos(relPos.z());
                        }
                        os << subCellIds;
                        forAll(cellParcels, parcelI)
                            writeParcelState(os,*cellParcels[parcelI],dlbCompactTaskPayload);
                    }
                }
            }
            else
            {
            UOPstream os(peerProc, taskBufs);
            os  << dlbCompactTaskPayload
                << token::SPACE
                << dlbPersistentTaskPayload
                << token::SPACE
                << dlbPersistentTaskExpire
                << token::SPACE
                << offloadCellLabels.size();

            DynamicList<dsmcParcel*> cellParcels;

            forAll(offloadCellLabels, taskI)
            {
                const label cellI = offloadCellLabels[taskI];
                const point& cC = mesh.cellCentres()[cellI];
                collectCellParcels(cellI, cellParcels);

                os  << cellI
                    << cloud_.nCandidatesPerCell()[cellI]
                    << cloud_.sigmaTcRMax()[cellI]
                    << cellParcels.size();

                labelList subCellIds(cellParcels.size(), 0);

                forAll(cellParcels, parcelI)
                {
                    const dsmcParcel& p = *cellParcels[parcelI];
                    const vector relPos = p.position() - cC;
                    subCellIds[parcelI] =
                        pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());
                }

                os << subCellIds;

                forAll(cellParcels, parcelI)
                {
                    writeParcelState(os, *cellParcels[parcelI], dlbCompactTaskPayload);
                }
            }
            } // end else
        }

        taskBufs.finishedSends();

        // Drain incoming data from ALL ranks to avoid unconsumed-data
        // errors in PstreamBuffers destructor. Only ranks that are
        // legitimate helpers (donorProc matches sender) execute remote
        // collision; other ranks just consume and discard.
        for (label fromRank = 0; fromRank < Pstream::nProcs(); ++fromRank)
        {
            if (fromRank == Pstream::myProcNo()) continue;
            if (!taskBufs.recvDataCount(fromRank)) continue;

            UIPstream is(fromRank, taskBufs);
            is >> dlbIncomingCompactPayload;
            bool incomingPersistentTask = false;
            label incomingPersistentExpire = -1;
            is >> incomingPersistentTask;
            is >> incomingPersistentExpire;
            label nIncomingTasks = 0;
            is >> nIncomingTasks;

            const bool isMyDonor =
                (donorProc >= 0) && (fromRank == donorProc);

            if (dlbOffloadPersistent && incomingPersistentTask && nIncomingTasks > 0)
            {
                dlbPersistentHelperPeer = fromRank;
                dlbPersistentHelperExpire = incomingPersistentExpire;
            }
            else if (dlbOffloadPersistent && dlbOffloadIntervalHit && nIncomingTasks == 0)
            {
                dlbPersistentHelperPeer = -1;
                dlbPersistentHelperExpire = -1;
            }

            // Only execute remote collision if this rank is a helper
            // and the sender is its assigned donor.
            if (nIncomingTasks > 0 && isMyDonor)
            {
            UOPstream os(fromRank, resultBufs);
            os << nIncomingTasks;

            const label helperCell = 0;
            const point helperPosition = mesh.cellCentres()[helperCell];

            for (label taskI = 0; taskI < nIncomingTasks; ++taskI)
            {
                label donorCellI = -1;
                label nCandidates = 0;
                scalar sigmaTcRMaxInitial = SMALL;
                label nParcels = 0;

                is  >> donorCellI
                    >> nCandidates
                    >> sigmaTcRMaxInitial
                    >> nParcels;

                DynamicList<dsmcParcel*> tempParcelPtrs;
                labelList subCellIds;
                tempParcelPtrs.setCapacity(nParcels);
                is >> subCellIds;

                if (subCellIds.size() != nParcels)
                {
                    FatalErrorInFunction
                        << "DLB offload sub-cell count mismatch for cell " << donorCellI
                        << ": expected " << nParcels
                        << ", received " << subCellIds.size()
                        << exit(FatalError);
                }

                for (label parcelI = 0; parcelI < nParcels; ++parcelI)
                {
                    vector U = Zero;
                    scalar RWF = 1.0;
                    scalar ERot = 0.0;
                    label ELevel = 0;
                    label typeId = -1;
                    label newParcel = -1;
                    label classification = 0;
                    labelList vibLevel;

                    readParcelState
                    (
                        is,
                        U,
                        RWF,
                        ERot,
                        ELevel,
                        typeId,
                        newParcel,
                        classification,
                        vibLevel,
                        dlbIncomingCompactPayload
                    );

                    tempParcelPtrs.append
                    (
                        new dsmcParcel
                        (
                            mesh,
                            helperPosition,
                            U,
                            RWF,
                            ERot,
                            ELevel,
                            helperCell,
                            0,
                            0,
                            typeId,
                            newParcel,
                            classification,
                            vibLevel
                        )
                    );
                }

                scalar sigmaTcRMaxUpdated = sigmaTcRMaxInitial;
                boolList dirtyParcels;
                const label accepted =
                    executeRemoteCell
                    (
                        donorCellI,
                        tempParcelPtrs,
                        subCellIds,
                        nCandidates,
                        sigmaTcRMaxInitial,
                        sigmaTcRMaxUpdated,
                        dirtyParcels
                    );

                label nDirtyParcels = 0;
                forAll(dirtyParcels, parcelI)
                {
                    if (dirtyParcels[parcelI])
                    {
                        ++nDirtyParcels;
                    }
                }

                os  << donorCellI
                    << accepted
                    << sigmaTcRMaxUpdated
                    << tempParcelPtrs.size()
                    << nDirtyParcels;

                forAll(tempParcelPtrs, parcelI)
                {
                    if (dirtyParcels[parcelI])
                    {
                        os  << parcelI
                            << token::SPACE;
                        writeParcelState(os, *tempParcelPtrs[parcelI], dlbIncomingCompactPayload);
                    }
                    delete tempParcelPtrs[parcelI];
                }
            }
            } // end if (nIncomingTasks > 0 && isMyDonor)
            else
            {
                // No tasks or not our donor: send 0 response to drain buffer
                UOPstream os(fromRank, resultBufs);
                os << label(0);
            }
        } // end for (fromRank)

        // Defer resultBufs.finishedSends() to after local collision loop
        // so donor's local collision overlaps with helper's remote execution.
        // All ranks in dlbActive must participate in finishedSends (collective).
        dlbResultsPending = true;

        if (dlbOffloadReport && dlbDonorCandidate)
        {
            Pout<< "DLB offload planner summary: local candidates "
                << dlbLocalCandidateTotal << ", peer candidates "
                << dlbPeerCandidateTotal << ", imbalance "
                << dlbCandidateImbalance << ", desired candidates "
                << dlbDesiredOffloadCandidates << ", selected candidates "
                << dlbSelectedCandidatesBeforeGate << ", selected cells "
                << dlbSelectedCells << ", predicted donor cost "
                << dlbPredictedDonorCost << ", predicted helper cost "
                << dlbPredictedHelperCost << ", predicted max cost "
                << dlbPredictedMaxCost << ", current max cost "
                << dlbCurrentMaxCost << ", gate cancelled "
                << Switch(dlbGateCancelled) << ", offloaded candidates "
                << offloadCandidateCount << ", selected parcels "
                << dlbSelectedParcels << ", offloaded cells "
                << offloadCellLabels.size() << ", compact payload "
                << Switch(dlbCompactTaskPayload) << ", persistent "
                << Switch(dlbPersistentTaskPayload) << ", lease reused "
                << Switch(dlbUsingPersistentLease) << ", lease expire "
                << dlbPersistentTaskExpire << endl;
        }

        dlbOffloadWallTime =
            std::chrono::duration<scalar>
            (
                dlb_clock_type::now() - dlbBegin
            ).count();
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

        // Unified RNG interface: FastRng (configurable) or cloud Random (fallback)
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

        if (offloadCells[cellI])
        {
            return 0;
        }

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
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                //label candidateP = rndGen_.position<label>(0, nC - 1);
                label candidateP = rngPos(nC);

                // Declare the second collision candidate
                label candidateQ = -1;

                const List<label>& subCellPs = subCells[whichSubCell[candidateP]];

                const label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        //candidateQ = subCellPs[rndGen_.position<label>(0, nSC - 1)]; OLD
                        candidateQ = subCellPs[rngPos(nSC)];

                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        //candidateQ = rndGen_.position<label>(0, nC - 1); OLD
                        candidateQ = rngPos(nC);

                    } while (candidateP == candidateQ);
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = cloud_.randomLabel(0, nC-1);

                // // Select a possible second collision candidate
                // label candidateQ = rngPos(nC);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = rngPos(nC);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                const label typeIdP = typeIds[candidateP];
                const label typeIdQ = typeIds[candidateQ];
                label chargeP = charges[candidateP];
                label chargeQ = charges[candidateQ];

                if (chargeP == -2)
                {
                    chargeP = cloud_.constProps(typeIdP).charge();
                }

                if (chargeQ == -2)
                {
                    chargeQ = cloud_.constProps(typeIdQ).charge();
                }

                //do not allow electron-electron collisions

                if(!(chargeP == -1 && chargeQ == -1))
                {

                    const scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        *parcelPtrs[candidateP],
                        *parcelPtrs[candidateQ]
                    );


                    // Update the maximum value of sigmaTcR stored, but use the
                    // initial value in the acceptance-rejection criteria because
                    // the number of collision candidates selected was based on this


                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    if ((sigmaTcR/sigmaTcRMax) > rng01())
                    {
                        // chemical reactions

                        // find which reaction model parcel p and q should use
                        const label rMId =
                            cloud_.reactionsActive()
                          ? cloud_.reactions().pairModelAddressing()[typeIdP][typeIdQ]
                          : -1;

                        dsmcParcel& parcelP = *parcelPtrs[candidateP];
                        dsmcParcel& parcelQ = *parcelPtrs[candidateQ];

    //                             Info << " parcelP id: " <<  parcelP.typeId()
    //                                 << " parcelQ id: " << parcelQ.typeId()
    //                                 << " reaction model: " << rMId
    //                                 << endl;

                        if(rMId != -1)
                        {
                            ++reactionHitCount;
                            // try to react molecules
    //                         if(cloud_.reactions().reactions()[rMId]->reactWithLists())
    //                         {
                                // so far for recombination only
    //                                     reactions_.reactions()[rMId]->reaction
    //                                     (
    //                                         parcelP,
    //                                         parcelQ,
    //                                         candidateList,
    //                                         candidateSubList,
    //                                         candidateP,
    //                                         whichSubCell
    //                                     );
    //                         }
    //                         else
    //                         {
                                cloud_.reactions().reactions()[rMId]->reaction
                                (
                                    parcelP,
                                    parcelQ
                                );
    //                         }
                            // if reaction unsuccessful use conventional collision model
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
                        else // if reaction model not found, use conventional collision model
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
        if (cloud_.openmpCollisionStrategy() == "partition")
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

    if (dlbResultsPending)
    {
        resultBufs.finishedSends();

        using dlb_clock_type = std::chrono::steady_clock;
        const auto dlbApplyBegin = dlb_clock_type::now();

        // Drain result data from ALL ranks
        for (label fromRank = 0; fromRank < Pstream::nProcs(); ++fromRank)
        {
            if (fromRank == Pstream::myProcNo()) continue;
            if (!resultBufs.recvDataCount(fromRank)) continue;

        UIPstream is(fromRank, resultBufs);
            label nResults = 0;
            is >> nResults;

            DynamicList<dsmcParcel*> cellParcels;

            for (label resultI = 0; resultI < nResults; ++resultI)
            {
                label cellI = -1;
                label accepted = 0;
                scalar sigmaTcRMaxUpdated = SMALL;
                label nParcels = 0;
                label nDirtyParcels = 0;

                is  >> cellI
                    >> accepted
                    >> sigmaTcRMaxUpdated
                    >> nParcels
                    >> nDirtyParcels;

                collectCellParcels(cellI, cellParcels);

                if (cellParcels.size() != nParcels)
                {
                    FatalErrorInFunction
                        << "DLB offload result size mismatch for cell " << cellI
                        << ": local parcels " << cellParcels.size()
                        << ", returned parcels " << nParcels
                        << exit(FatalError);
                }

                for (label dirtyI = 0; dirtyI < nDirtyParcels; ++dirtyI)
                {
                    label parcelI = -1;
                    is >> parcelI;

                    if (parcelI < 0 || parcelI >= cellParcels.size())
                    {
                        FatalErrorInFunction
                            << "DLB offload dirty parcel index out of range for cell "
                            << cellI << ": index " << parcelI
                            << ", local parcels " << cellParcels.size()
                            << exit(FatalError);
                    }

                    applyParcelState(is, *cellParcels[parcelI], dlbCompactTaskPayload);
                }

                cloud_.sigmaTcRMax()[cellI] =
                    max(cloud_.sigmaTcRMax()[cellI], sigmaTcRMaxUpdated);
                remoteAcceptedCount += accepted;
                ++remoteActiveCellCount;
            }
        }

        dlbOffloadWallTime +=
            std::chrono::duration<scalar>
            (
                dlb_clock_type::now() - dlbApplyBegin
            ).count();
    }

    if (remoteAcceptedCount || remoteActiveCellCount || dlbOffloadWallTime > 0)
    {
        threadAcceptedCounts[0] += remoteAcceptedCount;
        threadActiveCellCounts[0] += remoteActiveCellCount;
        threadWallTimes[0] += dlbOffloadWallTime;
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
        if (cloud_.nCandidatesPerCell()[cellI] > 0)
        {
            const scalar sigma = cloud_.sigmaTcRMax()[cellI];
            ++localCandidateCells;
            localSigmaSum += sigma;
            localSigmaMax = max(localSigmaMax, sigma);
            localSigmaMin = min(localSigmaMin, sigma);
        }
    }
    scalar localRankCollisionWall = 0.0;

    forAll(threadWallTimes, threadI)
    {
        localRankCollisionWall = max(localRankCollisionWall, threadWallTimes[threadI]);
    }

    dlbPrevRankCollisionWall = localRankCollisionWall;
    dlbPrevSecondsPerCandidate =
        collisionCandidates > 0
      ? localRankCollisionWall/scalar(collisionCandidates)
      : 0.0;

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

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
    //             << "    Acceptance rate                 = "
    //             << scalar(collisions)/scalar(collisionCandidates) << nl
                << endl;

            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No collisions" << endl;

            infoCounter_ = 0;
        }
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
