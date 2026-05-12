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
    dsmcDynamicLoadBalancing

Description

\*----------------------------------------------------------------------------*/

#include "dsmcDynamicLoadBalancing.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "wallPolyPatch.H"
#include "dsmcCloud.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "treeDataCell.H"

#include <algorithm>
#include <chrono>

namespace Foam
{

namespace
{
struct redistributedDsmcParcel
{
    point position;
    barycentric coordinates;
    vector U;
    scalar RWF = 1.0;
    scalar ERot = 0.0;
    label ELevel = 0;
    label face = -1;
    label tetFace = -1;
    label tetPt = -1;
    label typeId = -1;
    label newParcel = -1;
    label classification = 0;
    scalar stepFraction = 1.0;
    scalar behind = 0.0;
    label nBehind = 0;
    label origProc = -1;
    label origId = -1;
    labelList vibLevel;
    label stuckFlag = 0;
    scalarField stuckWallTemperature;
    vectorField stuckWallVectors;
};

Ostream& operator<<(Ostream& os, const redistributedDsmcParcel& p)
{
    os  << p.position << token::SPACE
        << p.coordinates << token::SPACE
        << p.U << token::SPACE
        << p.RWF << token::SPACE
        << p.ERot << token::SPACE
        << p.ELevel << token::SPACE
        << p.face << token::SPACE
        << p.tetFace << token::SPACE
        << p.tetPt << token::SPACE
        << p.typeId << token::SPACE
        << p.newParcel << token::SPACE
        << p.classification << token::SPACE
        << p.stepFraction << token::SPACE
        << p.behind << token::SPACE
        << p.nBehind << token::SPACE
        << p.origProc << token::SPACE
        << p.origId << token::SPACE
        << p.vibLevel << token::SPACE
        << p.stuckFlag;

    if (p.stuckFlag)
    {
        os  << token::SPACE << p.stuckWallTemperature
            << token::SPACE << p.stuckWallVectors;
    }

    return os;
}

Istream& operator>>(Istream& is, redistributedDsmcParcel& p)
{
    is  >> p.position
        >> p.coordinates
        >> p.U
        >> p.RWF
        >> p.ERot
        >> p.ELevel
        >> p.face
        >> p.tetFace
        >> p.tetPt
        >> p.typeId
        >> p.newParcel
        >> p.classification
        >> p.stepFraction
        >> p.behind
        >> p.nBehind
        >> p.origProc
        >> p.origId
        >> p.vibLevel
        >> p.stuckFlag;

    if (p.stuckFlag)
    {
        is >> p.stuckWallTemperature >> p.stuckWallVectors;
    }
    else
    {
        p.stuckWallTemperature.clear();
        p.stuckWallVectors.clear();
    }

    return is;
}

bool operator==
(
    const redistributedDsmcParcel& a,
    const redistributedDsmcParcel& b
)
{
    return &a == &b;
}

bool operator!=
(
    const redistributedDsmcParcel& a,
    const redistributedDsmcParcel& b
)
{
    return !(a == b);
}

scalar elapsedSeconds(const std::chrono::steady_clock::time_point& start)
{
    return std::chrono::duration_cast<std::chrono::duration<scalar>>
    (
        std::chrono::steady_clock::now() - start
    ).count();
}

scalar initialSigmaTcRMaxEstimate(const dsmcCloud& cloud)
{
    scalar sigmaTcRMax = SMALL;
    const List<dsmcParcel::constantProperties>& constProps = cloud.constProps();

    forAll(constProps, i)
    {
        const scalar estimate =
            constProps[i].sigmaT()
           *cloud.maxwellianMostProbableSpeed(300.0, constProps[i].mass());

        sigmaTcRMax = max(sigmaTcRMax, estimate);
    }

    return max(sigmaTcRMax, SMALL);
}

redistributedDsmcParcel makeRedistributedDsmcParcel(const dsmcParcel& p)
{
    redistributedDsmcParcel data;
    data.position = p.position();
    data.coordinates = p.coordinates();
    data.U = p.U();
    data.RWF = p.RWF();
    data.ERot = p.ERot();
    data.ELevel = p.ELevel();
    data.face = p.face();
    data.tetFace = p.tetFace();
    data.tetPt = p.tetPt();
    data.typeId = p.typeId();
    data.newParcel = p.newParcel();
    data.classification = p.classification();
    data.stepFraction = p.stepFraction();
    data.behind = p.behind();
    data.nBehind = p.nBehind();
    data.origProc = p.origProc();
    data.origId = p.origId();
    data.vibLevel = p.vibLevel();
    data.stuckFlag = p.isStuck() ? 1 : 0;

    if (data.stuckFlag)
    {
        data.stuckWallTemperature = p.stuck().wallTemperature();
        data.stuckWallVectors = p.stuck().wallVectors();
    }

    return data;
}

scalar lookupBalanceUntilTime(const dictionary& dict)
{
    if (dict.found("balanceUntilTime"))
    {
        return dict.lookupOrDefault<scalar>("balanceUntilTime", VGREAT);
    }

    return dict.lookupOrDefault<scalar>("loadBalancingUntilTime", VGREAT);
}
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar dsmcDynamicLoadBalancing::totalBalanceWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalPrepareWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalReconstructMeshWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalReconstructFieldsWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalDecomposeWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalProcessorMeshSyncWallTime_ = 0.0;
scalar dsmcDynamicLoadBalancing::totalBackupWallTime_ = 0.0;
label dsmcDynamicLoadBalancing::totalBalanceCount_ = 0;
bool dsmcDynamicLoadBalancing::reportTimingEnabled_ = false;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor
dsmcDynamicLoadBalancing::dsmcDynamicLoadBalancing
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud
)
:
    time_(t),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    dsmcLoadBalanceDict_
    (
        IOobject
        (
            "loadBalanceDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    performBalance_(false),
    enableBalancing_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "enableBalancing",
            false
        )
    ),
    reportTiming_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "reportBalancingTiming",
            enableBalancing_
        )
    ),
    inMemoryBalancing_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "inMemoryBalancing",
            false
        )
    ),
    inMemoryMigrationFraction_
    (
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<scalar>
            (
                "inMemoryMigrationFraction",
                0.50
            ),
            scalar(0)
        )
    ),
    inMemoryMinLocalCells_
    (
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<label>
            (
                "inMemoryMinLocalCells",
                1
            ),
            label(1)
        )
    ),
    inMemoryUseTimerGate_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "inMemoryUseTimerGate",
            false
        )
    ),
    inMemoryRequiredGainFactor_
    (
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<scalar>
            (
                "inMemoryRequiredGainFactor",
                1.0
            ),
            scalar(0)
        )
    ),
    inMemoryAssumedBalanceCost_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<scalar>
        (
            "inMemoryAssumedBalanceCost",
            -1.0
        )
    ),
    balanceUntilTime_(lookupBalanceUntilTime(dsmcLoadBalanceDict_)),
    originalEndTime_(time_.time().endTime().value()),
    maxImbalance_(Foam::hyCompat::toScalar
    (
        dsmcLoadBalanceDict_, "maximumAllowableImbalance"
    )),
    balanceCheckInterval_
    (
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<label>
            (
                "balanceCheckInterval",
                1
            ),
            label(1)
        )
    ),
    lastBalanceCheckTimeIndex_(-1),
    balanceCooldownSteps_
    (
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<label>
            (
                "balanceCooldownSteps",
                0
            ),
            label(0)
        )
    ),
    lastBalancePerformTimeIndex_(-1),
    lastCheckCollisionWallTime_(0.0),
    lastCheckMoveWallTime_(0.0),
    lastCheckEvolveProfileCalls_(0),
    limitTimeDirBackups_
    (
        dsmcLoadBalanceDict_.lookupOrDefault<label>
        (
            "limitTimeDirBackups",
            -1
        )
    )
{
    reportTimingEnabled_ = reportTimingEnabled_ || reportTiming_;

    Info<< "DSMC dynamic load balancing: "
        << (enableBalancing_ ? "enabled" : "disabled")
        << ", mode "
        << (inMemoryBalancing_ ? "in-memory" : "external-repartition")
        << ", check interval " << balanceCheckInterval_
        << ", cooldown " << balanceCooldownSteps_
        << ", max imbalance " << maxImbalance_
        << ", timer gate " << (inMemoryUseTimerGate_ ? "on" : "off")
        << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDynamicLoadBalancing::~dsmcDynamicLoadBalancing()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDynamicLoadBalancing::update()
{
    updateProperties();

    if (!enableBalancing_)
    {
        performBalance_ = false;
        return;
    }

    if (!Pstream::parRun())
    {
        return;
    }

    const label currentTimeIndex = time_.timeIndex();

    if
    (
        currentTimeIndex == lastBalanceCheckTimeIndex_
     || (currentTimeIndex % balanceCheckInterval_) != 0
    )
    {
        return;
    }

    lastBalanceCheckTimeIndex_ = currentTimeIndex;

    if (time_.time().value() >= originalEndTime_ - SMALL)
    {
        return;
    }

    if
    (
        balanceCooldownSteps_ > 0
     && lastBalancePerformTimeIndex_ >= 0
     && currentTimeIndex - lastBalancePerformTimeIndex_ < balanceCooldownSteps_
    )
    {
        Info<< "    DLB cooldown active at time " << time_.timeName()
            << " (timeIndex " << currentTimeIndex
            << ", last balance timeIndex " << lastBalancePerformTimeIndex_
            << ", cooldown " << balanceCooldownSteps_ << ")"
            << nl << endl;
        return;
    }

    const scalar& allowableImbalance = maxImbalance_;

    scalar projectedCollisionGain = -1.0;
    scalar estimatedBalanceCost = -1.0;
    bool timerGateAllows = true;

    if (inMemoryBalancing_ && inMemoryUseTimerGate_ && cloud_.evolveProfileEnabled())
    {
        const scalar currentCollisionWallTime = cloud_.evolveCollisionWallTime();
        const label currentEvolveProfileCalls = cloud_.evolveProfileCalls();
        const scalar localCollisionWindow =
            max(currentCollisionWallTime - lastCheckCollisionWallTime_, scalar(0));
        const label localEvolveWindowCalls =
            max(currentEvolveProfileCalls - lastCheckEvolveProfileCalls_, label(0));

        lastCheckCollisionWallTime_ = currentCollisionWallTime;
        lastCheckEvolveProfileCalls_ = currentEvolveProfileCalls;

        if (localEvolveWindowCalls > 0)
        {
            scalarList procCollisionWindows(Pstream::nProcs(), 0.0);
            procCollisionWindows[Pstream::myProcNo()] = localCollisionWindow;
            Pstream::allGatherList(procCollisionWindows);

            scalar maxCollisionWindow = 0.0;
            scalar totalCollisionWindow = 0.0;
            forAll(procCollisionWindows, proci)
            {
                totalCollisionWindow += procCollisionWindows[proci];
                maxCollisionWindow = max(maxCollisionWindow, procCollisionWindows[proci]);
            }

            const scalar meanCollisionWindow =
                totalCollisionWindow/max(scalar(Pstream::nProcs()), scalar(1));

            projectedCollisionGain =
                max(maxCollisionWindow - meanCollisionWindow, scalar(0));

            if (inMemoryAssumedBalanceCost_ > SMALL)
            {
                estimatedBalanceCost = inMemoryAssumedBalanceCost_;
            }
            else if (totalBalanceCount_ > 0)
            {
                estimatedBalanceCost =
                    totalBalanceWallTime_/max(scalar(totalBalanceCount_), scalar(1));
            }

            timerGateAllows =
                estimatedBalanceCost <= SMALL
             || projectedCollisionGain
                > inMemoryRequiredGainFactor_*estimatedBalanceCost;

            Info<< "    DLB timer gate: collision window max/avg gain [s] = "
                << projectedCollisionGain
                << ", estimated migration cost [s] = "
                << estimatedBalanceCost
                << ", required factor = "
                << inMemoryRequiredGainFactor_
                << ", decision = "
                << (timerGateAllows ? "pass" : "skip")
                << nl << endl;
        }
    }

    // Candidates-per-rank imbalance (true computational cost, no MPI barrier bias)
    scalar localCandidates = 0.0;
    forAll(cloud_.nCandidatesPerCell(), celli)
        localCandidates += scalar(cloud_.nCandidatesPerCell()[celli]);

    scalarList procCandidates(Pstream::nProcs(), 0.0);
    procCandidates[Pstream::myProcNo()] = localCandidates;
    Pstream::allGatherList(procCandidates);

    scalar globalCandidates = 0.0;
    forAll(procCandidates, pi) globalCandidates += procCandidates[pi];
    const scalar idealCandidates = globalCandidates/scalar(Pstream::nProcs());

    scalar candImbLocal = mag(localCandidates - idealCandidates);
    Foam::reduce(candImbLocal, maxOp<scalar>());
    scalar maxCandImbalance =
        idealCandidates > SMALL ? candImbLocal/idealCandidates : 0.0;

    // Particle count imbalance for reporting
    scalar nGlobalParticles = cloud_.size();
    Foam::reduce(nGlobalParticles, sumOp<scalar>());
    scalar idealNParticles = scalar(nGlobalParticles)/scalar(Pstream::nProcs());
    scalar nParticles = cloud_.size();
    scalar localImb = mag(nParticles - idealNParticles);
    Foam::reduce(localImb, maxOp<scalar>());
    scalar maxImbalance = localImb/idealNParticles;

    Info<< "    DLB imbalance check at time " << time_.timeName()
        << " (timeIndex " << currentTimeIndex
        << ", interval " << balanceCheckInterval_ << ")" << nl
        << "    Particle imbalance = " << 100*maxImbalance << "%" << nl
        << "    Candidates imbalance = " << 100*maxCandImbalance << "%" << nl
        << "    Rank candidates = " << procCandidates << nl
        << endl;

    if
    (
           enableBalancing_
        && time_.time().value() <= balanceUntilTime_
        && maxCandImbalance > allowableImbalance
        && timerGateAllows
    )
    {
        if (inMemoryBalancing_)
        {
            Info<< "    DLB trigger: in-memory ownership migration" << nl
                << endl;

            performBalance_ = true;
            return;
        }

        Info<< "    DLB trigger: forcing write of current time before "
            << "mesh repartition" << nl << endl;

        writeParticleWeightField();
        time_.writeNow();

        performBalance_ = true;

        originalEndTime_ = time_.time().endTime().value();

        scalar currentTime = time_.time().value();

        time_.setEndTime(currentTime);
    }
    else if
    (
           enableBalancing_
        && time_.time().value() <= balanceUntilTime_
        && maxImbalance > allowableImbalance
        && !timerGateAllows
    )
    {
        Info<< "    DLB trigger skipped: timer gate projected gain [s] = "
            << projectedCollisionGain
            << ", estimated migration cost [s] = "
            << estimatedBalanceCost
            << nl << endl;
    }
}


void dsmcDynamicLoadBalancing::copyPolyMeshToLatestTimeFolder() const
{
    for (label i=0; i<Pstream::nProcs(); i++)
    {
        const word processorName = "processor" + name(i) + "/";
        const word copyPolyMesh =
            word("latest=$(find ")
          + processorName
          + word(" -mindepth 1 -maxdepth 1 -type d | sed 's#.*/##' | ")
          + word("grep -E '^[0-9.+-eE]+$' | sort -g | tail -n 1);")
          + word("starting=$(find ")
          + processorName
          + word(" -mindepth 1 -maxdepth 1 -type d | sed 's#.*/##' | ")
          + word("grep -E '^[0-9.+-eE]+$' | sort -g | head -n 1);")
          + word("if [ -n \"$latest\" ] && [ -n \"$starting\" ] && ")
          + word("[ ! -d ")
          + processorName
          + word("$latest/polyMesh ]; then meshSource=")
          + processorName
          + word("$starting/polyMesh; ")
          + word("if [ ! -d \"$meshSource\" ]; then meshSource=")
          + processorName
          + word("constant/polyMesh; fi; mkdir -p ")
          + processorName
          + word("$latest; cp -r \"$meshSource\" ")
          + processorName
          + word("$latest/; fi");

        Foam::system(copyPolyMesh);
    }
}


void dsmcDynamicLoadBalancing::writeParticleWeightField() const
{
    volScalarField particleWeights
    (
        IOobject
        (
            "dsmcParticleCount",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const label celli = iter().cell();

        if (celli >= 0 && celli < particleWeights.size())
        {
            particleWeights[celli] += 1.0;
        }
    }

    particleWeights.correctBoundaryConditions();
    particleWeights.write();
}


bool dsmcDynamicLoadBalancing::performInMemory()
{
    if (!Pstream::parRun())
    {
        performBalance_ = false;
        return false;
    }

    const label myProc = Pstream::myProcNo();
    const label nProcs = Pstream::nProcs();
    const label oldNCells = mesh_.nCells();

    labelList procParticleCounts(nProcs, 0);
    procParticleCounts[myProc] = cloud_.size();
    Pstream::allGatherList(procParticleCounts);

    label globalParticles = 0;
    forAll(procParticleCounts, proci)
    {
        globalParticles += procParticleCounts[proci];
    }

    const label oldNFaces = mesh_.nFaces();
    labelList cellLoads(oldNCells, 0);
    scalarField cellSigmaTcRMax(cloud_.sigmaTcRMax().primitiveField());
    scalarField cellCollisionRemainder(cloud_.collisionSelectionRemainder());
    labelList cellMigratedFlag(oldNCells, 0);
    scalar prePairPotentialLocal = 0.0;
    scalar preCellSquareLoadLocal = 0.0;
    label preCollisionCellsLocal = 0;

    const auto loadCountStart = std::chrono::steady_clock::now();
    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const dsmcParcel& p = iter();
        const label celli = p.cell();

        if (celli < 0 || celli >= oldNCells)
        {
            continue;
        }

        ++cellLoads[celli];
    }
    scalar loadCountWallTime = elapsedSeconds(loadCountStart);

    // TACF-style cost: particle count as spatial resolution.
    // Per-particle cost = 1.0 (equivalently: balance particle count).
    // Can be augmented with per-rank wall time ratio later.
    scalarList cellTacfCost(oldNCells, 0.0);
    scalar localTacfCost = 0.0;

    forAll(cellLoads, celli)
    {
        const scalar nC = scalar(cellLoads[celli]);
        preCellSquareLoadLocal += nC*nC;
        if (nC > SMALL)
        {
            cellTacfCost[celli] = nC;
            localTacfCost += nC;
        }
        if (cellLoads[celli] >= 2)
        {
            prePairPotentialLocal += 0.5*nC*(nC - 1.0);
            ++preCollisionCellsLocal;
        }
    }

    scalarList procTacfCosts(nProcs, 0.0);
    procTacfCosts[myProc] = localTacfCost;
    Pstream::allGatherList(procTacfCosts);

    scalar globalTacfCost = 0.0;
    forAll(procTacfCosts, proci)
        globalTacfCost += procTacfCosts[proci];

    const scalar idealTacfCost =
        globalTacfCost/max(scalar(nProcs), scalar(1));
    const scalar localCostExcess =
        max(localTacfCost - idealTacfCost, scalar(0));

    const scalar idealParticles =
        scalar(globalParticles)/max(scalar(nProcs), scalar(1));
    // Cap: don't migrate if it would push particle imbalance beyond 15%
    const scalar maxParticleMigration =
        max(scalar(procParticleCounts[myProc]) - 0.85*idealParticles, scalar(0));

    labelList distribution(oldNCells, myProc);
    labelList procSelectedLoad(nProcs, 0);
    label localSelectedCells = 0;
    label localSelectedParticles = 0;
    scalar remainingExcess = inMemoryMigrationFraction_*localCostExcess;
    scalar remainingParticleBudget = maxParticleMigration;

    // Boundary penalty: don't migrate cells that would leave <50% local neighbors.
    // This prevents thin protrusions that cause transfer explosion.
    const labelListList& cellCells = mesh_.cellCells();
    auto boundaryPenalty = [&](const label celli) -> bool
    {
        const labelList& nbrs = cellCells[celli];
        label localNbrs = 0;
        forAll(nbrs, ni)
        {
            const label nbr = nbrs[ni];
            if (nbr >= 0 && nbr < oldNCells && distribution[nbr] == myProc)
                ++localNbrs;
        }
        // Require >= 50% local neighbors after migration
        return scalar(localNbrs) >= 0.5*scalar(nbrs.size());
    };

    const auto selectStart = std::chrono::steady_clock::now();
    const auto selectStart = std::chrono::steady_clock::now();
    if (remainingExcess > SMALL && oldNCells > inMemoryMinLocalCells_)
    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        const labelListList& cellCells = mesh_.cellCells();

        forAll(patches, patchi)
        {
            if (remainingExcess <= SMALL)
            {
                break;
            }

            if (!isA<processorPolyPatch>(patches[patchi]))
            {
                continue;
            }

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);
            const label nbrProc = procPatch.neighbProcNo();

            if
            (
                nbrProc < 0
             || nbrProc >= nProcs
             || procTacfCosts[nbrProc] >= idealTacfCost
            )
            {
                continue;
            }

            scalar nbrCapacity =
                min
                (
                    idealTacfCost - procTacfCosts[nbrProc],
                    remainingExcess
                );

            DynamicList<label> candidates;
            DynamicList<label> frontier;
            boolList seen(oldNCells, false);
            const labelUList& faceCells = procPatch.faceCells();

            forAll(faceCells, faceI)
            {
                const label celli = faceCells[faceI];

                if
                (
                    celli >= 0
                 && celli < oldNCells
                 && distribution[celli] == myProc
                 && cellTacfCost[celli] > SMALL
                )
                {
                    if (!seen[celli])
                    {
                        candidates.append(celli);
                        seen[celli] = true;
                    }
                }
            }

            std::sort
            (
                candidates.begin(),
                candidates.end(),
                [&](const label a, const label b)
                {
                    return cellTacfCost[a] > cellTacfCost[b];
                }
            );

            forAll(candidates, candidateI)
            {
                if
                (
                    nbrCapacity <= SMALL
                 || remainingParticleBudget <= SMALL
                 || oldNCells - localSelectedCells <= inMemoryMinLocalCells_
                )
                {
                    break;
                }

                const label celli = candidates[candidateI];
                if (scalar(cellLoads[celli]) > remainingParticleBudget)
                {
                    continue;
                }
                if (!boundaryPenalty(celli)) continue;
                distribution[celli] = nbrProc;
                cellMigratedFlag[celli] = 1;
                nbrCapacity -= cellTacfCost[celli];
                remainingExcess -= cellTacfCost[celli];
                remainingParticleBudget -= scalar(cellLoads[celli]);
                procSelectedLoad[nbrProc] += cellLoads[celli];
                frontier.append(celli);
                ++localSelectedCells;
                localSelectedParticles += cellLoads[celli];
            }

            while
            (
                   nbrCapacity > SMALL
                && remainingExcess > SMALL
                && remainingParticleBudget > SMALL
                && oldNCells - localSelectedCells > inMemoryMinLocalCells_
                && frontier.size() > 0
            )
            {
                DynamicList<label> ringCandidates;

                forAll(frontier, frontierI)
                {
                    const label seedCell = frontier[frontierI];
                    const labelList& neighbours = cellCells[seedCell];

                    forAll(neighbours, nbrI)
                    {
                        const label celli = neighbours[nbrI];

                        if
                        (
                            celli >= 0
                         && celli < oldNCells
                         && !seen[celli]
                         && distribution[celli] == myProc
                         && cellTacfCost[celli] > SMALL
                        )
                        {
                            ringCandidates.append(celli);
                            seen[celli] = true;
                        }
                    }
                }

                if (ringCandidates.empty())
                {
                    break;
                }

                std::sort
                (
                    ringCandidates.begin(),
                    ringCandidates.end(),
                    [&](const label a, const label b)
                    {
                        return cellTacfCost[a] > cellTacfCost[b];
                    }
                );

                frontier.clear();

                forAll(ringCandidates, candidateI)
                {
                    if
                    (
                        nbrCapacity <= SMALL
                     || remainingExcess <= SMALL
                     || remainingParticleBudget <= SMALL
                     || oldNCells - localSelectedCells <= inMemoryMinLocalCells_
                    )
                    {
                        break;
                    }

                    const label celli = ringCandidates[candidateI];
                    if (scalar(cellLoads[celli]) > remainingParticleBudget)
                        continue;
                    if (!boundaryPenalty(celli))
                        continue;
                    distribution[celli] = nbrProc;
                    cellMigratedFlag[celli] = 1;
                    nbrCapacity -= cellTacfCost[celli];
                    remainingExcess -= cellTacfCost[celli];
                    remainingParticleBudget -= scalar(cellLoads[celli]);
                    procSelectedLoad[nbrProc] += cellLoads[celli];
                    frontier.append(celli);
                    ++localSelectedCells;
                    localSelectedParticles += cellLoads[celli];
                }
            }
        }
    }
    scalar selectWallTime = elapsedSeconds(selectStart);

    label globalSelectedCells = localSelectedCells;
    label globalSelectedParticles = localSelectedParticles;
    reduce(globalSelectedCells, sumOp<label>());
    reduce(globalSelectedParticles, sumOp<label>());

    if (globalSelectedCells == 0)
    {
        if (Pstream::master())
        {
            Info<< "    In-memory DLB: no boundary cells selected for "
                << "local ownership migration" << nl << endl;
        }

        performBalance_ = false;
        return false;
    }

    const auto parcelPackStart = std::chrono::steady_clock::now();

    // Only serialize parcels from outgoing cells; detach staying parcels
    // into a local buffer to avoid full serialize/rebuild cycle.
    List<List<redistributedDsmcParcel>> cellParcels(oldNCells);
    forAll(cellParcels, celli)
    {
        if (distribution[celli] != myProc)
        {
            cellParcels[celli].setSize(cellLoads[celli]);
        }
    }

    labelList cellParcelOffsets(oldNCells, 0);
    label localPackedParcels = 0;
    label localPackedCells = 0;

    // Detach staying parcels from cloud without destroying them
    DynamicList<dsmcParcel*> stayingParcels(cloud_.size());
    DynamicList<label> stayingParcelOldCells(cloud_.size());
    DynamicList<point> stayingParcelPositions(cloud_.size());

    List<dsmcParcel*> allParcels(cloud_.size());
    {
        label idx = 0;
        for (dsmcParcel& p : cloud_)
        {
            allParcels[idx++] = &p;
        }
        allParcels.setSize(idx);
    }

    forAll(allParcels, i)
    {
        dsmcParcel* pPtr = allParcels[i];
        const label celli = pPtr->cell();

        if (celli < 0 || celli >= oldNCells)
        {
            continue;
        }

        if (distribution[celli] == myProc)
        {
            cloud_.remove(pPtr);
            stayingParcels.append(pPtr);
            stayingParcelOldCells.append(celli);
            stayingParcelPositions.append(pPtr->position());
        }
        else
        {
            const label offset = cellParcelOffsets[celli]++;
            cellParcels[celli][offset] = makeRedistributedDsmcParcel(*pPtr);
            ++localPackedParcels;
        }
    }

    forAll(cellLoads, celli)
    {
        if (cellLoads[celli] > 0 && distribution[celli] != myProc)
        {
            ++localPackedCells;
        }
    }
    scalar parcelPackWallTime = elapsedSeconds(parcelPackStart);
    label globalPackedParcels = localPackedParcels;
    label globalPackedCells = localPackedCells;
    reduce(globalPackedParcels, sumOp<label>());
    reduce(globalPackedCells, sumOp<label>());

    const auto balanceStart = std::chrono::steady_clock::now();
    scalar storePositionsWallTime = 0.0;
    scalar meshDistributeWallTime = 0.0;
    scalar parcelDistributeWallTime = 0.0;
    scalar parcelRebuildWallTime = 0.0;
    scalar occupancyRebuildWallTime = 0.0;
    scalar refreshWallTime = 0.0;

    const auto storePositionsStart = std::chrono::steady_clock::now();
    cloud_.storeGlobalPositions();
    cloud_.clearMoveOrderedParcels();
    cloud_.clearMoveAppendedParcels();
    cloud_.clearPendingMoveParcels();
    cloud_.clear();
    storePositionsWallTime = elapsedSeconds(storePositionsStart);

    fvMesh& mesh = const_cast<fvMesh&>(mesh_);
    fvMeshDistribute distributor(mesh);
    const auto meshDistributeStart = std::chrono::steady_clock::now();
    autoPtr<mapDistributePolyMesh> distMap = distributor.distribute(distribution);
    meshDistributeWallTime = elapsedSeconds(meshDistributeStart);

    const auto parcelDistributeStart = std::chrono::steady_clock::now();
    distMap().distributeCellData(cellParcels);
    distMap().distributeCellData(cellSigmaTcRMax);
    distMap().distributeCellData(cellCollisionRemainder);
    distMap().distributeCellData(cellMigratedFlag);
    parcelDistributeWallTime = elapsedSeconds(parcelDistributeStart);

    labelList oldFaceIds(oldNFaces, -1);
    for (label facei = 0; facei < oldNFaces; ++facei)
    {
        oldFaceIds[facei] = facei;
    }
    distMap().distributeFaceData(oldFaceIds);
    labelList oldFaceToNewFace(oldNFaces, -1);
    forAll(oldFaceIds, newFacei)
    {
        const label oldFacei = oldFaceIds[newFacei];

        if (oldFacei >= 0 && oldFacei < oldFaceToNewFace.size())
        {
            oldFaceToNewFace[oldFacei] = newFacei;
        }
    }

    cloud_.resetAfterMeshDistribution();

    if (cellSigmaTcRMax.size() == cloud_.sigmaTcRMax().primitiveField().size())
    {
        cloud_.sigmaTcRMax().primitiveFieldRef() = cellSigmaTcRMax;
        cloud_.sigmaTcRMax().correctBoundaryConditions();
    }

    if (cellCollisionRemainder.size() == cloud_.collisionSelectionRemainder().size())
    {
        cloud_.collisionSelectionRemainder() = cellCollisionRemainder;
    }

    // Rebuild the core tracking geometry/cache explicitly after distribution
    // so relocated parcels are located against the new mesh state.
    (void)mesh_.solutionD();
    (void)mesh_.tetBasePtIs();
    (void)mesh_.oldCellCentres();
    (void)mesh_.cellTree();

    const auto parcelRebuildStart = std::chrono::steady_clock::now();
    label rebuiltParcels = 0;
    label localCellLocateMismatchCount = 0;
    label localExactTopologyRebuildCount = 0;
    label localExactTopologyFallbackCount = 0;
    label localCellTreeLocateCount = 0;
    scalar localExactTopologyMaxPositionError = 0.0;
    label localStayingRemapped = 0;

    // --- Phase A: remap staying parcels using cell/face distribution maps ---
    // Build old-to-new cell index mapping for cells that stayed on this proc.
    const mapDistribute& cMap = distMap().cellMap();
    const labelList& cellSubLocal = cMap.subMap()[myProc];
    const labelList& cellConstructLocal = cMap.constructMap()[myProc];

    labelList oldCellToNewCell(oldNCells, -1);
    forAll(cellSubLocal, i)
    {
        const label oldIdx = cellSubLocal[i];
        const label newIdx = cellConstructLocal[i];
        if (oldIdx >= 0 && oldIdx < oldNCells)
        {
            oldCellToNewCell[oldIdx] = newIdx;
        }
    }

    label localStayingDirectRemap = 0;
    label localStayingRelocated = 0;

    forAll(stayingParcels, i)
    {
        dsmcParcel* pPtr = stayingParcels[i];
        const label oldCelli = stayingParcelOldCells[i];
        const label newCelli = oldCellToNewCell[oldCelli];

        if (newCelli < 0 || newCelli >= mesh_.nCells())
        {
            FatalErrorInFunction
                << "Staying parcel old cell " << oldCelli
                << " mapped to invalid new cell " << newCelli
                << " (nCells=" << mesh_.nCells() << ")"
                << exit(FatalError);
        }

        const point pos = stayingParcelPositions[i];

        // Try direct index remap first (avoids expensive locate)
        const label oldTetFace = pPtr->tetFace();
        const label mappedTetFace =
            (oldTetFace >= 0 && oldTetFace < oldFaceToNewFace.size())
          ? oldFaceToNewFace[oldTetFace]
          : -1;

        bool directRemapOk = false;
        if
        (
            mappedTetFace >= 0
         && mappedTetFace < mesh_.nFaces()
         && pPtr->tetPt() > 0
         && pPtr->tetPt() < mesh_.faces()[mappedTetFace].size() - 1
        )
        {
            pPtr->cell() = newCelli;
            pPtr->tetFace() = mappedTetFace;
            // tetPt unchanged

            const scalar posErr = mag(pPtr->position() - pos);
            if (posErr <= 1e-10 && pPtr->cell() == newCelli)
            {
                directRemapOk = true;
                ++localStayingDirectRemap;
            }
        }

        if (!directRemapOk)
        {
            pPtr->cell() = newCelli;
            pPtr->relocate(pos, newCelli);
            ++localStayingRelocated;
        }

        const label oldFace = pPtr->face();
        if (oldFace >= 0 && oldFace < oldFaceToNewFace.size())
        {
            pPtr->face() = oldFaceToNewFace[oldFace];
        }
        else
        {
            pPtr->face() = -1;
        }

        cloud_.addParticle(pPtr);
        ++localStayingRemapped;
    }

    // --- Phase B: rebuild only incoming parcels from other processors ---
    forAll(cellParcels, celli)
    {
        const List<redistributedDsmcParcel>& parcels = cellParcels[celli];

        forAll(parcels, parcelI)
        {
            const redistributedDsmcParcel& data = parcels[parcelI];
            dsmcParcel* pPtr = nullptr;
            const label mappedTetFace =
                data.tetFace >= 0 && data.tetFace < oldFaceToNewFace.size()
              ? oldFaceToNewFace[data.tetFace]
              : -1;

            if
            (
                mappedTetFace >= 0
             && mappedTetFace < mesh_.nFaces()
             && data.tetPt > 0
             && data.tetPt < mesh_.faces()[mappedTetFace].size() - 1
            )
            {
                pPtr =
                    new dsmcParcel
                    (
                        mesh_,
                        data.coordinates,
                        celli,
                        mappedTetFace,
                        data.tetPt,
                        data.U,
                        data.RWF,
                        data.ERot,
                        data.ELevel,
                        data.typeId,
                        data.newParcel,
                        data.classification,
                        data.vibLevel
                    );

                const scalar positionError = mag(pPtr->position() - data.position);
                localExactTopologyMaxPositionError =
                    max(localExactTopologyMaxPositionError, positionError);

                if
                (
                    pPtr->cell() == celli
                 && positionError <= 1e-12
                )
                {
                    ++localExactTopologyRebuildCount;
                }
                else
                {
                    delete pPtr;
                    pPtr = nullptr;
                    ++localExactTopologyFallbackCount;
                }
            }

            if (!pPtr)
            {
                ++localCellTreeLocateCount;
                const label locatedCell = mesh_.cellTree().findInside(data.position);

                if (locatedCell < 0)
                {
                    FatalErrorInFunction
                        << "Unable to locate redistributed parcel at position "
                        << data.position << " after in-memory ownership migration"
                        << exit(FatalError);
                }

                if (locatedCell != celli)
                {
                    ++localCellLocateMismatchCount;
                }

                pPtr =
                    new dsmcParcel
                    (
                        mesh_,
                        data.position,
                        locatedCell,
                        data.U,
                        data.RWF,
                        data.ERot,
                        data.ELevel,
                        data.typeId,
                        data.newParcel,
                        data.classification,
                        data.vibLevel
                    );
            }

            pPtr->stepFraction() = data.stepFraction;
            pPtr->behind() = data.behind;
            pPtr->nBehind() = data.nBehind;
            pPtr->origProc() = data.origProc;
            pPtr->origId() = data.origId;

            const label mappedFace =
                data.face >= 0 && data.face < oldFaceToNewFace.size()
              ? oldFaceToNewFace[data.face]
              : -1;
            pPtr->face() = mappedFace;

            if (data.stuckFlag)
            {
                pPtr->setStuck
                (
                    data.stuckWallTemperature,
                    data.stuckWallVectors
                );
            }

            cloud_.addParticle(pPtr);
            ++rebuiltParcels;
        }
    }
    parcelRebuildWallTime = elapsedSeconds(parcelRebuildStart);

    const auto occupancyStart = std::chrono::steady_clock::now();
    cloud_.rebuildCellOccupancyAfterMeshDistribution();
    occupancyRebuildWallTime = elapsedSeconds(occupancyStart);

    scalar postPairPotentialLocal = 0.0;
    scalar postCellSquareLoadLocal = 0.0;
    label postCollisionCellsLocal = 0;
    label postOccupancyParcelsLocal = 0;
    label postCloudValidCellParcelsLocal = 0;
    label postCloudInvalidCellParcelsLocal = 0;

    for (label celli = 0; celli < mesh_.nCells(); ++celli)
    {
        const label occCount = cloud_.occupancyCount(celli);
        const scalar nC = scalar(occCount);
        postOccupancyParcelsLocal += occCount;
        postCellSquareLoadLocal += nC*nC;

        if (nC >= 2.0)
        {
            postPairPotentialLocal += 0.5*nC*(nC - 1.0);
            ++postCollisionCellsLocal;
        }
    }

    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const label celli = iter().cell();

        if (celli >= 0 && celli < mesh_.nCells())
        {
            ++postCloudValidCellParcelsLocal;
        }
        else
        {
            ++postCloudInvalidCellParcelsLocal;
        }
    }

    label localRecomputedSigmaCells = 0;
    forAll(cellMigratedFlag, celli)
    {
        if (!cellMigratedFlag[celli])
        {
            continue;
        }

        const DynamicList<dsmcParcel*>& cellParcelsInCell =
            cloud_.cellOccupancy()[celli];
        const label nCellParcels = cellParcelsInCell.size();

        if (nCellParcels < 2)
        {
            continue;
        }

        scalar sigmaTcRMaxCell = SMALL;

        for (label i = 0; i < nCellParcels; ++i)
        {
            const dsmcParcel& pP = *cellParcelsInCell[i];

            for (label j = i + 1; j < nCellParcels; ++j)
            {
                const dsmcParcel& pQ = *cellParcelsInCell[j];
                const scalar sigmaTcR =
                    cloud_.binaryCollision().sigmaTcR(pP, pQ);
                sigmaTcRMaxCell = max(sigmaTcRMaxCell, sigmaTcR);
            }
        }

        cloud_.sigmaTcRMax()[celli] =
            max(sigmaTcRMaxCell, initialSigmaTcRMaxEstimate(cloud_));
        ++localRecomputedSigmaCells;
    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();

    const auto refreshStart = std::chrono::steady_clock::now();
    cloud_.refreshAfterMeshDistribution();
    refreshWallTime = elapsedSeconds(refreshStart);

    scalar balanceWallTime = elapsedSeconds(balanceStart);
    reduce(balanceWallTime, maxOp<scalar>());
    reduce(storePositionsWallTime, maxOp<scalar>());
    reduce(meshDistributeWallTime, maxOp<scalar>());
    reduce(parcelDistributeWallTime, maxOp<scalar>());
    reduce(parcelRebuildWallTime, maxOp<scalar>());
    reduce(occupancyRebuildWallTime, maxOp<scalar>());
    reduce(refreshWallTime, maxOp<scalar>());

    totalBalanceWallTime_ += balanceWallTime;
    totalBalanceCount_++;
    lastBalancePerformTimeIndex_ = time_.timeIndex();

    label localRebuiltParcels = rebuiltParcels;
    label globalRebuiltParcels = localRebuiltParcels;
    reduce(globalRebuiltParcels, sumOp<label>());
    label globalStayingRemapped = localStayingRemapped;
    reduce(globalStayingRemapped, sumOp<label>());
    label globalStayingDirectRemap = localStayingDirectRemap;
    reduce(globalStayingDirectRemap, sumOp<label>());
    label globalStayingRelocated = localStayingRelocated;
    reduce(globalStayingRelocated, sumOp<label>());
    label globalCellLocateMismatchCount = localCellLocateMismatchCount;
    reduce(globalCellLocateMismatchCount, sumOp<label>());
    label globalExactTopologyRebuildCount = localExactTopologyRebuildCount;
    reduce(globalExactTopologyRebuildCount, sumOp<label>());
    label globalExactTopologyFallbackCount = localExactTopologyFallbackCount;
    reduce(globalExactTopologyFallbackCount, sumOp<label>());
    label globalCellTreeLocateCount = localCellTreeLocateCount;
    reduce(globalCellTreeLocateCount, sumOp<label>());
    scalar globalExactTopologyMaxPositionError = localExactTopologyMaxPositionError;
    reduce(globalExactTopologyMaxPositionError, maxOp<scalar>());
    label globalRecomputedSigmaCells = localRecomputedSigmaCells;
    reduce(globalRecomputedSigmaCells, sumOp<label>());
    scalar globalPrePairPotential = prePairPotentialLocal;
    reduce(globalPrePairPotential, sumOp<scalar>());
    scalar globalPostPairPotential = postPairPotentialLocal;
    reduce(globalPostPairPotential, sumOp<scalar>());
    scalar globalPreCellSquareLoad = preCellSquareLoadLocal;
    reduce(globalPreCellSquareLoad, sumOp<scalar>());
    scalar globalPostCellSquareLoad = postCellSquareLoadLocal;
    reduce(globalPostCellSquareLoad, sumOp<scalar>());
    label globalPreCollisionCells = preCollisionCellsLocal;
    reduce(globalPreCollisionCells, sumOp<label>());
    label globalPostCollisionCells = postCollisionCellsLocal;
    reduce(globalPostCollisionCells, sumOp<label>());
    label globalPostOccupancyParcels = postOccupancyParcelsLocal;
    reduce(globalPostOccupancyParcels, sumOp<label>());
    label globalPostCloudValidCellParcels = postCloudValidCellParcelsLocal;
    reduce(globalPostCloudValidCellParcels, sumOp<label>());
    label globalPostCloudInvalidCellParcels = postCloudInvalidCellParcelsLocal;
    reduce(globalPostCloudInvalidCellParcels, sumOp<label>());

    labelList postParticleCounts(nProcs, 0);
    postParticleCounts[myProc] = cloud_.size();
    Pstream::allGatherList(postParticleCounts);

    scalar postGlobalParticles = 0;
    forAll(postParticleCounts, proci)
    {
        postGlobalParticles += postParticleCounts[proci];
    }

    const scalar postIdealParticles =
        scalar(postGlobalParticles)/max(scalar(nProcs), scalar(1));

    scalar postLocalImbalance = mag(scalar(postParticleCounts[myProc]) - postIdealParticles);
    reduce(postLocalImbalance, maxOp<scalar>());
    const scalar postMaxImbalance =
        postLocalImbalance/max(postIdealParticles, scalar(VSMALL));

    Info<< "    In-memory DLB ownership migration:" << nl
        << "        selected cells global      = "
        << globalSelectedCells << nl
        << "        selected parcel load       = "
        << globalSelectedParticles << nl
        << "        staying parcels remapped   = "
        << globalStayingRemapped << nl
        << "        staying direct remap       = "
        << globalStayingDirectRemap << nl
        << "        staying relocated          = "
        << globalStayingRelocated << nl
        << "        rebuilt parcels global     = "
        << globalRebuiltParcels << nl
        << "        parcel cell relabels       = "
        << globalCellLocateMismatchCount << nl
        << "        exact topology rebuilds    = "
        << globalExactTopologyRebuildCount << nl
        << "        exact topology fallbacks   = "
        << globalExactTopologyFallbackCount << nl
        << "        cellTree locates           = "
        << globalCellTreeLocateCount << nl
        << "        exact topology max |dx|    = "
        << globalExactTopologyMaxPositionError << nl
        << "        sigma recomputed cells     = "
        << globalRecomputedSigmaCells << nl
        << "        pre collision cells        = "
        << globalPreCollisionCells << nl
        << "        post collision cells       = "
        << globalPostCollisionCells << nl
        << "        pre pair potential         = "
        << globalPrePairPotential << nl
        << "        post pair potential        = "
        << globalPostPairPotential << nl
        << "        pre square-load sum        = "
        << globalPreCellSquareLoad << nl
        << "        post square-load sum       = "
        << globalPostCellSquareLoad << nl
        << "        post occupancy parcels     = "
        << globalPostOccupancyParcels << nl
        << "        post cloud valid parcels   = "
        << globalPostCloudValidCellParcels << nl
        << "        post cloud invalid parcels = "
        << globalPostCloudInvalidCellParcels << nl
        << "        local cells before/after   = "
        << oldNCells << " / " << mesh_.nCells() << nl
        << "        post particle counts        = "
        << postParticleCounts << nl
        << "        post maximum imbalance      = "
        << 100*postMaxImbalance << "%" << nl
        << "        store/clear wall [s]       = "
        << storePositionsWallTime << nl
        << "        parcel pack/detach wall [s]= "
        << parcelPackWallTime << nl
        << "        mesh distribute wall [s]   = "
        << meshDistributeWallTime << nl
        << "        parcel distribute wall [s] = "
        << parcelDistributeWallTime << nl
        << "        parcel rebuild wall [s]    = "
        << parcelRebuildWallTime << nl
        << "        occupancy rebuild wall [s] = "
        << occupancyRebuildWallTime << nl
        << "        refresh wall [s]           = "
        << refreshWallTime << nl
        << "        wall time [s]              = "
        << balanceWallTime << nl
        << endl;

    performBalance_ = false;
    return true;
}


void dsmcDynamicLoadBalancing::perform(const label noRefinement)
{
    if (enableBalancing_ && performBalance_)
    {
        if (inMemoryBalancing_)
        {
            performInMemory();
            return;
        }

        const auto balanceStart = std::chrono::steady_clock::now();
        scalar prepareWallTime = 0.0;
        scalar reconstructMeshWallTime = 0.0;
        scalar reconstructFieldsWallTime = 0.0;
        scalar decomposeWallTime = 0.0;
        scalar processorMeshSyncWallTime = 0.0;
        scalar backupWallTime = 0.0;

        if (Pstream::master())
        {
            if (noRefinement == 0)
            {
                const auto sectionStart = std::chrono::steady_clock::now();
                copyPolyMeshToLatestTimeFolder();
                prepareWallTime += elapsedSeconds(sectionStart);
            }

            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system("reconstructParMesh -latestTime");
                reconstructMeshWallTime = elapsedSeconds(sectionStart);
            }

            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system("reconstructPar -latestTime");
                reconstructFieldsWallTime = elapsedSeconds(sectionStart);
            }

            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system
                (
                    "latest=$(foamListTimes -latestTime);"
                    "if [ ! -d $latest/polyMesh ]; then cp -r constant/polyMesh $latest/; fi"
                );
                Foam::system("rm -r processor*");
                prepareWallTime += elapsedSeconds(sectionStart);
            }

            const word decomposePar =
                word("timeDirs=$(foamListTimes -noZero);")
                // check if there are any time dirs, if not this indicates a
                // fatal error
                + word("if [ -z \"$timeDirs\" ];")
                + word("then echo \"error: no time directories to decompose\";")
                // decompose the latest time dir
                + word("else decomposeTool=\"${FOAM_USER_APPBIN}/decomposeDSMCLoadBalancePar\";")
                + word("if [ ! -x \"$decomposeTool\" ]; then ")
                + word("decomposeTool=decomposeDSMCLoadBalancePar; fi;")
                + word("\"$decomposeTool\" -force -latestTime -copyUniform;")
                + word("fi");
            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system(decomposePar);
                decomposeWallTime = elapsedSeconds(sectionStart);
            }

            {
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system
                (
                    "latest=$(foamListTimes -processor -latestTime); "
                    "for procDir in processor*; do "
                    "if [ -d \"$procDir/$latest/polyMesh\" ]; then "
                    "mkdir -p \"$procDir/constant\"; "
                    "rm -rf \"$procDir/constant/polyMesh\"; "
                    "cp -r \"$procDir/$latest/polyMesh\" \"$procDir/constant/\"; "
                    "fi; "
                    "done"
                );
                processorMeshSyncWallTime = elapsedSeconds(sectionStart);
            }

            // backup time dirs must be stored in resultFolders to prevent them
            // from being cleared. They can be moved back when the simulation
            // has finished.
            mkDir("resultFolders");

            // respect the specified limit of max. concurrent time dir backups
            if (limitTimeDirBackups_ >= 0)
            {
                // impose the limit by moving back the time dirs that are
                // currently already backuped up and then only keeping the
                // most recent time dirs <= the imposed limit.
                const word backupTimeDirsWithLimit =
                    // move the time dirs currently backed up to the case dir
                    // so the foamListTimes utility can be used
                    // make sure directory is not empty to prevent mv from
                    // printing a warning
                    word("if [ \"$(ls -A resultFolders)\" ];")
                    + word("then mv resultFolders/* .; fi;")
                    // total number of time directories
                    + word("timeDirs=$(foamListTimes);")
                    + word("nTimeDirs=$(echo $timeDirs | tr -cd ' ' | wc -c);")
                    + word("nTimeDirs=$((nTimeDirs+1));")
                    // convert OpenFOAM label to shell variable for limit
                    + word("limitNTimeDirs=")
                    + name(limitTimeDirBackups_)
                    + word(";")
                    // if the total number of time directories is larger than
                    // limit remove the oldest time directories first
                    + word("diffNTimeDirs=$((nTimeDirs-limitNTimeDirs));")
                    + word("if [ $diffNTimeDirs -gt 0 ];")
                    + word("then timeDirs=$(echo $timeDirs | ")
                    + word("cut -d ' ' -f$((diffNTimeDirs+1))-$nTimeDirs);")
                    + word("fi;")
                    // move the time dirs that were chosen to be eligible for
                    // backup to the backup dir
                    + word("mv $timeDirs resultFolders/;")
                    // clear all other time dirs
                    + word("foamListTimes -rm");
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system(backupTimeDirsWithLimit);
                backupWallTime = elapsedSeconds(sectionStart);
            }
            else
            {
                // keep all time dirs in the backup dir
                const auto sectionStart = std::chrono::steady_clock::now();
                Foam::system
                (
                    "timeDirs=$(foamListTimes); mv $timeDirs resultFolders/"
                );
                backupWallTime = elapsedSeconds(sectionStart);
            }

            performBalance_ = false;
        }

        reduce(prepareWallTime, maxOp<scalar>());
        reduce(reconstructMeshWallTime, maxOp<scalar>());
        reduce(reconstructFieldsWallTime, maxOp<scalar>());
        reduce(decomposeWallTime, maxOp<scalar>());
        reduce(processorMeshSyncWallTime, maxOp<scalar>());
        reduce(backupWallTime, maxOp<scalar>());

        const scalar balanceWallTime = elapsedSeconds(balanceStart);

        totalBalanceWallTime_ += balanceWallTime;
        totalPrepareWallTime_ += prepareWallTime;
        totalReconstructMeshWallTime_ += reconstructMeshWallTime;
        totalReconstructFieldsWallTime_ += reconstructFieldsWallTime;
        totalDecomposeWallTime_ += decomposeWallTime;
        totalProcessorMeshSyncWallTime_ += processorMeshSyncWallTime;
        totalBackupWallTime_ += backupWallTime;
        totalBalanceCount_++;

        if (reportTiming_ && Pstream::master())
        {
            Info<< "DLB wall-time summary [s]:" << nl
                << "    total              = " << balanceWallTime << nl
                << "    prepare            = " << prepareWallTime << nl
                << "    reconstruct mesh   = " << reconstructMeshWallTime << nl
                << "    reconstruct fields = " << reconstructFieldsWallTime << nl
                << "    decompose          = " << decomposeWallTime << nl
                << "    processor mesh sync= " << processorMeshSyncWallTime << nl
                << "    backup             = " << backupWallTime << nl
                << endl;
        }

        time_.setEndTime(originalEndTime_);
    }
}


void dsmcDynamicLoadBalancing::resetTiming()
{
    totalBalanceWallTime_ = 0.0;
    totalPrepareWallTime_ = 0.0;
    totalReconstructMeshWallTime_ = 0.0;
    totalReconstructFieldsWallTime_ = 0.0;
    totalDecomposeWallTime_ = 0.0;
    totalProcessorMeshSyncWallTime_ = 0.0;
    totalBackupWallTime_ = 0.0;
    totalBalanceCount_ = 0;
    reportTimingEnabled_ = false;
}

void dsmcDynamicLoadBalancing::updateProperties()
{
    enableBalancing_ =
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "enableBalancing",
            false
        );
    reportTiming_ =
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "reportBalancingTiming",
            enableBalancing_
        );
    reportTimingEnabled_ = reportTimingEnabled_ || reportTiming_;

    inMemoryBalancing_ =
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "inMemoryBalancing",
            false
        );
    inMemoryMigrationFraction_ =
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<scalar>
            (
                "inMemoryMigrationFraction",
                0.50
            ),
            scalar(0)
        );
    inMemoryMinLocalCells_ =
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<label>
            (
                "inMemoryMinLocalCells",
                1
            ),
            label(1)
        );
    inMemoryUseTimerGate_ =
        dsmcLoadBalanceDict_.lookupOrDefault<Switch>
        (
            "inMemoryUseTimerGate",
            false
        );
    inMemoryRequiredGainFactor_ =
        max
        (
            dsmcLoadBalanceDict_.lookupOrDefault<scalar>
            (
                "inMemoryRequiredGainFactor",
                1.0
            ),
            scalar(0)
        );
    inMemoryAssumedBalanceCost_ =
        dsmcLoadBalanceDict_.lookupOrDefault<scalar>
        (
            "inMemoryAssumedBalanceCost",
            -1.0
        );

    // if balancing is active this additional option allows to specify a time
    // after which balancing is deactivated (this is useful in conjunction with
    // resetAtOutput / resetAtOutputUntilTime and averaging across solver
    // restarts)
    balanceUntilTime_ = lookupBalanceUntilTime(dsmcLoadBalanceDict_);
    maxImbalance_ = Foam::hyCompat::toScalar
    (
        dsmcLoadBalanceDict_, "maximumAllowableImbalance"
    );
    balanceCheckInterval_ = max
    (
        dsmcLoadBalanceDict_.lookupOrDefault<label>
        (
            "balanceCheckInterval",
            1
        ),
        label(1)
    );
    balanceCooldownSteps_ = max
    (
        dsmcLoadBalanceDict_.lookupOrDefault<label>
        (
            "balanceCooldownSteps",
            0
        ),
        label(0)
    );
    limitTimeDirBackups_ = dsmcLoadBalanceDict_.lookupOrDefault<label>
    (
        "limitTimeDirBackups",
        -1
    );
}


}  // End namespace Foam

// ************************************************************************* //
