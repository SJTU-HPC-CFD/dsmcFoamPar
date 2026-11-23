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
#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <queue>
#include <vector>
#include <mpi.h>
#include <malloc.h>

using namespace Foam::constant;

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<dsmcParcel>, 0);
};


namespace
{
    typedef std::chrono::steady_clock steadyWallClock;

    thread_local Foam::dsmcFastRng* collisionRngContext = nullptr;

    // buildCellOccupancy subphase wall accumulators (§6.4 diagnostics,
    // worklog dlb_further_opt.md).  OMP path only; serial fallback leaves
    // these unchanged.  Reported (max over ranks) by printProfileSummary.
    Foam::scalar buildOccResetWall = 0.0;
    Foam::scalar buildOccCountWall = 0.0;
    Foam::scalar buildOccReduceWall = 0.0;
    Foam::scalar buildOccFillWall = 0.0;
    Foam::label buildOccOrderedHits = 0;
    Foam::label buildOccFallbacks = 0;

    // Overflow-safe cumulative move-detail counters (worklog §11).  The
    // per-step trackingData label counters cannot overflow within one step,
    // but run-long label accumulators wrap at 13M-cell scale (200 steps x
    // 1.5e8 parcels), so the by-rank table reads these double accumulators
    // instead (exact to 2^53).  Zero header change (file-static pattern).
    // Index order: parcels, trackCalls, sameTet, internalTet, faceHits,
    // procHits, cyclicHits, patchHits, stuckHits, trackWall, trackerWall,
    // boundaryWall.
    double dsmcMoveDetailCum[12] = {0};

    // §15.4: cumulative collision/candidate counters (hooked from the
    // collision partner selection once per collide call; read by the
    // per-rank load tables at nTerminalOutputs and post-DLB cadences)
    double dsmcCumCollisions_ = 0;
    double dsmcCumCollisionCandidates_ = 0;

    // §28: evolve-residual attribution (per-step model info updates).
    double dsmcEvolveUpdateInfoWall_ = 0.0;
    double dsmcEvolveCoordWall_ = 0.0;
    double dsmcEvolveCtrlWall_ = 0.0;

    // §11.6 nearWall/interior track accounting.  Cumulative doubles
    // (wall-cell track calls, wall-cell track wall, interior track calls,
    // interior track wall) plus the lazily built nearWall flag table.
    double dsmcMoveDetailSplitCum[4] = {0};
    bool dsmcNearWallBuilt = false;
    const Foam::fvMesh* dsmcNearWallMeshPtr = nullptr;
    Foam::boolList dsmcNearWallCells;

    void buildNearWallCells(const Foam::fvMesh& mesh)
    {
        if (dsmcNearWallBuilt && dsmcNearWallMeshPtr == &mesh)
        {
            return;
        }

        Foam::boolList& flags = dsmcNearWallCells;
        flags.setSize(mesh.nCells());
        flags = false;

        const Foam::labelList& owner = mesh.faceOwner();
        const Foam::polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
        forAll(pbMesh, patchi)
        {
            const Foam::polyPatch& pp = pbMesh[patchi];
            for (Foam::label fi = pp.start(); fi < pp.start() + pp.size(); ++fi)
            {
                flags[owner[fi]] = true;
            }
        }

        dsmcNearWallMeshPtr = &mesh;
        dsmcNearWallBuilt = true;
    }


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

    Foam::IOobject::readOption dsmcSigmaReadOption
    (
        const Foam::fvMesh& mesh
    )
    {
        const Foam::dictionary& controlDict = mesh.time().controlDict();
        const bool processorRestart =
            controlDict.lookupOrDefault<bool>("replicatedMesh", false)
         && controlDict.lookupOrDefault<Foam::word>
            ("replicatedMeshWriteMode", "gathered") == "processor"
         && controlDict.lookupOrDefault<Foam::word>
            ("startFrom", "latestTime") == "latestTime";

        if (processorRestart)
        {
            const Foam::fileName checkpoint =
                mesh.time().path()
               /"processor0"
               /mesh.time().timeName()
               /"dsmcSigmaTcRMax";

        if (Foam::isFile(checkpoint))
        {
            // Checkpoint restart: value constructor starts from zero and
            // restoreProcessorCheckpoint() fills it in explicitly.
            return Foam::IOobject::READ_IF_PRESENT;
        }
    }

    // Fresh start: the value constructor calls readIfPresent() internally,
    // which only performs a real read for READ_IF_PRESENT (MUST_READ would
    // just emit a warning and leave the field zeroed).
    return Foam::IOobject::READ_IF_PRESENT;
}


// Tier-2 initial-read filtering (see doc/worklog/v2506/dsmcFoam_todolist.md):
// for replicated-mesh fresh starts (startTime branch, positions file
// present) read only the parcels owned by this rank; the checkpoint
// restart path, single-rank runs and non-replicated runs are excluded.
// Controlled by replicatedMeshFilterInitialRead (default off).
static bool dsmcFilteredInitialRead
(
    const Foam::fvMesh& mesh,
    const Foam::word& cloudName
)
{
    const Foam::dictionary& controlDict = mesh.time().controlDict();

    if (!controlDict.lookupOrDefault<bool>("replicatedMesh", false))
    {
        return false;
    }

    if
    (
        !controlDict.lookupOrDefault<bool>
        (
            "replicatedMeshFilterInitialRead",
            false
        )
    )
    {
        return false;
    }

    const Foam::word startFrom = controlDict.lookupOrDefault<Foam::word>
    (
        "startFrom",
        "startTime"
    );

    if (startFrom == "latestTime")
    {
        // Checkpoint restart reads through restoreProcessorCheckpoint()
        return false;
    }

    // NB: no Pstream::nProcs() check here - without the -parallel flag
    // Pstream still reports 1 at this point (MPI is initialised later in
    // dsmcReplicatedMesh::initialize()).  For a genuine single-rank run
    // the filter keeps everything (cellOwner all zero == myRank), which
    // is harmless.

    // Same path resolution as IOPosition::headerOk (top-level positions
    // file only)
    const Foam::fileName posFile
    (
        mesh.time().path()
       /mesh.time().timeName()
       /"lagrangian"
       /cloudName
       /"positions"
    );

    return Foam::isFile(posFile);
}
}


// In-solver parallel fill (optimisation 3b): for replicated-mesh fresh
// starts without an initial cloud file, fill owned cells directly from
// dsmcInitialiseDict (zero file round-trip).  Off by default.
static bool dsmcInSolverFill
(
    const Foam::fvMesh& mesh,
    const Foam::word& cloudName
)
{
    const Foam::dictionary& controlDict = mesh.time().controlDict();

    if (!controlDict.lookupOrDefault<bool>("replicatedMesh", false))
    {
        return false;
    }

    if (!controlDict.lookupOrDefault<bool>("dsmcInSolverFill", false))
    {
        return false;
    }

    const Foam::word startFrom = controlDict.lookupOrDefault<Foam::word>
    (
        "startFrom",
        "startTime"
    );

    if (startFrom == "latestTime")
    {
        // Checkpoint restart reads through restoreProcessorCheckpoint()
        return false;
    }

    // Only when there is no initial cloud file to read
    const Foam::fileName posFile
    (
        mesh.time().path()
       /mesh.time().timeName()
       /"lagrangian"
       /cloudName
       /"positions"
    );

    return !Foam::isFile(posFile);
}


void Foam::dsmcCloud::initialFillFromDict(const IOdictionary& dsmcInitialiseDict)
{
    dsmcAllConfigurations conf(dsmcInitialiseDict, *this);
    conf.setInitialConfigParallel
    (
        replicatedMesh().cellOwner(),
        replicatedMesh().myRank()
    );
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


void Foam::dsmcCloud::printInitialiseTimeStepCalculation
(
    Time& runTime,
    const IOdictionary& dsmcInitialiseDict
) const
{
    if (!dsmcInitialiseDict.found("configurations"))
    {
        WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
            << "No configurations entry found in dsmcInitialiseDict. "
            << "Time step calculation skipped." << endl;
        return;
    }

    const PtrList<entry> configurations
    (
        dsmcInitialiseDict.lookup("configurations")
    );

    if (configurations.size() == 0)
    {
        WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
            << "No configurations found in dsmcInitialiseDict. "
            << "Time step calculation skipped." << endl;
        return;
    }

    const dictionary& configurationDict = configurations[0].dict();

    if
    (
       !configurationDict.found("translationalTemperature")
     || !configurationDict.found("velocity")
     || !configurationDict.found("numberDensities")
    )
    {
        WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
            << "The first dsmcInitialiseDict configuration must contain "
            << "translationalTemperature, velocity and numberDensities. "
            << "Time step calculation skipped." << endl;
        return;
    }

    const scalar temperature =
        readScalar(configurationDict.lookup("translationalTemperature"));
    const vector flowVelocity(configurationDict.lookup("velocity"));
    const scalar flowSpeed = mag(flowVelocity);
    const dictionary& numberDensitiesDict =
        configurationDict.subDict("numberDensities");
    const dictionary& moleculePropertiesDict =
        particleProperties_.subDict("moleculeProperties");

    scalar totalNumberDensity = 0.0;
    scalar totalMassDensity = 0.0;

    forAllConstIter(dictionary, numberDensitiesDict, iter)
    {
        const word& speciesName = iter().keyword();
        const scalar numberDensity =
            readScalar(numberDensitiesDict.lookup(speciesName));

        totalNumberDensity += numberDensity;

        if (moleculePropertiesDict.found(speciesName))
        {
            const dictionary& speciesDict =
                moleculePropertiesDict.subDict(speciesName);
            const scalar molecularMass =
                readScalar(speciesDict.lookup("mass"));

            totalMassDensity += numberDensity*molecularMass;

            if (Pstream::master())
            {
                Info<< "  Species " << speciesName
                    << ": mass = " << molecularMass
                    << " kg, numberDensity = " << numberDensity
                    << " m^-3" << endl;
            }
        }
        else if (Pstream::master())
        {
            WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
                << "Species '" << speciesName
                << "' not found in moleculeProperties. "
                << "Skipping this species in time step calculation." << endl;
        }
    }

    if (totalNumberDensity <= VSMALL || totalMassDensity <= VSMALL)
    {
        WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
            << "Invalid mixture density for time-step calculation. "
            << "totalNumberDensity = " << totalNumberDensity
            << ", totalMassDensity = " << totalMassDensity
            << ". Time step calculation skipped." << endl;
        return;
    }

    const scalar averageMolecularMass =
        totalMassDensity/totalNumberDensity;
    const scalar vThermal =
        maxwellianMostProbableSpeed(temperature, averageMolecularMass);
    const scalar characteristicSpeed = flowSpeed + vThermal;

    if (vThermal <= VSMALL || characteristicSpeed <= VSMALL)
    {
        WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
            << "Invalid molecular speed for time-step calculation. "
            << "vThermal = " << vThermal
            << ", characteristicSpeed = " << characteristicSpeed
            << ". Time step calculation skipped." << endl;
        return;
    }

    const boundBox localBounds(mesh_.bounds());
    scalar minX = localBounds.min().x();
    scalar minY = localBounds.min().y();
    scalar minZ = localBounds.min().z();
    scalar maxX = localBounds.max().x();
    scalar maxY = localBounds.max().y();
    scalar maxZ = localBounds.max().z();

    const scalarField& cellVolumes = mesh_.cellVolumes();
    scalar totalVolume = 0.0;
    scalar minVolume = GREAT;
    scalar maxVolume = -GREAT;

    forAll(cellVolumes, cellI)
    {
        const scalar volume = cellVolumes[cellI];
        totalVolume += volume;
        minVolume = min(minVolume, volume);
        maxVolume = max(maxVolume, volume);
    }

    label totalCells = mesh_.nCells();

    if (Pstream::parRun())
    {
        reduce(minX, minOp<scalar>());
        reduce(minY, minOp<scalar>());
        reduce(minZ, minOp<scalar>());
        reduce(maxX, maxOp<scalar>());
        reduce(maxY, maxOp<scalar>());
        reduce(maxZ, maxOp<scalar>());
        reduce(totalVolume, sumOp<scalar>());
        reduce(minVolume, minOp<scalar>());
        reduce(maxVolume, maxOp<scalar>());
        reduce(totalCells, sumOp<label>());
    }

    if (totalCells <= 0 || minVolume <= VSMALL || maxVolume <= VSMALL)
    {
        WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
            << "Invalid mesh statistics for time-step calculation. "
            << "totalCells = " << totalCells
            << ", minVolume = " << minVolume
            << ", maxVolume = " << maxVolume
            << ". Time step calculation skipped." << endl;
        return;
    }

    const scalar xLen = maxX - minX;
    const scalar yLen = maxY - minY;
    const scalar zLen = maxZ - minZ;
    const scalar maxDomainLength = max(xLen, max(yLen, zLen));
    const scalar meanCellVolume = totalVolume/scalar(totalCells);
    const scalar characteristicLength = pow(meanCellVolume, 1.0/3.0);

    const scalar CAC =
        controlDict_.lookupOrDefault<scalar>("CAC", 1.0);
    const scalar CTC =
        controlDict_.lookupOrDefault<scalar>("CTC", 1.0);
    const label sampleSteps =
        controlDict_.lookupOrDefault<label>("sampleSteps", 500);

    if (CTC <= VSMALL)
    {
        WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
            << "Invalid CTC = " << CTC
            << ". Time step calculation skipped." << endl;
        return;
    }

    const scalar preCtcDeltaT = characteristicLength/characteristicSpeed;
    const scalar newDeltaT = preCtcDeltaT*CTC;

    if (newDeltaT <= VSMALL)
    {
        WarningIn("Foam::dsmcCloud::printInitialiseTimeStepCalculation")
            << "Calculated invalid time step: " << newDeltaT
            << ". Time step calculation skipped." << endl;
        return;
    }

    const label naver =
        label(CAC*maxDomainLength/vThermal/newDeltaT + 1.0);
    const scalar steadyStateTime = naver*newDeltaT;
    const scalar totalTime = steadyStateTime + sampleSteps*newDeltaT;

    if (Pstream::master())
    {
        Info<< nl << "========== DSMC Time Step Calculation =========="
            << endl;
        Info<< "Input parameters:" << endl;
        Info<< "  Translational temperature: " << temperature << " K"
            << endl;
        Info<< "  Flow velocity: " << flowSpeed << " m/s" << endl;
        Info<< "  Total number density: " << totalNumberDensity
            << " m^-3" << endl;
        Info<< "Domain dimensions:" << endl;
        Info<< "  x-length: " << xLen << " m" << endl;
        Info<< "  y-length: " << yLen << " m" << endl;
        Info<< "  z-length: " << zLen << " m" << endl;
        Info<< "  Max domain length: " << maxDomainLength << " m"
            << endl;
        Info<< "  Total number of cells: " << totalCells << endl;
        Info<< "Calculated values:" << endl;
        Info<< "  Average molecular mass: " << averageMolecularMass
            << " kg" << endl;
        Info<< "  Thermal velocity: " << vThermal << " m/s" << endl;
        Info<< "  Characteristic speed: " << characteristicSpeed
            << " m/s" << endl;
        Info<< "  Mean cell volume: " << meanCellVolume << " m^3"
            << endl;
        Info<< "  Max cell volume: " << maxVolume << " m^3" << endl;
        Info<< "  Min cell volume: " << minVolume << " m^3" << endl;
        Info<< "  Max/Min: " << maxVolume/minVolume << endl;
        Info<< "  Characteristic length: " << characteristicLength
            << " m" << endl;
        Info<< "Time step adjustment:" << endl;
        Info<< "  Original time step: " << runTime.deltaTValue()
            << " s" << endl;
        Info<< "  New time step: " << newDeltaT
            << " s before CTC-" << CTC
            << " modification is " << preCtcDeltaT << endl;
        Info<< "  Steady state steps (naver): " << naver << endl;
        Info<< "  Sample steps: " << sampleSteps << endl;
        Info<< "  Total steps: " << naver + sampleSteps << endl;
        Info<< "------------------------------------------------" << endl;
        Info<< "For control file parameter replace:" << endl;
        Info<< "  deltaT             " << newDeltaT << endl;
        Info<< "  steadyStateTime    " << steadyStateTime << endl;
        Info<< "  totalTime          " << totalTime << endl;
        Info<< "================================================"
            << nl << endl;
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
                    appendedThreadOffsets[threadI] = label
                    (
                        (
                            static_cast<long long>(threadI)
                          * static_cast<long long>(appendedParcels)
                        )
                      / static_cast<long long>(ompNumThreads_)
                    );
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
                    generatedThreadOffsets[threadI] = label
                    (
                        (
                            static_cast<long long>(threadI)
                          * static_cast<long long>(gatheredParcels.size())
                        )
                      / static_cast<long long>(ompNumThreads_)
                    );
                }

            generatedThreadOffsets[ompNumThreads_] = gatheredParcels.size();
            parcelsPtr = &gatheredParcels;
            parcelThreadOffsetsPtr = &generatedThreadOffsets;
        }

        const List<dsmcParcel*>& parcels = *parcelsPtr;
        const labelList& parcelThreadOffsets = *parcelThreadOffsetsPtr;
        const bool useAppendedParcels = appendedParcelsPtr != nullptr;

        if (useMoveOrderedParcels)
        {
            ++buildOccOrderedHits;
        }
        else
        {
            ++buildOccFallbacks;
        }

        // ---- subphase: per-thread count reset ----
        const auto tBuildOccReset0 = steadyWallClock::now();

        if (occupancyThreadCellCounts_.size() != ompNumThreads_)
        {
            occupancyThreadCellCounts_.setSize(ompNumThreads_);
        }

        if (occupancyThreadActiveCells_.size() != ompNumThreads_)
        {
            occupancyThreadActiveCells_.setSize(ompNumThreads_);
        }

        // Per-thread reset is thread-disjoint: each thread clears only its
        // own lists (§8.5).
        #pragma omp parallel for num_threads(ompNumThreads_) schedule(static)
        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
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

        buildOccResetWall += elapsedWallSeconds(tBuildOccReset0);

        // ---- subphase: parallel counting ----
        const auto tBuildOccCount0 = steadyWallClock::now();

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

        buildOccCountWall += elapsedWallSeconds(tBuildOccCount0);

        // ---- subphase: thread reduce + active/collision cell lists ----
        const auto tBuildOccReduce0 = steadyWallClock::now();

        // Persistent total counts (no per-step 52 MB allocation) + parallel
        // zero and per-thread atomic merge: the per-thread active-cell lists
        // are disjoint, cells may be shared across threads (§8.5).
        // §46 P1 verdict: a sparse clear (saved last-active set) plus a
        // parallelised classification measured SLOWER than this serial
        // form (reduce 2.75 s vs 1.20 s per 100 steps, same node) — the
        // barrier/atomic overhead exceeds the two serial scans at
        // nCells = 5.8M.  Reverted; see worklog §46.5.
        labelList& totalCounts = occupancyTotalCounts_;
        if (totalCounts.size() != nCells)
        {
            totalCounts.setSize(nCells, 0);
        }

        #pragma omp parallel num_threads(ompNumThreads_)
        {
            #pragma omp for schedule(static)
            for (label cellI = 0; cellI < nCells; ++cellI)
            {
                totalCounts[cellI] = 0;
            }

            #pragma omp for schedule(static)
            for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
            {
                const labelList& localCounts = occupancyThreadCellCounts_[threadI];
                const DynamicList<label>& activeCells =
                    occupancyThreadActiveCells_[threadI];

                forAll(activeCells, activeI)
                {
                    const label cellI = activeCells[activeI];
                    #pragma omp atomic
                    totalCounts[cellI] += localCounts[cellI];
                }
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

        // Build owned collision cells for replicated mesh
        if (replicatedMeshActive())
        {
            DynamicList<label> ownedCC(occupancyCollisionCells_.size());
            forAll(occupancyCollisionCells_, i)
            {
                const label cellI = occupancyCollisionCells_[i];
                if (replicatedMesh_->isMyCell(cellI))
                {
                    ownedCC.append(cellI);
                }
            }
            occupancyOwnedCollisionCells_.transfer(ownedCC);
        }
        else
        {
            occupancyOwnedCollisionCells_ = occupancyCollisionCells_;
        }

        buildOccReduceWall += elapsedWallSeconds(tBuildOccReduce0);

        // ---- subphase: offsets + flat ordered fill ----
        const auto tBuildOccFill0 = steadyWallClock::now();

        occupancyCellOffsets_.setSize(nCells + 1, 0);

        // Chunked parallel exclusive prefix sum over nCells (§8.5): the
        // serial 13M-cell pass was a measurable fill-subphase cost.
        const label nPrefixThreads = max(label(1), min(ompNumThreads_, nCells));
        const label chunkSize = (nCells + nPrefixThreads - 1)/nPrefixThreads;
        labelList chunkStart(nPrefixThreads, 0);
        #pragma omp parallel num_threads(nPrefixThreads)
        {
            #pragma omp for schedule(static)
            for (label c = 0; c < nPrefixThreads; ++c)
            {
                const label b = c*chunkSize;
                const label e = min(b + chunkSize, nCells);
                label s = 0;
                for (label i = b; i < e; ++i)
                {
                    s += totalCounts[i];
                }
                chunkStart[c] = s;
            }

            #pragma omp single
            {
                label run = 0;
                for (label c = 0; c < nPrefixThreads; ++c)
                {
                    const label v = chunkStart[c];
                    chunkStart[c] = run;
                    run += v;
                }
            }

            #pragma omp for schedule(static)
            for (label c = 0; c < nPrefixThreads; ++c)
            {
                const label b = c*chunkSize;
                const label e = min(b + chunkSize, nCells);
                label acc = chunkStart[c];
                for (label i = b; i < e; ++i)
                {
                    occupancyCellOffsets_[i + 1] = acc + totalCounts[i];
                    acc += totalCounts[i];
                }
            }
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

        buildOccFillWall += elapsedWallSeconds(tBuildOccFill0);

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

    // Build owned collision cells for replicated mesh
    if (replicatedMeshActive())
    {
        DynamicList<label> ownedCC(occupancyCollisionCells_.size());
        forAll(occupancyCollisionCells_, i)
        {
            const label cellI = occupancyCollisionCells_[i];
            if (replicatedMesh_->isMyCell(cellI))
            {
                ownedCC.append(cellI);
            }
        }
        occupancyOwnedCollisionCells_.transfer(ownedCC);
    }
    else
    {
        occupancyOwnedCollisionCells_ = occupancyCollisionCells_;
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



// comm_further_opt.md §2.3: occupancy "already current" flag for the
// post-migration build guard, and incremental-append switch caches.
// clean M6: real internal linkage (the previous definitions sat outside the
// anonymous namespace, so "File-static" was wrong and the symbols were
// exported).  Single cloud instance per process; no cross-TU consumer.
static bool dsmcOverlapOccCurrent_ = false;
static bool dsmcOverlapThresholdRead_ = false;
static Foam::scalar dsmcOverlapThresholdCache_ = 0.1;
static bool dsmcOverlapValidateRead_ = false;
static bool dsmcOverlapValidateCache_ = false;

// H4: cumulative wall time of the incremental-occupancy validation (only
// advanced while the debug validation is enabled); reported separately in
// the §35 migration stage table.
static double dsmcOccValidateWall_ = 0.0;

// H4 coverage counters: which input shapes the incremental append actually
// saw (empty receive, one affected cell, several affected cells, and a
// receive touching the last cell index).  Reported at the end of the run so
// the coverage claim is measured, not assumed.
static double dsmcOccAppendEmpty_ = 0.0;
static double dsmcOccAppendSingleCell_ = 0.0;
static double dsmcOccAppendMultiCell_ = 0.0;
static double dsmcOccAppendLastCell_ = 0.0;

// H4 self-test: while true, the debug full-rebuild validation is skipped —
// the synthetic received parcels are not part of the cloud, so a reference
// rebuild could not match them.
static bool dsmcOccSelfTestActive_ = false;

// ---- reactionOccupancyAppend / shared append core (O phase 2) ----
// Per-call scratch written by appendParcelsToOccupancy() and read by the
// reaction path (the async-migration path ignores it).
struct DsmcOccAppendLast
{
    double tSort = 0.0;
    double tShift = 0.0;
    double tOffsets = 0.0;
    double tDerived = 0.0;
    Foam::label nAffected = 0;
    Foam::label maxPerCell = 0;
    Foam::label nKeptOne = 0;
    Foam::label nKeptZero = 0;
    bool touchesLast = false;
};
static DsmcOccAppendLast dsmcOccLast_;

// Reaction-path counters (file-static, zero header change for the data).
struct DsmcRxnAppendStats
{
    double nCalls = 0;          // post-collision updates attempted (reaction steps)
    double nAppend = 0;         // successful incremental appends
    double nFbThreshold = 0;    // fallback: new fraction above threshold
    double nFbGuard = 0;        // fallback: a guard condition failed
    double nSkipZero = 0;       // reactions but no new parcel: nothing to do
    double nLate = 0;           // flat view invalidated after controlAfterCollisions
    double nParcels = 0;        // cumulative appended parcels
    double nShapeSingle = 0;    // appends touching exactly one cell
    double nShapeMulti = 0;     // appends touching several cells
    double nShapeSameCell = 0;  // appends with >=2 new parcels in one cell
    double nKeptOneCells = 0;   // cells that had 1 parcel (1 -> >=2 transition)
    double nKeptZeroCells = 0;  // cells that were empty (0 -> >=1 activation)
    double nShapeLast = 0;      // appends touching the last cell index
    double tAppend = 0.0;       // wall in append path (validation excluded)
    double tSort = 0.0;
    double tShift = 0.0;
    double tOffsets = 0.0;
    double tDerived = 0.0;
    double tValidate = 0.0;     // reference rebuild + comparison (separate)
    double nValidated = 0;
};
static DsmcRxnAppendStats dsmcRxn_;
static bool dsmcRxnAppendRead_ = false;
static bool dsmcRxnAppendEnabled_ = false;
static bool dsmcRxnAppendValidate_ = false;
static bool dsmcRxnAppendHasElectrons_ = false;

static void dsmcReadOverlapControls(const Foam::dictionary& controlDict)
{
    if (!dsmcOverlapThresholdRead_)
    {
        dsmcOverlapThresholdCache_ =
            controlDict.lookupOrDefault<Foam::scalar>
            (
                "replicatedMeshOverlapBuildThreshold",
                0.1
            );
        dsmcOverlapValidateCache_ =
            controlDict.lookupOrDefault<bool>
            (
                "replicatedMeshOverlapBuildValidate",
                false
            );
        if (dsmcOverlapValidateCache_ && Foam::Pstream::master())
        {
            Foam::Info<< "Replicated mesh: incremental occupancy validation"
                << " enabled (debug: full reference rebuild every append)"
                << Foam::endl;
        }
        dsmcOverlapThresholdRead_ = true;
    }
}

// Core of the incremental occupancy update (shared by the async-migration
// path and the post-reaction path).  Precondition (NOT checked here, it is
// the caller's guard): occupancyOrderedParcels_/occupancyCellOffsets_ hold
// exactly the cloud minus 'received'.  Returns false, with the flat state
// untouched, when the received fraction exceeds the threshold; the caller
// decides how to fall back (it must call buildCellOccupancy() itself and must
// NOT use markCurrent unless it is the async path).  Sets no validity flag.
bool Foam::dsmcCloud::appendParcelsToOccupancy
(
    const UList<dsmcParcel*>& received,
    const bool requireOwned
)
{
    typedef std::chrono::steady_clock clk;
    dsmcReadOverlapControls(controlDict_);
    dsmcOccLast_ = DsmcOccAppendLast();

    if (received.size() == 0)
    {
        return true;
    }

    const label nCells = mesh_.nCells();

    // ---- 0. Input validation (H4) ----
    // Runs before every consumer, including the threshold fallback below:
    // buildCellOccupancy() also classifies these parcels by cell, so an
    // out-of-range cell would corrupt memory there as well.  The check is
    // O(nRecv) and always on; its cost is accumulated in the validation
    // timer so it is never mixed into the incremental-append gain.
    {
        const auto tInput0 = clk::now();
        label nBadCell = 0;
        label firstBadCell = -1;
        label nForeign = 0;
        label firstForeign = -1;
        forAll(received, i)
        {
            const label c = received[i]->cell();
            if (c < 0 || c >= nCells)
            {
                if (nBadCell == 0)
                {
                    firstBadCell = c;
                }
                ++nBadCell;
            }
            else if
            (
                requireOwned
             && replicatedMeshActive()
             && !replicatedMesh_->isMyCell(c)
            )
            {
                if (nForeign == 0)
                {
                    firstForeign = c;
                }
                ++nForeign;
            }
        }
        if (nBadCell != 0 || nForeign != 0)
        {
            FatalErrorInFunction
                << "Incremental occupancy append received invalid parcels:"
                << " rank "
                << (replicatedMeshActive() ? replicatedMesh_->myRank() : 0)
                << " nRecv " << received.size()
                << " nCells " << nCells
                << " badCell " << nBadCell << " (first " << firstBadCell
                << ") foreignOwner " << nForeign << " (first "
                << firstForeign << ")"
                << abort(FatalError);
        }
        dsmcOccValidateWall_ += std::chrono::duration<double>
        (
            clk::now() - tInput0
        ).count();
    }

    // Fallback condition: received fraction too large (DLB repartition
    // windows).  The caller does the full rebuild.
    if (received.size() > dsmcOverlapThresholdCache_*this->size())
    {
        return false;
    }

    const auto tStep1 = clk::now();

    // ---- 1. Received per-cell counts (O(nRecv)) ----
    labelList recvCells(received.size());
    forAll(received, i)
    {
        recvCells[i] = received[i]->cell();
    }
    std::sort(recvCells.begin(), recvCells.end());

    label nAffected = 0;
    for (label i = 0; i < recvCells.size(); ++i)
    {
        if (i == 0 || recvCells[i] != recvCells[i-1])
        {
            ++nAffected;
        }
    }
    labelList affectedCells(nAffected);
    labelList affectedRecvCounts(nAffected);
    {
        label a = -1;
        for (label i = 0; i < recvCells.size(); ++i)
        {
            if (i == 0 || recvCells[i] != recvCells[i-1])
            {
                ++a;
                affectedCells[a] = recvCells[i];
                affectedRecvCounts[a] = 1;
            }
            else
            {
                ++affectedRecvCounts[a];
            }
        }
    }

    // Old offsets are needed only at the affected cells: capture them now
    // (replaces the former full copy of occupancyCellOffsets_ and lets the
    // offsets array be updated in place in step 4).
    //   oldStart[a] = offsets[cell_a],  oldAt[a] = offsets[cell_a + 1]
    labelList oldStart(nAffected);
    labelList oldAt(nAffected);
    for (label a = 0; a < nAffected; ++a)
    {
        oldStart[a] = occupancyCellOffsets_[affectedCells[a]];
        oldAt[a] = occupancyCellOffsets_[affectedCells[a] + 1];
    }
    const label oldTotal = occupancyCellOffsets_[nCells];
    const label recvTotal = received.size();

    // Shape of this receive (H4 coverage / reaction-path coverage).
    label maxPerCell = 0;
    label nKeptOne = 0;
    label nKeptZero = 0;
    for (label a = 0; a < nAffected; ++a)
    {
        maxPerCell = max(maxPerCell, affectedRecvCounts[a]);
        const label kept = oldAt[a] - oldStart[a];
        if (kept == 1)
        {
            ++nKeptOne;
        }
        else if (kept == 0)
        {
            ++nKeptZero;
        }
    }

    const auto tStep2 = clk::now();

    // ---- 2. Shift kept segments to their new positions ----
    // Block formulation: block a = cells (affected[a-1], affected[a]]
    // (block 0 = cells [0, affected[0]]; top block = cells above the
    // last affected).  Every cell in block a moves by the SAME shift
    // S(a) = sum(recv[a'] for a' < a): its own received parcels go after
    // its kept segment and do not move its start.  Process the top block
    // first, then affected blocks DESCENDING, so each memmove
    // destination never clobbers data that still needs moving (dest end
    // = srcEnd + S(a) <= the next-higher block's final start).
    occupancyOrderedParcels_.resize(oldTotal + recvTotal);

    {
        // Top block: cells (affectedCells[last], nCells) move by
        // recvTotal.
        if (nAffected > 0 && recvTotal > 0)
        {
            const label srcBegin = oldAt[nAffected - 1];
            const label srcEnd = oldTotal;
            if (srcEnd > srcBegin)
            {
                std::memmove
                (
                    occupancyOrderedParcels_.data() + srcBegin + recvTotal,
                    occupancyOrderedParcels_.data() + srcBegin,
                    sizeof(dsmcParcel*)*label(srcEnd - srcBegin)
                );
            }
        }

        // Affected blocks, descending.  shiftAcc holds
        // sum(recv[a'] for a' >= a) on entry of iteration a.
        label shiftAcc = recvTotal;
        for (label a = nAffected - 1; a >= 0; --a)
        {
            shiftAcc -= affectedRecvCounts[a];
            const label blockShift = shiftAcc;  // = S(a)
            if (blockShift == 0) continue;
            const label srcBegin = a > 0 ? oldAt[a - 1] : 0;
            const label srcEnd = oldAt[a];
            if (srcEnd > srcBegin)
            {
                std::memmove
                (
                    occupancyOrderedParcels_.data()
                  + srcBegin + blockShift,
                    occupancyOrderedParcels_.data() + srcBegin,
                    sizeof(dsmcParcel*)*label(srcEnd - srcBegin)
                );
            }
        }
        // Cells below the first affected cell: shift 0, no move.
    }

    // ---- 3. Place received parcels at their cell segment tails ----
    {
        // cursors[a] = new start of cell_a + kept count
        //            = oldStart[a] + S(a) + (oldAt[a] - oldStart[a])
        //            = oldAt[a] + S(a)
        labelList cursors(nAffected);
        label shift = 0;
        for (label a = 0; a < nAffected; ++a)
        {
            cursors[a] = oldAt[a] + shift;
            shift += affectedRecvCounts[a];
        }
        forAll(received, i)
        {
            const label cellI = received[i]->cell();
            const label a =
                std::lower_bound
                (
                    affectedCells.begin(),
                    affectedCells.end(),
                    cellI
                ) - affectedCells.begin();
            occupancyOrderedParcels_[cursors[a]++] = received[i];
        }
    }

    const auto tStep4 = clk::now();

    // ---- 4. Update the offsets array in place ----
    // offsets[c+1] += (number of received parcels in cells <= c).  Segment
    // a = cells [affected[a], affected[a+1]) all gain cum[a] = sum of the
    // first a+1 receive counts; cells below affected[0] are unchanged.
    {
        labelList cum(nAffected);
        label run = 0;
        for (label a = 0; a < nAffected; ++a)
        {
            run += affectedRecvCounts[a];
            cum[a] = run;
        }

        labelList& off = occupancyCellOffsets_;

        #ifdef _OPENMP
        #pragma omp parallel for num_threads(ompNumThreads_) \
            schedule(dynamic, 16) if (ompNumThreads_ > 1 && nAffected > 64)
        #endif
        for (label a = 0; a < nAffected; ++a)
        {
            const label cBegin = affectedCells[a];
            const label cEnd =
                (a + 1 < nAffected) ? affectedCells[a + 1] : nCells;
            const label add = cum[a];
            for (label c = cBegin; c < cEnd; ++c)
            {
                off[c + 1] += add;
            }
        }
    }

    const auto tStep5 = clk::now();

    // ---- 5. Derived lists (active / collision / owned) ----
    {
        label activeCellCount = 0;
        label collisionCellCount = 0;
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            const label count =
                occupancyCellOffsets_[cellI + 1]
              - occupancyCellOffsets_[cellI];
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
            const label count =
                occupancyCellOffsets_[cellI + 1]
              - occupancyCellOffsets_[cellI];
            if (count > 0)
            {
                occupancyActiveCells_[activeCellCount++] = cellI;
                if (count > 1)
                {
                    occupancyCollisionCells_[collisionCellCount++] = cellI;
                }
            }
        }

        if (replicatedMeshActive())
        {
            DynamicList<label> ownedCC(occupancyCollisionCells_.size());
            forAll(occupancyCollisionCells_, i)
            {
                const label cellI = occupancyCollisionCells_[i];
                if (replicatedMesh_->isMyCell(cellI))
                {
                    ownedCC.append(cellI);
                }
            }
            occupancyOwnedCollisionCells_.transfer(ownedCC);
        }
        else
        {
            occupancyOwnedCollisionCells_ = occupancyCollisionCells_;
        }
    }

    const auto tStep6 = clk::now();

    // ---- 6. Per-cell DynamicList occupancy ----
    if (cellOccupancyMaterialized_)
    {
        for (label a = 0; a < nAffected; ++a)
        {
            const label cellI = affectedCells[a];
            DynamicList<dsmcParcel*>& cellParcels = cellOccupancy_[cellI];
            cellParcels.clear();
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

    auto sec = [](const clk::time_point& a, const clk::time_point& b)
    {
        return std::chrono::duration<double>(b - a).count();
    };
    dsmcOccLast_.tSort = sec(tStep1, tStep2);
    dsmcOccLast_.tShift = sec(tStep2, tStep4);
    dsmcOccLast_.tOffsets = sec(tStep4, tStep5);
    dsmcOccLast_.tDerived = sec(tStep5, tStep6);
    dsmcOccLast_.nAffected = nAffected;
    dsmcOccLast_.maxPerCell = maxPerCell;
    dsmcOccLast_.nKeptOne = nKeptOne;
    dsmcOccLast_.nKeptZero = nKeptZero;
    dsmcOccLast_.touchesLast = (affectedCells[nAffected - 1] == nCells - 1);

    return true;
}
bool Foam::dsmcCloud::appendReceivedToOccupancy
(
    const List<dsmcParcel*>& received
)
{
    // Set the "current" flag on every exit path so the caller's
    // post-migration buildCellOccupancy() is skipped exactly once, and
    // restore the ordered-view validity: the occupancy data IS the
    // complete post-migration state on every exit path (no-received
    // early exit keeps the gap build; fallback rebuilds; incremental
    // updates).  Without this, collide's lazy materialize would clear
    // the per-cell lists and skip the rebuild (worklog
    // comm_further_opt.md §6 root cause).
    auto markCurrent = [&](bool ok)
    {
        dsmcOverlapOccCurrent_ = true;
        occupancyOrderedParcelsValid_ = true;
        return ok;
    };

    dsmcReadOverlapControls(controlDict_);

    // Nothing received: the gap build over kept is already complete.
    if (received.size() == 0)
    {
        dsmcOccAppendEmpty_ += 1.0;
        return markCurrent(true);
    }

    const label nCells = mesh_.nCells();

    // Input validation (H4), threshold fallback and the incremental update
    // itself live in the shared core.  NOTE: the entry validity is the
    // caller's gapBuildDone gate (passed via the mesh), NOT
    // occupancyOrderedParcelsValid_: the moveOrdered set/clear in
    // migrateFinish invalidates that flag while the ordered DATA from
    // the gap build remains intact.
    if (!appendParcelsToOccupancy(received, true))
    {
        buildCellOccupancy();
        return markCurrent(false);
    }

    // H4 coverage: record the shape of this receive.
    if (dsmcOccLast_.nAffected == 1)
    {
        dsmcOccAppendSingleCell_ += 1.0;
    }
    else
    {
        dsmcOccAppendMultiCell_ += 1.0;
    }
    if (dsmcOccLast_.touchesLast)
    {
        dsmcOccAppendLastCell_ += 1.0;
    }
    // ---- 7. Validation mode: recompute the reference and compare ----
    const auto tOccValidate0 = std::chrono::steady_clock::now();
    if (dsmcOverlapValidateCache_ && !dsmcOccSelfTestActive_)
    {
        labelList refCounts(nCells, 0);
        forAllIter(dsmcCloud, *this, iter)
        {
            const label cellI = iter().cell();
            if (cellI >= 0 && cellI < nCells)
            {
                ++refCounts[cellI];
            }
        }
        labelList refOff(nCells + 1, 0);
        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            refOff[cellI + 1] = refOff[cellI] + refCounts[cellI];
        }
        List<dsmcParcel*> refOrdered(refOff[nCells]);
        {
            labelList cursor(refOff);
            forAllIter(dsmcCloud, *this, iter)
            {
                const label cellI = iter().cell();
                if (cellI >= 0 && cellI < nCells)
                {
                    refOrdered[cursor[cellI]++] = &iter();
                }
            }
        }

        label nMismatch = 0;
        if (refOff[nCells] != occupancyCellOffsets_[nCells])
        {
            ++nMismatch;
        }
        for
        (
            label cellI = 0;
            cellI <= nCells && nMismatch == 0;
            ++cellI
        )
        {
            if (refOff[cellI] != occupancyCellOffsets_[cellI])
            {
                ++nMismatch;
                break;
            }
        }
        if (nMismatch == 0)
        {
            for
            (
                label cellI = 0;
                cellI < nCells && nMismatch == 0;
                ++cellI
            )
            {
                const label b = refOff[cellI];
                const label e = refOff[cellI + 1];
                if
                (
                    e - b
                 != occupancyCellOffsets_[cellI + 1]
                  - occupancyCellOffsets_[cellI]
                )
                {
                    ++nMismatch;
                    break;
                }
                // Membership scan (order-free); validation mode only
                for (label i = b; i < e; ++i)
                {
                    const dsmcParcel* p = refOrdered[i];
                    bool found = false;
                    for
                    (
                        label j = occupancyCellOffsets_[cellI];
                        j < occupancyCellOffsets_[cellI + 1];
                        ++j
                    )
                    {
                        if (occupancyOrderedParcels_[j] == p)
                        {
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        ++nMismatch;
                        break;
                    }
                }
            }
        }
        if (nMismatch != 0)
        {
            FatalErrorInFunction
                << "Incremental occupancy append validation failed"
                << " (mismatch count " << nMismatch << ")"
                << " at time " << mesh_.time().timeName()
                << abort(FatalError);
        }
    }
    dsmcOccValidateWall_ += std::chrono::duration<double>
    (
        std::chrono::steady_clock::now() - tOccValidate0
    ).count();

    static label nIncrementalRuns_ = 0;
    if (nIncrementalRuns_ == 0 && Pstream::master())
    {
        Info<< "Replicated mesh: overlap incremental occupancy active"
            << " (first append: " << received.size() << " parcels)"
            << endl;
    }
    ++nIncrementalRuns_;

    return markCurrent(true);
}

// ---- reactionOccupancyAppend (O phase 2): post-reaction occupancy update ----

bool Foam::dsmcCloud::reactionOccupancyAppendEnabled() const
{
    if (!dsmcRxnAppendRead_)
    {
        dsmcRxnAppendEnabled_ =
            controlDict_.lookupOrDefault<bool>("reactionOccupancyAppend", false);
        dsmcRxnAppendValidate_ =
            controlDict_.lookupOrDefault<bool>
            (
                "reactionOccupancyAppendValidate",
                false
            );
        dsmcRxnAppendHasElectrons_ = findIndex(typeIdList_, word("e-")) != -1;
        dsmcRxnAppendRead_ = true;
        if (dsmcRxnAppendEnabled_ && isOutputRank())
        {
            Info<< "Reaction occupancy append: enabled"
                << (dsmcRxnAppendValidate_ ? " + VALIDATE (debug)" : "")
                << (dsmcRxnAppendHasElectrons_
                  ? " (e- species present: always falls back)" : "")
                << endl;
        }
    }
    return dsmcRxnAppendEnabled_;
}


void Foam::dsmcCloud::noteReactionOccupancyLateInvalidate()
{
    dsmcRxn_.nLate += 1.0;
}


namespace
{
    void dsmcRxnValidateFail(const std::string& msg)
    {
        FatalErrorInFunction
            << "reactionOccupancyAppend validation failed: " << msg.c_str()
            << Foam::abort(Foam::FatalError);
    }
}


void Foam::dsmcCloud::validateReactionOccupancy
(
    const char* outcome,
    const label appendedSizeBefore
)
{
    // Independent reference: the flat occupancy view is recomputed from the
    // cloud's IDLList (ground truth, no dependence on moveOrdered/appended)
    // and compared with the state the update left.  Nothing the solver
    // reads is modified - the state left behind is the APPEND result (or
    // the fallback's full-rebuild result), never the reference, so the
    // following consumers are exercised on the state under test.
    const auto t0 = std::chrono::steady_clock::now();
    const label nCells = mesh_.nCells();
    const label nParcels = this->size();
    const std::string tag = std::string("[") + outcome + "] time "
      + mesh_.time().timeName().c_str() + ": ";

    // ---- exit flags ----
    if (!occupancyOrderedParcelsValid_)
    {
        dsmcRxnValidateFail(tag + "occupancyOrderedParcelsValid_ is false");
    }
    if (dsmcOverlapOccCurrent_)
    {
        dsmcRxnValidateFail(tag + "dsmcOverlapOccCurrent_ was set");
    }
    if (moveAppendedParcels_.size() != appendedSizeBefore)
    {
        dsmcRxnValidateFail
        (
            tag + "moveAppendedParcels_ size changed "
          + std::to_string(appendedSizeBefore) + " -> "
          + std::to_string(moveAppendedParcels_.size())
        );
    }
    if (occupancyCellOffsets_.size() != nCells + 1)
    {
        dsmcRxnValidateFail(tag + "offsets size != nCells+1");
    }
    if (occupancyCellOffsets_.last() != nParcels)
    {
        dsmcRxnValidateFail
        (
            tag + "offsets.last()=" + std::to_string(occupancyCellOffsets_.last())
          + " != cloud size " + std::to_string(nParcels)
        );
    }
    if (occupancyOrderedParcels_.size() != nParcels)
    {
        dsmcRxnValidateFail(tag + "orderedParcels size != cloud size");
    }

    // ---- reference counts / offsets ----
    labelList refCount(nCells, 0);
    forAllIter(dsmcCloud, *this, iter)
    {
        const label c = iter().cell();
        if (c >= 0 && c < nCells)
        {
            ++refCount[c];
        }
    }
    labelList refOff(nCells + 1, 0);
    for (label c = 0; c < nCells; ++c)
    {
        refOff[c + 1] = refOff[c] + refCount[c];
    }
    if (occupancyCellOffsets_[0] != 0)
    {
        dsmcRxnValidateFail(tag + "offsets[0] != 0");
    }
    for (label c = 0; c <= nCells; ++c)
    {
        if (refOff[c] != occupancyCellOffsets_[c])
        {
            // first differing offset: report the cell and counts
            const label cc = (c > 0 ? c - 1 : 0);
            dsmcRxnValidateFail
            (
                tag + "offsets[" + std::to_string(c) + "] ref "
              + std::to_string(refOff[c]) + " got "
              + std::to_string(occupancyCellOffsets_[c]) + " (first differing"
                " cell ~" + std::to_string(cc) + ", ref count "
              + std::to_string(refCount[cc]) + ")"
            );
        }
    }

    // ---- per-cell pointer multiset (sorted comparison, O(n log n)) ----
    std::vector<dsmcParcel*> refOrdered(nParcels);
    {
        labelList cursor(refOff);
        forAllIter(dsmcCloud, *this, iter)
        {
            const label c = iter().cell();
            if (c >= 0 && c < nCells)
            {
                refOrdered[cursor[c]++] = &iter();
            }
        }
    }
    std::vector<dsmcParcel*> tmp;
    for (label c = 0; c < nCells; ++c)
    {
        const label b = refOff[c];
        const label e = refOff[c + 1];
        if (e - b == 0)
        {
            continue;
        }
        std::sort(refOrdered.begin() + b, refOrdered.begin() + e);
        tmp.assign
        (
            occupancyOrderedParcels_.begin() + b,
            occupancyOrderedParcels_.begin() + e
        );
        std::sort(tmp.begin(), tmp.end());
        for (label i = 0; i < e - b; ++i)
        {
            if (tmp[i] != refOrdered[b + i])
            {
                dsmcRxnValidateFail
                (
                    tag + "parcel multiset differs in cell " + std::to_string(c)
                  + " (count " + std::to_string(e - b) + ", position "
                  + std::to_string(i) + ")"
                );
            }
        }
    }

    // ---- derived lists ----
    DynamicList<label> act;
    DynamicList<label> col;
    for (label c = 0; c < nCells; ++c)
    {
        if (refCount[c] > 0)
        {
            act.append(c);
            if (refCount[c] > 1)
            {
                col.append(c);
            }
        }
    }
    DynamicList<label> own;
    forAll(col, i)
    {
        if (!replicatedMeshActive() || replicatedMesh_->isMyCell(col[i]))
        {
            own.append(col[i]);
        }
    }
    auto cmpList = [&](const char* name, const labelList& got, const UList<label>& ref)
    {
        if (got.size() != ref.size())
        {
            dsmcRxnValidateFail
            (
                tag + name + " size ref " + std::to_string(ref.size())
              + " got " + std::to_string(got.size())
            );
        }
        forAll(got, i)
        {
            if (got[i] != ref[i])
            {
                dsmcRxnValidateFail
                (
                    tag + name + "[" + std::to_string(i) + "] ref "
                  + std::to_string(ref[i]) + " got " + std::to_string(got[i])
                );
            }
        }
    };
    cmpList("occupancyActiveCells", occupancyActiveCells_, act);
    cmpList("occupancyCollisionCells", occupancyCollisionCells_, col);
    cmpList("occupancyOwnedCollisionCells", occupancyOwnedCollisionCells_, own);

    dsmcRxn_.tValidate += std::chrono::duration<double>
    (
        std::chrono::steady_clock::now() - t0
    ).count();
    dsmcRxn_.nValidated += 1.0;
}


bool Foam::dsmcCloud::updateOccupancyAfterReactions
(
    const label preSize,
    const label mark,
    const bool preValid
)
{
    if (!reactionOccupancyAppendEnabled())
    {
        buildCellOccupancy();
        return false;
    }

    typedef std::chrono::steady_clock clk;
    const auto tEntry = clk::now();
    const scalar cpu0 = profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;
    const label nCells = mesh_.nCells();
    const label appSize0 = moveAppendedParcels_.size();
    const label nNew = appSize0 - mark;
    dsmcRxn_.nCalls += 1.0;

    // Guard (audit_1.3.md §5.2.3 / §5.2.5): every condition must hold or the
    // full rebuild runs.
    const bool structural =
        openmpEnabled_
     && openmpMoveEnabled_
     && ompNumThreads_ > 1
     && !dsmcRxnAppendHasElectrons_;
    const bool guard =
        structural
     && preValid
     && moveOrderedParcelsValid_
     && mark >= 0
     && nNew >= 0
     && this->size() == preSize + nNew
     && occupancyCellOffsets_.size() == nCells + 1
     && occupancyOrderedParcels_.size() == preSize;

    const char* outcome = "fallback-guard";
    bool current = false;

    if (!guard)
    {
        buildCellOccupancy();
        dsmcRxn_.nFbGuard += 1.0;
    }
    else if (nNew == 0)
    {
        // Only typeId-type reactions (or allowSplitting no): the flat view
        // was not invalidated and is still exactly the cloud.
        if (occupancyOrderedParcelsValid_)
        {
            dsmcRxn_.nSkipZero += 1.0;
            current = true;
            outcome = "skip-zero";
        }
        else
        {
            buildCellOccupancy();
            dsmcRxn_.nFbGuard += 1.0;
        }
    }
    else
    {
        UList<dsmcParcel*> newParcels(moveAppendedParcels_.data() + mark, nNew);
        const double inputValidate0 = dsmcOccValidateWall_;

        if (appendParcelsToOccupancy(newParcels, true))
        {
            // Full-rebuild exit state of the OMP path: valid=true,
            // materialized=false (addParticle already cleared the latter).
            // dsmcOverlapOccCurrent_ is deliberately NOT touched and
            // moveAppendedParcels_ is NOT cleared (Cloud::move reuses it).
            occupancyOrderedParcelsValid_ = true;

            const double wall =
                std::chrono::duration<double>(clk::now() - tEntry).count()
              - (dsmcOccValidateWall_ - inputValidate0);
            dsmcRxn_.nAppend += 1.0;
            dsmcRxn_.nParcels += double(nNew);
            if (dsmcOccLast_.nAffected == 1)
            {
                dsmcRxn_.nShapeSingle += 1.0;
            }
            else
            {
                dsmcRxn_.nShapeMulti += 1.0;
            }
            if (dsmcOccLast_.maxPerCell >= 2)
            {
                dsmcRxn_.nShapeSameCell += 1.0;
            }
            if (dsmcOccLast_.touchesLast)
            {
                dsmcRxn_.nShapeLast += 1.0;
            }
            dsmcRxn_.nKeptOneCells += double(dsmcOccLast_.nKeptOne);
            dsmcRxn_.nKeptZeroCells += double(dsmcOccLast_.nKeptZero);
            dsmcRxn_.tAppend += wall;
            dsmcRxn_.tSort += dsmcOccLast_.tSort;
            dsmcRxn_.tShift += dsmcOccLast_.tShift;
            dsmcRxn_.tOffsets += dsmcOccLast_.tOffsets;
            dsmcRxn_.tDerived += dsmcOccLast_.tDerived;

            // Counted in build_occupancy so the total stays comparable.
            if (profileSummary_ && profileTimingActive_)
            {
                profileBuildCellOccupancyWall_ += wall;
                profileBuildCellOccupancyCpu_ +=
                    mesh_.time().elapsedCpuTime() - cpu0;
            }
            current = true;
            outcome = "append";
        }
        else
        {
            buildCellOccupancy();
            dsmcRxn_.nFbThreshold += 1.0;
            outcome = "fallback-threshold";
        }
    }

    if (dsmcRxnAppendValidate_)
    {
        validateReactionOccupancy(outcome, appSize0);
    }

    return current;
}


double Foam::dsmcCloud::occupancyValidationWall() const
{
    return dsmcOccValidateWall_;
}


void Foam::dsmcCloud::occupancyAppendCoverage
(
    double& nEmpty,
    double& nSingleCell,
    double& nMultiCell,
    double& nLastCell
) const
{
    nEmpty = dsmcOccAppendEmpty_;
    nSingleCell = dsmcOccAppendSingleCell_;
    nMultiCell = dsmcOccAppendMultiCell_;
    nLastCell = dsmcOccAppendLastCell_;
}


bool Foam::dsmcCloud::occupancySelfTest
(
    const word& mode,
    label& nScratch,
    label& nErrors
)
{
    // H4: exercise the real incremental append with the input shapes that a
    // production schedule rarely produces (empty receive, all parcels in one
    // cell) and with the last cell index.  The scratch parcels are not in
    // the cloud; the occupancy is rebuilt from the cloud afterwards, so the
    // test leaves no trace beyond the counters.
    const label nCells = mesh_.nCells();
    if (nCells <= 0)
    {
        nScratch = 0;
        nErrors = 1;
        return false;
    }

    label celli = 0;
    if (mode == "lastCell" || mode == "singleCell")
    {
        if (replicatedMeshActive())
        {
            celli = -1;
            for (label c = nCells - 1; c >= 0; --c)
            {
                if (replicatedMesh_->isMyCell(c))
                {
                    celli = c;
                    break;
                }
            }
            if (celli < 0)
            {
                nScratch = 0;
                nErrors = 1;
                return false;
            }
        }
        else
        {
            celli = (mode == "lastCell") ? nCells - 1 : 0;
        }
    }

    List<dsmcParcel*> scratch;
    if (mode == "empty")
    {
        nScratch = 0;
    }
    else if (mode == "singleCell" || mode == "lastCell")
    {
        nScratch = 3;
        scratch.setSize(nScratch);
        for (label i = 0; i < nScratch; ++i)
        {
            auto* p = new dsmcParcel
            (
                mesh_,
                mesh_.cellCentres()[celli],
                celli,
                false
            );
            p->cell() = celli;
            scratch[i] = p;
        }
    }
    else
    {
        nScratch = 0;
        nErrors = 1;
        return false;
    }

    const labelList oldOff = occupancyCellOffsets_;
    const label oldSize = occupancyOrderedParcels_.size();
    List<dsmcParcel*> oldOrdered;
    if (oldSize > 0)
    {
        oldOrdered.setSize(oldSize);
        forAll(oldOrdered, i)
        {
            oldOrdered[i] = occupancyOrderedParcels_[i];
        }
    }

    dsmcOccSelfTestActive_ = true;
    const bool appended = appendReceivedToOccupancy(scratch);
    dsmcOccSelfTestActive_ = false;

    nErrors = 0;

    // 1. offsets: size, monotonic, end == total
    if (occupancyCellOffsets_.size() != nCells + 1)
    {
        ++nErrors;
    }
    else
    {
        for (label c = 0; c < nCells; ++c)
        {
            if (occupancyCellOffsets_[c] > occupancyCellOffsets_[c+1])
            {
                ++nErrors;
                break;
            }
        }
        if
        (
            occupancyCellOffsets_[nCells]
         != occupancyOrderedParcels_.size()
        )
        {
            ++nErrors;
        }
    }

    // 2. each scratch parcel appears exactly once
    forAll(scratch, i)
    {
        label nFound = 0;
        forAll(occupancyOrderedParcels_, k)
        {
            if (occupancyOrderedParcels_[k] == scratch[i])
            {
                ++nFound;
            }
        }
        if (nFound != 1)
        {
            ++nErrors;
        }
    }

    // 3. scratch parcels sit in the requested cell's segment
    if (nScratch > 0 && occupancyCellOffsets_.size() == nCells + 1)
    {
        forAll(scratch, i)
        {
            bool inCell = false;
            for
            (
                label k = occupancyCellOffsets_[celli];
                k < occupancyCellOffsets_[celli + 1];
                ++k
            )
            {
                if (occupancyOrderedParcels_[k] == scratch[i])
                {
                    inCell = true;
                    break;
                }
            }
            if (!inCell)
            {
                ++nErrors;
            }
        }
    }

    // 4. cells other than the target keep their old sequence, shifted by
    //    the inserted count
    if
    (
        occupancyCellOffsets_.size() == nCells + 1
     && oldOff.size() == nCells + 1
    )
    {
        for (label c = 0; c < nCells; ++c)
        {
            if (c == celli) continue;
            const label shift = (c > celli) ? nScratch : 0;
            const label oldLen = oldOff[c+1] - oldOff[c];
            const label newLen =
                occupancyCellOffsets_[c+1] - occupancyCellOffsets_[c];
            if (newLen != oldLen)
            {
                ++nErrors;
                break;
            }
            bool same = true;
            for (label k = 0; k < oldLen; ++k)
            {
                if
                (
                    occupancyOrderedParcels_[occupancyCellOffsets_[c] + k]
                 != oldOrdered[oldOff[c] + k]
                )
                {
                    same = false;
                    break;
                }
            }
            if (!same)
            {
                ++nErrors;
                break;
            }
            (void)shift;
        }
    }

    // Restore the occupancy from the cloud (also drops the scratch pointers).
    buildCellOccupancy();

    forAll(scratch, i)
    {
        delete scratch[i];
    }

    return appended;
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


void Foam::dsmcCloud::invalidateParcelTraversalCaches()
{
    // Keep the current arrays alive until the next explicit rebuild. Reaction
    // and split models can add parcels inside the OpenMP collision loop while
    // other threads are still reading the current-step occupancy snapshot.
    moveOrderedParcelsValid_ = false;
    moveAppendedParcels_.clear();
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;
}


void Foam::dsmcCloud::addParticle(dsmcParcel* pPtr)
{
    Cloud<dsmcParcel>::addParticle(pPtr);

    if (moveOrderedParcelsValid_)
    {
        moveAppendedParcels_.append(pPtr);
        occupancyOrderedParcelsValid_ = false;
        cellOccupancyMaterialized_ = false;
    }
    else if (!moveAppendCaptureActive_)
    {
        invalidateParcelTraversalCaches();
    }
    else
    {
        occupancyOrderedParcelsValid_ = false;
        cellOccupancyMaterialized_ = false;
    }
}


void Foam::dsmcCloud::deleteParticle(dsmcParcel& p)
{
    Cloud<dsmcParcel>::deleteParticle(p);
    invalidateParcelTraversalCaches();
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

    #ifdef _OPENMP
    if (openmpEnabled_ && openmpMoveEnabled_ && omp_in_parallel())
    {
        #pragma omp critical(dsmcAddParticle)
        {
            addParticle(pPtr);
        }

        return;
    }
    #endif

    addParticle(pPtr);
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

    #ifdef _OPENMP
    if (openmpEnabled_ && openmpMoveEnabled_ && omp_in_parallel())
    {
        #pragma omp critical(dsmcAddParticle)
        {
            addParticle(pPtr);
        }

        return;
    }
    #endif

    addParticle(pPtr);
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
        return collisionSample01();
    }

    scalar energyRatio;

    scalar P;

    do
    {
        P = 0;

        energyRatio = collisionSample01();

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
    } while (P < collisionSample01());

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
        return collisionSample01();
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
        rPSIm = collisionSample01();
        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);
    } while (prob < collisionSample01());

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
    Cloud<dsmcParcel>(mesh, cloudName, false, !dsmcFilteredInitialRead(mesh, cloudName)),
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
    offloadProbe_
    (
        controlDict_.lookupOrDefault<bool>("offloadProbe", false)
    ),
    offloadObserve_
    (
        controlDict_.lookupOrDefault<bool>("offloadObserve", false)
    ),
    offloadObserveDistMap_
    (
        controlDict_.lookupOrDefault<bool>("offloadObserveDistMap", true)
    ),
    offloadEnable_
    (
        controlDict_.lookupOrDefault<bool>("offloadEnable", false)
    ),
    offloadProbeInited_(false),
    offloadProbeRank_(0),
    offloadProbeEpoch_(0.0),
    tsMoveStart_(-1.0),
    tsMoveEnd_(-1.0),
    tsMigEnter_(-1.0),
    tsOccEnd_(-1.0),
    tsCollEnd_(-1.0),
    tsPostEnd_(-1.0),
    offHistP_(),
    offHistT_(),
    offHistB_(),
    offSafeP_(0),
    offSafeT_(0),
    offSafeB_(0),
    offMaxSpeed_(0.0),
    offDistBuilt_(false),
    offDistRef_(1.0),
    offDist_(),
    offSafeDist_(),
    offDistBuildWall_(0.0),
    nCandidatesPerCell_(mesh_.nCells(), 0),
    moveItersPerCell_(mesh_.nCells(), 0),
    moveItersPerCellCumulative_(mesh_.nCells(), 0),
    moveOrderedParcels_(),
    moveOrderedThreadOffsets_(),
    moveOrderedParcelsValid_(false),
    moveAppendCaptureActive_(false),
    moveAppendCaptureDepth_(0),
    moveAppendedParcels_(),
    profileFullEvolveWall_(0.0),
    profileMoveAndCollideWall_(0.0),
    profileMoveWall_(0.0),
    profileBuildCellOccupancyWall_(0.0),
    profileCollisionWall_(0.0),
    profilePostFieldsWall_(0.0),
    profilePostReactionsWall_(0.0),
    profilePostFieldCalcWall_(0.0),
    profilePostFieldWriteWall_(0.0),
    profilePostControllersWall_(0.0),
    profilePostBoundariesWall_(0.0),
    profilePostBoundaryMeasWall_(0.0),
    profilePostCleanWall_(0.0),
    profilePostTrackingCleanWall_(0.0),
    profilePostBoundaryMeasCleanWall_(0.0),
    profilePostCellMeasCleanWall_(0.0),
    profileFullEvolveCpu_(0.0),
    profileMoveAndCollideCpu_(0.0),
    profileMoveCpu_(0.0),
    profileBuildCellOccupancyCpu_(0.0),
    profileCollisionCpu_(0.0),
    profilePostFieldsCpu_(0.0),
    profilePostReactionsCpu_(0.0),
    profilePostFieldCalcCpu_(0.0),
    profilePostFieldWriteCpu_(0.0),
    profilePostControllersCpu_(0.0),
    profilePostBoundariesCpu_(0.0),
    profilePostBoundaryMeasCpu_(0.0),
    profilePostCleanCpu_(0.0),
    profilePostTrackingCleanCpu_(0.0),
    profilePostBoundaryMeasCleanCpu_(0.0),
    profilePostCellMeasCleanCpu_(0.0),
    occupancyOrderedParcels_(),
    occupancyCellOffsets_(),
    occupancyActiveCells_(),
    occupancyCollisionCells_(),
    occupancyOwnedCollisionCells_(),
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
            dsmcSigmaReadOption(mesh_),
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
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
    )
    // §47: "partition" removed from the valid set — it was accepted but
    // never implemented (noTimeCounter has no partition branch), so it
    // silently behaved as dynamic; validation now matches reality.
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

    buildConstProps();

    bool inSolverFilled = false;

    if (readFields)
    {
        if (dsmcFilteredInitialRead(mesh, cloudName))
        {
            // Tier-2: construct the replicated mesh early so cellOwner_ is
            // available, then read only owned parcels (positions file and
            // parcel fields) - see doc/worklog/v2506/dsmcFoam_todolist.md.
            replicatedMesh_.reset(new dsmcReplicatedMesh(*this, mesh_));
            replicatedMesh_->initialize();

            Cloud<dsmcParcel>::readCloudFiltered
            (
                replicatedMesh_->cellOwner(),
                replicatedMesh_->myRank()
            );

            dsmcParcel::readFieldsFiltered
            (
                *this,
                initialReadKeep(),
                initialReadNFull()
            );

            Info<< "Replicated mesh: rank " << replicatedMesh_->myRank()
                << " filtered initial read: kept " << this->size()
                << " / " << initialReadNFull() << " particles (skipped "
                << initialReadNFull() - this->size()
                << " non-owned during read)" << endl;
        }
        else if (dsmcInSolverFill(mesh, cloudName))
        {
            // Optimisation 3b: in-solver parallel fill from
            // dsmcInitialiseDict (zero file round-trip)
            replicatedMesh_.reset(new dsmcReplicatedMesh(*this, mesh_));
            replicatedMesh_->initialize();

            IOdictionary dsmcInitialiseDict
            (
                IOobject
                (
                    "dsmcInitialiseDict",
                    mesh_.time().system(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            initialFillFromDict(dsmcInitialiseDict);

            inSolverFilled = true;

            Info<< "Replicated mesh: rank " << replicatedMesh_->myRank()
                << " in-solver fill: " << this->size() << " particles"
                << endl;
        }
        else
        {
            dsmcParcel::readFields(*this);
        }
    }

    if (!inSolverFilled)
    {
    // v1706 value-constructor + readIfPresent() does NOT read the file when
    // the read option is MUST_READ (warning only).  Read dsmcSigmaTcRMax
    // explicitly when present so fresh starts get the initial max cross
    // section (processor checkpoint restart restores it separately in
    // dsmcReplicatedMesh::restoreProcessorCheckpoint()).
    // dsmcSigmaTcRMax read guard.  The value constructor with
    // READ_IF_PRESENT reads the field when the file exists; the processor
    // checkpoint restart restores it explicitly.  If neither is available
    // the field would stay zero, which silently disables all NTC collisions
    // (zero is an absorbing state) - fail loudly instead.
    {
        const fileName sigmaFile
        (
            mesh_.time().path()
           /mesh_.time().timeName()
           /(this->name() + "SigmaTcRMax")
        );

        if (isFile(sigmaFile))
        {
            Info<< "dsmcSigmaTcRMax: initial field found at " << sigmaFile
                << endl;
        }
        else
        {
            const fileName checkpoint
            (
                mesh_.time().path()
               /"processor0"
               /mesh_.time().timeName()
               /"dsmcSigmaTcRMax"
            );

            if (!isFile(checkpoint))
            {
                FatalErrorInFunction
                    << "Cannot find dsmcSigmaTcRMax at " << sigmaFile << nl
                    << "A zero initial field would silently disable all "
                    << "NTC collisions." << exit(FatalError);
            }

            Info<< "dsmcSigmaTcRMax: expecting processor checkpoint "
                << "restore from " << checkpoint << endl;
        }
    }
    }

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
        if (replicatedMesh_.empty())
        {
            replicatedMesh_.reset(new dsmcReplicatedMesh(*this, mesh_));
            replicatedMesh_->initialize();
        }
        // else: already constructed and initialized before the filtered
        // initial read (Tier-2 path)
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
    offloadProbe_(false),
    offloadObserve_(false),
    offloadObserveDistMap_(true),
    offloadEnable_(false),
    offloadProbeInited_(false),
    offloadProbeRank_(0),
    offloadProbeEpoch_(0.0),
    tsMoveStart_(-1.0),
    tsMoveEnd_(-1.0),
    tsMigEnter_(-1.0),
    tsOccEnd_(-1.0),
    tsCollEnd_(-1.0),
    tsPostEnd_(-1.0),
    offHistP_(),
    offHistT_(),
    offHistB_(),
    offSafeP_(0),
    offSafeT_(0),
    offSafeB_(0),
    offMaxSpeed_(0.0),
    offDistBuilt_(false),
    offDistRef_(1.0),
    offDist_(),
    offSafeDist_(),
    offDistBuildWall_(0.0),
    nCandidatesPerCell_(),
    moveItersPerCell_(),
    moveItersPerCellCumulative_(),
    moveOrderedParcels_(),
    moveOrderedThreadOffsets_(),
    moveOrderedParcelsValid_(false),
    moveAppendCaptureActive_(false),
    moveAppendCaptureDepth_(0),
    moveAppendedParcels_(),
    profileFullEvolveWall_(0.0),
    profileMoveAndCollideWall_(0.0),
    profileMoveWall_(0.0),
    profileBuildCellOccupancyWall_(0.0),
    profileCollisionWall_(0.0),
    profilePostFieldsWall_(0.0),
    profilePostReactionsWall_(0.0),
    profilePostFieldCalcWall_(0.0),
    profilePostFieldWriteWall_(0.0),
    profilePostControllersWall_(0.0),
    profilePostBoundariesWall_(0.0),
    profilePostBoundaryMeasWall_(0.0),
    profilePostCleanWall_(0.0),
    profilePostTrackingCleanWall_(0.0),
    profilePostBoundaryMeasCleanWall_(0.0),
    profilePostCellMeasCleanWall_(0.0),
    profileFullEvolveCpu_(0.0),
    profileMoveAndCollideCpu_(0.0),
    profileMoveCpu_(0.0),
    profileBuildCellOccupancyCpu_(0.0),
    profileCollisionCpu_(0.0),
    profilePostFieldsCpu_(0.0),
    profilePostReactionsCpu_(0.0),
    profilePostFieldCalcCpu_(0.0),
    profilePostFieldWriteCpu_(0.0),
    profilePostControllersCpu_(0.0),
    profilePostBoundariesCpu_(0.0),
    profilePostBoundaryMeasCpu_(0.0),
    profilePostCleanCpu_(0.0),
    profilePostTrackingCleanCpu_(0.0),
    profilePostBoundaryMeasCleanCpu_(0.0),
    profilePostCellMeasCleanCpu_(0.0),
    occupancyOrderedParcels_(),
    occupancyCellOffsets_(),
    occupancyActiveCells_(),
    occupancyCollisionCells_(),
    occupancyOwnedCollisionCells_(),
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

        printInitialiseTimeStepCalculation(t, dsmcInitialiseDict);
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

// * * * * * * * * * * dlb_offload.md P0/P1 probe  * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- One-time conservative boundary-clearance map (dlb_offload.md §3).
//  dist[c]: Dijkstra over cell adjacency with centre-to-centre edge
//  weights, sources = cells owning any boundary face at 0 (v1 treats
//  every patch as a side-effect boundary).  Lower-bound chain is valid
//  because distance-to-a-set is 1-Lipschitz.  safeDist[c] = dist[c] -
//  circumradius[c] (exact circumradius over the cell's vertices) then
//  bounds the clearance of ANY point of cell c from the boundary.
void offloadComputeClearance
(
    const dynamicFvMesh& mesh,
    scalarField& dist,
    scalarField& safeDist
)
{
    const label nCells = mesh.nCells();
    const label nIntFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const vectorField& cc = mesh.cellCentres();
    const vectorField& pts = mesh.points();
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();

    safeDist.setSize(nCells);

    // exact circumradius per cell (a convex function attains its max over
    // a polytope at a vertex; every cell vertex lies on one of its faces)
    for (label c = 0; c < nCells; ++c)
    {
        scalar s = 0.0;
        const labelList& fcs = cells[c];
        forAll(fcs, i)
        {
            const face& f = faces[fcs[i]];
            forAll(f, j)
            {
                s = max(s, mag(pts[f[j]] - cc[c]));
            }
        }
        safeDist[c] = s;
    }

    dist.setSize(nCells);
    dist = GREAT;
    typedef std::pair<scalar, label> clearItem;
    std::priority_queue
    <
        clearItem,
        std::vector<clearItem>,
        std::greater<clearItem>
    > pq;
    for (label f = nIntFaces; f < faceOwner.size(); ++f)
    {
        dist[faceOwner[f]] = 0.0;
    }
    for (label c = 0; c < nCells; ++c)
    {
        if (dist[c] == 0.0)
        {
            pq.push(clearItem(0.0, c));
        }
    }
    while (!pq.empty())
    {
        const clearItem top = pq.top();
        pq.pop();
        const label c = top.second;
        if (top.first > dist[c])
        {
            continue;
        }
        const labelList& fcs = cells[c];
        forAll(fcs, i)
        {
            const label f = fcs[i];
            if (f >= nIntFaces)
            {
                continue;
            }
            const label nb = (faceOwner[f] == c) ? faceNeighbour[f] : faceOwner[f];
            const scalar nd = dist[c] + mag(cc[nb] - cc[c]);
            if (nd < dist[nb])
            {
                dist[nb] = nd;
                pq.push(clearItem(nd, nb));
            }
        }
    }

    for (label c = 0; c < nCells; ++c)
    {
        safeDist[c] = dist[c] - safeDist[c];
    }
}

} // anonymous namespace
} // namespace Foam


double Foam::dsmcCloud::offloadProbeNow() const
{
    return std::chrono::duration<double>
    (
        std::chrono::steady_clock::now().time_since_epoch()
    ).count() - offloadProbeEpoch_;
}


void Foam::dsmcCloud::offloadProbeInit()
{
    if (offloadProbeInited_)
    {
        return;
    }
    offloadProbeInited_ = true;

    // MPI may not be initialised here: pure-OMP runs never start it and
    // the replicated mesh (which does) may not exist yet on some paths.
    int mpiIsInitialized = 0;
    MPI_Initialized(&mpiIsInitialized);
    if (mpiIsInitialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &offloadProbeRank_);

        // One shared clock reference: capture after a barrier so same-node
        // ranks agree within barrier-exit skew (CLOCK_MONOTONIC is per
        // boot, so cross-node absolute comparison is not claimed).
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        offloadProbeRank_ = 0;
    }
    offloadProbeEpoch_ = std::chrono::duration<double>
    (
        std::chrono::steady_clock::now().time_since_epoch()
    ).count();

    for (label b = 0; b < 12; ++b)
    {
        offHistP_[b] = 0;
        offHistT_[b] = 0;
        offHistB_[b] = 0;
    }
    offSafeP_ = 0;
    offSafeT_ = 0;
    offSafeB_ = 0;
    offMaxSpeed_ = 0.0;

    // Per-rank probe CSV in the run directory (CWD = case dir under the
    // job's mpirun).  Per-rank Info does not reach the job stdout in
    // replicated-mesh mode — only the cloud master's does — so all ranks
    // write their own file.
    offloadProbeFile_.open
    (
        "offloadProbe-r" + std::to_string(offloadProbeRank_) + ".csv"
    );
    offloadProbeFile_
        << "step,mv0,mv1,mig0,migPrep1,migSx1,mig1,migSent,migRecv,"
           "migSentB,migRecvB,occ1,coll1,post1";
    if (offloadObserve_)
    {
        offloadProbeFile_ << ",safeP,safeT,safeB,maxU,bytes";
        for (label b = 0; b < 12; ++b)
        {
            offloadProbeFile_ << ",bP" << b;
        }
        for (label b = 0; b < 12; ++b)
        {
            offloadProbeFile_ << ",bT" << b;
        }
        for (label b = 0; b < 12; ++b)
        {
            offloadProbeFile_ << ",bB" << b;
        }
    }
    offloadProbeFile_ << "\n";

    if ((offloadObserve_ || offloadEnable_) && offloadObserveDistMap_)
    {
        const auto t0 = steadyWallClock::now();
        offloadComputeClearance(mesh_, offDist_, offSafeDist_);
        offDistBuildWall_ =
            std::chrono::duration<scalar>(steadyWallClock::now() - t0).count();
        offDistBuilt_ = true;

        scalar maxDist = 0.0;
        forAll(offDist_, c)
        {
            maxDist = max(maxDist, offDist_[c]);
        }
        offDistRef_ = max(maxDist/12.0, ROOTVSMALL);

        // dlb_offload.md P1: connected-component sizes of the
        // positive-clearance set (the candidate task region structure).
        label nComp = 0;
        labelList compId(offDist_.size(), -1);
        DynamicList<label> stack(offDist_.size()/(1000*24) + 1024);
        DynamicList<label> compSizes;
        forAll(offSafeDist_, c)
        {
            if (offSafeDist_[c] <= 0 || compId[c] >= 0)
            {
                continue;
            }
            const label id = nComp++;
            compId[c] = id;
            stack.clear();
            stack.append(c);
            label size = 0;
            while (stack.size())
            {
                const label cc = stack.remove();
                ++size;
                const labelList& fcs = mesh_.cells()[cc];
                forAll(fcs, i)
                {
                    const label f = fcs[i];
                    if (f >= mesh_.nInternalFaces())
                    {
                        continue;
                    }
                    const label nb =
                        (mesh_.faceOwner()[f] == cc)
                      ? mesh_.faceNeighbour()[f]
                      : mesh_.faceOwner()[f];
                    if (offSafeDist_[nb] > 0 && compId[nb] < 0)
                    {
                        compId[nb] = id;
                        stack.append(nb);
                    }
                }
            }
            compSizes.append(size);
        }
        std::sort(compSizes.begin(), compSizes.end(),
                  std::greater<label>());

        if (offloadProbeRank_ == 0)
        {
            OStringStream os;
            os << "clearance>0 components: n=" << nComp
               << " top=";
            for (label i = 0; i < min(label(8), nComp); ++i)
            {
                os << compSizes[i] << (i < min(label(7), nComp - 1)
                                       ? "," : "");
            }
            Info<< "Offload observe: clearance map built in "
                << offDistBuildWall_ << " s"
                << " (nCells " << offDist_.size()
                << ", maxDist " << maxDist
                << ", bucket width " << offDistRef_
                << ", all patches treated as side-effect boundaries)"
                << nl << "Offload observe: " << os.str() << endl;
        }
    }

    if (offloadProbeRank_ == 0)
    {
        Info<< "Offload probe: enabled (observe=" << offloadObserve_
            << ", epoch aligned via MPI_Barrier)" << endl;
    }
}


void Foam::dsmcCloud::offloadProbeEmit()
{
    if (!offloadProbeFile_.is_open())
    {
        return;
    }
    std::ofstream& os = offloadProbeFile_;
    os << mesh_.time().timeIndex()
       << ',' << tsMoveStart_ << ',' << tsMoveEnd_
       << ',' << tsMigEnter_;
    if
    (
        replicatedMeshActive()
     && !replicatedMesh_->asyncOverlapEnabled()
    )
    {
        os << ',' << replicatedMesh_->offPrep1_
           << ',' << replicatedMesh_->offSx1_
           << ',' << replicatedMesh_->offMig1_
           << ',' << replicatedMesh_->offSent_
           << ',' << replicatedMesh_->offRecv_
           << ',' << replicatedMesh_->offSentB_
           << ',' << replicatedMesh_->offRecvB_;
    }
    else
    {
        // clean M1: the migration stamps are only produced by the
        // synchronous flat path.  Async overlap runs must report NA (-1)
        // instead of last step's values.
        os << ",-1,-1,-1,-1,-1,-1,-1";
    }
    os << ',' << tsOccEnd_
       << ',' << tsCollEnd_
       << ',' << tsPostEnd_;
    if (offloadObserve_)
    {
        // candidate task bytes: safe particles round-trip (task + result)
        // clean M14: widen before the multiply (int32 overflows above
        // ~9e6 safe parcels); the probe column is a diagnostic estimate.
        const double taskBytes =
            2.0*double(offSafeP_)*double(sizeof(dsmcParcel::TransferData));
        os << ',' << offSafeP_
           << ',' << offSafeT_
           << ',' << offSafeB_
           << ',' << offMaxSpeed_
           << ',' << taskBytes;
        for (label b = 0; b < 12; ++b)
        {
            os << ',' << offHistP_[b];
        }
        for (label b = 0; b < 12; ++b)
        {
            os << ',' << offHistT_[b];
        }
        for (label b = 0; b < 12; ++b)
        {
            os << ',' << offHistB_[b];
        }
    }
    os << '\n';
    if ((mesh_.time().timeIndex() % 50) == 0)
    {
        os.flush();
    }
}


void Foam::dsmcCloud::evolve()
{
    if
    (
        (offloadProbe_ || offloadObserve_ || offloadEnable_)
     && !offloadProbeInited_
    )
    {
        offloadProbeInit();
    }

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

    if (offloadProbe_)
    {
        tsPostEnd_ = offloadProbeNow();
        offloadProbeEmit();
    }
}


void Foam::dsmcCloud::evolve_moveAndCollide()
{
    const steadyWallClock::time_point moveAndCollideWallStart =
        steadyWallClock::now();
    const scalar moveAndCollideCpuStart =
        profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;

    // §28: evolve-residual attribution.  The four core phases (move/
    // build/coll/post) leave ~20% of fullEvolve unattributed; two known
    // residents are the per-step model-info updates (timed here) and
    // updateParticleCounts (already accumulated in the mesh, reported
    // next to the residual).  File-static: zero header change.
    const auto dsmcEvolveInfoStart = steadyWallClock::now();
    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();
    dsmcEvolveUpdateInfoWall_ +=
        elapsedWallSeconds(dsmcEvolveInfoStart);

    dsmcParcel::trackingData td(*this);
    td.moveDetailProfile =
        profileDetail_
     && controlDict_.lookupOrDefault<bool>("moveDetailProfile", false);
    if (td.moveDetailProfile)
    {
        buildNearWallCells(mesh_);
        td.nearWallCells = &dsmcNearWallCells;
    }

    if (debug)
    {
        this->dumpParticlePositions();
    }

    const bool debugReplicatedParticleState =
        replicatedMeshActive()
     && controlDict_.lookupOrDefault<bool>
        (
            "replicatedMeshDebugParticleState",
            false
        );

    auto reportReplicatedParticleState = [&](const char* tag)
    {
        if (!debugReplicatedParticleState)
        {
            return;
        }

        label nParcels = 0;
        label nFree = 0;
        label nStuck = 0;
        label badCell = 0;
        label nonOwner = 0;
        label badFace = 0;
        label badTetFace = 0;
        label badTetPt = 0;
        label badStepFraction = 0;
        label minCell = mesh_.nCells();
        label maxCell = -1;
        label firstBadCell = -1;
        label firstBadFace = -1;
        label firstBadTetFace = -1;
        label firstBadTetPt = -1;

        forAllConstIter(dsmcCloud, *this, iter)
        {
            const dsmcParcel& p = iter();
            ++nParcels;

            if (p.isFree())
            {
                ++nFree;
            }
            else
            {
                ++nStuck;
            }

            const label cellI = p.cell();
            if (cellI >= 0 && cellI < mesh_.nCells())
            {
                minCell = min(minCell, cellI);
                maxCell = max(maxCell, cellI);

                if (!replicatedMesh_->isMyCell(cellI))
                {
                    ++nonOwner;
                }
            }
            else
            {
                ++badCell;
                if (firstBadCell < 0)
                {
                    firstBadCell = cellI;
                }
            }

            const label faceI = p.face();
            if (faceI < -1 || faceI >= mesh_.nFaces())
            {
                ++badFace;
                if (firstBadFace < 0)
                {
                    firstBadFace = faceI;
                }
            }

            const label tetFaceI = p.tetFace();
            const label tetPtI = p.tetPt();
            if (tetFaceI < 0 || tetFaceI >= mesh_.nFaces())
            {
                ++badTetFace;
                if (firstBadTetFace < 0)
                {
                    firstBadTetFace = tetFaceI;
                }
            }
            else if (tetPtI < 0 || tetPtI >= mesh_.faces()[tetFaceI].size())
            {
                ++badTetPt;
                if (firstBadTetPt < 0)
                {
                    firstBadTetPt = tetPtI;
                }
            }

            if (p.stepFraction() < -SMALL || p.stepFraction() > 1 + SMALL)
            {
                ++badStepFraction;
            }
        }

        Pout<< "Replicated mesh particle state [" << tag << "]: "
            << "rank=" << replicatedMesh_->myRank()
            << " parcels=" << nParcels
            << " free=" << nFree
            << " stuck=" << nStuck
            << " cellRange=[" << minCell << "," << maxCell << "]"
            << " badCell=" << badCell
            << " nonOwner=" << nonOwner
            << " badFace=" << badFace
            << " badTetFace=" << badTetFace
            << " badTetPt=" << badTetPt
            << " badStepFraction=" << badStepFraction
            << " firstBad(cell face tetFace tetPt)=("
            << firstBadCell << " " << firstBadFace << " "
            << firstBadTetFace << " " << firstBadTetPt << ")"
            << endl;
    };

    if (openmpEnabled_ && openmpMoveEnabled_)
    {
        beginMoveAppendCapture();
    }

    controllers_.controlBeforeMove();
    boundaries_.controlBeforeMove();

    if (replicatedMeshActive() && replicatedMesh_->migrationCalls() == 0)
    {
        Info<< "Replicated mesh: initial particle distribution" << nl << endl;
        const label oldMigrationContext =
            replicatedMesh_->migrationProfileContext();
        replicatedMesh_->setMigrationProfileContext
        (
            dsmcReplicatedMesh::migrationProfileInitial
        );
        if (Pstream::parRun())
        {
            replicatedMesh_->distributeInitialParticles();
        }
        else
        {
            replicatedMesh_->migrateParticlesByCellOwner();
        }
        replicatedMesh_->updateParticleCounts();
        replicatedMesh_->setMigrationProfileContext(oldMigrationContext);
        if (!hasMoveOrderedParcels())
        {
            clearMoveOrderedParcels();
        }
        buildCellOccupancy();
        reportReplicatedParticleState("afterInitialDistribution");
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
        const label oldMigrationContext =
            replicatedMesh_->migrationProfileContext();
        replicatedMesh_->setMigrationProfileContext
        (
            dsmcReplicatedMesh::migrationProfileRegular
        );
        replicatedMesh_->migrateFinish();
        replicatedMesh_->updateParticleCounts();
        replicatedMesh_->setMigrationProfileContext(oldMigrationContext);
        if (!hasMoveOrderedParcels())
        {
            clearMoveOrderedParcels();
        }
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

    if (openmpEnabled_ && openmpMoveEnabled_ && ompNumThreads_ > 1)
    {
        prepareOpenMPMoveMeshData();
    }

    // dlb_offload.md P2-P5: front dispatch of offload tasks before the
    // move (donor packs/sends/deletes; helper posts receives)
    if
    (
        replicatedMeshActive()
     && replicatedMesh_->offloadEnabled()
    )
    {
        replicatedMesh_->offloadDispatch();
    }

    //scalar timer = mesh_.time().elapsedCpuTime();
    const steadyWallClock::time_point moveWallStart = steadyWallClock::now();
    const scalar moveCpuStart =
        profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;
    reportReplicatedParticleState("beforeMove");
    if (offloadProbe_)
    {
        tsMoveStart_ = offloadProbeNow();
    }
    Cloud<dsmcParcel>::move(td, deltaTValue());
    if (openmpEnabled_ && openmpMoveEnabled_)
    {
        endMoveAppendCapture();
    }
    if (offloadProbe_)
    {
        tsMoveEnd_ = offloadProbeNow();
    }
    if (offloadObserve_)
    {
        // merge thread-local observe counters (master td accumulated
        // inside Cloud::move via the per-thread critical sections)
        for (label b = 0; b < 12; ++b)
        {
            offHistP_[b] = td.offHistP[b];
            offHistT_[b] = td.offHistT[b];
            offHistB_[b] = td.offHistB[b];
        }
        offSafeP_ = td.offSafeP;
        offSafeT_ = td.offSafeT;
        offSafeB_ = td.offSafeB;
        offMaxSpeed_ = sqrt(td.offMaxMagSqrU);
    }

    // dlb_offload.md P2-P5: helper executes received tasks and inserts
    // results before the regular migration routes them to final owners
    if
    (
        replicatedMeshActive()
     && replicatedMesh_->offloadEnabled()
    )
    {
        replicatedMesh_->offloadExecuteAndInsert();
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

        dsmcMoveDetailCum[0] += double(td.moveParcels);
        dsmcMoveDetailCum[1] += double(td.moveTrackCalls);
        dsmcMoveDetailCum[2] += double(td.moveSameTetNoFaceHits);
        dsmcMoveDetailCum[3] += double(td.moveInternalTetNoFaceHits);
        dsmcMoveDetailCum[4] += double(td.moveFaceHits);
        dsmcMoveDetailCum[5] += double(td.moveProcessorHits);
        dsmcMoveDetailCum[6] += double(td.moveCyclicHits);
        dsmcMoveDetailCum[7] += double(td.movePatchHits);
        dsmcMoveDetailCum[8] += double(td.moveStuckHits);
        dsmcMoveDetailCum[9] += td.moveTrackWallTime;
        dsmcMoveDetailCum[10] += td.moveTrackerWallTime;
        dsmcMoveDetailCum[11] += td.moveBoundaryWallTime;

        dsmcMoveDetailSplitCum[0] += double(td.moveWallCellTrackCalls);
        dsmcMoveDetailSplitCum[1] += td.moveWallCellTrackWallTime;
        dsmcMoveDetailSplitCum[2] += double(td.moveInteriorTrackCalls);
        dsmcMoveDetailSplitCum[3] += td.moveInteriorTrackWallTime;
    }
    if (profileSummary_)
    {
        profileMoveWall_ += elapsedWallSeconds(moveWallStart);
        profileMoveCpu_ += mesh_.time().elapsedCpuTime() - moveCpuStart;
    }
    //Info<< "move" << tab << mesh_.time().elapsedCpuTime() - timer << " s" << endl;

    if (replicatedMeshActive())
    {
        if (offloadProbe_)
        {
            tsMigEnter_ = offloadProbeNow();
        }
        if
        (
            replicatedMesh_->stepCounter() == 0
         || (
                replicatedMesh_->migrateInterval() > 0
             && replicatedMesh_->stepCounter() % replicatedMesh_->migrateInterval() == 0
            )
        )
        {
            const label oldMigrationContext =
                replicatedMesh_->migrationProfileContext();
            replicatedMesh_->setMigrationProfileContext
            (
                dsmcReplicatedMesh::migrationProfileRegular
            );
            const bool asyncOverlap =
                replicatedMesh_->asyncOverlapEnabled()
             && replicatedMesh_->asyncOverlapCapable();
            bool asyncOverlapRan = false;
            if (asyncOverlap)
            {
                // comm_further_opt.md §2: begin -> gap build -> finish.
                // begin packs + posts Isends and deletes migrated-out;
                // the gap builds occupancy over the kept parcels while
                // transfers are in flight; finish deserializes received
                // and appends them into the occupancy incrementally.
                // begin/finish own the moveOrdered lifecycle on this
                // path - the clear-if below must NOT run here: it
                // invalidates occupancyOrderedParcelsValid_ after
                // finish has established the current occupancy, and
                // collide's lazy materialize would then destroy the
                // per-cell lists (clear-then-skip-rebuild semantics).
                replicatedMesh_->migrateBegin();
                if (replicatedMesh_->asyncMigrationPending())
                {
                    replicatedMesh_->markAsyncOverlapGapBuild();
                    buildCellOccupancy();
                    replicatedMesh_->migrateFinish();
                    asyncOverlapRan = true;
                }
            }
            else if (replicatedMeshDelayedReceive)
            {
                replicatedMesh_->migrateBegin();
                replicatedMesh_->migrateFinish();
            }
            else
            {
                replicatedMesh_->migrateParticlesByCellOwner();
            }
            replicatedMesh_->updateParticleCounts();
            replicatedMesh_->setMigrationProfileContext(oldMigrationContext);
            if (!asyncOverlapRan && !hasMoveOrderedParcels())
            {
                clearMoveOrderedParcels();
            }
        }

        replicatedMesh_->advanceStepCounter();
    }

    // dlb_offload.md P2: post-migration particle-multiset checksum at
    // the configured dump step (off/on bitwise comparison)
    if
    (
        replicatedMeshActive()
     && replicatedMesh_->offloadEnabled()
    )
    {
        replicatedMesh_->offloadDumpChecksum();
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
                    const label oldMigrationContext =
                        replicatedMesh_->migrationProfileContext();
                    replicatedMesh_->setMigrationProfileContext
                    (
                        dsmcReplicatedMesh::migrationProfileRegular
                    );
                    replicatedMesh_->migrateFinish();
                    replicatedMesh_->updateParticleCounts();
                    replicatedMesh_->setMigrationProfileContext
                    (
                        oldMigrationContext
                    );
                    if (!hasMoveOrderedParcels())
                    {
                        clearMoveOrderedParcels();
                    }
                }

                Info<< nl
                    << "Replicated mesh: manual cell owner reassignment at step "
                    << currentStep << nl << endl;
                // Manual rebalance re-migrates parcels after the async
                // overlap may have built the occupancy — invalidate it so
                // the unconditional build below runs (correctness first).
                dsmcOverlapOccCurrent_ = false;
                replicatedMesh_->reassignCellOwner();
                const label oldMigrationContext =
                    replicatedMesh_->migrationProfileContext();
                replicatedMesh_->setMigrationProfileContext
                (
                    dsmcReplicatedMesh::migrationProfileManual
                );
                replicatedMesh_->migrateParticlesByCellOwner();
                replicatedMesh_->updateParticleCounts();
                replicatedMesh_->setMigrationProfileContext(oldMigrationContext);
                if (!hasMoveOrderedParcels())
                {
                    clearMoveOrderedParcels();
                }
                break;
            }
        }
    }

    //- Update cell occupancy
    //timer = mesh_.time().elapsedCpuTime();
    // comm_further_opt.md §2: the async path already built occupancy over
    // kept in the gap and appended the received parcels incrementally in
    // migrateFinish — consume the flag instead of rebuilding.
    // clean M12 contract: exactly one consumer per migrateFinish()/
    // appendReceivedToOccupancy() pair.  The flag is set by
    // appendReceivedToOccupancy() on every exit path and consumed here
    // (reset) or in the equivalent post-migration branch; a migrateFinish()
    // that is not followed by one of these consumers must reset it.
    if (dsmcOverlapOccCurrent_)
    {
        dsmcOverlapOccCurrent_ = false;
    }
    else
    {
        buildCellOccupancy();
    }
    //Info<< "buildCellOccupancy" << tab << mesh_.time().elapsedCpuTime() - timer << " s " << endl;

    //- Add electrons back after the move function
    if (findIndex(typeIdList_, "e-") != -1)
    {
        // TODO VINCENT: there is a clever way than rebuilding entire cell occ.
        addElectrons();
        buildCellOccupancy();
    }

    if (offloadProbe_)
    {
        tsOccEnd_ = offloadProbeNow();
    }

    //- Radial weighting for non-Cartesian flows (e.g., axisymmetric). This is
    // where parcels will receive their new RWF and will possibly be cloned or
    // deleted.
    const auto dsmcCoordStart = steadyWallClock::now();
    coordSystem().evolve();
    dsmcEvolveCoordWall_ += elapsedWallSeconds(dsmcCoordStart);

    const auto dsmcCtrlStart = steadyWallClock::now();
    controllers_.controlBeforeCollisions();
    boundaries_.controlBeforeCollisions();
    dsmcEvolveCtrlWall_ += elapsedWallSeconds(dsmcCtrlStart);

    //- Calculate new velocities via stochastic collisions
    //timer = mesh_.time().elapsedCpuTime();
    const steadyWallClock::time_point collisionWallStart =
        steadyWallClock::now();
    const scalar collisionCpuStart =
        profileSummary_ ? mesh_.time().elapsedCpuTime() : 0.0;
    // reactionOccupancyAppend: pre-collision snapshot (audit_1.3.md §5.2.3).
    // Collisions only ever ADD parcels (addNewParcel); addParticle clears
    // occupancyOrderedParcelsValid_ during the loop, so the state the
    // incremental path needs must be captured here.
    // (read the switch on every rank before the first reaction so the
    // profile print does not depend on whether rank 0 ever reacted)
    reactionOccupancyAppendEnabled();
    const label occPreSize = this->size();
    const label occPreMark = moveAppendedParcels_.size();
    const bool occPreValid =
        occupancyOrderedParcelsValid_
     && occupancyCellOffsets_.size() == mesh_.nCells() + 1
     && occupancyCellOffsets_.last() == occPreSize;
    collisions();
    if (profileSummary_)
    {
        profileCollisionWall_ += elapsedWallSeconds(collisionWallStart);
        profileCollisionCpu_ +=
            mesh_.time().elapsedCpuTime() - collisionCpuStart;
    }
    if (offloadProbe_)
    {
        tsCollEnd_ = offloadProbeNow();
    }
    //Info<< "collisions" << tab << mesh_.time().elapsedCpuTime() - timer << " s" << endl;

    //- Reactions may have changed cell occupancy, update only if this step reacted
    bool occAppendCurrent = false;
    if (reactions_.nReactionsPerTimeStep() != 0)
    {
        // Switch off (default): exactly the former buildCellOccupancy().
        occAppendCurrent =
            updateOccupancyAfterReactions(occPreSize, occPreMark, occPreValid);
    }

    controllers_.controlAfterCollisions();
    boundaries_.controlAfterCollisions();

    // The incremental update left a current flat view; if a controller
    // changed the cloud afterwards (WIP density controllers add/delete
    // parcels), fall back to the full rebuild and count it.
    if (occAppendCurrent && !occupancyOrderedParcelsValid_)
    {
        noteReactionOccupancyLateInvalidate();
        buildCellOccupancy();
    }

    if (replicatedMeshActive())
    {
        if (replicatedMesh_->asyncMigrationPending())
        {
            const label oldMigrationContext =
                replicatedMesh_->migrationProfileContext();
            replicatedMesh_->setMigrationProfileContext
            (
                dsmcReplicatedMesh::migrationProfileRegular
            );
            replicatedMesh_->migrateFinish();
            replicatedMesh_->updateParticleCounts();
            replicatedMesh_->setMigrationProfileContext(oldMigrationContext);
            if (!hasMoveOrderedParcels())
            {
                clearMoveOrderedParcels();
            }
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
            if (!hasMoveOrderedParcels())
            {
                clearMoveOrderedParcels();
            }
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
    const bool detailTiming = profileSummary_ && profileDetail_;
    steadyWallClock::time_point detailWallStart = wallStart;
    scalar detailCpu0 = cpu0;

    auto resetDetailTimer = [&]()
    {
        if (detailTiming)
        {
            detailWallStart = steadyWallClock::now();
            detailCpu0 = mesh_.time().elapsedCpuTime();
        }
    };

    auto addDetailTimer = [&](scalar& wallTime, scalar& cpuTime)
    {
        if (detailTiming)
        {
            wallTime += elapsedWallSeconds(detailWallStart);
            cpuTime += mesh_.time().elapsedCpuTime() - detailCpu0;
            resetDetailTimer();
        }
    };

    reactions_.outputData();
    addDetailTimer(profilePostReactionsWall_, profilePostReactionsCpu_);

    if
    (
        !noScheduledFieldOutput()
     || collisionModelUsesMacroscopicFieldTemperature()
    )
    {
        fields_.calculateFields();
        addDetailTimer(profilePostFieldCalcWall_, profilePostFieldCalcCpu_);

        fields_.writeFields();
        addDetailTimer(profilePostFieldWriteWall_, profilePostFieldWriteCpu_);
    }

    controllers_.calculateProps();
    controllers_.outputResults();
    addDetailTimer(profilePostControllersWall_, profilePostControllersCpu_);

    boundaries_.calculateProps();
    boundaries_.outputResults();
    addDetailTimer(profilePostBoundariesWall_, profilePostBoundariesCpu_);

    boundaryMeas_.outputResults();
    addDetailTimer(profilePostBoundaryMeasWall_, profilePostBoundaryMeasCpu_);

    const steadyWallClock::time_point cleanWallStart = detailWallStart;
    const scalar cleanCpu0 = detailCpu0;
    steadyWallClock::time_point cleanStageWallStart = cleanWallStart;
    scalar cleanStageCpu0 = cleanCpu0;

    trackingInfo_.clean();
    if (detailTiming)
    {
        profilePostTrackingCleanWall_ +=
            elapsedWallSeconds(cleanStageWallStart);
        profilePostTrackingCleanCpu_ +=
            mesh_.time().elapsedCpuTime() - cleanStageCpu0;
        cleanStageWallStart = steadyWallClock::now();
        cleanStageCpu0 = mesh_.time().elapsedCpuTime();
    }

    boundaryMeas_.clean();
    if (detailTiming)
    {
        profilePostBoundaryMeasCleanWall_ +=
            elapsedWallSeconds(cleanStageWallStart);
        profilePostBoundaryMeasCleanCpu_ +=
            mesh_.time().elapsedCpuTime() - cleanStageCpu0;
        cleanStageWallStart = steadyWallClock::now();
        cleanStageCpu0 = mesh_.time().elapsedCpuTime();
    }

    cellMeas_.clean();

    if (detailTiming)
    {
        profilePostCellMeasCleanWall_ +=
            elapsedWallSeconds(cleanStageWallStart);
        profilePostCellMeasCleanCpu_ +=
            mesh_.time().elapsedCpuTime() - cleanStageCpu0;

        profilePostCleanWall_ += elapsedWallSeconds(cleanWallStart);
        profilePostCleanCpu_ += mesh_.time().elapsedCpuTime() - cleanCpu0;
        resetDetailTimer();
    }

    if (profileSummary_)
    {
        profilePostFieldsWall_ += elapsedWallSeconds(wallStart);
        profilePostFieldsCpu_ += mesh_.time().elapsedCpuTime() - cpu0;
    }
}


bool Foam::dsmcCloud::noScheduledFieldOutput() const
{
    const Time& runTime = mesh_.time();

    if (runTime.outputTime())
    {
        return false;
    }

    const dictionary& dict = runTime.controlDict();
    const word writeControl =
        dict.lookupOrDefault<word>("writeControl", "timeStep");

    if (writeControl != "runTime" && writeControl != "adjustableRunTime")
    {
        // §31 revert: the conservative legacy behaviour is restored —
        // timeStep-based writing assumes output may happen, so the
        // per-step field calculation always runs.  The lazyFieldCalculation
        // experiment (precise interval-vs-span prediction, −19.6% on the
        // no-output test config) was removed: production runs always have
        // scheduled field output, so the optimisation never applies there,
        // while as a default it would silently skip the per-step fields
        // for any run without output.
        return false;
    }

    const scalar writeInterval =
        dict.lookupOrDefault<scalar>("writeInterval", GREAT);
    const scalar runSpan =
        runTime.endTime().value() - runTime.startTime().value();

    return writeInterval > runSpan + 0.5*runTime.deltaT().value();
}


bool Foam::dsmcCloud::collisionModelUsesMacroscopicFieldTemperature() const
{
    if (!binaryCollisionModel_.valid())
    {
        return false;
    }

    const dictionary& coeffs = binaryCollisionModel_->coeffDict();
    const word inverseZvFormulation =
        coeffs.lookupOrDefault<word>("inverseZvFormulation", word::null);

    return inverseZvFormulation == "2008";
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
        moveOrderedThreadOffsets_[threadI] = label
        (
            (
                static_cast<long long>(threadI)
              * static_cast<long long>(i)
            )
          / static_cast<long long>(nThreads)
        );
    }
    moveOrderedParcelsValid_ = true;
    moveAppendedParcels_.clear();
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;
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
    moveAppendedParcels_.clear();
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;
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
    moveAppendedParcels_.clear();
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;
}


void Foam::dsmcCloud::beginMoveAppendCapture()
{
    ++moveAppendCaptureDepth_;
    moveAppendCaptureActive_ = true;
}


void Foam::dsmcCloud::endMoveAppendCapture()
{
    if (moveAppendCaptureDepth_ > 0)
    {
        --moveAppendCaptureDepth_;
    }

    moveAppendCaptureActive_ = moveAppendCaptureDepth_ > 0;
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
    // §34.2: ownership-transfer fill.  List::setSize() allocates a fresh
    // array AND copies min(old,new) old elements (List.C:283 — no
    // capacity retention), which the fill loop below overwrites
    // completely: up to ~136 MB of dead memcpy per step at 17M parcels.
    // Fill a fresh list and transfer() its storage (O(1) swap) instead.
    List<dsmcParcel*> tmp(parcels.size());
    forAll(parcels, i)
    {
        tmp[i] = parcels[i];
    }
    moveOrderedParcels_.transfer(tmp);

    const label nThreads = max(ompNumThreads_, label(1));
    moveOrderedThreadOffsets_.setSize(nThreads + 1);

    for (label threadI = 0; threadI <= nThreads; ++threadI)
    {
        moveOrderedThreadOffsets_[threadI] = label
        (
            (
                static_cast<long long>(threadI)
              * static_cast<long long>(moveOrderedParcels_.size())
            )
          / static_cast<long long>(nThreads)
        );
    }
    moveOrderedParcelsValid_ = true;
    moveAppendedParcels_.clear();
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;
}


void Foam::dsmcCloud::takeMoveOrderedParcels
(
    DynamicList<dsmcParcel*>& parcels
)
{
    // §36.1: zero-copy ownership transfer.  The caller's DynamicList was
    // filled by index during the migration pass; List::transfer() steals
    // its storage (v_/size_) instead of re-walking it into a fresh array
    // (the previous path copied ~136 MB of pointers per step at 17M
    // parcels, on top of the fill itself).  The DynamicList is left with
    // v_=nullptr and must not be reused by the caller.
    moveOrderedParcels_.transfer(parcels);

    const label nThreads = max(ompNumThreads_, label(1));
    moveOrderedThreadOffsets_.setSize(nThreads + 1);

    for (label threadI = 0; threadI <= nThreads; ++threadI)
    {
        moveOrderedThreadOffsets_[threadI] = label
        (
            (
                static_cast<long long>(threadI)
              * static_cast<long long>(moveOrderedParcels_.size())
            )
          / static_cast<long long>(nThreads)
        );
    }
    moveOrderedParcelsValid_ = true;
    moveAppendedParcels_.clear();
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;
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
        moveOrderedThreadOffsets_[threadI] = label
        (
            (
                static_cast<long long>(threadI)
              * static_cast<long long>(moveOrderedParcels_.size())
            )
          / static_cast<long long>(nThreads)
        );
    }

    moveOrderedParcelsValid_ = true;
    moveAppendedParcels_.clear();
    occupancyOrderedParcelsValid_ = false;
    cellOccupancyMaterialized_ = false;
}


void Foam::dsmcCloud::deleteParcel(dsmcParcel* p)
{
    if (p)
    {
        deleteParticle(*p);
    }
}


void Foam::dsmcCloud::setLastStepWallTime(const scalar t)
{
    if (replicatedMeshActive())
    {
        replicatedMesh_->setLastStepWallTime(t);
    }
}

// §45 P0: 碰撞相每 cell 候选快照 → 复制网格的窗口观测累加
void Foam::dsmcCloud::accumulateCellCostCandidates
(
    const labelList& counts,
    const labelList& touched
)
{
    if (replicatedMeshActive())
    {
        replicatedMesh_->accumulateCellCostCandidates(counts, touched);
    }
}


// comm_further_opt.md §15.4: cumulative counters accessors
double Foam::dsmcCloud::moveDetailCumulative(const label idx) const
{
    return (idx >= 0 && idx < 12) ? dsmcMoveDetailCum[idx] : 0.0;
}

void Foam::dsmcCloud::addToCumulativeCollisions
(
    const label nColl,
    const label nCand
)
{
    dsmcCumCollisions_ += nColl;
    dsmcCumCollisionCandidates_ += nCand;
}

double Foam::dsmcCloud::cumulativeCollisions() const
{
    return dsmcCumCollisions_;
}

double Foam::dsmcCloud::cumulativeCollisionCandidates() const
{
    return dsmcCumCollisionCandidates_;
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
    scalar postReactionsWall = profilePostReactionsWall_;
    scalar postFieldCalcWall = profilePostFieldCalcWall_;
    scalar postFieldWriteWall = profilePostFieldWriteWall_;
    scalar postControllersWall = profilePostControllersWall_;
    scalar postBoundariesWall = profilePostBoundariesWall_;
    scalar postBoundaryMeasWall = profilePostBoundaryMeasWall_;
    scalar postCleanWall = profilePostCleanWall_;
    scalar postTrackingCleanWall = profilePostTrackingCleanWall_;
    scalar postBoundaryMeasCleanWall = profilePostBoundaryMeasCleanWall_;
    scalar postCellMeasCleanWall = profilePostCellMeasCleanWall_;

    const bool replicatedRawMpi =
        replicatedMeshActive() && replicatedMesh_->nProcs() > 1;
    const bool replicatedRawMpiOutput =
        replicatedRawMpi && replicatedMesh_->myRank() == 0;
    bool replicatedRawMpiInitialized = false;
    int replicatedRawMpiSize = 1;

    const label localProfileSteps = profileSteps_;
    const scalar localFullEvolveWall = profileFullEvolveWall_;
    const scalar localMoveWall = profileMoveWall_;
    const scalar localBuildCellOccupancyWall = profileBuildCellOccupancyWall_;
    const scalar localCollisionWall = profileCollisionWall_;
    const scalar localPostFieldsWall = profilePostFieldsWall_;
    const label localParcels = this->size();
    const label localOwnedCollisionCells = occupancyOwnedCollisionCells_.size();
    label localFinalCandidates = 0;
    forAll(nCandidatesPerCell_, cellI)
    {
        localFinalCandidates += nCandidatesPerCell_[cellI];
    }

    // §31: residual accounting.  The phase values below are per-phase
    // MAXes over ranks (possibly from different ranks), so
    // max(evolve) - sum(max(phase)) is not any rank's residual.  Compute
    // each rank's own residual and reduce MAX, together with the
    // residual sub-items, so the whole section shares one convention.
    const scalar localCoreWall =
        localMoveWall
      + localBuildCellOccupancyWall
      + localCollisionWall
      + localPostFieldsWall;
    scalar maxResidualWall = localFullEvolveWall - localCoreWall;
    scalar maxResidualUpdateInfo = dsmcEvolveUpdateInfoWall_;
    scalar maxResidualCoord = dsmcEvolveCoordWall_;
    scalar maxResidualCtrl = dsmcEvolveCtrlWall_;
    // Guarded: pure-OMP runs have no replicated mesh allocated.
    scalar maxResidualUpc =
        replicatedMeshActive()
      ? replicatedMesh_->updateParticleCountsWallTime()
      : 0.0;

    // §33: per-rank residual table.  max(E) and max(R) may come from
    // different ranks, so max E - max R is not any rank's phase sum;
    // this gathers each rank's own (E, C, R) plus the residual's known
    // residents (DLB migration wall, updateParticleCounts) so the
    // attribution is per-rank and unambiguous.
    List<scalar> perRankResidualRows;

    // reactionOccupancyAppend profile values: times (max over rank), rank-max
    // counters and rank-sum counters.  Local values unless raw-MPI reduces.
    double rxnTime[6] =
    {
        dsmcRxn_.tAppend, dsmcRxn_.tSort, dsmcRxn_.tShift,
        dsmcRxn_.tOffsets, dsmcRxn_.tDerived, dsmcRxn_.tValidate
    };
    double rxnMaxCnt[6] =
    {
        dsmcRxn_.nCalls, dsmcRxn_.nAppend, dsmcRxn_.nFbThreshold,
        dsmcRxn_.nFbGuard, dsmcRxn_.nSkipZero, dsmcRxn_.nLate
    };
    double rxnSumCnt[8] =
    {
        dsmcRxn_.nParcels, dsmcRxn_.nShapeSingle, dsmcRxn_.nShapeMulti,
        dsmcRxn_.nShapeSameCell, dsmcRxn_.nKeptOneCells,
        dsmcRxn_.nKeptZeroCells, dsmcRxn_.nShapeLast, dsmcRxn_.nValidated
    };

    if (replicatedRawMpi)
    {
        int mpiInit = 0;
        MPI_Initialized(&mpiInit);

        if (mpiInit)
        {
            replicatedRawMpiInitialized = true;
            MPI_Comm_size(MPI_COMM_WORLD, &replicatedRawMpiSize);

            MPI_Allreduce
            (
                MPI_IN_PLACE,
                &steps,
                1,
                MPI_INT,
                MPI_MAX,
                MPI_COMM_WORLD
            );

            const int nProfileValues = 16;
            scalar profileValues[nProfileValues] =
            {
                fullEvolveWall,
                moveAndCollideWall,
                moveWall,
                buildCellOccupancyWall,
                collisionWall,
                postFieldsWall,
                postReactionsWall,
                postFieldCalcWall,
                postFieldWriteWall,
                postControllersWall,
                postBoundariesWall,
                postBoundaryMeasWall,
                postCleanWall,
                postTrackingCleanWall,
                postBoundaryMeasCleanWall,
                postCellMeasCleanWall
            };

            MPI_Allreduce
            (
                MPI_IN_PLACE,
                profileValues,
                nProfileValues,
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
            postReactionsWall = profileValues[6];
            postFieldCalcWall = profileValues[7];
            postFieldWriteWall = profileValues[8];
            postControllersWall = profileValues[9];
            postBoundariesWall = profileValues[10];
            postBoundaryMeasWall = profileValues[11];
            postCleanWall = profileValues[12];
            postTrackingCleanWall = profileValues[13];
            postBoundaryMeasCleanWall = profileValues[14];
            postCellMeasCleanWall = profileValues[15];

            // §31: residual + sub-items reduced together (one call).
            scalar residualValues[5] =
            {
                maxResidualWall,
                maxResidualUpdateInfo,
                maxResidualCoord,
                maxResidualCtrl,
                maxResidualUpc
            };
            MPI_Allreduce
            (
                MPI_IN_PLACE,
                residualValues,
                5,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD
            );
            maxResidualWall = residualValues[0];
            maxResidualUpdateInfo = residualValues[1];
            maxResidualCoord = residualValues[2];
            maxResidualCtrl = residualValues[3];
            maxResidualUpc = residualValues[4];

            // §33: gather the per-rank rows (E, C, R, DLB migration wall,
            // DLB count, updateParticleCounts).
            double myResidualRow[6] =
            {
                double(localFullEvolveWall),
                double(localCoreWall),
                double(localFullEvolveWall - localCoreWall),
                double(replicatedMesh_->autoRebalanceMigrationWallTime()),
                double(replicatedMesh_->autoRebalanceCount()),
                double(replicatedMesh_->updateParticleCountsWallTime())
            };
            perRankResidualRows.setSize(6*replicatedRawMpiSize, 0.0);
            MPI_Allgather
            (
                myResidualRow, 6, MPI_DOUBLE,
                perRankResidualRows.data(), 6, MPI_DOUBLE,
                MPI_COMM_WORLD
            );

            // buildCellOccupancy subphase diagnostics (§6.4)
            scalar buildOccValues[4] =
            {
                buildOccResetWall,
                buildOccCountWall,
                buildOccReduceWall,
                buildOccFillWall
            };
            MPI_Allreduce
            (
                MPI_IN_PLACE,
                buildOccValues,
                4,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD
            );
            buildOccResetWall = buildOccValues[0];
            buildOccCountWall = buildOccValues[1];
            buildOccReduceWall = buildOccValues[2];
            buildOccFillWall = buildOccValues[3];

            label buildOccHits[2] =
            {
                buildOccOrderedHits,
                buildOccFallbacks
            };
            MPI_Allreduce
            (
                MPI_IN_PLACE,
                buildOccHits,
                2,
                MPI_INT,
                MPI_MAX,
                MPI_COMM_WORLD
            );
            buildOccOrderedHits = buildOccHits[0];
            buildOccFallbacks = buildOccHits[1];

            MPI_Allreduce
            (
                MPI_IN_PLACE, rxnTime, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD
            );
            MPI_Allreduce
            (
                MPI_IN_PLACE, rxnMaxCnt, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD
            );
            MPI_Allreduce
            (
                MPI_IN_PLACE, rxnSumCnt, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD
            );
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
        reduce(postReactionsWall, maxOp<scalar>());
        reduce(postFieldCalcWall, maxOp<scalar>());
        reduce(postFieldWriteWall, maxOp<scalar>());
        reduce(postControllersWall, maxOp<scalar>());
        reduce(postBoundariesWall, maxOp<scalar>());
        reduce(postBoundaryMeasWall, maxOp<scalar>());
        reduce(postCleanWall, maxOp<scalar>());
        reduce(postTrackingCleanWall, maxOp<scalar>());
        reduce(postBoundaryMeasCleanWall, maxOp<scalar>());
        reduce(postCellMeasCleanWall, maxOp<scalar>());
    }

    if (profileDetail_ && replicatedRawMpi && replicatedRawMpiInitialized)
    {
        const int nScalarDetail = 5;
        const int nLabelDetail = 4;

        scalar localScalarDetail[nScalarDetail] =
        {
            localFullEvolveWall,
            localMoveWall,
            localBuildCellOccupancyWall,
            localCollisionWall,
            localPostFieldsWall
        };

        label localLabelDetail[nLabelDetail] =
        {
            localProfileSteps,
            localParcels,
            localOwnedCollisionCells,
            localFinalCandidates
        };

        List<scalar> allScalarDetail;
        List<label> allLabelDetail;

        if (replicatedRawMpiOutput)
        {
            allScalarDetail.setSize(replicatedRawMpiSize*nScalarDetail, 0.0);
            allLabelDetail.setSize(replicatedRawMpiSize*nLabelDetail, 0);
        }

        MPI_Gather
        (
            localScalarDetail,
            nScalarDetail,
            MPI_DOUBLE,
            replicatedRawMpiOutput ? allScalarDetail.data() : nullptr,
            nScalarDetail,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        MPI_Gather
        (
            localLabelDetail,
            nLabelDetail,
            MPI_INT,
            replicatedRawMpiOutput ? allLabelDetail.data() : nullptr,
            nLabelDetail,
            MPI_INT,
            0,
            MPI_COMM_WORLD
        );

        if (replicatedRawMpiOutput)
        {
            scalar minRankFull = GREAT;
            scalar maxRankFull = 0.0;
            scalar minRankCollision = GREAT;
            scalar maxRankCollision = 0.0;

            Info<< nl
                << "Replicated mesh profile detail by rank:" << nl
                << "    rank steps parcels ownedCollCells finalCandidates"
                << " full move build_occupancy coll post" << nl;

            for (int rankI = 0; rankI < replicatedRawMpiSize; ++rankI)
            {
                const label labelBase = rankI*nLabelDetail;
                const label scalarBase = rankI*nScalarDetail;
                const scalar rankFull = allScalarDetail[scalarBase + 0];
                const scalar rankCollision = allScalarDetail[scalarBase + 3];

                minRankFull = min(minRankFull, rankFull);
                maxRankFull = max(maxRankFull, rankFull);
                minRankCollision = min(minRankCollision, rankCollision);
                maxRankCollision = max(maxRankCollision, rankCollision);

                Info<< "    rank" << rankI
                    << " " << allLabelDetail[labelBase + 0]
                    << " " << allLabelDetail[labelBase + 1]
                    << " " << allLabelDetail[labelBase + 2]
                    << " " << allLabelDetail[labelBase + 3]
                    << " " << rankFull
                    << " " << allScalarDetail[scalarBase + 1]
                    << " " << allScalarDetail[scalarBase + 2]
                    << " " << rankCollision
                    << " " << allScalarDetail[scalarBase + 4]
                    << nl;
            }

            Info<< "    rank full max/min            = "
                << maxRankFull/max(minRankFull, SMALL) << nl
                << "    rank collision max/min       = "
                << maxRankCollision/max(minRankCollision, SMALL) << nl;
        }
    }

    if
    (
        (replicatedRawMpi && replicatedRawMpiOutput)
     || (!replicatedRawMpi && Pstream::master())
    )
    {
        const scalar coreAccountedWall =
            moveWall
          + buildCellOccupancyWall
          + collisionWall
          + postFieldsWall;
        const scalar postDetailWall =
            postReactionsWall
          + postFieldCalcWall
          + postFieldWriteWall
          + postControllersWall
          + postBoundariesWall
          + postBoundaryMeasWall
          + postCleanWall;

        Info<< nl
            << "DSMC cloud profile v2" << nl
            << "    profile steps                 = " << steps << nl
            << "    full evolve max [s]           = " << fullEvolveWall << nl
            << "    move+collide max [s]          = " << moveAndCollideWall << nl
            << "    move max [s]                  = " << moveWall << nl
            << "    build_occupancy max [s]       = "
            << buildCellOccupancyWall << nl
            << "    build_occupancy reset max [s] = "
            << buildOccResetWall << nl
            << "    build_occupancy count max [s] = "
            << buildOccCountWall << nl
            << "    build_occupancy reduce max [s] = "
            << buildOccReduceWall << nl
            << "    build_occupancy fill max [s]  = "
            << buildOccFillWall << nl
            << "    build_occupancy ordered/fallback = "
            << buildOccOrderedHits << "/" << buildOccFallbacks << nl
            << (
                dsmcRxnAppendRead_ && dsmcRxnAppendEnabled_
              ? "    build_occupancy reaction append max [s] = "
              + std::to_string(rxnTime[0]) + "  (included in build_occupancy)\n"
              + "      reaction append calls/append/fb-threshold/fb-guard/skip-zero/late (rank max) = "
              + std::to_string(label(rxnMaxCnt[0])) + "/"
              + std::to_string(label(rxnMaxCnt[1])) + "/"
              + std::to_string(label(rxnMaxCnt[2])) + "/"
              + std::to_string(label(rxnMaxCnt[3])) + "/"
              + std::to_string(label(rxnMaxCnt[4])) + "/"
              + std::to_string(label(rxnMaxCnt[5])) + "\n"
              + "      reaction append parcels (rank sum) = "
              + std::to_string(label(rxnSumCnt[0])) + "\n"
              + "      reaction append shapes (rank sum): single/multi/sameCellMulti/keptOne cells/keptZero cells/lastCell = "
              + std::to_string(label(rxnSumCnt[1])) + "/"
              + std::to_string(label(rxnSumCnt[2])) + "/"
              + std::to_string(label(rxnSumCnt[3])) + "/"
              + std::to_string(label(rxnSumCnt[4])) + "/"
              + std::to_string(label(rxnSumCnt[5])) + "/"
              + std::to_string(label(rxnSumCnt[6])) + "\n"
              + "      reaction append subphase max [s]: sort = "
              + std::to_string(rxnTime[1]) + " shift = "
              + std::to_string(rxnTime[2]) + " offsets = "
              + std::to_string(rxnTime[3]) + " derived = "
              + std::to_string(rxnTime[4]) + "\n"
              + "      reaction append validate: checked (rank sum) = "
              + std::to_string(label(rxnSumCnt[7]))
              + " mismatches = 0 (any mismatch aborts) wall max [s] = "
              + std::to_string(rxnTime[5]) + " (NOT in build_occupancy)\n"
              : std::string("")
              ).c_str()
            << "    coll max [s]                  = " << collisionWall << nl
            << "    post max [s]                  = " << postFieldsWall << nl
            << "    core phase sum [s]            = " << coreAccountedWall << nl
            << "    evolve residual (per-rank max) [s] = "
            << maxResidualWall << nl
            << "      residual: model info update (max) [s] = "
            << maxResidualUpdateInfo << nl
            << "      residual: coordSystem evolve (max) [s] = "
            << maxResidualCoord << nl
            << "      residual: controlBeforeCollisions (max) [s] = "
            << maxResidualCtrl << nl
            << "      residual: updateParticleCounts (max) [s] = "
            << maxResidualUpc << nl
            << "    per-rank residual table (E=evolve, C=4-phase sum,"
            << " R=E-C, dlbMig=DLB migration wall, nDLB=executed,"
            << " upc=updateParticleCounts):" << nl;
            if (perRankResidualRows.size() == 6*replicatedRawMpiSize)
            {
                label argmaxR = 0;
                for (label r = 0; r < replicatedRawMpiSize; ++r)
                {
                    if
                    (
                        perRankResidualRows[6*r + 2]
                      > perRankResidualRows[6*argmaxR + 2]
                    )
                    {
                        argmaxR = r;
                    }
                }
                for (label r = 0; r < replicatedRawMpiSize; ++r)
                {
                    Info<< "      r" << r
                        << " E=" << perRankResidualRows[6*r + 0]
                        << " C=" << perRankResidualRows[6*r + 1]
                        << " R=" << perRankResidualRows[6*r + 2]
                        << " dlbMig=" << perRankResidualRows[6*r + 3]
                        << " nDLB=" << label(perRankResidualRows[6*r + 4])
                        << " upc=" << perRankResidualRows[6*r + 5]
                        << nl;
                }
                Info<< "      argmax-R rank = r" << argmaxR
                    << " (E=" << perRankResidualRows[6*argmaxR + 0]
                    << ", C=" << perRankResidualRows[6*argmaxR + 1]
                    << ", R=" << perRankResidualRows[6*argmaxR + 2]
                    << ", dlbMig=" << perRankResidualRows[6*argmaxR + 3]
                    << ", upc=" << perRankResidualRows[6*argmaxR + 5]
                    << ", R - dlbMig - upc = "
                    << (perRankResidualRows[6*argmaxR + 2]
                      - perRankResidualRows[6*argmaxR + 3]
                      - perRankResidualRows[6*argmaxR + 5])
                    << ")" << nl;
            }
            Info<< "    post definition               = dsmcCloud::evolve_fields()" << nl
            << "    post reactions [s]            = "
            << postReactionsWall << nl
            << "    post field calculate [s]      = "
            << postFieldCalcWall << nl
            << "    post field write [s]          = "
            << postFieldWriteWall << nl
            << "    post controllers [s]          = "
            << postControllersWall << nl
            << "    post boundaries [s]           = "
            << postBoundariesWall << nl
            << "    post boundary meas [s]        = "
            << postBoundaryMeasWall << nl
            << "    post clean [s]                = "
            << postCleanWall << nl
            << "    post tracking clean [s]      = "
            << postTrackingCleanWall << nl
            << "    post boundary meas clean [s] = "
            << postBoundaryMeasCleanWall << nl
            << "    post cell meas clean [s]     = "
            << postCellMeasCleanWall << nl
            << "    post subphase sum [s]         = "
            << postDetailWall << nl
            << "    post residual [s]             = "
            << postFieldsWall - postDetailWall << nl
            << "    OpenMP enabled                = " << openmpEnabled_ << nl
            << "    OpenMP max threads            = " << ompNumThreads_ << nl
            << "    OpenMP move                   = " << openmpMoveEnabled_
            << " (" << openmpMoveSchedule_ << ", chunk "
            << openmpMoveChunk_ << ")" << nl
            << "    OpenMP collision              = "
            << openmpCollisionSchedule_ << ", chunk "
            << openmpCollisionChunk_ << nl
            << endl;

        if
        (
            profileDetail_
         && controlDict_.lookupOrDefault<bool>
            (
                "moveDetailProfile",
                false
            )
        )
        {
            Info<< "DSMC move detail profile (thread-summed diagnostics)" << nl
                << "    parcels [count]              = "
                << dsmcMoveDetailCum[0] << nl
                << "    track calls [count]          = "
                << dsmcMoveDetailCum[1] << nl
                << "    same-tet no-face [count]     = "
                << dsmcMoveDetailCum[2] << nl
                << "    internal-tet [count]         = "
                << dsmcMoveDetailCum[3] << nl
                << "    face hits [count]            = "
                << dsmcMoveDetailCum[4] << nl
                << "    processor hits [count]       = "
                << dsmcMoveDetailCum[5] << nl
                << "    cyclic hits [count]          = "
                << dsmcMoveDetailCum[6] << nl
                << "    patch hits [count]           = "
                << dsmcMoveDetailCum[7] << nl
                << "    stuck hits [count]           = "
                << dsmcMoveDetailCum[8] << nl
                << "    trackToFace sum [s]           = "
                << dsmcMoveDetailCum[9] << nl
                << "    tracker sum [s]              = "
                << dsmcMoveDetailCum[10] << nl
                << "    boundary sum [s]             = "
                << dsmcMoveDetailCum[11] << nl
                << endl;
        }
    }

    // Worklog dlb_further_opt.md §11: per-rank move sub-phase table.
    // Under replicated mesh each rank is a spatial region, so per-rank
    // counters locate the source of the measured ~1.5x per-parcel move
    // cost variation (tracking crossing rate vs boundary critical-section
    // cost vs thread imbalance).
    // NOTE: this block must be OUTSIDE the rank0-gated section above so
    // that every rank reaches the MPI_Gather (same collective-communication
    // placement rule as the §5.6 report() Allreduce fix).
    if
    (
        profileDetail_
     && controlDict_.lookupOrDefault<bool>
        (
            "moveDetailProfile",
            false
        )
     && replicatedRawMpi
     && replicatedRawMpiInitialized
    )
    {
        const int nMoveVals = 14;
        scalar localMoveVals[nMoveVals] =
        {
            dsmcMoveDetailCum[0],
            dsmcMoveDetailCum[1],
            dsmcMoveDetailCum[4],
            dsmcMoveDetailCum[6],
            dsmcMoveDetailCum[7],
            dsmcMoveDetailCum[5],
            dsmcMoveDetailCum[9],
            dsmcMoveDetailCum[10],
            dsmcMoveDetailCum[11],
            localMoveWall,
            dsmcMoveDetailSplitCum[0],
            dsmcMoveDetailSplitCum[1],
            dsmcMoveDetailSplitCum[2],
            dsmcMoveDetailSplitCum[3]
        };

        List<scalar> allMoveVals;
        if (replicatedRawMpiOutput)
        {
            allMoveVals.setSize(replicatedRawMpiSize*nMoveVals, 0.0);
        }
        MPI_Gather
        (
            localMoveVals,
            nMoveVals,
            MPI_DOUBLE,
            replicatedRawMpiOutput ? allMoveVals.data() : nullptr,
            nMoveVals,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        if (replicatedRawMpiOutput)
        {
            Info<< nl
                << "Move sub-phase by rank (spatial cost split):" << nl
                << "    rank parcels trackCalls tc/parcel"
                << " faceHits patchHits trackWall moveWall"
                << " | wallCellCalls wallCell_s wCall_us"
                << " | intCalls int_s iCall_us" << nl;

            scalar minTcPerParcel = GREAT;
            scalar maxTcPerParcel = 0.0;
            for (int rankI = 0; rankI < replicatedRawMpiSize; ++rankI)
            {
                const scalar base = rankI*nMoveVals;
                const scalar parcels = allMoveVals[base + 0];
                const scalar trackCalls = allMoveVals[base + 1];
                const scalar tcPerParcel =
                    parcels > 0.5 ? trackCalls/parcels : 0.0;
                minTcPerParcel = min(minTcPerParcel, tcPerParcel);
                maxTcPerParcel = max(maxTcPerParcel, tcPerParcel);

                const scalar wallCalls = allMoveVals[base + 10];
                const scalar wallTime = allMoveVals[base + 11];
                const scalar intCalls = allMoveVals[base + 12];
                const scalar intTime = allMoveVals[base + 13];
                const scalar wallPerCall_us =
                    wallCalls > 0.5 ? 1.0e6*wallTime/wallCalls : 0.0;
                const scalar intPerCall_us =
                    intCalls > 0.5 ? 1.0e6*intTime/intCalls : 0.0;

                // Counts printed as scalars: Info<< has no long overload
                // and a long here silently degrades to label (int32),
                // wrapping counts above 2^31 (worklog §11).
                Info<< "    rank" << rankI
                    << " " << parcels
                    << " " << trackCalls
                    << " " << tcPerParcel
                    << " " << allMoveVals[base + 2]
                    << " " << allMoveVals[base + 4]
                    << " " << allMoveVals[base + 6]
                    << " " << allMoveVals[base + 9]
                    << " | " << wallCalls
                    << " " << wallTime
                    << " " << wallPerCall_us
                    << " | " << intCalls
                    << " " << intTime
                    << " " << intPerCall_us
                    << nl;
            }

            Info<< "    rank trackCalls/parcel max/min = "
                << maxTcPerParcel/max(minTcPerParcel, SMALL) << nl
                << endl;
        }
    }

    if (replicatedMeshActive())
    {
        replicatedMesh_->report();
    }

    return;
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


void Foam::dsmcCloud::setCollisionRngContext(dsmcFastRng* rng)
{
    collisionRngContext = rng;
}


void Foam::dsmcCloud::clearCollisionRngContext()
{
    collisionRngContext = nullptr;
}


Foam::scalar Foam::dsmcCloud::collisionSample01()
{
    if (collisionRngContext)
    {
        return collisionRngContext->sample01();
    }

    #ifdef _OPENMP
    if (omp_in_parallel())
    {
        scalar value = 0.0;
        #pragma omp critical(dsmcCollisionRndGen)
        {
            value = rndGen_.sample01<scalar>();
        }
        return value;
    }
    #endif

    return rndGen_.sample01<scalar>();
}


Foam::label Foam::dsmcCloud::collisionRandomLabel
(
    const label valOne,
    const label valTwo
)
{
    if (valOne == valTwo)
    {
        return valOne;
    }

    const label start = Foam::min(valOne, valTwo);
    const label end = Foam::max(valOne, valTwo);
    const label n = end - start + 1;

    if (collisionRngContext)
    {
        return start + collisionRngContext->position(n);
    }

    label val = start + label(collisionSample01()*n);
    if (val > end)
    {
        val = end;
    }

    return val;
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


Foam::vector Foam::dsmcCloud::equipartitionLinearVelocity
(
    const scalar temperature,
    const scalar mass,
            Random& random
)
{
    return sqrt(physicoChemical::k.value()*temperature/mass)
        *random.GaussNormal<vector>();
}

Foam::scalar Foam::dsmcCloud::equipartitionRotationalEnergy
(
    const scalar temperature,
    const scalar rotationalDof,
            Random& random
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
        ERot = -log(random.sample01<scalar>())*physicoChemical::k.value()*temperature;
    }
    else
    {
        scalar a = 0.5*rotationalDof - 1;

        scalar energyRatio;

        scalar P = -1;

        do
        {
            energyRatio = 10*random.sample01<scalar>();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < random.sample01<scalar>());

        ERot = energyRatio*physicoChemical::k.value()*temperature;
    }

    return ERot;
}

Foam::labelList Foam::dsmcCloud::equipartitionVibrationalEnergyLevel
(
    const scalar temperature,
    const label nVibrationalModes,
    const label typeId,
            Random& random
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
            vibLevel[mode] = -log(random.sample01<scalar>())*temperature
                /constProps(typeId).thetaV_m(mode);
        }
    }

    return vibLevel;
}

Foam::label Foam::dsmcCloud::equipartitionElectronicLevel
(
    const scalar temperature,
    const labelList& electronicDegeneracyList,
    const scalarList& electronicEnergyList,
            Random& random
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
          //jDash = random.position<label>(0,jMax); OLD
            jDash = randomLabel(0, jMax);
            func =
                electronicDegeneracyList[jDash]
               *exp(-electronicEnergyList[jDash]/EMax)
               /expMax;
        } while(func < random.sample01<scalar>());
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
        energyRatio = 1.0 - pow(collisionSample01(), 1.0/ChiB);
    }
    else
    {
        const scalar ChiA = 0.5*rotationalDof;

        scalar ChiAMinusOne = ChiA - 1.;

        scalar ChiBMinusOne = ChiB - 1.;

        if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
        {
            return collisionSample01();
        }

        scalar P = 0.0;

        do
        {
            P = 0;

            energyRatio = collisionSample01();

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
        } while (P < collisionSample01());
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
            iDash = collisionRandomLabel(0, iMax);
            EVib = iDash*physicoChemical::k.value()*thetaV;

            // - equation 5.61, Bird
            func = pow(1.0 - EVib/Ec, 1.5 - omega);

        } while(func < collisionSample01());
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

        if (inverseVibrationalCollisionNumber > collisionSample01())
        {
            // post-collision quantum number
            scalar func = 0.0;
            scalar EVib = 0.0;

            do // acceptance - rejection
            {
                iDash = collisionRandomLabel(0, iMax);

                EVib = iDash*physicoChemical::k.value()*thetaV;

                // - equation 5.61, Bird
                func = pow(1.0 - EVib/Ec, 1.5 - omega);

            } while(func < collisionSample01());
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
        const label nState = collisionRandomLabel(1, nPossibleStates);
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

            if (prob > collisionSample01())
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
       jDash = collisionRandomLabel(0, jSelectA);
       prob = gList[jDash]*pow(Ec - EElist[jDash], 1.5 - omega)/denomMax;

    } while(prob < collisionSample01());

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
