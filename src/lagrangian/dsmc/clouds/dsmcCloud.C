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
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "meshTools.H"

using namespace Foam::constant;

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<dsmcParcel>, 0);
};


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
    const label nCells = cellOccupancy_.size();

    // If cell index is available and up-to-date, use it for parallel rebuild
    if (cellFirst_.size() == nCells && parcelPtrs_.size() > 0)
    {
        #pragma omp parallel for schedule(static)
        for (label cellI = 0; cellI < nCells; cellI++)
        {
            // Pre-size to exact count (avoids realloc during fill)
            cellOccupancy_[cellI].setSize(cellCount_[cellI]);

            label idx = 0;
            for (label ip = cellFirst_[cellI]; ip >= 0; ip = cellNext_[ip])
            {
                cellOccupancy_[cellI][idx++] = parcelPtrs_[ip];
            }
        }
    }
    else
    {
        // Fallback: serial rebuild from IDLList
        #pragma omp parallel for schedule(static)
        for (label i = 0; i < nCells; i++)
        {
            cellOccupancy_[i].clear();
        }

        forAllIter(dsmcCloud, *this, iter)
        {
            cellOccupancy_[iter().cell()].append(&iter());
        }
    }
}


void Foam::dsmcCloud::buildParcelPtrs()
{
    const label nPart = this->size();
    parcelPtrs_.setSize(nPart);

    // Serial IDLList traversal (unavoidable for linked list)
    label i = 0;
    forAllIter(dsmcCloud, *this, iter)
    {
        parcelPtrs_[i++] = &iter();
    }
}


void Foam::dsmcCloud::buildCellIndex()
{
    const label nCells = mesh_.nCells();
    const label nPart = parcelPtrs_.size();

    cellFirst_.setSize(nCells);
    cellCount_.setSize(nCells);
    cellNext_.setSize(nPart);

    // Step 1: parallel count particles per cell
    cellCount_ = 0;

    #pragma omp parallel
    {
        // Thread-local count
        labelList localCount(nCells, 0);

        #pragma omp for schedule(static)
        for (label i = 0; i < nPart; i++)
        {
            localCount[parcelPtrs_[i]->cell()]++;
        }

        // Merge into global count
        #pragma omp critical
        {
            for (label c = 0; c < nCells; c++)
            {
                cellCount_[c] += localCount[c];
            }
        }
    }

    // Step 2: serial build linked list (O(n), cache-friendly reverse scan)
    cellFirst_ = -1;
    cellNext_ = -1;

    for (label i = nPart - 1; i >= 0; --i)
    {
        const label cellI = parcelPtrs_[i]->cell();
        cellNext_[i] = cellFirst_[cellI];
        cellFirst_[cellI] = i;
    }
}

void Foam::dsmcCloud::validateCellIndex() const
{
    // Check total particle count
    label totalFromIndex = 0;
    forAll(cellCount_, c)
    {
        totalFromIndex += cellCount_[c];
    }

    if (totalFromIndex != parcelPtrs_.size())
    {
        FatalErrorInFunction
            << "Cell index total (" << totalFromIndex
            << ") != parcelPtrs size (" << parcelPtrs_.size() << ")"
            << exit(FatalError);
    }

    // Check against cellOccupancy
    forAll(cellOccupancy_, c)
    {
        if (cellCount_[c] != cellOccupancy_[c].size())
        {
            FatalErrorInFunction
                << "Cell " << c << ": cellCount=" << cellCount_[c]
                << " but cellOccupancy size=" << cellOccupancy_[c].size()
                << exit(FatalError);
        }
    }
}


Foam::Random& Foam::dsmcCloud::rng(const label tid)
{
    if (tid < 0)
    {
        return rndGen_;
    }
    return threadRng_[tid];
}


void Foam::dsmcCloud::buildBoundaryCellMarking()
{
    const label nCells = mesh_.nCells();
    boundaryCell_.setSize(nCells);
    boundaryCell_ = false;  // Force reset all to false

    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    Info<< "DSMC OpenMP move: scanning " << pbMesh.size()
        << " patches for boundary cells" << endl;

    forAll(pbMesh, patchI)
    {
        const polyPatch& pp = pbMesh[patchI];

        // Only mark cells adjacent to patches that cause particle interactions:
        // wall, processor, inlet/outlet (generic patch type).
        // Skip: empty, symmetry, symmetryPlane, wedge — these don't
        // cause boundary interactions in DSMC.
        const word& pType = pp.type();
        if
        (
            pp.size() > 0
         && pType != "empty"
         && pType != "symmetry"
         && pType != "symmetryPlane"
         && pType != "wedge"
        )
        {
            forAll(pp, faceI)
            {
                const label cellI = mesh_.faceOwner()[pp.start() + faceI];
                boundaryCell_[cellI] = true;
            }
        }
    }

    label nBoundary = 0;
    forAll(boundaryCell_, c)
    {
        if (boundaryCell_[c]) nBoundary++;
    }

    Info<< "DSMC OpenMP move: " << nBoundary << " boundary cells / "
        << nCells << " total ("
        << 100.0*scalar(nBoundary)/max(scalar(nCells), SMALL) << "%)"
        << endl;
}


void Foam::dsmcCloud::moveParallel(const scalar trackTime)
{
    // Phase 5: Parallel move with boundary cell classification
    // Internal particles (not in boundary cells) are moved in parallel.
    // Boundary particles are moved serially with full OF tracking.

    // Ensure boundary cell marking is initialized
    if (boundaryCell_.size() != mesh_.nCells())
    {
        buildBoundaryCellMarking();
    }

    // Force lazy-evaluated mesh data to be computed before entering
    // parallel region. Call trackToFace once in serial to trigger ALL
    // demand-driven mesh data that tracking might need.
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellCentres();
    (void)mesh_.cellVolumes();
    (void)mesh_.cells();
    (void)mesh_.faceOwner();
    (void)mesh_.faceNeighbour();
    (void)mesh_.faces();
    (void)mesh_.points();
    (void)mesh_.boundaryMesh();
    (void)this->cellHasWallFaces();

    // Warm up: do one serial tracking step to ensure all lazy data is computed
    if (parcelPtrs_.size() > 0)
    {
        dsmcParcel& firstP = *parcelPtrs_[0];
        dsmcParcel::trackingData warmupTd(*this);
        warmupTd.switchProcessor = false;
        warmupTd.keepParticle = true;
        // Track with zero displacement (no-op but triggers all lazy paths)
        firstP.trackToFace(firstP.position(), warmupTd, true);
    }

    // Pre-compute cell sizes for parallel use (avoid lazy mesh queries)
    const scalarField& cellVols = mesh_.cellVolumes();

    dsmcParcel::trackingData td(*this);

    buildParcelPtrs();
    const label nPart = parcelPtrs_.size();

    // Reset stepFraction for all particles (same as Cloud::move() does)
    // Parallelized since parcelPtrs_ provides random access
    #pragma omp parallel for schedule(static)
    for (label i = 0; i < nPart; i++)
    {
        parcelPtrs_[i]->stepFraction() = 0;
    }

    // Classify particles: internal vs boundary
    DynamicList<label> internalParcels(nPart);
    DynamicList<label> boundaryParcels;

    for (label i = 0; i < nPart; i++)
    {
        if (boundaryCell_[parcelPtrs_[i]->cell()])
        {
            boundaryParcels.append(i);
        }
        else
        {
            internalParcels.append(i);
        }
    }

    // --- Parallel move for internal particles ---
    // Use trackToFace for proper cell crossing, but only for internal particles.
    // cloud.labels() has been made local in particleTemplates.C (thread-safe).
    // Boundary measurements are skipped during parallel move (parallelMoveActive_).
    // Boundary particles are handled serially.

    parallelMoveActive_ = true;

    label nParallelMoved = 0;
    DynamicList<dsmcParcel*> parallelDeleteList;

    #pragma omp parallel reduction(+:nParallelMoved)
    {
        #ifdef _OPENMP
            const int tid = omp_get_thread_num();
        #else
            const int tid = 0;
        #endif

        DynamicList<dsmcParcel*> myDeleteList;

        #pragma omp for schedule(dynamic, 256)
        for (label ii = 0; ii < internalParcels.size(); ii++)
        {
            const label idx = internalParcels[ii];
            dsmcParcel& p = *parcelPtrs_[idx];

            dsmcParcel::trackingData localTd(*this);
            localTd.switchProcessor = false;
            localTd.keepParticle = true;

            if (p.isFree())
            {
                if (p.newParcel() != -1)
                {
                    p.stepFraction() = rng(tid).sample01<scalar>();
                    p.newParcel() = -1;
                }

                const label orgCell = p.cell();
                scalar tEnd =
                    (1.0 - p.stepFraction()) * deltaTValue(orgCell);
                vector Utracking = p.U();

                while (localTd.keepParticle && !localTd.switchProcessor
                       && tEnd > ROOTVSMALL)
                {
                    Utracking = p.U();

                    // 2D/reduced-D constraints (critical for 2D cases!)
                    if (coordSystem().type() == "dsmcCartesian")
                    {
                        meshTools::constrainToMeshCentre(mesh_, p.position());
                        meshTools::constrainDirection
                        (
                            mesh_, mesh_.solutionD(), Utracking
                        );
                    }

                    scalar dt = tEnd;
                    const label orgCellLocal = p.cell();
                    dt *= p.trackToFace
                    (
                        p.position() + dt*Utracking, localTd, true
                    );
                    tEnd -= dt;
                    p.stepFraction() = 1.0 - tEnd/deltaTValue(orgCellLocal);

                    if (p.onBoundary() && localTd.keepParticle)
                    {
                        if (isA<processorPolyPatch>
                            (mesh_.boundaryMesh()[p.patch(p.face())]))
                        {
                            localTd.switchProcessor = true;
                        }
                    }
                }

                nParallelMoved++;
            }

            if (!localTd.keepParticle)
            {
                myDeleteList.append(&p);
            }
        }

        #pragma omp critical
        {
            forAll(myDeleteList, i)
            {
                parallelDeleteList.append(myDeleteList[i]);
            }
        }
    } // end omp parallel

    parallelMoveActive_ = false;

    // Mark deleted particles with stepFraction = -1 (dead flag)
    // Actual deletion deferred to after boundary move
    forAll(parallelDeleteList, i)
    {
        parallelDeleteList[i]->stepFraction() = -1.0;
    }

    // Serial move for boundary particles
    forAll(boundaryParcels, bi)
    {
        dsmcParcel& p = *parcelPtrs_[boundaryParcels[bi]];
        if (p.stepFraction() < 0) continue; // skip if already dead
        td.switchProcessor = false;
        td.keepParticle = true;
        p.move(td, trackTime);
        if (!td.keepParticle)
        {
            p.stepFraction() = -1.0; // mark dead
        }
    }

    // Now delete all dead particles and compact parcelPtrs_
    label writeIdx = 0;
    for (label i = 0; i < parcelPtrs_.size(); i++)
    {
        if (parcelPtrs_[i]->stepFraction() < 0)
        {
            deleteParticle(*parcelPtrs_[i]);
        }
        else
        {
            parcelPtrs_[writeIdx++] = parcelPtrs_[i];
        }
    }
    parcelPtrs_.setSize(writeIdx);
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
    trackingInfo_.reset();
    fields_.resetFields();

    cellMeas_.reset();
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
    if (newParcel != -1)
    {
        tracker().trackFaceTransition(typeId, U, RWF, tetFaceI);
    }

    porousMeas().additionInteraction(*pPtr, newParcel);

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
        return rndGen().sample01<scalar>();
    }

    scalar energyRatio;

    scalar P;

    do
    {
        P = 0;

        energyRatio = rndGen().sample01<scalar>();

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
    } while (P < rndGen().sample01<scalar>());

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
        return rndGen().sample01<scalar>();
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
        rPSIm = rndGen().sample01<scalar>();
        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);
    } while (prob < rndGen().sample01<scalar>());

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
    porousMeasurements_(porousMeasurements::New(t, mesh, *this)),
    nTerminalOutputs_
    (
        controlDict_.lookupOrDefault<label>("nTerminalOutputs", 1)
    ),
    ompTimingEnabled_
    (
        controlDict_.lookupOrDefault<Switch>("ompTimingEnabled", true)
    ),
    parallelMoveActive_(false),
    timeMove_(0),
    timeBuildCellOccupancy_(0),
    timeCoordSystem_(0),
    timeCollisions_(0),
    timeFields_(0),
    timeControllersBeforeMove_(0),
    timeBoundariesBeforeMove_(0),
    timeControllersBeforeCollisions_(0),
    timeControllersAfterCollisions_(0),
    timeOther_(0),
    timingStepCounter_(0),
    timingOutputInterval_
    (
        controlDict_.lookupOrDefault<label>("timingOutputInterval", 100)
    ),
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

    // Initialize per-thread RNG
    {
        label nThreads = 1;
        #ifdef _OPENMP
            nThreads = omp_get_max_threads();
        #endif
        threadRng_.setSize(nThreads);
        for (label t = 0; t < nThreads; t++)
        {
            threadRng_.set
            (
                t,
                new Random(label(rndGen_.sample01<scalar>() * 1000000) + t + 1)
            );
        }
        Info<< "DSMC OpenMP: initialized " << nThreads
            << " thread RNGs" << endl;
    }

    // Build parallel data structures
    buildBoundaryCellMarking();
    buildParcelPtrs();
    buildCellIndex();

    collisionPartnerSelectionModel_ = autoPtr<collisionPartnerSelection>
    (
        collisionPartnerSelection::New(mesh, *this, particleProperties_)
    );

    collisionPartnerSelectionModel_->initialConfiguration();

    fields_.createFields();
    boundaryMeas_.setInitialConfig();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();
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
    porousMeasurements_(porousMeasurements::New(t, mesh, *this)),
    nTerminalOutputs_
    (
        controlDict_.lookupOrDefault<label>("nTerminalOutputs", 1)
    ),
    ompTimingEnabled_
    (
        controlDict_.lookupOrDefault<Switch>("ompTimingEnabled", true)
    ),
    parallelMoveActive_(false),
    timeMove_(0),
    timeBuildCellOccupancy_(0),
    timeCoordSystem_(0),
    timeCollisions_(0),
    timeFields_(0),
    timeControllersBeforeMove_(0),
    timeBoundariesBeforeMove_(0),
    timeControllersBeforeCollisions_(0),
    timeControllersAfterCollisions_(0),
    timeOther_(0),
    timingStepCounter_(0),
    timingOutputInterval_
    (
        controlDict_.lookupOrDefault<label>("timingOutputInterval", 100)
    ),
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
    evolve_moveAndCollide();
    evolve_fields();

    // --- Phase timing output ---
    if (ompTimingEnabled_)
    {
        timingStepCounter_++;

        if (timingStepCounter_ >= timingOutputInterval_)
        {
            const scalar N = scalar(timingStepCounter_);
            const scalar totalTime = timeMove_ + timeBuildCellOccupancy_
                + timeCoordSystem_ + timeCollisions_ + timeFields_
                + timeControllersBeforeMove_ + timeBoundariesBeforeMove_
                + timeControllersBeforeCollisions_
                + timeControllersAfterCollisions_;

            Info<< nl
                << "DSMC phase timings (last " << timingStepCounter_
                << " steps, wall-time seconds):" << nl
                << "  move                    " << timeMove_ << " ("
                << 100.0*timeMove_/max(totalTime, SMALL) << "%)" << nl
                << "  buildCellOccupancy      " << timeBuildCellOccupancy_ << " ("
                << 100.0*timeBuildCellOccupancy_/max(totalTime, SMALL) << "%)" << nl
                << "  coordSystem             " << timeCoordSystem_ << " ("
                << 100.0*timeCoordSystem_/max(totalTime, SMALL) << "%)" << nl
                << "  collisions              " << timeCollisions_ << " ("
                << 100.0*timeCollisions_/max(totalTime, SMALL) << "%)" << nl
                << "  fields                  " << timeFields_ << " ("
                << 100.0*timeFields_/max(totalTime, SMALL) << "%)" << nl
                << "  controllersBeforeMove   " << timeControllersBeforeMove_ << nl
                << "  boundariesBeforeMove    " << timeBoundariesBeforeMove_ << nl
                << "  controllersBefCollision " << timeControllersBeforeCollisions_ << nl
                << "  controllersAftCollision " << timeControllersAfterCollisions_ << nl
                << "  total                   " << totalTime << nl
                << "  per-step average        " << totalTime/N << nl
                << "DSMC diagnostics:" << nl
                << "  nParticles              " << this->size() << nl
                << "  nCells                  " << mesh_.nCells() << nl
                << endl;

            // Reset accumulators
            timeMove_ = 0;
            timeBuildCellOccupancy_ = 0;
            timeCoordSystem_ = 0;
            timeCollisions_ = 0;
            timeFields_ = 0;
            timeControllersBeforeMove_ = 0;
            timeBoundariesBeforeMove_ = 0;
            timeControllersBeforeCollisions_ = 0;
            timeControllersAfterCollisions_ = 0;
            timeOther_ = 0;
            timingStepCounter_ = 0;
        }
    }
}


void Foam::dsmcCloud::evolve_moveAndCollide()
{
    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    dsmcParcel::trackingData td(*this);

    if (debug)
    {
        this->dumpParticlePositions();
    }

    // --- Phase timing ---
    scalar t0 = 0, t1 = 0;
    #ifdef _OPENMP
        #define DSMC_WTIME() omp_get_wtime()
    #else
        #define DSMC_WTIME() mesh_.time().elapsedClockTime()
    #endif

    t0 = DSMC_WTIME();
    controllers_.controlBeforeMove();
    t1 = DSMC_WTIME();
    timeControllersBeforeMove_ += (t1 - t0);

    t0 = t1;
    boundaries_.controlBeforeMove();
    t1 = DSMC_WTIME();
    timeBoundariesBeforeMove_ += (t1 - t0);

    //- Remove electrons
    if (findIndex(typeIdList_, "e-") != -1)
    {
        removeElectrons();
        buildCellOccupancy();
    }

    //- Move the particles ballistically with their current velocities
    t0 = DSMC_WTIME();
    #ifdef _OPENMP
    if (omp_get_max_threads() > 1)
    {
        moveParallel(deltaTValue());
    }
    else
    #endif
    {
        Cloud<dsmcParcel>::move(td, deltaTValue());
    }
    t1 = DSMC_WTIME();
    timeMove_ += (t1 - t0);

    //- Update cell occupancy
    t0 = t1;
    #ifdef _OPENMP
    if (omp_get_max_threads() > 1)
    {
        // moveParallel already compacted parcelPtrs_, only rebuild cell index
        buildCellIndex();
    }
    else
    #endif
    {
        buildParcelPtrs();
        buildCellIndex();
    }
    // buildCellOccupancy only if legacy consumers need it
    if (collisionPartnerSelectionModel_->type() != "noTimeCounterOMP")
    {
        buildCellOccupancy();
    }
    t1 = DSMC_WTIME();
    timeBuildCellOccupancy_ += (t1 - t0);

    //- Add electrons back after the move function
    if (findIndex(typeIdList_, "e-") != -1)
    {
        addElectrons();
        buildCellOccupancy();
    }

    //- Radial weighting for non-Cartesian flows
    t0 = DSMC_WTIME();
    coordSystem().evolve();
    t1 = DSMC_WTIME();
    timeCoordSystem_ += (t1 - t0);

    t0 = t1;
    controllers_.controlBeforeCollisions();
    boundaries_.controlBeforeCollisions();
    t1 = DSMC_WTIME();
    timeControllersBeforeCollisions_ += (t1 - t0);

    // Validate cell index in debug mode
    if (debug)
    {
        validateCellIndex();
    }

    //- Calculate new velocities via stochastic collisions
    t0 = t1;
    collisions();
    t1 = DSMC_WTIME();
    timeCollisions_ += (t1 - t0);

    //- Reactions may have changed particle list, rebuild index
    if (reactions_.nReactions() != 0)
    {
        buildParcelPtrs();
        buildCellIndex();
    }

    t0 = DSMC_WTIME();
    controllers_.controlAfterCollisions();
    boundaries_.controlAfterCollisions();
    t1 = DSMC_WTIME();
    timeControllersAfterCollisions_ += (t1 - t0);

    #undef DSMC_WTIME
}


void Foam::dsmcCloud::evolve_fields()
{
    #ifdef _OPENMP
        #define DSMC_WTIME() omp_get_wtime()
    #else
        #define DSMC_WTIME() mesh_.time().elapsedClockTime()
    #endif

    scalar t0 = DSMC_WTIME();

    // Build sorted index once for fields (cache-friendly cell traversal)

    reactions_.outputData();

    fields_.calculateFields();

    fields_.writeFields();

    controllers_.calculateProps();
    controllers_.outputResults();

    boundaries_.calculateProps();
    boundaries_.outputResults();

    boundaryMeas_.outputResults();

    trackingInfo_.clean();
    boundaryMeas_.clean();
    cellMeas_.clean();

    scalar t1 = DSMC_WTIME();
    timeFields_ += (t1 - t0);

    #undef DSMC_WTIME
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
    buildParcelPtrs();
    buildCellIndex();
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

        label val = start + label(rndGen().sample01<scalar>()*(end - start + 1));

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
        *rndGen().GaussNormal<vector>();
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
        CTry = rndGen().GaussNormal<vector>()/sqrt(2.0);

        const scalar gammaTry = 1.0 + (q & CTry)*(0.4*(CTry & CTry) - 1.0)
            - (CTry & (tau & CTry));

        if (gammaTry >= A*rndGen().sample01<scalar>())
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

        CTry = rndGen().GaussNormal<vector>()/sqrt(2.0);

        const scalar gammaTry = 1.0 + 2.0*(D & CTry)
            + (qTra & CTry)*(0.4*(CTry & CTry) - 1.0)
            + (qRot & CTry)*(epsRot - epsRotAv)
            + (qVib & CTry)*(epsVib - epsVibAv)
            - (CTry & (tau & CTry));

        if (gammaTry >= A*rndGen().sample01<scalar>())
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
        ERot = -log(rndGen().sample01<scalar>())*physicoChemical::k.value()*temperature;
    }
    else
    {
        scalar a = 0.5*rotationalDof - 1;

        scalar energyRatio;

        scalar P = -1;

        do
        {
            energyRatio = 10*rndGen().sample01<scalar>();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < rndGen().sample01<scalar>());

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
            vibLevel[mode] = -log(rndGen().sample01<scalar>())*temperature
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
