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
#include "zeroGradientFvPatchFields.H"
#include <chrono>

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


void Foam::dsmcCloud::buildCellOccupancy()
{
    using clock_type = std::chrono::steady_clock;

    const auto t0 = clock_type::now();

    forAll(cellOccupancy_, celli)
    {
        cellOccupancy_[celli].clear();
    }

    const bool useMoveOrderedParcels =
        openmpEnabled_
     && openmpMoveEnabled_
     && moveOrderedParcelsValid_
     && moveOrderedParcels_.size() == this->size()
     && moveOrderedThreadOffsets_.size() == ompNumThreads_ + 1;

    const label nParcels =
        useMoveOrderedParcels ? moveOrderedParcels_.size() : this->size();
    const label nCells = cellOccupancy_.size();

    #ifdef _OPENMP
    if (openmpEnabled_ && ompNumThreads_ > 1 && nParcels > 0)
    {
        List<dsmcParcel*> parcels;
        labelList parcelThreadOffsets;

        if (useMoveOrderedParcels)
        {
            parcels = moveOrderedParcels_;
            parcelThreadOffsets = moveOrderedThreadOffsets_;
        }
        else
        {
            parcels.setSize(nParcels);
            label parcelI = 0;

            forAllIter(dsmcCloud, *this, iter)
            {
                parcels[parcelI++] = &iter();
            }
        }

        const auto t1 = clock_type::now();

        List<labelList> threadCellCounts(ompNumThreads_);

        forAll(threadCellCounts, threadI)
        {
            threadCellCounts[threadI].setSize(nCells, 0);
        }

        #pragma omp parallel
        {
            const label threadI = currentThreadId();
            labelList& localCounts = threadCellCounts[threadI];

            if (useMoveOrderedParcels)
            {
                for (label i = parcelThreadOffsets[threadI]; i < parcelThreadOffsets[threadI + 1]; ++i)
                {
                    const label celli = parcels[i]->cell();

                    if (celli >= 0 && celli < nCells)
                    {
                        ++localCounts[celli];
                    }
                }
            }
            else
            {
                #pragma omp for schedule(static)
                for (label i = 0; i < nParcels; ++i)
                {
                    const label celli = parcels[i]->cell();

                    if (celli >= 0 && celli < nCells)
                    {
                        ++localCounts[celli];
                    }
                }
            }
        }

        labelList totalCounts(nCells, 0);

        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            for (label celli = 0; celli < nCells; ++celli)
            {
                totalCounts[celli] += threadCellCounts[threadI][celli];
            }
        }

        const auto t2 = clock_type::now();

        for (label celli = 0; celli < nCells; ++celli)
        {
            cellOccupancy_[celli].setCapacity(totalCounts[celli]);
            cellOccupancy_[celli].setSize(totalCounts[celli]);
        }

        rebuildParticleLoadPartition();

        for (label celli = 0; celli < nCells; ++celli)
        {
            label offset = 0;

            for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
            {
                const label count = threadCellCounts[threadI][celli];
                threadCellCounts[threadI][celli] = offset;
                offset += count;
            }
        }

        #pragma omp parallel
        {
            const label threadI = currentThreadId();
            labelList& localOffsets = threadCellCounts[threadI];

            if (useMoveOrderedParcels)
            {
                for (label i = parcelThreadOffsets[threadI]; i < parcelThreadOffsets[threadI + 1]; ++i)
                {
                    dsmcParcel* pPtr = parcels[i];
                    const label celli = pPtr->cell();

                    if (celli >= 0 && celli < nCells)
                    {
                        cellOccupancy_[celli][localOffsets[celli]++] = pPtr;
                    }
                }
            }
            else
            {
                #pragma omp for schedule(static)
                for (label i = 0; i < nParcels; ++i)
                {
                    dsmcParcel* pPtr = parcels[i];
                    const label celli = pPtr->cell();

                    if (celli >= 0 && celli < nCells)
                    {
                        cellOccupancy_[celli][localOffsets[celli]++] = pPtr;
                    }
                }
            }
        }

        const auto t3 = clock_type::now();

        if (evolveProfileEnabled_)
        {
            buildOccupancyExtractWallTime_ += std::chrono::duration<scalar>(t1 - t0).count();
            buildOccupancyCountWallTime_ += std::chrono::duration<scalar>(t2 - t1).count();
            buildOccupancyAssembleWallTime_ += std::chrono::duration<scalar>(t3 - t2).count();
            ++buildOccupancyProfileCalls_;

            if (mesh_.time().writeTime() && Pstream::master())
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

                buildOccupancyExtractWallTime_ = 0.0;
                buildOccupancyCountWallTime_ = 0.0;
                buildOccupancyAssembleWallTime_ = 0.0;
                buildOccupancyProfileCalls_ = 0;
            }
        }

        return;
    }
    #endif

    forAllIter(dsmcCloud, *this, iter)
    {
        if (iter().cell() >= 0 && iter().cell() < cellOccupancy_.size())
        {
            cellOccupancy_[iter().cell()].append(&iter());
        }
    }

    rebuildParticleLoadPartition();

    const auto t1 = clock_type::now();

    if (evolveProfileEnabled_)
    {
        buildOccupancyExtractWallTime_ += std::chrono::duration<scalar>(t1 - t0).count();
        ++buildOccupancyProfileCalls_;

        if (mesh_.time().writeTime() && Pstream::master())
        {
            Info<< "BuildCellOccupancy profiling summary:" << nl
                << "    buildCellOccupancy calls      = " << buildOccupancyProfileCalls_ << nl
                << "    extract parcels [s]           = " << buildOccupancyExtractWallTime_ << nl
                << "    count/reduce [s]              = " << buildOccupancyCountWallTime_ << nl
                << "    allocate/fill [s]             = " << buildOccupancyAssembleWallTime_ << nl
                << "    total profiled [s]            = " << buildOccupancyExtractWallTime_ << nl
                << endl;

            buildOccupancyExtractWallTime_ = 0.0;
            buildOccupancyCountWallTime_ = 0.0;
            buildOccupancyAssembleWallTime_ = 0.0;
            buildOccupancyProfileCalls_ = 0;
        }
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
    if (openmpEnabled_ && openmpCollisionStrategy_ == "partition")
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

        if (mesh_.time().writeTime() && Pstream::master())
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

            collisionPrecomputeWallTime_ = 0.0;
            collisionPartitionWallTime_ = 0.0;
            collisionSelectionWallTime_ = 0.0;
            collisionProfileCalls_ = 0;
        }
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
    openmpMoveEnabled_ = controlDict.lookupOrDefault<bool>("openmpMove", false);
    ompNumThreads_ = controlDict.lookupOrDefault<label>("openmpThreads", 0);
    openmpCollisionStrategy_ =
        controlDict.lookupOrDefault<word>("openmpCollisionStrategy", "dynamic");
    openmpMoveSchedule_ =
        controlDict.lookupOrDefault<word>("openmpMoveSchedule", "static");
    openmpMoveChunk_ =
        controlDict.lookupOrDefault<label>("openmpMoveChunk", 64);
    collisionProfileEnabled_ =
        controlDict.lookupOrDefault<bool>("profileCollisionPhases", false);
    evolveProfileEnabled_ =
        controlDict.lookupOrDefault<bool>("profileEvolvePhases", false);

    #ifdef _OPENMP
    if (openmpEnabled_)
    {
        if
        (
            openmpCollisionStrategy_ != "dynamic"
         && openmpCollisionStrategy_ != "partition"
        )
        {
            WarningInFunction
                << "Unknown openmpCollisionStrategy '"
                << openmpCollisionStrategy_
                << "'. Falling back to 'dynamic'." << endl;

            openmpCollisionStrategy_ = "dynamic";
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

        forAll(ompRndGens_, threadI)
        {
            const label seed = 104729*Pstream::myProcNo() + threadI + 1;
            ompRndGens_[threadI].reset(seed);
        }

        Info<< "OpenMP enabled for dsmcCloud with "
            << ompNumThreads_ << " thread-local RNG streams"
            << " using collision strategy '" << openmpCollisionStrategy_
            << "', move kernel "
            << (openmpMoveEnabled_ ? "enabled" : "disabled")
            << " (" << openmpMoveSchedule_ << ", chunk "
            << openmpMoveChunk_ << ")"
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
    #endif
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

    #ifdef _OPENMP
    if (openmpEnabled_)
    {
        #pragma omp parallel for schedule(static)
        for (label celli = 0; celli < mesh_.nCells(); ++celli)
        {
            const label nC = cellOccupancy_[celli].size();
            const scalar sigmaTcRMaxCell = max(sigmaTcRMax_[celli], SMALL);
            const scalar selectedPairs =
                collisionSelectionRemainder_[celli]
              + 0.5*nC*(nC - 1)*nParticles(celli)*sigmaTcRMaxCell*deltaTValue(celli)
               /mesh_.cellVolumes()[celli];

            selectedPairsPerCell_[celli] = selectedPairs;
            nCandidatesPerCell_[celli] = label(selectedPairs);
            collisionSelectionRemainder_[celli] =
                selectedPairs - scalar(nCandidatesPerCell_[celli]);
        }
    }
    else
    #endif
    {
        forAll(selectedPairsPerCell_, celli)
        {
            const label nC = cellOccupancy_[celli].size();
            const scalar sigmaTcRMaxCell = max(sigmaTcRMax_[celli], SMALL);
            const scalar selectedPairs =
                collisionSelectionRemainder_[celli]
              + 0.5*nC*(nC - 1)*nParticles(celli)*sigmaTcRMaxCell*deltaTValue(celli)
               /mesh_.cellVolumes()[celli];

            selectedPairsPerCell_[celli] = selectedPairs;
            nCandidatesPerCell_[celli] = label(selectedPairs);
            collisionSelectionRemainder_[celli] =
                selectedPairs - scalar(nCandidatesPerCell_[celli]);
        }
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

    collisionLoadStart_.setSize(ompNumThreads_, mesh_.nCells());
    collisionLoadEnd_.setSize(ompNumThreads_, mesh_.nCells());

    label totalCandidates = 0;

    forAll(nCandidatesPerCell_, celli)
    {
        totalCandidates += nCandidatesPerCell_[celli];
    }

    if (totalCandidates <= 0)
    {
        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            collisionLoadStart_[threadI] = threadI*mesh_.nCells()/ompNumThreads_;
            collisionLoadEnd_[threadI] = (threadI + 1)*mesh_.nCells()/ompNumThreads_;
        }

        return;
    }

    const scalar avgCandidates = scalar(totalCandidates)/scalar(ompNumThreads_);
    scalar accumulatedCandidates = 0.0;
    label start = 0;

    for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
    {
        collisionLoadStart_[threadI] = start;

        if (threadI == ompNumThreads_ - 1)
        {
            collisionLoadEnd_[threadI] = mesh_.nCells();
            break;
        }

        const scalar targetCandidates = avgCandidates*scalar(threadI + 1);
        label end = start;

        while (end < mesh_.nCells() && accumulatedCandidates < targetCandidates)
        {
            accumulatedCandidates += nCandidatesPerCell_[end];
            ++end;
        }

        collisionLoadEnd_[threadI] = end;
        start = end;
    }
}


void Foam::dsmcCloud::rebuildParticleLoadPartition()
{
    if (!openmpEnabled_ || ompNumThreads_ <= 1)
    {
        particleLoadStart_.setSize(1, 0);
        particleLoadEnd_.setSize(1, mesh_.nCells());
        return;
    }

    particleLoadStart_.setSize(ompNumThreads_, mesh_.nCells());
    particleLoadEnd_.setSize(ompNumThreads_, mesh_.nCells());

    label totalParticles = 0;

    forAll(cellOccupancy_, celli)
    {
        totalParticles += cellOccupancy_[celli].size();
    }

    if (totalParticles <= 0)
    {
        for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
        {
            particleLoadStart_[threadI] = threadI*mesh_.nCells()/ompNumThreads_;
            particleLoadEnd_[threadI] = (threadI + 1)*mesh_.nCells()/ompNumThreads_;
        }

        return;
    }

    const scalar avgParticles = scalar(totalParticles)/scalar(ompNumThreads_);
    scalar accumulatedParticles = 0.0;
    label start = 0;

    for (label threadI = 0; threadI < ompNumThreads_; ++threadI)
    {
        particleLoadStart_[threadI] = start;

        if (threadI == ompNumThreads_ - 1)
        {
            particleLoadEnd_[threadI] = mesh_.nCells();
            break;
        }

        const scalar targetParticles = avgParticles*scalar(threadI + 1);
        label end = start;

        while (end < mesh_.nCells() && accumulatedParticles < targetParticles)
        {
            accumulatedParticles += cellOccupancy_[end].size();
            ++end;
        }

        particleLoadEnd_[threadI] = end;
        start = end;
    }
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
    openmpCollisionStrategy_("dynamic"),
    openmpMoveSchedule_("static"),
    openmpMoveChunk_(64),
    collisionProfileEnabled_(false),
    evolveProfileEnabled_(false),
    ompRndGens_(),
    particleLoadStart_(),
    particleLoadEnd_(),
    moveOrderedParcels_(),
    moveOrderedThreadOffsets_(),
    moveOrderedParcelsValid_(false),
    selectedPairsPerCell_(mesh_.nCells(), 0.0),
    nCandidatesPerCell_(mesh_.nCells(), 0),
    collisionPrecomputeWallTime_(0.0),
    collisionPartitionWallTime_(0.0),
    collisionSelectionWallTime_(0.0),
    collisionProfileCalls_(0),
    buildOccupancyExtractWallTime_(0.0),
    buildOccupancyCountWallTime_(0.0),
    buildOccupancyAssembleWallTime_(0.0),
    buildOccupancyProfileCalls_(0),
    evolveMoveWallTime_(0.0),
    evolveBuildWallTime_(0.0),
    evolveCoordWallTime_(0.0),
    evolveCollisionWallTime_(0.0),
    evolveReactionWallTime_(0.0),
    evolvePostWallTime_(0.0),
    evolveProfileCalls_(0),
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

    coordSystem().checkCoordinateSystemInputs();
    buildConstProps();

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

    buildCellOccupancy();
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
    this->addParticle(pPtr);
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
    this->addParticle(pPtr);
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
    this->addParticle(pPtr);
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

    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    dsmcParcel::trackingData td(*this);

    const auto t0 = clock_type::now();

    controllers_.controlBeforeMove();
    boundaries_.controlBeforeMove();

    if (openmpEnabled_ && openmpMoveEnabled_)
    {
        rebuildParticleLoadPartition();
    }

    Cloud<dsmcParcel>::move(*this, td, mesh_.time().deltaTValue());
    const auto t1 = clock_type::now();

    buildCellOccupancy();
    const auto t2 = clock_type::now();

    coordSystem().evolve();
    const auto t3 = clock_type::now();

    controllers_.controlBeforeCollisions();
    boundaries_.controlBeforeCollisions();

    collisions();
    const auto t4 = clock_type::now();

    if (reactionsActive())
    {
        buildCellOccupancy();
        reactions().outputData();
    }
    const auto t5 = clock_type::now();

    controllers_.controlAfterCollisions();
    boundaries_.controlAfterCollisions();

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

    const auto t6 = clock_type::now();

    if (evolveProfileEnabled_)
    {
        evolveMoveWallTime_ += std::chrono::duration<scalar>(t1 - t0).count();
        evolveBuildWallTime_ += std::chrono::duration<scalar>(t2 - t1).count();
        evolveCoordWallTime_ += std::chrono::duration<scalar>(t3 - t2).count();
        evolveCollisionWallTime_ += std::chrono::duration<scalar>(t4 - t3).count();
        evolveReactionWallTime_ += std::chrono::duration<scalar>(t5 - t4).count();
        evolvePostWallTime_ += std::chrono::duration<scalar>(t6 - t5).count();
        ++evolveProfileCalls_;

        if (mesh_.time().writeTime() && Pstream::master())
        {
            const scalar totalProfiled =
                evolveMoveWallTime_
              + evolveBuildWallTime_
              + evolveCoordWallTime_
              + evolveCollisionWallTime_
              + evolveReactionWallTime_
              + evolvePostWallTime_;

            Info<< "Evolve profiling summary:" << nl
                << "    evolve calls                  = " << evolveProfileCalls_ << nl
                << "    move only [s]                 = " << evolveMoveWallTime_ << nl
                << "    buildCellOccupancy [s]        = " << evolveBuildWallTime_ << nl
                << "    coordSystem [s]               = " << evolveCoordWallTime_ << nl
                << "    collision phase [s]           = " << evolveCollisionWallTime_ << nl
                << "    reaction/output [s]           = " << evolveReactionWallTime_ << nl
                << "    post fields/output [s]        = " << evolvePostWallTime_ << nl
                << "    total profiled [s]            = " << totalProfiled << nl
                << endl;

            evolveMoveWallTime_ = 0.0;
            evolveBuildWallTime_ = 0.0;
            evolveCoordWallTime_ = 0.0;
            evolveCollisionWallTime_ = 0.0;
            evolveReactionWallTime_ = 0.0;
            evolvePostWallTime_ = 0.0;
            evolveProfileCalls_ = 0;
        }
    }
}


void Foam::dsmcCloud::loadBalanceCheck()
{}


void Foam::dsmcCloud::info() const
{
    label nParcels = this->size();
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

    Info<< "Cloud name: " << this->name() << nl
        << "    Number of dsmc particles        = " << nParcels << nl;

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


void Foam::dsmcCloud::reportProfiling() const
{
    if (!Pstream::master())
    {
        return;
    }

    if (buildOccupancyProfileCalls_ > 0)
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
    }

    if (collisionProfileCalls_ > 0)
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

    if (evolveProfileCalls_ > 0)
    {
        const scalar totalProfiled =
            evolveMoveWallTime_
          + evolveBuildWallTime_
          + evolveCoordWallTime_
          + evolveCollisionWallTime_
          + evolveReactionWallTime_
          + evolvePostWallTime_;

        Info<< "Evolve profiling summary:" << nl
            << "    evolve calls                  = " << evolveProfileCalls_ << nl
            << "    move only [s]                 = " << evolveMoveWallTime_ << nl
            << "    buildCellOccupancy [s]        = " << evolveBuildWallTime_ << nl
            << "    coordSystem [s]               = " << evolveCoordWallTime_ << nl
            << "    collision phase [s]           = " << evolveCollisionWallTime_ << nl
            << "    reaction/output [s]           = " << evolveReactionWallTime_ << nl
            << "    post fields/output [s]        = " << evolvePostWallTime_ << nl
            << "    total profiled [s]            = " << totalProfiled << nl
            << endl;
    }
}

// ************************************************************************* //


















