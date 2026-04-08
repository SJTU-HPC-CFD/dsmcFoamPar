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

Description

Measures DSMC macroscopic fields for a single species or a gas mixture and
writes the results to a volume field that can be viewed in Paraview.

Translational, rotatational and vibrational temperature fields will also be
written automatically.

Boundary fields are measured in conjunction with the boundaryMeasurements class
and are also written.

\*---------------------------------------------------------------------------*/

#include "dsmcVolFields.H"
#include "addToRunTimeSelectionTable.H"
#include <chrono>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace
{

struct dsmcVolSharedSampleCache
{
    struct BuildProfile
    {
        scalar allocateWallTime = 0.0;
        scalar resetWallTime = 0.0;
        scalar parcelAccumWallTime = 0.0;
        scalar baseAccumWallTime = 0.0;
        scalar vibAccumWallTime = 0.0;
        scalar electronicAccumWallTime = 0.0;
        scalar classAccumWallTime = 0.0;
        label detailSampleCells = 0;
        label detailSampleParcels = 0;
    };

    const dsmcCloud* cloudPtr = nullptr;
    scalar timeValue = -GREAT;
    label nCells = 0;
    label nTypes = 0;
    bool built = false;

    List<scalarField> dsmcN;
    List<scalarField> dsmcM;
    List<scalarField> dsmcLinearKE;
    List<vectorField> dsmcMomentum;
    List<scalarField> dsmcErot;
    List<scalarField> dsmcZetaRot;
    List<scalarField> dsmcSpeciesEelec;
    List<scalarField> dsmcNElecLvl;
    List<scalarField> nGrndElecLvl;
    List<scalarField> n1stElecLvl;
    List<scalarField> nReal;
    List<scalarField> mReal;
    List<vectorField> momentumReal;
    List<scalarField> linearKEReal;
    List<scalarField> dsmcMuu;
    List<scalarField> dsmcMuv;
    List<scalarField> dsmcMuw;
    List<scalarField> dsmcMvv;
    List<scalarField> dsmcMvw;
    List<scalarField> dsmcMww;
    List<scalarField> dsmcMcc;
    List<scalarField> dsmcMccu;
    List<scalarField> dsmcMccv;
    List<scalarField> dsmcMccw;
    List<scalarField> dsmcEu;
    List<scalarField> dsmcEv;
    List<scalarField> dsmcEw;
    List<scalarField> dsmcECum;
    List<scalarField> dsmcNClassI;
    List<scalarField> dsmcNClassII;
    List<scalarField> dsmcNClassIII;
    List<List<scalarField>> dsmcSpeciesEvibMod;

    bool validFor(const dsmcCloud& cloud, const scalar currentTime) const
    {
        return cloudPtr == &cloud
            && mag(timeValue - currentTime) < SMALL
            && nCells == cloud.mesh().nCells()
            && nTypes == cloud.constProps().size();
    }

    void allocateFields(const dsmcCloud& cloud)
    {
        cloudPtr = &cloud;
        nCells = cloud.mesh().nCells();
        nTypes = cloud.constProps().size();
        built = false;

        auto initScalarFields =
            [this](List<scalarField>& fields)
            {
                fields.setSize(nTypes);
                for (label typei = 0; typei < nTypes; ++typei)
                {
                    fields[typei].setSize(nCells, 0.0);
                }
            };

        auto initVectorFields =
            [this](List<vectorField>& fields)
            {
                fields.setSize(nTypes);
                for (label typei = 0; typei < nTypes; ++typei)
                {
                    fields[typei].setSize(nCells, vector::zero);
                }
            };

        initScalarFields(dsmcN);
        initScalarFields(dsmcM);
        initScalarFields(dsmcLinearKE);
        initVectorFields(dsmcMomentum);
        initScalarFields(dsmcErot);
        initScalarFields(dsmcZetaRot);
        initScalarFields(dsmcSpeciesEelec);
        initScalarFields(dsmcNElecLvl);
        initScalarFields(nGrndElecLvl);
        initScalarFields(n1stElecLvl);
        initScalarFields(nReal);
        initScalarFields(mReal);
        initVectorFields(momentumReal);
        initScalarFields(linearKEReal);
        initScalarFields(dsmcMuu);
        initScalarFields(dsmcMuv);
        initScalarFields(dsmcMuw);
        initScalarFields(dsmcMvv);
        initScalarFields(dsmcMvw);
        initScalarFields(dsmcMww);
        initScalarFields(dsmcMcc);
        initScalarFields(dsmcMccu);
        initScalarFields(dsmcMccv);
        initScalarFields(dsmcMccw);
        initScalarFields(dsmcEu);
        initScalarFields(dsmcEv);
        initScalarFields(dsmcEw);
        initScalarFields(dsmcECum);
        initScalarFields(dsmcNClassI);
        initScalarFields(dsmcNClassII);
        initScalarFields(dsmcNClassIII);

        dsmcSpeciesEvibMod.setSize(nTypes);
        for (label typei = 0; typei < nTypes; ++typei)
        {
            const label nMods = cloud.constProps(typei).thetaV().size();
            dsmcSpeciesEvibMod[typei].setSize(nMods);

            for (label mod = 0; mod < nMods; ++mod)
            {
                dsmcSpeciesEvibMod[typei][mod].setSize(nCells, 0.0);
            }
        }
    }

    void build
    (
        const dsmcCloud& cloud,
        const List<DynamicList<dsmcParcel*>>& cellOccupancy,
        const List<dsmcParcel*>* occupancyOrderedParcelsPtr,
        const labelList* occupancyCellOffsetsPtr,
        const scalar currentTime,
        const bool needVibrational,
        const bool needElectronic,
        const bool needClassification,
        const bool needHeatFluxShearStress,
        BuildProfile* buildProfile = nullptr
    )
    {
        const scalar kBoltzmann = constant::physicoChemical::k.value();
        const bool doProfile = (buildProfile != nullptr);
        const bool useFlatOccupancy =
            occupancyOrderedParcelsPtr
         && occupancyCellOffsetsPtr
         && occupancyCellOffsetsPtr->size() == cellOccupancy.size() + 1
         && occupancyOrderedParcelsPtr->size() == occupancyCellOffsetsPtr->last();

        auto wallClockNow = []()
        {
            return std::chrono::steady_clock::now();
        };

        auto wallSeconds =
            [](const std::chrono::steady_clock::time_point& start,
               const std::chrono::steady_clock::time_point& end)
            {
                return std::chrono::duration<scalar>(end - start).count();
            };

        if (validFor(cloud, currentTime) && built)
        {
            return;
        }

        if (!validFor(cloud, currentTime))
        {
            const auto allocateStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            allocateFields(cloud);
            if (doProfile)
            {
                buildProfile->allocateWallTime +=
                    wallSeconds(allocateStart, wallClockNow());
            }
        }
        else
        {
            const auto resetStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            forAll(dsmcN, typei)
            {
                dsmcN[typei] = 0.0;
                dsmcM[typei] = 0.0;
                dsmcLinearKE[typei] = 0.0;
                dsmcMomentum[typei] = vector::zero;
                dsmcErot[typei] = 0.0;
                dsmcZetaRot[typei] = 0.0;
                dsmcSpeciesEelec[typei] = 0.0;
                dsmcNElecLvl[typei] = 0.0;
                nGrndElecLvl[typei] = 0.0;
                n1stElecLvl[typei] = 0.0;
                nReal[typei] = 0.0;
                mReal[typei] = 0.0;
                momentumReal[typei] = vector::zero;
                linearKEReal[typei] = 0.0;
                dsmcMuu[typei] = 0.0;
                dsmcMuv[typei] = 0.0;
                dsmcMuw[typei] = 0.0;
                dsmcMvv[typei] = 0.0;
                dsmcMvw[typei] = 0.0;
                dsmcMww[typei] = 0.0;
                dsmcMcc[typei] = 0.0;
                dsmcMccu[typei] = 0.0;
                dsmcMccv[typei] = 0.0;
                dsmcMccw[typei] = 0.0;
                dsmcEu[typei] = 0.0;
                dsmcEv[typei] = 0.0;
                dsmcEw[typei] = 0.0;
                dsmcECum[typei] = 0.0;
                dsmcNClassI[typei] = 0.0;
                dsmcNClassII[typei] = 0.0;
                dsmcNClassIII[typei] = 0.0;

                forAll(dsmcSpeciesEvibMod[typei], mod)
                {
                    dsmcSpeciesEvibMod[typei][mod] = 0.0;
                }
            }

            if (doProfile)
            {
                buildProfile->resetWallTime +=
                    wallSeconds(resetStart, wallClockNow());
            }
        }

        timeValue = currentTime;

        const auto parcelAccumStart =
            doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();

        auto accumulateCell =
            [&](const label celli)
            {
                scalar localBaseAccumWallTime = 0.0;
                scalar localVibAccumWallTime = 0.0;
                scalar localElectronicAccumWallTime = 0.0;
                scalar localClassAccumWallTime = 0.0;
                label localDetailSampleCells = 0;
                label localDetailSampleParcels = 0;
                const DynamicList<dsmcParcel*>& parcels = cellOccupancy[celli];
                const label parcelBegin =
                    useFlatOccupancy ? (*occupancyCellOffsetsPtr)[celli] : 0;
                const label parcelEnd =
                    useFlatOccupancy ? (*occupancyCellOffsetsPtr)[celli + 1] : parcels.size();
                const label parcelCount =
                    useFlatOccupancy ? (parcelEnd - parcelBegin) : parcels.size();
                const bool profileCellDetail =
                    doProfile && parcelCount > 0 && (celli % 32 == 0);

                if (profileCellDetail)
                {
                    localDetailSampleCells = 1;
                    localDetailSampleParcels = parcelCount;
                }

                for (label pi = 0; pi < parcelCount; ++pi)
                {
                    const dsmcParcel& p =
                        useFlatOccupancy
                      ? *(*occupancyOrderedParcelsPtr)[parcelBegin + pi]
                      : *parcels[pi];

                    if (!p.isFree())
                    {
                        continue;
                    }

                    const label typeId = p.typeId();
                    const dsmcParcel::constantProperties& cP = cloud.constProps(typeId);
                    const scalar nParticles = cloud.nParticles(celli);
                    const scalar mp = cP.mass();
                    const vector& Up = p.U();
                    const scalar linearKE = mp*(Up & Up);
                    const scalar Erotp = p.ERot();
                    const scalar zetaRotp = cP.rotationalDegreesOfFreedom();

                    scalar Evibp = 0.0;
                    std::chrono::steady_clock::time_point vibAccumStart;
                    if (profileCellDetail)
                    {
                        vibAccumStart = wallClockNow();
                    }
                    if (needVibrational)
                    {
                        const labelList& vibLevels = p.vibLevel();
                        const scalarList& thetaV = cP.thetaV();

                        forAll(thetaV, mod)
                        {
                            const scalar EvibMod =
                                kBoltzmann*thetaV[mod]*(vibLevels[mod] + 0.5);
                            dsmcSpeciesEvibMod[typeId][mod][celli] += EvibMod;
                            Evibp += EvibMod;
                        }
                    }
                    if (profileCellDetail)
                    {
                        localVibAccumWallTime +=
                            wallSeconds(vibAccumStart, wallClockNow());
                    }

                    const label nElecLevels = cP.nElectronicLevels();
                    const label eLevel = p.ELevel();

                    std::chrono::steady_clock::time_point baseAccumStart;
                    if (profileCellDetail)
                    {
                        baseAccumStart = wallClockNow();
                    }
                    dsmcN[typeId][celli] += 1.0;
                    dsmcM[typeId][celli] += mp;
                    dsmcLinearKE[typeId][celli] += linearKE;
                    dsmcMomentum[typeId][celli] += mp*Up;
                    dsmcErot[typeId][celli] += Erotp;
                    dsmcZetaRot[typeId][celli] += zetaRotp;

                    nReal[typeId][celli] += nParticles;
                    mReal[typeId][celli] += mp*nParticles;
                    momentumReal[typeId][celli] += mp*Up*nParticles;
                    linearKEReal[typeId][celli] += linearKE*nParticles;

                    if (needHeatFluxShearStress)
                    {
                        const scalar Eintp = Erotp + Evibp;

                        dsmcMuu[typeId][celli] += mp*sqr(Up.x());
                        dsmcMuv[typeId][celli] += mp*Up.x()*Up.y();
                        dsmcMuw[typeId][celli] += mp*Up.x()*Up.z();
                        dsmcMvv[typeId][celli] += mp*sqr(Up.y());
                        dsmcMvw[typeId][celli] += mp*Up.y()*Up.z();
                        dsmcMww[typeId][celli] += mp*sqr(Up.z());

                        dsmcMcc[typeId][celli] += linearKE;
                        dsmcMccu[typeId][celli] += linearKE*Up.x();
                        dsmcMccv[typeId][celli] += linearKE*Up.y();
                        dsmcMccw[typeId][celli] += linearKE*Up.z();

                        dsmcEu[typeId][celli] += Eintp*Up.x();
                        dsmcEv[typeId][celli] += Eintp*Up.y();
                        dsmcEw[typeId][celli] += Eintp*Up.z();
                        dsmcECum[typeId][celli] += Eintp;
                    }

                    if (needElectronic)
                    {
                        dsmcSpeciesEelec[typeId][celli] +=
                            cP.electronicEnergyList()[eLevel];
                    }
                    if (profileCellDetail)
                    {
                        localBaseAccumWallTime +=
                            wallSeconds(baseAccumStart, wallClockNow());
                    }

                    if (needElectronic && nElecLevels > 1)
                    {
                        std::chrono::steady_clock::time_point electronicAccumStart;
                        if (profileCellDetail)
                        {
                            electronicAccumStart = wallClockNow();
                        }
                        dsmcNElecLvl[typeId][celli] += 1.0;

                        if (eLevel == 0)
                        {
                            nGrndElecLvl[typeId][celli] += 1.0;
                        }
                        if (eLevel == 1)
                        {
                            n1stElecLvl[typeId][celli] += 1.0;
                        }
                        if (profileCellDetail)
                        {
                            localElectronicAccumWallTime +=
                                wallSeconds(electronicAccumStart, wallClockNow());
                        }
                    }

                    std::chrono::steady_clock::time_point classAccumStart;
                    if (profileCellDetail && needClassification)
                    {
                        classAccumStart = wallClockNow();
                    }

                    const label classification = p.classification();
                    if (needClassification)
                    {
                        if (classification == 0)
                        {
                            dsmcNClassI[typeId][celli] += 1.0;
                        }
                        else if (classification == 1)
                        {
                            dsmcNClassII[typeId][celli] += 1.0;
                        }
                        else if (classification == 2)
                        {
                            dsmcNClassIII[typeId][celli] += 1.0;
                        }
                    }
                    if (profileCellDetail && needClassification)
                    {
                        localClassAccumWallTime +=
                            wallSeconds(classAccumStart, wallClockNow());
                    }
                }

                if (doProfile)
                {
                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    buildProfile->baseAccumWallTime += localBaseAccumWallTime;

                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    buildProfile->vibAccumWallTime += localVibAccumWallTime;

                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    buildProfile->electronicAccumWallTime += localElectronicAccumWallTime;

                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    buildProfile->classAccumWallTime += localClassAccumWallTime;

                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    buildProfile->detailSampleCells += localDetailSampleCells;

                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                    buildProfile->detailSampleParcels += localDetailSampleParcels;
                }
            };

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) if (cloud.openmpEnabled())
        forAll(cellOccupancy, celli)
        {
            accumulateCell(celli);
        }
        #else
        forAll(cellOccupancy, celli)
        {
            accumulateCell(celli);
        }
        #endif

        if (doProfile)
        {
            buildProfile->parcelAccumWallTime +=
                wallSeconds(parcelAccumStart, wallClockNow());
        }

        built = true;
    }
};

static dsmcVolSharedSampleCache sharedSampleCache_;

}

defineTypeNameAndDebug(dsmcVolFields, 0);

addToRunTimeSelectionTable(dsmcField, dsmcVolFields, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void dsmcVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];

        if (isA<wallPolyPatch>(pPatch))
        {
            const vectorField& fC = pPatch.faceCentres();

            forAll(n_[patchi], facei)
            {
                n_[patchi][facei] = pPatch.faceAreas()[facei]
                    /mag(pPatch.faceAreas()[facei]);

                //- Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                t1_[patchi][facei] = fC[facei]
                    - mesh_.points()[mesh_.faces()[pPatch.start() + facei][0]];
                t1_[patchi][facei] /= mag(t1_[patchi][facei]);

                //- Other tangential unit vector.  Rescaling in case face is not
                //  flat and n and t1 aren't perfectly orthogonal
                t2_[patchi][facei] = n_[patchi][facei]^t1_[patchi][facei];
                t2_[patchi][facei] /= mag(t2_[patchi][facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcVolFields::dsmcVolFields
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    sampleInterval_(1),
    sampleCounter_(0),
    nTimeSteps_(0.0),
    isSample_(false),
    mfpTref_(273.0),
    fieldName_(Foam::hyCompat::lookup(Foam::hyCompat::lookup(propsDict_, "fieldName"))),
    speciesIds_(),
    dsmcN_
    (
        IOobject
        (
            "dsmcN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    dsmcNMean_
    (
        IOobject
        (
            "dsmcNMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimMass/dimVolume, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimPressure, 0.0)
    ),
    Ttra_
    (
        IOobject
        (
            "Ttra_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    Trot_
    (
        IOobject
        (
            "Trot_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    Tvib_
    (
        IOobject
        (
            "Tvib_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    Telec_
    (
        IOobject
        (
            "Telec_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    Tov_
    (
        IOobject
        (
            "Tov_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    q_
    (
        IOobject
        (
            "wallHeatFlux_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    tau_
    (
        IOobject
        (
            "wallShearStress_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimPressure, 0.0)
    ),
    mfp_
    (
        IOobject
        (
            "mfp_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    mfpToDx_
    (
        IOobject
        (
            "mfpToDx_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    DxToMfp_
    (
        IOobject
        (
            "DxToMfp_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    meanCollisionRate_
    (
        IOobject
        (
            "meanCollisionRate_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    measuredCollisionRate_
    (
        IOobject
        (
            "measuredCollisionRate",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    meanCollisionTime_
    (
        IOobject
        (
            "mct_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, 1, 0, 0), 0.0)
    ),
    mctToDt_
    (
        IOobject
        (
            "mctToDt_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimless, 0.0)
    ),
    meanCollisionSeparation_
    (
        IOobject
        (
            "mcs_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    SOF_
    (
        IOobject
        (
            "SOFP_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    Ma_
    (
        IOobject
        (
            "Ma_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIDistribution_
    (
        IOobject
        (
            "classIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIIDistribution_
    (
        IOobject
        (
            "classIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIIIDistribution_
    (
        IOobject
        (
            "classIIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    densityError_
    (
        IOobject
        (
            "rhoMError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    velocityError_
    (
        IOobject
        (
            "UError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    temperatureError_
    (
        IOobject
        (
            "TError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    pressureError_
    (
        IOobject
        (
            "pError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    UMean_
    (
        IOobject
        (
            "U_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero)
    ),
    fD_
    (
        IOobject
        (
            "fD_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVector_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 0, -3, 0, 0),
            vector::zero
        )
    ),
    pressureTensor_
    (
        IOobject
        (
            "pressureTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimPressure,
            tensor::zero
        )
    ),
    shearStressTensor_
    (
        IOobject
        (
            "shearStressTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimPressure,
            tensor::zero
        )
    ),
    Cq_
    (
        IOobject
        (
            "Cq_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0",  dimless, 0.0)
    ),
    Cp_
    (
        IOobject
        (
            "Cp_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    totalWallForce_(vector::zero),
    freeStreamVelocity_(vector::zero),
    freeStreamVelocityValue_(0.0),
    freeStreamDensity_(0.0),
    cqCoef_(0.0),
    cpCoef_(0.0),
    dsmcNCum_(mesh_.nCells(), 0.0),
    nCum_(mesh_.nCells(), 0.0),
    dsmcNElecLvlCum_(mesh_.nCells(), 0.0),
    dsmcMCum_(mesh_.nCells(), 0.0),
    mCum_(mesh_.nCells(), 0.0),
    dsmcLinearKECum_(mesh_.nCells(), 0.0),
    linearKECum_(mesh_.nCells(), 0.0),
    dsmcErotCum_(mesh_.nCells(), 0.0),
    dsmcZetaRotCum_(mesh_.nCells(), 0.0),
    dsmcMuuCum_(mesh_.nCells(), 0.0),
    dsmcMuvCum_(mesh_.nCells(), 0.0),
    dsmcMuwCum_(mesh_.nCells(), 0.0),
    dsmcMvvCum_(mesh_.nCells(), 0.0),
    dsmcMvwCum_(mesh_.nCells(), 0.0),
    dsmcMwwCum_(mesh_.nCells(), 0.0),
    dsmcMccCum_(mesh_.nCells(), 0.0),
    dsmcMccuCum_(mesh_.nCells(), 0.0),
    dsmcMccvCum_(mesh_.nCells(), 0.0),
    dsmcMccwCum_(mesh_.nCells(), 0.0),
    dsmcEuCum_(mesh_.nCells(), 0.0),
    dsmcEvCum_(mesh_.nCells(), 0.0),
    dsmcEwCum_(mesh_.nCells(), 0.0),
    dsmcECum_(mesh_.nCells(), 0.0),
    zetaVib_(mesh_.nCells(), 0.0),
    dsmcNClassICum_(mesh_.nCells(), 0.0),
    dsmcNClassIICum_(mesh_.nCells(), 0.0),
    dsmcNClassIIICum_(mesh_.nCells(), 0.0),
    collisionSeparation_(mesh_.nCells(), 0.0),
    dsmcNCollsCum_(mesh_.nCells(), 0.0),
    dsmcMomentumCum_(mesh.nCells(), vector::zero),
    momentumCum_(mesh.nCells(), vector::zero),
    boundaryCells_(),
    dsmcSpeciesEvibModCum_(),
    dsmcSpeciesEelecCum_(),
    dsmcNSpeciesCum_(),
    nSpeciesCum_(),
    dsmcMccSpeciesCum_(),
    speciesTvib_(),
    dsmcNGrndElecLvlSpeciesCum_(),
    dsmcN1stElecLvlSpeciesCum_(),
    speciesMfp_(),
    speciesMcr_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    ErotBF_(),
    zetaRotBF_(),
    qBF_(),
    zetaVibBF_(),
    rhoNIntBF_(),
    rhoNElecBF_(),
    momentumBF_(),
    fDBF_(),
    speciesEvibBF_(),
    speciesEelecBF_(),
    speciesRhoNBF_(),
    speciesMccBF_(),
    speciesTvibBF_(),
    speciesZetaVibBF_(),
    speciesEvibModBF_(),
    n_(),
    t1_(),
    t2_(),
    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false),
    densityOnly_(false),
    measureHeatFluxShearStress_(false),
    writeRotationalTemperature_(false),
    writeVibrationalTemperature_(false),
    writeElectronicTemperature_(false),
    profileSampleAccumWallTime_(0.0),
    profileSharedCacheBuildWallTime_(0.0),
    profileSharedCacheAllocateWallTime_(0.0),
    profileSharedCacheResetWallTime_(0.0),
    profileSharedCacheParcelAccumWallTime_(0.0),
    profileSharedCacheBaseAccumWallTime_(0.0),
    profileSharedCacheVibAccumWallTime_(0.0),
    profileSharedCacheElectronicAccumWallTime_(0.0),
    profileSharedCacheClassAccumWallTime_(0.0),
    profileSharedCacheDetailSampleCells_(0),
    profileSharedCacheDetailSampleParcels_(0),
    profileFieldCombineWallTime_(0.0),
    profileCellReduceWallTime_(0.0),
    profileBoundaryAccumWallTime_(0.0),
    profileOutputComputeWallTime_(0.0),
    profileFieldWriteWallTime_(0.0),
    profileOutputResetWallTime_(0.0),
    profileOutputTimeWallTime_(0.0),
    profileCalls_(0),
    finalProfilePrinted_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVolFields::~dsmcVolFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcVolFields::readIn()
{
    IOdictionary resumeSampling
    (
        IOobject
        (
            "resumeSampling_" + fieldName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );
    
    if (resumeSampling.size() > 0)
    {
        dictionary dict
        (
            resumeSampling.readStream(resumeSampling.filePath())
        );

        dict.readIfPresent("nTimeSteps", nTimeSteps_);

        // DSMC parcel related cumulative values
        dict.readIfPresent("dsmcNCum", dsmcNCum_);
        dict.readIfPresent("dsmcMCum", dsmcMCum_);

        dict.readIfPresent("dsmcLinearKECum", dsmcLinearKECum_);
        dict.readIfPresent("dsmcMomentumCum", dsmcMomentumCum_);
        dict.readIfPresent("dsmcErotCum", dsmcErotCum_);
        dict.readIfPresent("dsmcZetaRotCum", dsmcZetaRotCum_);
        dict.readIfPresent("dsmcSpeciesEelecCum", dsmcSpeciesEelecCum_);
        dict.readIfPresent("dsmcNSpeciesCum", dsmcNSpeciesCum_);
        dict.readIfPresent("dsmcMccSpeciesCum", dsmcMccSpeciesCum_);
        dict.readIfPresent("dsmcMuuCum", dsmcMuuCum_);
        dict.readIfPresent("dsmcMuvCum", dsmcMuvCum_);
        dict.readIfPresent("dsmcMuwCum", dsmcMuwCum_);
        dict.readIfPresent("dsmcMvvCum", dsmcMvvCum_);
        dict.readIfPresent("dsmcMvwCum", dsmcMvwCum_);
        dict.readIfPresent("dsmcMwwCum", dsmcMwwCum_);
        dict.readIfPresent("dsmcMccCum", dsmcMccCum_);
        dict.readIfPresent("dsmcMccuCum", dsmcMccuCum_);
        dict.readIfPresent("dsmcMccvCum", dsmcMccvCum_);
        dict.readIfPresent("dsmcMccwCum", dsmcMccwCum_);
        dict.readIfPresent("dsmcEuCum", dsmcEuCum_);
        dict.readIfPresent("dsmcEvCum", dsmcEvCum_);
        dict.readIfPresent("dsmcEwCum", dsmcEwCum_);
        dict.readIfPresent("dsmcECum", dsmcECum_);
        dict.readIfPresent("dsmcNElecLvlCum", dsmcNElecLvlCum_);
        dict.readIfPresent
        (
            "dsmcNGrndElecLvlSpeciesCum",
            dsmcNGrndElecLvlSpeciesCum_
        );
        dict.readIfPresent
        (
            "dsmcN1stElecLvlSpeciesCum",
            dsmcN1stElecLvlSpeciesCum_
        );
        if (measureClassifications_)
        {
          dict.readIfPresent("dsmcNClassICum", dsmcNClassICum_);
          dict.readIfPresent("dsmcNClassIICum", dsmcNClassIICum_);
          dict.readIfPresent("dsmcNClassIIICum", dsmcNClassIIICum_);
        }
        dict.readIfPresent("dsmcSpeciesEvibModCum", dsmcSpeciesEvibModCum_);
        dict.readIfPresent("dsmcNCollsCum", dsmcNCollsCum_);
        
        // boundary measurements
        dict.readIfPresent("rhoNBF", rhoNBF_);
        dict.readIfPresent("rhoMBF", rhoMBF_);
        dict.readIfPresent("linearKEBF", linearKEBF_);
        dict.readIfPresent("momentumBF", momentumBF_);
        dict.readIfPresent("ErotBF", ErotBF_);
        dict.readIfPresent("zetaRotBF", zetaRotBF_);
        dict.readIfPresent("rhoNIntBF", rhoNIntBF_);
        dict.readIfPresent("rhoNElecBF", rhoNElecBF_);
        dict.readIfPresent("qBF", qBF_);
        dict.readIfPresent("fDBF", fDBF_);
        dict.readIfPresent("speciesRhoNBF", speciesRhoNBF_);
        dict.readIfPresent("speciesEvibBF", speciesEvibBF_);
        dict.readIfPresent("speciesEelecBF",  speciesEelecBF_);
        dict.readIfPresent("speciesMccBF", speciesMccBF_);
        dict.readIfPresent("speciesEvibModBF", speciesEvibModBF_);

        // cumulative values
        dict.readIfPresent("nCum", nCum_);
        dict.readIfPresent("mCum", mCum_);
        dict.readIfPresent("nSpeciesCum", nSpeciesCum_);
        dict.readIfPresent("momentumCum", momentumCum_);
        dict.readIfPresent("linearKECum", linearKECum_);
        dict.readIfPresent("collisionSeparation", collisionSeparation_);
    }
}


void dsmcVolFields::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "resumeSampling_" + fieldName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("nTimeSteps", nTimeSteps_);

        // DSMC parcel related cumulative values
        dict.add("dsmcNCum", dsmcNCum_);
        dict.add("dsmcMCum", dsmcMCum_);

        dict.add("dsmcLinearKECum", dsmcLinearKECum_);
        dict.add("dsmcMomentumCum", dsmcMomentumCum_);
        dict.add("dsmcErotCum", dsmcErotCum_);
        dict.add("dsmcZetaRotCum", dsmcZetaRotCum_);
        dict.add("dsmcSpeciesEelecCum", dsmcSpeciesEelecCum_);
        dict.add("dsmcNSpeciesCum", dsmcNSpeciesCum_);
        dict.add("dsmcMccSpeciesCum", dsmcMccSpeciesCum_);
        dict.add("dsmcMuuCum", dsmcMuuCum_);
        dict.add("dsmcMuvCum", dsmcMuvCum_);
        dict.add("dsmcMuwCum", dsmcMuwCum_);
        dict.add("dsmcMvvCum", dsmcMvvCum_);
        dict.add("dsmcMvwCum", dsmcMvwCum_);
        dict.add("dsmcMwwCum", dsmcMwwCum_);
        dict.add("dsmcMccCum", dsmcMccCum_);
        dict.add("dsmcMccuCum", dsmcMccuCum_);
        dict.add("dsmcMccvCum", dsmcMccvCum_);
        dict.add("dsmcMccwCum", dsmcMccwCum_);
        dict.add("dsmcEuCum", dsmcEuCum_);
        dict.add("dsmcEvCum", dsmcEvCum_);
        dict.add("dsmcEwCum", dsmcEwCum_);
        dict.add("dsmcECum", dsmcECum_);
        dict.add("dsmcNElecLvlCum", dsmcNElecLvlCum_);
        dict.add("dsmcNGrndElecLvlSpeciesCum", dsmcNGrndElecLvlSpeciesCum_);
        dict.add("dsmcN1stElecLvlSpeciesCum", dsmcN1stElecLvlSpeciesCum_);
        if (measureClassifications_)
        {
            dict.add("dsmcNClassICum", dsmcNClassICum_);
            dict.add("dsmcNClassIICum", dsmcNClassIICum_);
            dict.add("dsmcNClassIIICum", dsmcNClassIIICum_);
        }
        dict.add("dsmcSpeciesEvibModCum", dsmcSpeciesEvibModCum_);
        dict.add("dsmcNCollsCum", dsmcNCollsCum_);

        // cumulative values
        dict.add("nCum", nCum_);
        dict.add("mCum", mCum_);
        dict.add("nSpeciesCum", nSpeciesCum_);
        dict.add("momentumCum", momentumCum_);
        dict.add("linearKECum", linearKECum_);
        dict.add("collisionSeparation", collisionSeparation_);
        
        // boundary measurements
        dict.add("rhoNBF", rhoNBF_);
        dict.add("rhoMBF", rhoMBF_);
        dict.add("linearKEBF", linearKEBF_);
        dict.add("momentumBF", momentumBF_);
        dict.add("ErotBF", ErotBF_);
        dict.add("zetaRotBF", zetaRotBF_);
        dict.add("rhoNIntBF", rhoNIntBF_);
        dict.add("rhoNElecBF", rhoNElecBF_);
        dict.add("qBF", qBF_);
        dict.add("fDBF", fDBF_);
        dict.add("speciesRhoNBF", speciesRhoNBF_);
        dict.add("speciesEvibBF", speciesEvibBF_);
        dict.add("speciesEelecBF",  speciesEelecBF_);
        dict.add("speciesMccBF", speciesMccBF_);
        dict.add("speciesEvibModBF", speciesEvibModBF_);

        IOstream::streamFormat fmt = time_.time().writeFormat();
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();

        dict.regIOobject::writeObject(fmt, ver, cmp, true);
    }
}


//- Initial configuration
void dsmcVolFields::createField()
{
    Info << "Initialising dsmcVolFields field" << endl;

    const List<word>& species(Foam::hyCompat::lookup(Foam::hyCompat::lookup(propsDict_, "typeIds")));

    DynamicList<word> speciesReduced(0);

    forAll(species, i)
    {
        const word& speciesName(species[i]);

        if (Foam::hyCompat::indexOf(speciesReduced, speciesName) == -1)
        {
            speciesReduced.append(speciesName);
        }
    }

    speciesReduced.shrink();

    speciesIds_.setSize(speciesReduced.size(), -1);

    forAll(speciesReduced, i)
    {
        const word& speciesName = speciesReduced[i];

        const label spId = Foam::hyCompat::indexOf(cloud_.typeIdList(), speciesName);

        if (spId == -1)
        {
            FatalErrorIn("dsmcVolFields::dsmcVolFields()")
                << "Cannot find typeId: " << speciesName << nl << "in: "
                << mesh_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        speciesIds_[i] = spId;
    }
    
    const label nCells = mesh_.nCells();
    const label nSpecies = speciesIds_.size();
    const label nPatches = mesh_.boundaryMesh().size();

    //- Volume fields initialisation
    dsmcNSpeciesCum_.setSize(nSpecies);
    nSpeciesCum_.setSize(nSpecies);
    dsmcMccSpeciesCum_.setSize(nSpecies);
    speciesMfp_.setSize(nSpecies);
    speciesMcr_.setSize(nSpecies);
    speciesTvib_.setSize(nSpecies);
    dsmcSpeciesEvibModCum_.setSize(nSpecies);
    dsmcSpeciesEelecCum_.setSize(nSpecies);
    dsmcNGrndElecLvlSpeciesCum_.setSize(nSpecies);
    dsmcN1stElecLvlSpeciesCum_.setSize(nSpecies);
    
    forAll(speciesIds_, i)
    {
        const label spId = speciesIds_[i];
        const label zetaRot =
            cloud_.constProps(spId).rotationalDegreesOfFreedom();
        const label nVibMod = cloud_.constProps(spId).nVibrationalModes();
        const label nElecLevels = cloud_.constProps(spId).nElectronicLevels();
        
        if (zetaRot > 0 and (not writeRotationalTemperature_))
        {
            writeRotationalTemperature_ = true;
        }
        
        if (nVibMod > 0 and (not writeVibrationalTemperature_))
        {
            writeVibrationalTemperature_ = true;
        }
        
        if (nElecLevels > 1 and (not writeElectronicTemperature_))
        {
            writeElectronicTemperature_ = true;
        }
        
        dsmcNSpeciesCum_[i].setSize(nCells, 0.0);
        nSpeciesCum_[i].setSize(nCells, 0.0);
        dsmcMccSpeciesCum_[i].setSize(nCells, 0.0);
        speciesMfp_[i].setSize(nCells, 0.0);
        speciesMcr_[i].setSize(nCells, 0.0);
        speciesTvib_[i].setSize(nCells);
        dsmcSpeciesEvibModCum_[i].setSize(nVibMod);

        forAll(dsmcSpeciesEvibModCum_[i], j)
        {
            dsmcSpeciesEvibModCum_[i][j].setSize(nCells, 0.0);
        }
        
        dsmcSpeciesEelecCum_[i].setSize(nCells, 0.0);
        dsmcNGrndElecLvlSpeciesCum_[i].setSize(nCells, 0.0);
        dsmcN1stElecLvlSpeciesCum_[i].setSize(nCells, 0.0);
    }

    //- Boundary fields initialisation
    boundaryCells_.setSize(nPatches);
    rhoNBF_.setSize(nPatches);
    rhoMBF_.setSize(nPatches);
    linearKEBF_.setSize(nPatches);
    momentumBF_.setSize(nPatches);
    ErotBF_.setSize(nPatches);
    zetaRotBF_.setSize(nPatches);
    qBF_.setSize(nPatches);
    fDBF_.setSize(nPatches);
    zetaVibBF_.setSize(nPatches);
    rhoNIntBF_.setSize(nPatches);
    rhoNElecBF_.setSize(nPatches);

    n_.setSize(nPatches);
    t1_.setSize(nPatches);
    t2_.setSize(nPatches);
    
    forAll(boundaryCells_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        const label nFaces = patch.size();

        boundaryCells_[j].setSize(nFaces);
        rhoNBF_[j].setSize(nFaces, 0.0);
        rhoMBF_[j].setSize(nFaces, 0.0);
        linearKEBF_[j].setSize(nFaces, 0.0);
        momentumBF_[j].setSize(nFaces, vector::zero);
        ErotBF_[j].setSize(nFaces, 0.0);
        zetaRotBF_[j].setSize(nFaces, 0.0);
        qBF_[j].setSize(nFaces, 0.0);
        fDBF_[j].setSize(nFaces, vector::zero);
        zetaVibBF_[j].setSize(nFaces, 0.0);
        rhoNIntBF_[j].setSize(nFaces, 0.0);
        rhoNElecBF_[j].setSize(nFaces, 0.0);

        n_[j].setSize(nFaces, vector::zero);
        t1_[j].setSize(nFaces, vector::zero);
        t2_[j].setSize(nFaces, vector::zero);

        forAll(boundaryCells_[j], k)
        {
            boundaryCells_[j][k] = patch.faceCells()[k];
        }
    }

    calculateWallUnitVectors();

    speciesRhoNBF_.setSize(nSpecies);
    speciesMccBF_.setSize(nSpecies);
    speciesEvibBF_.setSize(nSpecies);
    speciesTvibBF_.setSize(nSpecies);
    speciesZetaVibBF_.setSize(nSpecies);
    speciesEvibModBF_.setSize(nSpecies);
    speciesEelecBF_.setSize(nSpecies);

    forAll(speciesIds_, i)
    {
        const label spId = speciesIds_[i];
        
        speciesRhoNBF_[i].setSize(nPatches);
        speciesMccBF_[i].setSize(nPatches);
        
        speciesEvibBF_[i].setSize(nPatches);
        speciesTvibBF_[i].setSize(nPatches);
        speciesZetaVibBF_[i].setSize(nPatches);
        speciesEvibModBF_[i].setSize
        (
            cloud_.constProps(spId).nVibrationalModes()
        );
        speciesEelecBF_[i].setSize(nPatches);

        forAll(speciesEvibBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            const label nFaces = patch.size();

            speciesRhoNBF_[i][j].setSize(nFaces, 0.0);
            speciesMccBF_[i][j].setSize(nFaces, 0.0);
            speciesEvibBF_[i][j].setSize(nFaces, 0.0);
            speciesTvibBF_[i][j].setSize(nFaces, 0.0);
            speciesZetaVibBF_[i][j].setSize(nFaces, 0.0);
            speciesEelecBF_[i][j].setSize(nFaces, 0.0);
        }

        forAll(speciesEvibModBF_[i], mod)
        {
            speciesEvibModBF_[i][mod].setSize(nPatches);
            forAll(speciesEvibModBF_[i][mod], j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                const label nFaces = patch.size();
                speciesEvibModBF_[i][mod][j].setSize(nFaces, 0.0);
            }
        }
    }

    sampleInterval_ = propsDict_.lookupOrDefault("sampleInterval", 1);

    measureClassifications_ =
        propsDict_.lookupOrDefault<bool>("measureClassifications", false);

    measureErrors_ = propsDict_.lookupOrDefault<bool>("measureErrors", false);

    densityOnly_ = propsDict_.lookupOrDefault<bool>("densityOnly", false);

    measureHeatFluxShearStress_ =
        propsDict_.lookupOrDefault<bool>("measureHeatFluxShearStress", false);

    measureMeanFreePath_ =
        propsDict_.lookupOrDefault<bool>("measureMeanFreePath", false);

    mfpTref_ =
        propsDict_.lookupOrDefault<scalar>("mfpReferenceTemperature", 273.0);

    averagingAcrossManyRuns_ =
        propsDict_.lookupOrDefault<bool>("averagingAcrossManyRuns", false);

    //- read in stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        if (!time_.resetFieldsAtOutput())
        {
            Info<< "Averaging across many runs for field " << fieldName_
                << " is enabled. Sampled data will be read from file."
                << endl;
            readIn();
        }
        else
        {
            Info<< "Averaging across many runs for field " << fieldName_
                << " will be enabled as soon as resetAtOutput is turned off."
                << endl;
        }
    }
}


void dsmcVolFields::calculateField()
{
    sampleCounter_++;

    const scalar kB = physicoChemical::k.value();
    const scalar NAvo = physicoChemical::NA.value();
    constexpr bool doProfile = false;

    auto wallClockNow = []()
    {
        return std::chrono::steady_clock::now();
    };

    auto wallSeconds =
        [](const std::chrono::steady_clock::time_point& start,
           const std::chrono::steady_clock::time_point& end)
        {
            return std::chrono::duration<scalar>(end - start).count();
        };
    
    //- Reset instantaneous number of DSMC parcels
    dsmcN_ = 0.0;
    if (sampleInterval_ <= sampleCounter_)
    {
        nTimeSteps_ += 1.0;
        const scalar nAvTimeSteps = nTimeSteps_;
        const auto& cellOccupancy = cloud_.cellOccupancy();
        const bool useOpenMPSampling = cloud_.openmpEnabled();
        const auto sampleAccumStart =
            doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();

        if (densityOnly_)
        {
            if (useOpenMPSampling)
            {
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static)
                #endif
                forAll(cellOccupancy, celli)
                {
                    const DynamicList<dsmcParcel*>& parcels = cellOccupancy[celli];

                    scalar dsmcNLocal = 0.0;
                    scalar nLocal = 0.0;
                    scalar mLocal = 0.0;

                    forAll(parcels, pi)
                    {
                        const dsmcParcel& p = *parcels[pi];
                        const label spId = Foam::hyCompat::indexOf(speciesIds_, p.typeId());

                        if (spId != -1 && p.isFree())
                        {
                            const scalar nParticles = cloud_.nParticles(celli);
                            const scalar mass = cloud_.constProps(p.typeId()).mass();

                            dsmcNLocal += 1.0;
                            nLocal += nParticles;
                            mLocal += mass*nParticles;
                        }
                    }

                    dsmcNCum_[celli] += dsmcNLocal;
                    dsmcN_[celli] += dsmcNLocal;
                    nCum_[celli] += nLocal;
                    mCum_[celli] += mLocal;
                }
            }
            else
            {
                //- Loop over the the entire parcel cloud
                forAllConstIter(dsmcCloud, cloud_, iter)
                {
                    const dsmcParcel& p = iter();
                    const label spId = Foam::hyCompat::indexOf(speciesIds_, p.typeId());

                    //- Do not consider adsorbed parcels
                    if (spId != -1 && p.isFree())
                    {
                        const label cell = p.cell();
                        const scalar nParticles = cloud_.nParticles(cell);
                        const scalar mass = cloud_.constProps(p.typeId()).mass();

                        // cumulative number of DSMC parcels
                        dsmcNCum_[cell] += 1.0;
                        // instantaneous number of DSMC parcels in this time step
                        dsmcN_[cell] += 1.0;
                        // cumulative number of real particles
                        nCum_[cell] += nParticles;
                        // cumulative mass of real particles
                        mCum_[cell] += mass*nParticles;
                    }
                }
            }

            if (doProfile)
            {
                profileSampleAccumWallTime_ +=
                    wallSeconds(sampleAccumStart, wallClockNow());
            }
        }
        else
        {
            const bool needVibrational = writeVibrationalTemperature_;
            const bool needElectronic = writeElectronicTemperature_;
            const bool needClassification = measureClassifications_;
            const bool needHeatFluxShearStress = measureHeatFluxShearStress_;
            const auto sharedCacheBuildStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            dsmcVolSharedSampleCache::BuildProfile sharedCacheBuildProfile;
            sharedSampleCache_.build
            (
                cloud_,
                cellOccupancy,
                cloud_.hasOccupancyOrderedParcels()
                  ? &cloud_.occupancyOrderedParcels()
                  : nullptr,
                cloud_.hasOccupancyOrderedParcels()
                  ? &cloud_.occupancyCellOffsets()
                  : nullptr,
                time_.time().value(),
                needVibrational,
                needElectronic,
                needClassification,
                needHeatFluxShearStress,
                doProfile ? &sharedCacheBuildProfile : nullptr
            );

            if (doProfile)
            {
                profileSharedCacheBuildWallTime_ +=
                    wallSeconds(sharedCacheBuildStart, wallClockNow());
                profileSharedCacheAllocateWallTime_ +=
                    sharedCacheBuildProfile.allocateWallTime;
                profileSharedCacheResetWallTime_ +=
                    sharedCacheBuildProfile.resetWallTime;
                profileSharedCacheParcelAccumWallTime_ +=
                    sharedCacheBuildProfile.parcelAccumWallTime;
                profileSharedCacheBaseAccumWallTime_ +=
                    sharedCacheBuildProfile.baseAccumWallTime;
                profileSharedCacheVibAccumWallTime_ +=
                    sharedCacheBuildProfile.vibAccumWallTime;
                profileSharedCacheElectronicAccumWallTime_ +=
                    sharedCacheBuildProfile.electronicAccumWallTime;
                profileSharedCacheClassAccumWallTime_ +=
                    sharedCacheBuildProfile.classAccumWallTime;
                profileSharedCacheDetailSampleCells_ +=
                    sharedCacheBuildProfile.detailSampleCells;
                profileSharedCacheDetailSampleParcels_ +=
                    sharedCacheBuildProfile.detailSampleParcels;
            }

            const auto fieldCombineStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) if (useOpenMPSampling)
            #endif
            forAll(cellOccupancy, celli)
            {
                scalar dsmcNLocal = 0.0;
                scalar dsmcMLocal = 0.0;
                scalar dsmcLinearKELocal = 0.0;
                vector dsmcMomentumLocal = vector::zero;
                scalar dsmcErotLocal = 0.0;
                scalar dsmcZetaRotLocal = 0.0;
                scalar dsmcNElecLvlLocal = 0.0;
                scalar nLocal = 0.0;
                scalar mLocal = 0.0;
                vector momentumLocal = vector::zero;
                scalar linearKELocal = 0.0;
                scalar dsmcMuuLocal = 0.0;
                scalar dsmcMuvLocal = 0.0;
                scalar dsmcMuwLocal = 0.0;
                scalar dsmcMvvLocal = 0.0;
                scalar dsmcMvwLocal = 0.0;
                scalar dsmcMwwLocal = 0.0;
                scalar dsmcMccLocal = 0.0;
                scalar dsmcMccuLocal = 0.0;
                scalar dsmcMccvLocal = 0.0;
                scalar dsmcMccwLocal = 0.0;
                scalar dsmcEuLocal = 0.0;
                scalar dsmcEvLocal = 0.0;
                scalar dsmcEwLocal = 0.0;
                scalar dsmcEIntLocal = 0.0;
                scalar dsmcNClassILocal = 0.0;
                scalar dsmcNClassIILocal = 0.0;
                scalar dsmcNClassIIILocal = 0.0;

                forAll(speciesIds_, i)
                {
                    const label spId = speciesIds_[i];

                    dsmcNLocal += sharedSampleCache_.dsmcN[spId][celli];
                    dsmcMLocal += sharedSampleCache_.dsmcM[spId][celli];
                    dsmcLinearKELocal += sharedSampleCache_.dsmcLinearKE[spId][celli];
                    dsmcMomentumLocal += sharedSampleCache_.dsmcMomentum[spId][celli];
                    dsmcErotLocal += sharedSampleCache_.dsmcErot[spId][celli];
                    dsmcZetaRotLocal += sharedSampleCache_.dsmcZetaRot[spId][celli];
                    dsmcNElecLvlLocal += sharedSampleCache_.dsmcNElecLvl[spId][celli];
                    nLocal += sharedSampleCache_.nReal[spId][celli];
                    mLocal += sharedSampleCache_.mReal[spId][celli];
                    momentumLocal += sharedSampleCache_.momentumReal[spId][celli];
                    linearKELocal += sharedSampleCache_.linearKEReal[spId][celli];
                    dsmcMuuLocal += sharedSampleCache_.dsmcMuu[spId][celli];
                    dsmcMuvLocal += sharedSampleCache_.dsmcMuv[spId][celli];
                    dsmcMuwLocal += sharedSampleCache_.dsmcMuw[spId][celli];
                    dsmcMvvLocal += sharedSampleCache_.dsmcMvv[spId][celli];
                    dsmcMvwLocal += sharedSampleCache_.dsmcMvw[spId][celli];
                    dsmcMwwLocal += sharedSampleCache_.dsmcMww[spId][celli];
                    dsmcMccLocal += sharedSampleCache_.dsmcMcc[spId][celli];
                    dsmcMccuLocal += sharedSampleCache_.dsmcMccu[spId][celli];
                    dsmcMccvLocal += sharedSampleCache_.dsmcMccv[spId][celli];
                    dsmcMccwLocal += sharedSampleCache_.dsmcMccw[spId][celli];
                    dsmcEuLocal += sharedSampleCache_.dsmcEu[spId][celli];
                    dsmcEvLocal += sharedSampleCache_.dsmcEv[spId][celli];
                    dsmcEwLocal += sharedSampleCache_.dsmcEw[spId][celli];
                    dsmcEIntLocal += sharedSampleCache_.dsmcECum[spId][celli];

                    if (needElectronic)
                    {
                        dsmcSpeciesEelecCum_[i][celli] +=
                            sharedSampleCache_.dsmcSpeciesEelec[spId][celli];
                    }
                    dsmcNSpeciesCum_[i][celli] +=
                        sharedSampleCache_.dsmcN[spId][celli];
                    dsmcMccSpeciesCum_[i][celli] +=
                        sharedSampleCache_.dsmcLinearKE[spId][celli];
                    nSpeciesCum_[i][celli] +=
                        sharedSampleCache_.nReal[spId][celli];

                    if (needElectronic)
                    {
                        dsmcNGrndElecLvlSpeciesCum_[i][celli] +=
                            sharedSampleCache_.nGrndElecLvl[spId][celli];
                        dsmcN1stElecLvlSpeciesCum_[i][celli] +=
                            sharedSampleCache_.n1stElecLvl[spId][celli];
                    }

                    if (needVibrational)
                    {
                        forAll(dsmcSpeciesEvibModCum_[i], mod)
                        {
                            dsmcSpeciesEvibModCum_[i][mod][celli] +=
                                sharedSampleCache_.dsmcSpeciesEvibMod[spId][mod][celli];
                        }
                    }

                    if (needClassification)
                    {
                        dsmcNClassILocal += sharedSampleCache_.dsmcNClassI[spId][celli];
                        dsmcNClassIILocal += sharedSampleCache_.dsmcNClassII[spId][celli];
                        dsmcNClassIIILocal += sharedSampleCache_.dsmcNClassIII[spId][celli];
                    }
                }

                dsmcNCum_[celli] += dsmcNLocal;
                dsmcN_[celli] += dsmcNLocal;
                dsmcMCum_[celli] += dsmcMLocal;
                dsmcLinearKECum_[celli] += dsmcLinearKELocal;
                dsmcMomentumCum_[celli] += dsmcMomentumLocal;
                dsmcErotCum_[celli] += dsmcErotLocal;
                dsmcZetaRotCum_[celli] += dsmcZetaRotLocal;
                dsmcNElecLvlCum_[celli] += dsmcNElecLvlLocal;
                nCum_[celli] += nLocal;
                mCum_[celli] += mLocal;
                momentumCum_[celli] += momentumLocal;
                linearKECum_[celli] += linearKELocal;
                dsmcMuuCum_[celli] += dsmcMuuLocal;
                dsmcMuvCum_[celli] += dsmcMuvLocal;
                dsmcMuwCum_[celli] += dsmcMuwLocal;
                dsmcMvvCum_[celli] += dsmcMvvLocal;
                dsmcMvwCum_[celli] += dsmcMvwLocal;
                dsmcMwwCum_[celli] += dsmcMwwLocal;
                dsmcMccCum_[celli] += dsmcMccLocal;
                dsmcMccuCum_[celli] += dsmcMccuLocal;
                dsmcMccvCum_[celli] += dsmcMccvLocal;
                dsmcMccwCum_[celli] += dsmcMccwLocal;
                dsmcEuCum_[celli] += dsmcEuLocal;
                dsmcEvCum_[celli] += dsmcEvLocal;
                dsmcEwCum_[celli] += dsmcEwLocal;
                dsmcECum_[celli] += dsmcEIntLocal;

                if (needClassification)
                {
                    dsmcNClassICum_[celli] += dsmcNClassILocal;
                    dsmcNClassIICum_[celli] += dsmcNClassIILocal;
                    dsmcNClassIIICum_[celli] += dsmcNClassIIILocal;
                }
            }

            if (doProfile)
            {
                profileFieldCombineWallTime_ +=
                    wallSeconds(fieldCombineStart, wallClockNow());
                profileSampleAccumWallTime_ +=
                    wallSeconds(sampleAccumStart, wallClockNow());
            }

            //- Loop over all cells
            const auto cellReduceStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            #ifdef _OPENMP
            if (useOpenMPSampling)
            {
                #pragma omp parallel for schedule(static)
                for (label celli = 0; celli < dsmcNCum_.size(); ++celli)
                {
                    collisionSeparation_[celli] +=
                        cloud_.cellPropMeasurements().collisionSeparation()[celli];

                    dsmcNCollsCum_[celli] +=
                        cloud_.cellPropMeasurements().nColls()[celli];

                    if (dsmcNCum_[celli] > 1e-3)
                    {
                        const scalar cellVolume = mesh_.cellVolumes()[celli];

                        dsmcNMean_[celli] = dsmcNCum_[celli]/nAvTimeSteps;

                        const scalar rhoNMean = nCum_[celli]
                            /(nAvTimeSteps*cellVolume);
                        const scalar rhoMMean = mCum_[celli]
                            /(nAvTimeSteps*cellVolume);

                        rhoN_[celli] = rhoNMean;
                        rhoM_[celli] = rhoMMean;

                        UMean_[celli] = momentumCum_[celli]/mCum_[celli];

                        const scalar linearKEMean = 0.5*linearKECum_[celli]
                            /(cellVolume*nAvTimeSteps);

                        Ttra_[celli] =
                            2.0/(3.0*kB*rhoNMean)
                           *(
                                linearKEMean - 0.5*rhoMMean
                               *(
                                    UMean_[celli] & UMean_[celli]
                                )
                            );

                        p_[celli] = rhoNMean*kB*Ttra_[celli];
                    }
                    else
                    {
                        dsmcNMean_[celli] = 0.001;
                        rhoN_[celli] = 0.0;
                        rhoM_[celli] = 0.0;
                        UMean_[celli] = vector::zero;
                        Ttra_[celli] = 0.0;
                        p_[celli] = 0.0;
                    }
                }
            }
            else
            #endif
            forAll(dsmcNCum_, celli)
            {
                collisionSeparation_[celli] +=
                    cloud_.cellPropMeasurements().collisionSeparation()[celli];
                    
                dsmcNCollsCum_[celli] +=
                    cloud_.cellPropMeasurements().nColls()[celli];

                if (dsmcNCum_[celli] > 1e-3)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[celli];

                    dsmcNMean_[celli] = dsmcNCum_[celli]/nAvTimeSteps;

                    const scalar rhoNMean = nCum_[celli]
                        /(nAvTimeSteps*cellVolume);
                    const scalar rhoMMean = mCum_[celli]
                        /(nAvTimeSteps*cellVolume);

                    rhoN_[celli] = rhoNMean;
                    rhoM_[celli] = rhoMMean;
                    
                    UMean_[celli] = momentumCum_[celli]/mCum_[celli];

                    const scalar linearKEMean = 0.5*linearKECum_[celli]
                        /(cellVolume*nAvTimeSteps);

                    Ttra_[celli] =
                        2.0/(3.0*kB*rhoNMean)
                       *(
                            linearKEMean - 0.5*rhoMMean
                           *(
                                UMean_[celli] & UMean_[celli]
                            )
                        );

                    p_[celli] = rhoNMean*kB*Ttra_[celli];
                }
                else
                {
                    dsmcNMean_[celli] = 0.001;
                    rhoN_[celli] = 0.0;
                    rhoM_[celli] = 0.0;
                    UMean_[celli] = vector::zero;
                    Ttra_[celli] = 0.0;
                    p_[celli] = 0.0;
                }
            }

            if (doProfile)
            {
                profileCellReduceWallTime_ +=
                    wallSeconds(cellReduceStart, wallClockNow());
            }

            //- Obtain boundary measurements
            const auto boundaryAccumStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            forAll(speciesIds_, i)
            {
                const label spId = speciesIds_[i];
                
                forAll(mesh_.boundaryMesh(), j)
                {
                    forAll(mesh_.boundaryMesh()[j], k)
                    {
                        rhoNBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoNBF(spId, j, k);
                        rhoMBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoMBF(spId, j, k);
                        linearKEBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesLinearKEBF(spId, j, k);
                        momentumBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesMomentumBF(spId, j, k);
                        ErotBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesErotBF(spId, j, k);
                        zetaRotBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesZetaRotBF(spId, j, k);
                        rhoNIntBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoNIntBF(spId, j, k); 
                        rhoNElecBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoNElecBF(spId, j, k);
                        qBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesqBF(spId, j, k);
                        fDBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesfDBF(spId, j, k);
                        
                        speciesRhoNBF_[i][j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoNBF(spId, j, k);
                        speciesEvibBF_[i][j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesEvibBF(spId, j, k);
                        speciesEelecBF_[i][j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesEelecBF(spId, j, k);
                        speciesMccBF_[i][j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesMccBF(spId, j, k);
                    }
                }

                forAll(speciesEvibModBF_[i], mod)
                {
                    forAll(mesh_.boundaryMesh(), j)
                    {
                        forAll(mesh_.boundaryMesh()[j], k)
                        {
                            speciesEvibModBF_[i][mod][j][k] +=
                                cloud_.boundaryFluxMeasurements()
                                  .speciesEvibModBF(spId, mod, j, k);
                        }
                    }
                }
            }

            if (doProfile)
            {
                profileBoundaryAccumWallTime_ +=
                    wallSeconds(boundaryAccumStart, wallClockNow());
            }
        }

        sampleCounter_ = 0;
    }

    // XCX: if time > resetAtOutputUntilTime_, means reach steady state
    // clear field and start sample
    if
    (
          time_.time().value()
        + time_.time().deltaT().value()
        > time_.resetAtOutputUntilTime()
        && !isSample_
    )
    {
        // Info << "resetFieldsAtOutput is set to false for the rest of the simulation" << endl;
        // Info << "time " << time_.time().value() << " > resetAtOutputUntilTime " << resetAtOutputUntilTime << endl;
        time_.resetFieldsAtOutput() = false;
        isSample_ = true;

        //- Reset fields after printing the instantaneous solution ... or continue sampling
        Info << endl;
        Info << "Start sample at time " << time_.time().value()<< " and step " << nTimeSteps_ << endl;
        Info << endl;

        nTimeSteps_ = 0.0;
        
        forAll(dsmcNCum_, celli)
        {
            dsmcNCum_[celli] = 0.0;
            dsmcMCum_[celli] = 0.0;
            dsmcLinearKECum_[celli] = 0.0;
            dsmcMomentumCum_[celli] = vector::zero;
            dsmcErotCum_[celli] = 0.0;
            dsmcZetaRotCum_[celli] = 0.0;
            dsmcNElecLvlCum_[celli] = 0.0,
            dsmcNClassICum_[celli] = 0.0;
            dsmcNClassIICum_[celli] = 0.0;
            dsmcNClassIIICum_[celli] = 0.0;
            collisionSeparation_[celli] = 0.0;
            dsmcNCollsCum_[celli] = 0.0;
            measuredCollisionRate_[celli] = 0.0;
            dsmcMuuCum_[celli] = 0.0;
            dsmcMuvCum_[celli] = 0.0;
            dsmcMuwCum_[celli] = 0.0;
            dsmcMvvCum_[celli] = 0.0;
            dsmcMvwCum_[celli] = 0.0;
            dsmcMwwCum_[celli] = 0.0;
            dsmcMccCum_[celli] = 0.0;
            dsmcMccuCum_[celli] = 0.0;
            dsmcMccvCum_[celli] = 0.0;
            dsmcMccwCum_[celli] = 0.0;
            dsmcEuCum_[celli] = 0.0;
            dsmcEvCum_[celli] = 0.0;
            dsmcEwCum_[celli] = 0.0;
            dsmcECum_[celli] = 0.0;
            zetaVib_[celli] = 0.0;
            nCum_[celli] = 0.0;
            mCum_[celli] = 0.0;
            momentumCum_[celli] = vector::zero;
            linearKECum_[celli] = 0.0;
        }
        forAll(speciesIds_, i)
        {
            forAll(speciesTvib_[i], celli)
            {
                dsmcNSpeciesCum_[i][celli] = 0.0;
                nSpeciesCum_[i][celli] = 0.0;
                dsmcMccSpeciesCum_[i][celli] = 0.0;
                speciesMfp_[i][celli] = 0.0;
                speciesMcr_[i][celli] = 0.0;
                speciesTvib_[i][celli] = 0.0;
                dsmcSpeciesEelecCum_[i][celli] = 0.0;
                dsmcNGrndElecLvlSpeciesCum_[i][celli] = 0.0;
                dsmcN1stElecLvlSpeciesCum_[i][celli] = 0.0;
            }
            forAll(dsmcSpeciesEvibModCum_[i], mod)
            {
                forAll(dsmcSpeciesEvibModCum_[i][mod], celli)
                {
                    dsmcSpeciesEvibModCum_[i][mod][celli] = 0.0;
                }
            }
        }
        //- Reset boundary information
        forAll(rhoNBF_, j)
        {
            rhoNBF_[j] = 0.0;
            rhoMBF_[j] = 0.0;
            linearKEBF_[j] = 0.0;
            rhoNIntBF_[j] = 0.0;
            
            ErotBF_[j] = 0.0;
            zetaRotBF_[j] = 0.0;
            zetaVibBF_[j] = 0.0;
            rhoNElecBF_[j] = 0.0;
            
            qBF_[j] = 0.0;
            fDBF_[j] = vector::zero;
            momentumBF_[j] = vector::zero;
        }
        forAll(speciesIds_, i)
        {
            forAll(speciesTvibBF_[i], j)
            {
                speciesRhoNBF_[i][j] = 0.0;
                speciesMccBF_[i][j] = 0.0;
                speciesEvibBF_[i][j] = 0.0;
                speciesTvibBF_[i][j] = 0.0;
                speciesZetaVibBF_[i][j] = 0.0;
                speciesEelecBF_[i][j] = 0.0;
            }
            
            forAll(speciesEvibModBF_[i], mod)
            {
                forAll(speciesEvibModBF_[i][mod], j)
                {
                    speciesEvibModBF_[i][mod][j] = 0.0;
                }
            }
        }
    }

    

    if (time_.time().outputTime())
    {
        const auto outputTimeStart =
            doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
        const scalar nAvTimeSteps = nTimeSteps_;

        if (Pstream::myProcNo() == 0)
            Info << endl << "Sample average steps = " << nAvTimeSteps << endl;

        //- XCX: get freestreaminflow velocity and density for heat and pressure coeffcient
        IOdictionary dsmcInitialiseDict
        (
            IOobject
            (
                "dsmcInitialiseDict",
                mesh_.time().system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
        IOdictionary particleProperties_
        (
            IOobject
            (
                "dsmcProperties",
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        const PtrList<entry> configurations(Foam::hyCompat::lookup(Foam::hyCompat::lookup(dsmcInitialiseDict, "configurations")));

        if (configurations.size() > 0)
        {
            const dictionary& configurationDict = configurations[0].dict();

            freeStreamVelocity_ = configurationDict.get<vector>("velocity");
            freeStreamVelocityValue_ = mag(freeStreamVelocity_);
            freeStreamDensity_ = 0.0;
            totalWallForce_ = vector::zero;

            const dictionary& moleculePropertiesDict = particleProperties_.subDict("moleculeProperties");
            const dictionary& numberDensitiesDict = configurationDict.subDict("numberDensities");

            scalar totalNumberDensity = 0.0;

            forAllConstIters(numberDensitiesDict, iter)
            {
                const word& speciesName = iter().keyword();
                const scalar numberDensity = numberDensitiesDict.get<scalar>(speciesName);
                totalNumberDensity += numberDensity;
                freeStreamDensity_ +=
                    numberDensity*moleculePropertiesDict.subDict(speciesName).get<scalar>("mass");
            }

            if
            (
                totalNumberDensity > SMALL
             && freeStreamDensity_ > SMALL
             && freeStreamVelocityValue_ > SMALL
            )
            {
                const scalar u2 = sqr(freeStreamVelocityValue_);
                cqCoef_ = 2.0/(freeStreamDensity_*u2*freeStreamVelocityValue_);
                cpCoef_ = 2.0/(freeStreamDensity_*u2);
            }
            else
            {
                cqCoef_ = 0.0;
                cpCoef_ = 0.0;

                if (Pstream::myProcNo() == 0)
                {
                    Info<< endl
                        << "Vacuum or zero-velocity freestream condition. "
                        << "Cq and Cp coefficients are set to zero." << endl;
                }
            }
        }

        if (densityOnly_)
        {
            forAll(dsmcNCum_, celli)
            {
                if (dsmcNCum_[celli] > SMALL)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[celli];

                    dsmcNMean_[celli] = dsmcNCum_[celli]/nAvTimeSteps;

                    rhoN_[celli] = nCum_[celli]/(nAvTimeSteps*cellVolume);
                    rhoM_[celli] = mCum_[celli]/(nAvTimeSteps*cellVolume);
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcNMean_[celli] = 0.001;
                    rhoN_[celli] = 0.0;
                    rhoM_[celli] = 0.0;
                }

                if (dsmcN_[celli] < SMALL)
                {
                    // not zero so that weighted decomposition still works
                    dsmcN_[celli] = 0.001;
                }
            }
        }
        else
        {
            scalarField& MaInternal = Ma_.primitiveFieldRef();

            forAll(dsmcNCum_, celli)
            {
                //- Fields initialisation 
                scalar moleculesRhoN = 0.0;
                Tvib_[celli] = 0.0;
                zetaVib_[celli] = 0.0;
                
                scalar molarCv_trarot = 0.0;
                scalar molarCp_trarot = 0.0;
                scalar molecularMass = 0.0;
                scalar particleCv = 0.0;
                scalar gamma = 0.0;
            
                //- Rotational energy mode
                const scalar zetaRotTot
                (
                    dsmcNCum_[celli] > SMALL
                  ? dsmcZetaRotCum_[celli]/dsmcNCum_[celli]
                  : 0.0
                );

                Trot_[celli] =
                (
                    dsmcZetaRotCum_[celli] > SMALL
                  ? 2.0*dsmcErotCum_[celli]/(kB*dsmcZetaRotCum_[celli])
                  : 0.0
                );

                //- Vibrational energy mode
                forAll(speciesIds_, i)
                {
                    const label spId = speciesIds_[i];
                    const scalar speciesCount = nSpeciesCum_[i][celli];

                    if (speciesCount <= SMALL)
                    {
                        continue;
                    }

                    scalar speciesZetaVib = 0.0;
                    scalar zetaByTvibMod = 0.0;
                    const scalarList& thetaV = cloud_.constProps(spId).thetaV();

                    forAll(dsmcSpeciesEvibModCum_[i], mod)
                    {
                        const scalar evibCum =
                            dsmcSpeciesEvibModCum_[i][mod][celli];

                        if
                        (
                            evibCum > VSMALL
                         && thetaV[mod] > SMALL
                        )
                        {
                            const scalar iMean =
                                evibCum/(kB*thetaV[mod]*speciesCount);
                               
                            if (iMean > SMALL)
                            {
                                const scalar logFactor = log(1.0 + 1.0/iMean);
                                const scalar speciesTvibMod =
                                    thetaV[mod]/logFactor;
                                const scalar speciesZetaVibMod =
                                    2.0*iMean*logFactor;

                                speciesZetaVib += speciesZetaVibMod;
                                    
                                // XCX: = should be +=, same as line 2251
                                zetaByTvibMod +=
                                    speciesZetaVibMod*speciesTvibMod;
                            }
                        }
                    }

                    if (speciesZetaVib > SMALL)
                    {
                        const scalar speciesTvib =
                            zetaByTvibMod/speciesZetaVib;

                        moleculesRhoN += speciesCount;
                        Tvib_[celli] += speciesCount*speciesTvib;
                            
                        zetaVib_[celli] += speciesCount*speciesZetaVib;    
                    }
                    
                } //- end species loop

                if (moleculesRhoN > SMALL)
                {
                    Tvib_[celli] /= moleculesRhoN;
                    zetaVib_[celli] /= moleculesRhoN;
                }

                //- Electronic energy mode // TODO Vincent
                //  To reintroduce - I do not trust this part
                scalar zetaElecTot = 0.0;
                Telec_[celli] = 0.0;

                //- Overall temperature
                Tov_[celli] =
                    (
                        3.0*Ttra_[celli]
                      + zetaRotTot*Trot_[celli]
                      + zetaVib_[celli]*Tvib_[celli]
                      + zetaElecTot*Telec_[celli]
                    ) /
                    (3.0 + zetaRotTot + zetaVib_[celli] + zetaElecTot);


                if (measureHeatFluxShearStress_)
                {
                    if (dsmcNCum_[celli] > SMALL)
                    {
                        pressureTensor_[celli].xx() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                dsmcMuuCum_[celli]
                              - dsmcMCum_[celli]*sqr(UMean_[celli].x())
                            );
                        pressureTensor_[celli].xy() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                dsmcMuvCum_[celli]
                              - dsmcMCum_[celli]*UMean_[celli].x()
                              * UMean_[celli].y()
                            );
                        pressureTensor_[celli].xz() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                dsmcMuwCum_[celli]
                              - dsmcMCum_[celli]*UMean_[celli].x()
                              * UMean_[celli].z()
                            );

                        pressureTensor_[celli].yx() =
                            pressureTensor_[celli].xy();
                        pressureTensor_[celli].yy() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                dsmcMvvCum_[celli]
                              - dsmcMCum_[celli]*sqr(UMean_[celli].y())
                            );
                        pressureTensor_[celli].yz() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                dsmcMvwCum_[celli]
                              - dsmcMCum_[celli]*UMean_[celli].y()
                              * UMean_[celli].z()
                            );

                        pressureTensor_[celli].zx() =
                            pressureTensor_[celli].xz();
                        pressureTensor_[celli].zy() =
                            pressureTensor_[celli].yz();
                        pressureTensor_[celli].zz() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                dsmcMwwCum_[celli]
                              - dsmcMCum_[celli]*sqr(UMean_[celli].z())
                            );

                        const scalar scalarPressure =
                            1.0/3.0
                           *(
                                pressureTensor_[celli].xx()
                              + pressureTensor_[celli].yy()
                              + pressureTensor_[celli].zz()
                            );

                        shearStressTensor_[celli] = -pressureTensor_[celli];
                        shearStressTensor_[celli].xx() += scalarPressure;
                        shearStressTensor_[celli].yy() += scalarPressure;
                        shearStressTensor_[celli].zz() += scalarPressure;

                        //- terms involving pressure tensor should not be
                        //  multiplied by the number density
                        //  (see Bird corrigendum)

                        heatFluxVector_[celli].x() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                0.5*dsmcMccuCum_[celli]
                              - 0.5*dsmcMccCum_[celli]*UMean_[celli].x()
                              + dsmcEuCum_[celli]
                              - dsmcECum_[celli]*UMean_[celli].x()
                            )
                          - pressureTensor_[celli].xx()*UMean_[celli].x()
                          - pressureTensor_[celli].xy()*UMean_[celli].y()
                          - pressureTensor_[celli].xz()*UMean_[celli].z();

                        heatFluxVector_[celli].y() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                0.5*dsmcMccvCum_[celli]
                              - 0.5*dsmcMccCum_[celli]*UMean_[celli].y()
                              + dsmcEvCum_[celli]
                              - dsmcECum_[celli]*UMean_[celli].y()
                            )
                          - pressureTensor_[celli].yx()*UMean_[celli].x()
                          - pressureTensor_[celli].yy()*UMean_[celli].y()
                          - pressureTensor_[celli].yz()*UMean_[celli].z();

                        heatFluxVector_[celli].z() =
                            rhoN_[celli]/dsmcNCum_[celli]
                           *(
                                0.5*dsmcMccwCum_[celli]
                              - 0.5*dsmcMccCum_[celli]*UMean_[celli].z()
                              + dsmcEwCum_[celli]
                              - dsmcECum_[celli]*UMean_[celli].z()
                            )
                          - pressureTensor_[celli].zx()*UMean_[celli].x()
                          - pressureTensor_[celli].zy()*UMean_[celli].y()
                          - pressureTensor_[celli].zz()*UMean_[celli].z();
                    }
                    else
                    {
                        pressureTensor_[celli] = tensor::zero;
                        shearStressTensor_[celli] = tensor::zero;
                        heatFluxVector_[celli] = vector::zero;
                    }
                }                if (dsmcNCum_[celli] > SMALL && Ttra_[celli] > SMALL)
                {
                    forAll(speciesIds_, i)
                    {
                        const label spId = speciesIds_[i];
                        const scalar speciesZetaRot =
                            cloud_.constProps(spId)
                              .rotationalDegreesOfFreedom();
                        
                        const scalar Xs = nSpeciesCum_[i][celli]
                            /nCum_[celli];

                        molecularMass += Xs*cloud_.constProps(spId).mass();
                            
                         //- Heat capacity at constant volume/(0.5*kB)
                         //  trans-rotational
                         molarCv_trarot += Xs*(3.0 + speciesZetaRot);
                            
                        //- Heat capacity at constant pressure/(0.5*kB)
                        //  trans-rotational
                        molarCp_trarot += Xs*(5.0 + speciesZetaRot);
                    }

                    particleCv = molarCv_trarot/NAvo;

                    if (molarCv_trarot > SMALL && molarCp_trarot > SMALL && molecularMass > VSMALL)
                    {
                        gamma = molarCp_trarot/molarCv_trarot;

                        const scalar speedOfSound = sqrt
                        (
                            max(gamma*kB/molecularMass*Ttra_[celli], scalar(0.0))
                        );

                        if (speedOfSound > SMALL)
                        {
                            MaInternal[celli] = mag(UMean_[celli])/speedOfSound;
                        }
                        else
                        {
                            MaInternal[celli] = 0.0;
                        }
                    }
                    else
                    {
                        MaInternal[celli] = 0.0;
                    }
                }
                else
                {
                    MaInternal[celli] = 0.0;
                }

                if (measureMeanFreePath_ && Ttra_[celli] > 1.0)
                {
                    const scalar deltaT = cloud_.deltaTValue(celli);
                    
                    mfp_[celli] = 0.0;
                    meanCollisionRate_[celli] = 0.0;
                    
                    forAll(speciesIds_, s)
                    {
                        const label spIdp = speciesIds_[s];
                        
                        speciesMfp_[s][celli] = 0.0;
                        speciesMcr_[s][celli] = 0.0;

                        forAll(speciesIds_, r)
                        {
                            const label spIdq = speciesIds_[r];
                            
                            const scalar dPQ =
                                0.5*
                                (
                                    cloud_.constProps(spIdp).d()
                                  + cloud_.constProps(spIdq).d()
                                );
                            const scalar omegaPQ =
                                0.5*
                                (
                                    cloud_.constProps(spIdp).omega()
                                  + cloud_.constProps(spIdq).omega()
                                );
                            const scalar massRatio =
                                cloud_.constProps(spIdp).mass()
                               /cloud_.constProps(spIdq).mass();

                            if
                            (
                                dsmcNSpeciesCum_[r][celli] > SMALL
                             && Ttra_[celli] > SMALL
                            )
                            {
                                const scalar nDensQ =
                                    nSpeciesCum_[r][celli]
                                   /(mesh_.cellVolumes()[celli]*nAvTimeSteps);
                                const scalar reducedMass =
                                    cloud_.constProps(spIdp).mass()
                                   *cloud_.constProps(spIdq).mass()
                                   /
                                    (
                                       cloud_.constProps(spIdp).mass()
                                     + cloud_.constProps(spIdq).mass()
                                    );

                                // Bird 1994, eq (4.76)
                                speciesMfp_[s][celli] += pi*sqr(dPQ)*nDensQ
                                   *pow
                                    (
                                        mfpTref_/Ttra_[celli], omegaPQ - 0.5
                                    )*sqrt(1.0+massRatio);

                                // Bird 1994, eq (4.74)
                                speciesMcr_[s][celli] +=
                                    2.0*sqrt(pi)*sqr(dPQ)*nDensQ
                                   *pow
                                    (
                                        Ttra_[celli]/mfpTref_, 1.0 - omegaPQ
                                    )
                                   *sqrt
                                    (
                                        2.0*kB*mfpTref_/reducedMass
                                    );
                            }
                        }

                        if (speciesMfp_[s][celli] > SMALL)
                        {
                            speciesMfp_[s][celli] = 1.0/speciesMfp_[s][celli];
                        }
                    }

                    meanCollisionSeparation_[celli] =
                    (
                        dsmcNCollsCum_[celli] > SMALL
                      ? collisionSeparation_[celli]/dsmcNCollsCum_[celli]
                      : GREAT
                    );

                    if (nCum_[celli] > SMALL)
                    {
                        // const scalar symmFactor = 2.0;
                        // TODO (s == r ? 1.0 : 2.0);
                        measuredCollisionRate_[celli] = dsmcNCollsCum_[celli]
                            *cloud_.nParticles(celli)/(nCum_[celli]*deltaT);
                    }

                    if (rhoN_[celli] > SMALL)
                    {
                        forAll(speciesIds_, i)
                        {
                            const scalar rhoNi = nSpeciesCum_[i][celli];

                            // Bird 1994, eq (4.77)
                            mfp_[celli] += speciesMfp_[i][celli]
                                *rhoNi/nCum_[celli];

                            // Bird 1994, eq (1.38)
                            meanCollisionRate_[celli] +=
                                speciesMcr_[i][celli]*rhoNi/nCum_[celli];
                        }
                    }

                    if (mfp_[celli] < SMALL)
                    {
                        mfp_[celli] = GREAT;
                    }

                    if (meanCollisionRate_[celli] > SMALL)
                    {
                        meanCollisionTime_[celli] =
                            1.0/meanCollisionRate_[celli];
                        mctToDt_[celli] = meanCollisionTime_[celli]/deltaT;
                    }
                    else
                    {
                        meanCollisionTime_[celli] = GREAT;
                        mctToDt_[celli] = GREAT;
                    }

                    if (mfp_[celli] != GREAT)
                    {
                        scalar maxCellDx = 0.0;
                        scalarField cellDx(3, 0.0);

                        const labelList& pLabels
                        (
                            mesh_.cells()[celli].labels(mesh_.faces())
                        );
                        pointField pLocal(pLabels.size(), vector::zero);

                        forAll (pLabels, pointi)
                        {
                            pLocal[pointi] = mesh_.points()[pLabels[pointi]];
                        }

                        cellDx[0] = Foam::max(pLocal & vector(1,0,0))
                            - Foam::min(pLocal & vector(1,0,0));
                        cellDx[1] = Foam::max(pLocal & vector(0,1,0))
                            - Foam::min(pLocal & vector(0,1,0));
                        cellDx[2] = Foam::max(pLocal & vector(0,0,1))
                            - Foam::min(pLocal & vector(0,0,1));

                        maxCellDx = cellDx[0];

                        forAll(cellDx, dim)
                        {
                            if (cellDx[dim] > maxCellDx)
                            {
                                maxCellDx = cellDx[dim];
                            }
                        }

                        mfpToDx_[celli] = mfp_[celli]/maxCellDx;

                        SOF_[celli] =
                        (
                            mfp_[celli] > SMALL
                          ? meanCollisionSeparation_[celli]/mfp_[celli]
                          : 0.0
                        );
                    }
                    else
                    {
                        mfpToDx_[celli] = GREAT;
                        SOF_[celli] = GREAT;
                    }

                    // when few particles in cell, undesired refinement
                    // this condition should eliminates this problem
                    if (dsmcN_[celli] >= 4.0)
                    {
                        DxToMfp_[celli] = (mfpToDx_[celli] > SMALL ? 1.0/mfpToDx_[celli] : GREAT);
                    }
                }

                if (measureClassifications_)
                {
                    if (dsmcNCum_[celli] > SMALL)
                    {
                        classIDistribution_[celli] = dsmcNClassICum_[celli]
                            /dsmcNCum_[celli];
                        classIIDistribution_[celli] = dsmcNClassIICum_[celli]
                            /dsmcNCum_[celli];
                        classIIIDistribution_[celli] = dsmcNClassIIICum_[celli]
                            /dsmcNCum_[celli];
                    }
                }

                if (measureErrors_)
                {
                    if
                    (
                         dsmcNMean_[celli] > SMALL && Ma_[celli] > SMALL
                      && gamma > SMALL && particleCv > SMALL
                    )
                    {
                        const scalar deno = sqrt(dsmcNMean_[celli]*nAvTimeSteps);
                        
                        densityError_[celli] = 1.0/deno;
                        velocityError_[celli] = 1.0/(deno*Ma_[celli]*sqrt(gamma));
                        temperatureError_[celli] = sqrt(kB/particleCv)/deno;
                        pressureError_[celli] = sqrt(gamma)/deno;
                    }

                }
            } //- end loop over cells

            //- Computing boundary measurements: loop over all boundary patches
            forAll(rhoNBF_, j)
            {
                //- Determine of the type of patch: patch, wall, cyclic, ...
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                
                const bool isWall = isA<wallPolyPatch>(patch);
                
                const bool isNonEmptyNonCyclic = isA<polyPatch>(patch)
                    && !isA<emptyPolyPatch>(patch)
                    && !isA<cyclicPolyPatch>(patch);

                if (isWall)
                {
                    //- Loop over all wall boundary faces
                    forAll(patch, k)
                    {
                        const label celli = boundaryCells_[j][k];
                        
                        //- Initialise face fields
                        Tvib_.boundaryFieldRef()[j][k] = 0.0;
                        zetaVibBF_[j][k] = 0.0;
                        scalar molecularMassBF = 0.0;
                        scalar molarCvBF_trarot = 0.0;
                        scalar molarCpBF_trarot = 0.0;
                        
                        // Note: do not use the nParticles value that includes
                        // the RWF here. This is wrong because boundary
                        // measurements are performend during move steps. Hence
                        // radial weighting (if simulation is axisymmetric or
                        // spherical) is performed after the measurement. That
                        // is why the parcel RWFs are already included during
                        // the measurement step and we only need the FNUM value
                        // here.
                        const scalar nParticles = cloud_.coordSystem().dtModel()
                            .nParticles(j, k);

                        const scalar rhoNMean =
                            rhoNBF_[j][k]*nParticles/nAvTimeSteps;
                        const scalar rhoMMean =
                            rhoMBF_[j][k]*nParticles/nAvTimeSteps;
                        const scalar linearKEMean =
                            linearKEBF_[j][k]*nParticles/nAvTimeSteps;

                        rhoN_.boundaryFieldRef()[j][k] = rhoNMean;
                        rhoM_.boundaryFieldRef()[j][k] = rhoMMean;
                        
                        //- Instantaneous and sampled numbers of DSMC parcels
                        //  are that of the neighbouring cell
                        dsmcN_.boundaryFieldRef()[j][k] = dsmcN_[celli];
                        dsmcNMean_.boundaryFieldRef()[j][k] = dsmcNMean_[celli];

                        //- Translational energy mode and velocity
                        if (rhoMMean > VSMALL)
                        {
                            UMean_.boundaryFieldRef()[j][k] = momentumBF_[j][k]
                                /rhoMBF_[j][k];

                            Ttra_.boundaryFieldRef()[j][k] =
                                2.0/(3.0*kB*rhoNMean)*
                                (
                                    linearKEMean - 0.5*rhoMMean*
                                    (
                                        UMean_.boundaryField()[j][k]
                                      & UMean_.boundaryField()[j][k]
                                    )
                                );
                        }
                        else
                        {
                            UMean_.boundaryFieldRef()[j][k] = vector::zero;
                            Ttra_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Rotational energy mode
                        const scalar zetaRotTot =
                        (
                            rhoNBF_[j][k] > SMALL
                          ? zetaRotBF_[j][k]/rhoNBF_[j][k]
                          : 0.0
                        );

                        Trot_.boundaryFieldRef()[j][k] =
                        (
                            zetaRotBF_[j][k] > SMALL
                          ? 2.0*ErotBF_[j][k]/(kB*zetaRotBF_[j][k])
                          : 0.0
                        );

                        //- Vibrational energy mode: loop over all species
                        scalar moleculesRhoN = 0.0;
                        
                        forAll(speciesIds_, i)
                        {
                            const label spId = speciesIds_[i];
                            const scalar speciesRhoN = speciesRhoNBF_[i][j][k];
                            const scalarList& thetaV =
                                cloud_.constProps(spId).thetaV();
                            
                            speciesZetaVibBF_[i][j][k] = 0.0;
                            speciesTvibBF_[i][j][k] = 0.0;
                            
                            scalar zetaByTvibMod = 0.0;

                            if (speciesRhoN > SMALL)
                            {
                                forAll(thetaV, mod)
                                {
                                    const scalar evibMod =
                                        speciesEvibModBF_[i][mod][j][k];

                                    if (evibMod <= VSMALL || thetaV[mod] <= SMALL)
                                    {
                                        continue;
                                    }

                                    const scalar iMean =
                                        evibMod/(kB*thetaV[mod]*speciesRhoN);

                                    if (iMean > SMALL)
                                    {
                                        const scalar logFactor =
                                            log(1.0 + 1.0/iMean);
                                        const scalar speciesTvibMod =
                                            thetaV[mod]/logFactor;
                                        const scalar speciesZetaVibMod =
                                            2.0*iMean*logFactor;

                                        speciesZetaVibBF_[i][j][k] +=
                                            speciesZetaVibMod;
                                            
                                        zetaByTvibMod +=
                                            speciesZetaVibMod*speciesTvibMod;
                                    }
                                }
                            }

                            if (speciesZetaVibBF_[i][j][k] > SMALL)
                            {
                                moleculesRhoN += speciesRhoN;
                                
                                speciesTvibBF_[i][j][k] = zetaByTvibMod
                                    /speciesZetaVibBF_[i][j][k];
                                    
                                Tvib_.boundaryFieldRef()[j][k] +=
                                    speciesRhoN
                                   *speciesTvibBF_[i][j][k];

                                zetaVibBF_[j][k] +=
                                    speciesRhoN
                                   *speciesZetaVibBF_[i][j][k];
                            }
                        }

                        if (moleculesRhoN > SMALL)
                        {
                            Tvib_.boundaryFieldRef()[j][k] /= moleculesRhoN;
                            zetaVibBF_[j][k] /= moleculesRhoN;
                        }

                        //- Electronic energy mode // TODO Vincent
                        //  Removed temporarily - I don't trust this part
                        scalar zetaElecTot = 0.0;
                        Telec_.boundaryFieldRef()[j][k] = 0.0;

                        Tov_.boundaryFieldRef()[j][k] =
                            (
                                (3.0*Ttra_.boundaryField()[j][k])
                              + (zetaRotTot*Trot_.boundaryField()[j][k])
                              + (
                                    zetaVibBF_[j][k]
                                   *Tvib_.boundaryField()[j][k]
                                )
                              + (
                                    zetaElecTot
                                   *Telec_.boundaryFieldRef()[j][k]
                                )
                            )
                           /(
                                3.0 + zetaRotTot + zetaVibBF_[j][k]
                              + zetaElecTot
                            );                        if (rhoNBF_[j][k] > SMALL)
                        {
                            //- Loop over all species
                            forAll(speciesIds_, i)
                            {
                                const label spId = speciesIds_[i];
                                
                                const scalar speciesZetaRotBF =
                                    cloud_.constProps(spId)
                                      .rotationalDegreesOfFreedom();

                                const scalar Xs =
                                    speciesRhoNBF_[i][j][k]/rhoNBF_[j][k];

                                molecularMassBF += Xs
                                    *cloud_.constProps(spId).mass();

                                //- Heat capacity at constant volume/(0.5*kB)
                                //  trans-rotational energy mode
                                molarCvBF_trarot += Xs*(3.0 + speciesZetaRotBF);

                                //- Heat capacity at constant pressure/(0.5*kB)
                                //  trans-rotational energy mode
                                molarCpBF_trarot += Xs*(5.0 + speciesZetaRotBF);
                            }

                            scalar gasConstant = 0.0;
                            scalar gamma = 0.0;
                            scalar speedOfSound = 0.0;

                            if
                            (
                                Ttra_.boundaryFieldRef()[j][k] > SMALL
                             && molecularMassBF > VSMALL
                             && molarCvBF_trarot > SMALL
                             && molarCpBF_trarot > SMALL
                            )
                            {
                                gasConstant = kB/molecularMassBF;
                                gamma = molarCpBF_trarot/molarCvBF_trarot;

                                if (gamma > SMALL && gasConstant > SMALL)
                                {
                                    speedOfSound = sqrt
                                    (
                                        max
                                        (
                                            gamma*gasConstant*Ttra_.boundaryField()[j][k],
                                            scalar(0.0)
                                        )
                                    );
                                }
                            }

                            if (speedOfSound > SMALL)
                            {
                                Ma_.boundaryFieldRef()[j][k] =
                                    mag(UMean_.boundaryField()[j][k])/speedOfSound;
                            }
                            else
                            {
                                Ma_.boundaryFieldRef()[j][k] = 0.0;
                            }
                        }
                        else
                        {
                            Ma_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Force density
                        fD_.boundaryFieldRef()[j][k] = fDBF_[j][k]/nAvTimeSteps;

                        //- XCX: Total Force
                        //- XCX: boundaryFieldRef-used for modify, boundaryField-used for read
                        scalar faceArea = mag(patch.faceAreas()[k]);
                        vector faceForce = fD_.boundaryField()[j][k] * faceArea;
                        totalWallForce_ += faceForce;

                        //- Surface pressure
                        p_.boundaryFieldRef()[j][k] =
                            fD_.boundaryField()[j][k] & n_[j][k];

                        //- Pressure coefficient
                        Cp_.boundaryFieldRef()[j][k] = p_.boundaryFieldRef()[j][k] * cpCoef_;
                            
                        //- Wall shear stress
                        tau_.boundaryFieldRef()[j][k] =
                            sqrt
                            (
                                sqr(fD_.boundaryField()[j][k] & t1_[j][k])
                              + sqr(fD_.boundaryField()[j][k] & t2_[j][k])
                            );
                            
                        //- Heat flux
                        q_.boundaryFieldRef()[j][k] = qBF_[j][k]/nAvTimeSteps;

                        //- Pressure coefficient
                        Cq_.boundaryFieldRef()[j][k] = q_.boundaryField()[j][k] * cqCoef_;
                        
                        //- ZeroGradient condition assumed for Optional fields
                        if (measureMeanFreePath_)
                        {
                            mfp_.boundaryFieldRef()[j][k] = mfp_[celli];
                            SOF_.boundaryFieldRef()[j][k] = SOF_[celli];
                            mfpToDx_.boundaryFieldRef()[j][k] = mfpToDx_[celli];
                            meanCollisionRate_.boundaryFieldRef()[j][k] =
                                meanCollisionRate_[celli];
                            meanCollisionTime_.boundaryFieldRef()[j][k] =
                                meanCollisionTime_[celli];
                            mctToDt_.boundaryFieldRef()[j][k] = mctToDt_[celli];
                        }
                    }
                }
                else if (isNonEmptyNonCyclic)
                {
                    //- Loop over all boundary faces and set zeroGradient
                    //  conditions
                    forAll(boundaryCells_[j], k)
                    {
                        const label celli = boundaryCells_[j][k];

                        //- Instantaneous and sampled numbers of DSMC parcels
                        //  are that of the neighbouring cell
                        dsmcN_.boundaryFieldRef()[j][k] = dsmcN_[celli];
                        dsmcNMean_.boundaryFieldRef()[j][k] =
                            dsmcNMean_[celli];
                            
                        //- Number density and mass density fields
                        rhoN_.boundaryFieldRef()[j][k] = rhoN_[celli];
                        rhoM_.boundaryFieldRef()[j][k] = rhoM_[celli];
                        
                        //- Temperature fields
                        Ttra_.boundaryFieldRef()[j][k] = Ttra_[celli];
                        Trot_.boundaryFieldRef()[j][k] = Trot_[celli];
                        Tvib_.boundaryFieldRef()[j][k] = Tvib_[celli];
                        Tov_.boundaryFieldRef()[j][k] = Tov_[celli];
                        
                        //- Pressure, Mach and velocity fields
                        p_.boundaryFieldRef()[j][k] = p_[celli];
                        Ma_.boundaryFieldRef()[j][k] = Ma_[celli];
                        UMean_.boundaryFieldRef()[j][k] = UMean_[celli];
                        
                        //- Optional fields
                        if (measureMeanFreePath_)
                        {
                            mfp_.boundaryFieldRef()[j][k] = mfp_[celli];
                            SOF_.boundaryFieldRef()[j][k] = SOF_[celli];
                            mfpToDx_.boundaryFieldRef()[j][k] = mfpToDx_[celli];
                            meanCollisionRate_.boundaryFieldRef()[j][k] =
                                meanCollisionRate_[celli];
                            meanCollisionTime_.boundaryFieldRef()[j][k] =
                                meanCollisionTime_[celli];
                            mctToDt_.boundaryFieldRef()[j][k] = mctToDt_[celli];
                        }
                        
                        if (measureHeatFluxShearStress_)
                        {
                            shearStressTensor_.boundaryFieldRef()[j][k] =
                                shearStressTensor_[celli];
                            heatFluxVector_.boundaryFieldRef()[j][k] =
                                heatFluxVector_[celli];
                            pressureTensor_.boundaryFieldRef()[j][k] =
                                pressureTensor_[celli];
                        }
                        
                        if (measureClassifications_)
                        {
                            classIDistribution_.boundaryFieldRef()[j][k] =
                                classIDistribution_[celli];
                            classIIDistribution_.boundaryFieldRef()[j][k] =
                                classIIDistribution_[celli];
                            classIIIDistribution_.boundaryFieldRef()[j][k] =
                                classIIIDistribution_[celli];
                        }
                    }
                }
            }
            
            //- XCX: write force coefficient
            reduce(totalWallForce_, sumOp<vector>());
            if (Pstream::myProcNo() == 0)
            {
                OFstream forceFile
                (
                    "wallForce.dat"
                );
                forceFile << "wallForce_x " << totalWallForce_.x() << endl;
                forceFile << "wallForce_y " << totalWallForce_.y() << endl;
                forceFile << "wallForce_z " << totalWallForce_.z() << endl;

                Info << endl << "Total wall force: " << endl;
                Info << " wallForce_x " << totalWallForce_.x() << endl;
                Info << " wallForce_y " << totalWallForce_.y() << endl;
                Info << " wallForce_z " << totalWallForce_.z() << endl;
            }

            if (doProfile)
            {
                profileOutputComputeWallTime_ +=
                    wallSeconds(outputTimeStart, wallClockNow());
            }

            const auto fieldWriteStart =
                doProfile
              ? wallClockNow()
              : std::chrono::steady_clock::time_point();

            //- Write solution fields
            p_.write();
            Ttra_.write();
            UMean_.write();
            Ma_.write();
            q_.write();
            fD_.write();
            tau_.write();
            Cp_.write();
            Cq_.write();
            
            if (writeRotationalTemperature_)
            {
                Trot_.write();
            }
            if (writeVibrationalTemperature_)
            {
                Tvib_.write();
            }
            if (writeElectronicTemperature_)
            {
                Telec_.write();
            }
            if
            (
                  writeRotationalTemperature_ or writeVibrationalTemperature_
                or writeElectronicTemperature_
            )
            {
                Tov_.write();
            }
            
            if (measureMeanFreePath_)
            {
                mfp_.write();
                mfpToDx_.write();
                meanCollisionTime_.write();
                mctToDt_.write();
                SOF_.write();
            }

            if (measureClassifications_)
            {
                classIDistribution_.write();
                classIIDistribution_.write();
                classIIIDistribution_.write();
            }

            if (measureErrors_)
            {
                densityError_.write();
                velocityError_.write();
                temperatureError_.write();
                pressureError_.write();
            }

            if (measureHeatFluxShearStress_)
            {
                heatFluxVector_.write();
                pressureTensor_.write();
                shearStressTensor_.write();
            }

            if (doProfile)
            {
                profileFieldWriteWallTime_ +=
                    wallSeconds(fieldWriteStart, wallClockNow());
            }
        }

        const auto outputResetStart =
            doProfile
          ? wallClockNow()
          : std::chrono::steady_clock::time_point();

        //- Reset fields after printing the instantaneous solution ... or
        //  continue sampling
        // if (time_.resetFieldsAtOutput())
        if(isSample_)
        {
            nTimeSteps_ = 0.0;
            
            forAll(dsmcNCum_, celli)
            {
                dsmcNCum_[celli] = 0.0;
                dsmcMCum_[celli] = 0.0;
                dsmcLinearKECum_[celli] = 0.0;
                dsmcMomentumCum_[celli] = vector::zero;
                dsmcErotCum_[celli] = 0.0;
                dsmcZetaRotCum_[celli] = 0.0;
                dsmcNElecLvlCum_[celli] = 0.0,
                dsmcNClassICum_[celli] = 0.0;
                dsmcNClassIICum_[celli] = 0.0;
                dsmcNClassIIICum_[celli] = 0.0;
                collisionSeparation_[celli] = 0.0;
                dsmcNCollsCum_[celli] = 0.0;
                measuredCollisionRate_[celli] = 0.0;
                dsmcMuuCum_[celli] = 0.0;
                dsmcMuvCum_[celli] = 0.0;
                dsmcMuwCum_[celli] = 0.0;
                dsmcMvvCum_[celli] = 0.0;
                dsmcMvwCum_[celli] = 0.0;
                dsmcMwwCum_[celli] = 0.0;
                dsmcMccCum_[celli] = 0.0;
                dsmcMccuCum_[celli] = 0.0;
                dsmcMccvCum_[celli] = 0.0;
                dsmcMccwCum_[celli] = 0.0;
                dsmcEuCum_[celli] = 0.0;
                dsmcEvCum_[celli] = 0.0;
                dsmcEwCum_[celli] = 0.0;
                dsmcECum_[celli] = 0.0;
                zetaVib_[celli] = 0.0;
                nCum_[celli] = 0.0;
                mCum_[celli] = 0.0;
                momentumCum_[celli] = vector::zero;
                linearKECum_[celli] = 0.0;
            }

            forAll(speciesIds_, i)
            {
                forAll(speciesTvib_[i], celli)
                {
                    dsmcNSpeciesCum_[i][celli] = 0.0;
                    nSpeciesCum_[i][celli] = 0.0;
                    dsmcMccSpeciesCum_[i][celli] = 0.0;
                    speciesMfp_[i][celli] = 0.0;
                    speciesMcr_[i][celli] = 0.0;
                    speciesTvib_[i][celli] = 0.0;
                    dsmcSpeciesEelecCum_[i][celli] = 0.0;
                    dsmcNGrndElecLvlSpeciesCum_[i][celli] = 0.0;
                    dsmcN1stElecLvlSpeciesCum_[i][celli] = 0.0;
                }

                forAll(dsmcSpeciesEvibModCum_[i], mod)
                {
                    forAll(dsmcSpeciesEvibModCum_[i][mod], celli)
                    {
                        dsmcSpeciesEvibModCum_[i][mod][celli] = 0.0;
                    }
                }
            }

            //- Reset boundary information
            forAll(rhoNBF_, j)
            {
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                linearKEBF_[j] = 0.0;
                rhoNIntBF_[j] = 0.0;
                
                ErotBF_[j] = 0.0;
                zetaRotBF_[j] = 0.0;
                zetaVibBF_[j] = 0.0;
                rhoNElecBF_[j] = 0.0;
                
                qBF_[j] = 0.0;
                fDBF_[j] = vector::zero;
                momentumBF_[j] = vector::zero;
            }

            forAll(speciesIds_, i)
            {
                forAll(speciesTvibBF_[i], j)
                {
                    speciesRhoNBF_[i][j] = 0.0;
                    speciesMccBF_[i][j] = 0.0;
                    speciesEvibBF_[i][j] = 0.0;
                    speciesTvibBF_[i][j] = 0.0;
                    speciesZetaVibBF_[i][j] = 0.0;
                    speciesEelecBF_[i][j] = 0.0;
                }
                
                forAll(speciesEvibModBF_[i], mod)
                {
                    forAll(speciesEvibModBF_[i][mod], j)
                    {
                        speciesEvibModBF_[i][mod][j] = 0.0;
                    }
                }
            }
        }

        if (doProfile)
        {
            profileOutputResetWallTime_ +=
                wallSeconds(outputResetStart, wallClockNow());
        }

        if (averagingAcrossManyRuns_ && !time_.resetFieldsAtOutput())
        {
            const auto writeOutStart =
                doProfile
              ? wallClockNow()
              : std::chrono::steady_clock::time_point();

            writeOut();

            if (doProfile)
            {
                profileFieldWriteWallTime_ +=
                    wallSeconds(writeOutStart, wallClockNow());
            }
        }

        if (doProfile)
        {
            profileOutputTimeWallTime_ +=
                wallSeconds(outputTimeStart, wallClockNow());
        }
    }

    if
    (
        doProfile
     && Pstream::master()
     && !finalProfilePrinted_
     && time_.time().value() + time_.time().deltaT().value()
        >= time_.time().endTime().value() - SMALL
    )
    {
        finalProfilePrinted_ = true;
        Info<< "dsmcVolFields final profiling [" << fieldName_ << "]" << nl
            << "    calls                 = " << profileCalls_ << nl
            << "    sample accumulation   = " << profileSampleAccumWallTime_ << " s" << nl
            << "    shared cache build    = " << profileSharedCacheBuildWallTime_ << " s" << nl
            << "      cache allocate      = " << profileSharedCacheAllocateWallTime_ << " s" << nl
            << "      cache reset         = " << profileSharedCacheResetWallTime_ << " s" << nl
            << "      parcel accumulate   = " << profileSharedCacheParcelAccumWallTime_ << " s" << nl
            << "      detail sample cells = " << profileSharedCacheDetailSampleCells_ << nl
            << "      detail sample parcels = " << profileSharedCacheDetailSampleParcels_ << nl
            << "      base accum (sampled)= " << profileSharedCacheBaseAccumWallTime_ << " s" << nl
            << "      vib accum (sampled) = " << profileSharedCacheVibAccumWallTime_ << " s" << nl
            << "      electronic (sampled)= " << profileSharedCacheElectronicAccumWallTime_ << " s" << nl
            << "      class accum (sampled)= " << profileSharedCacheClassAccumWallTime_ << " s" << nl
            << "    field combine         = " << profileFieldCombineWallTime_ << " s" << nl
            << "    cell reduction        = " << profileCellReduceWallTime_ << " s" << nl
            << "    boundary accumulation = " << profileBoundaryAccumWallTime_ << " s" << nl
            << "    output compute        = " << profileOutputComputeWallTime_ << " s" << nl
            << "    field writes          = " << profileFieldWriteWallTime_ << " s" << nl
            << "    output reset          = " << profileOutputResetWallTime_ << " s" << nl
            << "    output-time block     = " << profileOutputTimeWallTime_ << " s" << nl
            << endl;
    }
}


//- reset fields when mesh is edited
void dsmcVolFields::resetField()
{
    const label nCells = mesh_.nCells();
    
    nTimeSteps_ = 0.0;

    //- Reset volume information
    dsmcNCum_.clear();
    dsmcMCum_.clear();
    dsmcLinearKECum_.clear();
    dsmcMomentumCum_.clear();
    dsmcErotCum_.clear();
    dsmcZetaRotCum_.clear();
    dsmcNElecLvlCum_.clear();
    dsmcNClassICum_.clear();
    dsmcNClassIICum_.clear();
    dsmcNClassIIICum_.clear();
    collisionSeparation_.clear();
    dsmcNCollsCum_.clear();
    dsmcMuuCum_.clear();
    dsmcMuvCum_.clear();
    dsmcMuwCum_.clear();
    dsmcMvvCum_.clear();
    dsmcMvwCum_.clear();
    dsmcMwwCum_.clear();
    dsmcMccCum_.clear();
    dsmcMccuCum_.clear();
    dsmcMccvCum_.clear();
    dsmcMccwCum_.clear();
    dsmcEuCum_.clear();
    dsmcEvCum_.clear();
    dsmcEwCum_.clear();
    dsmcECum_.clear();
    zetaVib_.clear();
    nCum_.clear();
    mCum_.clear();
    momentumCum_.clear();
    linearKECum_.clear();

    dsmcNCum_.setSize(nCells, 0.0);
    dsmcMCum_.setSize(nCells, 0.0);
    dsmcLinearKECum_.setSize(nCells, 0.0);
    dsmcMomentumCum_.setSize(nCells, vector::zero);
    dsmcErotCum_.setSize(nCells, 0.0);
    dsmcZetaRotCum_.setSize(nCells, 0.0);
    dsmcNElecLvlCum_.setSize(nCells, 0.0);
    dsmcNClassICum_.setSize(nCells, 0.0);
    dsmcNClassIICum_.setSize(nCells, 0.0);
    dsmcNClassIIICum_.setSize(nCells, 0.0);
    collisionSeparation_.setSize(nCells, 0.0);
    dsmcNCollsCum_.setSize(nCells, 0.0);
    measuredCollisionRate_.setSize(nCells, 0.0);
    dsmcMuuCum_.setSize(nCells, 0.0);
    dsmcMuvCum_.setSize(nCells, 0.0);
    dsmcMuwCum_.setSize(nCells, 0.0);
    dsmcMvvCum_.setSize(nCells, 0.0);
    dsmcMvwCum_.setSize(nCells, 0.0);
    dsmcMwwCum_.setSize(nCells, 0.0);
    dsmcMccCum_.setSize(nCells, 0.0);
    dsmcMccuCum_.setSize(nCells, 0.0);
    dsmcMccvCum_.setSize(nCells, 0.0);
    dsmcMccwCum_.setSize(nCells, 0.0);
    dsmcEuCum_.setSize(nCells, 0.0);
    dsmcEvCum_.setSize(nCells, 0.0);
    dsmcEwCum_.setSize(nCells, 0.0);
    dsmcECum_.setSize(nCells, 0.0);
    zetaVib_.setSize(nCells, 0.0);
    nCum_.setSize(nCells, 0.0);
    mCum_.setSize(nCells, 0.0);
    momentumCum_.setSize(nCells, vector::zero);
    linearKECum_.setSize(nCells, 0.0);

    forAll(speciesIds_, i)
    {
        dsmcNSpeciesCum_[i].clear();
        nSpeciesCum_[i].clear();
        dsmcMccSpeciesCum_[i].clear();
        speciesMfp_[i].clear();
        speciesMcr_[i].clear();
        speciesTvib_[i].clear();
        dsmcSpeciesEelecCum_[i].clear();
        dsmcNGrndElecLvlSpeciesCum_[i].clear();
        dsmcN1stElecLvlSpeciesCum_[i].clear();

        dsmcNSpeciesCum_[i].setSize(nCells, 0.0);
        nSpeciesCum_[i].setSize(nCells, 0.0);
        dsmcMccSpeciesCum_[i].setSize(nCells, 0.0);
        speciesMfp_[i].setSize(nCells, 0.0);
        speciesMcr_[i].setSize(nCells, 0.0);
        speciesTvib_[i].setSize(nCells, 0.0);
        dsmcSpeciesEelecCum_[i].setSize(nCells, 0.0);
        dsmcNGrndElecLvlSpeciesCum_[i].setSize(nCells, 0.0);
        dsmcN1stElecLvlSpeciesCum_[i].setSize(nCells, 0.0);
        
        forAll(dsmcSpeciesEvibModCum_[i], mod)
        {
           dsmcSpeciesEvibModCum_[i][mod].clear();
           dsmcSpeciesEvibModCum_[i][mod].setSize(nCells, 0.0);
        }
    }

    //- Reset boundary information
    forAll(mesh_.boundaryMesh(), j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        const label nFaces = patch.size();

        rhoNBF_[j].clear();
        rhoMBF_[j].clear();
        linearKEBF_[j].clear();
        momentumBF_[j].clear();
        ErotBF_[j].clear();
        zetaRotBF_[j].clear();
        qBF_[j].clear();
        fDBF_[j].clear();
        zetaVibBF_[j].clear();
        rhoNIntBF_[j].clear();
        rhoNElecBF_[j].clear();

        n_[j].clear();
        t1_[j].clear();
        t2_[j].clear();

        rhoNBF_[j].setSize(nFaces, 0.0);
        rhoMBF_[j].setSize(nFaces, 0.0);
        linearKEBF_[j].setSize(nFaces, 0.0);
        momentumBF_[j].setSize(nFaces, vector::zero);
        ErotBF_[j].setSize(nFaces, 0.0);
        zetaRotBF_[j].setSize(nFaces, 0.0);
        qBF_[j].setSize(nFaces, 0.0);
        fDBF_[j].setSize(nFaces, vector::zero);
        zetaVibBF_[j].setSize(nFaces, 0.0);
        rhoNIntBF_[j].setSize(nFaces, 0.0);
        rhoNElecBF_[j].setSize(nFaces, 0.0);

        n_[j].setSize(nFaces, vector::zero);
        t1_[j].setSize(nFaces, vector::zero);
        t2_[j].setSize(nFaces, vector::zero);
    }

    forAll(speciesIds_, i)
    {
        const label nPatches = mesh_.boundaryMesh().size();
        
        speciesEvibBF_[i].clear();
        speciesEelecBF_[i].clear();
        speciesRhoNBF_[i].clear();
        speciesMccBF_[i].clear();
        speciesTvibBF_[i].clear();
        speciesZetaVibBF_[i].clear();
        speciesEvibModBF_[i].clear();

        speciesRhoNBF_[i].setSize(nPatches);
        speciesMccBF_[i].setSize(nPatches);
        speciesEvibBF_[i].setSize(nPatches);
        speciesTvibBF_[i].setSize(nPatches);
        speciesZetaVibBF_[i].setSize(nPatches);
        speciesEelecBF_[i].setSize(nPatches);

        forAll(mesh_.boundaryMesh(), j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            const label nFaces = patch.size();

            speciesRhoNBF_[i][j].clear();
            speciesMccBF_[i][j].clear();
            speciesEvibBF_[i][j].clear();
            speciesTvibBF_[i][j].clear();
            speciesZetaVibBF_[i][j].clear();
            speciesEelecBF_[i][j].clear();

            speciesRhoNBF_[i][j].setSize(nFaces, 0.0);
            speciesMccBF_[i][j].setSize(nFaces, 0.0);
            speciesEvibBF_[i][j].setSize(nFaces, 0.0);
            speciesTvibBF_[i][j].setSize(nFaces, 0.0);
            speciesZetaVibBF_[i][j].setSize(nFaces, 0.0);
            speciesEelecBF_[i][j].setSize(nFaces, 0.0);
        }

        forAll(speciesEvibModBF_[i], mod)
        {
            const label nPatches = mesh_.boundaryMesh().size();
            speciesEvibModBF_[i][mod].setSize(nPatches);

            forAll(speciesEvibModBF_[i][mod], j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                const label nFaces = patch.size();
                speciesEvibModBF_[i][mod][j].setSize(nFaces, 0.0);
            }
        }
    }

    forAll(boundaryCells_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        const label nFaces = patch.size();

        boundaryCells_[j].clear();
        boundaryCells_[j].setSize(nFaces);

        forAll(boundaryCells_[j], k)
        {
            boundaryCells_[j][k] = patch.faceCells()[k];
        }
    }
}


void dsmcVolFields::writeField()
{}


void dsmcVolFields::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);
}

} // End namespace Foam

// ************************************************************************** //




