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
#include <mpi.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcVolFields, 0);

addToRunTimeSelectionTable(dsmcField, dsmcVolFields, dictionary);


namespace
{
    template<class Type>
    void sumReduceList(Foam::List<Type>& values)
    {
        Foam::Pstream::listCombineGather(values, Foam::plusEqOp<Type>());
        Foam::Pstream::listCombineScatter(values);
    }

    template<class Type>
    void sumReduceField(Foam::Field<Type>& values)
    {
        if (!values.size())
        {
            return;
        }

        int mpiInit = 0;
        MPI_Initialized(&mpiInit);
        if (!mpiInit)
        {
            return;
        }

        const int nScalars =
            values.size()*sizeof(Type)/sizeof(Foam::scalar);

        MPI_Allreduce
        (
            MPI_IN_PLACE,
            reinterpret_cast<Foam::scalar*>(values.begin()),
            nScalars,
            MPI_DOUBLE,
            MPI_SUM,
            MPI_COMM_WORLD
        );
    }

    template<class Type>
    void sumReduceFieldList(Foam::List<Foam::Field<Type>>& values)
    {
        forAll(values, i)
        {
            sumReduceField(values[i]);
        }
    }

    template<class Type>
    void sumReduceFieldListList
    (
        Foam::List<Foam::List<Foam::Field<Type>>>& values
    )
    {
        forAll(values, i)
        {
            sumReduceFieldList(values[i]);
        }
    }

    template<class Type>
    void sumReduceFieldListListList
    (
        Foam::List<Foam::List<Foam::List<Foam::Field<Type>>>>& values
    )
    {
        forAll(values, i)
        {
            sumReduceFieldListList(values[i]);
        }
    }

    struct dsmcVolSharedSampleCache
    {
        struct BuildProfile
        {
            scalar allocateWallTime;
            scalar resetWallTime;
            scalar parcelAccumWallTime;
            scalar baseAccumWallTime;
            scalar vibAccumWallTime;
            scalar electronicAccumWallTime;
            scalar classAccumWallTime;
            label detailSampleCells;
            label detailSampleParcels;

            BuildProfile()
            :
                allocateWallTime(0.0),
                resetWallTime(0.0),
                parcelAccumWallTime(0.0),
                baseAccumWallTime(0.0),
                vibAccumWallTime(0.0),
                electronicAccumWallTime(0.0),
                classAccumWallTime(0.0),
                detailSampleCells(0),
                detailSampleParcels(0)
            {}
        };

        const dsmcCloud* cloudPtr;
        scalar timeValue;
        label nCells;
        label nTypes;
        bool built;
        bool supportsVibrational;
        bool supportsElectronic;
        bool supportsClassification;
        bool supportsHeatFluxShearStress;

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
        labelList touchedCells;

        dsmcVolSharedSampleCache()
        :
            cloudPtr(nullptr),
            timeValue(-GREAT),
            nCells(0),
            nTypes(0),
            built(false),
            supportsVibrational(false),
            supportsElectronic(false),
            supportsClassification(false),
            supportsHeatFluxShearStress(false),
            touchedCells()
        {}

        bool validFor
        (
            const dsmcCloud& cloud,
            const scalar currentTime,
            const bool needVibrational,
            const bool needElectronic,
            const bool needClassification,
            const bool needHeatFluxShearStress
        ) const
        {
            return cloudPtr == &cloud
                && mag(timeValue - currentTime) < SMALL
                && nCells == cloud.mesh().nCells()
                && nTypes == cloud.constProps().size()
                && (!needVibrational || supportsVibrational)
                && (!needElectronic || supportsElectronic)
                && (!needClassification || supportsClassification)
                && (!needHeatFluxShearStress || supportsHeatFluxShearStress);
        }

        void allocateFields
        (
            const dsmcCloud& cloud,
            const bool needVibrational,
            const bool needElectronic,
            const bool needClassification,
            const bool needHeatFluxShearStress
        )
        {
            cloudPtr = &cloud;
            nCells = cloud.mesh().nCells();
            nTypes = cloud.constProps().size();
            built = false;
            supportsVibrational = needVibrational;
            supportsElectronic = needElectronic;
            supportsClassification = needClassification;
            supportsHeatFluxShearStress = needHeatFluxShearStress;
            touchedCells.clear();

            allocateAllFields(cloud);
        }

        void initScalarFields(List<scalarField>& fields)
        {
            fields.setSize(nTypes);
            for (label typei = 0; typei < nTypes; ++typei)
            {
                fields[typei].setSize(nCells, 0.0);
            }
        }

        void initVectorFields(List<vectorField>& fields)
        {
            fields.setSize(nTypes);
            for (label typei = 0; typei < nTypes; ++typei)
            {
                fields[typei].setSize(nCells, vector::zero);
            }
        }

        void allocateAllFields(const dsmcCloud& cloud)
        {
            initScalarFields(dsmcN);
            initScalarFields(dsmcM);
            initScalarFields(dsmcLinearKE);
            initVectorFields(dsmcMomentum);
            initScalarFields(dsmcErot);
            initScalarFields(dsmcZetaRot);
            initScalarFields(nReal);
            initScalarFields(mReal);
            initVectorFields(momentumReal);
            initScalarFields(linearKEReal);

            if (supportsElectronic)
            {
                initScalarFields(dsmcSpeciesEelec);
                initScalarFields(dsmcNElecLvl);
                initScalarFields(nGrndElecLvl);
                initScalarFields(n1stElecLvl);
            }
            else
            {
                dsmcSpeciesEelec.clear();
                dsmcNElecLvl.clear();
                nGrndElecLvl.clear();
                n1stElecLvl.clear();
            }

            if (supportsHeatFluxShearStress)
            {
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
            }
            else
            {
                dsmcMuu.clear();
                dsmcMuv.clear();
                dsmcMuw.clear();
                dsmcMvv.clear();
                dsmcMvw.clear();
                dsmcMww.clear();
                dsmcMcc.clear();
                dsmcMccu.clear();
                dsmcMccv.clear();
                dsmcMccw.clear();
                dsmcEu.clear();
                dsmcEv.clear();
                dsmcEw.clear();
                dsmcECum.clear();
            }

            if (supportsClassification)
            {
                initScalarFields(dsmcNClassI);
                initScalarFields(dsmcNClassII);
                initScalarFields(dsmcNClassIII);
            }
            else
            {
                dsmcNClassI.clear();
                dsmcNClassII.clear();
                dsmcNClassIII.clear();
            }

            if (supportsVibrational)
            {
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
            else
            {
                dsmcSpeciesEvibMod.clear();
            }
        }

        void resetFields()
        {
            forAll(dsmcN, typei)
            {
                if (!touchedCells.size())
                {
                    dsmcN[typei] = 0.0;
                    dsmcM[typei] = 0.0;
                    dsmcLinearKE[typei] = 0.0;
                    dsmcMomentum[typei] = vector::zero;
                    dsmcErot[typei] = 0.0;
                    dsmcZetaRot[typei] = 0.0;
                    nReal[typei] = 0.0;
                    mReal[typei] = 0.0;
                    momentumReal[typei] = vector::zero;
                    linearKEReal[typei] = 0.0;
                }
                else
                {
                    forAll(touchedCells, ci)
                    {
                        const label celli = touchedCells[ci];
                        dsmcN[typei][celli] = 0.0;
                        dsmcM[typei][celli] = 0.0;
                        dsmcLinearKE[typei][celli] = 0.0;
                        dsmcMomentum[typei][celli] = vector::zero;
                        dsmcErot[typei][celli] = 0.0;
                        dsmcZetaRot[typei][celli] = 0.0;
                        nReal[typei][celli] = 0.0;
                        mReal[typei][celli] = 0.0;
                        momentumReal[typei][celli] = vector::zero;
                        linearKEReal[typei][celli] = 0.0;
                    }
                }

                if (supportsElectronic)
                {
                    if (!touchedCells.size())
                    {
                        dsmcSpeciesEelec[typei] = 0.0;
                        dsmcNElecLvl[typei] = 0.0;
                        nGrndElecLvl[typei] = 0.0;
                        n1stElecLvl[typei] = 0.0;
                    }
                    else
                    {
                        forAll(touchedCells, ci)
                        {
                            const label celli = touchedCells[ci];
                            dsmcSpeciesEelec[typei][celli] = 0.0;
                            dsmcNElecLvl[typei][celli] = 0.0;
                            nGrndElecLvl[typei][celli] = 0.0;
                            n1stElecLvl[typei][celli] = 0.0;
                        }
                    }
                }

                if (supportsHeatFluxShearStress)
                {
                    if (!touchedCells.size())
                    {
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
                    }
                    else
                    {
                        forAll(touchedCells, ci)
                        {
                            const label celli = touchedCells[ci];
                            dsmcMuu[typei][celli] = 0.0;
                            dsmcMuv[typei][celli] = 0.0;
                            dsmcMuw[typei][celli] = 0.0;
                            dsmcMvv[typei][celli] = 0.0;
                            dsmcMvw[typei][celli] = 0.0;
                            dsmcMww[typei][celli] = 0.0;
                            dsmcMcc[typei][celli] = 0.0;
                            dsmcMccu[typei][celli] = 0.0;
                            dsmcMccv[typei][celli] = 0.0;
                            dsmcMccw[typei][celli] = 0.0;
                            dsmcEu[typei][celli] = 0.0;
                            dsmcEv[typei][celli] = 0.0;
                            dsmcEw[typei][celli] = 0.0;
                            dsmcECum[typei][celli] = 0.0;
                        }
                    }
                }

                if (supportsClassification)
                {
                    if (!touchedCells.size())
                    {
                        dsmcNClassI[typei] = 0.0;
                        dsmcNClassII[typei] = 0.0;
                        dsmcNClassIII[typei] = 0.0;
                    }
                    else
                    {
                        forAll(touchedCells, ci)
                        {
                            const label celli = touchedCells[ci];
                            dsmcNClassI[typei][celli] = 0.0;
                            dsmcNClassII[typei][celli] = 0.0;
                            dsmcNClassIII[typei][celli] = 0.0;
                        }
                    }
                }

                if (supportsVibrational)
                {
                    forAll(dsmcSpeciesEvibMod[typei], mod)
                    {
                        if (!touchedCells.size())
                        {
                            dsmcSpeciesEvibMod[typei][mod] = 0.0;
                        }
                        else
                        {
                            forAll(touchedCells, ci)
                            {
                                dsmcSpeciesEvibMod[typei][mod][touchedCells[ci]] =
                                    0.0;
                            }
                        }
                    }
                }
            }

            touchedCells.clear();
        }

        void build
        (
            const dsmcCloud& cloud,
            const List<DynamicList<dsmcParcel*>>* cellOccupancyPtr,
            const List<dsmcParcel*>* occupancyOrderedParcelsPtr,
            const labelList* occupancyCellOffsetsPtr,
            const scalar currentTime,
            const bool needVibrational,
            const bool needElectronic,
            const bool needClassification,
            const bool needHeatFluxShearStress,
            BuildProfile* buildProfile = nullptr,
            const bool doDetailProfile = false
        )
        {
            const bool doProfile = (buildProfile != nullptr);
            const bool doDetailedProfile = doProfile && doDetailProfile;
            const bool useFlatOccupancy =
                occupancyOrderedParcelsPtr
             && occupancyCellOffsetsPtr
             && occupancyCellOffsetsPtr->size() == cloud.mesh().nCells() + 1
             && occupancyOrderedParcelsPtr->size() == occupancyCellOffsetsPtr->last();
            const UList<label>* activeCellsPtr =
                cloud.replicatedMeshActive()
              ? static_cast<const UList<label>*>(&cloud.replicatedMesh().myCells())
              : (
                    useFlatOccupancy
                  ? static_cast<const UList<label>*>(&cloud.occupancyActiveCells())
                  : nullptr
                );

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

            if
            (
                validFor
                (
                    cloud,
                    currentTime,
                    needVibrational,
                    needElectronic,
                    needClassification,
                    needHeatFluxShearStress
                )
             && built
            )
            {
                return;
            }

            if
            (
                !validFor
                (
                    cloud,
                    currentTime,
                    needVibrational,
                    needElectronic,
                    needClassification,
                    needHeatFluxShearStress
                )
            )
            {
                const auto allocateStart =
                    doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
                allocateFields
                (
                    cloud,
                    needVibrational,
                    needElectronic,
                    needClassification,
                    needHeatFluxShearStress
                );
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
                resetFields();
                if (doProfile)
                {
                    buildProfile->resetWallTime +=
                        wallSeconds(resetStart, wallClockNow());
                }
            }

            timeValue = currentTime;

            if (activeCellsPtr)
            {
                touchedCells.setSize(activeCellsPtr->size());
                forAll(touchedCells, i)
                {
                    touchedCells[i] = (*activeCellsPtr)[i];
                }
            }
            else
            {
                touchedCells.setSize(nCells);
                for (label celli = 0; celli < nCells; ++celli)
                {
                    touchedCells[celli] = celli;
                }
            }

            const auto parcelAccumStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();

            auto accumulateIntoCache =
                [&]
                (
                    const label celli,
                    const dsmcParcel& p,
                    const scalar nParticles,
                    scalarField& localDsmcN,
                    scalarField& localDsmcM,
                    scalarField& localDsmcLinearKE,
                    vectorField& localDsmcMomentum,
                    scalarField& localDsmcErot,
                    scalarField& localDsmcZetaRot,
                    scalarField& localNReal,
                    scalarField& localMReal,
                    vectorField& localMomentumReal,
                    scalarField& localLinearKEReal,
                    scalarField& localSpeciesEelec,
                    scalarField& localNElecLvl,
                    scalarField& localNGrndElecLvl,
                    scalarField& localN1stElecLvl,
                    scalarField& localDsmcMuu,
                    scalarField& localDsmcMuv,
                    scalarField& localDsmcMuw,
                    scalarField& localDsmcMvv,
                    scalarField& localDsmcMvw,
                    scalarField& localDsmcMww,
                    scalarField& localDsmcMcc,
                    scalarField& localDsmcMccu,
                    scalarField& localDsmcMccv,
                    scalarField& localDsmcMccw,
                    scalarField& localDsmcEu,
                    scalarField& localDsmcEv,
                    scalarField& localDsmcEw,
                    scalarField& localDsmcECum,
                    scalarField& localDsmcNClassI,
                    scalarField& localDsmcNClassII,
                    scalarField& localDsmcNClassIII,
                    List<scalarField>& localSpeciesEvibMod,
                    DynamicList<label>& activeTypes
                )
            {
                if (!p.isFree())
                {
                    return;
                }

                const label typeId = p.typeId();
                if (typeId < 0 || typeId >= nTypes)
                {
                    return;
                }

                const dsmcParcel::constantProperties& cP =
                    cloud.constProps(typeId);

                if (localDsmcN[typeId] == 0.0)
                {
                    activeTypes.append(typeId);
                }

                const scalar mp = cP.mass();
                const vector& Up = p.U();
                const scalar Upx = Up.x();
                const scalar Upy = Up.y();
                const scalar Upz = Up.z();
                const vector mpUp = mp*Up;
                const scalar mpNParticles = mp*nParticles;
                const vector mpUpNParticles = mpUp*nParticles;
                const scalar linearKE = mp*(Up & Up);
                const scalar Erotp = p.ERot();
                const scalar zetaRotp = cP.rotationalDegreesOfFreedom();

                scalar Evibp = 0.0;
                if (needVibrational)
                {
                    const labelList& vibLevels = p.vibLevel();
                    const scalarList& thetaV = cP.thetaV();
                    scalarField& localEvibMods = localSpeciesEvibMod[typeId];

                    forAll(thetaV, mod)
                    {
                        const scalar EvibMod = cP.eVib_m(mod, vibLevels[mod]);
                        localEvibMods[mod] += EvibMod;
                        Evibp += EvibMod;
                    }
                }

                localDsmcN[typeId] += 1.0;
                localDsmcM[typeId] += mp;
                localDsmcLinearKE[typeId] += linearKE;
                localDsmcMomentum[typeId] += mpUp;
        localDsmcErot[typeId] += Erotp;
        localDsmcZetaRot[typeId] += zetaRotp;
        localNReal[typeId] += nParticles;
        localMReal[typeId] += mpNParticles;
        localMomentumReal[typeId] += mpUpNParticles;
        localLinearKEReal[typeId] += linearKE*nParticles;

                if (needHeatFluxShearStress)
                {
                    const scalar Eintp = Erotp + Evibp;

                    localDsmcMuu[typeId] += mp*sqr(Upx);
                    localDsmcMuv[typeId] += mp*Upx*Upy;
                    localDsmcMuw[typeId] += mp*Upx*Upz;
                    localDsmcMvv[typeId] += mp*sqr(Upy);
                    localDsmcMvw[typeId] += mp*Upy*Upz;
                    localDsmcMww[typeId] += mp*sqr(Upz);

                    localDsmcMcc[typeId] += linearKE;
                    localDsmcMccu[typeId] += linearKE*Upx;
                    localDsmcMccv[typeId] += linearKE*Upy;
                    localDsmcMccw[typeId] += linearKE*Upz;
                    localDsmcEu[typeId] += Eintp*Upx;
                    localDsmcEv[typeId] += Eintp*Upy;
                    localDsmcEw[typeId] += Eintp*Upz;
                    localDsmcECum[typeId] += Eintp;
                }

                if (needElectronic)
                {
                    const label eLevel = p.ELevel();
                    localSpeciesEelec[typeId] += cP.electronicEnergyList()[eLevel];

                    if (cP.nElectronicLevels() > 1)
                    {
                        localNElecLvl[typeId] += 1.0;

                        if (eLevel == 0)
                        {
                            localNGrndElecLvl[typeId] += 1.0;
                        }
                        if (eLevel == 1)
                        {
                            localN1stElecLvl[typeId] += 1.0;
                        }
                    }
                }

                if (needClassification)
                {
                    const label classification = p.classification();
                    if (classification == 0)
                    {
                        localDsmcNClassI[typeId] += 1.0;
                    }
                    else if (classification == 1)
                    {
                        localDsmcNClassII[typeId] += 1.0;
                    }
                    else if (classification == 2)
                    {
                        localDsmcNClassIII[typeId] += 1.0;
                    }
                }
            };

            auto flushCell =
                [&]
                (
                    const label celli,
                    scalarField& localDsmcN,
                    scalarField& localDsmcM,
                    scalarField& localDsmcLinearKE,
                    vectorField& localDsmcMomentum,
                    scalarField& localDsmcErot,
                    scalarField& localDsmcZetaRot,
                    scalarField& localNReal,
                    scalarField& localMReal,
                    vectorField& localMomentumReal,
                    scalarField& localLinearKEReal,
                    scalarField& localSpeciesEelec,
                    scalarField& localNElecLvl,
                    scalarField& localNGrndElecLvl,
                    scalarField& localN1stElecLvl,
                    scalarField& localDsmcMuu,
                    scalarField& localDsmcMuv,
                    scalarField& localDsmcMuw,
                    scalarField& localDsmcMvv,
                    scalarField& localDsmcMvw,
                    scalarField& localDsmcMww,
                    scalarField& localDsmcMcc,
                    scalarField& localDsmcMccu,
                    scalarField& localDsmcMccv,
                    scalarField& localDsmcMccw,
                    scalarField& localDsmcEu,
                    scalarField& localDsmcEv,
                    scalarField& localDsmcEw,
                    scalarField& localDsmcECum,
                    scalarField& localDsmcNClassI,
                    scalarField& localDsmcNClassII,
                    scalarField& localDsmcNClassIII,
                    List<scalarField>& localSpeciesEvibMod,
                    DynamicList<label>& activeTypes
                )
            {
                forAll(activeTypes, activeI)
                {
                    const label typeId = activeTypes[activeI];
                    dsmcN[typeId][celli] = localDsmcN[typeId];
                    dsmcM[typeId][celli] = localDsmcM[typeId];
                    dsmcLinearKE[typeId][celli] = localDsmcLinearKE[typeId];
                    dsmcMomentum[typeId][celli] = localDsmcMomentum[typeId];
                    dsmcErot[typeId][celli] = localDsmcErot[typeId];
                    dsmcZetaRot[typeId][celli] = localDsmcZetaRot[typeId];
                    nReal[typeId][celli] = localNReal[typeId];
                    mReal[typeId][celli] = localMReal[typeId];
                    momentumReal[typeId][celli] = localMomentumReal[typeId];
                    linearKEReal[typeId][celli] = localLinearKEReal[typeId];

                    if (needElectronic)
                    {
                        dsmcSpeciesEelec[typeId][celli] =
                            localSpeciesEelec[typeId];
                        dsmcNElecLvl[typeId][celli] = localNElecLvl[typeId];
                        nGrndElecLvl[typeId][celli] =
                            localNGrndElecLvl[typeId];
                        n1stElecLvl[typeId][celli] =
                            localN1stElecLvl[typeId];
                    }

                    if (needHeatFluxShearStress)
                    {
                        dsmcMuu[typeId][celli] = localDsmcMuu[typeId];
                        dsmcMuv[typeId][celli] = localDsmcMuv[typeId];
                        dsmcMuw[typeId][celli] = localDsmcMuw[typeId];
                        dsmcMvv[typeId][celli] = localDsmcMvv[typeId];
                        dsmcMvw[typeId][celli] = localDsmcMvw[typeId];
                        dsmcMww[typeId][celli] = localDsmcMww[typeId];
                        dsmcMcc[typeId][celli] = localDsmcMcc[typeId];
                        dsmcMccu[typeId][celli] = localDsmcMccu[typeId];
                        dsmcMccv[typeId][celli] = localDsmcMccv[typeId];
                        dsmcMccw[typeId][celli] = localDsmcMccw[typeId];
                        dsmcEu[typeId][celli] = localDsmcEu[typeId];
                        dsmcEv[typeId][celli] = localDsmcEv[typeId];
                        dsmcEw[typeId][celli] = localDsmcEw[typeId];
                        dsmcECum[typeId][celli] = localDsmcECum[typeId];
                    }

                    if (needClassification)
                    {
                        dsmcNClassI[typeId][celli] =
                            localDsmcNClassI[typeId];
                        dsmcNClassII[typeId][celli] =
                            localDsmcNClassII[typeId];
                        dsmcNClassIII[typeId][celli] =
                            localDsmcNClassIII[typeId];
                    }

                    if (needVibrational)
                    {
                        scalarField& localEvibMods =
                            localSpeciesEvibMod[typeId];
                        forAll(localEvibMods, mod)
                        {
                            dsmcSpeciesEvibMod[typeId][mod][celli] =
                                localEvibMods[mod];
                            localEvibMods[mod] = 0.0;
                        }
                    }

                    localDsmcN[typeId] = 0.0;
                    localDsmcM[typeId] = 0.0;
                    localDsmcLinearKE[typeId] = 0.0;
                    localDsmcMomentum[typeId] = vector::zero;
                    localDsmcErot[typeId] = 0.0;
                    localDsmcZetaRot[typeId] = 0.0;
                    localNReal[typeId] = 0.0;
                    localMReal[typeId] = 0.0;
                    localMomentumReal[typeId] = vector::zero;
                    localLinearKEReal[typeId] = 0.0;

                    if (needElectronic)
                    {
                        localSpeciesEelec[typeId] = 0.0;
                        localNElecLvl[typeId] = 0.0;
                        localNGrndElecLvl[typeId] = 0.0;
                        localN1stElecLvl[typeId] = 0.0;
                    }

                    if (needHeatFluxShearStress)
                    {
                        localDsmcMuu[typeId] = 0.0;
                        localDsmcMuv[typeId] = 0.0;
                        localDsmcMuw[typeId] = 0.0;
                        localDsmcMvv[typeId] = 0.0;
                        localDsmcMvw[typeId] = 0.0;
                        localDsmcMww[typeId] = 0.0;
                        localDsmcMcc[typeId] = 0.0;
                        localDsmcMccu[typeId] = 0.0;
                        localDsmcMccv[typeId] = 0.0;
                        localDsmcMccw[typeId] = 0.0;
                        localDsmcEu[typeId] = 0.0;
                        localDsmcEv[typeId] = 0.0;
                        localDsmcEw[typeId] = 0.0;
                        localDsmcECum[typeId] = 0.0;
                    }

                    if (needClassification)
                    {
                        localDsmcNClassI[typeId] = 0.0;
                        localDsmcNClassII[typeId] = 0.0;
                        localDsmcNClassIII[typeId] = 0.0;
                    }
                }

                activeTypes.clear();
            };

            #ifdef _OPENMP
            if (cloud.openmpEnabled())
            {
                #pragma omp parallel
                {
                    scalarField localDsmcN(nTypes, 0.0);
                    scalarField localDsmcM(nTypes, 0.0);
                    scalarField localDsmcLinearKE(nTypes, 0.0);
                    vectorField localDsmcMomentum(nTypes, vector::zero);
                    scalarField localDsmcErot(nTypes, 0.0);
                    scalarField localDsmcZetaRot(nTypes, 0.0);
                    scalarField localNReal(nTypes, 0.0);
                    scalarField localMReal(nTypes, 0.0);
                    vectorField localMomentumReal(nTypes, vector::zero);
                    scalarField localLinearKEReal(nTypes, 0.0);
                    scalarField localSpeciesEelec;
                    scalarField localNElecLvl;
                    scalarField localNGrndElecLvl;
                    scalarField localN1stElecLvl;
                    scalarField localDsmcMuu;
                    scalarField localDsmcMuv;
                    scalarField localDsmcMuw;
                    scalarField localDsmcMvv;
                    scalarField localDsmcMvw;
                    scalarField localDsmcMww;
                    scalarField localDsmcMcc;
                    scalarField localDsmcMccu;
                    scalarField localDsmcMccv;
                    scalarField localDsmcMccw;
                    scalarField localDsmcEu;
                    scalarField localDsmcEv;
                    scalarField localDsmcEw;
                    scalarField localDsmcECum;
                    scalarField localDsmcNClassI;
                    scalarField localDsmcNClassII;
                    scalarField localDsmcNClassIII;
                    List<scalarField> localSpeciesEvibMod;
                    DynamicList<label> activeTypes;

                    if (needElectronic)
                    {
                        localSpeciesEelec.setSize(nTypes, 0.0);
                        localNElecLvl.setSize(nTypes, 0.0);
                        localNGrndElecLvl.setSize(nTypes, 0.0);
                        localN1stElecLvl.setSize(nTypes, 0.0);
                    }

                    if (needHeatFluxShearStress)
                    {
                        localDsmcMuu.setSize(nTypes, 0.0);
                        localDsmcMuv.setSize(nTypes, 0.0);
                        localDsmcMuw.setSize(nTypes, 0.0);
                        localDsmcMvv.setSize(nTypes, 0.0);
                        localDsmcMvw.setSize(nTypes, 0.0);
                        localDsmcMww.setSize(nTypes, 0.0);
                        localDsmcMcc.setSize(nTypes, 0.0);
                        localDsmcMccu.setSize(nTypes, 0.0);
                        localDsmcMccv.setSize(nTypes, 0.0);
                        localDsmcMccw.setSize(nTypes, 0.0);
                        localDsmcEu.setSize(nTypes, 0.0);
                        localDsmcEv.setSize(nTypes, 0.0);
                        localDsmcEw.setSize(nTypes, 0.0);
                        localDsmcECum.setSize(nTypes, 0.0);
                    }

                    if (needClassification)
                    {
                        localDsmcNClassI.setSize(nTypes, 0.0);
                        localDsmcNClassII.setSize(nTypes, 0.0);
                        localDsmcNClassIII.setSize(nTypes, 0.0);
                    }

                    if (needVibrational)
                    {
                        localSpeciesEvibMod.setSize(nTypes);
                        for (label typei = 0; typei < nTypes; ++typei)
                        {
                            localSpeciesEvibMod[typei].setSize
                            (
                                cloud.constProps(typei).thetaV().size(),
                                0.0
                            );
                        }
                    }

                    activeTypes.setCapacity(min(nTypes, label(64)));

                    const label loopSize =
                        activeCellsPtr ? activeCellsPtr->size() : nCells;

                    #pragma omp for schedule(static)
                    for (label activeI = 0; activeI < loopSize; ++activeI)
                    {
                        const label celli =
                            activeCellsPtr ? (*activeCellsPtr)[activeI] : activeI;
                        const DynamicList<dsmcParcel*>* parcelsPtr =
                            useFlatOccupancy ? nullptr : &(*cellOccupancyPtr)[celli];
                        const label parcelBegin =
                            useFlatOccupancy ? (*occupancyCellOffsetsPtr)[celli] : 0;
                        const label parcelEnd =
                            useFlatOccupancy
                          ? (*occupancyCellOffsetsPtr)[celli + 1]
                          : parcelsPtr->size();
                        const label parcelCount =
                            useFlatOccupancy ? parcelEnd - parcelBegin : parcelsPtr->size();
                        const scalar nParticles = cloud.nParticles(celli);

                        if (doDetailedProfile && parcelCount > 0 && (celli % 32 == 0))
                        {
                            #pragma omp atomic
                            buildProfile->detailSampleCells++;
                            #pragma omp atomic
                            buildProfile->detailSampleParcels += parcelCount;
                        }

                        for (label pi = 0; pi < parcelCount; ++pi)
                        {
                            const dsmcParcel& p =
                                useFlatOccupancy
                              ? *(*occupancyOrderedParcelsPtr)[parcelBegin + pi]
                              : *(*parcelsPtr)[pi];

                            if (useFlatOccupancy && pi + 1 < parcelCount)
                            {
                                __builtin_prefetch
                                (
                                    (*occupancyOrderedParcelsPtr)[parcelBegin + pi + 1],
                                    0,
                                    1
                                );
                            }

                            accumulateIntoCache
                            (
                                celli,
                                p,
                                nParticles,
                                localDsmcN,
                                localDsmcM,
                                localDsmcLinearKE,
                                localDsmcMomentum,
                                localDsmcErot,
                                localDsmcZetaRot,
                                localNReal,
                                localMReal,
                                localMomentumReal,
                                localLinearKEReal,
                                localSpeciesEelec,
                                localNElecLvl,
                                localNGrndElecLvl,
                                localN1stElecLvl,
                                localDsmcMuu,
                                localDsmcMuv,
                                localDsmcMuw,
                                localDsmcMvv,
                                localDsmcMvw,
                                localDsmcMww,
                                localDsmcMcc,
                                localDsmcMccu,
                                localDsmcMccv,
                                localDsmcMccw,
                                localDsmcEu,
                                localDsmcEv,
                                localDsmcEw,
                                localDsmcECum,
                                localDsmcNClassI,
                                localDsmcNClassII,
                                localDsmcNClassIII,
                                localSpeciesEvibMod,
                                activeTypes
                            );
                        }

                        flushCell
                        (
                            celli,
                            localDsmcN,
                            localDsmcM,
                            localDsmcLinearKE,
                            localDsmcMomentum,
                            localDsmcErot,
                            localDsmcZetaRot,
                            localNReal,
                            localMReal,
                            localMomentumReal,
                            localLinearKEReal,
                            localSpeciesEelec,
                            localNElecLvl,
                            localNGrndElecLvl,
                            localN1stElecLvl,
                            localDsmcMuu,
                            localDsmcMuv,
                            localDsmcMuw,
                            localDsmcMvv,
                            localDsmcMvw,
                            localDsmcMww,
                            localDsmcMcc,
                            localDsmcMccu,
                            localDsmcMccv,
                            localDsmcMccw,
                            localDsmcEu,
                            localDsmcEv,
                            localDsmcEw,
                            localDsmcECum,
                            localDsmcNClassI,
                            localDsmcNClassII,
                            localDsmcNClassIII,
                            localSpeciesEvibMod,
                            activeTypes
                        );
                    }
                }
            }
            else
            #endif
            {
                forAllConstIter(dsmcCloud, cloud, iter)
                {
                    const dsmcParcel& p = iter();
                    if (!p.isFree())
                    {
                        continue;
                    }

                    const label celli = p.cell();
                    const label typeId = p.typeId();
                    if (celli < 0 || celli >= nCells || typeId < 0 || typeId >= nTypes)
                    {
                        continue;
                    }

                    const dsmcParcel::constantProperties& cP =
                        cloud.constProps(typeId);
                    const scalar nParticles = cloud.nParticles(celli);
                    const scalar mp = cP.mass();
                    const vector& Up = p.U();
                    const scalar linearKE = mp*(Up & Up);
                    const scalar Erotp = p.ERot();
                    const scalar zetaRotp = cP.rotationalDegreesOfFreedom();

                    scalar Evibp = 0.0;
                    if (needVibrational)
                    {
                        const labelList& vibLevels = p.vibLevel();
                        forAll(cP.thetaV(), mod)
                        {
                            const scalar EvibMod = cP.eVib_m(mod, vibLevels[mod]);
                            dsmcSpeciesEvibMod[typeId][mod][celli] += EvibMod;
                            Evibp += EvibMod;
                        }
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
                        const label eLevel = p.ELevel();
                        dsmcSpeciesEelec[typeId][celli] +=
                            cP.electronicEnergyList()[eLevel];
                        if (cP.nElectronicLevels() > 1)
                        {
                            dsmcNElecLvl[typeId][celli] += 1.0;
                            if (eLevel == 0)
                            {
                                nGrndElecLvl[typeId][celli] += 1.0;
                            }
                            if (eLevel == 1)
                            {
                                n1stElecLvl[typeId][celli] += 1.0;
                            }
                        }
                    }

                    if (needClassification)
                    {
                        const label classification = p.classification();
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
                }
            }

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
    mfpTref_(273.0),
    nMinParcelsTvib_(1),
    iMeanMinTvib_(0.0),
    fieldName_(propsDict_.lookup("fieldName")),
    speciesIds_(),
    typeIdToSpeciesIndex_(),
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
    sampledBoundaryPatches_(),
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
    openmpFieldSampling_(true),
    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false),
    densityOnly_(false),
    measureHeatFluxShearStress_(false),
    writeRotationalTemperature_(false),
    writeVibrationalTemperature_(false),
    writeElectronicTemperature_(false),
    profileSummaryEnabled_(false),
    profileDetailEnabled_(false),
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
{
    const dictionary& controlDict = mesh_.time().controlDict();
    const bool globalProfileSummary =
        controlDict.lookupOrDefault<bool>("profileSummary", true);
    const bool globalProfileDetail =
        controlDict.lookupOrDefault<bool>("profileDetail", false);

    profileSummaryEnabled_ =
        controlDict.lookupOrDefault<bool>
        (
            "profilePostSummary",
            globalProfileSummary
        );
    profileDetailEnabled_ =
        controlDict.lookupOrDefault<bool>
        (
            "profilePostDetail",
            globalProfileDetail
        );
}


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

        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}


//- Initial configuration
void dsmcVolFields::createField()
{
    Info << "Initialising dsmcVolFields field" << endl;

    const List<word>& species (propsDict_.lookup("typeIds"));

    DynamicList<word> speciesReduced(0);

    forAll(species, i)
    {
        const word& speciesName(species[i]);

        if (findIndex(speciesReduced, speciesName) == -1)
        {
            speciesReduced.append(speciesName);
        }
    }

    speciesReduced.shrink();

    speciesIds_.setSize(speciesReduced.size(), -1);

    forAll(speciesReduced, i)
    {
        const word& speciesName = speciesReduced[i];

        const label spId = findIndex(cloud_.typeIdList(), speciesName);

        if (spId == -1)
        {
            FatalErrorIn("dsmcVolFields::dsmcVolFields()")
                << "Cannot find typeId: " << speciesName << nl << "in: "
                << mesh_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        speciesIds_[i] = spId;
    }

    typeIdToSpeciesIndex_.setSize(cloud_.typeIdList().size(), -1);
    forAll(speciesIds_, i)
    {
        const label typeId = speciesIds_[i];

        if (typeId >= 0 && typeId < typeIdToSpeciesIndex_.size())
        {
            typeIdToSpeciesIndex_[typeId] = i;
        }
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

    DynamicList<label> sampledPatchIds(nPatches);
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchi];
        const word& patchType = patch.type();

        if
        (
            patch.size()
         && patchType != "processor"
         && patchType != "empty"
         && patchType != "symmetry"
         && patchType != "wedge"
        )
        {
            sampledPatchIds.append(patchi);
        }
    }

    sampledPatchIds.shrink();
    sampledBoundaryPatches_.setSize(sampledPatchIds.size());
    forAll(sampledPatchIds, i)
    {
        sampledBoundaryPatches_[i] = sampledPatchIds[i];
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

    openmpFieldSampling_ =
        mesh_.time().controlDict().lookupOrDefault<bool>
        (
            "openmpFieldSampling",
            true
        );

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

    nMinParcelsTvib_ =
        mesh_.time().controlDict().lookupOrDefault<label>("nMinParcelsTvib", 1);

    iMeanMinTvib_ =
        mesh_.time().controlDict().lookupOrDefault<scalar>("iMeanMinTvib", 0.0);

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

    if (cloud_.replicatedMeshActive())
    {
        ownedBoundaryFaces_.setSize(mesh_.boundaryMesh().size());

        forAll(mesh_.boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            const label startFace = pp.start();
            DynamicList<label> owned(pp.size()/8 + 1);

            forAll(pp, facei)
            {
                if
                (
                    cloud_.replicatedMesh().isMyCell
                    (
                        mesh_.faceOwner()[startFace + facei]
                    )
                )
                {
                    owned.append(facei);
                }
            }

            ownedBoundaryFaces_[patchi].transfer(owned);
        }
    }
}


void dsmcVolFields::calculateField()
{
    sampleCounter_++;

    const scalar kB = physicoChemical::k.value();
    const scalar NAvo = physicoChemical::NA.value();
    const bool doProfile = profileDetailEnabled_;

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

    ++profileCalls_;
    
    //- Reset instantaneous number of DSMC parcels
    dsmcN_ = 0.0;

    if (sampleInterval_ <= sampleCounter_)
    {
        nTimeSteps_ += 1.0;
        const scalar nAvTimeSteps = nTimeSteps_;
        const auto sampleAccumStart =
            doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
        const bool useOpenMPSampling =
            openmpFieldSampling_
         && cloud_.openmpEnabled()
         && cloud_.hasOccupancyOrderedParcels();
        const label nMappedTypes = typeIdToSpeciesIndex_.size();

        if (densityOnly_)
        {
            const UList<label>* sampleCellsPtr =
                cloud_.replicatedMeshActive()
              ? static_cast<const UList<label>*>(&cloud_.replicatedMesh().myCells())
              : nullptr;
            const label sampleLoopSize =
                sampleCellsPtr ? sampleCellsPtr->size() : mesh_.nCells();

            if (useOpenMPSampling)
            {
                #ifdef _OPENMP
                #pragma omp parallel for schedule(static)
                #endif
                for (label sampleI = 0; sampleI < sampleLoopSize; ++sampleI)
                {
                    const label cell =
                        sampleCellsPtr ? (*sampleCellsPtr)[sampleI] : sampleI;
                    const label occStart = cloud_.occupancyStart(cell);
                    const label occEnd = cloud_.occupancyEnd(cell);

                    if (occStart == occEnd)
                    {
                        continue;
                    }

                    scalar dsmcNLocal = 0.0;
                    scalar nLocal = 0.0;
                    scalar mLocal = 0.0;

                    for (label occI = occStart; occI < occEnd; ++occI)
                    {
                        const dsmcParcel& p = *cloud_.occupancyParcel(occI);
                        const label typeId = p.typeId();
                        const label spId =
                        (
                            typeId >= 0 && typeId < nMappedTypes
                          ? typeIdToSpeciesIndex_[typeId]
                          : -1
                        );

                        //- Do not consider adsorbed parcels
                        if (spId != -1 && p.isFree())
                        {
                            const scalar nParticles = cloud_.nParticles(cell);
                            const scalar mass = cloud_.constProps(typeId).mass();

                            dsmcNLocal += 1.0;
                            nLocal += nParticles;
                            mLocal += mass*nParticles;
                        }
                    }

                    dsmcNCum_[cell] += dsmcNLocal;
                    dsmcN_[cell] += dsmcNLocal;
                    nCum_[cell] += nLocal;
                    mCum_[cell] += mLocal;
                }
            }
            else
            {
                //- Loop over the the entire parcel cloud
                forAllConstIter(dsmcCloud, cloud_, iter)
                {
                    const dsmcParcel& p = iter();
                    const label typeId = p.typeId();
                    const label spId =
                    (
                        typeId >= 0 && typeId < nMappedTypes
                      ? typeIdToSpeciesIndex_[typeId]
                      : -1
                    );

                    //- Do not consider adsorbed parcels
                    if (spId != -1 && p.isFree())
                    {
                        const label cell = p.cell();
                        if
                        (
                            cloud_.replicatedMeshActive()
                         && !cloud_.replicatedMesh().isMyCell(cell)
                        )
                        {
                            continue;
                        }
                        const scalar nParticles = cloud_.nParticles(cell);
                        const scalar mass = cloud_.constProps(typeId).mass();

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
        }
        else
        {
            bool needVibrational = false;
            forAll(speciesIds_, i)
            {
                if (cloud_.constProps(speciesIds_[i]).thetaV().size() > 0)
                {
                    needVibrational = true;
                    break;
                }
            }

            const bool needElectronic = false;
            const bool needClassification = measureClassifications_;
            const bool needHeatFluxShearStress = measureHeatFluxShearStress_;
            const auto sharedCacheBuildStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            dsmcVolSharedSampleCache::BuildProfile sharedCacheBuildProfile;

            sharedSampleCache_.build
            (
                cloud_,
                cloud_.hasOccupancyOrderedParcels() ? nullptr : &cloud_.cellOccupancy(),
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
                doProfile ? &sharedCacheBuildProfile : nullptr,
                profileDetailEnabled_
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
            const bool singleSpeciesField = (speciesIds_.size() == 1);
            const label onlyTypeId = singleSpeciesField ? speciesIds_[0] : -1;
            const UList<label>* combineCellsPtr =
                cloud_.replicatedMeshActive()
              ? static_cast<const UList<label>*>(&cloud_.replicatedMesh().myCells())
              : (
                    cloud_.hasOccupancyOrderedParcels()
                  ? static_cast<const UList<label>*>(&cloud_.occupancyActiveCells())
                  : nullptr
                );
            const label combineLoopSize =
                combineCellsPtr ? combineCellsPtr->size() : dsmcNCum_.size();

            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) if (useOpenMPSampling)
            #endif
            for (label combineI = 0; combineI < combineLoopSize; ++combineI)
            {
                const label cell =
                    combineCellsPtr ? (*combineCellsPtr)[combineI] : combineI;
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
                scalar dsmcECumLocal = 0.0;
                scalar dsmcNClassILocal = 0.0;
                scalar dsmcNClassIILocal = 0.0;
                scalar dsmcNClassIIILocal = 0.0;

                if (singleSpeciesField)
                {
                    dsmcNLocal = sharedSampleCache_.dsmcN[onlyTypeId][cell];
                    dsmcMLocal = sharedSampleCache_.dsmcM[onlyTypeId][cell];
                    dsmcLinearKELocal =
                        sharedSampleCache_.dsmcLinearKE[onlyTypeId][cell];
                    dsmcMomentumLocal =
                        sharedSampleCache_.dsmcMomentum[onlyTypeId][cell];
                    dsmcErotLocal = sharedSampleCache_.dsmcErot[onlyTypeId][cell];
                    dsmcZetaRotLocal =
                        sharedSampleCache_.dsmcZetaRot[onlyTypeId][cell];
                    nLocal = sharedSampleCache_.nReal[onlyTypeId][cell];
                    mLocal = sharedSampleCache_.mReal[onlyTypeId][cell];
                    momentumLocal =
                        sharedSampleCache_.momentumReal[onlyTypeId][cell];
                    linearKELocal =
                        sharedSampleCache_.linearKEReal[onlyTypeId][cell];

                    dsmcNSpeciesCum_[0][cell] += dsmcNLocal;
                    dsmcMccSpeciesCum_[0][cell] += dsmcLinearKELocal;
                    nSpeciesCum_[0][cell] += nLocal;

                    if (needVibrational)
                    {
                        forAll(dsmcSpeciesEvibModCum_[0], mod)
                        {
                            dsmcSpeciesEvibModCum_[0][mod][cell] +=
                                sharedSampleCache_.dsmcSpeciesEvibMod[onlyTypeId][mod][cell];
                        }
                    }

                    if (needHeatFluxShearStress)
                    {
                        dsmcMuuLocal = sharedSampleCache_.dsmcMuu[onlyTypeId][cell];
                        dsmcMuvLocal = sharedSampleCache_.dsmcMuv[onlyTypeId][cell];
                        dsmcMuwLocal = sharedSampleCache_.dsmcMuw[onlyTypeId][cell];
                        dsmcMvvLocal = sharedSampleCache_.dsmcMvv[onlyTypeId][cell];
                        dsmcMvwLocal = sharedSampleCache_.dsmcMvw[onlyTypeId][cell];
                        dsmcMwwLocal = sharedSampleCache_.dsmcMww[onlyTypeId][cell];
                        dsmcMccLocal = sharedSampleCache_.dsmcMcc[onlyTypeId][cell];
                        dsmcMccuLocal = sharedSampleCache_.dsmcMccu[onlyTypeId][cell];
                        dsmcMccvLocal = sharedSampleCache_.dsmcMccv[onlyTypeId][cell];
                        dsmcMccwLocal = sharedSampleCache_.dsmcMccw[onlyTypeId][cell];
                        dsmcEuLocal = sharedSampleCache_.dsmcEu[onlyTypeId][cell];
                        dsmcEvLocal = sharedSampleCache_.dsmcEv[onlyTypeId][cell];
                        dsmcEwLocal = sharedSampleCache_.dsmcEw[onlyTypeId][cell];
                        dsmcECumLocal = sharedSampleCache_.dsmcECum[onlyTypeId][cell];
                    }

                    if (needClassification)
                    {
                        dsmcNClassILocal =
                            sharedSampleCache_.dsmcNClassI[onlyTypeId][cell];
                        dsmcNClassIILocal =
                            sharedSampleCache_.dsmcNClassII[onlyTypeId][cell];
                        dsmcNClassIIILocal =
                            sharedSampleCache_.dsmcNClassIII[onlyTypeId][cell];
                    }
                }
                else
                {
                    forAll(speciesIds_, i)
                    {
                        const label typeId = speciesIds_[i];

                        const scalar speciesDsmcN =
                            sharedSampleCache_.dsmcN[typeId][cell];
                        const scalar speciesDsmcLinearKE =
                            sharedSampleCache_.dsmcLinearKE[typeId][cell];
                        const scalar speciesNReal =
                            sharedSampleCache_.nReal[typeId][cell];

                        dsmcNLocal += speciesDsmcN;
                        dsmcMLocal += sharedSampleCache_.dsmcM[typeId][cell];
                        dsmcLinearKELocal += speciesDsmcLinearKE;
                        dsmcMomentumLocal +=
                            sharedSampleCache_.dsmcMomentum[typeId][cell];
                        dsmcErotLocal += sharedSampleCache_.dsmcErot[typeId][cell];
                        dsmcZetaRotLocal +=
                            sharedSampleCache_.dsmcZetaRot[typeId][cell];
                        nLocal += speciesNReal;
                        mLocal += sharedSampleCache_.mReal[typeId][cell];
                        momentumLocal +=
                            sharedSampleCache_.momentumReal[typeId][cell];
                        linearKELocal +=
                            sharedSampleCache_.linearKEReal[typeId][cell];

                        dsmcNSpeciesCum_[i][cell] += speciesDsmcN;
                        dsmcMccSpeciesCum_[i][cell] += speciesDsmcLinearKE;
                        nSpeciesCum_[i][cell] += speciesNReal;

                        if (needVibrational)
                        {
                            forAll(dsmcSpeciesEvibModCum_[i], mod)
                            {
                                dsmcSpeciesEvibModCum_[i][mod][cell] +=
                                    sharedSampleCache_.dsmcSpeciesEvibMod[typeId][mod][cell];
                            }
                        }

                        if (needHeatFluxShearStress)
                        {
                            dsmcMuuLocal += sharedSampleCache_.dsmcMuu[typeId][cell];
                            dsmcMuvLocal += sharedSampleCache_.dsmcMuv[typeId][cell];
                            dsmcMuwLocal += sharedSampleCache_.dsmcMuw[typeId][cell];
                            dsmcMvvLocal += sharedSampleCache_.dsmcMvv[typeId][cell];
                            dsmcMvwLocal += sharedSampleCache_.dsmcMvw[typeId][cell];
                            dsmcMwwLocal += sharedSampleCache_.dsmcMww[typeId][cell];
                            dsmcMccLocal += sharedSampleCache_.dsmcMcc[typeId][cell];
                            dsmcMccuLocal += sharedSampleCache_.dsmcMccu[typeId][cell];
                            dsmcMccvLocal += sharedSampleCache_.dsmcMccv[typeId][cell];
                            dsmcMccwLocal += sharedSampleCache_.dsmcMccw[typeId][cell];
                            dsmcEuLocal += sharedSampleCache_.dsmcEu[typeId][cell];
                            dsmcEvLocal += sharedSampleCache_.dsmcEv[typeId][cell];
                            dsmcEwLocal += sharedSampleCache_.dsmcEw[typeId][cell];
                            dsmcECumLocal += sharedSampleCache_.dsmcECum[typeId][cell];
                        }

                        if (needClassification)
                        {
                            dsmcNClassILocal +=
                                sharedSampleCache_.dsmcNClassI[typeId][cell];
                            dsmcNClassIILocal +=
                                sharedSampleCache_.dsmcNClassII[typeId][cell];
                            dsmcNClassIIILocal +=
                                sharedSampleCache_.dsmcNClassIII[typeId][cell];
                        }
                    }
                }

                dsmcNCum_[cell] += dsmcNLocal;
                dsmcN_[cell] += dsmcNLocal;
                dsmcMCum_[cell] += dsmcMLocal;
                dsmcLinearKECum_[cell] += dsmcLinearKELocal;
                dsmcMomentumCum_[cell] += dsmcMomentumLocal;
                dsmcErotCum_[cell] += dsmcErotLocal;
                dsmcZetaRotCum_[cell] += dsmcZetaRotLocal;
                dsmcNElecLvlCum_[cell] += dsmcNElecLvlLocal;
                nCum_[cell] += nLocal;
                mCum_[cell] += mLocal;
                momentumCum_[cell] += momentumLocal;
                linearKECum_[cell] += linearKELocal;

                if (needHeatFluxShearStress)
                {
                    dsmcMuuCum_[cell] += dsmcMuuLocal;
                    dsmcMuvCum_[cell] += dsmcMuvLocal;
                    dsmcMuwCum_[cell] += dsmcMuwLocal;
                    dsmcMvvCum_[cell] += dsmcMvvLocal;
                    dsmcMvwCum_[cell] += dsmcMvwLocal;
                    dsmcMwwCum_[cell] += dsmcMwwLocal;
                    dsmcMccCum_[cell] += dsmcMccLocal;
                    dsmcMccuCum_[cell] += dsmcMccuLocal;
                    dsmcMccvCum_[cell] += dsmcMccvLocal;
                    dsmcMccwCum_[cell] += dsmcMccwLocal;
                    dsmcEuCum_[cell] += dsmcEuLocal;
                    dsmcEvCum_[cell] += dsmcEvLocal;
                    dsmcEwCum_[cell] += dsmcEwLocal;
                    dsmcECum_[cell] += dsmcECumLocal;
                }

                if (needClassification)
                {
                    dsmcNClassICum_[cell] += dsmcNClassILocal;
                    dsmcNClassIICum_[cell] += dsmcNClassIILocal;
                    dsmcNClassIIICum_[cell] += dsmcNClassIIILocal;
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
            const UList<label>* reduceCellsPtr =
                cloud_.replicatedMeshActive()
              ? static_cast<const UList<label>*>(&cloud_.replicatedMesh().myCells())
              : nullptr;
            const label reduceLoopSize =
                reduceCellsPtr ? reduceCellsPtr->size() : dsmcNCum_.size();
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) if (useOpenMPSampling)
            #endif
            for (label reduceI = 0; reduceI < reduceLoopSize; ++reduceI)
            {
                const label celli =
                    reduceCellsPtr ? (*reduceCellsPtr)[reduceI] : reduceI;
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
                    // not zero so that weighted decomposition still works
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
            const boundaryMeasurements& boundaryFlux =
                cloud_.boundaryFluxMeasurements();
            const bool useOwnedBoundaryFaces = ownedBoundaryFaces_.size() > 0;
            forAll(speciesIds_, i)
            {
                const label spId = speciesIds_[i];

                forAll(sampledBoundaryPatches_, patchi)
                {
                    const label j = sampledBoundaryPatches_[patchi];
                    const label nFaces =
                        useOwnedBoundaryFaces
                      ? ownedBoundaryFaces_[j].size()
                      : mesh_.boundaryMesh()[j].size();
                    for (label facei = 0; facei < nFaces; ++facei)
                    {
                        const label k =
                            useOwnedBoundaryFaces
                          ? ownedBoundaryFaces_[j][facei]
                          : facei;
                        rhoNBF_[j][k] +=
                            boundaryFlux.speciesRhoNBF(spId, j, k);
                        rhoMBF_[j][k] +=
                            boundaryFlux.speciesRhoMBF(spId, j, k);
                        linearKEBF_[j][k] +=
                            boundaryFlux.speciesLinearKEBF(spId, j, k);
                        momentumBF_[j][k] +=
                            boundaryFlux.speciesMomentumBF(spId, j, k);
                        ErotBF_[j][k] +=
                            boundaryFlux.speciesErotBF(spId, j, k);
                        zetaRotBF_[j][k] +=
                            boundaryFlux.speciesZetaRotBF(spId, j, k);
                        rhoNIntBF_[j][k] +=
                            boundaryFlux.speciesRhoNIntBF(spId, j, k);
                        rhoNElecBF_[j][k] +=
                            boundaryFlux.speciesRhoNElecBF(spId, j, k);
                        qBF_[j][k] +=
                            boundaryFlux.speciesqBF(spId, j, k);
                        fDBF_[j][k] +=
                            boundaryFlux.speciesfDBF(spId, j, k);
                        
                        speciesRhoNBF_[i][j][k] +=
                            boundaryFlux.speciesRhoNBF(spId, j, k);
                        speciesEvibBF_[i][j][k] +=
                            boundaryFlux.speciesEvibBF(spId, j, k);
                        speciesEelecBF_[i][j][k] +=
                            boundaryFlux.speciesEelecBF(spId, j, k);
                        speciesMccBF_[i][j][k] +=
                            boundaryFlux.speciesMccBF(spId, j, k);
                    }
                }

                forAll(speciesEvibModBF_[i], mod)
                {
                    forAll(sampledBoundaryPatches_, patchi)
                    {
                        const label j = sampledBoundaryPatches_[patchi];
                        const label nFaces =
                            useOwnedBoundaryFaces
                          ? ownedBoundaryFaces_[j].size()
                          : mesh_.boundaryMesh()[j].size();
                        for (label facei = 0; facei < nFaces; ++facei)
                        {
                            const label k =
                                useOwnedBoundaryFaces
                              ? ownedBoundaryFaces_[j][facei]
                              : facei;
                            speciesEvibModBF_[i][mod][j][k] +=
                                boundaryFlux.speciesEvibModBF(spId, mod, j, k);
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

    if (time_.time().outputTime())
    {
        const bool computeOutputFields =
            !cloud_.replicatedMeshActive() || cloud_.isOutputRank();

        if (cloud_.replicatedMeshActive())
        {
            sumReduceField(dsmcN_.primitiveFieldRef());
            sumReduceField(dsmcNCum_);
            sumReduceField(nCum_);
            sumReduceField(dsmcNElecLvlCum_);
            sumReduceField(dsmcMCum_);
            sumReduceField(mCum_);
            sumReduceField(dsmcLinearKECum_);
            sumReduceField(linearKECum_);
            sumReduceField(dsmcErotCum_);
            sumReduceField(dsmcZetaRotCum_);
            sumReduceField(dsmcMuuCum_);
            sumReduceField(dsmcMuvCum_);
            sumReduceField(dsmcMuwCum_);
            sumReduceField(dsmcMvvCum_);
            sumReduceField(dsmcMvwCum_);
            sumReduceField(dsmcMwwCum_);
            sumReduceField(dsmcMccCum_);
            sumReduceField(dsmcMccuCum_);
            sumReduceField(dsmcMccvCum_);
            sumReduceField(dsmcMccwCum_);
            sumReduceField(dsmcEuCum_);
            sumReduceField(dsmcEvCum_);
            sumReduceField(dsmcEwCum_);
            sumReduceField(dsmcECum_);
            sumReduceField(zetaVib_);
            sumReduceField(dsmcNClassICum_);
            sumReduceField(dsmcNClassIICum_);
            sumReduceField(dsmcNClassIIICum_);
            sumReduceField(collisionSeparation_);
            sumReduceField(dsmcNCollsCum_);
            sumReduceField(dsmcMomentumCum_);
            sumReduceField(momentumCum_);

            sumReduceFieldList(dsmcSpeciesEelecCum_);
            sumReduceFieldList(dsmcNSpeciesCum_);
            sumReduceFieldList(nSpeciesCum_);
            sumReduceFieldList(dsmcMccSpeciesCum_);
            sumReduceFieldList(dsmcNGrndElecLvlSpeciesCum_);
            sumReduceFieldList(dsmcN1stElecLvlSpeciesCum_);
            sumReduceFieldListList(dsmcSpeciesEvibModCum_);

            sumReduceFieldList(rhoNBF_);
            sumReduceFieldList(rhoMBF_);
            sumReduceFieldList(linearKEBF_);
            sumReduceFieldList(ErotBF_);
            sumReduceFieldList(zetaRotBF_);
            sumReduceFieldList(qBF_);
            sumReduceFieldList(zetaVibBF_);
            sumReduceFieldList(rhoNIntBF_);
            sumReduceFieldList(rhoNElecBF_);
            sumReduceFieldList(momentumBF_);
            sumReduceFieldList(fDBF_);

            sumReduceFieldListList(speciesEvibBF_);
            sumReduceFieldListList(speciesEelecBF_);
            sumReduceFieldListList(speciesRhoNBF_);
            sumReduceFieldListList(speciesMccBF_);
            sumReduceFieldListListList(speciesEvibModBF_);
        }

        if (computeOutputFields)
        {
            const scalar nAvTimeSteps = nTimeSteps_;

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
                const label nSpecies = speciesIds_.size();

            forAll(dsmcNCum_, celli)
            {
                //- Fields initialisation 
                scalar moleculesRhoN = 0.0;
                Tvib_[celli] = 0.0;
                scalarList speciesTvib(nSpecies, 0.0);
                List<scalarList> speciesTvibMod(nSpecies);
                scalarList speciesZetaVib(nSpecies, 0.0);
                List<scalarList> speciesZetaVibMod(nSpecies);
                
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
                    const label nVibMod =
                        cloud_.constProps(spId).nVibrationalModes();
                        
                    speciesZetaVibMod[i].setSize(nVibMod, 0.0);
                    speciesTvibMod[i].setSize(nVibMod, 0.0);
                    scalar zetaByTvibMod = 0.0;

                    forAll(dsmcSpeciesEvibModCum_[i], mod)
                    {
                        if
                        (
                            dsmcSpeciesEvibModCum_[i][mod][celli] > VSMALL
                         && dsmcNSpeciesCum_[i][celli] >= nMinParcelsTvib_
                         && speciesZetaVibMod.size() > SMALL
                        )
                        {
                            const scalar thetaV =
                                cloud_.constProps(spId).thetaV_m(mod);

                            const scalar iMean =
                                dsmcSpeciesEvibModCum_[i][mod][celli]
                               /(kB*thetaV*dsmcNSpeciesCum_[i][celli]);
                               
                            if (iMean > iMeanMinTvib_)
                            {
                                const scalar logFactor = log(1.0 + 1.0/iMean);

                                speciesTvibMod[i][mod] = thetaV/logFactor;

                                speciesZetaVibMod[i][mod] = 2.0*iMean*logFactor;

                                speciesZetaVib[i] += speciesZetaVibMod[i][mod];
                                    
                                zetaByTvibMod = speciesZetaVibMod[i][mod]
                                    *speciesTvibMod[i][mod];
                            }
                        }
                    }

                    if (speciesZetaVib[i] > SMALL)
                    {
                        moleculesRhoN += nSpeciesCum_[i][celli];
                        
                        speciesTvib[i] = zetaByTvibMod/speciesZetaVib[i];
                        
                        Tvib_[celli] += nSpeciesCum_[i][celli]*speciesTvib[i];
                            
                        zetaVib_[celli] += nSpeciesCum_[i][celli]
                            *speciesZetaVib[i];    
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
                }
                
                if (dsmcNCum_[celli] > SMALL and Ttra_[celli] > SMALL)
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

                    gamma = molarCp_trarot/molarCv_trarot;

                    const scalar speedOfSound = sqrt
                        (
                            gamma*kB/molecularMass*Ttra_[celli]
                        );

                    Ma_[celli] = mag(UMean_[celli])/speedOfSound;
                }
                else
                {
                    Ma_[celli] = 0.0;
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
                        DxToMfp_[celli] = 1.0/mfpToDx_[celli];
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
                            const label nVibMod =
                                cloud_.constProps(spId).nVibrationalModes();
                            
                            speciesZetaVibBF_[i][j][k] = 0.0;
                            speciesTvibBF_[i][j][k] = 0.0;
                            
                            scalar zetaByTvibMod = 0.0;
                            scalarList speciesZetaVibMod(nVibMod, 0.0);
                            scalarList speciesTvibMod(nVibMod, 0.0);

                            if (speciesRhoNBF_[i][j][k] > SMALL)
                            {
                                forAll(speciesZetaVibMod, mod)
                                {
                                    const scalar thetaV =
                                        cloud_.constProps(spId).thetaV()[mod];

                                    const scalar iMean =
                                        speciesEvibModBF_[i][mod][j][k]
                                       /(kB*thetaV*speciesRhoNBF_[i][j][k]);

                                    if (iMean > iMeanMinTvib_)
                                    {
                                        const scalar logFactor =
                                            log(1.0 + 1.0/iMean);
                                        
                                        speciesTvibMod[mod] = thetaV/logFactor;

                                        speciesZetaVibMod[mod] =
                                            2.0*iMean*logFactor;

                                        speciesZetaVibBF_[i][j][k] +=
                                            speciesZetaVibMod[mod];
                                            
                                        zetaByTvibMod += speciesZetaVibMod[mod]
                                            *speciesTvibMod[mod];
                                    }
                                }
                            }

                            if (speciesZetaVibBF_[i][j][k] > SMALL)
                            {
                                moleculesRhoN += speciesRhoNBF_[i][j][k];
                                
                                speciesTvibBF_[i][j][k] = zetaByTvibMod
                                    /speciesZetaVibBF_[i][j][k];
                                    
                                Tvib_.boundaryFieldRef()[j][k] +=
                                    speciesRhoNBF_[i][j][k]
                                   *speciesTvibBF_[i][j][k];

                                zetaVibBF_[j][k] +=
                                    speciesRhoNBF_[i][j][k]
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
                            );

                        if (rhoNBF_[j][k] > SMALL)
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

                            //- Mach number calculation
                            const scalar gasConstant = 
                                (
                                    rhoNBF_[j][k] > SMALL
                                 && Ttra_.boundaryFieldRef()[j][k] > SMALL
                                  ? kB/molecularMassBF
                                  : 0.0
                                );

                            const scalar gamma = molarCpBF_trarot/molarCvBF_trarot;

                            const scalar speedOfSound =
                                sqrt
                                (
                                    gamma*gasConstant*Ttra_.boundaryField()[j][k]
                                );

                            Ma_.boundaryFieldRef()[j][k] =
                                mag(UMean_.boundaryField()[j][k])/speedOfSound;
                        }
                        else
                        {
                            Ma_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Force density
                        fD_.boundaryFieldRef()[j][k] = fDBF_[j][k]/nAvTimeSteps;

                        //- Surface pressure
                        p_.boundaryFieldRef()[j][k] =
                            fD_.boundaryField()[j][k] & n_[j][k];
                            
                        //- Wall shear stress
                        tau_.boundaryFieldRef()[j][k] =
                            sqrt
                            (
                                sqr(fD_.boundaryField()[j][k] & t1_[j][k])
                              + sqr(fD_.boundaryField()[j][k] & t2_[j][k])
                            );
                            
                        //- Heat flux
                        q_.boundaryFieldRef()[j][k] = qBF_[j][k]/nAvTimeSteps;
                        
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
            
            //- Write solution fields
            p_.write();
            Ttra_.write();
            UMean_.write();
            Ma_.write();
            q_.write();
            fD_.write();
            tau_.write();
            
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
        }
        }
        
        //- Reset fields after printing the instantaneous solution ... or
        //  continue sampling
        if (time_.resetFieldsAtOutput())
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

        if (averagingAcrossManyRuns_ && !time_.resetFieldsAtOutput())
        {
            writeOut();
        }
    }

    if
    (
        profileDetailEnabled_
     && cloud_.isOutputRank()
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
            << "    field combine         = " << profileFieldCombineWallTime_ << " s" << nl
            << "    cell reduction        = " << profileCellReduceWallTime_ << " s" << nl
            << "    boundary accumulation = " << profileBoundaryAccumWallTime_ << " s" << nl
            << "    output compute        = " << profileOutputComputeWallTime_ << " s" << nl
            << "    field writes          = " << profileFieldWriteWallTime_ << " s" << nl
            << "    output reset          = " << profileOutputResetWallTime_ << " s" << nl
            << "    output-time block     = " << profileOutputTimeWallTime_ << " s" << nl
            << "      detail sample cells = " << profileSharedCacheDetailSampleCells_ << nl
            << "      detail sample parcels = " << profileSharedCacheDetailSampleParcels_ << nl
            << "      base accum (sampled)= " << profileSharedCacheBaseAccumWallTime_ << " s" << nl
            << "      vib accum (sampled) = " << profileSharedCacheVibAccumWallTime_ << " s" << nl
            << "      electronic (sampled)= " << profileSharedCacheElectronicAccumWallTime_ << " s" << nl
            << "      class accum (sampled)= " << profileSharedCacheClassAccumWallTime_ << " s" << nl
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
