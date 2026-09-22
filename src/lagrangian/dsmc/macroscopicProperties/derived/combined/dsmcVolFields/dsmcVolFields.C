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
#include "dsmcMasterInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
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

    bool readMixtureFreestream
    (
        const Foam::dictionary& properties,
        const Foam::dictionary& moleculeProperties,
        Foam::scalar& rhoInf,
        Foam::vector& UInf,
        Foam::scalar& pInf
    )
    {
        if
        (
           !properties.found("velocity")
         || !properties.found("numberDensities")
         || !properties.found("translationalTemperature")
        )
        {
            return false;
        }

        UInf = Foam::vector(properties.lookup("velocity"));
        rhoInf = 0.0;
        pInf = 0.0;
        const Foam::scalar translationalTemperature = Foam::readScalar
        (
            properties.lookup("translationalTemperature")
        );

        if (translationalTemperature <= Foam::VSMALL)
        {
            return false;
        }

        const Foam::dictionary& numberDensities =
            properties.subDict("numberDensities");

        forAllConstIter(Foam::dictionary, numberDensities, iter)
        {
            const Foam::word& speciesName = iter().keyword();

            if (!moleculeProperties.found(speciesName))
            {
                return false;
            }

            const Foam::scalar numberDensity = Foam::readScalar
            (
                numberDensities.lookup(speciesName)
            );
            const Foam::scalar mass = Foam::readScalar
            (
                moleculeProperties.subDict(speciesName).lookup("mass")
            );

            if (numberDensity < 0.0 || mass <= Foam::VSMALL)
            {
                return false;
            }

            rhoInf += numberDensity*mass;
            pInf +=
                numberDensity
               *physicoChemical::k.value()
               *translationalTemperature;
        }

        return
            rhoInf > Foam::VSMALL
         && Foam::mag(UInf) > Foam::VSMALL
         && pInf > Foam::VSMALL;
    }


    bool readAutomaticFreestream
    (
        const Foam::fvMesh& mesh,
        const Foam::dsmcCloud& cloud,
        Foam::scalar& rhoInf,
        Foam::vector& UInf,
        Foam::scalar& pInf,
        Foam::word& source
    )
    {
        const Foam::dictionary& moleculeProperties =
            cloud.particleProperties().subDict("moleculeProperties");

        Foam::IOdictionary boundariesDict
        (
            Foam::IOobject
            (
                "boundariesDict",
                mesh.time().system(),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            )
        );

        if (boundariesDict.found("dsmcGeneralBoundaries"))
        {
            const Foam::PtrList<Foam::entry> generalBoundaries
            (
                boundariesDict.lookup("dsmcGeneralBoundaries")
            );

            forAll(generalBoundaries, boundaryI)
            {
                const Foam::dictionary& boundary =
                    generalBoundaries[boundaryI].dict();
                const Foam::word boundaryModel =
                    boundary.lookupOrDefault<Foam::word>
                    (
                        "boundaryModel",
                        Foam::word::null
                    );

                if
                (
                    boundaryModel != "dsmcFreeStreamInflowPatch"
                 && boundaryModel != "dsmcChapmanEnskogFreeStreamInflowPatch"
                )
                {
                    continue;
                }

                const Foam::word propertiesName =
                    boundaryModel + "Properties";

                if
                (
                    boundary.found(propertiesName)
                 && readMixtureFreestream
                    (
                        boundary.subDict(propertiesName),
                        moleculeProperties,
                        rhoInf,
                        UInf,
                        pInf
                    )
                )
                {
                    Foam::word patchName("unknownPatch");
                    if (boundary.found("generalBoundaryProperties"))
                    {
                        patchName = boundary.subDict
                        (
                            "generalBoundaryProperties"
                        ).lookupOrDefault<Foam::word>("patchName", patchName);
                    }

                    source = boundaryModel + "_" + patchName;
                    return true;
                }
            }
        }

        Foam::IOdictionary dsmcInitialiseDict
        (
            Foam::IOobject
            (
                "dsmcInitialiseDict",
                mesh.time().system(),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            )
        );

        if (!dsmcInitialiseDict.found("configurations"))
        {
            return false;
        }

        const Foam::PtrList<Foam::entry> configurations
        (
            dsmcInitialiseDict.lookup("configurations")
        );

        if (!configurations.size())
        {
            return false;
        }

        if
        (
            readMixtureFreestream
            (
                configurations[0].dict(),
                moleculeProperties,
                rhoInf,
                UInf,
                pInf
            )
        )
        {
            source = "dsmcInitialiseDict_configuration0";
            return true;
        }

        return false;
    }


    void calculateAutomaticReferenceGeometry
    (
        const Foam::fvMesh& mesh,
        const Foam::labelList& wallPatchIds,
        const Foam::vector& dragDirection,
        const bool replicatedMesh,
        Foam::scalar& referenceArea,
        Foam::scalar& referenceLength,
        Foam::vector& referencePoint,
        Foam::word& areaDefinition,
        Foam::word& lengthDefinition
    )
    {
        int mpiInitialised = 0;
        MPI_Initialized(&mpiInitialised);

        int rank = 0;
        if (mpiInitialised)
        {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        }

        const bool contributes = !replicatedMesh || rank == 0;
        Foam::scalar projectedAreaSum = 0.0;
        Foam::scalar wettedAreaSum = 0.0;
        Foam::vector surfaceAreaSum(Foam::vector::zero);
        Foam::vector weightedCentreSum(Foam::vector::zero);
        Foam::scalar minProjection = Foam::GREAT;
        Foam::scalar maxProjection = -Foam::GREAT;

        if (contributes)
        {
            forAll(wallPatchIds, wallPatchI)
            {
                const Foam::polyPatch& patch =
                    mesh.boundaryMesh()[wallPatchIds[wallPatchI]];

                forAll(patch, faceI)
                {
                    const Foam::label meshFaceI = patch.start() + faceI;
                    const Foam::vector& faceArea =
                        mesh.faceAreas()[meshFaceI];
                    const Foam::scalar faceAreaMagnitude = Foam::mag(faceArea);

                    projectedAreaSum +=
                        Foam::mag(faceArea & dragDirection);
                    wettedAreaSum += faceAreaMagnitude;
                    surfaceAreaSum += faceArea;
                    weightedCentreSum +=
                        faceAreaMagnitude*mesh.faceCentres()[meshFaceI];

                    const Foam::face& face = mesh.faces()[meshFaceI];
                    forAll(face, pointI)
                    {
                        const Foam::scalar projection =
                            mesh.points()[face[pointI]] & dragDirection;
                        minProjection = Foam::min(minProjection, projection);
                        maxProjection = Foam::max(maxProjection, projection);
                    }
                }
            }
        }

        if (mpiInitialised)
        {
            Foam::scalar scalarSums[2] =
            {
                projectedAreaSum,
                wettedAreaSum
            };
            Foam::scalar vectorSums[6] =
            {
                surfaceAreaSum.x(),
                surfaceAreaSum.y(),
                surfaceAreaSum.z(),
                weightedCentreSum.x(),
                weightedCentreSum.y(),
                weightedCentreSum.z()
            };
            Foam::scalar projectionBounds[2] =
            {
                minProjection,
                maxProjection
            };

            MPI_Allreduce
            (
                MPI_IN_PLACE,
                scalarSums,
                2,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD
            );
            MPI_Allreduce
            (
                MPI_IN_PLACE,
                vectorSums,
                6,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD
            );
            MPI_Allreduce
            (
                MPI_IN_PLACE,
                projectionBounds,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                MPI_COMM_WORLD
            );
            MPI_Allreduce
            (
                MPI_IN_PLACE,
                projectionBounds + 1,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD
            );

            projectedAreaSum = scalarSums[0];
            wettedAreaSum = scalarSums[1];
            surfaceAreaSum = Foam::vector
            (
                vectorSums[0],
                vectorSums[1],
                vectorSums[2]
            );
            weightedCentreSum = Foam::vector
            (
                vectorSums[3],
                vectorSums[4],
                vectorSums[5]
            );
            minProjection = projectionBounds[0];
            maxProjection = projectionBounds[1];
        }

        if (wettedAreaSum <= Foam::VSMALL)
        {
            FatalErrorInFunction
                << "No wall-face area is available for automatic force-moment "
                << "reference geometry." << exit(FatalError);
        }

        const bool closedSurface =
            Foam::mag(surfaceAreaSum) <= 1.0e-8*wettedAreaSum;

        referenceArea =
            (closedSurface ? 0.5 : 1.0)*projectedAreaSum;
        areaDefinition =
            closedSurface
          ? "closedWallProjectedArea"
          : "wallProjectedArea";

        if (referenceArea <= Foam::VSMALL)
        {
            FatalErrorInFunction
                << "Automatic reference area is zero in the free-stream "
                << "direction " << dragDirection << "."
                << exit(FatalError);
        }

        const Foam::scalar streamwiseLength = maxProjection - minProjection;
        if (streamwiseLength > Foam::VSMALL)
        {
            referenceLength = streamwiseLength;
            lengthDefinition = "wallPointStreamwiseSpan";
        }
        else
        {
            referenceLength = Foam::sqrt(referenceArea);
            lengthDefinition = "sqrtProjectedAreaFallback";
        }

        referencePoint = weightedCentreSum/wettedAreaSum;
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

    template<class GeoField>
    void setProcessorWriteOpt(GeoField& field, const bool enabled)
    {
        field.writeOpt() =
            enabled ? Foam::IOobject::AUTO_WRITE : Foam::IOobject::NO_WRITE;
    }

    template<class Type>
    void snapshotField
    (
        const Foam::Field<Type>& src,
        Foam::Field<Type>& dst
    )
    {
        dst.setSize(src.size());
        forAll(src, i)
        {
            dst[i] = src[i];
        }
    }

    template<class Type>
    void restoreField
    (
        Foam::Field<Type>& dst,
        const Foam::Field<Type>& src
    )
    {
        dst.setSize(src.size());
        forAll(src, i)
        {
            dst[i] = src[i];
        }
    }

    template<class Type>
    void snapshotFieldList
    (
        const Foam::List<Foam::Field<Type>>& src,
        Foam::List<Foam::Field<Type>>& dst
    )
    {
        dst.setSize(src.size());
        forAll(src, i)
        {
            snapshotField(src[i], dst[i]);
        }
    }

    template<class Type>
    void restoreFieldList
    (
        Foam::List<Foam::Field<Type>>& dst,
        const Foam::List<Foam::Field<Type>>& src
    )
    {
        dst.setSize(src.size());
        forAll(src, i)
        {
            restoreField(dst[i], src[i]);
        }
    }

    template<class Type>
    void snapshotFieldListList
    (
        const Foam::List<Foam::List<Foam::Field<Type>>>& src,
        Foam::List<Foam::List<Foam::Field<Type>>>& dst
    )
    {
        dst.setSize(src.size());
        forAll(src, i)
        {
            snapshotFieldList(src[i], dst[i]);
        }
    }

    template<class Type>
    void restoreFieldListList
    (
        Foam::List<Foam::List<Foam::Field<Type>>>& dst,
        const Foam::List<Foam::List<Foam::Field<Type>>>& src
    )
    {
        dst.setSize(src.size());
        forAll(src, i)
        {
            restoreFieldList(dst[i], src[i]);
        }
    }

    template<class Type>
    void snapshotFieldListListList
    (
        const Foam::List<Foam::List<Foam::List<Foam::Field<Type>>>>& src,
        Foam::List<Foam::List<Foam::List<Foam::Field<Type>>>>& dst
    )
    {
        dst.setSize(src.size());
        forAll(src, i)
        {
            snapshotFieldListList(src[i], dst[i]);
        }
    }

    template<class Type>
    void restoreFieldListListList
    (
        Foam::List<Foam::List<Foam::List<Foam::Field<Type>>>>& dst,
        const Foam::List<Foam::List<Foam::List<Foam::Field<Type>>>>& src
    )
    {
        dst.setSize(src.size());
        forAll(src, i)
        {
            restoreFieldListList(dst[i], src[i]);
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
        scalarField totalDsmcN;
        scalarField totalDsmcM;
        scalarField totalDsmcLinearKE;
        vectorField totalDsmcMomentum;
        scalarField totalDsmcErot;
        scalarField totalDsmcZetaRot;
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

        // Capacity-retained initialisation: the arrays persist across steps
        // and are cleared by resetFields() (sparse over touchedCells) before
        // each accumulate pass.  Fresh allocations are zeroed by setSize;
        // the previous full-operator= zeroing pass per step is removed
        // (worklog dlb_further_opt.md §8.4, cache allocate 25.7 s/200 steps).
        void initScalarFields(List<scalarField>& fields)
        {
            if (fields.size() != nTypes)
            {
                fields.setSize(nTypes);
            }
            for (label typei = 0; typei < nTypes; ++typei)
            {
                if (fields[typei].size() != nCells)
                {
                    fields[typei].setSize(nCells, 0.0);
                }
            }
        }

        void initVectorFields(List<vectorField>& fields)
        {
            if (fields.size() != nTypes)
            {
                fields.setSize(nTypes);
            }
            for (label typei = 0; typei < nTypes; ++typei)
            {
                if (fields[typei].size() != nCells)
                {
                    fields[typei].setSize(nCells, vector::zero);
                }
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
            // Totals persist as well; resetFields() zeroes them sparsely via
            // touchedCells before each accumulate pass.
            if (totalDsmcN.size() != nCells) totalDsmcN.setSize(nCells, 0.0);
            if (totalDsmcM.size() != nCells) totalDsmcM.setSize(nCells, 0.0);
            if (totalDsmcLinearKE.size() != nCells) totalDsmcLinearKE.setSize(nCells, 0.0);
            if (totalDsmcMomentum.size() != nCells) totalDsmcMomentum.setSize(nCells, vector::zero);
            if (totalDsmcErot.size() != nCells) totalDsmcErot.setSize(nCells, 0.0);
            if (totalDsmcZetaRot.size() != nCells) totalDsmcZetaRot.setSize(nCells, 0.0);
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
                        if (dsmcSpeciesEvibMod[typei][mod].size() != nCells)
                        {
                            dsmcSpeciesEvibMod[typei][mod].setSize(nCells, 0.0);
                        }
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

                    if (typei == 0)
                    {
                        totalDsmcN = 0.0;
                        totalDsmcM = 0.0;
                        totalDsmcLinearKE = 0.0;
                        totalDsmcMomentum = vector::zero;
                        totalDsmcErot = 0.0;
                        totalDsmcZetaRot = 0.0;
                    }
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

                        if (typei == 0)
                        {
                            totalDsmcN[celli] = 0.0;
                            totalDsmcM[celli] = 0.0;
                            totalDsmcLinearKE[celli] = 0.0;
                            totalDsmcMomentum[celli] = vector::zero;
                            totalDsmcErot[celli] = 0.0;
                            totalDsmcZetaRot[celli] = 0.0;
                        }
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
                useFlatOccupancy
              ? static_cast<const UList<label>*>(&cloud.occupancyActiveCells())
              : nullptr;

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

                // The arrays persist across steps: clear only the cells this
                // build will write (the current occupancy active set — every
                // parcel accumulates into its cell) — the same set the field
                // derivation reads.  Sparse reset via touchedCells.
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

                const auto resetStart =
                    doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
                resetFields();
                if (doProfile)
                {
                    buildProfile->resetWallTime +=
                        wallSeconds(resetStart, wallClockNow());
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

	            scalarField typeMass(nTypes, 0.0);
	            scalarField typeZetaRot(nTypes, 0.0);
	            labelList typeNVibModes(nTypes, 0);
	            List<scalarField> typeVibEnergyQuantum;

            if (needVibrational)
            {
                typeVibEnergyQuantum.setSize(nTypes);
            }

            const scalar kB = physicoChemical::k.value();
            for (label typei = 0; typei < nTypes; ++typei)
            {
                const dsmcParcel::constantProperties& cP =
                    cloud.constProps(typei);

                typeMass[typei] = cP.mass();
                typeZetaRot[typei] = cP.rotationalDegreesOfFreedom();

	                if (needVibrational)
	                {
	                    const scalarList& thetaV = cP.thetaV();
	                    typeNVibModes[typei] = thetaV.size();
	                    typeVibEnergyQuantum[typei].setSize(thetaV.size(), 0.0);
	                    forAll(thetaV, mod)
                    {
                        typeVibEnergyQuantum[typei][mod] = kB*thetaV[mod];
                    }
                }
            }

            const auto parcelAccumStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();

            auto accumulateIntoCache =
                [&]
                (
                    const label celli,
                    const dsmcParcel& p,
                    scalarField& localDsmcN,
                    scalarField& localDsmcM,
                    scalarField& localDsmcLinearKE,
                    vectorField& localDsmcMomentum,
                    scalarField& localDsmcErot,
                    scalarField& localDsmcZetaRot,
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

                if (localDsmcN[typeId] == 0.0)
                {
                    activeTypes.append(typeId);
                }

                const scalar mp = typeMass[typeId];
                const vector& Up = p.U();
                const scalar Upx = Up.x();
                const scalar Upy = Up.y();
                const scalar Upz = Up.z();
                const vector mpUp = mp*Up;
                const scalar linearKE = mp*(Up & Up);
                const scalar Erotp = p.ERot();

                scalar Evibp = 0.0;
	                if (needVibrational && typeNVibModes[typeId])
	                {
	                    const labelList& vibLevels = p.vibLevel();
	                    const scalarField& vibEnergyQuantum =
                        typeVibEnergyQuantum[typeId];
                    scalarField& localEvibMods = localSpeciesEvibMod[typeId];

                    forAll(vibEnergyQuantum, mod)
                    {
                        const scalar EvibMod =
                            vibLevels[mod]*vibEnergyQuantum[mod];
                        localEvibMods[mod] += EvibMod;
                        if (needHeatFluxShearStress)
                        {
                            Evibp += EvibMod;
                        }
                    }
                }

                localDsmcN[typeId] += 1.0;
                localDsmcLinearKE[typeId] += linearKE;
                localDsmcMomentum[typeId] += mpUp;
                localDsmcErot[typeId] += Erotp;

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
                    const dsmcParcel::constantProperties& cP =
                        cloud.constProps(typeId);
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
                    const scalar localN = localDsmcN[typeId];
                    const scalar localM = localN*typeMass[typeId];
                    const scalar localLinearKE = localDsmcLinearKE[typeId];
                    const vector localMomentum = localDsmcMomentum[typeId];
                    const scalar localErot = localDsmcErot[typeId];
                    const scalar localZetaRot = localN*typeZetaRot[typeId];

                    dsmcN[typeId][celli] = localN;
                    dsmcM[typeId][celli] = localM;
                    dsmcLinearKE[typeId][celli] = localLinearKE;
                    dsmcMomentum[typeId][celli] = localMomentum;
                    dsmcErot[typeId][celli] = localErot;
                    dsmcZetaRot[typeId][celli] = localZetaRot;

                    totalDsmcN[celli] += localN;
                    totalDsmcM[celli] += localM;
                    totalDsmcLinearKE[celli] += localLinearKE;
                    totalDsmcMomentum[celli] += localMomentum;
                    totalDsmcErot[celli] += localErot;
                    totalDsmcZetaRot[celli] += localZetaRot;

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
                const bool baseVibFastPath =
                    needVibrational
                 && !needElectronic
                 && !needClassification
                 && !needHeatFluxShearStress;

                if (baseVibFastPath)
                {
                    #pragma omp parallel
                    {
                        scalarField localDsmcN(nTypes, 0.0);
                        scalarField localUSqr(nTypes, 0.0);
                        vectorField localU(nTypes, vector::zero);
                        scalarField localDsmcErot(nTypes, 0.0);
                        List<scalarField> localSpeciesEvibMod(nTypes);
	                        labelList activeTypes(nTypes, -1);
	                        label nActiveTypes = 0;
	                        label localDetailSampleCells = 0;
                        label localDetailSampleParcels = 0;

                        for (label typei = 0; typei < nTypes; ++typei)
                        {
                            localSpeciesEvibMod[typei].setSize
                            (
                                typeVibEnergyQuantum[typei].size(),
                                0.0
                            );
                        }

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

                            if
                            (
                                doDetailedProfile
                             && parcelCount > 0
                             && (celli % 32 == 0)
                            )
                            {
                                ++localDetailSampleCells;
                                localDetailSampleParcels += parcelCount;
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

                                if (!p.isFree())
                                {
                                    continue;
                                }

                                const label typeId = p.typeId();
                                if (typeId < 0 || typeId >= nTypes)
                                {
                                    continue;
                                }

	                                if (localDsmcN[typeId] == 0.0)
	                                {
	                                    activeTypes[nActiveTypes++] = typeId;
	                                }

                                const vector& Up = p.U();
                                const scalar Erotp = p.ERot();

	                                if (typeNVibModes[typeId])
	                                {
	                                    const labelList& vibLevels = p.vibLevel();
	                                    const scalarField& vibEnergyQuantum =
	                                        typeVibEnergyQuantum[typeId];
	                                    scalarField& localEvibMods =
	                                        localSpeciesEvibMod[typeId];

	                                    for (label mod = 0; mod < typeNVibModes[typeId]; ++mod)
	                                    {
	                                        const scalar EvibMod =
	                                            vibLevels[mod]*vibEnergyQuantum[mod];
	                                        localEvibMods[mod] += EvibMod;
	                                    }
	                                }

	                                localDsmcN[typeId] += 1.0;
	                                localUSqr[typeId] += (Up & Up);
	                                localU[typeId] += Up;
	                                localDsmcErot[typeId] += Erotp;
                            }

		                            for (label activeI = 0; activeI < nActiveTypes; ++activeI)
		                            {
		                                const label typeId = activeTypes[activeI];
		                                const scalar localN = localDsmcN[typeId];
		                                const scalar typeM =
		                                    typeMass[typeId];
		                                const scalar localM =
		                                    localN*typeM;
		                                const scalar localLinearKE =
		                                    typeM*localUSqr[typeId];
		                                const vector localMomentum =
		                                    typeM*localU[typeId];
		                                const scalar localErot =
		                                    localDsmcErot[typeId];
		                                const scalar localZetaRot =
		                                    localN*typeZetaRot[typeId];

	                                dsmcN[typeId][celli] = localN;
	                                dsmcM[typeId][celli] = localM;
	                                dsmcLinearKE[typeId][celli] =
	                                    localLinearKE;
	                                dsmcMomentum[typeId][celli] =
	                                    localMomentum;
	                                dsmcErot[typeId][celli] = localErot;
	                                dsmcZetaRot[typeId][celli] =
	                                    localZetaRot;

	                                totalDsmcN[celli] += localN;
	                                totalDsmcM[celli] += localM;
	                                totalDsmcLinearKE[celli] += localLinearKE;
	                                totalDsmcMomentum[celli] += localMomentum;
	                                totalDsmcErot[celli] += localErot;
	                                totalDsmcZetaRot[celli] += localZetaRot;

                                scalarField& localEvibMods =
                                    localSpeciesEvibMod[typeId];
                                forAll(localEvibMods, mod)
                                {
                                    dsmcSpeciesEvibMod[typeId][mod][celli] =
                                        localEvibMods[mod];
                                    localEvibMods[mod] = 0.0;
                                }

	                                localDsmcN[typeId] = 0.0;
	                                localUSqr[typeId] = 0.0;
	                                localU[typeId] = vector::zero;
	                                localDsmcErot[typeId] = 0.0;
                            }

	                            nActiveTypes = 0;
	                        }

	                        if (doDetailedProfile)
	                        {
	                            #pragma omp atomic
	                            buildProfile->detailSampleCells += localDetailSampleCells;
	                            #pragma omp atomic
	                            buildProfile->detailSampleParcels += localDetailSampleParcels;
	                        }
                    }
                }
                else
                {
                #pragma omp parallel
                {
                    scalarField localDsmcN(nTypes, 0.0);
                    scalarField localDsmcM(nTypes, 0.0);
                    scalarField localDsmcLinearKE(nTypes, 0.0);
                    vectorField localDsmcMomentum(nTypes, vector::zero);
                    scalarField localDsmcErot(nTypes, 0.0);
                    scalarField localDsmcZetaRot(nTypes, 0.0);
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
                    label localDetailSampleCells = 0;
                    label localDetailSampleParcels = 0;

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

                        if (doDetailedProfile && parcelCount > 0 && (celli % 32 == 0))
                        {
                            ++localDetailSampleCells;
                            localDetailSampleParcels += parcelCount;
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
                                localDsmcN,
                                localDsmcM,
                                localDsmcLinearKE,
                                localDsmcMomentum,
                                localDsmcErot,
                                localDsmcZetaRot,
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

                    if (doDetailedProfile)
                    {
                        #pragma omp atomic
                        buildProfile->detailSampleCells += localDetailSampleCells;
                        #pragma omp atomic
                        buildProfile->detailSampleParcels += localDetailSampleParcels;
                    }
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

                    const scalar mp = typeMass[typeId];
                    const vector& Up = p.U();
                    const scalar linearKE = mp*(Up & Up);
                    const scalar Erotp = p.ERot();
                    const scalar zetaRotp = typeZetaRot[typeId];

                    scalar Evibp = 0.0;
	                    if (needVibrational && typeNVibModes[typeId])
	                    {
	                        const labelList& vibLevels = p.vibLevel();
	                        const scalarField& vibEnergyQuantum =
                            typeVibEnergyQuantum[typeId];
                        forAll(vibEnergyQuantum, mod)
                        {
                            const scalar EvibMod =
                                vibLevels[mod]*vibEnergyQuantum[mod];
                            dsmcSpeciesEvibMod[typeId][mod][celli] += EvibMod;
                            if (needHeatFluxShearStress)
                            {
                                Evibp += EvibMod;
                            }
                        }
                    }

                    dsmcN[typeId][celli] += 1.0;
                    dsmcM[typeId][celli] += mp;
                    dsmcLinearKE[typeId][celli] += linearKE;
                    dsmcMomentum[typeId][celli] += mp*Up;
                    dsmcErot[typeId][celli] += Erotp;
                    dsmcZetaRot[typeId][celli] += zetaRotp;

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
                        const dsmcParcel::constantProperties& cP =
                            cloud.constProps(typeId);
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


bool dsmcVolFields::outputFieldEnabled
(
    const word& fieldKey,
    const word& objectName,
    const bool defaultEnabled
) const
{
    if (!defaultEnabled)
    {
        return false;
    }

    if (!restrictOutputFields_)
    {
        return true;
    }

    if
    (
        outputFieldNames_.found("all")
     || outputFieldNames_.found(fieldKey)
     || outputFieldNames_.found(objectName)
    )
    {
        return true;
    }

    if
    (
        outputFieldNames_.found("basic")
     && (
            fieldKey == "dsmcN"
         || fieldKey == "dsmcNMean"
         || fieldKey == "rhoN"
         || fieldKey == "rhoM"
         || fieldKey == "p"
         || fieldKey == "Ttra"
         || fieldKey == "U"
         || fieldKey == "Ma"
        )
    )
    {
        return true;
    }

    if
    (
        outputFieldNames_.found("thermal")
     && (
            fieldKey == "Ttra"
         || fieldKey == "Trot"
         || fieldKey == "Tvib"
         || fieldKey == "Telec"
         || fieldKey == "Tov"
        )
    )
    {
        return true;
    }

    if
    (
        outputFieldNames_.found("wall")
     && (
            fieldKey == "wallHeatFlux"
         || fieldKey == "wallShearStress"
         || fieldKey == "fD"
         || fieldKey == "wallPressureCoefficient"
         || fieldKey == "wallHeatFluxCoefficient"
        )
    )
    {
        return true;
    }

    if
    (
        outputFieldNames_.found("transport")
     && (
            fieldKey == "mfp"
         || fieldKey == "mfpToDx"
         || fieldKey == "mct"
         || fieldKey == "mctToDt"
         || fieldKey == "SOFP"
        )
    )
    {
        return true;
    }

    if
    (
        outputFieldNames_.found("errors")
     && (
            fieldKey == "rhoMError"
         || fieldKey == "UError"
         || fieldKey == "TError"
         || fieldKey == "pError"
        )
    )
    {
        return true;
    }

    if
    (
        outputFieldNames_.found("heatFluxShearStress")
     && (
            fieldKey == "heatFluxVector"
         || fieldKey == "pressureTensor"
         || fieldKey == "shearStressTensor"
        )
    )
    {
        return true;
    }

    if
    (
        outputFieldNames_.found("classifications")
     && (
            fieldKey == "classIDistribution"
         || fieldKey == "classIIDistribution"
         || fieldKey == "classIIIDistribution"
        )
    )
    {
        return true;
    }

    return false;
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
    iMeanMinTvib_(0.01),
    fieldName_(propsDict_.lookup("fieldName")),
    speciesIds_(),
    typeIdToSpeciesIndex_(),
    Cp_(),
    Ch_(),
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
    restrictOutputFields_(false),
    outputFieldNames_(),
    writeForceMoment_(false),
    forceMomentReferencePoint_(vector::zero),
    forceMomentWallPatchIds_(),
    forceMomentWallPatchIndices_(),
    writeForceMomentCoefficients_(false),
    writeWallPressureCoefficient_(false),
    forceMomentQInf_(0.0),
    forceMomentPInf_(0.0),
    forceMomentHeatFluxInf_(0.0),
    forceMomentReferenceArea_(0.0),
    forceMomentReferenceLength_(0.0),
    forceMomentDragDirection_(vector::zero),
    forceMomentLiftDirection_(vector::zero),
    forceMomentSideDirection_(vector::zero),
    forceMomentFreestreamSource_(word::null),
    forceMomentReferenceAreaDefinition_(word::null),
    forceMomentReferenceLengthDefinition_(word::null),
    forceMomentReferencePointAutomatic_(false),
    forceMomentAutomaticReferencePoint_(vector::zero),
    forceMomentReferenceGeometryDirection_(vector::zero),
    forceMomentReferenceGeometryCached_(false),
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
    profileOutputSnapshotWallTime_(0.0),
    profileOutputReduceWallTime_(0.0),
    profileOutputComputeWallTime_(0.0),
    profileFieldWriteWallTime_(0.0),
    profileOutputResetWallTime_(0.0),
    profileOutputRestoreWallTime_(0.0),
    profileOutputTimeWallTime_(0.0),
    profileCalls_(0),
    finalProfilePrinted_(false),
    ownedBoundaryFaces_(),
    ownedBoundaryFacesOwnerVersion_(-1),
    ownedStorageActive_(false),
    ownedStorageOwnerVersion_(-1),
    toOwnedCell_(),
    pendingResumeReadIn_(false)
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

// M4: copy the owned entries (positions 0..nOwned-1 in the current storage)
// into a full-mesh field at the global cell positions.  Used when the
// resume-sampling file must keep its full-mesh layout in owned mode.
template<class Type>
static void scatterOwnedToFull
(
    const Field<Type>& owned,
    Field<Type>& full,
    const UList<label>& myCells
)
{
    full.setSize(myCells.size(), pTraits<Type>::zero);
    full = pTraits<Type>::zero;
    forAll(myCells, i)
    {
        full[myCells[i]] = owned[i];
    }
}


// M4: keep only the owned entries of a field whose size reflects another
// layout (full mesh, or a pre-rebalance owned set).
template<class Type>
static void gatherOwnedEntries
(
    Field<Type>& f,
    const UList<label>& myCells
)
{
    if (f.size() == myCells.size())
    {
        return;
    }

    Field<Type> owned(myCells.size(), pTraits<Type>::zero);
    forAll(myCells, i)
    {
        const label globalCell = myCells[i];
        if (globalCell >= 0 && globalCell < f.size())
        {
            owned[i] = f[globalCell];
        }
    }
    f.transfer(owned);
}


// M4: temporarily expand an owned-size cumulative array to the full-mesh
// layout in place; the owned values are kept in `backup` for the restore.
template<class Type>
static void expandOwnedToFull
(
    Field<Type>& f,
    Field<Type>& backup,
    const UList<label>& myCells
)
{
    backup = f;
    Field<Type> full(myCells.size(), pTraits<Type>::zero);
    forAll(myCells, i)
    {
        full[myCells[i]] = backup[i];
    }
    f.transfer(full);
}


template<class Type>
static void expandOwnedListToFull
(
    List<Field<Type>>& fs,
    List<Field<Type>>& backups,
    const UList<label>& myCells
)
{
    backups.setSize(fs.size());
    forAll(fs, i)
    {
        expandOwnedToFull(fs[i], backups[i], myCells);
    }
}


template<class Type>
static void restoreOwnedList
(
    List<Field<Type>>& fs,
    const List<Field<Type>>& backups
)
{
    forAll(backups, i)
    {
        restoreField(fs[i], backups[i]);
    }
}


void dsmcVolFields::initOwnedStorage()
{
    const bool ownedMode =
        cloud_.replicatedMeshActive()
     && cloud_.replicatedMesh().processorWriteEnabled();

    if (!ownedMode)
    {
        if (pendingResumeReadIn_)
        {
            pendingResumeReadIn_ = false;
            readIn();
        }
        return;
    }

    const dsmcReplicatedMesh& replMesh = cloud_.replicatedMesh();
    const UList<label>& myCells = replMesh.myCells();
    const label ownerVersion =
        replMesh.rebalanceCount() + replMesh.autoRebalanceCount();

    if
    (
        ownedStorageOwnerVersion_ == ownerVersion
     && dsmcNCum_.size() == myCells.size()
     && !pendingResumeReadIn_
    )
    {
        return;
    }

    ownedStorageActive_ = true;
    ownedStorageOwnerVersion_ = ownerVersion;

    // Rebuild the global-cell -> owned-index mapping.
    toOwnedCell_.setSize(mesh_.nCells(), -1);
    forAll(myCells, i)
    {
        toOwnedCell_[myCells[i]] = i;
    }

    // First owned allocation (ctor storage is full-size and zeroed) or
    // re-allocation after a DLB owner reassignment: per-cell history cannot
    // be remapped onto the new owned layout, so the averaging window
    // restarts here (memory_opt.md M4 step 1, accepted cost).
    nTimeSteps_ = 0.0;

    dsmcNCum_ = scalarField(myCells.size(), 0.0);
    nCum_ = scalarField(myCells.size(), 0.0);
    dsmcNElecLvlCum_ = scalarField(myCells.size(), 0.0);
    dsmcMCum_ = scalarField(myCells.size(), 0.0);
    mCum_ = scalarField(myCells.size(), 0.0);
    dsmcLinearKECum_ = scalarField(myCells.size(), 0.0);
    linearKECum_ = scalarField(myCells.size(), 0.0);
    dsmcErotCum_ = scalarField(myCells.size(), 0.0);
    dsmcZetaRotCum_ = scalarField(myCells.size(), 0.0);
    dsmcMuuCum_ = scalarField(myCells.size(), 0.0);
    dsmcMuvCum_ = scalarField(myCells.size(), 0.0);
    dsmcMuwCum_ = scalarField(myCells.size(), 0.0);
    dsmcMvvCum_ = scalarField(myCells.size(), 0.0);
    dsmcMvwCum_ = scalarField(myCells.size(), 0.0);
    dsmcMwwCum_ = scalarField(myCells.size(), 0.0);
    dsmcMccCum_ = scalarField(myCells.size(), 0.0);
    dsmcMccuCum_ = scalarField(myCells.size(), 0.0);
    dsmcMccvCum_ = scalarField(myCells.size(), 0.0);
    dsmcMccwCum_ = scalarField(myCells.size(), 0.0);
    dsmcEuCum_ = scalarField(myCells.size(), 0.0);
    dsmcEvCum_ = scalarField(myCells.size(), 0.0);
    dsmcEwCum_ = scalarField(myCells.size(), 0.0);
    dsmcECum_ = scalarField(myCells.size(), 0.0);
    zetaVib_ = scalarField(myCells.size(), 0.0);
    dsmcNClassICum_ = scalarField(myCells.size(), 0.0);
    dsmcNClassIICum_ = scalarField(myCells.size(), 0.0);
    dsmcNClassIIICum_ = scalarField(myCells.size(), 0.0);
    collisionSeparation_ = scalarField(myCells.size(), 0.0);
    dsmcNCollsCum_ = scalarField(myCells.size(), 0.0);
    dsmcMomentumCum_ = vectorField(myCells.size(), vector::zero);
    momentumCum_ = vectorField(myCells.size(), vector::zero);

    forAll(dsmcNSpeciesCum_, i)
    {
        dsmcNSpeciesCum_[i] = scalarField(myCells.size(), 0.0);
        nSpeciesCum_[i] = scalarField(myCells.size(), 0.0);
        dsmcMccSpeciesCum_[i] = scalarField(myCells.size(), 0.0);
        speciesMfp_[i] = scalarField(myCells.size(), 0.0);
        speciesMcr_[i] = scalarField(myCells.size(), 0.0);
        speciesTvib_[i] = scalarField(myCells.size(), 0.0);
        dsmcSpeciesEelecCum_[i] = scalarField(myCells.size(), 0.0);
        dsmcNGrndElecLvlSpeciesCum_[i] = scalarField(myCells.size(), 0.0);
        dsmcN1stElecLvlSpeciesCum_[i] = scalarField(myCells.size(), 0.0);

        forAll(dsmcSpeciesEvibModCum_[i], mod)
        {
            dsmcSpeciesEvibModCum_[i][mod] =
                scalarField(myCells.size(), 0.0);
        }
    }

    if (pendingResumeReadIn_)
    {
        pendingResumeReadIn_ = false;
        readIn();

        // A resume file written by an earlier full-size layout holds
        // full-mesh lists: keep only the owned entries.
        if (dsmcNCum_.size() != myCells.size())
        {
            gatherOwnedEntries(dsmcNCum_, myCells);
            gatherOwnedEntries(nCum_, myCells);
            gatherOwnedEntries(dsmcNElecLvlCum_, myCells);
            gatherOwnedEntries(dsmcMCum_, myCells);
            gatherOwnedEntries(mCum_, myCells);
            gatherOwnedEntries(dsmcLinearKECum_, myCells);
            gatherOwnedEntries(linearKECum_, myCells);
            gatherOwnedEntries(dsmcErotCum_, myCells);
            gatherOwnedEntries(dsmcZetaRotCum_, myCells);
            gatherOwnedEntries(dsmcMuuCum_, myCells);
            gatherOwnedEntries(dsmcMuvCum_, myCells);
            gatherOwnedEntries(dsmcMuwCum_, myCells);
            gatherOwnedEntries(dsmcMvvCum_, myCells);
            gatherOwnedEntries(dsmcMvwCum_, myCells);
            gatherOwnedEntries(dsmcMwwCum_, myCells);
            gatherOwnedEntries(dsmcMccCum_, myCells);
            gatherOwnedEntries(dsmcMccuCum_, myCells);
            gatherOwnedEntries(dsmcMccvCum_, myCells);
            gatherOwnedEntries(dsmcMccwCum_, myCells);
            gatherOwnedEntries(dsmcEuCum_, myCells);
            gatherOwnedEntries(dsmcEvCum_, myCells);
            gatherOwnedEntries(dsmcEwCum_, myCells);
            gatherOwnedEntries(dsmcECum_, myCells);
            gatherOwnedEntries(zetaVib_, myCells);
            gatherOwnedEntries(dsmcNClassICum_, myCells);
            gatherOwnedEntries(dsmcNClassIICum_, myCells);
            gatherOwnedEntries(dsmcNClassIIICum_, myCells);
            gatherOwnedEntries(collisionSeparation_, myCells);
            gatherOwnedEntries(dsmcNCollsCum_, myCells);
            gatherOwnedEntries(dsmcMomentumCum_, myCells);
            gatherOwnedEntries(momentumCum_, myCells);

            forAll(dsmcNSpeciesCum_, i)
            {
                gatherOwnedEntries(dsmcNSpeciesCum_[i], myCells);
                gatherOwnedEntries(nSpeciesCum_[i], myCells);
                gatherOwnedEntries(dsmcMccSpeciesCum_[i], myCells);
                gatherOwnedEntries(speciesMfp_[i], myCells);
                gatherOwnedEntries(speciesMcr_[i], myCells);
                gatherOwnedEntries(speciesTvib_[i], myCells);
                gatherOwnedEntries(dsmcSpeciesEelecCum_[i], myCells);
                gatherOwnedEntries(dsmcNGrndElecLvlSpeciesCum_[i], myCells);
                gatherOwnedEntries(dsmcN1stElecLvlSpeciesCum_[i], myCells);

                forAll(dsmcSpeciesEvibModCum_[i], mod)
                {
                    gatherOwnedEntries
                    (
                        dsmcSpeciesEvibModCum_[i][mod],
                        myCells
                    );
                }
            }
        }
    }
}


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

        // M4: the resume file keeps the full-mesh layout regardless of the
        // storage mode.  In owned mode, temporarily expand the cumulative
        // arrays in place, write, then restore.
        bool expandedOwned = false;
        const UList<label>* writeCellsPtr = nullptr;
        Field<scalar> bDsmcNCum, bNCum, bDsmcNElecLvlCum, bDsmcMCum, bMCum,
            bDsmcLinearKECum, bLinearKECum, bDsmcErotCum, bDsmcZetaRotCum,
            bDsmcMuuCum, bDsmcMuvCum, bDsmcMuwCum, bDsmcMvvCum, bDsmcMvwCum,
            bDsmcMwwCum, bDsmcMccCum, bDsmcMccuCum, bDsmcMccvCum, bDsmcMccwCum,
            bDsmcEuCum, bDsmcEvCum, bDsmcEwCum, bDsmcECum, bZetaVib,
            bDsmcNClassICum, bDsmcNClassIICum, bDsmcNClassIIICum,
            bCollisionSeparation, bDsmcNCollsCum;
        vectorField bDsmcMomentumCum, bMomentumCum;
        List<Field<scalar>> bDsmcSpeciesEelecCum, bDsmcNSpeciesCum,
            bNSpeciesCum, bDsmcMccSpeciesCum, bDsmcNGrndElecLvlSpeciesCum,
            bDsmcN1stElecLvlSpeciesCum;
        List<List<Field<scalar>>> bDsmcSpeciesEvibModCum;

        if (ownedStorageActive_)
        {
            expandedOwned = true;
            writeCellsPtr = &cloud_.replicatedMesh().myCells();
            const UList<label>& writeCells = *writeCellsPtr;

            expandOwnedToFull(dsmcNCum_, bDsmcNCum, writeCells);
            expandOwnedToFull(nCum_, bNCum, writeCells);
            expandOwnedToFull(dsmcNElecLvlCum_, bDsmcNElecLvlCum, writeCells);
            expandOwnedToFull(dsmcMCum_, bDsmcMCum, writeCells);
            expandOwnedToFull(mCum_, bMCum, writeCells);
            expandOwnedToFull(dsmcLinearKECum_, bDsmcLinearKECum, writeCells);
            expandOwnedToFull(linearKECum_, bLinearKECum, writeCells);
            expandOwnedToFull(dsmcErotCum_, bDsmcErotCum, writeCells);
            expandOwnedToFull(dsmcZetaRotCum_, bDsmcZetaRotCum, writeCells);
            expandOwnedToFull(dsmcMuuCum_, bDsmcMuuCum, writeCells);
            expandOwnedToFull(dsmcMuvCum_, bDsmcMuvCum, writeCells);
            expandOwnedToFull(dsmcMuwCum_, bDsmcMuwCum, writeCells);
            expandOwnedToFull(dsmcMvvCum_, bDsmcMvvCum, writeCells);
            expandOwnedToFull(dsmcMvwCum_, bDsmcMvwCum, writeCells);
            expandOwnedToFull(dsmcMwwCum_, bDsmcMwwCum, writeCells);
            expandOwnedToFull(dsmcMccCum_, bDsmcMccCum, writeCells);
            expandOwnedToFull(dsmcMccuCum_, bDsmcMccuCum, writeCells);
            expandOwnedToFull(dsmcMccvCum_, bDsmcMccvCum, writeCells);
            expandOwnedToFull(dsmcMccwCum_, bDsmcMccwCum, writeCells);
            expandOwnedToFull(dsmcEuCum_, bDsmcEuCum, writeCells);
            expandOwnedToFull(dsmcEvCum_, bDsmcEvCum, writeCells);
            expandOwnedToFull(dsmcEwCum_, bDsmcEwCum, writeCells);
            expandOwnedToFull(dsmcECum_, bDsmcECum, writeCells);
            expandOwnedToFull(zetaVib_, bZetaVib, writeCells);
            expandOwnedToFull(dsmcNClassICum_, bDsmcNClassICum, writeCells);
            expandOwnedToFull(dsmcNClassIICum_, bDsmcNClassIICum, writeCells);
            expandOwnedToFull(dsmcNClassIIICum_, bDsmcNClassIIICum, writeCells);
            expandOwnedToFull(collisionSeparation_, bCollisionSeparation, writeCells);
            expandOwnedToFull(dsmcNCollsCum_, bDsmcNCollsCum, writeCells);
            expandOwnedToFull(dsmcMomentumCum_, bDsmcMomentumCum, writeCells);
            expandOwnedToFull(momentumCum_, bMomentumCum, writeCells);

            expandOwnedListToFull(dsmcSpeciesEelecCum_, bDsmcSpeciesEelecCum, writeCells);
            expandOwnedListToFull(dsmcNSpeciesCum_, bDsmcNSpeciesCum, writeCells);
            expandOwnedListToFull(nSpeciesCum_, bNSpeciesCum, writeCells);
            expandOwnedListToFull(dsmcMccSpeciesCum_, bDsmcMccSpeciesCum, writeCells);
            expandOwnedListToFull(dsmcNGrndElecLvlSpeciesCum_, bDsmcNGrndElecLvlSpeciesCum, writeCells);
            expandOwnedListToFull(dsmcN1stElecLvlSpeciesCum_, bDsmcN1stElecLvlSpeciesCum, writeCells);

            bDsmcSpeciesEvibModCum.setSize(dsmcSpeciesEvibModCum_.size());
            forAll(dsmcSpeciesEvibModCum_, i)
            {
                bDsmcSpeciesEvibModCum[i].setSize
                (
                    dsmcSpeciesEvibModCum_[i].size()
                );
                forAll(dsmcSpeciesEvibModCum_[i], mod)
                {
                    expandOwnedToFull
                    (
                        dsmcSpeciesEvibModCum_[i][mod],
                        bDsmcSpeciesEvibModCum[i][mod],
                        writeCells
                    );
                }
            }
        }

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

        if (expandedOwned)
        {
            restoreField(dsmcNCum_, bDsmcNCum);
            restoreField(nCum_, bNCum);
            restoreField(dsmcNElecLvlCum_, bDsmcNElecLvlCum);
            restoreField(dsmcMCum_, bDsmcMCum);
            restoreField(mCum_, bMCum);
            restoreField(dsmcLinearKECum_, bDsmcLinearKECum);
            restoreField(linearKECum_, bLinearKECum);
            restoreField(dsmcErotCum_, bDsmcErotCum);
            restoreField(dsmcZetaRotCum_, bDsmcZetaRotCum);
            restoreField(dsmcMuuCum_, bDsmcMuuCum);
            restoreField(dsmcMuvCum_, bDsmcMuvCum);
            restoreField(dsmcMuwCum_, bDsmcMuwCum);
            restoreField(dsmcMvvCum_, bDsmcMvvCum);
            restoreField(dsmcMvwCum_, bDsmcMvwCum);
            restoreField(dsmcMwwCum_, bDsmcMwwCum);
            restoreField(dsmcMccCum_, bDsmcMccCum);
            restoreField(dsmcMccuCum_, bDsmcMccuCum);
            restoreField(dsmcMccvCum_, bDsmcMccvCum);
            restoreField(dsmcMccwCum_, bDsmcMccwCum);
            restoreField(dsmcEuCum_, bDsmcEuCum);
            restoreField(dsmcEvCum_, bDsmcEvCum);
            restoreField(dsmcEwCum_, bDsmcEwCum);
            restoreField(dsmcECum_, bDsmcECum);
            restoreField(zetaVib_, bZetaVib);
            restoreField(dsmcNClassICum_, bDsmcNClassICum);
            restoreField(dsmcNClassIICum_, bDsmcNClassIICum);
            restoreField(dsmcNClassIIICum_, bDsmcNClassIIICum);
            restoreField(collisionSeparation_, bCollisionSeparation);
            restoreField(dsmcNCollsCum_, bDsmcNCollsCum);
            restoreField(dsmcMomentumCum_, bDsmcMomentumCum);
            restoreField(momentumCum_, bMomentumCum);

            restoreOwnedList(dsmcSpeciesEelecCum_, bDsmcSpeciesEelecCum);
            restoreOwnedList(dsmcNSpeciesCum_, bDsmcNSpeciesCum);
            restoreOwnedList(nSpeciesCum_, bNSpeciesCum);
            restoreOwnedList(dsmcMccSpeciesCum_, bDsmcMccSpeciesCum);
            restoreOwnedList(dsmcNGrndElecLvlSpeciesCum_, bDsmcNGrndElecLvlSpeciesCum);
            restoreOwnedList(dsmcN1stElecLvlSpeciesCum_, bDsmcN1stElecLvlSpeciesCum);

            forAll(dsmcSpeciesEvibModCum_, i)
            {
                forAll(dsmcSpeciesEvibModCum_[i], mod)
                {
                    restoreField
                    (
                        dsmcSpeciesEvibModCum_[i][mod],
                        bDsmcSpeciesEvibModCum[i][mod]
                    );
                }
            }
        }
    }
}


//- Initial configuration
// M2: construct the GeoFields on demand.  Core fields (written by the
// derive pass unconditionally or consumed through the virtual accessors)
// are always built; the rest follow their measure/write switches so that
// disabled families never pay their memory.
void dsmcVolFields::constructOutputFields()
{
    // Core set: derive pass writes these unconditionally; Tov_ backs the
    // virtual overallT() accessor used by the collision model.
    #define M2_MAKE_FIELD(NAME, TYPE, DIMS)                          \
        if (NAME##_.empty())                                         \
        {                                                            \
            NAME##_.reset                                            \
            (                                                        \
                new TYPE                                             \
                (                                                    \
                    IOobject                                         \
                    (                                                \
                        #NAME + fieldName_,                          \
                        time_.time().timeName(),                     \
                        mesh_,                                       \
                        IOobject::NO_READ,                           \
                        IOobject::AUTO_WRITE                         \
                    ),                                               \
                    mesh_,                                           \
                    dimensionedScalar("0.0", DIMS, 0.0)              \
                )                                                    \
            );                                                       \
        }

    M2_MAKE_FIELD(dsmcN, volScalarField, dimless)
    M2_MAKE_FIELD(dsmcNMean, volScalarField, dimless)
    M2_MAKE_FIELD(rhoN, volScalarField, dimless/dimVolume)
    M2_MAKE_FIELD(rhoM, volScalarField, dimMass/dimVolume)
    M2_MAKE_FIELD(p, volScalarField, dimPressure)
    M2_MAKE_FIELD(Ttra, volScalarField, dimTemperature)
    M2_MAKE_FIELD(Trot, volScalarField, dimTemperature)
    M2_MAKE_FIELD(Tvib, volScalarField, dimTemperature)
    M2_MAKE_FIELD(Telec, volScalarField, dimTemperature)
    M2_MAKE_FIELD(Ma, volScalarField, dimless)
    M2_MAKE_FIELD(q, volScalarField, dimensionSet(1, 0, -3, 0, 0))
    if (fD_.empty())
    {
        fD_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "fD_"+ fieldName_,
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "0.0",
                    dimensionSet(1, -1, -2, 0, 0),
                    vector::zero
                )
            )
        );
    }
    M2_MAKE_FIELD(tau, volScalarField, dimPressure)
    if (UMean_.empty())
    {
        UMean_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "UMean_"+ fieldName_,
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector("0.0", dimensionSet(0, 1, -1, 0, 0), vector::zero)
            )
        );
    }
    M2_MAKE_FIELD(Tov, volScalarField, dimTemperature)

    #undef M2_MAKE_FIELD

}


void dsmcVolFields::createField()
{
    if (Foam::dsmcIsPrintingRank())
    {
        Info << "Initialising dsmcVolFields field" << endl;
    }

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

    restrictOutputFields_ = false;
    outputFieldNames_.clear();
    if (propsDict_.found("writeFields"))
    {
        wordList requestedOutputFields(propsDict_.lookup("writeFields"));
        restrictOutputFields_ = true;
        forAll(requestedOutputFields, i)
        {
            outputFieldNames_.insert(requestedOutputFields[i]);
        }

        Info<< "dsmcVolFields [" << fieldName_
            << "]: restricting output fields to "
            << requestedOutputFields << endl;
    }

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
        mesh_.time().controlDict().lookupOrDefault<scalar>("iMeanMinTvib", 0.01);

    averagingAcrossManyRuns_ =
        propsDict_.lookupOrDefault<bool>("averagingAcrossManyRuns", false);

    //- read in stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        if (!time_.resetFieldsAtOutput())
        {
            // M4: the load is deferred to initOwnedStorage() (first
            // sampling call) so that, in owned mode, the resume file can
            // be gathered onto the owned-cell layout.
            pendingResumeReadIn_ = true;
            Info<< "Averaging across many runs for field " << fieldName_
                << " is enabled. Sampled data will be read from file at the"
                << " first sampling step."
                << endl;
        }
        else
        {
            Info<< "Averaging across many runs for field " << fieldName_
                << " will be enabled as soon as resetAtOutput is turned off."
                << endl;
        }
    }

    updateOwnedBoundaryFaces();

    // M2: GeoFields are built here, after all gate flags are known.
    constructOutputFields();
}


void dsmcVolFields::updateOwnedBoundaryFaces()
{
    if (!cloud_.replicatedMeshActive())
    {
        ownedBoundaryFaces_.clear();
        ownedBoundaryFacesOwnerVersion_ = -1;
        return;
    }

    const dsmcReplicatedMesh& replicatedMesh = cloud_.replicatedMesh();
    const label ownerVersion =
        replicatedMesh.rebalanceCount() + replicatedMesh.autoRebalanceCount();

    if
    (
        ownedBoundaryFacesOwnerVersion_ == ownerVersion
     && ownedBoundaryFaces_.size() == mesh_.boundaryMesh().size()
    )
    {
        return;
    }

    ownedBoundaryFaces_.setSize(mesh_.boundaryMesh().size());

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const label startFace = pp.start();
        DynamicList<label> owned(pp.size()/8 + 1);

        forAll(pp, facei)
        {
            if (replicatedMesh.isMyCell(mesh_.faceOwner()[startFace + facei]))
            {
                owned.append(facei);
            }
        }

        ownedBoundaryFaces_[patchi].transfer(owned);
    }

    ownedBoundaryFacesOwnerVersion_ = ownerVersion;
}


void dsmcVolFields::calculateField()
{
    sampleCounter_++;

    updateOwnedBoundaryFaces();

    // M4: ensure the cumulative storage matches the current owned-cell set
    // (first call: allocate owned-size and apply any deferred resume load;
    // after a DLB owner reassignment: restart the averaging window).
    initOwnedStorage();

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
    
    // dsmcN is an instantaneous field.  In a replicated mesh, a cell can
    // change owner between samples, so an active-cell-only reset can expose
    // an older value left on its new owner.
    dsmcN_() = 0.0;

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
            // M4: sample over all owned cells (ascending order) rather than
            // occupancyActiveCells so the same cell list is used by every
            // pass; empty cells contribute nothing and keep the RNG stream
            // untouched.
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

                    dsmcNCum_[toOwned(cell)] += dsmcNLocal;
                    dsmcN_()[cell] += dsmcNLocal;
                    nCum_[toOwned(cell)] += nLocal;
                    mCum_[toOwned(cell)] += mLocal;
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
                        const scalar nParticles = cloud_.nParticles(cell);
                        const scalar mass = cloud_.constProps(typeId).mass();

                        // cumulative number of DSMC parcels
                        dsmcNCum_[toOwned(cell)] += 1.0;
                        // instantaneous number of DSMC parcels in this time step
                        dsmcN_()[cell] += 1.0;
                        // cumulative number of real particles
                        nCum_[toOwned(cell)] += nParticles;
                        // cumulative mass of real particles
                        mCum_[toOwned(cell)] += mass*nParticles;
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
	            bool allSpeciesField =
	                !singleSpeciesField
	             && speciesIds_.size() == cloud_.constProps().size();
	            if (allSpeciesField)
	            {
	                for (label typeId = 0; typeId < typeIdToSpeciesIndex_.size(); ++typeId)
	                {
	                    if (typeIdToSpeciesIndex_[typeId] < 0)
	                    {
	                        allSpeciesField = false;
	                        break;
	                    }
	                }
	            }
            // M4: combine over all owned cells (ascending order) rather than
            // occupancyActiveCells so the same cell list is used by every
            // pass; empty cells contribute nothing and keep the RNG stream
            // untouched.
            const UList<label>* combineCellsPtr =
                cloud_.replicatedMeshActive()
              ? static_cast<const UList<label>*>(&cloud_.replicatedMesh().myCells())
              : nullptr;
            const label combineLoopSize =
                combineCellsPtr ? combineCellsPtr->size() : dsmcNCum_.size();
            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) if (useOpenMPSampling)
            #endif
            for (label combineI = 0; combineI < combineLoopSize; ++combineI)
            {
                const label cell =
                    combineCellsPtr ? (*combineCellsPtr)[combineI] : combineI;
                const scalar cellNParticles = cloud_.nParticles(cell);
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
                    nLocal = dsmcNLocal*cellNParticles;
                    mLocal = dsmcMLocal*cellNParticles;
                    momentumLocal = dsmcMomentumLocal*cellNParticles;
                    linearKELocal = dsmcLinearKELocal*cellNParticles;

                    dsmcNSpeciesCum_[0][toOwned(cell)] += dsmcNLocal;
                    dsmcMccSpeciesCum_[0][toOwned(cell)] += dsmcLinearKELocal;
                    nSpeciesCum_[0][toOwned(cell)] += nLocal;

                    if (needVibrational)
                    {
                        forAll(dsmcSpeciesEvibModCum_[0], mod)
                        {
                            dsmcSpeciesEvibModCum_[0][mod][toOwned(cell)] +=
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
	                    if (allSpeciesField)
	                    {
	                        dsmcNLocal = sharedSampleCache_.totalDsmcN[cell];
	                        dsmcMLocal = sharedSampleCache_.totalDsmcM[cell];
	                        dsmcLinearKELocal =
	                            sharedSampleCache_.totalDsmcLinearKE[cell];
	                        dsmcMomentumLocal =
	                            sharedSampleCache_.totalDsmcMomentum[cell];
	                        dsmcErotLocal = sharedSampleCache_.totalDsmcErot[cell];
	                        dsmcZetaRotLocal =
	                            sharedSampleCache_.totalDsmcZetaRot[cell];
	                        nLocal = dsmcNLocal*cellNParticles;
	                        mLocal = dsmcMLocal*cellNParticles;
	                        momentumLocal = dsmcMomentumLocal*cellNParticles;
	                        linearKELocal = dsmcLinearKELocal*cellNParticles;
	                    }

	                    forAll(speciesIds_, i)
	                    {
	                        const label typeId = speciesIds_[i];

	                        const scalar speciesDsmcN =
	                            sharedSampleCache_.dsmcN[typeId][cell];
	                        const scalar speciesDsmcLinearKE =
	                            sharedSampleCache_.dsmcLinearKE[typeId][cell];
	                        const scalar speciesNReal =
	                            speciesDsmcN*cellNParticles;

	                        if (!allSpeciesField)
	                        {
	                            dsmcNLocal += speciesDsmcN;
	                            dsmcMLocal +=
	                                sharedSampleCache_.dsmcM[typeId][cell];
	                            dsmcLinearKELocal += speciesDsmcLinearKE;
	                            dsmcMomentumLocal +=
	                                sharedSampleCache_.dsmcMomentum[typeId][cell];
	                            dsmcErotLocal +=
	                                sharedSampleCache_.dsmcErot[typeId][cell];
	                            dsmcZetaRotLocal +=
	                                sharedSampleCache_.dsmcZetaRot[typeId][cell];
	                            nLocal += speciesNReal;
	                            mLocal +=
	                                sharedSampleCache_.dsmcM[typeId][cell]
	                               *cellNParticles;
	                            momentumLocal +=
	                                sharedSampleCache_.dsmcMomentum[typeId][cell]
	                               *cellNParticles;
	                            linearKELocal +=
	                                speciesDsmcLinearKE*cellNParticles;
	                        }

                        dsmcNSpeciesCum_[i][toOwned(cell)] += speciesDsmcN;
                        dsmcMccSpeciesCum_[i][toOwned(cell)] += speciesDsmcLinearKE;
                        nSpeciesCum_[i][toOwned(cell)] += speciesNReal;

                        if (needVibrational)
                        {
                            forAll(dsmcSpeciesEvibModCum_[i], mod)
                            {
                                dsmcSpeciesEvibModCum_[i][mod][toOwned(cell)] +=
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

                dsmcNCum_[toOwned(cell)] += dsmcNLocal;
                dsmcN_()[cell] += dsmcNLocal;
                dsmcMCum_[toOwned(cell)] += dsmcMLocal;
                dsmcLinearKECum_[toOwned(cell)] += dsmcLinearKELocal;
                dsmcMomentumCum_[toOwned(cell)] += dsmcMomentumLocal;
                dsmcErotCum_[toOwned(cell)] += dsmcErotLocal;
                dsmcZetaRotCum_[toOwned(cell)] += dsmcZetaRotLocal;
                dsmcNElecLvlCum_[toOwned(cell)] += dsmcNElecLvlLocal;
                nCum_[toOwned(cell)] += nLocal;
                mCum_[toOwned(cell)] += mLocal;
                momentumCum_[toOwned(cell)] += momentumLocal;
                linearKECum_[toOwned(cell)] += linearKELocal;

                if (needHeatFluxShearStress)
                {
                    dsmcMuuCum_[toOwned(cell)] += dsmcMuuLocal;
                    dsmcMuvCum_[toOwned(cell)] += dsmcMuvLocal;
                    dsmcMuwCum_[toOwned(cell)] += dsmcMuwLocal;
                    dsmcMvvCum_[toOwned(cell)] += dsmcMvvLocal;
                    dsmcMvwCum_[toOwned(cell)] += dsmcMvwLocal;
                    dsmcMwwCum_[toOwned(cell)] += dsmcMwwLocal;
                    dsmcMccCum_[toOwned(cell)] += dsmcMccLocal;
                    dsmcMccuCum_[toOwned(cell)] += dsmcMccuLocal;
                    dsmcMccvCum_[toOwned(cell)] += dsmcMccvLocal;
                    dsmcMccwCum_[toOwned(cell)] += dsmcMccwLocal;
                    dsmcEuCum_[toOwned(cell)] += dsmcEuLocal;
                    dsmcEvCum_[toOwned(cell)] += dsmcEvLocal;
                    dsmcEwCum_[toOwned(cell)] += dsmcEwLocal;
                    dsmcECum_[toOwned(cell)] += dsmcECumLocal;
                }

                if (needClassification)
                {
                    dsmcNClassICum_[toOwned(cell)] += dsmcNClassILocal;
                    dsmcNClassIICum_[toOwned(cell)] += dsmcNClassIILocal;
                    dsmcNClassIIICum_[toOwned(cell)] += dsmcNClassIIILocal;
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
                collisionSeparation_[toOwned(celli)] +=
                    cloud_.cellPropMeasurements().collisionSeparation()[celli];
                    
                dsmcNCollsCum_[toOwned(celli)] +=
                    cloud_.cellPropMeasurements().nColls()[celli];

                if (dsmcNCum_[toOwned(celli)] > 1e-3)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[celli];

                    dsmcNMean_()[celli] = dsmcNCum_[toOwned(celli)]/nAvTimeSteps;

                    const scalar rhoNMean = nCum_[toOwned(celli)]
                        /(nAvTimeSteps*cellVolume);
                    const scalar rhoMMean = mCum_[toOwned(celli)]
                        /(nAvTimeSteps*cellVolume);

                    rhoN_()[celli] = rhoNMean;
                    rhoM_()[celli] = rhoMMean;
                    
                    UMean_()[celli] = momentumCum_[toOwned(celli)]/mCum_[toOwned(celli)];

                    const scalar linearKEMean = 0.5*linearKECum_[toOwned(celli)]
                        /(cellVolume*nAvTimeSteps);

                    Ttra_()[celli] =
                        2.0/(3.0*kB*rhoNMean)
                       *(
                            linearKEMean - 0.5*rhoMMean
                           *(
                                UMean_()[celli] & UMean_()[celli]
                            )
                        );

                    p_()[celli] = rhoNMean*kB*Ttra_()[celli];
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcNMean_()[celli] = 0.001;
                    rhoN_()[celli] = 0.0;
                    rhoM_()[celli] = 0.0;
                    UMean_()[celli] = vector::zero;
                    Ttra_()[celli] = 0.0;
                    p_()[celli] = 0.0;
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
            // Wall hits are recorded on the rank where the move occurs;
            // migration to the cell owner happens only after the move.
            // Accumulate every local face here and apply ownership only when
            // processor output is selected below.
            forAll(speciesIds_, i)
            {
                const label spId = speciesIds_[i];

                // Boundary accumulation is per-face disjoint (each face's
                // slots are written only here) and the boundaryFlux
                // accessors are const reads — OMP over the sampled patches
                // (§8.6).
                const label nSampledPatches = sampledBoundaryPatches_.size();
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic)
                #endif
                for (label patchi = 0; patchi < nSampledPatches; ++patchi)
                {
                    const label j = sampledBoundaryPatches_[patchi];
                    const label nFaces = mesh_.boundaryMesh()[j].size();
                    for (label facei = 0; facei < nFaces; ++facei)
                    {
                        const label k = facei;
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
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic)
                    #endif
                    for (label patchi = 0; patchi < nSampledPatches; ++patchi)
                    {
                        const label j = sampledBoundaryPatches_[patchi];
                        const label nFaces = mesh_.boundaryMesh()[j].size();
                        for (label facei = 0; facei < nFaces; ++facei)
                        {
                            speciesEvibModBF_[i][mod][j][facei] +=
                                boundaryFlux.speciesEvibModBF(spId, mod, j, facei);
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
        const auto outputTimeStart =
            doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
        const scalar nAvTimeSteps = nTimeSteps_;

        if (nAvTimeSteps <= SMALL)
        {
            if (cloud_.isOutputRank())
            {
                WarningInFunction
                    << "No valid DSMC samples are available at output time "
                    << time_.time().value()
                    << ". Skipping dsmcVolFields output calculation for field "
                    << fieldName_ << '.' << endl;
            }

            return;
        }

        const bool processorWrite =
            cloud_.replicatedMeshActive()
         && cloud_.replicatedMesh().processorWriteEnabled();
        const bool computeOutputFields =
            !cloud_.replicatedMeshActive()
         || cloud_.isOutputRank()
         || processorWrite;
        const bool restoreLocalAfterOutput =
            cloud_.replicatedMeshActive()
         && !time_.resetFieldsAtOutput();
        // DLB leaves sampled history on prior ranks when cell ownership moves.
        // Sum it before applying the current output owner map in every mode.
        const bool reduceForOutput = cloud_.replicatedMeshActive();
        const UList<label>* outputCellsPtr =
            (
                cloud_.replicatedMeshActive()
             && processorWrite
            )
          ? static_cast<const UList<label>*>(&cloud_.replicatedMesh().myCells())
          : nullptr;
        const label outputLoopSize =
            outputCellsPtr ? outputCellsPtr->size() : dsmcNCum_.size();

        scalarField localDsmcNCum;
        scalarField localNCum;
        scalarField localDsmcNElecLvlCum;
        scalarField localDsmcMCum;
        scalarField localMCum;
        scalarField localDsmcLinearKECum;
        scalarField localLinearKECum;
        scalarField localDsmcErotCum;
        scalarField localDsmcZetaRotCum;
        scalarField localDsmcMuuCum;
        scalarField localDsmcMuvCum;
        scalarField localDsmcMuwCum;
        scalarField localDsmcMvvCum;
        scalarField localDsmcMvwCum;
        scalarField localDsmcMwwCum;
        scalarField localDsmcMccCum;
        scalarField localDsmcMccuCum;
        scalarField localDsmcMccvCum;
        scalarField localDsmcMccwCum;
        scalarField localDsmcEuCum;
        scalarField localDsmcEvCum;
        scalarField localDsmcEwCum;
        scalarField localDsmcECum;
        scalarField localZetaVib;
        scalarField localDsmcNClassICum;
        scalarField localDsmcNClassIICum;
        scalarField localDsmcNClassIIICum;
        scalarField localCollisionSeparation;
        scalarField localDsmcNCollsCum;
        vectorField localDsmcMomentumCum;
        vectorField localMomentumCum;

        List<scalarField> localDsmcSpeciesEelecCum;
        List<scalarField> localDsmcNSpeciesCum;
        List<scalarField> localNSpeciesCum;
        List<scalarField> localDsmcMccSpeciesCum;
        List<scalarField> localDsmcNGrndElecLvlSpeciesCum;
        List<scalarField> localDsmcN1stElecLvlSpeciesCum;
        List<List<scalarField>> localDsmcSpeciesEvibModCum;

        List<scalarField> localRhoNBF;
        List<scalarField> localRhoMBF;
        List<scalarField> localLinearKEBF;
        List<scalarField> localErotBF;
        List<scalarField> localZetaRotBF;
        List<scalarField> localQBF;
        List<scalarField> localZetaVibBF;
        List<scalarField> localRhoNIntBF;
        List<scalarField> localRhoNElecBF;
        List<vectorField> localMomentumBF;
        List<vectorField> localFDBF;

        List<List<scalarField>> localSpeciesEvibBF;
        List<List<scalarField>> localSpeciesEelecBF;
        List<List<scalarField>> localSpeciesRhoNBF;
        List<List<scalarField>> localSpeciesMccBF;
        List<List<List<scalarField>>> localSpeciesEvibModBF;

        if (restoreLocalAfterOutput)
        {
            const auto outputSnapshotStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            snapshotField(dsmcNCum_, localDsmcNCum);
            snapshotField(nCum_, localNCum);
            snapshotField(dsmcNElecLvlCum_, localDsmcNElecLvlCum);
            snapshotField(dsmcMCum_, localDsmcMCum);
            snapshotField(mCum_, localMCum);
            snapshotField(dsmcLinearKECum_, localDsmcLinearKECum);
            snapshotField(linearKECum_, localLinearKECum);
            snapshotField(dsmcErotCum_, localDsmcErotCum);
            snapshotField(dsmcZetaRotCum_, localDsmcZetaRotCum);
            snapshotField(dsmcMuuCum_, localDsmcMuuCum);
            snapshotField(dsmcMuvCum_, localDsmcMuvCum);
            snapshotField(dsmcMuwCum_, localDsmcMuwCum);
            snapshotField(dsmcMvvCum_, localDsmcMvvCum);
            snapshotField(dsmcMvwCum_, localDsmcMvwCum);
            snapshotField(dsmcMwwCum_, localDsmcMwwCum);
            snapshotField(dsmcMccCum_, localDsmcMccCum);
            snapshotField(dsmcMccuCum_, localDsmcMccuCum);
            snapshotField(dsmcMccvCum_, localDsmcMccvCum);
            snapshotField(dsmcMccwCum_, localDsmcMccwCum);
            snapshotField(dsmcEuCum_, localDsmcEuCum);
            snapshotField(dsmcEvCum_, localDsmcEvCum);
            snapshotField(dsmcEwCum_, localDsmcEwCum);
            snapshotField(dsmcECum_, localDsmcECum);
            snapshotField(zetaVib_, localZetaVib);
            snapshotField(dsmcNClassICum_, localDsmcNClassICum);
            snapshotField(dsmcNClassIICum_, localDsmcNClassIICum);
            snapshotField(dsmcNClassIIICum_, localDsmcNClassIIICum);
            snapshotField(collisionSeparation_, localCollisionSeparation);
            snapshotField(dsmcNCollsCum_, localDsmcNCollsCum);
            snapshotField(dsmcMomentumCum_, localDsmcMomentumCum);
            snapshotField(momentumCum_, localMomentumCum);

            snapshotFieldList(dsmcSpeciesEelecCum_, localDsmcSpeciesEelecCum);
            snapshotFieldList(dsmcNSpeciesCum_, localDsmcNSpeciesCum);
            snapshotFieldList(nSpeciesCum_, localNSpeciesCum);
            snapshotFieldList(dsmcMccSpeciesCum_, localDsmcMccSpeciesCum);
            snapshotFieldList
            (
                dsmcNGrndElecLvlSpeciesCum_,
                localDsmcNGrndElecLvlSpeciesCum
            );
            snapshotFieldList
            (
                dsmcN1stElecLvlSpeciesCum_,
                localDsmcN1stElecLvlSpeciesCum
            );
            snapshotFieldListList
            (
                dsmcSpeciesEvibModCum_,
                localDsmcSpeciesEvibModCum
            );

            snapshotFieldList(rhoNBF_, localRhoNBF);
            snapshotFieldList(rhoMBF_, localRhoMBF);
            snapshotFieldList(linearKEBF_, localLinearKEBF);
            snapshotFieldList(ErotBF_, localErotBF);
            snapshotFieldList(zetaRotBF_, localZetaRotBF);
            snapshotFieldList(qBF_, localQBF);
            snapshotFieldList(zetaVibBF_, localZetaVibBF);
            snapshotFieldList(rhoNIntBF_, localRhoNIntBF);
            snapshotFieldList(rhoNElecBF_, localRhoNElecBF);
            snapshotFieldList(momentumBF_, localMomentumBF);
            snapshotFieldList(fDBF_, localFDBF);

            snapshotFieldListList(speciesEvibBF_, localSpeciesEvibBF);
            snapshotFieldListList(speciesEelecBF_, localSpeciesEelecBF);
            snapshotFieldListList(speciesRhoNBF_, localSpeciesRhoNBF);
            snapshotFieldListList(speciesMccBF_, localSpeciesMccBF);
            snapshotFieldListListList
            (
                speciesEvibModBF_,
                localSpeciesEvibModBF
            );

            if (doProfile)
            {
                profileOutputSnapshotWallTime_ +=
                    wallSeconds(outputSnapshotStart, wallClockNow());
            }
        }

        if (reduceForOutput)
        {
            const auto outputReduceStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            sumReduceField(dsmcN_().primitiveFieldRef());
            // M4: in owned mode the per-cell cumulative arrays are
            // owned-size and each cell has exactly one owner rank, so the
            // cross-rank sum is unnecessary (and would mismatch sizes
            // across ranks).
            if (!ownedStorageActive_)
            {
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
            }

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

            if (doProfile)
            {
                profileOutputReduceWallTime_ +=
                    wallSeconds(outputReduceStart, wallClockNow());
            }
        }

        const auto outputComputeStart =
            doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();

        if
        (
            cloud_.replicatedMeshActive()
         && computeOutputFields
         && !densityOnly_
        )
        {
            const scalar kBLocal = physicoChemical::k.value();

            for (label outputI = 0; outputI < outputLoopSize; ++outputI)
            {
                const label celli =
                    outputCellsPtr ? (*outputCellsPtr)[outputI] : outputI;
                if (dsmcNCum_[toOwned(celli)] > 1e-3)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[celli];

                    dsmcNMean_()[celli] = dsmcNCum_[toOwned(celli)]/nAvTimeSteps;

                    const scalar rhoNMean =
                        nCum_[toOwned(celli)]/(nAvTimeSteps*cellVolume);
                    const scalar rhoMMean =
                        mCum_[toOwned(celli)]/(nAvTimeSteps*cellVolume);

                    rhoN_()[celli] = rhoNMean;
                    rhoM_()[celli] = rhoMMean;
                    UMean_()[celli] = momentumCum_[toOwned(celli)]/mCum_[toOwned(celli)];

                    const scalar linearKEMean =
                        0.5*linearKECum_[toOwned(celli)]/(cellVolume*nAvTimeSteps);

                    Ttra_()[celli] =
                        2.0/(3.0*kBLocal*rhoNMean)
                       *(
                            linearKEMean
                          - 0.5*rhoMMean*(UMean_()[celli] & UMean_()[celli])
                        );

                    p_()[celli] = rhoNMean*kBLocal*Ttra_()[celli];
                }
                else
                {
                    dsmcNMean_()[celli] = 0.001;
                    rhoN_()[celli] = 0.0;
                    rhoM_()[celli] = 0.0;
                    UMean_()[celli] = vector::zero;
                    Ttra_()[celli] = 0.0;
                    p_()[celli] = 0.0;
                }
            }
        }

        if (computeOutputFields)
        {
            if (densityOnly_)
            {
                for (label outputI = 0; outputI < outputLoopSize; ++outputI)
                {
                    const label celli =
                        outputCellsPtr ? (*outputCellsPtr)[outputI] : outputI;
                    if (dsmcNCum_[toOwned(celli)] > SMALL)
                    {
                        const scalar cellVolume = mesh_.cellVolumes()[celli];

                        dsmcNMean_()[celli] = dsmcNCum_[toOwned(celli)]/nAvTimeSteps;

                        rhoN_()[celli] = nCum_[toOwned(celli)]/(nAvTimeSteps*cellVolume);
                        rhoM_()[celli] = mCum_[toOwned(celli)]/(nAvTimeSteps*cellVolume);
                    }
                    else
                    {
                        // not zero so that weighted decomposition still works
                        dsmcNMean_()[celli] = 0.001;
                        rhoN_()[celli] = 0.0;
                        rhoM_()[celli] = 0.0;
                    }

                    if (dsmcN_()[celli] < SMALL)
                    {
                        // not zero so that weighted decomposition still works
                        dsmcN_()[celli] = 0.001;
                    }
                }
            }
            else
            {
                const label nSpecies = speciesIds_.size();

            for (label outputI = 0; outputI < outputLoopSize; ++outputI)
            {
                const label celli =
                    outputCellsPtr ? (*outputCellsPtr)[outputI] : outputI;
                //- Fields initialisation 
                scalar moleculesRhoN = 0.0;
                Tvib_()[celli] = 0.0;
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
                    dsmcNCum_[toOwned(celli)] > SMALL
                  ? dsmcZetaRotCum_[toOwned(celli)]/dsmcNCum_[toOwned(celli)]
                  : 0.0
                );

                Trot_()[celli] =
                (
                    dsmcZetaRotCum_[toOwned(celli)] > SMALL
                  ? 2.0*dsmcErotCum_[toOwned(celli)]/(kB*dsmcZetaRotCum_[toOwned(celli)])
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
                    const bool enoughTvibSamples =
                        dsmcNSpeciesCum_[i][toOwned(celli)] >= nMinParcelsTvib_;

                    if (nVibMod > 0 && enoughTvibSamples)
                    {
                        moleculesRhoN += nSpeciesCum_[i][toOwned(celli)];
                    }

                    forAll(dsmcSpeciesEvibModCum_[i], mod)
                    {
                        if
                        (
                            dsmcSpeciesEvibModCum_[i][mod][toOwned(celli)] > VSMALL
                         && enoughTvibSamples
                         && nVibMod > 0
                        )
                        {
                            const scalar thetaV =
                                cloud_.constProps(spId).thetaV_m(mod);

                            const scalar iMean =
                                dsmcSpeciesEvibModCum_[i][mod][toOwned(celli)]
                               /(kB*thetaV*dsmcNSpeciesCum_[i][toOwned(celli)]);
                               
                            if (iMean > iMeanMinTvib_)
                            {
                                const scalar logFactor = log(1.0 + 1.0/iMean);

                                speciesTvibMod[i][mod] = thetaV/logFactor;

                                speciesZetaVibMod[i][mod] = 2.0*iMean*logFactor;

                                speciesZetaVib[i] += speciesZetaVibMod[i][mod];
                                    
                                zetaByTvibMod += speciesZetaVibMod[i][mod]
                                    *speciesTvibMod[i][mod];
                            }
                        }
                    }

                    if (speciesZetaVib[i] > SMALL)
                    {
                        speciesTvib[i] = zetaByTvibMod/speciesZetaVib[i];
                        
                        Tvib_()[celli] += nSpeciesCum_[i][toOwned(celli)]*speciesTvib[i];
                            
                        zetaVib_[toOwned(celli)] += nSpeciesCum_[i][toOwned(celli)]
                            *speciesZetaVib[i];    
                    }
                    
                } //- end species loop

                if (moleculesRhoN > SMALL)
                {
                    Tvib_()[celli] /= moleculesRhoN;
                    zetaVib_[toOwned(celli)] /= moleculesRhoN;
                }

                //- Electronic energy mode // TODO Vincent
                //  To reintroduce - I do not trust this part
                scalar zetaElecTot = 0.0;
                Telec_()[celli] = 0.0;

                //- Overall temperature
                Tov_()[celli] =
                    (
                        3.0*Ttra_()[celli]
                      + zetaRotTot*Trot_()[celli]
                      + zetaVib_[toOwned(celli)]*Tvib_()[celli]
                      + zetaElecTot*Telec_()[celli]
                    ) /
                    (3.0 + zetaRotTot + zetaVib_[toOwned(celli)] + zetaElecTot);


                if (measureHeatFluxShearStress_)
                {
                    if (dsmcNCum_[toOwned(celli)] > SMALL)
                    {
                        pressureTensor_()[celli].xx() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                dsmcMuuCum_[toOwned(celli)]
                              - dsmcMCum_[toOwned(celli)]*sqr(UMean_()[celli].x())
                            );
                        pressureTensor_()[celli].xy() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                dsmcMuvCum_[toOwned(celli)]
                              - dsmcMCum_[toOwned(celli)]*UMean_()[celli].x()
                              * UMean_()[celli].y()
                            );
                        pressureTensor_()[celli].xz() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                dsmcMuwCum_[toOwned(celli)]
                              - dsmcMCum_[toOwned(celli)]*UMean_()[celli].x()
                              * UMean_()[celli].z()
                            );

                        pressureTensor_()[celli].yx() =
                            pressureTensor_()[celli].xy();
                        pressureTensor_()[celli].yy() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                dsmcMvvCum_[toOwned(celli)]
                              - dsmcMCum_[toOwned(celli)]*sqr(UMean_()[celli].y())
                            );
                        pressureTensor_()[celli].yz() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                dsmcMvwCum_[toOwned(celli)]
                              - dsmcMCum_[toOwned(celli)]*UMean_()[celli].y()
                              * UMean_()[celli].z()
                            );

                        pressureTensor_()[celli].zx() =
                            pressureTensor_()[celli].xz();
                        pressureTensor_()[celli].zy() =
                            pressureTensor_()[celli].yz();
                        pressureTensor_()[celli].zz() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                dsmcMwwCum_[toOwned(celli)]
                              - dsmcMCum_[toOwned(celli)]*sqr(UMean_()[celli].z())
                            );

                        const scalar scalarPressure =
                            1.0/3.0
                           *(
                                pressureTensor_()[celli].xx()
                              + pressureTensor_()[celli].yy()
                              + pressureTensor_()[celli].zz()
                            );

                        shearStressTensor_()[celli] = -pressureTensor_()[celli];
                        shearStressTensor_()[celli].xx() += scalarPressure;
                        shearStressTensor_()[celli].yy() += scalarPressure;
                        shearStressTensor_()[celli].zz() += scalarPressure;

                        //- terms involving pressure tensor should not be
                        //  multiplied by the number density
                        //  (see Bird corrigendum)

                        heatFluxVector_()[celli].x() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                0.5*dsmcMccuCum_[toOwned(celli)]
                              - 0.5*dsmcMccCum_[toOwned(celli)]*UMean_()[celli].x()
                              + dsmcEuCum_[toOwned(celli)]
                              - dsmcECum_[toOwned(celli)]*UMean_()[celli].x()
                            )
                          - pressureTensor_()[celli].xx()*UMean_()[celli].x()
                          - pressureTensor_()[celli].xy()*UMean_()[celli].y()
                          - pressureTensor_()[celli].xz()*UMean_()[celli].z();

                        heatFluxVector_()[celli].y() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                0.5*dsmcMccvCum_[toOwned(celli)]
                              - 0.5*dsmcMccCum_[toOwned(celli)]*UMean_()[celli].y()
                              + dsmcEvCum_[toOwned(celli)]
                              - dsmcECum_[toOwned(celli)]*UMean_()[celli].y()
                            )
                          - pressureTensor_()[celli].yx()*UMean_()[celli].x()
                          - pressureTensor_()[celli].yy()*UMean_()[celli].y()
                          - pressureTensor_()[celli].yz()*UMean_()[celli].z();

                        heatFluxVector_()[celli].z() =
                            rhoN_()[celli]/dsmcNCum_[toOwned(celli)]
                           *(
                                0.5*dsmcMccwCum_[toOwned(celli)]
                              - 0.5*dsmcMccCum_[toOwned(celli)]*UMean_()[celli].z()
                              + dsmcEwCum_[toOwned(celli)]
                              - dsmcECum_[toOwned(celli)]*UMean_()[celli].z()
                            )
                          - pressureTensor_()[celli].zx()*UMean_()[celli].x()
                          - pressureTensor_()[celli].zy()*UMean_()[celli].y()
                          - pressureTensor_()[celli].zz()*UMean_()[celli].z();
                    }
                    else
                    {
                        pressureTensor_()[celli] = tensor::zero;
                        shearStressTensor_()[celli] = tensor::zero;
                        heatFluxVector_()[celli] = vector::zero;
                    }
                }
                
                if (dsmcNCum_[toOwned(celli)] > SMALL and Ttra_()[celli] > SMALL)
                {
                    forAll(speciesIds_, i)
                    {
                        const label spId = speciesIds_[i];
                        const scalar speciesZetaRot =
                            cloud_.constProps(spId)
                              .rotationalDegreesOfFreedom();
                        
                        const scalar Xs = nSpeciesCum_[i][toOwned(celli)]
                            /nCum_[toOwned(celli)];

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
                            gamma*kB/molecularMass*Ttra_()[celli]
                        );

                    Ma_()[celli] = mag(UMean_()[celli])/speedOfSound;
                }
                else
                {
                    Ma_()[celli] = 0.0;
                }

                if (measureMeanFreePath_ && Ttra_()[celli] > 1.0)
                {
                    const scalar deltaT = cloud_.deltaTValue(celli);
                    
                    mfp_()[celli] = 0.0;
                    meanCollisionRate_()[celli] = 0.0;
                    
                    forAll(speciesIds_, s)
                    {
                        const label spIdp = speciesIds_[s];
                        
                        speciesMfp_[s][toOwned(celli)] = 0.0;
                        speciesMcr_[s][toOwned(celli)] = 0.0;

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
                                dsmcNSpeciesCum_[r][toOwned(celli)] > SMALL
                             && Ttra_()[celli] > SMALL
                            )
                            {
                                const scalar nDensQ =
                                    nSpeciesCum_[r][toOwned(celli)]
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
                                speciesMfp_[s][toOwned(celli)] += pi*sqr(dPQ)*nDensQ
                                   *pow
                                    (
                                        mfpTref_/Ttra_()[celli], omegaPQ - 0.5
                                    )*sqrt(1.0+massRatio);

                                // Bird 1994, eq (4.74)
                                speciesMcr_[s][toOwned(celli)] +=
                                    2.0*sqrt(pi)*sqr(dPQ)*nDensQ
                                   *pow
                                    (
                                        Ttra_()[celli]/mfpTref_, 1.0 - omegaPQ
                                    )
                                   *sqrt
                                    (
                                        2.0*kB*mfpTref_/reducedMass
                                    );
                            }
                        }

                        if (speciesMfp_[s][toOwned(celli)] > SMALL)
                        {
                            speciesMfp_[s][toOwned(celli)] = 1.0/speciesMfp_[s][toOwned(celli)];
                        }
                    }

                    meanCollisionSeparation_()[celli] =
                    (
                        dsmcNCollsCum_[toOwned(celli)] > SMALL
                      ? collisionSeparation_[toOwned(celli)]/dsmcNCollsCum_[toOwned(celli)]
                      : GREAT
                    );

                    if (nCum_[toOwned(celli)] > SMALL)
                    {
                        // const scalar symmFactor = 2.0;
                        // TODO (s == r ? 1.0 : 2.0);
                        measuredCollisionRate_()[celli] = dsmcNCollsCum_[toOwned(celli)]
                            *cloud_.nParticles(celli)/(nCum_[toOwned(celli)]*deltaT);
                    }

                    if (rhoN_()[celli] > SMALL)
                    {
                        forAll(speciesIds_, i)
                        {
                            const scalar rhoNi = nSpeciesCum_[i][toOwned(celli)];

                            // Bird 1994, eq (4.77)
                            mfp_()[celli] += speciesMfp_[i][toOwned(celli)]
                                *rhoNi/nCum_[toOwned(celli)];

                            // Bird 1994, eq (1.38)
                            meanCollisionRate_()[celli] +=
                                speciesMcr_[i][toOwned(celli)]*rhoNi/nCum_[toOwned(celli)];
                        }
                    }

                    if (mfp_()[celli] < SMALL)
                    {
                        mfp_()[celli] = GREAT;
                    }

                    if (meanCollisionRate_()[celli] > SMALL)
                    {
                        meanCollisionTime_()[celli] =
                            1.0/meanCollisionRate_()[celli];
                        mctToDt_()[celli] = meanCollisionTime_()[celli]/deltaT;
                    }
                    else
                    {
                        meanCollisionTime_()[celli] = GREAT;
                        mctToDt_()[celli] = GREAT;
                    }

                    if (mfp_()[celli] != GREAT)
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

                        mfpToDx_()[celli] = mfp_()[celli]/maxCellDx;

                        SOF_()[celli] =
                        (
                            mfp_()[celli] > SMALL
                          ? meanCollisionSeparation_()[celli]/mfp_()[celli]
                          : 0.0
                        );
                    }
                    else
                    {
                        mfpToDx_()[celli] = GREAT;
                        SOF_()[celli] = GREAT;
                    }

                    // when few particles in cell, undesired refinement
                    // this condition should eliminates this problem
                    if (dsmcN_()[celli] >= 4.0)
                    {
                        DxToMfp_()[celli] = 1.0/mfpToDx_()[celli];
                    }
                }

                if (measureClassifications_)
                {
                    if (dsmcNCum_[toOwned(celli)] > SMALL)
                    {
                        classIDistribution_()[celli] = dsmcNClassICum_[toOwned(celli)]
                            /dsmcNCum_[toOwned(celli)];
                        classIIDistribution_()[celli] = dsmcNClassIICum_[toOwned(celli)]
                            /dsmcNCum_[toOwned(celli)];
                        classIIIDistribution_()[celli] = dsmcNClassIIICum_[toOwned(celli)]
                            /dsmcNCum_[toOwned(celli)];
                    }
                }

                if (measureErrors_)
                {
                    if
                    (
                         dsmcNMean_()[celli] > SMALL && Ma_()[celli] > SMALL
                      && gamma > SMALL && particleCv > SMALL
                    )
                    {
                        const scalar deno = sqrt(dsmcNMean_()[celli]*nAvTimeSteps);
                        
                        densityError_()[celli] = 1.0/deno;
                        velocityError_()[celli] = 1.0/(deno*Ma_()[celli]*sqrt(gamma));
                        temperatureError_()[celli] = sqrt(kB/particleCv)/deno;
                        pressureError_()[celli] = sqrt(gamma)/deno;
                    }

                }
            } //- end loop over cells

            //- Computing boundary measurements: loop over all boundary patches
            forAll(rhoNBF_, j)
            {
                //- Determine of the type of patch: patch, wall, cyclic, ...
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                const bool useOwnedOutputBoundaryFaces =
                    processorWrite
                 && cloud_.replicatedMeshActive()
                 && ownedBoundaryFaces_.size() == mesh_.boundaryMesh().size();
                
                const bool isWall = isA<wallPolyPatch>(patch);
                
                const bool isNonEmptyNonCyclic = isA<polyPatch>(patch)
                    && !isA<emptyPolyPatch>(patch)
                    && !isA<cyclicPolyPatch>(patch);

                if (isWall)
                {
                    //- Loop over all wall boundary faces
                    const label nFaces =
                        useOwnedOutputBoundaryFaces
                      ? ownedBoundaryFaces_[j].size()
                      : patch.size();

                    for (label facei = 0; facei < nFaces; ++facei)
                    {
                        const label k =
                            useOwnedOutputBoundaryFaces
                          ? ownedBoundaryFaces_[j][facei]
                          : facei;
                        const label celli = boundaryCells_[j][k];
                        
                        //- Initialise face fields
                        Tvib_().boundaryFieldRef()[j][k] = 0.0;
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

                        rhoN_().boundaryFieldRef()[j][k] = rhoNMean;
                        rhoM_().boundaryFieldRef()[j][k] = rhoMMean;
                        
                        //- Instantaneous and sampled numbers of DSMC parcels
                        //  are that of the neighbouring cell
                        dsmcN_().boundaryFieldRef()[j][k] = dsmcN_()[celli];
                        dsmcNMean_().boundaryFieldRef()[j][k] = dsmcNMean_()[celli];

                        //- Translational energy mode and velocity
                        if (rhoMMean > VSMALL)
                        {
                            UMean_().boundaryFieldRef()[j][k] = momentumBF_[j][k]
                                /rhoMBF_[j][k];

                            Ttra_().boundaryFieldRef()[j][k] =
                                2.0/(3.0*kB*rhoNMean)*
                                (
                                    linearKEMean - 0.5*rhoMMean*
                                    (
                                        UMean_().boundaryField()[j][k]
                                      & UMean_().boundaryField()[j][k]
                                    )
                                );
                        }
                        else
                        {
                            UMean_().boundaryFieldRef()[j][k] = vector::zero;
                            Ttra_().boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Rotational energy mode
                        const scalar zetaRotTot =
                        (
                            rhoNBF_[j][k] > SMALL
                          ? zetaRotBF_[j][k]/rhoNBF_[j][k]
                          : 0.0
                        );

                        Trot_().boundaryFieldRef()[j][k] =
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

                            const label tvibCellI = toOwned(celli);
                            const bool enoughTvibSamples =
                                tvibCellI >= 0
                             && tvibCellI < dsmcNSpeciesCum_[i].size()
                             && dsmcNSpeciesCum_[i][tvibCellI] >= nMinParcelsTvib_;

                            if
                            (
                                nVibMod > 0
                             && enoughTvibSamples
                             && speciesRhoNBF_[i][j][k] > SMALL
                            )
                            {
                                moleculesRhoN += speciesRhoNBF_[i][j][k];
                            }

                            if
                            (
                                speciesRhoNBF_[i][j][k] > SMALL
                             && enoughTvibSamples
                            )
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
                                speciesTvibBF_[i][j][k] = zetaByTvibMod
                                    /speciesZetaVibBF_[i][j][k];
                                    
                                Tvib_().boundaryFieldRef()[j][k] +=
                                    speciesRhoNBF_[i][j][k]
                                   *speciesTvibBF_[i][j][k];

                                zetaVibBF_[j][k] +=
                                    speciesRhoNBF_[i][j][k]
                                   *speciesZetaVibBF_[i][j][k];
                            }
                        }

                        if (moleculesRhoN > SMALL)
                        {
                            Tvib_().boundaryFieldRef()[j][k] /= moleculesRhoN;
                            zetaVibBF_[j][k] /= moleculesRhoN;
                        }

                        //- Electronic energy mode // TODO Vincent
                        //  Removed temporarily - I don't trust this part
                        scalar zetaElecTot = 0.0;
                        Telec_().boundaryFieldRef()[j][k] = 0.0;

                        Tov_().boundaryFieldRef()[j][k] =
                            (
                                (3.0*Ttra_().boundaryField()[j][k])
                              + (zetaRotTot*Trot_().boundaryField()[j][k])
                              + (
                                    zetaVibBF_[j][k]
                                   *Tvib_().boundaryField()[j][k]
                                )
                              + (
                                    zetaElecTot
                                   *Telec_().boundaryFieldRef()[j][k]
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
                                 && Ttra_().boundaryFieldRef()[j][k] > SMALL
                                  ? kB/molecularMassBF
                                  : 0.0
                                );

                            const scalar gamma = molarCpBF_trarot/molarCvBF_trarot;

                            const scalar speedOfSound =
                                sqrt
                                (
                                    gamma*gasConstant*Ttra_().boundaryField()[j][k]
                                );

                            Ma_().boundaryFieldRef()[j][k] =
                                mag(UMean_().boundaryField()[j][k])/speedOfSound;
                        }
                        else
                        {
                            Ma_().boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Force density
                        fD_().boundaryFieldRef()[j][k] = fDBF_[j][k]/nAvTimeSteps;

                        //- Surface pressure
                        p_().boundaryFieldRef()[j][k] =
                            fD_().boundaryField()[j][k] & n_[j][k];
                            
                        //- Wall shear stress
                        tau_().boundaryFieldRef()[j][k] =
                            sqrt
                            (
                                sqr(fD_().boundaryField()[j][k] & t1_[j][k])
                              + sqr(fD_().boundaryField()[j][k] & t2_[j][k])
                            );
                            
                        //- Heat flux
                        q_().boundaryFieldRef()[j][k] = qBF_[j][k]/nAvTimeSteps;

                        if (writeWallPressureCoefficient_)
                        {
                            Cp_->boundaryFieldRef()[j][k] =
                                (
                                    p_().boundaryField()[j][k]
                                  - forceMomentPInf_
                                )/forceMomentQInf_;
                        }

                        if (writeForceMomentCoefficients_)
                        {
                            Ch_->boundaryFieldRef()[j][k] =
                                q_().boundaryField()[j][k]
                               /forceMomentHeatFluxInf_;
                        }
                        
                        //- ZeroGradient condition assumed for Optional fields
                        if (measureMeanFreePath_)
                        {
                            mfp_().boundaryFieldRef()[j][k] = mfp_()[celli];
                            SOF_().boundaryFieldRef()[j][k] = SOF_()[celli];
                            mfpToDx_().boundaryFieldRef()[j][k] = mfpToDx_()[celli];
                            meanCollisionRate_().boundaryFieldRef()[j][k] =
                                meanCollisionRate_()[celli];
                            meanCollisionTime_().boundaryFieldRef()[j][k] =
                                meanCollisionTime_()[celli];
                            mctToDt_().boundaryFieldRef()[j][k] = mctToDt_()[celli];
                        }
                    }
                }
                else if (isNonEmptyNonCyclic)
                {
                    //- Loop over all boundary faces and set zeroGradient
                    //  conditions
                    const label nFaces =
                        useOwnedOutputBoundaryFaces
                      ? ownedBoundaryFaces_[j].size()
                      : boundaryCells_[j].size();

                    for (label facei = 0; facei < nFaces; ++facei)
                    {
                        const label k =
                            useOwnedOutputBoundaryFaces
                          ? ownedBoundaryFaces_[j][facei]
                          : facei;
                        const label celli = boundaryCells_[j][k];

                        //- Instantaneous and sampled numbers of DSMC parcels
                        //  are that of the neighbouring cell
                        dsmcN_().boundaryFieldRef()[j][k] = dsmcN_()[celli];
                        dsmcNMean_().boundaryFieldRef()[j][k] =
                            dsmcNMean_()[celli];
                            
                        //- Number density and mass density fields
                        rhoN_().boundaryFieldRef()[j][k] = rhoN_()[celli];
                        rhoM_().boundaryFieldRef()[j][k] = rhoM_()[celli];
                        
                        //- Temperature fields
                        Ttra_().boundaryFieldRef()[j][k] = Ttra_()[celli];
                        Trot_().boundaryFieldRef()[j][k] = Trot_()[celli];
                        Tvib_().boundaryFieldRef()[j][k] = Tvib_()[celli];
                        Tov_().boundaryFieldRef()[j][k] = Tov_()[celli];
                        
                        //- Pressure, Mach and velocity fields
                        p_().boundaryFieldRef()[j][k] = p_()[celli];
                        Ma_().boundaryFieldRef()[j][k] = Ma_()[celli];
                        UMean_().boundaryFieldRef()[j][k] = UMean_()[celli];
                        
                        //- Optional fields
                        if (measureMeanFreePath_)
                        {
                            mfp_().boundaryFieldRef()[j][k] = mfp_()[celli];
                            SOF_().boundaryFieldRef()[j][k] = SOF_()[celli];
                            mfpToDx_().boundaryFieldRef()[j][k] = mfpToDx_()[celli];
                            meanCollisionRate_().boundaryFieldRef()[j][k] =
                                meanCollisionRate_()[celli];
                            meanCollisionTime_().boundaryFieldRef()[j][k] =
                                meanCollisionTime_()[celli];
                            mctToDt_().boundaryFieldRef()[j][k] = mctToDt_()[celli];
                        }
                        
                        if (measureHeatFluxShearStress_)
                        {
                            shearStressTensor_().boundaryFieldRef()[j][k] =
                                shearStressTensor_()[celli];
                            heatFluxVector_().boundaryFieldRef()[j][k] =
                                heatFluxVector_()[celli];
                            pressureTensor_().boundaryFieldRef()[j][k] =
                                pressureTensor_()[celli];
                        }
                        
                        if (measureClassifications_)
                        {
                            classIDistribution_().boundaryFieldRef()[j][k] =
                                classIDistribution_()[celli];
                            classIIDistribution_().boundaryFieldRef()[j][k] =
                                classIIDistribution_()[celli];
                            classIIIDistribution_().boundaryFieldRef()[j][k] =
                                classIIIDistribution_()[celli];
                        }
                    }
                }
            }

            if (doProfile)
            {
                profileOutputComputeWallTime_ +=
                    wallSeconds(outputComputeStart, wallClockNow());
            }

            const auto fieldWriteStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();

            if (writeForceMoment_)
            {
                vector totalForce(vector::zero);
                vector totalMoment(vector::zero);
                vectorField wallForces
                (
                    forceMomentWallPatchIds_.size(),
                    vector::zero
                );
                vectorField wallMoments
                (
                    forceMomentWallPatchIds_.size(),
                    vector::zero
                );
                const bool useOwnedForceMomentFaces =
                    processorWrite
                 && cloud_.replicatedMeshActive()
                 && ownedBoundaryFaces_.size() == mesh_.boundaryMesh().size();

                forAll(mesh_.boundaryMesh(), patchi)
                {
                    const polyPatch& patch = mesh_.boundaryMesh()[patchi];
                    const label wallPatchI = forceMomentWallPatchIndices_[patchi];

                    if (wallPatchI < 0)
                    {
                        continue;
                    }

                    const label nFaces =
                        useOwnedForceMomentFaces
                      ? ownedBoundaryFaces_[patchi].size()
                      : patch.size();

                    for (label localFaceI = 0; localFaceI < nFaces; ++localFaceI)
                    {
                        const label facei =
                            useOwnedForceMomentFaces
                          ? ownedBoundaryFaces_[patchi][localFaceI]
                          : localFaceI;
                        const label meshFacei = patch.start() + facei;
                        const vector faceForce =
                            fD_().boundaryField()[patchi][facei]
                           *mag(mesh_.faceAreas()[meshFacei]);

                        const vector faceMoment =
                            (mesh_.faceCentres()[meshFacei]
                           - forceMomentReferencePoint_) ^ faceForce;

                        totalForce += faceForce;
                        totalMoment += faceMoment;
                        wallForces[wallPatchI] += faceForce;
                        wallMoments[wallPatchI] += faceMoment;
                    }
                }

                if (processorWrite || !cloud_.replicatedMeshActive())
                {
                    vectorField forceMoment
                    (
                        2 + 2*forceMomentWallPatchIds_.size(),
                        vector::zero
                    );
                    forceMoment[0] = totalForce;
                    forceMoment[1] = totalMoment;
                    forAll(forceMomentWallPatchIds_, wallPatchI)
                    {
                        forceMoment[2 + 2*wallPatchI] = wallForces[wallPatchI];
                        forceMoment[3 + 2*wallPatchI] = wallMoments[wallPatchI];
                    }
                    sumReduceField(forceMoment);
                    totalForce = forceMoment[0];
                    totalMoment = forceMoment[1];
                    forAll(forceMomentWallPatchIds_, wallPatchI)
                    {
                        wallForces[wallPatchI] = forceMoment[2 + 2*wallPatchI];
                        wallMoments[wallPatchI] = forceMoment[3 + 2*wallPatchI];
                    }
                }

                if (cloud_.isOutputRank())
                {
                    const fileName timePath
                    (
                        time_.time().path()/time_.time().timeName()
                    );
                    const fileName uniformPath(timePath/"uniform");
                    mkDir(timePath);
                    mkDir(uniformPath);

                    OFstream forceMomentFile
                    (
                        uniformPath/("forceMoment_" + fieldName_)
                    );

                    if (!forceMomentFile.good())
                    {
                        FatalErrorInFunction
                            << "Cannot write force and moment output for field "
                            << fieldName_ << " to " << uniformPath
                            << exit(FatalError);
                    }

                    forceMomentFile
                        << "FoamFile" << nl
                        << "{" << nl
                        << "    version     2.0;" << nl
                        << "    format      ascii;" << nl
                        << "    class       dictionary;" << nl
                        << "    location    \"" << time_.time().timeName()
                        << "/uniform\";" << nl
                        << "    object      forceMoment_" << fieldName_ << ";" << nl
                        << "}" << nl << nl
                        << "referencePoint " << forceMomentReferencePoint_ << ";" << nl
                        << "referencePointAutomatic "
                        << forceMomentReferencePointAutomatic_ << ";" << nl
                        << "patchSelection allWallPatches;" << nl
                        << "force " << totalForce << "; // N" << nl
                        << "moment " << totalMoment << "; // N m" << nl;

                    vector totalForceCoefficients(vector::zero);
                    vector totalMomentCoefficients(vector::zero);
                    if (writeForceMomentCoefficients_)
                    {
                        const scalar forceScale =
                            1.0/(forceMomentQInf_*forceMomentReferenceArea_);
                        const scalar momentScale =
                            forceScale/forceMomentReferenceLength_;

                        totalForceCoefficients = vector
                        (
                            (totalForce & forceMomentDragDirection_)*forceScale,
                            (totalForce & forceMomentLiftDirection_)*forceScale,
                            (totalForce & forceMomentSideDirection_)*forceScale
                        );
                        totalMomentCoefficients = vector
                        (
                            (totalMoment & forceMomentDragDirection_)*momentScale,
                            (totalMoment & forceMomentLiftDirection_)*momentScale,
                            (totalMoment & forceMomentSideDirection_)*momentScale
                        );

                        forceMomentFile
                            << nl << "coefficients" << nl << "{" << nl
                            << "    qInf " << forceMomentQInf_ << "; // Pa" << nl
                            << "    pInf " << forceMomentPInf_ << "; // Pa" << nl
                            << "    heatFluxInf " << forceMomentHeatFluxInf_
                            << "; // W/m2" << nl
                            << "    referenceArea " << forceMomentReferenceArea_
                            << "; // m2" << nl
                            << "    referenceLength "
                            << forceMomentReferenceLength_ << "; // m" << nl
                            << "    freestreamSource "
                            << forceMomentFreestreamSource_ << ";" << nl
                            << "    referenceAreaDefinition "
                            << forceMomentReferenceAreaDefinition_ << ";" << nl
                            << "    referenceLengthDefinition "
                            << forceMomentReferenceLengthDefinition_ << ";" << nl
                            << "    dragDirection " << forceMomentDragDirection_
                            << ";" << nl
                            << "    liftDirection " << forceMomentLiftDirection_
                            << ";" << nl
                            << "    sideDirection " << forceMomentSideDirection_
                            << ";" << nl
                            << "    CD " << totalForceCoefficients.x() << ";" << nl
                            << "    CL " << totalForceCoefficients.y() << ";" << nl
                            << "    CS " << totalForceCoefficients.z() << ";" << nl
                            << "    CMdrag " << totalMomentCoefficients.x()
                            << ";" << nl
                            << "    CMlift " << totalMomentCoefficients.y()
                            << ";" << nl
                            << "    CMside " << totalMomentCoefficients.z()
                            << ";" << nl
                            << "    CMglobal " << totalMoment*momentScale
                            << ";" << nl
                            << "}" << nl;
                    }

                    if (!forceMomentWallPatchIds_.empty())
                    {
                        forceMomentFile << nl << "wallPatches" << nl << "{" << nl;

                        forAll(forceMomentWallPatchIds_, wallPatchI)
                        {
                            const label patchi = forceMomentWallPatchIds_[wallPatchI];
                            forceMomentFile
                                << "    " << mesh_.boundaryMesh()[patchi].name()
                                << nl << "    {" << nl
                                << "        force " << wallForces[wallPatchI]
                                << "; // N" << nl
                                << "        moment " << wallMoments[wallPatchI]
                                << "; // N m" << nl;

                            if (writeForceMomentCoefficients_)
                            {
                                const scalar forceScale =
                                    1.0/
                                    (
                                        forceMomentQInf_
                                       *forceMomentReferenceArea_
                                    );
                                const scalar momentScale =
                                    forceScale/forceMomentReferenceLength_;
                                const vector wallForceCoefficients
                                (
                                    (wallForces[wallPatchI]
                                   & forceMomentDragDirection_)*forceScale,
                                    (wallForces[wallPatchI]
                                   & forceMomentLiftDirection_)*forceScale,
                                    (wallForces[wallPatchI]
                                   & forceMomentSideDirection_)*forceScale
                                );
                                const vector wallMomentCoefficients
                                (
                                    (wallMoments[wallPatchI]
                                   & forceMomentDragDirection_)*momentScale,
                                    (wallMoments[wallPatchI]
                                   & forceMomentLiftDirection_)*momentScale,
                                    (wallMoments[wallPatchI]
                                   & forceMomentSideDirection_)*momentScale
                                );

                                forceMomentFile
                                    << "        coefficients" << nl
                                    << "        {" << nl
                                    << "            CD "
                                    << wallForceCoefficients.x() << ";" << nl
                                    << "            CL "
                                    << wallForceCoefficients.y() << ";" << nl
                                    << "            CS "
                                    << wallForceCoefficients.z() << ";" << nl
                                    << "            CMdrag "
                                    << wallMomentCoefficients.x() << ";" << nl
                                    << "            CMlift "
                                    << wallMomentCoefficients.y() << ";" << nl
                                    << "            CMside "
                                    << wallMomentCoefficients.z() << ";" << nl
                                    << "            CMglobal "
                                    << wallMoments[wallPatchI]*momentScale
                                    << ";" << nl
                                    << "        }" << nl;
                            }

                            forceMomentFile << "    }" << nl;
                        }

                        forceMomentFile << "}" << nl;
                    }

                    Info<< "dsmcVolFields [" << fieldName_
                        << "]: force " << totalForce << " N, moment about "
                        << forceMomentReferencePoint_ << " " << totalMoment
                        << " N m";

                    if (writeForceMomentCoefficients_)
                    {
                        Info<< ", CD " << totalForceCoefficients.x()
                            << ", CL " << totalForceCoefficients.y()
                            << ", CS " << totalForceCoefficients.z()
                            << ", CMside " << totalMomentCoefficients.z();
                    }

                    Info<< endl;
                }
            }

            const bool writeDsmcN =
                outputFieldEnabled("dsmcN", dsmcN_().name());
            const bool writeDsmcNMean =
                outputFieldEnabled("dsmcNMean", dsmcNMean_().name());
            const bool writeRhoN =
                outputFieldEnabled("rhoN", rhoN_().name());
            const bool writeRhoM =
                outputFieldEnabled("rhoM", rhoM_().name());
            const bool writeP =
                outputFieldEnabled("p", p_().name());
            const bool writeTtra =
                outputFieldEnabled("Ttra", Ttra_().name());
            const bool writeU =
                outputFieldEnabled("U", UMean_().name());
            const bool writeMa =
                outputFieldEnabled("Ma", Ma_().name());
            const bool writeQ =
                outputFieldEnabled("wallHeatFlux", q_().name());
            const bool writeFD =
                outputFieldEnabled("fD", fD_().name());
            const bool writeTau =
                outputFieldEnabled("wallShearStress", tau_().name());
            const bool writeCp =
                writeWallPressureCoefficient_
             && outputFieldEnabled
                (
                    "wallPressureCoefficient",
                    Cp_->name(),
                    writeWallPressureCoefficient_
                );
            const bool writeCh =
                writeForceMomentCoefficients_
             && outputFieldEnabled
                (
                    "wallHeatFluxCoefficient",
                    Ch_->name(),
                    writeForceMomentCoefficients_
                );
            const bool writeTrot =
                outputFieldEnabled
                (
                    "Trot",
                    Trot_().name(),
                    writeRotationalTemperature_
                );
            const bool writeTvib =
                outputFieldEnabled
                (
                    "Tvib",
                    Tvib_().name(),
                    writeVibrationalTemperature_
                );
            const bool writeTelec =
                outputFieldEnabled
                (
                    "Telec",
                    Telec_().name(),
                    writeElectronicTemperature_
                );
            const bool writeTov =
                outputFieldEnabled
                (
                    "Tov",
                    Tov_().name(),
                    writeRotationalTemperature_
                 || writeVibrationalTemperature_
                 || writeElectronicTemperature_
                );
            const bool writeMfp =
                measureMeanFreePath_
             && outputFieldEnabled("mfp", mfp_().name(), measureMeanFreePath_);
            const bool writeMfpToDx =
                measureMeanFreePath_
             && outputFieldEnabled
                (
                    "mfpToDx",
                    mfpToDx_().name(),
                    measureMeanFreePath_
                );
            const bool writeMct =
                measureMeanFreePath_
             && outputFieldEnabled
                (
                    "mct",
                    meanCollisionTime_().name(),
                    measureMeanFreePath_
                );
            const bool writeMctToDt =
                measureMeanFreePath_
             && outputFieldEnabled
                (
                    "mctToDt",
                    mctToDt_().name(),
                    measureMeanFreePath_
                );
            const bool writeSOF =
                measureMeanFreePath_
             && outputFieldEnabled("SOFP", SOF_().name(), measureMeanFreePath_);
            const bool writeClassI =
                measureClassifications_
             && outputFieldEnabled
                (
                    "classIDistribution",
                    classIDistribution_().name(),
                    measureClassifications_
                );
            const bool writeClassII =
                measureClassifications_
             && outputFieldEnabled
                (
                    "classIIDistribution",
                    classIIDistribution_().name(),
                    measureClassifications_
                );
            const bool writeClassIII =
                measureClassifications_
             && outputFieldEnabled
                (
                    "classIIIDistribution",
                    classIIIDistribution_().name(),
                    measureClassifications_
                );
            const bool writeDensityError =
                measureErrors_
             && outputFieldEnabled
                (
                    "rhoMError",
                    densityError_().name(),
                    measureErrors_
                );
            const bool writeVelocityError =
                measureErrors_
             && outputFieldEnabled
                (
                    "UError",
                    velocityError_().name(),
                    measureErrors_
                );
            const bool writeTemperatureError =
                measureErrors_
             && outputFieldEnabled
                (
                    "TError",
                    temperatureError_().name(),
                    measureErrors_
                );
            const bool writePressureError =
                measureErrors_
             && outputFieldEnabled
                (
                    "pError",
                    pressureError_().name(),
                    measureErrors_
                );
            const bool writeHeatFluxVector =
                measureHeatFluxShearStress_
             && outputFieldEnabled
                (
                    "heatFluxVector",
                    heatFluxVector_().name(),
                    measureHeatFluxShearStress_
                );
            const bool writePressureTensor =
                measureHeatFluxShearStress_
             && outputFieldEnabled
                (
                    "pressureTensor",
                    pressureTensor_().name(),
                    measureHeatFluxShearStress_
                );
            const bool writeShearStressTensor =
                measureHeatFluxShearStress_
             && outputFieldEnabled
                (
                    "shearStressTensor",
                    shearStressTensor_().name(),
                    measureHeatFluxShearStress_
                );

            setProcessorWriteOpt(dsmcN_(), writeDsmcN);
            setProcessorWriteOpt(dsmcNMean_(), writeDsmcNMean);
            setProcessorWriteOpt(rhoN_(), writeRhoN);
            setProcessorWriteOpt(rhoM_(), writeRhoM);
            setProcessorWriteOpt(p_(), writeP);
            setProcessorWriteOpt(Ttra_(), writeTtra);
            setProcessorWriteOpt(UMean_(), writeU);
            if (Ma_.valid()) setProcessorWriteOpt(Ma_(), writeMa);
            if (q_.valid()) setProcessorWriteOpt(q_(), writeQ);
            if (fD_.valid()) setProcessorWriteOpt(fD_(), writeFD);
            if (tau_.valid()) setProcessorWriteOpt(tau_(), writeTau);
            if (Cp_.valid()) setProcessorWriteOpt(Cp_(), writeCp);
            if (Ch_.valid()) setProcessorWriteOpt(Ch_(), writeCh);
            if (Trot_.valid()) setProcessorWriteOpt(Trot_(), writeTrot);
            if (Tvib_.valid()) setProcessorWriteOpt(Tvib_(), writeTvib);
            if (Telec_.valid()) setProcessorWriteOpt(Telec_(), writeTelec);
            if (Tov_.valid()) setProcessorWriteOpt(Tov_(), writeTov);
            if (mfp_.valid()) setProcessorWriteOpt(mfp_(), writeMfp);
            if (mfpToDx_.valid()) setProcessorWriteOpt(mfpToDx_(), writeMfpToDx);
            if (meanCollisionTime_.valid()) setProcessorWriteOpt(meanCollisionTime_(), writeMct);
            if (mctToDt_.valid()) setProcessorWriteOpt(mctToDt_(), writeMctToDt);
            if (SOF_.valid()) setProcessorWriteOpt(SOF_(), writeSOF);
            if (classIDistribution_.valid()) setProcessorWriteOpt(classIDistribution_(), writeClassI);
            if (classIIDistribution_.valid()) setProcessorWriteOpt(classIIDistribution_(), writeClassII);
            if (classIIIDistribution_.valid()) setProcessorWriteOpt(classIIIDistribution_(), writeClassIII);
            if (densityError_.valid()) setProcessorWriteOpt(densityError_(), writeDensityError);
            if (velocityError_.valid()) setProcessorWriteOpt(velocityError_(), writeVelocityError);
            if (temperatureError_.valid()) setProcessorWriteOpt(temperatureError_(), writeTemperatureError);
            if (pressureError_.valid()) setProcessorWriteOpt(pressureError_(), writePressureError);
            if (heatFluxVector_.valid()) setProcessorWriteOpt(heatFluxVector_(), writeHeatFluxVector);
            if (pressureTensor_.valid()) setProcessorWriteOpt(pressureTensor_(), writePressureTensor);
            if (shearStressTensor_.valid()) setProcessorWriteOpt(shearStressTensor_(), writeShearStressTensor);

            if (processorWrite)
            {
                // Processor output consumes writeOpt as a decompose/write filter.
            }
            else
            {
                //- Write solution fields
                if (writeP) p_().write();
                if (writeTtra) Ttra_().write();
                if (writeU) UMean_().write();
                if (writeMa) Ma_().write();
                if (writeQ) q_().write();
                if (writeFD) fD_().write();
                if (writeTau) tau_().write();
                if (writeCp) Cp_->write();
                if (writeCh) Ch_->write();
                if (writeTrot) Trot_().write();
                if (writeTvib) Tvib_().write();
                if (writeTelec) Telec_().write();
                if (writeTov) Tov_().write();
                if (writeMfp) mfp_().write();
                if (writeMfpToDx) mfpToDx_().write();
                if (writeMct) meanCollisionTime_().write();
                if (writeMctToDt) mctToDt_().write();
                if (writeSOF) SOF_().write();
                if (writeClassI) classIDistribution_().write();
                if (writeClassII) classIIDistribution_().write();
                if (writeClassIII) classIIIDistribution_().write();
                if (writeDensityError) densityError_().write();
                if (writeVelocityError) velocityError_().write();
                if (writeTemperatureError) temperatureError_().write();
                if (writePressureError) pressureError_().write();
                if (writeHeatFluxVector) heatFluxVector_().write();
                if (writePressureTensor) pressureTensor_().write();
                if (writeShearStressTensor) shearStressTensor_().write();
            }

            if (doProfile)
            {
                profileFieldWriteWallTime_ +=
                    wallSeconds(fieldWriteStart, wallClockNow());
            }
        }
        }
        
        //- Reset fields after printing the instantaneous solution ... or
        //  continue sampling
        const bool resetAtOutput = time_.resetFieldsAtOutput();
        const auto outputResetStart =
            doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
        if (resetAtOutput)
        {
            nTimeSteps_ = 0.0;

            // M4: measuredCollisionRate_ is a GeoField (full mesh size);
            // zero it over all cells independently of the cumulative
            // storage layout.  It only exists when measureMeanFreePath_ is
            // enabled (M2: lazily constructed).
            if (measuredCollisionRate_.valid())
            {
                forAll(measuredCollisionRate_(), gCelli)
                {
                    measuredCollisionRate_()[gCelli] = 0.0;
                }
            }

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

        if (doProfile && resetAtOutput)
        {
            profileOutputResetWallTime_ +=
                wallSeconds(outputResetStart, wallClockNow());
        }

        if (averagingAcrossManyRuns_ && !time_.resetFieldsAtOutput())
        {
            writeOut();
        }

        if (restoreLocalAfterOutput)
        {
            const auto outputRestoreStart =
                doProfile ? wallClockNow() : std::chrono::steady_clock::time_point();
            restoreField(dsmcNCum_, localDsmcNCum);
            restoreField(nCum_, localNCum);
            restoreField(dsmcNElecLvlCum_, localDsmcNElecLvlCum);
            restoreField(dsmcMCum_, localDsmcMCum);
            restoreField(mCum_, localMCum);
            restoreField(dsmcLinearKECum_, localDsmcLinearKECum);
            restoreField(linearKECum_, localLinearKECum);
            restoreField(dsmcErotCum_, localDsmcErotCum);
            restoreField(dsmcZetaRotCum_, localDsmcZetaRotCum);
            restoreField(dsmcMuuCum_, localDsmcMuuCum);
            restoreField(dsmcMuvCum_, localDsmcMuvCum);
            restoreField(dsmcMuwCum_, localDsmcMuwCum);
            restoreField(dsmcMvvCum_, localDsmcMvvCum);
            restoreField(dsmcMvwCum_, localDsmcMvwCum);
            restoreField(dsmcMwwCum_, localDsmcMwwCum);
            restoreField(dsmcMccCum_, localDsmcMccCum);
            restoreField(dsmcMccuCum_, localDsmcMccuCum);
            restoreField(dsmcMccvCum_, localDsmcMccvCum);
            restoreField(dsmcMccwCum_, localDsmcMccwCum);
            restoreField(dsmcEuCum_, localDsmcEuCum);
            restoreField(dsmcEvCum_, localDsmcEvCum);
            restoreField(dsmcEwCum_, localDsmcEwCum);
            restoreField(dsmcECum_, localDsmcECum);
            restoreField(zetaVib_, localZetaVib);
            restoreField(dsmcNClassICum_, localDsmcNClassICum);
            restoreField(dsmcNClassIICum_, localDsmcNClassIICum);
            restoreField(dsmcNClassIIICum_, localDsmcNClassIIICum);
            restoreField(collisionSeparation_, localCollisionSeparation);
            restoreField(dsmcNCollsCum_, localDsmcNCollsCum);
            restoreField(dsmcMomentumCum_, localDsmcMomentumCum);
            restoreField(momentumCum_, localMomentumCum);

            restoreFieldList(dsmcSpeciesEelecCum_, localDsmcSpeciesEelecCum);
            restoreFieldList(dsmcNSpeciesCum_, localDsmcNSpeciesCum);
            restoreFieldList(nSpeciesCum_, localNSpeciesCum);
            restoreFieldList(dsmcMccSpeciesCum_, localDsmcMccSpeciesCum);
            restoreFieldList
            (
                dsmcNGrndElecLvlSpeciesCum_,
                localDsmcNGrndElecLvlSpeciesCum
            );
            restoreFieldList
            (
                dsmcN1stElecLvlSpeciesCum_,
                localDsmcN1stElecLvlSpeciesCum
            );
            restoreFieldListList
            (
                dsmcSpeciesEvibModCum_,
                localDsmcSpeciesEvibModCum
            );

            restoreFieldList(rhoNBF_, localRhoNBF);
            restoreFieldList(rhoMBF_, localRhoMBF);
            restoreFieldList(linearKEBF_, localLinearKEBF);
            restoreFieldList(ErotBF_, localErotBF);
            restoreFieldList(zetaRotBF_, localZetaRotBF);
            restoreFieldList(qBF_, localQBF);
            restoreFieldList(zetaVibBF_, localZetaVibBF);
            restoreFieldList(rhoNIntBF_, localRhoNIntBF);
            restoreFieldList(rhoNElecBF_, localRhoNElecBF);
            restoreFieldList(momentumBF_, localMomentumBF);
            restoreFieldList(fDBF_, localFDBF);

            restoreFieldListList(speciesEvibBF_, localSpeciesEvibBF);
            restoreFieldListList(speciesEelecBF_, localSpeciesEelecBF);
            restoreFieldListList(speciesRhoNBF_, localSpeciesRhoNBF);
            restoreFieldListList(speciesMccBF_, localSpeciesMccBF);
            restoreFieldListListList
            (
                speciesEvibModBF_,
                localSpeciesEvibModBF
            );

            if (doProfile)
            {
                profileOutputRestoreWallTime_ +=
                    wallSeconds(outputRestoreStart, wallClockNow());
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
        profileDetailEnabled_
	     && cloud_.isOutputRank()
	     && !finalProfilePrinted_
	     && time_.time().value() >= time_.time().endTime().value() - SMALL
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
            << "    output snapshot       = " << profileOutputSnapshotWallTime_ << " s" << nl
            << "    output reduce         = " << profileOutputReduceWallTime_ << " s" << nl
            << "    output compute        = " << profileOutputComputeWallTime_ << " s" << nl
            << "    field writes          = " << profileFieldWriteWallTime_ << " s" << nl
            << "    output reset          = " << profileOutputResetWallTime_ << " s" << nl
            << "    output restore        = " << profileOutputRestoreWallTime_ << " s" << nl
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
    // M4: align the cumulative storage with the current owned-cell set
    // before resizing (also invalidates any stale owner version).
    initOwnedStorage();

    const label nCells =
        ownedStorageActive_
      ? cloud_.replicatedMesh().myCells().size()
      : mesh_.nCells();

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
    if (measuredCollisionRate_.empty())
    {
        measuredCollisionRate_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "measuredCollisionRate_"+ fieldName_,
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0), 0.0)
            )
        );
    }
    measuredCollisionRate_() = 0.0;
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

    const dictionary& controlDict = mesh_.time().controlDict();

    const bool hasManualForceMomentCoefficients =
        controlDict.found("forceMomentCoefficients");
    const bool writeAutomaticForceMomentCoefficients =
        controlDict.lookupOrDefault<bool>
        (
            "writeForceMomentCoefficients",
            true
        );

    // Force, moment, and coefficients are enabled by default.  Explicit
    // switches remain available for cases that need to disable the output.
    writeForceMoment_ =
        controlDict.lookupOrDefault<bool>("writeForceMoment", true)
     || hasManualForceMomentCoefficients;

    writeForceMoment_ = writeForceMoment_ && fieldName_ == "mixture";

    const bool hasForceMomentReferencePoint =
        controlDict.found("forceMomentReferencePoint");
    const bool hasForceMomentReferenceArea = controlDict.found("Aref");
    const bool hasForceMomentReferenceLength = controlDict.found("lRef");

    forceMomentReferencePoint_ = controlDict.lookupOrDefault<vector>
    (
        "forceMomentReferencePoint",
        vector::zero
    );

    writeForceMomentCoefficients_ = false;
    writeWallPressureCoefficient_ = false;
    forceMomentQInf_ = 0.0;
    forceMomentPInf_ = 0.0;
    forceMomentHeatFluxInf_ = 0.0;
    forceMomentDragDirection_ = vector::zero;
    forceMomentLiftDirection_ = vector::zero;
    forceMomentSideDirection_ = vector::zero;
    forceMomentReferencePointAutomatic_ = false;

    if
    (
        hasManualForceMomentCoefficients
     || !writeAutomaticForceMomentCoefficients
    )
    {
        forceMomentReferenceArea_ = 0.0;
        forceMomentReferenceLength_ = 0.0;
        forceMomentFreestreamSource_ = word::null;
        forceMomentReferenceAreaDefinition_ = word::null;
        forceMomentReferenceLengthDefinition_ = word::null;
        forceMomentAutomaticReferencePoint_ = vector::zero;
        forceMomentReferenceGeometryDirection_ = vector::zero;
        forceMomentReferenceGeometryCached_ = false;
    }

    const label nBoundaryPatches = mesh_.boundaryMesh().size();
    forceMomentWallPatchIds_.clear();
    forceMomentWallPatchIndices_.setSize(nBoundaryPatches, -1);

    if (writeForceMoment_)
    {
        label nWallPatches = 0;

        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (isA<wallPolyPatch>(mesh_.boundaryMesh()[patchi]))
            {
                ++nWallPatches;
            }
        }

        forceMomentWallPatchIds_.setSize(nWallPatches);
        nWallPatches = 0;

        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (isA<wallPolyPatch>(mesh_.boundaryMesh()[patchi]))
            {
                forceMomentWallPatchIds_[nWallPatches] = patchi;
                forceMomentWallPatchIndices_[patchi] = nWallPatches;
                ++nWallPatches;
            }
        }
    }

    if
    (
        writeForceMoment_
     &&
        (
            hasManualForceMomentCoefficients
         || writeAutomaticForceMomentCoefficients
        )
    )
    {
        scalar rhoInf = 0.0;
        vector UInf(vector::zero);
        scalar pInf = 0.0;
        bool hasLiftDirection = false;
        vector liftDirection(vector(0, 1, 0));

        if (hasManualForceMomentCoefficients)
        {
            const dictionary& coefficientDict =
                controlDict.subDict("forceMomentCoefficients");

            rhoInf = readScalar(coefficientDict.lookup("rhoInf"));
            UInf = vector(coefficientDict.lookup("UInf"));
            pInf = coefficientDict.lookupOrDefault<scalar>("pInf", 0.0);
            forceMomentReferenceArea_ =
                readScalar(coefficientDict.lookup("Aref"));
            forceMomentReferenceLength_ =
                readScalar(coefficientDict.lookup("lRef"));
            hasLiftDirection = coefficientDict.found("liftDir");
            liftDirection = coefficientDict.lookupOrDefault<vector>
            (
                "liftDir",
                liftDirection
            );
            forceMomentFreestreamSource_ = "userSpecified";
            forceMomentReferenceAreaDefinition_ = "userSpecified";
            forceMomentReferenceLengthDefinition_ = "userSpecified";
        }
        else
        {
            if
            (
               !readAutomaticFreestream
                (
                    mesh_,
                    cloud_,
                    rhoInf,
                    UInf,
                    pInf,
                    forceMomentFreestreamSource_
                )
            )
            {
                FatalErrorInFunction
                    << "Cannot determine free-stream density and velocity. "
                    << "Define dsmcFreeStreamInflowPatch or "
                    << "dsmcChapmanEnskogFreeStreamInflowPatch in "
                    << "boundariesDict, or provide velocity and "
                    << "numberDensities in the first dsmcInitialiseDict "
                    << "configuration." << exit(FatalError);
            }

            hasLiftDirection = controlDict.found("forceMomentLiftDirection");
            liftDirection = controlDict.lookupOrDefault<vector>
            (
                "forceMomentLiftDirection",
                liftDirection
            );
        }

        const scalar magUInf = mag(UInf);
        if (rhoInf <= VSMALL || magUInf <= VSMALL)
        {
            FatalErrorInFunction
                << "Force-moment coefficients require positive free-stream "
                << "density and velocity." << exit(FatalError);
        }

        forceMomentDragDirection_ = UInf/magUInf;

        if (!hasManualForceMomentCoefficients)
        {
            const bool needsAutomaticReferenceGeometry =
                   !hasForceMomentReferencePoint
                || !hasForceMomentReferenceArea
                || !hasForceMomentReferenceLength;

            if
            (
                needsAutomaticReferenceGeometry
             &&
                (
                   !forceMomentReferenceGeometryCached_
                 || mag
                    (
                        forceMomentReferenceGeometryDirection_
                      - forceMomentDragDirection_
                    ) > SMALL
                )
            )
            {
                calculateAutomaticReferenceGeometry
                (
                    mesh_,
                    forceMomentWallPatchIds_,
                    forceMomentDragDirection_,
                    controlDict.lookupOrDefault<bool>("replicatedMesh", false),
                    forceMomentReferenceArea_,
                    forceMomentReferenceLength_,
                    forceMomentAutomaticReferencePoint_,
                    forceMomentReferenceAreaDefinition_,
                    forceMomentReferenceLengthDefinition_
                );
                forceMomentReferenceGeometryDirection_ =
                    forceMomentDragDirection_;
                forceMomentReferenceGeometryCached_ = true;
            }

            if (!hasForceMomentReferencePoint)
            {
                forceMomentReferencePoint_ = forceMomentAutomaticReferencePoint_;
                forceMomentReferencePointAutomatic_ = true;
            }

        }

        if (hasForceMomentReferenceArea)
        {
            forceMomentReferenceArea_ = readScalar(controlDict.lookup("Aref"));
            forceMomentReferenceAreaDefinition_ = "userSpecified";
        }

        if (hasForceMomentReferenceLength)
        {
            forceMomentReferenceLength_ = readScalar(controlDict.lookup("lRef"));
            forceMomentReferenceLengthDefinition_ = "userSpecified";
        }

        if
        (
            forceMomentReferenceArea_ <= VSMALL
         || forceMomentReferenceLength_ <= VSMALL
        )
        {
            FatalErrorInFunction
                << "Force-moment coefficients require positive reference area "
                << "and reference length." << exit(FatalError);
        }

        if (!hasLiftDirection)
        {
            const vector& dragDirection = forceMomentDragDirection_;

            if
            (
                   mag(dragDirection.y()) <= mag(dragDirection.x())
                && mag(dragDirection.y()) <= mag(dragDirection.z())
            )
            {
                liftDirection = vector(0, 1, 0);
            }
            else if (mag(dragDirection.z()) <= mag(dragDirection.x()))
            {
                liftDirection = vector(0, 0, 1);
            }
            else
            {
                liftDirection = vector(1, 0, 0);
            }
        }

        liftDirection -=
            (liftDirection & forceMomentDragDirection_)
           *forceMomentDragDirection_;

        if (mag(liftDirection) <= VSMALL)
        {
            FatalErrorInFunction
                << "forceMomentLiftDirection must not be parallel to the "
                << "free-stream velocity." << exit(FatalError);
        }

        forceMomentLiftDirection_ = liftDirection/mag(liftDirection);
        forceMomentSideDirection_ =
            forceMomentDragDirection_ ^ forceMomentLiftDirection_;
        forceMomentQInf_ = 0.5*rhoInf*sqr(magUInf);
        forceMomentPInf_ = pInf;
        forceMomentHeatFluxInf_ = forceMomentQInf_*magUInf;
        writeForceMomentCoefficients_ = true;
        writeWallPressureCoefficient_ = pInf > VSMALL;
    }

    if (writeWallPressureCoefficient_ && !Cp_.valid())
    {
        Cp_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "wallPressureCoefficient_"+ fieldName_,
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
    else if (!writeWallPressureCoefficient_ && Cp_.valid())
    {
        Cp_.clear();
    }

    if (writeForceMomentCoefficients_ && !Ch_.valid())
    {
        Ch_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "wallHeatFluxCoefficient_"+ fieldName_,
                    time_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
    else if (!writeForceMomentCoefficients_ && Ch_.valid())
    {
        Ch_.clear();
    }
}

} // End namespace Foam

// ************************************************************************** //
