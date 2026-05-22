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

\*---------------------------------------------------------------------------*/

#include "dsmcFreeStreamInflowPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcFreeStreamInflowPatch, 0);

addToRunTimeSelectionTable
(
    dsmcGeneralBoundary, dsmcFreeStreamInflowPatch, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcFreeStreamInflowPatch::dsmcFreeStreamInflowPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    translationalTemperature_(),
    rotationalTemperature_(),
    vibrationalTemperature_(),
    electronicTemperature_(),
    numberDensities_(),
    accumulatedParcelsToInsert_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcFreeStreamInflowPatch::~dsmcFreeStreamInflowPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcFreeStreamInflowPatch::initialConfiguration()
{}

void dsmcFreeStreamInflowPatch::calculateProperties()
{

}

void dsmcFreeStreamInflowPatch::controlParcelsBeforeMove()
{
    struct InflowInsertSlot
    {
        label faceLocalI = -1;
        label typeId = -1;
        label nInsert = 0;
        scalar mass = 0.0;
        scalar mostProbableSpeed = 0.0;
        scalar sCosTheta = 0.0;
        scalar coeffA = 0.0;
        scalar coeffB = 0.0;
        scalar randomScaling = 0.0;
        bool hasNormalVelocity = false;
    };

    Random& rndGen = cloud_.rndGen();
    const bool useOpenMPInflow = false;
    const label inflowThreads =
        useOpenMPInflow ? cloud_.ompNumThreads() : label(1);

    const scalar sqrtPi = sqrt(pi);

    // compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label typeId = typeIds_[i];
        scalar mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            // Replicated mesh: only the rank that owns the adjacent cell
            // processes this boundary face (avoids duplicate injection)
            if (cloud_.replicatedMeshActive())
            {
                const label adjCell =
                    mesh_.boundaryMesh()[patchId_].faceCells()[f];
                if (!cloud_.replicatedMesh().isMyCell(adjCell))
                    continue;
            }

            const label faceI = faces_[f];
            const vector sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            const scalar deltaT =
                cloud_.deltaTValue
                (
                    mesh_.boundaryMesh()[patchId_].faceCells()[f]
                );

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    translationalTemperature_,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity_ & -sF/fA )/mostProbableSpeed;

            // From Bird eqn 4.22
            accumulatedParcelsToInsert_[i][f] +=
                (
                    fA*numberDensities_[i]*deltaT*mostProbableSpeed
                    *
                    (
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticles(patchId_, f));
        }
    }

    // print stage
    // Info <<"compute parcels to insert end "<<" procs "<< Pstream::myProcNo()<<endl;

    const vector faceVelocity = velocity_;
    const scalar faceTranslationalTemperature = translationalTemperature_;
    const scalar faceRotationalTemperature = rotationalTemperature_;
    const scalar faceVibrationalTemperature = vibrationalTemperature_;
    const scalar faceElectronicTemperature = electronicTemperature_;

    List<List<tetIndices>> faceTets(faces_.size());
    List<scalarField> faceCTriAFracs(faces_.size());
    vectorField faceNormals(faces_.size(), Zero);
    vectorField faceTangential1(faces_.size(), Zero);
    vectorField faceTangential2(faces_.size(), Zero);
    scalarField faceRWF(faces_.size(), 0.0);

    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector fC = mesh_.faceCentres()[faceI];
        const vector sF = mesh_.faceAreas()[faces_[f]];
        const scalar fA = mag(sF);

        faceTets[f] = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

        scalarField& cTriAFracs = faceCTriAFracs[f];
        cTriAFracs.setSize(faceTets[f].size(), 0.0);

        scalar previousCumulativeSum = 0.0;

        forAll(faceTets[f], triI)
        {
            const tetIndices& faceTetIs = faceTets[f][triI];

            cTriAFracs[triI] =
                faceTetIs.faceTri(mesh_).mag()/fA
              + previousCumulativeSum;

            previousCumulativeSum = cTriAFracs[triI];
        }

        if (cTriAFracs.size())
        {
            cTriAFracs.last() = 1.0;
        }

        vector n = sF;
        n /= -mag(n);

        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

        vector t2 = n^t1;
        t2 /= mag(t2);

        faceNormals[f] = n;
        faceTangential1[f] = t1;
        faceTangential2[f] = t2;
        faceRWF[f] = cloud_.coordSystem().RWF(cellI);
    }

    DynamicList<InflowInsertSlot> insertionSlots;
    label totalInsertedParcels = 0;

    // insert parcels
    forAll(faces_, f)
    {
        const vector& n = faceNormals[f];

        forAll(typeIds_, m)
        {
            const label typeId = typeIds_[m];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[m][f];

            // Number of whole particles to insert
            label nI = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nI) > rndGen.sample01<scalar>())
            {
                nI++;
            }

            faceAccumulator -= nI;

            if (nI > 0)
            {
                const scalar mass = cloud_.constProps(typeId).mass();
                const scalar mostProbableSpeed =
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTranslationalTemperature,
                        mass
                    );
                const scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;
                const scalar coeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);
                const scalar coeffB =
                    0.5*
                    (
                        1.0
                      + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                insertionSlots.append(InflowInsertSlot());
                InflowInsertSlot& slot =
                    insertionSlots[insertionSlots.size() - 1];
                slot.faceLocalI = f;
                slot.typeId = typeId;
                slot.nInsert = nI;
                slot.mass = mass;
                slot.mostProbableSpeed = mostProbableSpeed;
                slot.sCosTheta = sCosTheta;
                slot.coeffA = coeffA;
                slot.coeffB = coeffB;
                slot.randomScaling = sCosTheta < -3 ? mag(sCosTheta) + 1 : 3.0;
                slot.hasNormalVelocity = mag(faceVelocity & n) > VSMALL;
                totalInsertedParcels += nI;
            }
        }
    }

    if (!insertionSlots.size())
    {
        return;
    }

    struct GeneratedParcel
    {
        point position = Zero;
        vector U = Zero;
        scalar RWF = 0.0;
        scalar ERot = 0.0;
        label ELevel = 0;
        label cellI = -1;
        label faceI = -1;
        label typeId = -1;
        label newParcel = -1;
        labelList vibLevel;
    };

    List<DynamicList<GeneratedParcel>> generatedByThread(inflowThreads);

    forAll(generatedByThread, threadI)
    {
        generatedByThread[threadI].reserve
        (
            max(totalInsertedParcels/max(inflowThreads, label(1)), label(16))
        );
    }

    #ifdef _OPENMP
    #pragma omp parallel for if(useOpenMPInflow) num_threads(inflowThreads) schedule(dynamic, 1)
    #endif
    for (label slotI = 0; slotI < insertionSlots.size(); ++slotI)
    {
        const label threadI = useOpenMPInflow ? cloud_.currentThreadId() : 0;
        Random& threadRndGen = cloud_.rndGen();
        DynamicList<GeneratedParcel>& localGenerated = generatedByThread[threadI];
        const InflowInsertSlot& slot = insertionSlots[slotI];
        const label faceLocalI = slot.faceLocalI;
        const label faceI = faces_[faceLocalI];
        const label cellI = cells_[faceLocalI];
        const vector& n = faceNormals[faceLocalI];
        const vector& t1 = faceTangential1[faceLocalI];
        const vector& t2 = faceTangential2[faceLocalI];
        const List<tetIndices>& localFaceTets = faceTets[faceLocalI];
        const scalarField& cTriAFracs = faceCTriAFracs[faceLocalI];
        const scalar tangentialV1 = t1 & faceVelocity;
        const scalar tangentialV2 = t2 & faceVelocity;

        for (label i = 0; i < slot.nInsert; ++i)
        {
            const scalar triSelection = threadRndGen.sample01<scalar>();
            label selectedTriI = -1;

            forAll(cTriAFracs, triI)
            {
                selectedTriI = triI;

                if (cTriAFracs[triI] >= triSelection)
                {
                    break;
                }
            }

            const tetIndices& faceTetIs = localFaceTets[selectedTriI];
            point p = faceTetIs.faceTri(mesh_).randomPoint(threadRndGen);
            scalar uNormal = 0.0;

            if (slot.hasNormalVelocity)
            {
                scalar P = -1.0;

                do
                {
                    const scalar uNormalThermal =
                        slot.randomScaling
                       *(2.0*threadRndGen.sample01<scalar>() - 1.0);

                    uNormal = uNormalThermal + slot.sCosTheta;

                    if (uNormal < 0.0)
                    {
                        P = -1.0;
                    }
                    else
                    {
                        P = 2.0*uNormal/slot.coeffA
                          *exp(slot.coeffB - sqr(uNormalThermal));
                    }
                } while (P < threadRndGen.sample01<scalar>());
            }
            else
            {
                uNormal = sqrt(-log(threadRndGen.sample01<scalar>()));
            }

            const vector U =
                sqrt(physicoChemical::k.value()*faceTranslationalTemperature/slot.mass)
               *(
                    threadRndGen.GaussNormal<scalar>()*t1
                  + threadRndGen.GaussNormal<scalar>()*t2
                )
              + tangentialV1*t1
              + tangentialV2*t2
              + slot.mostProbableSpeed*uNormal*n;
            const scalar ERot = cloud_.equipartitionRotationalEnergy
            (
                faceRotationalTemperature,
                cloud_.constProps(slot.typeId).rotationalDegreesOfFreedom()
            );
            const labelList vibLevel =
                cloud_.equipartitionVibrationalEnergyLevel
            (
                faceVibrationalTemperature,
                cloud_.constProps(slot.typeId).nVibrationalModes(),
                slot.typeId
            );
            const label ELevel = cloud_.equipartitionElectronicLevel
            (
                faceElectronicTemperature,
                cloud_.constProps(slot.typeId).electronicDegeneracyList(),
                cloud_.constProps(slot.typeId).electronicEnergyList()
            );

            GeneratedParcel parcel;
            parcel.position = p;
            parcel.U = U;
            parcel.RWF = faceRWF[faceLocalI];
            parcel.ERot = ERot;
            parcel.ELevel = ELevel;
            parcel.cellI = cellI;
            parcel.faceI = faceI;
            parcel.typeId = slot.typeId;
            parcel.newParcel = patchId();
            parcel.vibLevel = vibLevel;
            localGenerated.append(parcel);
        }
    }

    forAll(generatedByThread, threadI)
    {
        const DynamicList<GeneratedParcel>& localGenerated = generatedByThread[threadI];

        forAll(localGenerated, i)
        {
            const GeneratedParcel& parcel = localGenerated[i];

            cloud_.addNewParcel
            (
                parcel.position,
                parcel.U,
                parcel.RWF,
                parcel.ERot,
                parcel.ELevel,
                parcel.cellI,
                parcel.faceI,
                0,
                parcel.typeId,
                parcel.newParcel,
                0,
                parcel.vibLevel
            );
        }
    }
}

void dsmcFreeStreamInflowPatch::controlParcelsBeforeCollisions()
{

}

void dsmcFreeStreamInflowPatch::controlParcelsAfterCollisions()
{
}

void dsmcFreeStreamInflowPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcFreeStreamInflowPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

//     setProperties();
}



void dsmcFreeStreamInflowPatch::setProperties()
{
    velocity_ = propsDict_.get<vector>("velocity");
    translationalTemperature_ = Foam::hyCompat::toScalar(Foam::hyCompat::lookup(propsDict_, "translationalTemperature"));
    rotationalTemperature_ = propsDict_.lookupOrDefault<scalar>("rotationalTemperature", 0.0);
    vibrationalTemperature_ = propsDict_.lookupOrDefault<scalar>("vibrationalTemperature", 0.0);
    electronicTemperature_ = propsDict_.lookupOrDefault<scalar>("electronicTemperature", 0.0);

    //  read in the type ids

    const List<word> molecules(Foam::hyCompat::lookup(Foam::hyCompat::lookup(propsDict_, "typeIds")));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcFreeStreamInflowPatch::dsmcFreeStreamInflowPatch()")
            << "Cannot have zero typeIds being inserted." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word moleculeName(molecules[i]);

        if(Foam::hyCompat::indexOf(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word moleculeName(moleculesReduced[i]);

        label typeId(Foam::hyCompat::indexOf(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcFreeStreamInflowPatch::dsmcFreeStreamInflowPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    // read in the mass density per specie

    const dictionary& numberDensitiesDict
    (
        propsDict_.subDict("numberDensities")
    );

    numberDensities_.clear();

    numberDensities_.setSize(typeIds_.size(), 0.0);

    forAll(numberDensities_, i)
    {
        numberDensities_[i] = Foam::hyCompat::toScalar
        (
            numberDensitiesDict, moleculesReduced[i]
        );
    }

    // set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }
}


void dsmcFreeStreamInflowPatch::setNewBoundaryFields()
{
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    nFaces_ = 0;
    patchSurfaceArea_ = 0.0;

    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        nFaces_++;
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if(Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }

   forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }
}

// scalar dsmcFreeStreamInflowPatch::calculateThermalSpeed()
// {
//     const scalar kB = physicoChemical::k.value();  // 玻尔兹曼常数
//     const scalar T = translationalTemperature_;    // 平动温度
    
//     // 计算加权平均质量
//     scalar weightedMass = 0.0;
//     scalar totalNumberDensity = 0.0;
    
//     forAll(typeIds_, i)
//     {
//         const label typeId = typeIds_[i];
//         const scalar mass = cloud_.constProps(typeId).mass();
//         const scalar nDens = numberDensities_[i];
        
//         weightedMass += mass * nDens;
//         totalNumberDensity += nDens;
//     }
    
//     // 计算平均质量
//     weightedMass /= totalNumberDensity;
    
//     // 计算�速度
//     const scalar vThermal = sqrt(8.0 * kB * T / (constant::mathematical::pi * weightedMass));
    
//     return vThermal;
// }

} // End namespace Foam

// ************************************************************************* //

