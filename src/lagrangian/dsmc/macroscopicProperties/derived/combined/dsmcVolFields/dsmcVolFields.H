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
    dsmcVolFields

Description
    Geometric field for temperature, density, velocity, pressure.

SourceFiles
    dsmcVolFields.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcVolFields_H
#define dsmcVolFields_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "HashSet.H"
#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dsmcVolFields Declaration
\*---------------------------------------------------------------------------*/

class dsmcVolFields
:
	public dsmcField
{

protected:

    // Protected data

        dictionary propsDict_;

        label sampleInterval_;
        label sampleCounter_;
        scalar nTimeSteps_;
        
        scalar mfpTref_;
        label nMinParcelsTvib_;
        scalar iMeanMinTvib_;
        word fieldName_;
        
        labelList speciesIds_;
        labelList typeIdToSpeciesIndex_;
        
        autoPtr<volScalarField> dsmcN_;
        autoPtr<volScalarField> dsmcNMean_;
        autoPtr<volScalarField> rhoN_;
        autoPtr<volScalarField> rhoM_;
        autoPtr<volScalarField> p_;
        autoPtr<volScalarField> Ttra_;
        autoPtr<volScalarField> Trot_;
        autoPtr<volScalarField> Tvib_;
        autoPtr<volScalarField> Telec_;
        autoPtr<volScalarField> Tov_;
        autoPtr<volScalarField> q_;
        autoPtr<volScalarField> tau_;
        autoPtr<volScalarField> Cp_;
        autoPtr<volScalarField> Ch_;
        autoPtr<volScalarField> mfp_;
        autoPtr<volScalarField> mfpToDx_;
        autoPtr<volScalarField> DxToMfp_;
        autoPtr<volScalarField> meanCollisionRate_;
        autoPtr<volScalarField> measuredCollisionRate_;
        autoPtr<volScalarField> meanCollisionTime_;
        autoPtr<volScalarField> mctToDt_;
        autoPtr<volScalarField> meanCollisionSeparation_;
        autoPtr<volScalarField> SOF_;
        autoPtr<volScalarField> Ma_;
        autoPtr<volScalarField> classIDistribution_;
        autoPtr<volScalarField> classIIDistribution_;
        autoPtr<volScalarField> classIIIDistribution_;
        autoPtr<volScalarField> densityError_;
        autoPtr<volScalarField> velocityError_;
        autoPtr<volScalarField> temperatureError_;
        autoPtr<volScalarField> pressureError_;
        autoPtr<volVectorField> UMean_;
        autoPtr<volVectorField> fD_;
        autoPtr<volVectorField> heatFluxVector_;
        autoPtr<volTensorField> pressureTensor_;
        autoPtr<volTensorField> shearStressTensor_;

        scalarField dsmcNCum_;
        scalarField nCum_;
        scalarField dsmcNElecLvlCum_;
        scalarField dsmcMCum_;
        scalarField mCum_;
        scalarField dsmcLinearKECum_;
        scalarField linearKECum_;
        scalarField dsmcErotCum_;
        scalarField dsmcZetaRotCum_;
        scalarField dsmcMuuCum_;
        scalarField dsmcMuvCum_;
        scalarField dsmcMuwCum_;
        scalarField dsmcMvvCum_;
        scalarField dsmcMvwCum_;
        scalarField dsmcMwwCum_;
        scalarField dsmcMccCum_;
        scalarField dsmcMccuCum_;
        scalarField dsmcMccvCum_;
        scalarField dsmcMccwCum_;
        scalarField dsmcEuCum_;
        scalarField dsmcEvCum_;
        scalarField dsmcEwCum_;
        scalarField dsmcECum_;
        scalarField zetaVib_;
        scalarField dsmcNClassICum_;
        scalarField dsmcNClassIICum_;
        scalarField dsmcNClassIIICum_;
        scalarField collisionSeparation_;
        scalarField dsmcNCollsCum_;
        vectorField dsmcMomentumCum_;
        vectorField momentumCum_;
        
        List<labelList> boundaryCells_;
        labelList sampledBoundaryPatches_;
        List<List<scalarField>> dsmcSpeciesEvibModCum_;
        List<scalarField> dsmcSpeciesEelecCum_;
        List<scalarField> dsmcNSpeciesCum_;
        List<scalarField> nSpeciesCum_;
        List<scalarField> dsmcMccSpeciesCum_;
        List<scalarField> speciesTvib_;
        List<scalarField> dsmcNGrndElecLvlSpeciesCum_;
        List<scalarField> dsmcN1stElecLvlSpeciesCum_;
        List<scalarField> speciesMfp_;
        List<scalarField> speciesMcr_;
        //- Boundary measurements
      
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> linearKEBF_;
        List<scalarField> ErotBF_;
        List<scalarField> zetaRotBF_;
        List<scalarField> qBF_;
        List<scalarField> zetaVibBF_;
        List<scalarField> rhoNIntBF_;
        List<scalarField> rhoNElecBF_;
        List<vectorField> momentumBF_;
        List<vectorField> fDBF_;

        List<List<scalarField>> speciesEvibBF_;
        List<List<scalarField>> speciesEelecBF_;
        List<List<scalarField>> speciesRhoNBF_;
        List<List<scalarField>> speciesMccBF_;
        List<List<scalarField>> speciesTvibBF_;
        List<List<scalarField>> speciesZetaVibBF_;
        List<List<List<scalarField>>> speciesEvibModBF_;

        //- Wall unit vectors
        //  temporary storage, should be done in boundaries/patch
        List<vectorField> n_;
        List<vectorField> t1_;
        List<vectorField> t2_;

        bool openmpFieldSampling_;
        bool averagingAcrossManyRuns_;
        bool measureClassifications_;
        bool measureMeanFreePath_;
        bool measureErrors_;
        bool densityOnly_;
        bool measureHeatFluxShearStress_;
        bool writeRotationalTemperature_;
        bool writeVibrationalTemperature_;
        bool writeElectronicTemperature_;
        bool restrictOutputFields_;
        HashSet<word> outputFieldNames_;
        bool writeForceMoment_;
        vector forceMomentReferencePoint_;
        labelList forceMomentWallPatchIds_;
        labelList forceMomentWallPatchIndices_;
        bool writeForceMomentCoefficients_;
        bool writeWallPressureCoefficient_;
        scalar forceMomentQInf_;
        scalar forceMomentPInf_;
        scalar forceMomentHeatFluxInf_;
        scalar forceMomentReferenceArea_;
        scalar forceMomentReferenceLength_;
        vector forceMomentDragDirection_;
        vector forceMomentLiftDirection_;
        vector forceMomentSideDirection_;
        word forceMomentFreestreamSource_;
        word forceMomentReferenceAreaDefinition_;
        word forceMomentReferenceLengthDefinition_;
        bool forceMomentReferencePointAutomatic_;
        vector forceMomentAutomaticReferencePoint_;
        vector forceMomentReferenceGeometryDirection_;
        bool forceMomentReferenceGeometryCached_;
        bool profileSummaryEnabled_;
        bool profileDetailEnabled_;
        scalar profileSampleAccumWallTime_;
        scalar profileSharedCacheBuildWallTime_;
        scalar profileSharedCacheAllocateWallTime_;
        scalar profileSharedCacheResetWallTime_;
        scalar profileSharedCacheParcelAccumWallTime_;
        scalar profileSharedCacheBaseAccumWallTime_;
        scalar profileSharedCacheVibAccumWallTime_;
        scalar profileSharedCacheElectronicAccumWallTime_;
        scalar profileSharedCacheClassAccumWallTime_;
        label profileSharedCacheDetailSampleCells_;
        label profileSharedCacheDetailSampleParcels_;
        scalar profileFieldCombineWallTime_;
        scalar profileCellReduceWallTime_;
        scalar profileBoundaryAccumWallTime_;
        scalar profileOutputSnapshotWallTime_;
        scalar profileOutputReduceWallTime_;
        scalar profileOutputComputeWallTime_;
        scalar profileFieldWriteWallTime_;
        scalar profileOutputResetWallTime_;
        scalar profileOutputRestoreWallTime_;
        scalar profileOutputTimeWallTime_;
        label profileCalls_;
        bool finalProfilePrinted_;

        //- Precomputed owned boundary face indices for replicated mesh sampling
        List<labelList> ownedBoundaryFaces_;

        //- Repartition version used to construct ownedBoundaryFaces_
        label ownedBoundaryFacesOwnerVersion_;

        //- M4 owned-only cumulative storage: when the replicated mesh runs in
        //  processor-write mode the per-cell cumulative arrays are allocated
        //  at owned-cell size and indexed through toOwnedCell_ instead of
        //  the global cell index.  GeoFields stay full-size for output.
        bool ownedStorageActive_;
        label ownedStorageOwnerVersion_;
        labelList toOwnedCell_;
        bool pendingResumeReadIn_;

        //- Allocate the cumulative arrays at owned-cell size (or full size in
        //  legacy modes) and rebuild the global-to-owned mapping.  A version
        //  mismatch after a DLB rebalance reallocates and clears the
        //  accumulators (accepted loss of one averaging window).
        void initOwnedStorage();

        //- Global cell index to cumulative-array index.  In owned mode the
        //  mapping is precomputed; otherwise it is the identity.
        inline label toOwned(const label globalCell) const
        {
            return ownedStorageActive_ ? toOwnedCell_[globalCell] : globalCell;
        }

        //- True when per-cell cumulative arrays are owned-size
        inline bool ownedStorageActive() const
        {
            return ownedStorageActive_;
        }

    // Protected Member Functions

        void readIn();

        void writeOut();

        //- M2: build the GeoFields on demand (core set unconditionally,
        //  the rest behind their measure/write switches)
        void constructOutputFields();

        void calculateWallUnitVectors();

        void updateOwnedBoundaryFaces();

        bool outputFieldEnabled
        (
            const word& fieldKey,
            const word& objectName,
            const bool defaultEnabled = true
        ) const;


public:

    //- Runtime type information
    TypeName("dsmcVolFields");

    // Constructors

        //- Construct from components
        dsmcVolFields
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

    ~dsmcVolFields();


    // Member Functions

      //- create initial field
      virtual void createField();

      //- calculate field
      virtual void calculateField();

      //- reset field if mesh has changed
      virtual void resetField();

      //- write field
      virtual void writeField();

      virtual void updateProperties(const dictionary& dict);

    // Access

      virtual scalar translationalT(const label cellI)
      {
          return Ttra_()[cellI];
      }

      virtual scalar overallT(const label cellI)
      {
          return Tov_()[cellI];
      }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
