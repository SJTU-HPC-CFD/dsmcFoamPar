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

Class
    Foam::dsmcParcel

Description
    Minimal stage-1 DSMC parcel used for the v2506 migration.

SourceFiles
    dsmcParcelI.H
    dsmcParcel.C
    dsmcParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcParcel_H
#define dsmcParcel_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"

namespace Foam
{

class dsmcCloud;
class dsmcParcel;

Ostream& operator<<(Ostream&, const dsmcParcel&);

class dsmcParcel
:
    public particle
{
public:

    static const std::size_t sizeofFields;

    class constantProperties
    {
        label type_;
        scalar mass_;
        scalar d_;
        scalar rotationalDegreesOfFreedom_;
        label nVibrationalModes_;
        labelList vibrationalDegeneracyList_;
        scalar omega_;
        scalar alpha_;
        scalarList thetaV_;
        scalar thetaD_;
        scalarList Zref_;
        scalarList TrefZv_;
        labelList charDissQuantumLevel_;
        scalar ionisationTemperature_;
        label charge_;
        label nElectronicLevels_;
        scalarList electronicEnergyList_;
        labelList electronicDegeneracyList_;

    public:
        inline constantProperties();
        inline constantProperties(const dictionary& dict);

        inline label type() const;
        inline scalar mass() const;
        inline scalar d() const;
        inline scalar sigmaT() const;
        inline scalar rotationalDegreesOfFreedom() const;
        inline label nVibrationalModes() const;
        inline label vibrationalDegeneracy_m(const label mode) const;
        inline const labelList& vibrationalDegeneracyList() const;
        inline scalar omega() const;
        inline scalar alpha() const;
        inline scalar thetaV_m(const label mode) const;
        inline const scalarList& thetaV() const;
        inline scalar thetaD() const;
        inline scalar Zref_m(const label mode) const;
        inline const scalarList& Zref() const;
        inline scalar TrefZv_m(const label mode) const;
        inline const scalarList& TrefZv() const;
        inline label charDissQuantumLevel_m(const label mode) const;
        inline const labelList& charDissQuantumLevel() const;
        inline scalar eVib_m(const label mode, const label vibLevel) const;
        inline scalar eVib_tot(const labelList& vibLevel) const;
        inline scalarList eVib(const labelList& vibLevel) const;
        inline scalar ionisationTemperature() const;
        inline label charge() const;
        inline label nElectronicLevels() const;
        inline const scalarList& electronicEnergyList() const;
        inline const labelList& electronicDegeneracyList() const;
    };

    class StuckParcel
    {
        scalarField wallTemperature_;
        vectorField wallVectors_;

    public:
        StuckParcel()
        :
            wallTemperature_(4, 0.0),
            wallVectors_(4, vector::zero)
        {}

        StuckParcel(const scalarField& wallTemperature, const vectorField& wallVectors)
        :
            wallTemperature_(wallTemperature),
            wallVectors_(wallVectors)
        {}

        const scalarField& wallTemperature() const { return wallTemperature_; }
        const vectorField& wallVectors() const { return wallVectors_; }
        scalarField& wallTemperature() { return wallTemperature_; }
        vectorField& wallVectors() { return wallVectors_; }
    };


    class TrackedParcel
    {
    public:
        static label nDELETED;
    private:
        bool tracked_;
        label inPatchId_;
        bool storePositions_;
        scalar t0_;
        vector initialPosition_;
        vector currentPosition_;
        vector distanceTravelled_;
        vector msd_;

    public:
        TrackedParcel()
        : tracked_(false), inPatchId_(-1), storePositions_(false), t0_(0.0),
          initialPosition_(vector::zero), currentPosition_(vector::zero),
          distanceTravelled_(vector::zero), msd_(vector::zero)
        {}

        TrackedParcel(const bool& tracked, const label& inPatchId, const scalar& initialTime, const vector& initialPosition)
        : tracked_(tracked), inPatchId_(inPatchId), storePositions_(false), t0_(initialTime),
          initialPosition_(initialPosition), currentPosition_(initialPosition),
          distanceTravelled_(vector::zero), msd_(vector::zero)
        {}

        TrackedParcel(const bool& tracked, const label& inPatchId, const bool& storePositions, const scalar& initialTime, const vector& initialPosition)
        : tracked_(tracked), inPatchId_(inPatchId), storePositions_(storePositions), t0_(initialTime),
          initialPosition_(initialPosition), currentPosition_(initialPosition),
          distanceTravelled_(vector::zero), msd_(vector::zero)
        {}

        TrackedParcel(const bool& tracked, const label& inPatchId, const scalar& initialTime, const vector& initialPosition, const vector& distanceTravelled)
        : tracked_(tracked), inPatchId_(inPatchId), storePositions_(false), t0_(initialTime),
          initialPosition_(initialPosition), currentPosition_(initialPosition),
          distanceTravelled_(distanceTravelled), msd_(vector::zero)
        {}

        TrackedParcel(const bool& tracked, const label& inPatchId, const bool& storePositions, const scalar& initialTime, const vector& initialPosition, const vector& distanceTravelled)
        : tracked_(tracked), inPatchId_(inPatchId), storePositions_(storePositions), t0_(initialTime),
          initialPosition_(initialPosition), currentPosition_(initialPosition),
          distanceTravelled_(distanceTravelled), msd_(vector::zero)
        {}

        const bool& tracked() const { return tracked_; }
        bool& tracked() { return tracked_; }
        label inPatchId() const { return inPatchId_; }
        label& inPatchId() { return inPatchId_; }
        const bool& storePositions() const { return storePositions_; }
        bool& storePositions() { return storePositions_; }
        const scalar& initialTime() const { return t0_; }
        scalar& initialTime() { return t0_; }
        const vector& initialPosition() const { return initialPosition_; }
        vector& initialPosition() { return initialPosition_; }
        const vector& currentPosition() const { return currentPosition_; }
        vector& currentPosition() { return currentPosition_; }
        scalar distanceTravelled() const { return mag(distanceTravelled_); }
        scalar meanSquaredDisplacement() const { return msd_[0] + msd_[1] + msd_[2]; }
        const vector& distanceTravelledVector() const { return distanceTravelled_; }
        vector& distanceTravelledVector() { return distanceTravelled_; }
        const vector& meanSquaredDisplacementVector() const { return msd_; }
        vector& meanSquaredDisplacementVector() { return msd_; }

        void updateCurrentPosition(const vector& currentPos) { currentPosition_ = currentPos; }
        void updateDistanceTravelled(const vector& newPos)
        {
            distanceTravelled_ += newPos - currentPosition_;
            updateCurrentPosition(newPos);
        }
        void updateTotalDistanceTravelled(const vector& newPos)
        {
            distanceTravelled_ += vector(mag(newPos[0] - currentPosition_[0]), mag(newPos[1] - currentPosition_[1]), mag(newPos[2] - currentPosition_[2]));
            updateCurrentPosition(newPos);
        }
        void performSpecularReflectionOnDistanceTravelled(const vector& nw)
        {
            forAll(nw, dim)
            {
                if (mag(nw[dim]) > SMALL)
                {
                    distanceTravelled_[dim] *= -1;
                }
            }
        }
        void updateMeanSquaredDisplacement()
        {
            msd_ = vector(sqr(distanceTravelled_[0]), sqr(distanceTravelled_[1]), sqr(distanceTravelled_[2]));
        }
    };
    class trackingData
    :
        public particle::trackingData
    {
    public:
        trackingData(const dsmcCloud& cloud)
        :
            particle::trackingData(cloud)
        {}
    };

protected:

    vector U_;
    scalar RWF_;
    scalar ERot_;
    label ELevel_;
    label typeId_;
    label newParcel_;
    label classification_;
    TrackedParcel tracked_;
    StuckParcel* stuck_;
    labelList vibLevel_;

public:

    TypeName("dsmcParcel");

    friend class Cloud<dsmcParcel>;

    inline dsmcParcel
    (
        const polyMesh& mesh,
        const barycentric& coordinates,
        const label celli,
        const label tetFacei,
        const label tetPti,
        const vector& U,
        const scalar RWF,
        const scalar ERot,
        const label ELevel,
        const label typeId,
        const label newParcel,
        const label classification,
        const labelList& vibLevel
    );

    inline dsmcParcel
    (
        const polyMesh& mesh,
        const vector& position,
        const label celli,
        const vector& U,
        const scalar RWF,
        const scalar ERot,
        const label ELevel,
        const label typeId,
        const label newParcel,
        const label classification,
        const labelList& vibLevel
    );

    inline dsmcParcel
    (
        const polyMesh& mesh,
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
    );

    dsmcParcel
    (
        const polyMesh& mesh,
        Istream& is,
        bool readFields = true,
        bool newFormat = true
    );

    dsmcParcel(const dsmcParcel& dP, const polyMesh& mesh);
    dsmcParcel(const dsmcParcel& dP);
    ~dsmcParcel();

    virtual autoPtr<particle> clone() const
    {
        return particle::Clone(*this);
    }

    class iNew
    {
        const polyMesh& mesh_;

    public:
        iNew(const polyMesh& mesh)
        :
            mesh_(mesh)
        {}

        autoPtr<dsmcParcel> operator()(Istream& is) const
        {
            return autoPtr<dsmcParcel>(new dsmcParcel(mesh_, is, true));
        }
    };

    inline label typeId() const;
    inline const vector& U() const;
    inline scalar RWF() const;
    inline scalar ERot() const;
    inline label ELevel() const;
    inline label newParcel() const;
    inline const TrackedParcel& tracked() const;
    inline bool isTracked() const;
    inline label classification() const;
    inline const StuckParcel& stuck() const;
    inline bool isStuck() const;
    inline bool isFree() const;
    inline const labelList& vibLevel() const;

    inline label& typeId();
    inline vector& U();
    inline scalar& RWF();
    inline scalar& ERot();
    inline label& ELevel();
    inline label& newParcel();
    inline TrackedParcel& tracked();
    inline void setTracked();
    inline void setTracked(const bool& tracked, const label& inPatchId, const scalar& initialTime, const vector& initialPosition);
    inline void setTracked(const bool& tracked, const label& inPatchId, const bool& storePositions, const scalar& initialTime, const vector& initialPosition);
    inline void setTracked(const bool& tracked, const label& inPatchId, const scalar& initialTime, const vector& initialPosition, const vector& distanceTravelled);
    inline void setTracked(const bool& tracked, const label& inPatchId, const bool& storePositions, const scalar& initialTime, const vector& initialPosition, const vector& distanceTravelled);
    inline void deleteTracked();
    inline label& classification();
    inline StuckParcel& stuck();
    inline void setStuck();
    inline void setStuck(const scalarField& wallTemperature, const vectorField& wallVectors);
    inline void deleteStuck();
    inline labelList& vibLevel();

    bool move(dsmcCloud& cloud, trackingData& td, const scalar trackTime);
    bool hitPatch(dsmcCloud& cloud, trackingData& td);
    void hitProcessorPatch(dsmcCloud& cloud, trackingData& td);
    void hitWallPatch(dsmcCloud& cloud, trackingData& td);

    void transformProperties(const tensor& T);
    void transformProperties(const vector& separation);

    static void readFields(Cloud<dsmcParcel>& c);
    static void writeFields(const Cloud<dsmcParcel>& c);

private:
    dsmcParcel& operator=(const dsmcParcel&) = delete;
};

} // End namespace Foam

#include "dsmcParcelI.H"

#endif

// ************************************************************************* //


