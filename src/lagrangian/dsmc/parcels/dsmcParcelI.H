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

#include "constants.H"

inline Foam::dsmcParcel::constantProperties::constantProperties()
:
    type_(-1),
    mass_(0),
    d_(0),
    rotationalDegreesOfFreedom_(0),
    nVibrationalModes_(0),
    vibrationalDegeneracyList_(),
    omega_(0),
    alpha_(1.0),
    thetaV_(),
    thetaD_(0),
    Zref_(),
    TrefZv_(),
    charDissQuantumLevel_(),
    ionisationTemperature_(1e10),
    charge_(0),
    nElectronicLevels_(1),
    electronicEnergyList_(1, 0.0),
    electronicDegeneracyList_(1, 1)
{}

inline Foam::dsmcParcel::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    type_(-1),
    mass_(dict.get<scalar>("mass")),
    d_(dict.get<scalar>("diameter")),
    rotationalDegreesOfFreedom_(dict.lookupOrDefault<scalar>("rotationalDegreesOfFreedom", 0)),
    nVibrationalModes_(dict.lookupOrDefault<label>("nVibrationalModes", 0)),
    vibrationalDegeneracyList_(dict.lookupOrDefault<labelList>("vibrationalDegeneracyList", labelList(nVibrationalModes_, 1))),
    omega_(dict.get<scalar>("omega")),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1.0)),
    thetaV_(dict.lookupOrDefault<scalarList>("characteristicVibrationalTemperature", scalarList())),
    thetaD_(dict.lookupOrDefault<scalar>("dissociationTemperature", 0.0)),
    Zref_(dict.lookupOrDefault<scalarList>("Zref", scalarList())),
    TrefZv_(dict.lookupOrDefault<scalarList>("referenceTempForZref", scalarList())),
    charDissQuantumLevel_(),
    ionisationTemperature_(dict.lookupOrDefault<scalar>("ionisationTemperature", 1e10)),
    charge_(dict.lookupOrDefault<label>("charge", 0)),
    nElectronicLevels_(dict.lookupOrDefault<label>("nElectronicLevels", 1)),
    electronicEnergyList_(dict.lookupOrDefault<scalarList>("electronicEnergyList", scalarList(1, 0.0))),
    electronicDegeneracyList_(dict.lookupOrDefault<labelList>("electronicDegeneracyList", labelList(1, 1)))
{
    if (thetaV_.size() != nVibrationalModes_)
    {
        FatalErrorInFunction << "Number of characteristic vibrational temperatures is "
            << thetaV_.size() << ", instead of " << nVibrationalModes_ << nl
            << exit(FatalError);
    }

    if (Zref_.size() != nVibrationalModes_)
    {
        FatalErrorInFunction << "Number of reference vibrational relaxation numbers is "
            << Zref_.size() << ", instead of " << nVibrationalModes_ << nl
            << exit(FatalError);
    }

    if (TrefZv_.size() != nVibrationalModes_)
    {
        FatalErrorInFunction << "Number of reference temperature for vibrational relaxation is "
            << TrefZv_.size() << ", instead of " << nVibrationalModes_ << nl
            << exit(FatalError);
    }

    if (electronicDegeneracyList_.size() != nElectronicLevels_)
    {
        FatalErrorInFunction << "Number of degeneracy levels should be "
            << nElectronicLevels_ << ", instead of " << electronicDegeneracyList_.size() << nl
            << exit(FatalError);
    }

    if (electronicEnergyList_.size() != nElectronicLevels_)
    {
        FatalErrorInFunction << "Number of electronic energy levels should be "
            << nElectronicLevels_ << ", instead of " << electronicEnergyList_.size() << nl
            << exit(FatalError);
    }

    charDissQuantumLevel_.setSize(nVibrationalModes_, 0);
    forAll(charDissQuantumLevel_, mode)
    {
        charDissQuantumLevel_[mode] = thetaD_/thetaV_[mode];
    }

    if (charge_ == -1)
    {
        type_ = 0;
    }
    else if (charge_ == 1)
    {
        type_ = (nVibrationalModes_ == 0 ? 11 : (nVibrationalModes_ == 1 ? 21 : 31));
    }
    else
    {
        type_ = (nVibrationalModes_ == 0 ? 10 : (nVibrationalModes_ == 1 ? 20 : 30));
    }
}

inline Foam::dsmcParcel::dsmcParcel
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
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    U_(U),
    RWF_(RWF),
    ERot_(ERot),
    ELevel_(ELevel),
    typeId_(typeId),
    newParcel_(newParcel),
    classification_(classification),
    tracked_(),
    stuck_(nullptr),
    vibLevel_(vibLevel)
{}

inline Foam::dsmcParcel::dsmcParcel
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
)
:
    particle(mesh, position, celli),
    U_(U),
    RWF_(RWF),
    ERot_(ERot),
    ELevel_(ELevel),
    typeId_(typeId),
    newParcel_(newParcel),
    classification_(classification),
    tracked_(),
    stuck_(nullptr),
    vibLevel_(vibLevel)
{}

inline Foam::dsmcParcel::dsmcParcel
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
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI, false),
    U_(U),
    RWF_(RWF),
    ERot_(ERot),
    ELevel_(ELevel),
    typeId_(typeId),
    newParcel_(newParcel),
    classification_(classification),
    tracked_(),
    stuck_(nullptr),
    vibLevel_(vibLevel)
{}

inline Foam::label Foam::dsmcParcel::constantProperties::type() const { return type_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::mass() const { return mass_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::d() const { return d_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::sigmaT() const { return constant::mathematical::pi*d_*d_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::rotationalDegreesOfFreedom() const { return rotationalDegreesOfFreedom_; }
inline Foam::label Foam::dsmcParcel::constantProperties::nVibrationalModes() const { return nVibrationalModes_; }
inline Foam::label Foam::dsmcParcel::constantProperties::vibrationalDegeneracy_m(const label mode) const { return vibrationalDegeneracyList_[mode]; }
inline const Foam::labelList& Foam::dsmcParcel::constantProperties::vibrationalDegeneracyList() const { return vibrationalDegeneracyList_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::omega() const { return omega_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::alpha() const { return alpha_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::thetaV_m(const label mode) const { return thetaV_[mode]; }
inline const Foam::scalarList& Foam::dsmcParcel::constantProperties::thetaV() const { return thetaV_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::thetaD() const { return thetaD_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::Zref_m(const label mode) const { return Zref_[mode]; }
inline const Foam::scalarList& Foam::dsmcParcel::constantProperties::Zref() const { return Zref_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::TrefZv_m(const label mode) const { return TrefZv_[mode]; }
inline const Foam::scalarList& Foam::dsmcParcel::constantProperties::TrefZv() const { return TrefZv_; }
inline Foam::label Foam::dsmcParcel::constantProperties::charDissQuantumLevel_m(const label mode) const { return charDissQuantumLevel_[mode]; }
inline const Foam::labelList& Foam::dsmcParcel::constantProperties::charDissQuantumLevel() const { return charDissQuantumLevel_; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::eVib_m(const label mode, const label vibLevel) const { return constant::physicoChemical::k.value()*thetaV_[mode]*(vibLevel + 0.5); }
inline Foam::scalar Foam::dsmcParcel::constantProperties::eVib_tot(const labelList& vibLevel) const { scalar e = 0.0; forAll(vibLevel, mode) { e += eVib_m(mode, vibLevel[mode]); } return e; }
inline Foam::scalarList Foam::dsmcParcel::constantProperties::eVib(const labelList& vibLevel) const { scalarList e(vibLevel.size(), 0.0); forAll(vibLevel, mode) { e[mode] = eVib_m(mode, vibLevel[mode]); } return e; }
inline Foam::scalar Foam::dsmcParcel::constantProperties::ionisationTemperature() const { return ionisationTemperature_; }
inline Foam::label Foam::dsmcParcel::constantProperties::charge() const { return charge_; }
inline Foam::label Foam::dsmcParcel::constantProperties::nElectronicLevels() const { return nElectronicLevels_; }
inline const Foam::scalarList& Foam::dsmcParcel::constantProperties::electronicEnergyList() const { return electronicEnergyList_; }
inline const Foam::labelList& Foam::dsmcParcel::constantProperties::electronicDegeneracyList() const { return electronicDegeneracyList_; }

inline Foam::label Foam::dsmcParcel::typeId() const { return typeId_; }
inline const Foam::vector& Foam::dsmcParcel::U() const { return U_; }
inline Foam::scalar Foam::dsmcParcel::RWF() const { return RWF_; }
inline Foam::scalar Foam::dsmcParcel::ERot() const { return ERot_; }
inline Foam::label Foam::dsmcParcel::ELevel() const { return ELevel_; }
inline Foam::label Foam::dsmcParcel::newParcel() const { return newParcel_; }
inline const Foam::dsmcParcel::TrackedParcel& Foam::dsmcParcel::tracked() const { return tracked_; }
inline bool Foam::dsmcParcel::isTracked() const { return tracked_.tracked(); }
inline Foam::label Foam::dsmcParcel::classification() const { return classification_; }
inline const Foam::dsmcParcel::StuckParcel& Foam::dsmcParcel::stuck() const { return *stuck_; }
inline bool Foam::dsmcParcel::isStuck() const { return stuck_ != nullptr; }
inline bool Foam::dsmcParcel::isFree() const { return stuck_ == nullptr; }
inline const Foam::labelList& Foam::dsmcParcel::vibLevel() const { return vibLevel_; }

inline Foam::label& Foam::dsmcParcel::typeId() { return typeId_; }
inline Foam::vector& Foam::dsmcParcel::U() { return U_; }
inline Foam::scalar& Foam::dsmcParcel::RWF() { return RWF_; }
inline Foam::scalar& Foam::dsmcParcel::ERot() { return ERot_; }
inline Foam::label& Foam::dsmcParcel::ELevel() { return ELevel_; }
inline Foam::label& Foam::dsmcParcel::newParcel() { return newParcel_; }
inline Foam::dsmcParcel::TrackedParcel& Foam::dsmcParcel::tracked() { return tracked_; }
inline void Foam::dsmcParcel::setTracked() { tracked_ = dsmcParcel::TrackedParcel(); }
inline void Foam::dsmcParcel::setTracked(const bool& tracked, const label& inPatchId, const scalar& initialTime, const vector& initialPosition) { tracked_ = dsmcParcel::TrackedParcel(tracked, inPatchId, initialTime, initialPosition); }
inline void Foam::dsmcParcel::setTracked(const bool& tracked, const label& inPatchId, const bool& storePositions, const scalar& initialTime, const vector& initialPosition) { tracked_ = dsmcParcel::TrackedParcel(tracked, inPatchId, storePositions, initialTime, initialPosition); }
inline void Foam::dsmcParcel::setTracked(const bool& tracked, const label& inPatchId, const scalar& initialTime, const vector& initialPosition, const vector& distanceTravelled) { tracked_ = dsmcParcel::TrackedParcel(tracked, inPatchId, initialTime, initialPosition, distanceTravelled); }
inline void Foam::dsmcParcel::setTracked(const bool& tracked, const label& inPatchId, const bool& storePositions, const scalar& initialTime, const vector& initialPosition, const vector& distanceTravelled) { tracked_ = dsmcParcel::TrackedParcel(tracked, inPatchId, storePositions, initialTime, initialPosition, distanceTravelled); }
inline void Foam::dsmcParcel::deleteTracked() {}
inline Foam::label& Foam::dsmcParcel::classification() { return classification_; }
inline Foam::dsmcParcel::StuckParcel& Foam::dsmcParcel::stuck() { return *stuck_; }
inline void Foam::dsmcParcel::setStuck() { if (!stuck_) { stuck_ = new StuckParcel(); } }
inline void Foam::dsmcParcel::setStuck(const scalarField& wallTemperature, const vectorField& wallVectors) { deleteStuck(); stuck_ = new StuckParcel(wallTemperature, wallVectors); }
inline void Foam::dsmcParcel::deleteStuck() { if (stuck_) { delete stuck_; stuck_ = nullptr; } }
inline Foam::labelList& Foam::dsmcParcel::vibLevel() { return vibLevel_; }

// ************************************************************************* //


