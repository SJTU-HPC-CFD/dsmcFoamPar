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

#include "dsmcParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

const std::size_t Foam::dsmcParcel::sizeofFields
(
    sizeof(vector)
  + sizeof(scalar)
  + sizeof(scalar)
  + sizeof(label)
  + sizeof(label)
  + sizeof(label)
  + sizeof(label)
  + sizeof(label)
);

Foam::dsmcParcel::dsmcParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    U_(Zero),
    RWF_(1.0),
    ERot_(0.0),
    ELevel_(0),
    typeId_(-1),
    newParcel_(-1),
    classification_(0),
    localCellI_(-1),
    tracked_(),
    stuck_(nullptr),
    vibLevel_()
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            label stuckFlag = 0;
            scalarField wallTemperature;
            vectorField wallVectors;
            is >> U_ >> RWF_ >> ERot_ >> ELevel_ >> typeId_ >> newParcel_ >> classification_ >> vibLevel_ >> stuckFlag;
            if (stuckFlag)
            {
                is >> wallTemperature >> wallVectors;
                setStuck(wallTemperature, wallVectors);
            }
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            is.beginRawRead();
            readRawScalar(is, U_.data(), vector::nComponents);
            readRawScalar(is, &RWF_);
            readRawScalar(is, &ERot_);
            readRawLabel(is, &ELevel_);
            readRawLabel(is, &typeId_);
            readRawLabel(is, &newParcel_);
            readRawLabel(is, &classification_);
            is.endRawRead();
            label stuckFlag = 0;
            scalarField wallTemperature;
            vectorField wallVectors;
            is >> vibLevel_ >> stuckFlag;
            if (stuckFlag)
            {
                is >> wallTemperature >> wallVectors;
                setStuck(wallTemperature, wallVectors);
            }
        }
        else
        {
            is.read(reinterpret_cast<char*>(&U_), sizeof(U_));
            is.read(reinterpret_cast<char*>(&RWF_), sizeof(RWF_));
            is.read(reinterpret_cast<char*>(&ERot_), sizeof(ERot_));
            is.read(reinterpret_cast<char*>(&ELevel_), sizeof(ELevel_));
            is.read(reinterpret_cast<char*>(&typeId_), sizeof(typeId_));
            is.read(reinterpret_cast<char*>(&newParcel_), sizeof(newParcel_));
            is.read(reinterpret_cast<char*>(&classification_), sizeof(classification_));
            label stuckFlag = 0;
            scalarField wallTemperature;
            vectorField wallVectors;
            is >> vibLevel_ >> stuckFlag;
            if (stuckFlag)
            {
                is >> wallTemperature >> wallVectors;
                setStuck(wallTemperature, wallVectors);
            }
        }
    }

    is.check(FUNCTION_NAME);
}

Foam::dsmcParcel::dsmcParcel(const dsmcParcel& dP, const polyMesh& mesh)
:
    particle(dP, mesh),
    U_(dP.U_),
    RWF_(dP.RWF_),
    ERot_(dP.ERot_),
    ELevel_(dP.ELevel_),
    typeId_(dP.typeId_),
    newParcel_(dP.newParcel_),
    classification_(dP.classification_),
    localCellI_(-1),
    tracked_(dP.tracked_),
    stuck_(nullptr),
    vibLevel_(dP.vibLevel_)
{
    if (dP.stuck_)
    {
        setStuck(dP.stuck().wallTemperature(), dP.stuck().wallVectors());
    }
}

Foam::dsmcParcel::dsmcParcel(const dsmcParcel& dP)
:
    dsmcParcel(dP, dP.mesh())
{}

Foam::dsmcParcel::~dsmcParcel()
{
    deleteStuck();
}

void Foam::dsmcParcel::readFields(Cloud<dsmcParcel>& c)
{
    const bool readOnProc = c.size();

    particle::readFields(c);

    IOField<vector> U(c.newIOobject("U", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, U);

    IOField<scalar> RWF(c.newIOobject("radialWeight", IOobject::READ_IF_PRESENT), c.size());
    if (RWF.size() != c.size())
    {
        RWF.setSize(c.size());
        RWF = scalarField(c.size(), 1.0);
    }
    c.checkFieldIOobject(c, RWF);

    IOField<scalar> ERot(c.newIOobject("ERot", IOobject::READ_IF_PRESENT), c.size());
    if (ERot.size() != c.size())
    {
        ERot.setSize(c.size());
        ERot = scalarField(c.size(), 0.0);
    }
    c.checkFieldIOobject(c, ERot);

    IOField<label> ELevel(c.newIOobject("ELevel", IOobject::READ_IF_PRESENT), c.size());
    if (ELevel.size() != c.size())
    {
        ELevel.setSize(c.size());
        ELevel = labelField(c.size(), 0);
    }
    c.checkFieldIOobject(c, ELevel);

    IOField<label> typeId(c.newIOobject("typeId", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, typeId);

    IOField<label> newParcel(c.newIOobject("newParcel", IOobject::READ_IF_PRESENT), c.size());
    if (newParcel.size() != c.size())
    {
        newParcel.setSize(c.size());
        newParcel = labelField(c.size(), -1);
    }
    c.checkFieldIOobject(c, newParcel);

    IOField<label> classification(c.newIOobject("classification", IOobject::READ_IF_PRESENT), c.size());
    if (classification.size() != c.size())
    {
        classification.setSize(c.size());
        classification = labelField(c.size(), 0);
    }
    c.checkFieldIOobject(c, classification);

    IOField<labelField> vibLevel(c.newIOobject("vibLevel", IOobject::READ_IF_PRESENT), c.size());
    IOField<label> stuckFlag(c.newIOobject("stuck", IOobject::READ_IF_PRESENT), c.size());
    IOField<scalarField> stuckWallTemperature(c.newIOobject("stuckWallTemperature", IOobject::READ_IF_PRESENT), c.size());
    IOField<vectorField> stuckWallVectors(c.newIOobject("stuckWallVectors", IOobject::READ_IF_PRESENT), c.size());
    if (vibLevel.size() != c.size())
    {
        vibLevel.setSize(c.size());
        forAll(vibLevel, i)
        {
            vibLevel[i].clear();
        }
    }
    if (stuckFlag.size() != c.size()) { stuckFlag.setSize(c.size()); stuckFlag = labelField(c.size(), 0); }
    if (stuckWallTemperature.size() != c.size()) { stuckWallTemperature.setSize(c.size()); forAll(stuckWallTemperature, i) { stuckWallTemperature[i] = scalarField(); } }
    if (stuckWallVectors.size() != c.size()) { stuckWallVectors.setSize(c.size()); forAll(stuckWallVectors, i) { stuckWallVectors[i] = vectorField(); } }
    c.checkFieldIOobject(c, vibLevel);

    label i = 0;
    for (dsmcParcel& p : c)
    {
        p.U_ = U[i];
        p.RWF_ = RWF[i];
        p.ERot_ = ERot[i];
        p.ELevel_ = ELevel[i];
        p.typeId_ = typeId[i];
        p.newParcel_ = newParcel[i];
        p.classification_ = classification[i];
        p.vibLevel_ = vibLevel[i];
        p.deleteStuck();
        if (stuckFlag[i])
        {
            p.setStuck(stuckWallTemperature[i], stuckWallVectors[i]);
        }
        ++i;
    }
}

void Foam::dsmcParcel::writeFields(const Cloud<dsmcParcel>& c)
{
    particle::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<vector> U(c.newIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> RWF(c.newIOobject("radialWeight", IOobject::NO_READ), np);
    IOField<scalar> ERot(c.newIOobject("ERot", IOobject::NO_READ), np);
    IOField<label> ELevel(c.newIOobject("ELevel", IOobject::NO_READ), np);
    IOField<label> typeId(c.newIOobject("typeId", IOobject::NO_READ), np);
    IOField<label> newParcel(c.newIOobject("newParcel", IOobject::NO_READ), np);
    IOField<label> classification(c.newIOobject("classification", IOobject::NO_READ), np);
    IOField<labelField> vibLevel(c.newIOobject("vibLevel", IOobject::NO_READ), np);
    IOField<label> stuckFlag(c.newIOobject("stuck", IOobject::NO_READ), np);
    IOField<scalarField> stuckWallTemperature(c.newIOobject("stuckWallTemperature", IOobject::NO_READ), np);
    IOField<vectorField> stuckWallVectors(c.newIOobject("stuckWallVectors", IOobject::NO_READ), np);

    label i = 0;
    for (const dsmcParcel& p : c)
    {
        U[i] = p.U();
        RWF[i] = p.RWF();
        ERot[i] = p.ERot();
        ELevel[i] = p.ELevel();
        typeId[i] = p.typeId();
        newParcel[i] = p.newParcel();
        classification[i] = p.classification();
        vibLevel[i] = p.vibLevel();
        stuckFlag[i] = p.isStuck() ? 1 : 0;
        if (p.isStuck())
        {
            stuckWallTemperature[i] = p.stuck().wallTemperature();
            stuckWallVectors[i] = p.stuck().wallVectors();
        }
        ++i;
    }

    U.write(writeOnProc);
    RWF.write(writeOnProc);
    ERot.write(writeOnProc);
    ELevel.write(writeOnProc);
    typeId.write(writeOnProc);
    newParcel.write(writeOnProc);
    classification.write(writeOnProc);
    vibLevel.write(writeOnProc);
    stuckFlag.write(writeOnProc);
    stuckWallTemperature.write(writeOnProc);
    stuckWallVectors.write(writeOnProc);
}

Foam::Ostream& Foam::operator<<(Ostream& os, const dsmcParcel& p)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.U()
            << token::SPACE << p.RWF()
            << token::SPACE << p.ERot()
            << token::SPACE << p.ELevel()
            << token::SPACE << p.typeId()
            << token::SPACE << p.newParcel()
            << token::SPACE << p.classification()
            << token::SPACE << p.vibLevel()
            << token::SPACE << label(p.isStuck());
        if (p.isStuck())
        {
            os << token::SPACE << p.stuck().wallTemperature()
               << token::SPACE << p.stuck().wallVectors();
        }
    }
    else
    {
        const vector U = p.U();
        const scalar RWF = p.RWF();
        const scalar ERot = p.ERot();
        const label ELevel = p.ELevel();
        const label typeId = p.typeId();
        const label newParcel = p.newParcel();
        const label classification = p.classification();
        const label stuckFlag = p.isStuck() ? 1 : 0;

        os << static_cast<const particle&>(p);
        os.write(reinterpret_cast<const char*>(&U), sizeof(U));
        os.write(reinterpret_cast<const char*>(&RWF), sizeof(RWF));
        os.write(reinterpret_cast<const char*>(&ERot), sizeof(ERot));
        os.write(reinterpret_cast<const char*>(&ELevel), sizeof(ELevel));
        os.write(reinterpret_cast<const char*>(&typeId), sizeof(typeId));
        os.write(reinterpret_cast<const char*>(&newParcel), sizeof(newParcel));
        os.write(reinterpret_cast<const char*>(&classification), sizeof(classification));
        os << p.vibLevel() << stuckFlag;
        if (p.isStuck())
        {
            os << p.stuck().wallTemperature() << p.stuck().wallVectors();
        }
    }

    os.check(FUNCTION_NAME);
    return os;
}

// ************************************************************************* //

