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

#include "dsmcParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcParcel::dsmcParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    U_(vector::zero),
    RWF_(1.0),
    ERot_(0.0),
    ELevel_(0),
    typeId_(-1),
    newParcel_(-1),
    //tracked_(nullptr), // TODO TO BE REINTRODUCED
    //tracked_(0,0.0,vector::zero),
    classification_(0),
    stuck_(nullptr),
    vibLevel_(0)
{
    //dsmcParcel::TrackedParcel tP = dsmcParcel::TrackedParcel(); // TODO uncomment vincent  11/05/2018
    dsmcParcel::StuckParcel sP = dsmcParcel::StuckParcel();

    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
            RWF_ = readScalar(is);
            ERot_ = readScalar(is);
            ELevel_ = readLabel(is);
            typeId_ = readLabel(is);
            newParcel_ = readLabel(is);
            is >> tracked_; //tP; // TODO VINCENT  11/05/2018
            classification_ = readLabel(is);
            is >> sP;
            is >> vibLevel_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&U_),
                sizeof(U_)
                + sizeof(RWF_)
                + sizeof(ERot_)
                + sizeof(ELevel_)
                + sizeof(typeId_)
                + sizeof(newParcel_)
                + sizeof(classification_)
            );
            is >> tracked_; //tP;  // TODO VINCENT  11/05/2018
            is >> sP;
            is >> vibLevel_;
        }
    }

    if (false) // tP.tracked())
    {
        /*tracked_ = new dsmcParcel::TrackedParcel
            (
                tP.tracked(),
                tP.inPatchId(),
                tP.storePositions(),
                tP.initialTime(),
                tP.initialPosition(),
                tP.distanceTravelledVector()//,
                //tP.parcelTrajectory()
            );*/ // TODO TO BE REINTRODUCED

        /*tracked_ = dsmcParcel::TrackedParcel
            (
                tP.tracked(),
                tP.inPatchId(),
                tP.storePositions(),
                tP.initialTime(),
                tP.initialPosition(),
                tP.distanceTravelledVector()//,
                //tP.parcelTrajectory()
            );  */ // TODO uncomment vincent  11/05/2018
    }

    if (sP.wallTemperature()[0] != 0.0)
    {
        stuck_ = new dsmcParcel::StuckParcel
            (
                sP.wallTemperature(),
                sP.wallVectors()
            );
    }

    // Check state of Istream
    is.check
    (
        "Foam::dsmcParcel::dsmcParcel"
        "(const Cloud<dsmcParcel>& cloud, Foam::Istream&), bool"
    );
}


void Foam::dsmcParcel::readFields(Cloud<dsmcParcel>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<scalar> RWF
    (
        c.fieldIOobject
        (
            "radialWeight",
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 1.0)
    );
    c.checkFieldIOobject(c, RWF);

    IOField<scalar> ERot
    (
        c.fieldIOobject
        (
            "ERot",
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 0.0)
    );
    c.checkFieldIOobject(c, ERot);

    IOField<label> ELevel
    (
        c.fieldIOobject
        (
            "ELevel",
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, ELevel);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::MUST_READ));
    c.checkFieldIOobject(c, newParcel);

    IOField<label> classification(c.fieldIOobject("classification", IOobject::MUST_READ));
    c.checkFieldIOobject(c, classification);

    IOField<label> stuckToWall
    (
        c.fieldIOobject
        (
            "stuckToWall",
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, stuckToWall);

    IOField<scalarField> wallTemperature
    (
        c.fieldIOobject
        (
            "wallTemperature",
            IOobject::READ_IF_PRESENT
        )
    );

    if (wallTemperature.size() != c.size())
    {
        wallTemperature.setSize(c.size());
        forAll(wallTemperature, i)
        {
            wallTemperature[i] = scalarField(4, 0.0);
        }
    }

    c.checkFieldIOobject(c, wallTemperature);

    IOField<vectorField> wallVectors
    (
        c.fieldIOobject
        (
            "wallVectors",
            IOobject::READ_IF_PRESENT
        )
    );

    if (wallVectors.size() != c.size())
    {
        wallVectors.setSize(c.size());
        forAll(wallVectors, i)
        {
            wallVectors[i] = vectorField(4, vector::zero);
        }
    }

    c.checkFieldIOobject(c, wallVectors);

    IOField<label> isTracked
    (
        c.fieldIOobject
        (
            "isTracked",
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, isTracked);

    IOField<label> inPatchId
    (
        c.fieldIOobject
        (
            "inPatchId",
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), -1)
    );
    c.checkFieldIOobject(c, inPatchId);

    IOField<scalar> tracerInitialTime
    (
        c.fieldIOobject
        (
            "tracerInitialTime",
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 0.0)
    );
    c.checkFieldIOobject(c, tracerInitialTime);

    IOField<vector> tracerInitialPosition
    (
        c.fieldIOobject
        (
            "tracerInitialPosition",
            IOobject::READ_IF_PRESENT
        ),
        vectorField(c.size(), vector::zero)
    );
    c.checkFieldIOobject(c, tracerInitialPosition);

    IOField<labelField> vibLevel
    (
        c.fieldIOobject
        (
            "vibLevel",
            IOobject::READ_IF_PRESENT
        )
    );

    if (vibLevel.size() != c.size())
    {
        vibLevel.setSize(c.size());
        forAll(vibLevel, i)
        {
            vibLevel[i].setSize(0);
        }
    }

    c.checkFieldIOobject(c, vibLevel);

    label i = 0;
    forAllIter(dsmcCloud, c, iter)
    {
        dsmcParcel& p = iter();

        p.U_ = U[i];
        p.RWF_ = RWF[i];
        p.ERot_ = ERot[i];
        p.ELevel_ = ELevel[i];
        p.typeId_ = typeId[i];
        p.newParcel_ = newParcel[i];
        p.classification_ = classification[i];

        if (stuckToWall[i])
        {
            p.setStuck(wallTemperature[i], wallVectors[i]);
        }

        if (isTracked[i])
        {
            p.setTracked
            (
                isTracked[i],
                inPatchId[i],
                tracerInitialTime[i],
                tracerInitialPosition[i]
            );
        }

        p.vibLevel_ = vibLevel[i];

        i++;
    }
}


void Foam::dsmcParcel::readFieldsFiltered
(
    Cloud<dsmcParcel>& c,
    const labelList& keep,
    const label nFull
)
{
    if (!c.size())
    {
        return;
    }

    if (c.size() != keep.size())
    {
        FatalErrorInFunction
            << "Cloud size " << c.size()
            << " does not match the number of kept indices "
            << keep.size()
            << " from the filtered positions read."
            << exit(FatalIOError);
    }

    // A field read from file has nFull entries and is remapped through
    // 'keep'; a field whose file is absent keeps its default with
    // keep.size() entries.
    auto checkSize = [&](const Foam::word& name, const label size)
    {
        if (size != nFull && size != keep.size())
        {
            FatalErrorIn
            (
                "Foam::dsmcParcel::readFieldsFiltered"
            )   << "Field " << name << " has size " << size
                << " which matches neither the positions file count "
                << nFull << " nor the kept count " << keep.size()
                << exit(FatalError);
        }
    };

    // origProcId / origId (mirror particle::readFields with remap)
    {
        IOobject procIO(c.fieldIOobject("origProcId", IOobject::MUST_READ));

        if (procIO.typeHeaderOk<IOList<label>>(false))
        {
            IOField<label> origProcIdFull(procIO);
            checkSize("origProcId", origProcIdFull.size());
            IOField<label> origIdFull
            (
                c.fieldIOobject("origId", IOobject::MUST_READ)
            );
            checkSize("origId", origIdFull.size());

            label i = 0;
            forAllIter(dsmcCloud, c, iter)
            {
                particle& p = iter();

                if (origProcIdFull.size() == nFull)
                {
                    p.origProc() = origProcIdFull[keep[i]];
                }
                else
                {
                    p.origProc() = origProcIdFull[i];
                }

                if (origIdFull.size() == nFull)
                {
                    p.origId() = origIdFull[keep[i]];
                }
                else
                {
                    p.origId() = origIdFull[i];
                }

                i++;
            }
        }
    }

    // Mandatory and optional parcel fields, each remapped through keep
    IOField<vector> UFull(c.fieldIOobject("U", IOobject::MUST_READ));
    checkSize("U", UFull.size());

    IOField<scalar> RWFRead
    (
        c.fieldIOobject("radialWeight", IOobject::READ_IF_PRESENT),
        scalarField(keep.size(), 1.0)
    );
    checkSize("radialWeight", RWFRead.size());

    IOField<scalar> ERotRead
    (
        c.fieldIOobject
        (
            "ERot",
            IOobject::READ_IF_PRESENT
        ),
        scalarField(keep.size(), 0.0)
    );
    checkSize("ERot", ERotRead.size());

    IOField<label> ELevelRead
    (
        c.fieldIOobject
        (
            "ELevel",
            IOobject::READ_IF_PRESENT
        ),
        labelField(keep.size(), 0)
    );
    checkSize("ELevel", ELevelRead.size());

    IOField<label> typeIdFull(c.fieldIOobject("typeId", IOobject::MUST_READ));
    checkSize("typeId", typeIdFull.size());

    IOField<label> newParcelFull
    (
        c.fieldIOobject("newParcel", IOobject::MUST_READ)
    );
    checkSize("newParcel", newParcelFull.size());

    IOField<label> classificationFull
    (
        c.fieldIOobject("classification", IOobject::MUST_READ)
    );
    checkSize("classification", classificationFull.size());

    IOField<label> stuckToWallRead
    (
        c.fieldIOobject
        (
            "stuckToWall",
            IOobject::READ_IF_PRESENT
        ),
        labelField(keep.size(), 0)
    );
    checkSize("stuckToWall", stuckToWallRead.size());

    IOField<scalarField> wallTemperatureRead
    (
        c.fieldIOobject
        (
            "wallTemperature",
            IOobject::READ_IF_PRESENT
        )
    );

    if (wallTemperatureRead.size() != nFull)
    {
        wallTemperatureRead.setSize(keep.size());
        forAll(wallTemperatureRead, i)
        {
            wallTemperatureRead[i] = scalarField(4, 0.0);
        }
    }

    IOField<vectorField> wallVectorsRead
    (
        c.fieldIOobject
        (
            "wallVectors",
            IOobject::READ_IF_PRESENT
        )
    );

    if (wallVectorsRead.size() != nFull)
    {
        wallVectorsRead.setSize(keep.size());
        forAll(wallVectorsRead, i)
        {
            wallVectorsRead[i] = vectorField(4, vector::zero);
        }
    }

    IOField<label> isTrackedRead
    (
        c.fieldIOobject
        (
            "isTracked",
            IOobject::READ_IF_PRESENT
        ),
        labelField(keep.size(), 0)
    );
    checkSize("isTracked", isTrackedRead.size());

    IOField<label> inPatchIdRead
    (
        c.fieldIOobject
        (
            "inPatchId",
            IOobject::READ_IF_PRESENT
        ),
        labelField(keep.size(), -1)
    );
    checkSize("inPatchId", inPatchIdRead.size());

    IOField<scalar> tracerInitialTimeRead
    (
        c.fieldIOobject
        (
            "tracerInitialTime",
            IOobject::READ_IF_PRESENT
        ),
        scalarField(keep.size(), 0.0)
    );
    checkSize("tracerInitialTime", tracerInitialTimeRead.size());

    IOField<vector> tracerInitialPositionRead
    (
        c.fieldIOobject
        (
            "tracerInitialPosition",
            IOobject::READ_IF_PRESENT
        ),
        vectorField(keep.size(), vector::zero)
    );
    checkSize("tracerInitialPosition", tracerInitialPositionRead.size());

    IOField<labelField> vibLevelRead
    (
        c.fieldIOobject
        (
            "vibLevel",
            IOobject::READ_IF_PRESENT
        )
    );

    if (vibLevelRead.size() != nFull)
    {
        vibLevelRead.setSize(keep.size());
        forAll(vibLevelRead, i)
        {
            vibLevelRead[i].setSize(0);
        }
    }

    // Assignment with remap through keep
    label i = 0;
    forAllIter(dsmcCloud, c, iter)
    {
        dsmcParcel& p = iter();
        const label ki = keep[i];

        p.U_ =
            (UFull.size() == nFull) ? UFull[ki] : UFull[i];
        p.RWF_ =
            (RWFRead.size() == nFull) ? RWFRead[ki] : RWFRead[i];
        p.ERot_ =
            (ERotRead.size() == nFull) ? ERotRead[ki] : ERotRead[i];
        p.ELevel_ =
            (ELevelRead.size() == nFull) ? ELevelRead[ki] : ELevelRead[i];
        p.typeId_ =
            (typeIdFull.size() == nFull) ? typeIdFull[ki] : typeIdFull[i];
        p.newParcel_ =
            (newParcelFull.size() == nFull)
          ? newParcelFull[ki] : newParcelFull[i];
        p.classification_ =
            (classificationFull.size() == nFull)
          ? classificationFull[ki] : classificationFull[i];

        const label stuck =
            (stuckToWallRead.size() == nFull)
          ? stuckToWallRead[ki] : stuckToWallRead[i];

        if (stuck)
        {
            const scalarField wt =
                (wallTemperatureRead.size() == nFull)
              ? wallTemperatureRead[ki] : wallTemperatureRead[i];
            const vectorField wv =
                (wallVectorsRead.size() == nFull)
              ? wallVectorsRead[ki] : wallVectorsRead[i];

            p.setStuck(wt, wv);
        }

        const label tracked =
            (isTrackedRead.size() == nFull)
          ? isTrackedRead[ki] : isTrackedRead[i];

        if (tracked)
        {
            p.setTracked
            (
                tracked,
                (inPatchIdRead.size() == nFull)
              ? inPatchIdRead[ki] : inPatchIdRead[i],
                (tracerInitialTimeRead.size() == nFull)
              ? tracerInitialTimeRead[ki] : tracerInitialTimeRead[i],
                (tracerInitialPositionRead.size() == nFull)
              ? tracerInitialPositionRead[ki] : tracerInitialPositionRead[i]
            );
        }

        p.vibLevel_ =
            (vibLevelRead.size() == nFull)
          ? vibLevelRead[ki] : vibLevelRead[i];

        i++;
    }
}


void Foam::dsmcParcel::writeFields(const Cloud<dsmcParcel>& c)
{
    particle::writeFields(c);

    const label np = c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> RWF(c.fieldIOobject("radialWeight", IOobject::NO_READ), np);
    IOField<scalar> ERot(c.fieldIOobject("ERot", IOobject::NO_READ), np);
    IOField<labelField> vibLevel(c.fieldIOobject("vibLevel", IOobject::NO_READ), np);
    IOField<label> ELevel(c.fieldIOobject("ELevel", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::NO_READ), np);
    IOField<label> classification(c.fieldIOobject("classification", IOobject::NO_READ), np);

    IOField<label> stuckToWall(c.fieldIOobject("stuckToWall", IOobject::NO_READ), np);
    IOField<scalarField> wallTemperature(c.fieldIOobject("wallTemperature", IOobject::NO_READ), np);
    IOField<vectorField> wallVectors(c.fieldIOobject("wallVectors", IOobject::NO_READ), np);

    IOField<label> isTracked(c.fieldIOobject("isTracked", IOobject::NO_READ), np);
    IOField<label> inPatchId(c.fieldIOobject("inPatchId", IOobject::NO_READ), np);
    IOField<scalar> tracerInitialTime(c.fieldIOobject("tracerInitialTime", IOobject::NO_READ), np);
    IOField<vector> tracerInitialPosition(c.fieldIOobject("tracerInitialPosition", IOobject::NO_READ), np);
    IOField<vector> tracerCurrentPosition(c.fieldIOobject("tracerCurrentPosition", IOobject::NO_READ), np);
    IOField<vector> tracerDistanceTravelled(c.fieldIOobject("tracerDistanceTravelled", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(dsmcCloud, c, iter)
    {
        const dsmcParcel& p = iter();

        U[i] = p.U();
        RWF[i] = p.RWF();
        ERot[i] = p.ERot();
        vibLevel[i] = p.vibLevel();
        ELevel[i] = p.ELevel();
        typeId[i] = p.typeId();
        newParcel[i] = p.newParcel();
        classification[i] = p.classification();

        stuckToWall[i] = p.isStuck();
        if (stuckToWall[i])
        {
            wallTemperature[i] = p.stuck().wallTemperature();
            wallVectors[i] = p.stuck().wallVectors();
        }

        isTracked[i] = p.isTracked();
        if (isTracked[i])
        {
            inPatchId[i] = p.tracked().inPatchId();
            tracerInitialTime[i] = p.tracked().initialTime();
            tracerInitialPosition[i] = p.tracked().initialPosition();
            tracerCurrentPosition[i] = p.tracked().currentPosition();
            tracerDistanceTravelled[i] = p.tracked().distanceTravelledVector();
        }
        else
        {
            inPatchId[i] = -1;
            tracerInitialTime[i] = 0;
            tracerInitialPosition[i] = vector::zero;
            tracerCurrentPosition[i] = vector::zero;
            tracerDistanceTravelled[i] = vector::zero;
        }

        i++;
    }

    U.write();

    if (gMax(RWF) > 1.0)
    {
        //- this is an axi/spherically -symmetric simulation
        RWF.write();
    }

    if (gMax(ERot) > 0.0)
    {
        //- there is at least one molecule
        ERot.write();
    }

    if (gMax(ELevel) > 0)
    {
        //- the electronic mode is activated
        ELevel.write();
    }

    typeId.write();
    newParcel.write();
    classification.write();

    if (gSum(stuckToWall) > 0)
    {
        //- there is at least one stickingWallPatch with a particle stuck on it
        stuckToWall.write();
        wallTemperature.write();
        wallVectors.write();
    }

    if (gSum(isTracked) > 0)
    {
        //- there is at least one tracked parcel in the domain
        isTracked.write();
        inPatchId.write();
        tracerInitialTime.write();
        tracerInitialPosition.write();
        tracerCurrentPosition.write();
        tracerDistanceTravelled.write();
    }

    vibLevel.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const dsmcParcel& p
)
{
    /*dsmcParcel::TrackedParcel tP = dsmcParcel::TrackedParcel();
    if (p.isTracked())
    {
        tP = p.tracked();
    }*/ //TODO vincent 11/05/2018

    dsmcParcel::StuckParcel sP = dsmcParcel::StuckParcel();
    if (p.isStuck())
    {
        sP = p.stuck();
    }

    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.U()
            << token::SPACE << p.RWF()
            << token::SPACE << p.ERot()
            << token::SPACE << p.ELevel()
            << token::SPACE << p.typeId()
            << token::SPACE << p.newParcel()
            << token::SPACE << p.tracked() //tP; // TODO vincent  11/05/2018
            << token::SPACE << p.classification()
            << token::SPACE << sP
            << token::SPACE << p.vibLevel();
    }
    else
    {
        os  << static_cast<const particle&>(p);

        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            sizeof(p.U())
            + sizeof(p.RWF())
            + sizeof(p.ERot())
            + sizeof(p.ELevel())
            + sizeof(p.typeId())
            + sizeof(p.newParcel())
            + sizeof(p.classification())
        );

        os << p.tracked(); //tP; // TODO vincent  11/05/2018
        os << sP;
        os << p.vibLevel();
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::dsmcParcel&)"
    );

    return os;
}


// ************************************************************************* //
