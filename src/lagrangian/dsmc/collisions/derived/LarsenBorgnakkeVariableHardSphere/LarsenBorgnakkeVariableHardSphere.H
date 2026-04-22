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
    Foam::LarsenBorgnakkeVariableHardSphere

Description
    Variable Hard Sphere BinaryCollision Model with Larsen Borgnakke internal
    energy redistribution.

\*---------------------------------------------------------------------------*/

#ifndef LarsenBorgnakkeVariableHardSphere_H
#define LarsenBorgnakkeVariableHardSphere_H

#include "VariableHardSphere.H"

namespace Foam
{

class LarsenBorgnakkeVariableHardSphere
:
    public VariableHardSphere
{
    const dictionary coeffDictLB_;
    const scalar rotationalRelaxationCollisionNumber_;
    const scalar vibrationalRelaxationCollisionNumber_;
    label invZvFormulation_;
    const scalar electronicRelaxationCollisionNumber_;
    mutable scalar detailRedistributePWallTime_;
    mutable scalar detailRedistributeQWallTime_;
    mutable scalar detailScatterWallTime_;
    mutable label detailCollideCalls_;

public:

    TypeName("LarsenBorgnakkeVariableHardSphere");

    LarsenBorgnakkeVariableHardSphere
    (
        const dictionary& dict,
        dsmcCloud& cloud
    );

    ~LarsenBorgnakkeVariableHardSphere();

    inline const scalar& vibrationalRelaxationNumber() const
    {
        return vibrationalRelaxationCollisionNumber_;
    }

    virtual void collide
    (
        dsmcParcel& pP,
        dsmcParcel& pQ,
        const label cellI,
        scalar cR = -1
    );

    virtual void redistribute
    (
        dsmcParcel& p,
        scalar& translationalEnergy,
        const scalar omegaPQ,
        const bool postReaction = false
    );

    virtual const dictionary& coeffDict() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
