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

#include "BinaryCollisionModel.H"

namespace Foam
{
    defineTypeNameAndDebug(BinaryCollisionModel, 0);
    defineRunTimeSelectionTable(BinaryCollisionModel, dictionary);
}

Foam::BinaryCollisionModel::BinaryCollisionModel(dsmcCloud& owner)
:
    dict_(dictionary::null),
    cloud_(owner)
{}

Foam::BinaryCollisionModel::BinaryCollisionModel
(
    const dictionary& dict,
    dsmcCloud& owner
)
:
    dict_(dict),
    cloud_(owner)
{}

Foam::autoPtr<Foam::BinaryCollisionModel> Foam::BinaryCollisionModel::New
(
    const dictionary& dict,
    dsmcCloud& owner
)
{
    const word modelType(dict.get<word>("BinaryCollisionModel"));

    Info<< "Selecting BinaryCollisionModel " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "BinaryCollisionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<BinaryCollisionModel>(ctorPtr(dict, owner));
}

Foam::BinaryCollisionModel::~BinaryCollisionModel()
{}

const Foam::dictionary& Foam::BinaryCollisionModel::dict() const
{
    return dict_;
}

// ************************************************************************* //
