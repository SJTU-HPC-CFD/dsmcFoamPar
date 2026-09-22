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
    Foam::IOPosition

Description
    Helper IO class to read and write particle positions

SourceFiles
    IOPosition.C

\*---------------------------------------------------------------------------*/

#ifndef IOPosition_H
#define IOPosition_H

#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IOPosition Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IOPosition
:
    public regIOobject
{

    // Private data

        //- Reference to the cloud
        const CloudType& cloud_;


public:

    // Static data

        //- Runtime type name information. Use cloud type.
        virtual const word& type() const
        {
            return Cloud<typename CloudType::particleType>::typeName;
            //cloud_.type();
        }


    // Constructors

        //- Construct from cloud
        IOPosition(const CloudType&);


    // Member functions

        //- Read the positions header, return the particle count and set
        //  entriesStart to the byte offset of the first entry
        label readHeader(off_t& entriesStart);

        //- Read positions without filtering (all particles)
        virtual void readData(CloudType& c, bool checkClass);

        //- Read positions keeping only parcels whose cell is owned by
        //  filterRank according to cellOwner.  Appends the kept file
        //  indices to keep and sets nFull to the total entry count.
        //  Non-owned parcels are constructed and immediately deleted so
        //  that the cell/tet assignment semantics are identical to the
        //  unfiltered read.
        void readDataFiltered
        (
            CloudType& c,
            bool checkClass,
            const labelUList& cellOwner,
            const label filterRank,
            labelList& keep,
            label& nFull
        );

        virtual bool write() const;

        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOPosition.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
