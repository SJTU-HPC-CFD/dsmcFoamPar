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


Description

    no-time counter with reactions

SourceFiles
    noTimeCounterI.H
    noTimeCounter.C
    noTimeCounterIO.C

\*---------------------------------------------------------------------------*/

#ifndef noTimeCounter_H
#define noTimeCounter_H


#include "collisionPartnerSelection.H"
#include "dsmcCloud.H"
#include "../../../dsmcFastRng.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noTimeCounter Declaration
\*---------------------------------------------------------------------------*/

class noTimeCounter
:
    public collisionPartnerSelection
{

private:

    // Private data

//         dictionary propsDict_;
    label infoCounter_;
    List<DynamicList<label>> threadWhichSubCell_;
    List<List<DynamicList<label>>> threadSubCells_;
    List<DynamicList<dsmcParcel*>> threadParcelPtrs_;
    List<DynamicList<label>> threadTypeIds_;
    List<DynamicList<label>> threadCharges_;
    List<DynamicList<label>> threadCandidateCells_;
    List<dsmcPaddedFastRng*> collisionThreadFastRng_;
    DynamicList<label> candidateCellsToClear_;
    label profileLocalCollisions_;
    label profileLocalCollisionCandidates_;
    scalar profileCollisionLocalLoopWall_;
    scalar profileCollisionReduceWall_;
    scalar profileCollisionSigmaWall_;
    scalar profileCollisionPreOutputWall_;

    bool fastRng_;
    bool collisionReduceOnlyOnOutput_;
    bool collisionOutputGlobalReduce_;

    void readControlDictParams();

    void ensureCollisionThreadFastRngs(const label nThreads);

    void clearCandidateCounts(labelList& nCandidatesPerCell, const label nCells);

    void printReplicatedRankCollisionDetail
    (
        const label localCollisions,
        const label localCollisionCandidates,
        const label globalCollisions,
        const label globalCollisionCandidates
    ) const;




public:

    //- Runtime type information
    TypeName("noTimeCounter");

    // Static data members

    // Constructors

        //- construct from components
        noTimeCounter
        (
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~noTimeCounter();

    // Member Functions

        // Access

        // Check

        // Edit

            void initialConfiguration();

            void collide();

            bool hasCollisionSubphaseProfile() const override
            {
                return true;
            }

            scalar collisionLocalLoopWallTime() const override
            {
                return profileCollisionLocalLoopWall_;
            }

            scalar collisionReduceWallTime() const override
            {
                return profileCollisionReduceWall_;
            }

            scalar collisionSigmaWallTime() const override
            {
                return profileCollisionSigmaWall_;
            }

            scalar collisionTotalWallTime() const override
            {
                return profileCollisionPreOutputWall_;
            }

            label collisionLocalAcceptedCount() const override
            {
                return profileLocalCollisions_;
            }

            label collisionLocalCandidateCount() const override
            {
                return profileLocalCollisionCandidates_;
            }

        // Write





    // Member Operators

//         void operator=(const noTimeCounter&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const noTimeCounter&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
