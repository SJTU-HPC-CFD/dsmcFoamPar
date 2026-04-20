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
    noTimeCounter

Description

\*----------------------------------------------------------------------------*/

#include "noTimeCounter.H"
#include "addToRunTimeSelectionTable.H"
#include <chrono>

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounter, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
noTimeCounter::noTimeCounter
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    infoCounter_(0),
    threadWhichSubCell_(),
    threadSubCells_(),
    threadParcelPtrs_(),
    threadVelocities_(),
    threadTypeIds_(),
    threadCharges_()
//     propsDict_(dict.subDict(typeName + "Properties"))
{
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noTimeCounter::~noTimeCounter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void noTimeCounter::initialConfiguration()
{

}

void noTimeCounter::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    const label statsThreads =
        cloud_.openmpEnabled() ? max(cloud_.ompNumThreads(), label(1)) : label(1);
    labelList threadCandidateCounts(statsThreads, 0);
    labelList threadAcceptedCounts(statsThreads, 0);
    labelList threadActiveCellCounts(statsThreads, 0);
    labelList threadReactionHitCounts(statsThreads, 0);
    scalarField threadWallTimes(statsThreads, 0.0);

    const bool useFlatOccupancy = cloud_.hasOccupancyOrderedParcels();
    const List<DynamicList<dsmcParcel*>>* cellOccupancyPtr =
        useFlatOccupancy ? nullptr : &cloud_.cellOccupancy();

    const polyMesh& mesh = cloud_.mesh();
    const label nCells = mesh.nCells();
    const label collisionChunk = max(cloud_.openmpCollisionChunk(), label(1));

    if (threadWhichSubCell_.size() != statsThreads)
    {
        threadWhichSubCell_.setSize(statsThreads);
        threadSubCells_.setSize(statsThreads);
        threadParcelPtrs_.setSize(statsThreads);
        threadVelocities_.setSize(statsThreads);
        threadTypeIds_.setSize(statsThreads);
        threadCharges_.setSize(statsThreads);
    }

    for (label threadI = 0; threadI < statsThreads; ++threadI)
    {
        if (threadSubCells_[threadI].size() != 8)
        {
            threadSubCells_[threadI].setSize(8);
        }
    }

    auto processCell =
    [&]
    (
        const label cellI,
        const label threadI,
        label& activeCellCount,
        label& reactionHitCount
    )
    -> label
    {
        const DynamicList<dsmcParcel*>* cellParcelsPtr =
            useFlatOccupancy ? nullptr : &(*cellOccupancyPtr)[cellI];
        const label occStart = useFlatOccupancy ? cloud_.occupancyStart(cellI) : 0;
        const label nC =
            useFlatOccupancy ? cloud_.occupancyCount(cellI) : cellParcelsPtr->size();
        const label nCandidates = cloud_.nCandidatesPerCell()[cellI];
        label acceptedCollisions = 0;

        if (nC > 1 && nCandidates > 0)
        {
            ++activeCellCount;
            DynamicList<label>& whichSubCell = threadWhichSubCell_[threadI];
            List<DynamicList<label>>& subCells = threadSubCells_[threadI];
            DynamicList<dsmcParcel*>& parcelPtrs = threadParcelPtrs_[threadI];
            DynamicList<vector>& velocities = threadVelocities_[threadI];
            DynamicList<label>& typeIds = threadTypeIds_[threadI];
            DynamicList<label>& charges = threadCharges_[threadI];

            const point& cC = mesh.cellCentres()[cellI];
            label subCellCounts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            label subCellOffsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};

            whichSubCell.setSize(nC);
            parcelPtrs.setSize(nC);
            velocities.setSize(nC);
            typeIds.setSize(nC);
            charges.setSize(nC);

            for (label i = 0; i < nC; ++i)
            {
                dsmcParcel* pPtr =
                    useFlatOccupancy
                  ? cloud_.occupancyParcel(occStart + i)
                  : (*cellParcelsPtr)[i];
                const label typeId = pPtr->typeId();
                const vector relPos = pPtr->position() - cC;
                const label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                parcelPtrs[i] = pPtr;
                velocities[i] = pPtr->U();
                typeIds[i] = typeId;
                charges[i] = cloud_.constProps(typeId).charge();
                whichSubCell[i] = subCell;
                ++subCellCounts[subCell];
            }

            for (label subCellI = 0; subCellI < 8; ++subCellI)
            {
                subCells[subCellI].setSize(subCellCounts[subCellI]);
            }

            for (label i = 0; i < nC; ++i)
            {
                const label subCell = whichSubCell[i];
                subCells[subCell][subCellOffsets[subCell]++] = i;
            }

            scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];

            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                //label candidateP = rndGen_.position<label>(0, nC - 1);
                label candidateP = cloud_.randomLabel(0, nC-1);

                // Declare the second collision candidate
                label candidateQ = -1;

                const List<label>& subCellPs = subCells[whichSubCell[candidateP]];

                const label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        //candidateQ = subCellPs[rndGen_.position<label>(0, nSC - 1)]; OLD
                        candidateQ = subCellPs[cloud_.randomLabel(0, nSC-1)];

                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        //candidateQ = rndGen_.position<label>(0, nC - 1); OLD
                        candidateQ = cloud_.randomLabel(0, nC-1);

                    } while (candidateP == candidateQ);
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = cloud_.randomLabel(0, nC-1);

                // // Select a possible second collision candidate
                // label candidateQ = cloud_.randomLabel(0, nC-1);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = cloud_.randomLabel(0, nC-1);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                const label typeIdP = typeIds[candidateP];
                const label typeIdQ = typeIds[candidateQ];
                label chargeP = charges[candidateP];
                label chargeQ = charges[candidateQ];

                if (chargeP == -2)
                {
                    chargeP = cloud_.constProps(typeIdP).charge();
                }

                if (chargeQ == -2)
                {
                    chargeQ = cloud_.constProps(typeIdQ).charge();
                }

                //do not allow electron-electron collisions

                if(!(chargeP == -1 && chargeQ == -1))
                {

                    const scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        *parcelPtrs[candidateP],
                        *parcelPtrs[candidateQ]
                    );


                    // Update the maximum value of sigmaTcR stored, but use the
                    // initial value in the acceptance-rejection criteria because
                    // the number of collision candidates selected was based on this


                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    if ((sigmaTcR/sigmaTcRMax) > cloud_.rndGen().sample01<scalar>())
                    {
                        // chemical reactions

                        // find which reaction model parcel p and q should use
                        const label rMId =
                            cloud_.reactionsActive()
                          ? cloud_.reactions().pairModelAddressing()[typeIdP][typeIdQ]
                          : -1;

                        dsmcParcel& parcelP = *parcelPtrs[candidateP];
                        dsmcParcel& parcelQ = *parcelPtrs[candidateQ];

    //                             Info << " parcelP id: " <<  parcelP.typeId()
    //                                 << " parcelQ id: " << parcelQ.typeId()
    //                                 << " reaction model: " << rMId
    //                                 << endl;

                        if(rMId != -1)
                        {
                            ++reactionHitCount;
                            // try to react molecules
    //                         if(cloud_.reactions().reactions()[rMId]->reactWithLists())
    //                         {
                                // so far for recombination only
    //                                     reactions_.reactions()[rMId]->reaction
    //                                     (
    //                                         parcelP,
    //                                         parcelQ,
    //                                         candidateList,
    //                                         candidateSubList,
    //                                         candidateP,
    //                                         whichSubCell
    //                                     );
    //                         }
    //                         else
    //                         {
                                cloud_.reactions().reactions()[rMId]->reaction
                                (
                                    parcelP,
                                    parcelQ
                                );
    //                         }
                            // if reaction unsuccessful use conventional collision model
                            if(cloud_.reactions().reactions()[rMId]->relax())
                            {
                                cloud_.binaryCollision().collide
                                (
                                    parcelP,
                                    parcelQ,
                                    cellI
                                );
                            }
                        }
                        else // if reaction model not found, use conventional collision model
                        {
                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }

                        acceptedCollisions++;

                        velocities[candidateP] = parcelP.U();
                        velocities[candidateQ] = parcelQ.U();
                        typeIds[candidateP] = parcelP.typeId();
                        typeIds[candidateQ] = parcelQ.typeId();
                        charges[candidateP] =
                            cloud_.constProps(typeIds[candidateP]).charge();
                        charges[candidateQ] =
                            cloud_.constProps(typeIds[candidateQ]).charge();
                    }
                }
            }
        }
        return acceptedCollisions;
    };

    #ifdef _OPENMP
    if (cloud_.openmpEnabled())
    {
        if (cloud_.openmpCollisionStrategy() == "partition")
        {
            #pragma omp parallel
            {
                using clock_type = std::chrono::steady_clock;
                const auto tBegin = clock_type::now();
                const label threadI = cloud_.currentThreadId();
                const label startCell = cloud_.collisionLoadStart()[threadI];
                const label endCell = cloud_.collisionLoadEnd()[threadI];
                label localCandidates = 0;
                label localCollisions = 0;
                label localActiveCells = 0;
                label localReactionHits = 0;

                for (label cellI = startCell; cellI < endCell; ++cellI)
                {
                    localCandidates += cloud_.nCandidatesPerCell()[cellI];
                    localCollisions += processCell
                    (
                        cellI,
                        threadI,
                        localActiveCells,
                        localReactionHits
                    );
                }

                threadCandidateCounts[threadI] = localCandidates;
                threadAcceptedCounts[threadI] = localCollisions;
                threadActiveCellCounts[threadI] = localActiveCells;
                threadReactionHitCounts[threadI] = localReactionHits;
                threadWallTimes[threadI] =
                    std::chrono::duration<scalar>(clock_type::now() - tBegin).count();
            }
        }
        else
        {
            #pragma omp parallel
            {
                using clock_type = std::chrono::steady_clock;
                const auto tBegin = clock_type::now();
                const label threadI = cloud_.currentThreadId();
                label localCandidates = 0;
                label localCollisions = 0;
                label localActiveCells = 0;
                label localReactionHits = 0;

                #pragma omp for schedule(dynamic, collisionChunk)
                for (label cellI = 0; cellI < nCells; ++cellI)
                {
                    localCandidates += cloud_.nCandidatesPerCell()[cellI];
                    localCollisions += processCell
                    (
                        cellI,
                        threadI,
                        localActiveCells,
                        localReactionHits
                    );
                }

                threadCandidateCounts[threadI] = localCandidates;
                threadAcceptedCounts[threadI] = localCollisions;
                threadActiveCellCounts[threadI] = localActiveCells;
                threadReactionHitCounts[threadI] = localReactionHits;
                threadWallTimes[threadI] =
                    std::chrono::duration<scalar>(clock_type::now() - tBegin).count();
            }
        }
    }
    else
    #endif
    {
        using clock_type = std::chrono::steady_clock;
        const auto tBegin = clock_type::now();
        label localActiveCells = 0;
        label localReactionHits = 0;

        for (label cellI = 0; cellI < nCells; ++cellI)
        {
            threadCandidateCounts[0] += cloud_.nCandidatesPerCell()[cellI];
            threadAcceptedCounts[0] += processCell
            (
                cellI,
                0,
                localActiveCells,
                localReactionHits
            );
        }

        threadActiveCellCounts[0] = localActiveCells;
        threadReactionHitCounts[0] = localReactionHits;
        threadWallTimes[0] =
            std::chrono::duration<scalar>(clock_type::now() - tBegin).count();
    }

    cloud_.recordCollisionThreadProfile
    (
        threadCandidateCounts,
        threadAcceptedCounts,
        threadActiveCellCounts,
        threadReactionHitCounts,
        threadWallTimes
    );

    label collisionCandidates = sum(threadCandidateCounts);
    label collisions = sum(threadAcceptedCounts);

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    cloud_.sigmaTcRMax().correctBoundaryConditions();

    infoCounter_++;

    if(infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisionCandidates)
        {
            Info<< "    Collisions                      = "
                << collisions << nl
    //             << "    Acceptance rate                 = "
    //             << scalar(collisions)/scalar(collisionCandidates) << nl
                << endl;

            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No collisions" << endl;

            infoCounter_ = 0;
        }
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
