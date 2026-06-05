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
#include <cstdint>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noTimeCounter, 0);

addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounter, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace
{

struct FastRng
{
    uint64_t s0;
    uint64_t s1;

    FastRng(uint64_t seed = 1)
    {
        uint64_t z = seed + 0x9e3779b97f4a7c15ULL;
        z = (z ^ (z >> 30))*0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27))*0x94d049bb133111ebULL;
        s0 = z ^ (z >> 31);

        z = seed + 0x9e3779b97f4a7c15ULL + 1;
        z = (z ^ (z >> 30))*0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27))*0x94d049bb133111ebULL;
        s1 = z ^ (z >> 31);
    }

    inline uint64_t next()
    {
        const uint64_t oldS0 = s0;
        uint64_t oldS1 = s1;
        oldS1 ^= oldS1 << 23;
        s1 = oldS1 ^ oldS0 ^ (oldS1 >> 18) ^ (oldS0 >> 5);
        s0 = oldS1;
        return s1 + oldS0;
    }

    inline scalar sample01()
    {
        return scalar((next() >> 11)*0x1.0p-53);
    }

    inline label position(const label n)
    {
        return label(next()%uint64_t(n));
    }
};

}


void noTimeCounter::readControlDictParams()
{
    const dictionary& controlDict = cloud_.mesh().time().controlDict();

    collisionFastRng_ =
        controlDict.lookupOrDefault<bool>("collisionFastRng", false);
}


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
    threadCharges_(),
    collisionFastRng_(false)
//     propsDict_(dict.subDict(typeName + "Properties"))
{
    readControlDictParams();
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

    #ifdef _OPENMP
    if (cloud_.openmpEnabled())
    {
        const label statsThreads = max(cloud_.ompNumThreads(), label(1));
        const label collisionChunk = max(cloud_.openmpCollisionChunk(), label(1));
        const bool useFlatOccupancy = cloud_.hasOccupancyOrderedParcels();
        const DynamicList<DynamicList<dsmcParcel*>>* cellOccupancyPtr =
            useFlatOccupancy ? nullptr : &cloud_.cellOccupancy();
        const polyMesh& mesh = cloud_.mesh();
        const label nCells = mesh.nCells();

        labelList& nCandidatesPerCell = cloud_.nCandidatesPerCell();
        if (nCandidatesPerCell.size() != nCells)
        {
            nCandidatesPerCell.setSize(nCells, 0);
        }

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

        List<FastRng> threadFastRng(statsThreads);
        for (label threadI = 0; threadI < statsThreads; ++threadI)
        {
            threadFastRng[threadI] = FastRng
            (
                uint64_t(Pstream::myProcNo())*1000000ULL
              + uint64_t(threadI)*10000ULL
              + uint64_t(mesh.time().timeIndex() + 1)
            );
        }

        labelList threadCandidateCounts(statsThreads, 0);
        labelList threadAcceptedCounts(statsThreads, 0);

        auto processCell =
        [&]
        (
            const label cellI,
            const label threadI
        )
        -> label
        {
            FastRng& fastRng = threadFastRng[threadI];

            auto randomIndex = [&](const label n) -> label
            {
                return fastRng.position(n);
            };

            auto random01 = [&]() -> scalar
            {
                return fastRng.sample01();
            };

            const DynamicList<dsmcParcel*>* cellParcelsPtr =
                useFlatOccupancy ? nullptr : &(*cellOccupancyPtr)[cellI];
            const label occStart =
                useFlatOccupancy ? cloud_.occupancyStart(cellI) : 0;
            const label nC =
                useFlatOccupancy
              ? cloud_.occupancyCount(cellI)
              : cellParcelsPtr->size();

            nCandidatesPerCell[cellI] = 0;

            if (nC <= 1)
            {
                return 0;
            }

            const scalar selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
              + 0.5*nC*(nC - 1)
               *cloud_.nParticles(cellI)
               *cloud_.sigmaTcRMax()[cellI]
               *cloud_.deltaTValue(cellI)
               /mesh.cellVolumes()[cellI];

            const label nCandidates(selectedPairs);
            cloud_.collisionSelectionRemainder()[cellI] =
                selectedPairs - nCandidates;
            nCandidatesPerCell[cellI] = nCandidates;

            if (nCandidates <= 0)
            {
                return 0;
            }

            DynamicList<label>& whichSubCell = threadWhichSubCell_[threadI];
            List<DynamicList<label>>& subCells = threadSubCells_[threadI];
            DynamicList<dsmcParcel*>& parcelPtrs = threadParcelPtrs_[threadI];
            DynamicList<vector>& velocities = threadVelocities_[threadI];
            DynamicList<label>& typeIds = threadTypeIds_[threadI];
            DynamicList<label>& charges = threadCharges_[threadI];

            whichSubCell.setSize(nC);
            parcelPtrs.setSize(nC);
            velocities.setSize(nC);
            typeIds.setSize(nC);
            charges.setSize(nC);

            label subCellCounts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            label subCellOffsets[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            const point& cC = mesh.cellCentres()[cellI];

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
            label collisions = 0;

            for (label c = 0; c < nCandidates; ++c)
            {
                const label candidateP = randomIndex(nC);
                label candidateQ = -1;

                const List<label>& subCellPs =
                    subCells[whichSubCell[candidateP]];
                const label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    do
                    {
                        candidateQ = subCellPs[randomIndex(nSC)];
                    } while (candidateP == candidateQ);
                }
                else
                {
                    do
                    {
                        candidateQ = randomIndex(nC);
                    } while (candidateP == candidateQ);
                }

                const label typeIdP = typeIds[candidateP];
                const label typeIdQ = typeIds[candidateQ];
                const label chargeP = charges[candidateP];
                const label chargeQ = charges[candidateQ];

                if (chargeP == -1 && chargeQ == -1)
                {
                    continue;
                }

                dsmcParcel& parcelP = *parcelPtrs[candidateP];
                dsmcParcel& parcelQ = *parcelPtrs[candidateQ];

                const scalar sigmaTcR =
                    cloud_.binaryCollision().sigmaTcR(parcelP, parcelQ);

                if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                {
                    cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                }

                if ((sigmaTcR/sigmaTcRMax) > random01())
                {
                    const label rMId =
                        cloud_.reactions().nReactions() > 0
                      ? cloud_.reactions().pairModelAddressing()[typeIdP][typeIdQ]
                      : -1;

                    if (rMId != -1)
                    {
                        cloud_.reactions().reactions()[rMId]->reaction
                        (
                            parcelP,
                            parcelQ
                        );

                        if (cloud_.reactions().reactions()[rMId]->relax())
                        {
                            cloud_.binaryCollision().collide
                            (
                                parcelP,
                                parcelQ,
                                cellI
                            );
                        }
                    }
                    else
                    {
                        cloud_.binaryCollision().collide
                        (
                            parcelP,
                            parcelQ,
                            cellI
                        );
                    }

                    ++collisions;

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

            return collisions;
        };

        #pragma omp parallel num_threads(statsThreads)
        {
            const label threadI = omp_get_thread_num();
            label localCandidates = 0;
            label localCollisions = 0;

            if (cloud_.openmpCollisionSchedule() == "static")
            {
                #pragma omp for schedule(static, collisionChunk)
                for (label cellI = 0; cellI < nCells; ++cellI)
                {
                    localCollisions += processCell(cellI, threadI);
                    localCandidates += nCandidatesPerCell[cellI];
                }
            }
            else if (cloud_.openmpCollisionSchedule() == "guided")
            {
                #pragma omp for schedule(guided, collisionChunk)
                for (label cellI = 0; cellI < nCells; ++cellI)
                {
                    localCollisions += processCell(cellI, threadI);
                    localCandidates += nCandidatesPerCell[cellI];
                }
            }
            else
            {
                #pragma omp for schedule(dynamic, collisionChunk)
                for (label cellI = 0; cellI < nCells; ++cellI)
                {
                    localCollisions += processCell(cellI, threadI);
                    localCandidates += nCandidatesPerCell[cellI];
                }
            }

            threadCandidateCounts[threadI] = localCandidates;
            threadAcceptedCounts[threadI] = localCollisions;
        }

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
                    << "    Collision candidates           = "
                    << collisionCandidates << nl
                    << "    Collision acceptance rate      = "
                    << scalar(collisions)/scalar(collisionCandidates) << nl
                    << endl;
            }
            else
            {
                Info<< "    No collisions" << endl;
            }

            infoCounter_ = 0;
        }

        return;
    }
    #endif

    // Temporary storage for subCells
    List<DynamicList<label>> subCells(8);

    label collisionCandidates = 0;

    label collisions = 0;

    const polyMesh& mesh = cloud_.mesh();
    const label nCells = mesh.nCells();
    const bool useFlatOccupancy = cloud_.hasOccupancyOrderedParcels();
    const DynamicList<DynamicList<dsmcParcel*>>* cellOccupancyPtr =
        useFlatOccupancy ? nullptr : &cloud_.cellOccupancy();

    labelList& nCandidatesPerCell = cloud_.nCandidatesPerCell();
    if (nCandidatesPerCell.size() != nCells)
    {
        nCandidatesPerCell.setSize(nCells, 0);
    }
    forAll(nCandidatesPerCell, cellI)
    {
        nCandidatesPerCell[cellI] = 0;
    }

    FastRng fastRng
    (
        uint64_t(Pstream::myProcNo())*1000000ULL
      + uint64_t(mesh.time().timeIndex() + 1)
    );

    auto randomLabel = [&](const label minValue, const label maxValue) -> label
    {
        if (!collisionFastRng_)
        {
            return cloud_.randomLabel(minValue, maxValue);
        }

        return minValue + fastRng.position(maxValue - minValue + 1);
    };

    auto random01 = [&]() -> scalar
    {
        return collisionFastRng_
             ? fastRng.sample01()
             : rndGen_.sample01<scalar>();
    };

    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        const scalar deltaT = cloud_.deltaTValue(cellI);

        const DynamicList<dsmcParcel*>* cellParcelsPtr =
            useFlatOccupancy ? nullptr : &(*cellOccupancyPtr)[cellI];
        const label occStart =
            useFlatOccupancy ? cloud_.occupancyStart(cellI) : 0;

        const scalar& cellVolume = mesh.cellVolumes()[cellI];

        const label nC =
            useFlatOccupancy
          ? cloud_.occupancyCount(cellI)
          : cellParcelsPtr->size();

        if (nC > 1)
        {

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(nC);

            const point& cC = mesh.cellCentres()[cellI];

            for (label i = 0; i < nC; ++i)
            {
                const dsmcParcel& p =
                    useFlatOccupancy
                  ? *cloud_.occupancyParcel(occStart + i)
                  : *(*cellParcelsPtr)[i];

                vector relPos = p.position() - cC;

                label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                subCells[subCell].append(i);

                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];

            //scalar selectedPairs = 0.0;

            scalar selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*cloud_.nParticles(cellI)*sigmaTcRMax*deltaT
                /cellVolume;

            const label nCandidates(selectedPairs);

            cloud_.collisionSelectionRemainder()[cellI] = selectedPairs - nCandidates;

            nCandidatesPerCell[cellI] = nCandidates;

            collisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                //label candidateP = rndGen_.position<label>(0, nC - 1);
                label candidateP = randomLabel(0, nC-1);

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
                        candidateQ = subCellPs[randomLabel(0, nSC-1)];

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
                        candidateQ = randomLabel(0, nC-1);

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

                dsmcParcel& parcelP =
                    useFlatOccupancy
                  ? *cloud_.occupancyParcel(occStart + candidateP)
                  : *(*cellParcelsPtr)[candidateP];
                dsmcParcel& parcelQ =
                    useFlatOccupancy
                  ? *cloud_.occupancyParcel(occStart + candidateQ)
                  : *(*cellParcelsPtr)[candidateQ];

                label chargeP = -2;
                label chargeQ = -2;

                chargeP = cloud_.constProps(parcelP.typeId()).charge();
                chargeQ = cloud_.constProps(parcelQ.typeId()).charge();

                //do not allow electron-electron collisions

                if(!(chargeP == -1 && chargeQ == -1))
                {

                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP,
                        parcelQ
                    );


                    // Update the maximum value of sigmaTcR stored, but use the
                    // initial value in the acceptance-rejection criteria because
                    // the number of collision candidates selected was based on this


                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    if ((sigmaTcR/sigmaTcRMax) > random01())
                    {
                        // chemical reactions

                        // find which reaction model parcel p and q should use
                        label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);

    //                             Info << " parcelP id: " <<  parcelP.typeId()
    //                                 << " parcelQ id: " << parcelQ.typeId()
    //                                 << " reaction model: " << rMId
    //                                 << endl;

                        if(rMId != -1)
                        {
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

                        collisions++;
                    }
                }
            }
        }
    }

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
