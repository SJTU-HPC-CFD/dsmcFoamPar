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

#include "noTimeCounterOMP.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(noTimeCounterOMP, 0);
addToRunTimeSelectionTable(collisionPartnerSelection, noTimeCounterOMP, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noTimeCounterOMP::noTimeCounterOMP
(
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    collisionPartnerSelection(mesh, cloud, dict),
    infoCounter_(0)
{}

noTimeCounterOMP::~noTimeCounterOMP()
{}

void noTimeCounterOMP::initialConfiguration()
{}

void noTimeCounterOMP::collide()
{
    if (!cloud_.binaryCollision().active())
    {
        return;
    }

    // --- OpenMP parallel collision (zero-copy backend) ---
    // Reactions are handled inside the parallel loop with #pragma omp critical
    // for addNewParcel (rare events, acceptable overhead).
    collideParallel();
}

void noTimeCounterOMP::collideSerial()
{
    // Serial fallback: identical to noTimeCounter::collide()
    List<DynamicList<label>> subCells(8);
    label collisionCandidates = 0;
    label collisions = 0;

    const List<DynamicList<dsmcParcel*>>& cellOccupancy = cloud_.cellOccupancy();
    const polyMesh& mesh = cloud_.mesh();

    forAll(cellOccupancy, cellI)
    {
        const scalar deltaT = cloud_.deltaTValue(cellI);
        const DynamicList<dsmcParcel*>& cellParcels(cellOccupancy[cellI]);
        const scalar& cellVolume = mesh.cellVolumes()[cellI];
        const label nC(cellParcels.size());

        if (nC > 1)
        {
            forAll(subCells, i) { subCells[i].clear(); }
            List<label> whichSubCell(nC);
            const point& cC = mesh.cellCentres()[cellI];

            forAll(cellParcels, i)
            {
                const dsmcParcel& p = *cellParcels[i];
                vector relPos = p.position() - cC;
                label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());
                subCells[subCell].append(i);
                whichSubCell[i] = subCell;
            }

            scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];
            scalar selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*cloud_.nParticles(cellI)*sigmaTcRMax*deltaT
                /cellVolume;
            const label nCandidates(selectedPairs);
            cloud_.collisionSelectionRemainder()[cellI] = selectedPairs - nCandidates;
            collisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                label candidateP = cloud_.randomLabel(0, nC-1);
                label candidateQ = -1;
                const List<label>& subCellPs = subCells[whichSubCell[candidateP]];
                const label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    do
                    {
                        candidateQ = subCellPs[cloud_.randomLabel(0, nSC-1)];
                    } while (candidateP == candidateQ);
                }
                else
                {
                    do
                    {
                        candidateQ = cloud_.randomLabel(0, nC-1);
                    } while (candidateP == candidateQ);
                }

                dsmcParcel& parcelP = *cellParcels[candidateP];
                dsmcParcel& parcelQ = *cellParcels[candidateQ];

                label chargeP = cloud_.constProps(parcelP.typeId()).charge();
                label chargeQ = cloud_.constProps(parcelQ.typeId()).charge();

                if (!(chargeP == -1 && chargeQ == -1))
                {
                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP, parcelQ
                    );

                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    if ((sigmaTcR/sigmaTcRMax) > rndGen_.sample01<scalar>())
                    {
                        label rMId = cloud_.reactions().returnModelId(parcelP, parcelQ);
                        if (rMId != -1)
                        {
                            cloud_.reactions().reactions()[rMId]->reaction(parcelP, parcelQ);
                            if (cloud_.reactions().reactions()[rMId]->relax())
                            {
                                cloud_.binaryCollision().collide(parcelP, parcelQ, cellI);
                            }
                        }
                        else
                        {
                            cloud_.binaryCollision().collide(parcelP, parcelQ, cellI);
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
    if (infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisionCandidates)
        {
            Info<< "    Collisions                      = "
                << collisions << nl << endl;
        }
        else
        {
            Info<< "    No collisions" << endl;
        }
        infoCounter_ = 0;
    }
}


void noTimeCounterOMP::collideParallel()
{
    // OpenMP parallel collision using cellFirst_/cellNext_ index
    const label nCells = mesh_.nCells();
    const polyMesh& mesh = cloud_.mesh();
    const DynamicList<dsmcParcel*>& parcelPtrs = cloud_.parcelPtrs();
    const labelList& cellFirst = cloud_.cellFirst();
    const labelList& cellNext = cloud_.cellNext();
    const labelList& cellCount = cloud_.cellCount();

    label totalCollisions = 0;
    label totalCandidates = 0;

    #pragma omp parallel reduction(+:totalCollisions, totalCandidates)
    {
        #ifdef _OPENMP
            const int tid = omp_get_thread_num();
        #else
            const int tid = 0;
        #endif
        Random& rng = cloud_.rng(tid);

        // Thread-local scratch for subcells
        List<DynamicList<label>> subCells(8);

        #pragma omp for schedule(dynamic, 64)
        for (label cellI = 0; cellI < nCells; cellI++)
        {
            const label nC = cellCount[cellI];
            if (nC < 2) continue;

            const scalar deltaT = cloud_.deltaTValue(cellI);
            const scalar& cellVolume = mesh.cellVolumes()[cellI];

            // Build particle list for this cell
            DynamicList<label> plist(nC);
            for (label ip = cellFirst[cellI]; ip >= 0; ip = cellNext[ip])
            {
                plist.append(ip);
            }

            // Assign to 8 subcells
            forAll(subCells, i) { subCells[i].clear(); }
            List<label> whichSubCell(nC);
            const point& cC = mesh.cellCentres()[cellI];

            forAll(plist, i)
            {
                const dsmcParcel& p = *parcelPtrs[plist[i]];
                vector relPos = p.position() - cC;
                label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());
                subCells[subCell].append(i);
                whichSubCell[i] = subCell;
            }

            // NTC collision selection
            scalar sigmaTcRMax = cloud_.sigmaTcRMax()[cellI];
            scalar selectedPairs =
                cloud_.collisionSelectionRemainder()[cellI]
                + 0.5*nC*(nC - 1)*cloud_.nParticles(cellI)*sigmaTcRMax*deltaT
                /cellVolume;
            const label nCandidates(selectedPairs);
            cloud_.collisionSelectionRemainder()[cellI] = selectedPairs - nCandidates;
            totalCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                // Select first candidate
                label candidateP = label(rng.sample01<scalar>() * nC);
                if (candidateP >= nC) candidateP = nC - 1;

                label candidateQ = -1;
                const List<label>& subCellPs = subCells[whichSubCell[candidateP]];
                const label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    do
                    {
                        candidateQ = subCellPs
                        [
                            label(rng.sample01<scalar>() * nSC)
                            % nSC
                        ];
                    } while (candidateP == candidateQ);
                }
                else
                {
                    do
                    {
                        candidateQ = label(rng.sample01<scalar>() * nC);
                        if (candidateQ >= nC) candidateQ = nC - 1;
                    } while (candidateP == candidateQ);
                }

                dsmcParcel& parcelP = *parcelPtrs[plist[candidateP]];
                dsmcParcel& parcelQ = *parcelPtrs[plist[candidateQ]];

                label chargeP = cloud_.constProps(parcelP.typeId()).charge();
                label chargeQ = cloud_.constProps(parcelQ.typeId()).charge();

                if (!(chargeP == -1 && chargeQ == -1))
                {
                    scalar sigmaTcR = cloud_.binaryCollision().sigmaTcR
                    (
                        parcelP, parcelQ
                    );

                    if (sigmaTcR > cloud_.sigmaTcRMax()[cellI])
                    {
                        cloud_.sigmaTcRMax()[cellI] = sigmaTcR;
                    }

                    if ((sigmaTcR/sigmaTcRMax) > rng.sample01<scalar>())
                    {
                        // Check for chemical reactions
                        if (cloud_.reactions().nReactions() > 0)
                        {
                            label rMId = cloud_.reactions().returnModelId
                            (
                                parcelP, parcelQ
                            );

                            if (rMId != -1)
                            {
                                // Reaction may call addNewParcel — protect
                                #pragma omp critical(reactionCritical)
                                {
                                    cloud_.reactions().reactions()[rMId]
                                        ->reaction(parcelP, parcelQ);
                                }

                                if (cloud_.reactions().reactions()[rMId]
                                    ->relax())
                                {
                                    cloud_.binaryCollision().collide
                                    (
                                        parcelP, parcelQ, cellI
                                    );
                                }
                            }
                            else
                            {
                                cloud_.binaryCollision().collide
                                (
                                    parcelP, parcelQ, cellI
                                );
                            }
                        }
                        else
                        {
                            cloud_.binaryCollision().collide
                            (
                                parcelP, parcelQ, cellI
                            );
                        }
                        totalCollisions++;
                    }
                }
            }
        }
    } // end omp parallel

    reduce(totalCollisions, sumOp<label>());
    reduce(totalCandidates, sumOp<label>());
    cloud_.sigmaTcRMax().correctBoundaryConditions();

    infoCounter_++;
    if (infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (totalCandidates)
        {
            Info<< "    Collisions (OMP)                = "
                << totalCollisions << nl << endl;
        }
        else
        {
            Info<< "    No collisions" << endl;
        }
        infoCounter_ = 0;
    }
}

} // End namespace Foam
