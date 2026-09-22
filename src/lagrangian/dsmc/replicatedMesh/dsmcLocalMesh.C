/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------*/

#include "dsmcLocalMesh.H"

namespace Foam
{

dsmcLocalMesh::dsmcLocalMesh()
:
    globalMesh_(nullptr),
    globalToLocalCell_(),
    localToGlobalCell_(),
    isHaloCell_()
{}


dsmcLocalMesh::dsmcLocalMesh(const fvMesh& globalMesh)
:
    globalMesh_(&globalMesh),
    globalToLocalCell_(),
    localToGlobalCell_(),
    isHaloCell_()
{}


dsmcLocalMesh::~dsmcLocalMesh()
{}


void dsmcLocalMesh::build(const labelList& myCells)
{
    const label nGlobalCells = globalMesh_->nCells();

    // ---- collect myCells + 1-ring halo --------------------------------
    List<bool> includeCell(nGlobalCells, false);
    forAll(myCells, i) includeCell[myCells[i]] = true;

    forAll(myCells, i)
    {
        const label cellI = myCells[i];
        const cell& c = globalMesh_->cells()[cellI];
        forAll(c, j)
        {
            const label faceI = c[j];
            const label own = globalMesh_->faceOwner()[faceI];
            const label nei =
                globalMesh_->isInternalFace(faceI)
              ? globalMesh_->faceNeighbour()[faceI]
              : -1;
            const label other = (own == cellI) ? nei : own;
            if (other >= 0) includeCell[other] = true;
        }
    }

    // ---- build local→global mapping -----------------------------------
    label nTotal = 0;
    forAll(includeCell, cellI)
    {
        if (includeCell[cellI])
        {
            ++nTotal;
        }
    }
    localToGlobalCell_.setSize(nTotal, -1);
    globalToLocalCell_.setSize(nGlobalCells, -1);

    label localIdx = 0;
    for (label globalI = 0; globalI < nGlobalCells; ++globalI)
    {
        if (includeCell[globalI])
        {
            localToGlobalCell_[localIdx] = globalI;
            globalToLocalCell_[globalI] = localIdx;
            ++localIdx;
        }
    }

    // ---- identify halo cells -------------------------------------------
    // Ownership membership test by bitmap.  The previous per-cell
    // findIndex over myCells is O(nTotal*nOwned) and dominated DLB
    // repartition wall time on ~13M-cell meshes.
    List<bool> isOwnedCell(nGlobalCells, false);
    forAll(myCells, i) isOwnedCell[myCells[i]] = true;

    isHaloCell_.setSize(nTotal, false);
    for (label localI = 0; localI < nTotal; ++localI)
    {
        isHaloCell_[localI] = !isOwnedCell[localToGlobalCell_[localI]];
    }

    const label nOwned = myCells.size();
    const label nHalo = nTotal - nOwned;

    Info<< "Local mesh: " << nOwned << " owned + " << nHalo
        << " halo = " << nTotal << " / " << nGlobalCells << " cells"
        << endl;
}


void dsmcLocalMesh::report() const
{
    if (localToGlobalCell_.empty()) return;
    Info<< "Local mesh: " << localToGlobalCell_.size() << " cells "
        << "(of " << globalMesh_->nCells() << " global)"
        << endl;
}

} // End namespace Foam
