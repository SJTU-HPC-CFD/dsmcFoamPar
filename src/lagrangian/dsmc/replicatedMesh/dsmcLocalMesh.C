/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------*/

#include "dsmcLocalMesh.H"
#include "bitSet.H"

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
    bitSet includeCell(nGlobalCells, false);
    forAll(myCells, i) includeCell.set(myCells[i]);

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
            if (other >= 0) includeCell.set(other);
        }
    }

    // ---- build local→global mapping -----------------------------------
    const label nTotal = includeCell.count();
    localToGlobalCell_.setSize(nTotal, -1);
    globalToLocalCell_.setSize(nGlobalCells, -1);

    label localIdx = 0;
    for (label globalI = 0; globalI < nGlobalCells; ++globalI)
    {
        if (includeCell.test(globalI))
        {
            localToGlobalCell_[localIdx] = globalI;
            globalToLocalCell_[globalI] = localIdx;
            ++localIdx;
        }
    }

    // ---- identify halo cells -------------------------------------------
    isHaloCell_.reset();
    isHaloCell_.setSize(nTotal, false);
    for (label localI = 0; localI < nTotal; ++localI)
    {
        const label globalI = localToGlobalCell_[localI];
        if (myCells.found(globalI) < 0)
        {
            isHaloCell_.set(localI);
        }
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
