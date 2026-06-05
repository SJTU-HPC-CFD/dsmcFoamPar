/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------*/

#ifndef dsmcLocalMesh_H
#define dsmcLocalMesh_H

#include "fvMesh.H"
#include "labelList.H"

namespace Foam
{

class dsmcLocalMesh
{
    const fvMesh* globalMesh_;
    labelList globalToLocalCell_;       // global → local cell index (-1 = absent)
    labelList localToGlobalCell_;       // local → global cell index
    List<bool> isHaloCell_;             // local cell is halo (not owned)

public:
    dsmcLocalMesh();
    explicit dsmcLocalMesh(const fvMesh& globalMesh);
    ~dsmcLocalMesh();

    void build(const labelList& myCells);

    inline bool valid() const;
    inline label toLocal(label globalCell) const;
    inline label toGlobal(label localCell) const;
    inline bool isHalo(label localCell) const;
    inline label nCells() const;

    void report() const;
};

} // End namespace Foam

#include "dsmcLocalMeshI.H"

#endif
