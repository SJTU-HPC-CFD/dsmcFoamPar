namespace Foam
{

inline bool dsmcLocalMesh::valid() const
{
    return globalToLocalCell_.size() > 0;
}

inline label dsmcLocalMesh::toLocal(label globalCell) const
{
    return
        globalCell >= 0 && globalCell < globalToLocalCell_.size()
      ? globalToLocalCell_[globalCell]
      : -1;
}

inline label dsmcLocalMesh::toGlobal(label localCell) const
{
    return localToGlobalCell_[localCell];
}

inline bool dsmcLocalMesh::isHalo(label localCell) const
{
    return localCell >= 0 && isHaloCell_.test(localCell);
}

inline label dsmcLocalMesh::nCells() const
{
    return localToGlobalCell_.size();
}

} // End namespace Foam
