/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | www.cfmesh.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 Creative Fields, Ltd.
-------------------------------------------------------------------------------
Author
     Franjo Juretic (franjo.juretic@c-fields.com)

License
    This file is part of OpenFOAM.

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::Module::cellListPMG&
Foam::Module::polyMeshGenCells::cells() const
{
    return cells_;
}


inline void Foam::Module::polyMeshGenCells::addCellToSubset
(
    const label selID,
    const label cellI
)
{
    auto it = cellSubsets_.find(selID);
    if (it == cellSubsets_.end())
    {
        return;
    }

    it().addElement(cellI);
}


inline void Foam::Module::polyMeshGenCells::removeCellFromSubset
(
    const label selID,
    const label cellI
)
{
    auto it = cellSubsets_.find(selID);
    if (it == cellSubsets_.end())
    {
        return;
    }

    it().removeElement(cellI);
}


inline void Foam::Module::polyMeshGenCells::cellInSubsets
(
    const label cellI,
    DynamicList<label>& cellSubsets
) const
{
    cellSubsets.clear();

    forAllConstIter(Map<meshSubset>, cellSubsets_, it)
    {
        if (it().found(cellI))
        {
            cellSubsets.append(it.key());
        }
    }
}


inline void Foam::Module::polyMeshGenCells::cellSubsetIndices
(
    DynamicList<label>& indices
) const
{
    indices.clear();

    forAllConstIter(Map<meshSubset>, cellSubsets_, it)
    {
        indices.append(it.key());
    }
}


template<class ListType>
inline void Foam::Module::polyMeshGenCells::cellsInSubset
(
    const label selID,
    ListType& cellLabels
) const
{
    cellLabels.clear();

    auto it = cellSubsets_.find(selID);
    if (it == cellSubsets_.end())
    {
        return;
    }

    it().containedElements(cellLabels);
}


template<class ListType>
inline void Foam::Module::polyMeshGenCells::updateCellSubsets
(
    const ListType& newCellLabels
)
{
    forAllIter(Map<meshSubset>, cellSubsets_, it)
    {
        it().updateSubset(newCellLabels);
    }
}


inline void Foam::Module::polyMeshGenCells::updateCellSubsets
(
    const labelListList2D& newCellsForCell
)
{
    forAllIter(Map<meshSubset>, cellSubsets_, it)
    {
        it().updateSubset(newCellsForCell);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
