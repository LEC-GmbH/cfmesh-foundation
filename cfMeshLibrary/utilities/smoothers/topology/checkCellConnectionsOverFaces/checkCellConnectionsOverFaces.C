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

#include "checkCellConnectionsOverFaces.H"

#include <set>
#include <map>

#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

namespace meshConnectionsHelper
{

class meshConnectionsNeighbourOperator
{
    const polyMeshGen& mesh_;

public:

    meshConnectionsNeighbourOperator(const polyMeshGen& mesh)
    :
        mesh_(mesh)
    {}


    label size() const
    {
        return mesh_.cells().size();
    }

    void operator()(const label cellI, DynamicList<label>& neighbourCells) const
    {
        neighbourCells.clear();

        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        const cell& c = mesh_.cells()[cellI];

        forAll(c, fI)
        {
            label nei = owner[c[fI]];

            if (nei == cellI)
            {
                nei = neighbour[c[fI]];
            }

            if (nei >= 0)
            {
                neighbourCells.append(nei);
            }
        }
    }
};


class meshConnectionsSelectorOperator
{
public:

    meshConnectionsSelectorOperator()
    {}


    bool operator()(const label /*cellI*/) const
    {
        return true;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshConnectionsHelper
} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::Module::checkCellConnectionsOverFaces::findCellGroups()
{
    Info<< "Checking cell connections" << endl;

    mesh_.owner();
    nGroups_ =
        help::groupMarking
        (
            cellGroup_,
            meshConnectionsHelper::meshConnectionsNeighbourOperator(mesh_),
            meshConnectionsHelper::meshConnectionsSelectorOperator()
        );

    Info<< "Finished checking cell connections" << endl;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::Module::checkCellConnectionsOverFaces::checkCellConnectionsOverFaces
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    cellGroup_(mesh.cells().size(), -1),
    nGroups_(0)
{
    findCellGroups();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Module::checkCellConnectionsOverFaces::checkCellGroups
(
    const bool haveInsidePoint,
    const point& insidePoint
)
{
    if (nGroups_ == 1)
    {
        return false;
    }

    Warning << "Mesh has " << nGroups_ << " unconnected regions" << endl;

    label keepGroup(-1);

    if (haveInsidePoint)
    {
        // keep the group containing the cell whose centre is nearest
        // to insidePoint
        const vectorField& cellCentres = mesh_.addressingData().cellCentres();

        label nearestCell(-1);
        scalar nearestDistSqr(GREAT);

        forAll(cellCentres, cellI)
        {
            const scalar distSqr = magSqr(cellCentres[cellI] - insidePoint);

            if (distSqr < nearestDistSqr)
            {
                nearestDistSqr = distSqr;
                nearestCell = cellI;
            }
        }

        keepGroup = cellGroup_[nearestCell];
    }
    else
    {
        // keep the group which has the most cells
        labelList nCellsInGroup(nGroups_, 0);

        forAll(cellGroup_, cI)
        {
            ++nCellsInGroup[cellGroup_[cI]];
        }

        label maxNumCells(-1);
        forAll(nCellsInGroup, groupI)
        {
            if (nCellsInGroup[groupI] > maxNumCells)
            {
                maxNumCells = nCellsInGroup[groupI];
                keepGroup = groupI;
            }
        }
    }

    // remove cells which are not in the group which shall be kept
    boolList removeCell(mesh_.cells().size(), false);
    forAll(cellGroup_, cellI)
    {
        if (cellGroup_[cellI] != keepGroup)
        {
            removeCell[cellI] = true;
        }
    }

    polyMeshGenModifier(mesh_).removeCells(removeCell);

    return true;
}


// ************************************************************************* //
