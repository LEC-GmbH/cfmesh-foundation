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
#include "meshOctree.H"
#include "boundBox.H"

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
    const pointField& insidePoints,
    const meshOctree& octree
)
{
    if (nGroups_ == 1)
    {
        return false;
    }

    Warning << "Mesh has " << nGroups_ << " unconnected regions" << endl;

    labelHashSet keepGroups;

    if (insidePoints.size())
    {
        // keep the union of groups containing the cell that actually
        // contains each insidePoint
        const labelList resolvedCells =
            mesh_.addressingData().findCells(insidePoints, octree);

        forAll(insidePoints, pI)
        {
            const label cellI = resolvedCells[pI];

            if (cellI == -1)
            {
                FatalErrorInFunction
                    << "insidePoint " << insidePoints[pI]
                    << " is not inside the mesh." << nl
                    << "Mesh bounding box: " << boundBox(mesh_.points())
                    << exit(FatalError);
            }

            keepGroups.insert(cellGroup_[cellI]);
        }
    }
    else
    {
        // keep the group which has the most cells
        labelList nCellsInGroup(nGroups_, 0);

        forAll(cellGroup_, cI)
        {
            ++nCellsInGroup[cellGroup_[cI]];
        }

        label keepGroup(-1);
        label maxNumCells(-1);
        forAll(nCellsInGroup, groupI)
        {
            if (nCellsInGroup[groupI] > maxNumCells)
            {
                maxNumCells = nCellsInGroup[groupI];
                keepGroup = groupI;
            }
        }

        keepGroups.insert(keepGroup);
    }

    // remove cells which are not in a group which shall be kept. Keeping
    // more than one group (multiple insidePoints resolving to different
    // groups) is a legitimate, stable end state - nGroups_ > 1 no longer
    // implies more cells still need removing, so whether anything was
    // actually marked for removal (not just whether nGroups_ > 1) is
    // what the caller's do-while(changed) loop needs to terminate
    boolList removeCell(mesh_.cells().size(), false);
    bool anyRemoved(false);
    forAll(cellGroup_, cellI)
    {
        if (!keepGroups.found(cellGroup_[cellI]))
        {
            removeCell[cellI] = true;
            anyRemoved = true;
        }
    }

    if (!anyRemoved)
    {
        return false;
    }

    polyMeshGenModifier(mesh_).removeCells(removeCell);

    return true;
}


// ************************************************************************* //
