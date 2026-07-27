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

#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"
#include "polyMeshGenAddressing.H"
#include "renumberMethod.H"
#include "ListOps.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

// Determine the upper-triangular face order for a given new cell order.
// Ported from OpenFOAM-dev's applications/utilities/mesh/manipulation/
// renumberMesh/renumberMesh.C (getFaceOrder()) - that function is
// application-private (not part of a linkable library) and takes a
// Foam::primitiveMesh, so it is adapted here to read directly from
// polyMeshGen's own addressing. Ordering of boundary faces is left
// unchanged, matching the upstream function's documented behaviour.
Foam::labelList getFaceOrder
(
    const Foam::cellList& newCells,           // cell face-lists, new cell
                                               // order, old face numbers
    const Foam::labelList& oldOwner,          // old face owner
    const Foam::labelList& oldNeighbour,      // old face neighbour
    const Foam::DynamicList<Foam::label>& reverseCellOrder, // old cell ->
                                                             // new cell
    const Foam::label nInternalFaces,
    const Foam::label nFaces
)
{
    using namespace Foam;

    labelList oldToNewFace(nFaces, -1);

    label newFaceI = 0;

    labelList nbr;
    labelList order;

    forAll(newCells, newCellI)
    {
        const cell& cFaces = newCells[newCellI];

        nbr.setSize(cFaces.size());

        forAll(cFaces, i)
        {
            const label faceI = cFaces[i];

            if (faceI < nInternalFaces)
            {
                label nbrCellI = reverseCellOrder[oldNeighbour[faceI]];
                if (nbrCellI == newCellI)
                {
                    nbrCellI = reverseCellOrder[oldOwner[faceI]];
                }

                nbr[i] = (newCellI < nbrCellI) ? nbrCellI : -1;
            }
            else
            {
                // boundary face - handled below
                nbr[i] = -1;
            }
        }

        order.setSize(nbr.size());
        sortedOrder(nbr, order);

        forAll(order, i)
        {
            const label index = order[i];
            if (nbr[index] != -1)
            {
                oldToNewFace[cFaces[index]] = newFaceI++;
            }
        }
    }

    // Leave boundary faces intact (same patch grouping/order as before)
    for (label faceI = newFaceI; faceI < nFaces; ++faceI)
    {
        oldToNewFace[faceI] = faceI;
    }

    forAll(oldToNewFace, faceI)
    {
        if (oldToNewFace[faceI] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position for face " << faceI
                << abort(FatalError);
        }
    }

    return invert(nFaces, oldToNewFace);
}

} // End anonymous namespace


void Foam::Module::polyMeshGenModifier::renumberMesh()
{
    Info<< "Renumbering the mesh" << endl;

    // Determine the new cell order using OpenFOAM's Cuthill-McKee
    // bandwidth-reduction renumbering, operating directly on cfMesh's own
    // cell-adjacency addressing (no OpenFOAM polyMesh/fvMesh needed - the
    // renumberMethod interface used here only requires a cell-adjacency
    // graph and cell centres).
    const labelListList2D& cCells = mesh_.addressingData().cellCells();
    const vectorField& cellCentres = mesh_.addressingData().cellCentres();

    labelListList cellCells(cCells.size());
    forAll(cellCells, cellI)
    {
        labelList& cc = cellCells[cellI];
        cc.setSize(cCells[cellI].size());

        forAll(cc, i)
        {
            cc[i] = cCells[cellI][i];
        }
    }

    dictionary renumberDict;
    renumberDict.add("method", word("CuthillMcKee"));
    autoPtr<renumberMethod> renumberPtr = renumberMethod::New(renumberDict);

    // Order in which cells need to be visited, i.e. new cell -> old cell
    const labelList newOrder = renumberPtr().renumber(cellCells, cellCentres);

    cellListPMG& oldCells = this->cellsAccess();
    const labelList& oldOwner = mesh_.owner();
    const labelList& oldNeighbour = mesh_.neighbour();

    cellList newCells(oldCells.size());

    // The reverse order list gives the new cell label for every old cell
    DynamicList<label> reverseOrder;
    reverseOrder.setSize(newOrder.size());

    forAll(newOrder, cellI)
    {
        newCells[cellI].transfer(oldCells[newOrder[cellI]]);

        reverseOrder[newOrder[cellI]] = cellI;
    }

    // Face order gives the old face label for every new face
    const label nIntFaces = mesh_.nInternalFaces();
    const labelList faceOrder = getFaceOrder
    (
        newCells,
        oldOwner,
        oldNeighbour,
        reverseOrder,
        nIntFaces,
        oldOwner.size()
    );

    // Renumber the cells
    labelList reverseFaceOrder(invert(faceOrder.size(), faceOrder));
    forAll(newCells, cellI)
    {
        cell& c = newCells[cellI];

        forAll(c, fI)
        {
            c[fI] = reverseFaceOrder[c[fI]];
        }
    }

    faceListPMG& oldFaces = this->facesAccess();
    faceList newFaces(oldFaces.size());

    forAll(newFaces, faceI)
    {
        newFaces[faceI].transfer(oldFaces[faceOrder[faceI]]);
    }

    // Turn the faces that need to be turned
    // Only loop through internal faces
    forAll(oldNeighbour, faceI)
    {
        const label oldFaceI = faceOrder[faceI];
        if (oldNeighbour[oldFaceI] < 0)
        {
            continue;
        }

        if
        (
            reverseOrder[oldNeighbour[oldFaceI]]
          < reverseOrder[oldOwner[oldFaceI]]
        )
        {
            newFaces[faceI] = newFaces[faceI].reverseFace();
        }
    }

    // transfer faces and cells back to the original lists
    forAll(newCells, cellI)
    {
        oldCells[cellI].transfer(newCells[cellI]);
    }
    forAll(newFaces, faceI)
    {
        oldFaces[faceI].transfer(newFaces[faceI]);
    }

    mesh_.updateFaceSubsets(reverseFaceOrder);
    mesh_.updateCellSubsets(reverseOrder);
    this->clearOut();
    mesh_.clearOut();

    Info<< "Finished renumbering the mesh" << endl;
}


// ************************************************************************* //
