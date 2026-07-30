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

#include "decomposeCells.H"
#include "demandDrivenData.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"
#include "decomposeFaces.H"
#include "DynamicList.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::decomposeCells::decomposeMesh(const boolList& decomposeCell)
{
    checkFaceConnections(decomposeCell);

    createPointsAndCellFaces(decomposeCell);

    storeBoundaryFaces(decomposeCell);

    removeDecomposedCells(decomposeCell);

    addNewCells();
}


void Foam::Module::decomposeCells::checkFaceConnections
(
    const boolList& decomposeCell
)
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    boolList decomposeFace(faces.size(), false);
    forAll(cells, cellI)
    {
        if (decomposeCell[cellI])
        {
            DynamicList<label> vrt;
            DynamicList<edge> edges;
            labelListList2D faceEdges;
            DynamicList<DynamicList<label>> edgeFaces;

            findAddressingForCell(cellI, vrt, edges, faceEdges, edgeFaces);

            forAll(faceEdges, fI)
            {
                const DynamicList<label>& fEdges = faceEdges[fI];

                labelHashSet neiFaces;
                forAll(fEdges, feI)
                {
                    label neiFace = edgeFaces[fEdges[feI]][0];
                    if (neiFace == fI)
                    {
                        neiFace = edgeFaces[fEdges[feI]][1];
                    }

                    if (neiFaces.found(neiFace))
                    {
                        decomposeFace[cells[cellI][fI]] = true;
                    }
                    else
                    {
                        neiFaces.insert(neiFace);
                    }
                }
            }
        }
    }

    // decompose faces which would cause invalid connections
    decomposeFaces(mesh_).decomposeMeshFaces(decomposeFace);
}


void Foam::Module::decomposeCells::createPointsAndCellFaces
(
    const boolList& decomposeCell
)
{
    if (omp_in_parallel())
    {
        FatalErrorInFunction
            << "Building facesOfNewCells_ inside a parallel region."
            << " This is not thread safe" << exit(FatalError);
    }

    facesOfNewCells_.clear();

    forAll(decomposeCell, cI)
    {
        if (decomposeCell[cI])
        {
            decomposeCellIntoPyramids(cI);
        }
    }
}


void Foam::Module::decomposeCells::storeBoundaryFaces
(
    const boolList& /*decomposeCell*/
)
{
    meshSurfaceEngine mse(mesh_);
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch = mse.boundaryFacePatches();

    forAll(bFaces, bfI)
    {
        appendList(newBoundaryFaces_, bFaces[bfI]);
        newBoundaryPatches_.append(facePatch[bfI]);
    }
}


void Foam::Module::decomposeCells::removeDecomposedCells
(
    const boolList& decomposeCell
)
{
    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.removeCells(decomposeCell);
}


void Foam::Module::decomposeCells::addNewCells()
{
    Info<< "Adding new cells " << endl;
    polyMeshGenModifier(mesh_).addCells(facesOfNewCells_);
    facesOfNewCells_.clear();
    Info<< "Reordering bnd faces" << endl;
    polyMeshGenModifier(mesh_).reorderBoundaryFaces();

    Info<< "Finding bnd faces" << endl;
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const labelListList2D& pointFaces = mesh_.addressingData().pointFaces();

    DynamicList<label> newBoundaryOwners;

    forAll(newBoundaryFaces_, faceI)
    {
        face bf(newBoundaryFaces_[faceI].size());
        forAllRow(newBoundaryFaces_, faceI, pI)
        {
            bf[pI] = newBoundaryFaces_[faceI][pI];
        }

        forAllRow(pointFaces, bf[0], pfI)
        {
            const label fLabel = pointFaces[bf[0]][pfI];
            if ((mesh_.faceIsInPatch(fLabel) != -1) && (bf == faces[fLabel]))
            {
                newBoundaryOwners.append(owner[fLabel]);
            }
        }
    }

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames_,
        newBoundaryFaces_,
        newBoundaryOwners,
        newBoundaryPatches_
    );

    polyMeshGenModifier(mesh_).removeUnusedVertices();
    polyMeshGenModifier(mesh_).clearAll();

    PtrList<boundaryPatch>& boundaries =
        polyMeshGenModifier(mesh_).boundariesAccess();
    forAll(boundaries, patchI)
    {
        boundaries[patchI].patchType() = patchTypes_[patchI];
    }
}


// ************************************************************************* //
