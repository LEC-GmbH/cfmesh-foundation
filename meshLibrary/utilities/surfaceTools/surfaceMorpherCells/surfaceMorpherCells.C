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

#include "surfaceMorpherCells.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::surfaceMorpherCells::replaceMeshBoundary()
{
    wordList patchNames(1);
    patchNames[0] = "defaultFaces";

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames,
        newBoundaryFaces_,
        newBoundaryOwners_,
        newBoundaryPatches_
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Module::surfaceMorpherCells::surfaceMorpherCells
(
    polyMeshGen& mesh
)
    :
    mesh_(mesh),
    nIntFaces_(0),
    boundaryVertex_(mesh.points().size()),
    cellFlags_(mesh.cells().size()),
    newBoundaryFaces_(),
    newBoundaryOwners_(),
    newBoundaryPatches_()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::surfaceMorpherCells::morphMesh()
{
    // perform surface morphing
    bool changed;

    do
    {
        changed = false;

        findBoundaryVertices();

        findBoundaryCells();

        if (removeCellsWithAllVerticesAtTheBoundary())
        {
            changed = true;
            continue;
        }

        if (morphInternalFaces())
        {
            changed = true;
            continue;
        }

        if (morphBoundaryFaces())
        {
            changed = true;
            continue;
        }

    } while (changed);

    polyMeshGenModifier(mesh_).removeUnusedVertices();
}


// ************************************************************************* //
