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

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"
#include "DynamicList.H"

#include <map>

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::boundaryLayers::createNewFacesAndCells
(
    const boolList& treatPatches
)
{
    if (omp_in_parallel())
    {
        FatalErrorInFunction
            << "Building cellsToAdd inside a parallel region."
            << " This is not thread safe" << exit(FatalError);
    }

    Info<< "Starting creating layer cells" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelListList2D& faceEdges = mse.faceEdges();
    const labelListList2D& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const labelList& faceOwners = mse.faceOwners();

    // create lists for new boundary faces
    DynamicList<DynamicList<label>> newBoundaryFaces;
    DynamicList<label> newBoundaryOwners;
    DynamicList<label> newBoundaryPatches;

    // create storage for new cells
    DynamicList<faceList> cellsToAdd;

    // create layer cells and store boundary faces
    const label nOldCells = mesh_.cells().size();
    forAll(bFaces, bfI)
    {
        if (treatPatches[boundaryFacePatches[bfI]])
        {
            const face& f = bFaces[bfI];

            faceList cellFaces(f.size() + 2);

            label fI(0);

            // store boundary face
            cellFaces[fI++] = f.reverseFace();

            // create parallel face
            face newF(f.size());
            forAll(f, pI)
                newF[pI] = newLabelForVertex_[f[pI]];
            cellFaces[fI++] = newF;

            appendList(newBoundaryFaces,newF);
            newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);

            // create quad faces
            newF.setSize(4);
            forAll(f, pI)
            {
                newF[0] = f[pI];
                newF[1] = f.nextLabel(pI);
                newF[2] = newLabelForVertex_[f.nextLabel(pI)];
                newF[3] = newLabelForVertex_[f[pI]];

                cellFaces[fI++] = newF;

                // check if the face is at the boundary
                // of the treated partitions
                const label edgeI = faceEdges[bfI][pI];
                if (edgeFaces[edgeI].size() == 2)
                {
                    label neiFace = edgeFaces[edgeI][0];
                    if (neiFace == bfI)
                    {
                        neiFace = edgeFaces[edgeI][1];
                    }

                    if (!treatPatches[boundaryFacePatches[neiFace]])
                    {
                        appendList(newBoundaryFaces,newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(boundaryFacePatches[neiFace]);
                    }
                }
            }

            cellsToAdd.append(cellFaces);
        }
        else
        {
            appendList(newBoundaryFaces,bFaces[bfI]);
            newBoundaryOwners.append(faceOwners[bfI]);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);
        }
    }

    // create mesh modifier
    polyMeshGenModifier meshModifier(mesh_);

    meshModifier.addCells(cellsToAdd);
    cellsToAdd.clear();
    meshModifier.reorderBoundaryFaces();
    meshModifier.replaceBoundary
    (
        patchNames_,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        boundaries[patchI].patchType() = patchTypes_[patchI];
    }

    // delete meshSurfaceEngine
    this->clearOut();

    Info<< "Finished creating layer cells" << endl;
}


// ************************************************************************* //
