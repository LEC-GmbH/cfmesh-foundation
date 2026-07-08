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

#include "cartesianMeshExtractor.H"
#include "demandDrivenData.H"
#include "meshOctree.H"

#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "decomposeFaces.H"
#include "decomposeCells.H"

#include <map>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::cartesianMeshExtractor::createPolyMesh()
{
    Info<< "Creating polyMesh from octree" << endl;

    const meshOctree& octree = octreeCheck_.octree();

    // give labels to cubes which will be used as mesh cells
    const List<direction>& cType = octreeCheck_.boxType();

    labelList& leafCellLabel = *leafCellLabelPtr_;
    label nCells(0);
    forAll(cType, leafI)
    {
        if (cType[leafI] & meshOctreeAddressing::MESHCELL)
        {
            leafCellLabel[leafI] = nCells++;
        }
    }

    // access to mesh data
    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    // start creating octree mesh
    cells.setSize(nCells);
    List<direction> nFacesInCell(nCells, direction(0));
    label nFaces(0);

    const VRWGraph& octreeFaces = octreeCheck_.octreeFaces();
    const labelLongList& owner = octreeCheck_.octreeFaceOwner();
    const labelLongList& neighbour = octreeCheck_.octreeFaceNeighbour();

    forAll(octreeFaces, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const label ownLabel = leafCellLabel[own];
        label neiLabel(-1);
        if (nei != -1)
            neiLabel = leafCellLabel[nei];

        if ((ownLabel != -1) && (neiLabel != -1))
        {
            ++nFaces;
            ++nFacesInCell[ownLabel];
            ++nFacesInCell[neiLabel];
        }
        else if (ownLabel != -1)
        {
            ++nFaces;
            ++nFacesInCell[ownLabel];
        }
        else if (neiLabel != -1)
        {
            ++nFaces;
            ++nFacesInCell[neiLabel];
        }
    }

    faces.setSize(nFaces);
    forAll(cells, cI)
        cells[cI].setSize(nFacesInCell[cI]);
    nFacesInCell = 0;

    nFaces = 0;

    forAll(octreeFaces, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const label ownLabel = leafCellLabel[own];
        label neiLabel(-1);
        if (nei != -1)
            neiLabel = leafCellLabel[nei];

        if ((ownLabel != -1) && (neiLabel != -1))
        {
            // internal face
            faces[nFaces].setSize(octreeFaces.sizeOfRow(faceI));
            forAllRow(octreeFaces, faceI, pI)
                faces[nFaces][pI] = octreeFaces(faceI, pI);

            cells[ownLabel][nFacesInCell[ownLabel]++] = nFaces;
            cells[neiLabel][nFacesInCell[neiLabel]++] = nFaces;
            ++nFaces;
        }
        else if (ownLabel != -1)
        {
            // boundary face
            faces[nFaces].setSize(octreeFaces.sizeOfRow(faceI));
            forAllRow(octreeFaces, faceI, pI)
                faces[nFaces][pI] = octreeFaces(faceI, pI);

            cells[ownLabel][nFacesInCell[ownLabel]++] = nFaces;
            ++nFaces;
        }
        else if (neiLabel != -1)
        {
            // boundary face
            faces[nFaces].setSize(octreeFaces.sizeOfRow(faceI));
            faces[nFaces][0] = octreeFaces(faceI, 0);
            for (label pI = octreeFaces.sizeOfRow(faceI)-1; pI > 0; --pI)
                faces[nFaces][octreeFaces.sizeOfRow(faceI)-pI] =
                    octreeFaces(faceI, pI);

            cells[neiLabel][nFacesInCell[neiLabel]++] = nFaces;
            ++nFaces;
        }
    }

    meshModifier.reorderBoundaryFaces();

    if (octree.isQuadtree())
    {
        // generate empty patches
        // search for faces with a dominant z coordinate and store them
        // into an empty patch
        meshSurfaceEngine mse(mesh_);
        const vectorField& fNormals = mse.faceNormals();
        const faceList::subList& bFaces = mse.boundaryFaces();
        const labelList& fOwner = mse.faceOwners();
        const vectorField& fCentres = mse.faceCentres();

        const boundBox& bb = octree.rootBox();
        const scalar tZ = 0.05*(bb.max().z() - bb.min().z());

        wordList patchNames(3);
        patchNames[0] = "defaultFaces";
        patchNames[1] = "unusedFacesBottom";
        patchNames[2] = "unusedFacesTop";

        VRWGraph boundaryFaces;
        labelLongList newFaceOwner;
        labelLongList newFacePatch;

        forAll(fNormals, bfI)
        {
            // store the face and its owner
            boundaryFaces.appendList(bFaces[bfI]);
            newFaceOwner.append(fOwner[bfI]);

            const vector& fNormal = fNormals[bfI];

            if (Foam::mag(fNormal.z()) > Foam::mag(fNormal.x() + fNormal.y()))
            {
                if (Foam::mag(fCentres[bfI].z() - bb.min().z()) < tZ)
                {
                    newFacePatch.append(1);
                }
                else if (Foam::mag(fCentres[bfI].z() - bb.max().z()) < tZ)
                {
                    newFacePatch.append(2);
                }
                else
                {
                    FatalErrorInFunction
                        << "Cannot distribute the face!!" << exit(FatalError);
                }
            }
            else
            {
                newFacePatch.append(0);
            }
        }

        // replace the boundary with faces in correct patches
        meshModifier.replaceBoundary
        (
            patchNames,
            boundaryFaces,
            newFaceOwner,
            newFacePatch
        );

        meshModifier.boundariesAccess()[1].patchType() = "empty";
        meshModifier.boundariesAccess()[2].patchType() = "empty";
    }

    Info<< "Finished creating polyMesh" << endl;
}


// ************************************************************************* //
