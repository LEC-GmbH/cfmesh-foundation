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
#include "helperFunctions.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::surfaceMorpherCells::findBoundaryVertices()
{
    const faceListPMG& faces = mesh_.faces();

    boundaryVertex_.setSize(mesh_.points().size());
    boundaryVertex_ = false;

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for (label faceI = start; faceI < end; ++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                boundaryVertex_[f[pI]] = true;
            }
        }
    }
}


void Foam::Module::surfaceMorpherCells::findBoundaryCells()
{
    const labelList& owner = mesh_.owner();

    cellFlags_.setSize(mesh_.cells().size());
    cellFlags_ = NONE;

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for (label faceI = start; faceI < end; ++faceI)
        {
            cellFlags_[owner[faceI]] = BOUNDARY;
        }
    }
}


bool Foam::Module::surfaceMorpherCells::morphInternalFaces()
{
    Info<< "Morphing internal faces" << endl;

    bool changed(false);

    newBoundaryFaces_.setSize(0);
    newBoundaryOwners_.setSize(0);
    newBoundaryPatches_.setSize(0);

    const label nIntFaces = mesh_.nInternalFaces();
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    // copy boundary faces
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for (label faceI = start; faceI < end; ++faceI)
        {
            appendList(newBoundaryFaces_, faces[faceI]);
            newBoundaryOwners_.append(owner[faceI]);
            newBoundaryPatches_.append(0);
        }
    }

    // start morphing internal faces
    for (label faceI = 0; faceI < nIntFaces; ++faceI)
    {
        const face& f = faces[faceI];

        DynamicList<bool> removeFaceVertex(f.size(), false);

        face newF(f.size());
        label i(0);

        forAll(f, pI)
        {
            if (
                boundaryVertex_[f.prevLabel(pI)] &&
                boundaryVertex_[f[pI]] &&
                boundaryVertex_[f.nextLabel(pI)]
            )
            {
                removeFaceVertex[pI] = true;
            }
            else
            {
                newF[i++] = f[pI];
            }
        }

        if (i < f.size())
        {
            changed = true;

            // store shrinked face
            newF.setSize(i);

            // create new boundary faces from the removed part
            label mat(1);
            DynamicList<direction> nodeMaterial(f.size(), direction(0));
            DynamicList<DynamicList<edge>> edgeMats;
            forAll(nodeMaterial, nI)
            {
                if (!nodeMaterial[nI] && removeFaceVertex[nI])
                {
                    edgeMats.append(DynamicList<edge>());
                    DynamicList<label> front;
                    front.append(nI);

                    do
                    {
                        DynamicList<label> newFront;
                        forAll(front, pI)
                        {
                            const label fLabel = front[pI];
                            if (nodeMaterial[fLabel])
                            {
                                continue;
                            }
                            nodeMaterial[fLabel] = mat;

                            const edge& e0 = f.faceEdge(fLabel);
                            if (findIndex(edgeMats[mat - 1], e0) == -1)
                            {
                                edgeMats[mat - 1].append(e0);
                            }

                            const edge& e1 = f.faceEdge(f.rcIndex(fLabel));
                            if (findIndex(edgeMats[mat - 1], e1) == -1)
                            {
                                edgeMats[mat - 1].append(e1);
                            }

                            if (removeFaceVertex[f.rcIndex(fLabel)])
                            {
                                newFront.append(f.rcIndex(fLabel));
                            }
                            if (removeFaceVertex[f.fcIndex(fLabel)])
                            {
                                newFront.append(f.fcIndex(fLabel));
                            }
                        }
                        front = newFront;
                    }
                    while (front.size() != 0);

                    ++mat;
                }
            }

            forAll(edgeMats, matI)
            {
                edgeMats[matI].shrink();
                const face rf = help::removeEdgesFromFace(f, edgeMats[matI]);
                const face newBf = help::createFaceFromRemovedPart(f, rf);

                appendList(newBoundaryFaces_, newBf);
                newBoundaryOwners_.append(owner[faceI]);
                newBoundaryPatches_.append(0);

                appendList(newBoundaryFaces_, newBf.reverseFace());
                newBoundaryOwners_.append(neighbour[faceI]);
                newBoundaryPatches_.append(0);
            }

            const_cast<face&>(faces[faceI]) = newF;
        }
    }

    polyMeshGenModifier meshModifier(mesh_);

    if (changed)
    {
        // replace boundary of the mesh
        replaceMeshBoundary();

        // remove faces which do not exist any more
        boolList removeFace(faces.size(), false);
        bool removeFaces(false);

        for (label faceI = 0; faceI < nIntFaces; ++faceI)
        {
            if (faces[faceI].size() < 3)
            {
                removeFace[faceI] = true;
                removeFaces = true;
            }
        }

        if (removeFaces)
        {
            meshModifier.removeFaces(removeFace);
        }
    }

    Info<< "Finished morphing internal faces" << endl;

    return changed;
}


// ************************************************************************* //
