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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::polyMeshGenModifier::removeFaces(const boolList& removeFace)
{
    //mesh_.clearOut();

    Info<< "Removing faces" << endl;
    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    label nFaces(0);
    labelLongList newFaceLabel(faces.size(), -1);

    // copy internal faces
    const label nInternalFaces = mesh_.nInternalFaces();
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        if (!removeFace[faceI])
        {
            if (nFaces < faceI)
                faces[nFaces].transfer(faces[faceI]);

            newFaceLabel[faceI] = nFaces;
            ++nFaces;
        }
    }

    // copy boundary faces
    PtrList<boundaryPatch>& boundaries = mesh_.boundaries_;
    labelList patchStart(boundaries.size());
    labelList nFacesInPatch(boundaries.size());
    label npI(0);
    forAll(boundaries, patchI)
    {
        const label oldStart = boundaries[patchI].patchStart();
        const label oldNumFacesInPatch = boundaries[patchI].patchSize();

        patchStart[npI] = nFaces;
        nFacesInPatch[npI] = 0;

        for (label faceI = 0; faceI < oldNumFacesInPatch; ++faceI)
        {
            if (!removeFace[oldStart + faceI])
            {
                ++nFacesInPatch[npI];
                if (nFaces <(oldStart + faceI))
                {
                    faces[nFaces].transfer(faces[oldStart + faceI]);
                }
                newFaceLabel[oldStart + faceI] = nFaces++;
            }
        }

        ++npI;
    }

    forAll(boundaries, patchI)
    {
        boundaries[patchI].patchStart() = patchStart[patchI];
        boundaries[patchI].patchSize() = nFacesInPatch[patchI];
    }

    faces.setSize(nFaces);

    // update face subsets in the mesh
    mesh_.updateFaceSubsets(newFaceLabel);

    // change cells
    # pragma omp parallel for if (cells.size() > 1000) schedule(dynamic, 40)
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        DynList<label> newC;

        forAll(c, fI)
        {
            if (newFaceLabel[c[fI]] != -1)
            {
                newC.append(newFaceLabel[c[fI]]);
            }
        }

        c.setSize(newC.size());

        forAll(c, fI)
        {
            c[fI] = newC[fI];
        }
    }

    mesh_.clearOut();

    Info<< "Finished removing faces" << endl;
}


void Foam::Module::polyMeshGenModifier::removeDuplicateFaces()
{
    VRWGraph& pointFaces = this->pointFaces();

    faceListPMG& faces = mesh_.faces_;

    labelLongList newFaceLabel(faces.size(), -1);

    label nFaces(0);
    forAll(faces, faceI)
    {
        if (newFaceLabel[faceI] != -1)
        {
            continue;
        }

        const face& f = faces[faceI];
        DynList<label> duplicates;
        forAllRow(pointFaces, f[0], pfI)
        {
            const label faceJ = pointFaces(f[0], pfI);

            if (faceI == faceJ)
            {
                continue;
            }

            if (faces[faceJ] == f)
            {
                duplicates.append(faceJ);
            }
        }

        newFaceLabel[faceI] = nFaces;
        forAll(duplicates, i)
        {
            newFaceLabel[duplicates[i]] = nFaces;
        }
        ++nFaces;
    }

    label nInternalFaces(faces.size());
    if (mesh_.boundaries_.size() != 0)
    {
        nInternalFaces = mesh_.boundaries_[0].patchStart();
    }

    // remove internal faces
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        if (newFaceLabel[faceI] >= faceI)
        {
            continue;
        }

        faces[newFaceLabel[faceI]].transfer(faces[faceI]);
    }

    // update boundary faces (they cannot be duplicated)
    forAll(mesh_.boundaries_, patchI)
    {
        boundaryPatch& patch = mesh_.boundaries_[patchI];

        const label start = patch.patchStart();
        const label end = start + patch.patchSize();

        patch.patchStart() = newFaceLabel[start];
        for (label faceI = start; faceI < end; ++faceI)
        {
            if (newFaceLabel[faceI] <  faceI)
            {
                faces[newFaceLabel[faceI]].transfer(faces[faceI]);
            }
        }
    }

    faces.setSize(nFaces);
    mesh_.updateFaceSubsets(newFaceLabel);

    // change cells
    cellListPMG& cells = mesh_.cells_;
    # pragma omp parallel for if (cells.size() > 1000) schedule(dynamic, 40)
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        DynList<label> newC;

        forAll(c, fI)
        {
            if (newFaceLabel[c[fI]] != -1)
            {
                newC.append(newFaceLabel[c[fI]]);
            }
        }

        c.setSize(newC.size());

        forAll(c, fI)
        {
            c[fI] = newC[fI];
        }
    }

    mesh_.clearOut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
