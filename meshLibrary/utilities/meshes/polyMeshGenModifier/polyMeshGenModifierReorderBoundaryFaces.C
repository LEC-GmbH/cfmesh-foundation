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
#include "DynamicList.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::polyMeshGenModifier::reorderBoundaryFaces()
{
    Info<< "Reordering boundary faces " << endl;

    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    const labelList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();

    // count internal and boundary faces
    const label numBFaces = faces.size() - nInternalFaces;

    DynamicList<label> newFaceLabel(faces.size(), -1);

    // find faces which should be repositioned
    label nReplaced(0);
    labelList internalToChange;
    labelList boundaryToChange;

    const label nThreads = 3*omp_get_num_procs();
    labelList nInternalToChangeThread(nThreads);
    labelList nBoundaryToChangeThread(nThreads);

    # pragma omp parallel num_threads(nThreads)
    {
        const label threadI = omp_get_thread_num();

        label& nItc = nInternalToChangeThread[threadI];
        label& nBtc = nBoundaryToChangeThread[threadI];

        DynamicList<label> internalToChangeLocal, boundaryToChangeLocal;

        // find the boundary faces within the range of internal faces
        # pragma omp for schedule(static)
        for (label faceI = 0; faceI < nInternalFaces; ++faceI)
        {
            if (neighbour[faceI] == -1)
            {
                internalToChangeLocal.append(faceI);
            }
        }

        nItc = internalToChangeLocal.size();

        // find the internal faces within the range of boundary faces
        # pragma omp for schedule(static)
        for (label faceI = nInternalFaces; faceI < faces.size(); ++faceI)
        {
            if (neighbour[faceI] != -1)
            {
                boundaryToChangeLocal.append(faceI);
            }
        }

        nBtc = boundaryToChangeLocal.size();

        // perform reduction such that all threads know how many faces
        // need to be swapped
        # pragma omp critical
        nReplaced += nBtc;

        # pragma omp barrier

        # pragma omp master
        {
            internalToChange.setSize(nReplaced);
            boundaryToChange.setSize(nReplaced);
        }

        # pragma omp barrier

        label localStart(0);
        for (label i = 0; i < threadI; ++i)
        {
            localStart += nInternalToChangeThread[i];
        }

        forAll(internalToChangeLocal, i)
        {
            internalToChange[localStart++] = internalToChangeLocal[i];
        }

        localStart = 0;
        for (label i = 0; i < threadI; ++i)
        {
            localStart += nBoundaryToChangeThread[i];
        }

        forAll(boundaryToChangeLocal, i)
        {
            boundaryToChange[localStart++] = boundaryToChangeLocal[i];
        }

        # pragma omp barrier

        // start moving positions of faces
        # pragma omp for schedule(static)
        forAll(internalToChange, fI)
        {
            // swap with the face at the location the face should be
            face f;
            f.transfer(faces[internalToChange[fI]]);
            faces[internalToChange[fI]].transfer(faces[boundaryToChange[fI]]);
            faces[boundaryToChange[fI]].transfer(f);
            newFaceLabel[internalToChange[fI]] = boundaryToChange[fI];
            newFaceLabel[boundaryToChange[fI]] = internalToChange[fI];
        }

        # pragma omp barrier

        // renumber cells
        # pragma omp for schedule(dynamic, 40)
        forAll(cells, cellI)
        {
            cell& c = cells[cellI];

            forAll(c, fI)
            {
                if (newFaceLabel[c[fI]] != -1)
                {
                    c[fI] = newFaceLabel[c[fI]];
                }
            }
        }
    }

    // re-create boundary data
    PtrList<boundaryPatch>& boundaries = mesh_.boundaries_;
    if (boundaries.size() != 1)
    {
        boundaries.clear();
        boundaries.setSize(1);
        boundaries.set
        (
            0,
            new boundaryPatch
            (
                "defaultFaces",
                "patch",
                numBFaces,
                nInternalFaces
            )
        );
    }
    else
    {
        boundaries[0].patchStart() = nInternalFaces;
        boundaries[0].patchSize() = numBFaces;
    }

    // update face subsets
    mesh_.updateFaceSubsets(newFaceLabel);

    // delete invalid data
    mesh_.clearOut();
    this->clearOut();

    Info<< "Finished reordering boundary faces" << endl;
}


// ************************************************************************* //
