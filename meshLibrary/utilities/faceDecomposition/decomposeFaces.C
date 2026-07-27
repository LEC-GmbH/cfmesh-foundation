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

#include "decomposeFaces.H"
#include "boolList.H"
#include "helperFunctionsGeometryQueries.H"

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

Foam::Module::decomposeFaces::decomposeFaces(polyMeshGen& mesh)
:
    mesh_(mesh),
    newFacesForFace_(),
    done_(false)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::decomposeFaces::decomposeMeshFaces
(
    const boolList& decomposeFace
)
{
    done_ = false;

    newFacesForFace_.setSize(mesh_.faces().size());
    forAll(newFacesForFace_, fI)
    {
        newFacesForFace_[fI].setSize(0);
    }

    const label nIntFaces = mesh_.nInternalFaces();
    label nFaces(0), nPoints = mesh_.points().size();
    point p;

    pointFieldPMG& points = mesh_.points();
    forAll(decomposeFace, fI)
    {
        if (decomposeFace[fI])
        {
            ++nFaces;
        }
    }

    points.setSize(nPoints + nFaces);

    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();

    if (decomposeFace.size() != faces.size())
    {
        FatalErrorInFunction
            << "Incorrect size of the decomposeFace list!" << abort(FatalError);
    }

    nFaces = 0;
    DynamicList<DynamicList<label>> newFaces;

    // decompose internal faces
    for (label faceI = 0; faceI < nIntFaces; ++faceI)
    {
        const face& f = faces[faceI];

        if (decomposeFace[faceI])
        {

            FixedList<label, 3> newF;

            forAll(f, pI)
            {
                newF[0] = f[pI];
                newF[1] = f.nextLabel(pI);
                newF[2] = nPoints;


                appendList(newFaces, newF);
                newFacesForFace_[faceI].append(nFaces++);
            }

            p = help::stableFaceCentre(f, points);
            points[nPoints] = p;
            ++nPoints;
        }
        else
        {

            appendList(newFaces, f);
            newFacesForFace_[faceI].append(nFaces++);
        }
    }

    // decompose boundary faces
    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        boundaries[patchI].patchStart() = nFaces;

        for (label bfI = start; bfI < end; ++bfI)
        {
            const face& f = faces[bfI];
            if (decomposeFace[bfI])
            {

                FixedList<label, 3> newF;

                forAll(f, pI)
                {
                    newF[0] = f[pI];
                    newF[1] = f.nextLabel(pI);
                    newF[2] = nPoints;


                    appendList(newFaces, newF);
                    newFacesForFace_[bfI].append(nFaces++);
                }

                p = help::stableFaceCentre(f, points);
                points[nPoints++] = p;
            }
            else
            {

                appendList(newFaces, f);
                newFacesForFace_[bfI].append(nFaces++);
            }
        }

        boundaries[patchI].patchSize() =
            nFaces - boundaries[patchI].patchStart();
    }

    // store the faces back into their list
    faces.setSize(nFaces);
    forAll(faces, faceI)
    {
        face& f = faces[faceI];
        f.setSize(newFaces[faceI].size());

        forAll(f, pI)
        {
            f[pI] = newFaces[faceI][pI];
        }
    }
    newFaces.setSize(0);

    // update subsets
    mesh_.updateFaceSubsets(newFacesForFace_);

    // change the mesh
    cellListPMG& cells = meshModifier.cellsAccess();

    # pragma omp parallel for schedule(dynamic, 40)
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        DynamicList<label> newC;

        forAll(c, fJ)
        {
            const label faceI = c[fJ];
            forAllRow(newFacesForFace_, faceI, nfI)
            newC.append(newFacesForFace_[faceI][nfI]);
        }


        c.setSize(newC.size());
        forAll(newC, fJ)
        c[fJ] = newC[fJ];
    }

    meshModifier.clearAll();

    done_ = true;

}


void Foam::Module::decomposeFaces::decomposeConcaveInternalFaces
(
    const boolList& concaveVertex
)
{
    done_ = false;

    newFacesForFace_.setSize(mesh_.faces().size());
    forAll(newFacesForFace_, fI)
        newFacesForFace_[fI].setSize(0);

    const label nIntFaces = mesh_.nInternalFaces();

    polyMeshGenModifier meshModifier(mesh_);
    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();

    if (concaveVertex.size() != mesh_.points().size())
    {
        FatalErrorInFunction
            << "Incorrect size of the concaveVertex list!" << abort(FatalError);
    }

    DynamicList<DynamicList<label>> newFaces;
    DynamicList<label> newF;
    newF.setSize(3);


    // decompose internal faces
    for (label faceI = 0; faceI < nIntFaces; ++faceI)
    {
        const face& f = faces[faceI];

        DynamicList<label> concavePos;
        forAll(f, pI)
        {
            if (concaveVertex[f[pI]])
            {
                concavePos.append(pI);
            }
        }

        if (concavePos.size() == 1)
        {

            newF[0] = f[concavePos[0]];

            for (label pI = 1; pI<(f.size()-1); ++pI)
            {
                const label pJ = (concavePos[0] + pI) % f.size();
                newF[1] = f[pJ];
                newF[2] = f.nextLabel(pJ);


                newFacesForFace_[faceI].append(newFaces.size());
                appendList(newFaces, newF);
            }
        }
        else if (concavePos.size() > 1)
        {

            newF[0] = points.size();
            forAll(f, pI)
            {
                newF[1] = f[pI];
                newF[2] = f.nextLabel(pI);


                newFacesForFace_[faceI].append(newFaces.size());
                appendList(newFaces, newF);
            }

            const point fCent = help::stableFaceCentre(f, points);
            points.append(fCent);
        }
        else
        {

            newFacesForFace_[faceI].append(newFaces.size());
            appendList(newFaces, f);
        }
    }

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        // set new patch start
        boundaries[patchI].patchStart() = newFaces.size();

        // store faces into newFaces
        for (label bfI = start; bfI < end; ++bfI)
        {
            newFacesForFace_[bfI].append(newFaces.size());
            appendList(newFaces, faces[bfI]);
        }
    }

    // copy new faces into the faceListPMG
    faces.setSize(newFaces.size());
    forAll(newFaces, faceI)
    {
        faces[faceI].setSize(newFaces[faceI].size());

        forAllRow(newFaces, faceI, pI)
        {
            faces[faceI][pI] = newFaces[faceI][pI];
        }
    }

    newFaces.setSize(0);

    // update cells
    cellListPMG& cells = meshModifier.cellsAccess();

    # pragma omp parallel for schedule(dynamic, 40)
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        DynamicList<label> newC;
        newC.reserve(24);

        forAll(c, fJ)
        {
            const label faceI = c[fJ];
            forAllRow(newFacesForFace_, faceI, nfI)
            newC.append(newFacesForFace_[faceI][nfI]);
        }


        c.setSize(newC.size());
        forAll(newC, fJ)
        {
            c[fJ] = newC[fJ];
        }
    }

    meshModifier.clearAll();

    // update subsets
    mesh_.updateFaceSubsets(newFacesForFace_);


    meshModifier.removeUnusedVertices();

    done_ = true;
}


const Foam::Module::labelListList2D&
Foam::Module::decomposeFaces::newFacesForFace() const
{
    if (!done_)
    {
        WarningInFunction
            << "Decomposition is not yet performed!" << endl;
    }

    return newFacesForFace_;
}


// ************************************************************************* //
