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

#include "checkBoundaryFacesSharingTwoEdges.H"
#include "demandDrivenData.H"
#include "meshSurfaceEngine.H"
#include "decomposeFaces.H"
#include "decomposeCells.H"

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void Foam::Module::checkBoundaryFacesSharingTwoEdges::createMeshSurface() const
{
    meshSurfacePtr_ = new meshSurfaceEngine(mesh_);
}


void Foam::Module::checkBoundaryFacesSharingTwoEdges::findFacesAtBndEdge()
{
    const meshSurfaceEngine& mse = meshSurface();

    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const labelListList2D& pointEdges = mse.boundaryPointEdges();

    const label nIntFaces = mesh_.nInternalFaces();
    const faceListPMG& faces = mesh_.faces();

    // find the internal faces attached to the boundary points
    removeBndPoint_.setSize(pointEdges.size());
    removeBndPoint_ = true;

    # pragma omp parallel for if (nIntFaces > 100) schedule(dynamic, 20)
    for (label fI = 0; fI < nIntFaces; ++fI)
    {
        const face& f = faces[fI];

        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if (bpI < 0)
            {
                continue;
            }

            if (nBndFacesAtBndPoint_[bpI] == 2)
            {
                const edge ePrev = f.faceEdge(f.rcIndex(pI));
                const edge eNext = f.faceEdge(pI);

                bool foundNext(false), foundPrev(false);
                forAllRow(pointEdges, bpI, peI)
                {
                    const label beI = pointEdges[bpI][peI];

                    if (edges[beI] == ePrev)
                    {
                        foundPrev = true;
                    }
                    else if (edges[beI] == eNext)
                    {
                        foundNext = true;
                    }
                }

                if (!(foundPrev && foundNext))
                {
                    removeBndPoint_[bpI] = false;
                }
            }
            else
            {
                removeBndPoint_[bpI] = false;
            }
        }
    }

}


void Foam::Module::checkBoundaryFacesSharingTwoEdges::findBndFacesAtBndVertex()
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelListList2D& pointFaces = mse.pointFaces();

    nBndFacesAtBndPoint_.setSize(pointFaces.size());
    nBndFacesAtBndPoint_ = 0;

    forAll(nBndFacesAtBndPoint_, bpI)
    {
        nBndFacesAtBndPoint_[bpI] = pointFaces[bpI].size();
    }
}


void Foam::Module::checkBoundaryFacesSharingTwoEdges::removeExcessiveVertices()
{
    const labelList& bp = meshSurface().bp();
    const faceListPMG& faces = mesh_.faces();

    // remove points which can be safely be removed
    // internal faces
    const label nIntFaces = mesh_.nInternalFaces();

    # pragma omp parallel for if (nIntFaces > 100) schedule(dynamic, 10)
    for (label faceI = 0; faceI < nIntFaces; ++faceI)
    {
        const face& f = faces[faceI];

        DynamicList<label> newF;
        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if (
                (bpI >= 0) && removeBndPoint_[bpI] &&
                (nBndFacesAtBndPoint_[bpI] == 2)
            )
                continue;

            newF.append(f[pI]);
        }

        if (newF.size() < f.size())
        {
            face& mf = const_cast<face&>(f);
            mf.setSize(newF.size());
            forAll(mf, i)
            {
                mf[i] = newF[i];
            }
        }
    }

    // boundary faces
    forAll(mesh_.boundaries(), patchI)
    {
        const label start = mesh_.boundaries()[patchI].patchStart();
        const label end = start + mesh_.boundaries()[patchI].patchSize();

        # pragma omp parallel for if (end - start > 100) \
        schedule(dynamic, 10)
        for (label faceI = start; faceI < end; ++faceI)
        {
            const face& f = faces[faceI];

            DynamicList<label> newF;
            forAll(f, pI)
            {
                const label bpI = bp[f[pI]];

                if (removeBndPoint_[bpI] && (nBndFacesAtBndPoint_[bpI] == 2))
                {
                    continue;
                }

                newF.append(f[pI]);
            }

            if (newF.size() < f.size())
            {
                face& mf = const_cast<face&>(f);
                mf.setSize(newF.size());
                forAll(mf, i)
                {
                    mf[i] = newF[i];
                }
            }
        }
    }
}


Foam::label
Foam::Module::checkBoundaryFacesSharingTwoEdges::findBndFacesForDecomposition
(
    boolList& decomposeFace
)
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const faceList::subList& bFaces = mse.boundaryFaces();

    label nDecomposed(0);
    const label nIntFaces = mesh_.nInternalFaces();

    # pragma omp parallel for if (bFaces.size() > 100) \
    schedule(dynamic, 10) reduction(+ : nDecomposed)
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
        {
            const label bpI = bp[bf[pI]];

            if (nBndFacesAtBndPoint_[bpI] == 2)
            {
                ++nDecomposed;
                decomposeFace[nIntFaces + bfI] = true;
            }
        }
    }

    return nDecomposed;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::Module::checkBoundaryFacesSharingTwoEdges::
checkBoundaryFacesSharingTwoEdges
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    meshSurfacePtr_(nullptr),
    nBndFacesAtBndPoint_(),
    removeBndPoint_()
{}


// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

Foam::Module::checkBoundaryFacesSharingTwoEdges::
~checkBoundaryFacesSharingTwoEdges()
{
    deleteDemandDrivenData(meshSurfacePtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::checkBoundaryFacesSharingTwoEdges::findPoints
(
    labelHashSet& badPoints
)
{
    badPoints.clear();

    findBndFacesAtBndVertex();

    const labelList& bPoints = meshSurface().boundaryPoints();
    forAll(nBndFacesAtBndPoint_, bpI)
    {
        if (nBndFacesAtBndPoint_[bpI] != 2)
        {
            continue;
        }

        badPoints.insert(bPoints[bpI]);
    }
}


bool Foam::Module::checkBoundaryFacesSharingTwoEdges::improveTopology()
{
    bool changed(false);

    findBndFacesAtBndVertex();

    findFacesAtBndEdge();

    removeExcessiveVertices();

    boolList decomposeFace(mesh_.faces().size(), false);
    const label nDecomposed = findBndFacesForDecomposition(decomposeFace);

    Info<< "Marked " << nDecomposed << " faces for decomposition" << endl;

    if (nDecomposed != 0)
    {
        // delete the mesh surface engine
        deleteDemandDrivenData(meshSurfacePtr_);

        // find cells which will be decomposed
        boolList decomposeCell(mesh_.cells().size(), false);
        const labelList& owner = mesh_.owner();
        forAll(decomposeFace, faceI)
        {
            if (decomposeFace[faceI])
            {
                decomposeCell[owner[faceI]];
            }
        }

        // decompose marked faces
        decomposeFaces(mesh_).decomposeMeshFaces(decomposeFace);

        // decompose cells
        decomposeCells dc(mesh_);
        dc.decomposeMesh(decomposeCell);

        changed = true;
    }

    polyMeshGenModifier(mesh_).removeUnusedVertices();

    return changed;
}


// ************************************************************************* //
