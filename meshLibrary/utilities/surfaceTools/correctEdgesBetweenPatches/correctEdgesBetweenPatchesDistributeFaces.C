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

#include "demandDrivenData.H"
#include "correctEdgesBetweenPatches.H"
#include "decomposeFaces.H"
#include "triSurface.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceCheckEdgeTypes.H"
#include "meshSurfacePartitioner.H"
#include "helperFunctions.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::correctEdgesBetweenPatches::decomposeConcaveFaces()
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const labelListList2D& edgeFaces = mse.edgeFaces();
    const labelListList2D& bpEdges = mse.boundaryPointEdges();
    const labelList& facePatch = mse.boundaryFacePatches();

    // classify edges at the surface
    meshSurfaceCheckEdgeTypes edgeChecker(mse);
    const List<direction>& edgeType = edgeChecker.edgeTypes();

    // find concave points
    boolList concavePoint(bPoints.size(), false);

    labelList edgeInPatch(edges.size(), -1);

    direction problematicTypes = 0;
    problematicTypes |= meshSurfaceCheckEdgeTypes::CONCAVEEDGE;
    problematicTypes |= meshSurfaceCheckEdgeTypes::UNDETERMINED;

    # pragma omp parallel for schedule(dynamic, 100)
    forAll(edgeType, eI)
    {
        if (edgeType[eI] & meshSurfaceCheckEdgeTypes::PATCHEDGE)
        {
            if (edgeFaces[eI].size())
            {
                edgeInPatch[eI] = facePatch[edgeFaces[eI][0]];
            }

            continue;
        }

        if (edgeType[eI] & problematicTypes)
        {
            const edge& e = edges[eI];

            concavePoint[bp[e.start()]] = true;
            concavePoint[bp[e.end()]] = true;
        }
    }

    // decompose internal faces attached to concave vertices which have two
    // or more edges at the boundary
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    boolList decomposeFace(faces.size(), false);
    label nDecomposed(0);

    # pragma omp parallel for schedule(dynamic, 100) \
    reduction(+ : nDecomposed)
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); ++faceI)
    {
        const face& f = faces[faceI];

        bool hasConcave(false);
        label nBndEdges(0);
        DynamicList<label> bndEdgePatches;

        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if (bpI < 0)
            {
                continue;
            }

            if (concavePoint[bpI])
            {
                hasConcave = true;
            }

            // points is at a concave edge
            // count the number of boundary edge
            const edge e = f.faceEdge(pI);

            forAllRow(bpEdges, bpI, bpeI)
            {
                const label beI = bpEdges[bpI][bpeI];
                const edge& ee = edges[beI];

                if (e == ee)
                {
                    ++nBndEdges;
                    appendUniq(bndEdgePatches, edgeInPatch[beI]);
                    break;
                }
            }
        }

        if (hasConcave && (nBndEdges > 1) && (bndEdgePatches.size() > 1))
        {
            // the face has two or more edges at the boundary
            // Hence, it is marked for decomposition
            decomposeFace[faceI] = true;

            decomposeCell_[owner[faceI]] = true;
            decomposeCell_[neighbour[faceI]] = true;

            ++nDecomposed;
        }
    }

    // finally, perform decomposition of marked faces
    if (nDecomposed != 0)
    {
        Info<< "Decomposing " << nDecomposed << " internal faces" << endl;
        decomposeFaces(mesh_).decomposeMeshFaces(decomposeFace);

        decompose_ = true;

        clearMeshSurface();
        mesh_.clearAddressingData();
    }
}


void Foam::Module::correctEdgesBetweenPatches::patchCorrection()
{
    Info<< "Performing patch correction" << endl;

    const meshSurfaceEngine& mse = meshSurface();

    meshSurfacePartitioner surfacePartitioner(mse);

    const labelList& bPoints = mse.boundaryPoints();
    const labelListList2D& faceEdges = mse.faceEdges();
    const labelListList2D& edgeFaces = mse.edgeFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const labelList& boundaryFaceOwners = mse.faceOwners();
    const labelList& facePatches = mse.boundaryFacePatches();

    // set flag 1 to corner vertices, flag 2 to edge vertices
    List<direction> nodeType(bPoints.size(), direction(0));

    // set corner flags
    const labelHashSet& corners = surfacePartitioner.corners();
    for (const label corneri : corners.toc())
    {
        nodeType[corneri] |= 1;
    }

    // set flgs to edge vertices
    const labelHashSet& edgePoints = surfacePartitioner.edgePoints();
    for (const label pointi : edgePoints.toc())
    {
        nodeType[pointi] |= 2;
    }

    // set flags for feature edges
    boolList featureEdge(edgeFaces.size(), false);

    # pragma omp parallel for schedule(guided)
    forAll(edgeFaces, eI)
    {
        if (edgeFaces[eI].size() != 2)
        {
            continue;
        }

        if (facePatches[edgeFaces[eI][0]] != facePatches[edgeFaces[eI][1]])
        {
            featureEdge[eI] = true;
        }
    }

    // decompose bad faces into triangles
    newBoundaryFaces_.clear();
    newBoundaryOwners_.clear();
    newBoundaryPatches_.clear();
    face triF(3);

    label nDecomposedFaces(0);
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        bool store(true);

        forAll(bf, i)
        {
            if
            (
                (nodeType[bp[bf[i]]] == direction(2))
             && featureEdge[faceEdges[bfI][i]]
             && featureEdge[faceEdges[bfI][bf.rcIndex(i)]]
            )
            {
                if (bf.size() > 4)
                {
                    store = false;
                    ++nDecomposedFaces;
                    decomposeCell_[boundaryFaceOwners[bfI]] = true;
                    decompose_ = true;

                    // decompose into triangles
                    const point p =
                        help::stableFaceCentre(bf, mesh_.points());
                    triF[2] = mesh_.points().size();
                    mesh_.points().append(p);

                    forAll(bf, j)
                    {
                        triF[0] = bf[j];
                        triF[1] = bf.nextLabel(j);

                        appendList(newBoundaryFaces_, triF);
                        newBoundaryOwners_.append(boundaryFaceOwners[bfI]);
                        newBoundaryPatches_.append(facePatches[bfI]);
                    }

                    break;
                }
                else if (bf.size() == 4)
                {
                    store = false;
                    ++nDecomposedFaces;
                    decomposeCell_[boundaryFaceOwners[bfI]] = true;
                    decompose_ = true;

                    // decompose the quad into 2 triangles
                    triF[0] = bf[i];

                    triF[1] = bf.nextLabel(i);
                    triF[2] = bf[(i + 2)%4];
                    appendList(newBoundaryFaces_, triF);
                    newBoundaryOwners_.append(boundaryFaceOwners[bfI]);
                    newBoundaryPatches_.append(facePatches[bfI]);

                    triF[1] = bf[(i + 2)%4];
                    triF[2] = bf.prevLabel(i);
                    appendList(newBoundaryFaces_, triF);
                    newBoundaryOwners_.append(boundaryFaceOwners[bfI]);
                    newBoundaryPatches_.append(facePatches[bfI]);

                    break;
                }
            }
        }

        if (store)
        {
            // face has not been altered
            appendList(newBoundaryFaces_, bf);
            newBoundaryOwners_.append(boundaryFaceOwners[bfI]);
            newBoundaryPatches_.append(facePatches[bfI]);
        }
    }

    if (nDecomposedFaces != 0)
    {
        replaceBoundary();
        clearMeshSurface();
        mesh_.clearAddressingData();
    }

    Info<< "Finished with patch correction" << endl;
}


// ************************************************************************* //
