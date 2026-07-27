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

#include "meshSurfaceCheckEdgeTypes.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "boolList.H"
#include "demandDrivenData.H"
#include "triangle.H"
#include "tetrahedron.H"

#include <map>

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceCheckEdgeTypes::classifyEdges()
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bp = surfaceEngine_.bp();
    const labelListList2D& pointFaces = surfaceEngine_.pointFaces();
    const edgeList& edges = surfaceEngine_.edges();
    const labelListList2D& edgeFaces = surfaceEngine_.edgeFaces();
    const labelList& facePatch = surfaceEngine_.boundaryFacePatches();
    const vectorField& fCentres = surfaceEngine_.faceCentres();

    // check if there exist tangled parts of mesh surface where
    // classification is not reliable
    boolList problematicPoint(pointFaces.size(), false);

    meshSurfacePartitioner mPart(surfaceEngine_);
    meshSurfaceCheckInvertedVertices checkInverted(mPart);
    const labelHashSet& invertedPoints = checkInverted.invertedVertices();
    for (const label pointi : invertedPoints.toc())
    {
        problematicPoint[bp[pointi]] = true;
    }

    // classify edges
    edgeType_.setSize(edges.size());

    label nThreads = 3*omp_get_num_procs();
    if (edges.size() < 1000)
        nThreads = 1;

    # pragma omp parallel num_threads(nThreads)
    {
        // TODO: this is not valid for non-manifold meshes
        // start checking feature edges
        # pragma omp for schedule(static, 1)
        forAll(edgeFaces, edgeI)
        {
            edgeType_[edgeI] = NONE;

            if (edgeFaces[edgeI].size() == 2)
            {
                const label f0 = edgeFaces[edgeI][0];
                const label f1 = edgeFaces[edgeI][1];

                if (facePatch[f0] == facePatch[f1])
                {
                    edgeType_[edgeI] |= PATCHEDGE;
                }
                else
                {
                    edgeType_[edgeI] |= FEATUREEDGE;
                }

                const edge e = edges[edgeI];

                // check if the surface is tangled there
                if (problematicPoint[bp[e.start()]])
                {
                    edgeType_[edgeI] |= UNDETERMINED;
                    continue;
                }

                if (problematicPoint[bp[e.end()]])
                {
                    edgeType_[edgeI] |= UNDETERMINED;
                    continue;
                }

                // check the volumes pof tets which can be formed at the edge
                const tetrahedron<point, point> tet0
                (
                    points[e.start()],
                    points[e.end()],
                    fCentres[f0],
                    fCentres[f1]
                );

                if (tet0.mag() > -VSMALL)
                {
                    edgeType_[edgeI] |= CONCAVEEDGE;
                    continue;
                }

                const tetrahedron<point, point> tet1
                (
                    points[e.end()],
                    points[e.start()],
                    fCentres[f1],
                    fCentres[f0]
                );

                if (tet1.mag() > -VSMALL)
                {
                    edgeType_[edgeI] |= CONCAVEEDGE;
                    continue;
                }

                edgeType_[edgeI] |= CONVEXEDGE;
            }
        }
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshSurfaceCheckEdgeTypes::meshSurfaceCheckEdgeTypes
(
    const meshSurfaceEngine& mse
)
:
    surfaceEngine_(mse),
    edgeType_()
{
    classifyEdges();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceCheckEdgeTypes::convexEdges
(
    DynamicList<label>& convexEdges
) const
{
    convexEdges.clear();

    forAll(edgeType_, eI)
    {
        if (edgeType_[eI] & CONVEXEDGE)
            convexEdges.append(eI);
    }
}


void Foam::Module::meshSurfaceCheckEdgeTypes::concaveEdges
(
    DynamicList<label>& concaveEdges
) const
{
    concaveEdges.clear();

    forAll(edgeType_, eI)
    {
        if (edgeType_[eI] & CONCAVEEDGE)
            concaveEdges.append(eI);
    }
}


// ************************************************************************* //
