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
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceMapper.H"
#include "meshSurfacePartitioner.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "helperFunctions.H"

#include <map>

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceMapper::mapNodeToPatch
(
    const label bpI,
    const label patchI
)
{
    label patch, nt;
    point mapPoint;
    scalar dSq;

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const point p = points[bPoints[bpI]];

    if (patchI < 0)
    {
        meshOctree_.findNearestSurfacePoint(mapPoint, dSq, nt, patch, p);
    }
    else
    {
        meshOctree_.findNearestSurfacePointInRegion
        (
            mapPoint,
            dSq,
            nt,
            patchI,
            p
        );
    }

    meshSurfaceEngineModifier surfModifier(surfaceEngine_);
    surfModifier.moveBoundaryVertex(bpI, mapPoint);
}


void Foam::Module::meshSurfaceMapper::mapVerticesOntoSurface()
{
    Info<< "Mapping vertices onto surface" << endl;

    DynamicList<label> nodesToMap;
    nodesToMap.setSize(surfaceEngine_.boundaryPoints().size());
    forAll(nodesToMap, i)
    {
        nodesToMap[i] = i;
    }

    mapVerticesOntoSurface(nodesToMap);

    Info<< "Finished mapping vertices onto surface" << endl;
}


void Foam::Module::meshSurfaceMapper::mapVerticesOntoSurface
(
    const DynamicList<label>& nodesToMap
)
{
    const labelList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    const label size = nodesToMap.size();
    # pragma omp parallel for if (size > 1000) \
    schedule(dynamic, Foam::max(1, size/(3*omp_get_max_threads())))
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];

        label patch, nt;
        point mapPoint;
        scalar dSq;

        meshOctree_.findNearestSurfacePoint
        (
            mapPoint,
            dSq,
            nt,
            patch,
            points[boundaryPoints[bpI]]
        );

        surfaceModifier.moveBoundaryVertexNoUpdate(bpI, mapPoint);
    }

    // re-calculate face normals, point normals, etc.
    surfaceModifier.updateGeometry(nodesToMap);
}


void Foam::Module::meshSurfaceMapper::mapVerticesOntoSurfacePatches()
{
    Info<< "Mapping vertices with respect to surface patches" << endl;

    DynamicList<label> nodesToMap;
    nodesToMap.setSize(surfaceEngine_.boundaryPoints().size());
    forAll(nodesToMap, i)
    {
        nodesToMap[i] = i;
    }

    mapVerticesOntoSurfacePatches(nodesToMap);

    Info<< "Finished mapping vertices with respect to surface patches" << endl;
}


void Foam::Module::meshSurfaceMapper::mapVerticesOntoSurfacePatches
(
    const DynamicList<label>& nodesToMap
)
{
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPoints = mPart.corners();
    const labelHashSet& edgePoints = mPart.edgePoints();
    const labelListList2D& pointPatches = mPart.pointPatches();

    boolList treatedPoint(surfaceEngine_.boundaryPoints().size(), false);

    // find corner and edge points
    DynamicList<label> selectedCorners, selectedEdges;
    forAll(nodesToMap, i)
    {
        if (cornerPoints.found(nodesToMap[i]))
        {
            treatedPoint[nodesToMap[i]] = true;
            selectedCorners.append(nodesToMap[i]);
        }
        else if (edgePoints.found(nodesToMap[i]))
        {
            treatedPoint[nodesToMap[i]] = true;
            selectedEdges.append(nodesToMap[i]);
        }
    }

    // map the remaining selected points
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    const label size = nodesToMap.size();
    # pragma omp parallel for if (size > 1000) \
    schedule(dynamic, Foam::max(1, size /(3*omp_get_max_threads())))
    forAll(nodesToMap, nI)
    {
        const label bpI = nodesToMap[nI];

        if (treatedPoint[bpI])
        {
            continue;
        }

        const point& p = points[bPoints[bpI]];
        point mapPoint;
        scalar dSq;
        label nt;
        meshOctree_.findNearestSurfacePointInRegion
        (
            mapPoint,
            dSq,
            nt,
            pointPatches[bpI][0],
            p
        );

        surfaceModifier.moveBoundaryVertexNoUpdate(bpI, mapPoint);
    }

    // update face normals, point normals, etc.
    surfaceModifier.updateGeometry(nodesToMap);

    // map edge nodes
    mapEdgeNodes(selectedEdges);

    // map corner vertices
    mapCorners(selectedCorners);
}


// ************************************************************************* //
