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

#include "meshOctree.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "meshSurfaceMapper2D.H"
#include "meshSurfaceEngine.H"
#include "polyMeshGen2DEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfacePartitioner.H"
#include "labelledScalar.H"
#include "polyMeshGenAddressing.H"
#include "boundBox.H"


#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceMapper2D::findMappingDistance
(
    const DynamicList<label>& edgesToMap,
    Map<scalar>& mappingDistance
) const
{
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const labelListList2D& eFaces = surfaceEngine_.edgeFaces();
    const edgeList& edges = surfaceEngine_.edges();
    const pointFieldPMG& points = surfaceEngine_.points();

    // generate search distance for corner edges
    mappingDistance.clear();
    forAll(edgesToMap, i)
    {
        const label beI = edgesToMap[i];

        mappingDistance.insert(beI, 0.0);
        Map<scalar>::iterator mIter = mappingDistance.find(beI);

        const point p = edges[beI].centre(points);
        forAllRow(eFaces, beI, efI)
        {
            const scalar d = magSqr(faceCentres[eFaces[beI][efI]] - p);
            mIter() = Foam::max(mIter(), d);
        }

        // safety factor
        mIter() *= 16.0;
    }

}


void Foam::Module::meshSurfaceMapper2D::adjustZCoordinates()
{
    // create the bounding box of the surface mesh
    const boundBox sbb(meshOctree_.surface().points());

    const labelList& bp = surfaceEngine_.bp();
    const pointFieldPMG& points = surfaceEngine_.points();

    meshSurfaceEngineModifier modifier(surfaceEngine_);

    const polyMeshGen2DEngine& mesh2DEngine = this->mesh2DEngine();
    const boolList& minZPoint = mesh2DEngine.zMinPoints();
    const boolList& maxZPoint = mesh2DEngine.zMaxPoints();

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(minZPoint, pointI)
    {
        point newP = points[pointI];

        if (minZPoint[pointI])
        {
            newP.z() = sbb.min().z();
        }
        else if (maxZPoint[pointI])
        {
            newP.z() = sbb.max().z();
        }
        else
        {
            FatalErrorInFunction
                << "This mesh is not in the x - y plane!" << exit(FatalError);
        }

        modifier.moveBoundaryVertexNoUpdate(bp[pointI], newP);
    }

    modifier.updateGeometry();
}


void Foam::Module::meshSurfaceMapper2D::mapVerticesOntoSurface()
{
    DynamicList<label> edgesToMap;

    forAll(activeBoundaryEdges_, beI)
    {
        edgesToMap.append(activeBoundaryEdges_[beI]);
    }

    mapVerticesOntoSurface(edgesToMap);
}


void Foam::Module::meshSurfaceMapper2D::mapVerticesOntoSurface
(
    const DynamicList<label>& edgesToMap
)
{
    const edgeList& edges = surfaceEngine_.edges();
    const labelList& bp = surfaceEngine_.bp();
    const pointFieldPMG& points = surfaceEngine_.points();

    DynamicList<label> nodesToMap;
    forAll(edgesToMap, i)
    {
        const edge& e = edges[edgesToMap[i]];
        nodesToMap.append(bp[e.start()]);
        nodesToMap.append(bp[e.end()]);
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    const label size = edgesToMap.size();
    # pragma omp parallel for if (size > 1000) \
    schedule(dynamic, Foam::max(1, size/(3*omp_get_max_threads())))
    forAll(edgesToMap, i)
    {
        const label beI = edgesToMap[i];
        const edge& e = edges[beI];

        label patch, nt;
        point mapPoint;
        scalar dSq;

        meshOctree_.findNearestSurfacePoint
        (
            mapPoint,
            dSq,
            nt,
            patch,
            points[e.start()]
        );

        mapPoint.z() = points[e.start()].z();
        surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.start()], mapPoint);

        mapPoint.z() = points[e.end()].z();
        surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.end()], mapPoint);
    }

    surfaceModifier.updateGeometry(nodesToMap);
}


void Foam::Module::meshSurfaceMapper2D::mapCorners()
{
    const edgeList& edges = surfaceEngine_.edges();

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPoints = mPart.corners();

    DynamicList<label> selectedEdges;
    forAll(activeBoundaryEdges_, eI)
    {
        const edge& e = edges[activeBoundaryEdges_[eI]];

        if (cornerPoints.found(e.start()) || cornerPoints.found(e.end()))
        {
            selectedEdges.append(activeBoundaryEdges_[eI]);
        }
    }

    mapCorners(selectedEdges);
}


void Foam::Module::meshSurfaceMapper2D::mapCorners
(
    const DynamicList<label>& edgesToMap
)
{
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& corners = mPart.corners();
    const labelListList2D& pPatches = mPart.pointPatches();

    const pointFieldPMG& points = surfaceEngine_.points();
    const edgeList& edges = surfaceEngine_.edges();
    const labelList& bp = surfaceEngine_.bp();

    Map<scalar> mappingDistance;
    findMappingDistance(edgesToMap, mappingDistance);

    // for every corner in the mesh surface find the nearest corner in the
    // triSurface
    meshSurfaceEngineModifier sMod(surfaceEngine_);

    # pragma omp parallel for if (edgesToMap.size() > 10)
    forAll(edgesToMap, eI)
    {
        const label beI = edgesToMap[eI];
        const edge& be = edges[beI];

        const label bps = bp[be.start()];
        const label bpe = bp[be.end()];

        DynamicList<label> ePatches;
        forAllRow(pPatches, bps, i)
        {
            if (found(pPatches, bpe, pPatches[bps][i]))
            {
                ePatches.append(pPatches[bps][i]);
            }
        }

        if (!corners.found(bps) || !corners.found(bpe))
        {
            FatalErrorInFunction
                << "Trying to map a point that is not a corner"
                << abort(FatalError);
        }

        const point& p = points[be.start()];
        const scalar maxDist = mappingDistance[beI];

        // find the nearest position to the given point patches
        point mapPointApprox(p);
        scalar distSqApprox;

        label iter(0);
        while (iter++ < 40)
        {
            point newP(vector::zero);
            forAll(ePatches, epI)
            {
                point np;
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    distSqApprox,
                    nt,
                    ePatches[epI],
                    mapPointApprox
                );

                newP += np;
            }

            newP /= ePatches.size();
            if (magSqr(newP - mapPointApprox) < 1e-8*maxDist)
            {
                break;
            }

            mapPointApprox = newP;
        }

        distSqApprox = magSqr(mapPointApprox - p);

        // find the nearest triSurface corner for the given corner
        scalar distSq;
        point mapPoint;
        label nse;
        meshOctree_.findNearestEdgePoint(mapPoint, distSq, nse, p, ePatches);

        if (distSq > mappingDistance[beI])
        {
            const vector disp = mapPoint - p;

            mapPoint = Foam::sqrt(mappingDistance[beI] / distSq)*disp + p;
            distSq = mappingDistance[beI];
        }

        if (distSq > 1.2*distSqApprox)
        {
            mapPoint = mapPointApprox;
        }

        // move the point to the nearest corner
        mapPoint.z() = points[be.start()].z();
        sMod.moveBoundaryVertexNoUpdate(bps, mapPoint);

        mapPoint.z() = points[be.end()].z();
        sMod.moveBoundaryVertexNoUpdate(bpe, mapPoint);
    }

    DynamicList<label> nodesToMap;
    forAll(edgesToMap, eI)
    {
        const edge& be = edges[edgesToMap[eI]];

        nodesToMap.append(bp[be.start()]);
        nodesToMap.append(bp[be.end()]);
    }

    sMod.updateGeometry(nodesToMap);
}


void Foam::Module::meshSurfaceMapper2D::mapVerticesOntoSurfacePatches()
{
    DynamicList<label> edgesToMap;

    forAll(activeBoundaryEdges_, beI)
    {
        edgesToMap.append(activeBoundaryEdges_[beI]);
    }

    mapVerticesOntoSurfacePatches(edgesToMap);
}


void Foam::Module::meshSurfaceMapper2D::mapVerticesOntoSurfacePatches
(
    const DynamicList<label>& edgesToMap
)
{
    // map the selected edges on their patches
    const edgeList& edges = surfaceEngine_.edges();
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelListList2D& edgePatches = surfaceEngine_.edgePatches();
    const labelList& bp = surfaceEngine_.bp();

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    DynamicList<label> selectedCorners;

    const label size = edgesToMap.size();
    # pragma omp parallel for if (size > 1000) \
    schedule(dynamic, Foam::max(50, size/(3*omp_get_max_threads())))
    forAll(edgesToMap, nI)
    {
        const label beI = edgesToMap[nI];

        if (edgePatches[beI].size() == 2)
        {
            # pragma omp critical
            selectedCorners.append(beI);

            continue;
        }

        const edge& e = edges[beI];

        const point& p = points[e.start()];
        point mapPoint;
        scalar dSq;
        label nt;

        meshOctree_.findNearestSurfacePointInRegion
        (
            mapPoint,
            dSq,
            nt,
            edgePatches[beI][0],
            p
        );

        mapPoint.z() = points[e.start()].z();
        surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.start()], mapPoint);

        mapPoint.z() = points[e.end()].z();
        surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.end()], mapPoint);
    }

    mapCorners(selectedCorners);

    // update geometry at moved vertices
    selectedCorners.clear();
    forAll(edgesToMap, eI)
    {
        const edge& e = edges[edgesToMap[eI]];

        selectedCorners.append(bp[e.start()]);
        selectedCorners.append(bp[e.end()]);
    }

    surfaceModifier.updateGeometry(selectedCorners);
}


// ************************************************************************* //
