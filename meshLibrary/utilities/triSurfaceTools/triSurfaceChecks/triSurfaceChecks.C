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

#include "triSurfaceChecks.H"
#include "triSurf.H"
#include "boundBox.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::Module::triSurfaceChecks::checkAngles
(
    const triSurf& surf,
    DynamicList<label>& badTriangles,
    const scalar angleTol
)
{
    badTriangles.clear();

    const scalar tol = Foam::cos(angleTol*M_PI/180.0);

    const pointField& pts = surf.points();

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(surf, tI)
    {
        const labelledTri& tri = surf[tI];

        forAll(tri, pI)
        {
            vector vn = pts[tri[tri.fcIndex(pI)]] - pts[tri[pI]];
            vn /= Foam::max(mag(vn), VSMALL);
            vector vp = pts[tri[pI]] - pts[tri[tri.rcIndex(pI)]];
            vp /= Foam::max(mag(vp), VSMALL);

            if ((vp & vn) >= tol)
            {
                # pragma omp critical
                badTriangles.append(tI);

                break;
            }
        }
    }

    return badTriangles.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkAngles
(
    triSurf& surf,
    const word subsetName,
    const scalar angleTol
)
{
    DynamicList<label> badTriangles;

    if (checkAngles(surf, badTriangles, angleTol))
    {
        label setId = surf.facetSubsetIndex(subsetName);
        if (setId >= 0)
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetName);

        forAll(badTriangles, i)
            surf.addFacetToSubset(setId, badTriangles[i]);
    }

    return badTriangles.size();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
namespace manifoldOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class surfaceNeiOp
{
    // Private data

        // const reference to face-edges addressing
        const labelListList2D& faceEdges_;

        // const reference to edge-faces addressing
        const labelListList2D& edgeFaces_;


public:

    // Constructors

        surfaceNeiOp
        (
            const labelListList2D& faceEdges,
            const labelListList2D& edgeFaces
        )
        :
            faceEdges_(faceEdges),
            edgeFaces_(edgeFaces)
        {}


    // Public member functions

        label size() const
        {
            return faceEdges_.size();
        }

        void operator()(const label tI, DynamicList<label>& neiTriangles) const
        {
            neiTriangles.clear();

            forAllRow(faceEdges_, tI, teI)
            {
                const label eI = faceEdges_[tI][teI];

                if (edgeFaces_[eI].size() != 2)
                    continue;

                // ind the neighbour triangle over the edge
                label nei = edgeFaces_[eI][0];

                if (nei == tI)
                    nei = edgeFaces_[eI][1];

                neiTriangles.append(nei);
            }
        }

        template<class labelListType>
        void collectGroups
        (
            std::map<label, DynamicList<label>>& /*neiGroups*/,
            const labelListType& /*elementInGroup*/,
            const DynamicList<label>& /*localGroupLabel*/
        ) const
        {
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace manifoldOps
} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
namespace selectorOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class selectOp
{
public:

    selectOp()
    {}


    bool operator()(const label /*tI*/) const
    {
        return true;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace selectorOps
} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


Foam::label Foam::Module::triSurfaceChecks::checkSurfaceManifolds
(
    const triSurf& surf,
    DynamicList<label>& triangleInManifold
)
{
    const labelListList2D& faceEdges = surf.facetEdges();
    const labelListList2D& edgeFaces = surf.edgeFacets();

    const label nManifolds =
        help::groupMarking
        (
            triangleInManifold,
            manifoldOps::surfaceNeiOp(faceEdges, edgeFaces),
            selectorOps::selectOp()
        );

    return nManifolds;
}


Foam::label Foam::Module::triSurfaceChecks::checkSurfaceManifolds
(
    triSurf& surf,
    const word subsetPrefix
)
{
    DynamicList<label> facetInManifold;

    const label nManifolds = checkSurfaceManifolds(surf, facetInManifold);

    if (nManifolds > 1)
    {
        labelList groupIds(nManifolds);
        forAll(groupIds, i)
        {
            const word sName = subsetPrefix + Foam::name(i);
            label setId = surf.facetSubsetIndex(sName);
            if (setId >= 0)
                surf.removeFacetSubset(setId);
            groupIds[i] = surf.addFacetSubset(sName);
        }

        forAll(facetInManifold, tI)
            surf.addFacetToSubset(groupIds[facetInManifold[tI]], tI);
    }

    return nManifolds;
}


Foam::label Foam::Module::triSurfaceChecks::checkForHoles
(
    const triSurf& surf,
    DynamicList<label>& badTriangles
)
{
    badTriangles.clear();

    const labelListList2D& edgeFacets = surf.edgeFacets();

    # pragma omp parallel for schedule(dynamic, 100)
    forAll(edgeFacets, eI)
    {
        if (edgeFacets[eI].size() == 1)
        {
            # pragma omp critical
            badTriangles.append(edgeFacets[eI][0]);
        }
    }

    return badTriangles.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkForHoles
(
    triSurf& surf,
    const word subsetName
)
{
    DynamicList<label> trianglesNearHoles;

    if (checkForHoles(surf, trianglesNearHoles))
    {
        label setId = surf.facetSubsetIndex(subsetName);
        if (setId >= 0)
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetName);

        forAll(trianglesNearHoles, i)
            surf.addFacetToSubset(setId, trianglesNearHoles[i]);
    }

    return trianglesNearHoles.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkForNonManifoldEdges
(
    const triSurf& surf,
    DynamicList<label>& badTriangles
)
{
    badTriangles.clear();

    const labelListList2D& edgeFacets = surf.edgeFacets();

    # pragma omp parallel for schedule(dynamic, 100)
    forAll(edgeFacets, eI)
    {
        if (edgeFacets[eI].size() > 2)
        {
            # pragma omp critical
            {
                forAllRow(edgeFacets, eI, efI)
                    badTriangles.append(edgeFacets[eI][efI]);
            }
        }
    }

    return badTriangles.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkForNonManifoldEdges
(
    triSurf& surf,
    const word subsetPrefix
)
{
    DynamicList<label> trianglesNearHoles;

    if (checkForHoles(surf, trianglesNearHoles))
    {
        label setId = surf.facetSubsetIndex(subsetPrefix);
        if (setId >= 0)
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetPrefix);

        forAll(trianglesNearHoles, i)
            surf.addFacetToSubset(setId, trianglesNearHoles[i]);
    }

    return trianglesNearHoles.size();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
namespace orientationOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class surfaceNeiOp
{
    // Private data
        // const reference to face-edges addressing
        const labelListList2D& faceEdges_;

        // const reference to edge-faces addressing
        const labelListList2D& edgeFaces_;

        // const reference to triangles
        const DynamicList<labelledTri>& triangles_;

public:

    // Constructors

        surfaceNeiOp
        (
            const labelListList2D& faceEdges,
            const labelListList2D& edgeFaces,
            const DynamicList<labelledTri>& triangles
        )
        :
            faceEdges_(faceEdges),
            edgeFaces_(edgeFaces),
            triangles_(triangles)
        {}


    // Public member functions
        label size() const
        {
            return triangles_.size();
        }

        void operator()(const label tI, DynamicList<label>& neiTriangles) const
        {
            neiTriangles.clear();

            const labelledTri& tri = triangles_[tI];

            forAllRow(faceEdges_, tI, teI)
            {
                const label eI = faceEdges_[tI][teI];

                if (edgeFaces_[eI].size() != 2)
                    continue;

                // ind the neighbour triangle over the edge
                label nei = edgeFaces_[eI][0];

                if (nei == tI)
                    nei = edgeFaces_[eI][1];

                const labelledTri& neiTri = triangles_[nei];

                // check the orientation
                label pos(-1);
                forAll(neiTri, i)
                    if (neiTri[i] == tri[teI])
                    {
                        pos = i;
                        break;
                    }

                if (tri[tri.fcIndex(teI)] == neiTri[neiTri.rcIndex(pos)])
                {
                    // triangles are of the same orientation
                    neiTriangles.append(nei);
                }
            }
        }

        template<class labelListType>
        void collectGroups
        (
            std::map<label, DynamicList<label>>& /*neiGroups*/,
            const labelListType& /*elementInGroup*/,
            const DynamicList<label>& /*localGroupLabel*/
        ) const
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace orientationOps
} // End namespace Module
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::Module::triSurfaceChecks::checkOrientation
(
    const triSurf& surf,
    DynamicList<label>& triangleInGroup
)
{
    triangleInGroup.setSize(surf.size());
    triangleInGroup = -1;
    label nGroups(0);

    const labelListList2D& edgeFacets = surf.edgeFacets();
    const labelListList2D& faceEdges = surf.facetEdges();
    const DynamicList<labelledTri>& triangles = surf.facets();

    nGroups =
        help::groupMarking
        (
            triangleInGroup,
            orientationOps::surfaceNeiOp(faceEdges, edgeFacets, triangles),
            selectorOps::selectOp()
        );

    return nGroups;
}


Foam::label Foam::Module::triSurfaceChecks::checkOrientation
(
    triSurf& surf,
    const word subsetPrefix
)
{
    DynamicList<label> triangleInGroup;

    const label nGroups = checkOrientation(surf, triangleInGroup);

    if (nGroups > 1)
    {
        labelList groupIds(nGroups);
        forAll(groupIds, i)
        {
            const word sName = subsetPrefix + Foam::name(i);
            label setId = surf.facetSubsetIndex(sName);
            if (setId >= 0)
                surf.removeFacetSubset(setId);

            groupIds[i] = surf.addFacetSubset(sName);
        }

        forAll(triangleInGroup, tI)
            surf.addFacetToSubset(groupIds[triangleInGroup[tI]], tI);
    }

    return nGroups;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
namespace connectionOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class surfaceNeiOp
{
    // Private data

        // const reference to face-edges addressing
        const labelListList2D& faceEdges_;

        // const reference to edge-faces addressing
        const labelListList2D& edgeFaces_;


public:

    // Constructors

        surfaceNeiOp
        (
            const labelListList2D& faceEdges,
            const labelListList2D& edgeFaces
        )
        :
            faceEdges_(faceEdges),
            edgeFaces_(edgeFaces)
        {}


    // Public member functions

        label size() const
        {
            return faceEdges_.size();
        }

        void operator()(const label tI, DynamicList<label>& neiTriangles) const
        {
            neiTriangles.clear();

            forAllRow(faceEdges_, tI, teI)
            {
                const label eI = faceEdges_[tI][teI];

                forAllRow(edgeFaces_, eI, efI)
                {
                    const label tJ = edgeFaces_[eI][efI];

                    if (tJ == tI)
                        continue;

                    neiTriangles.append(tJ);
                }
            }
        }

        template<class labelListType>
        void collectGroups
        (
            std::map<label, DynamicList<label>>& /*neiGroups*/,
            const labelListType& /*elementInGroup*/,
            const DynamicList<label>& /*localGroupLabel*/
        ) const
        {
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace orientationOps
} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::Module::triSurfaceChecks::checkDisconnectedParts
(
    const triSurf& surf,
    DynamicList<label>& triangleInRegion
)
{
    triangleInRegion.setSize(surf.size());
    triangleInRegion = -1;
    label nGroups(0);

    const labelListList2D& edgeFacets = surf.edgeFacets();
    const labelListList2D& faceEdges = surf.facetEdges();

    nGroups =
        help::groupMarking
        (
            triangleInRegion,
            connectionOps::surfaceNeiOp(faceEdges, edgeFacets),
            selectorOps::selectOp()
        );

    return nGroups;
}


Foam::label Foam::Module::triSurfaceChecks::checkDisconnectedParts
(
    triSurf& surf,
    const word subsetPrefix
)
{
    DynamicList<label> triangleInRegion;

    const label nGroups = checkDisconnectedParts(surf, triangleInRegion);

    if (nGroups > 1)
    {
        labelList groupIds(nGroups);
        forAll(groupIds, i)
        {
            const word sName = subsetPrefix + Foam::name(i);
            label setId = surf.facetSubsetIndex(sName);
            if (setId >= 0)
                surf.removeFacetSubset(setId);

            groupIds[i] = surf.addFacetSubset(sName);
        }

        forAll(triangleInRegion, tI)
            surf.addFacetToSubset(groupIds[triangleInRegion[tI]], tI);
    }

    return nGroups;
}


void Foam::Module::triSurfaceChecks::calculateBoundingBox
(
    const triSurf& surf,
    boundBox& bb
)
{
    bb.min() = Foam::min(surf.points());
    bb.max() = Foam::max(surf.points());
}


Foam::label Foam::Module::triSurfaceChecks::checkCollocatedPoints
(
    const triSurf& surf,
    DynamicList<label>& collocatedPoints,
    const scalar distTol
)
{
    collocatedPoints.clear();

    meshOctree octree(surf);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 30);

    const pointField& pts = surf.points();

    boolList collocated(pts.size(), false);

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(collocated, pI)
    {
        const point& p = pts[pI];

        boundBox bb(p, p);

        bb.min() -= point(distTol, distTol, distTol);
        bb.max() += point(distTol, distTol, distTol);

        DynamicList<label> leavesInBox;
        octree.findLeavesContainedInBox(bb, leavesInBox);

        forAll(leavesInBox, i)
        {
            const label leafI = leavesInBox[i];

            DynamicList<label> trianglesInBox;

            octree.containedTriangles(leafI, trianglesInBox);

            forAll(trianglesInBox, j)
            {
                const label triJ = trianglesInBox[j];
                const labelledTri& nt = surf[triJ];

                forAll(nt, tpI)
                {
                    if (nt[tpI] == pI)
                        continue;

                    if (magSqr(pts[nt[tpI]] - p) < sqr(distTol))
                    {
                        collocated[pI] = true;
                        collocated[nt[tpI]] = true;
                    }
                }
            }
        }
    }

    forAll(collocated, pI)
        if (collocated[pI])
            collocatedPoints.append(pI);

    return collocatedPoints.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkCollocatedPoints
(
    triSurf& surf,
    const word subsetName,
    const scalar distTol
)
{
    DynamicList<label> collocatedPoints;

    if (checkCollocatedPoints(surf, collocatedPoints, distTol))
    {
        label setId = surf.pointSubsetIndex(subsetName);
        if (setId >= 0)
            surf.removePointSubset(setId);
        setId = surf.addPointSubset(subsetName);

        forAll(collocatedPoints, i)
            surf.addPointToSubset(setId, collocatedPoints[i]);
    }

    return collocatedPoints.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkSelfIntersections
(
    const triSurf& surf,
    DynamicList<label>& badFaces,
    const scalar tol
)
{
    badFaces.clear();

    meshOctree octree(surf);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 30);

    const pointField& pts = surf.points();

    boolList intersected(surf.size(), false);

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(intersected, tI)
    {
        const labelledTri& tri = surf[tI];

        const triangle<point, point> currTri
        (
            pts[tri[0]],
            pts[tri[1]],
            pts[tri[2]]
        );

        boundBox bb(pts[tri[0]], pts[tri[0]]);
        for (label i = 1; i < 3; ++i)
        {
            bb.min() = Foam::min(bb.min(), pts[tri[i]]);
            bb.max() = Foam::max(bb.max(), pts[tri[i]]);
        }

        bb.min() -= point(tol, tol, tol);
        bb.max() += point(tol, tol, tol);

        DynamicList<label> leavesInBox;
        octree.findLeavesContainedInBox(bb, leavesInBox);

        forAll(leavesInBox, i)
        {
            const label leafI = leavesInBox[i];

            DynamicList<label> trianglesInBox;

            octree.containedTriangles(leafI, trianglesInBox);

            forAll(trianglesInBox, j)
            {
                const label triJ = trianglesInBox[j];
                const labelledTri& nt = surf[triJ];

                if (tI >= triJ)
                    continue;

                DynamicList<label> sharedPoints;
                sharedPoints.reserve(3);
                forAll(nt, pJ)
                {
                    forAll(tri, pI)
                    if (tri[pI] == nt[pJ])
                        sharedPoints.append(pJ);
                }

                if (sharedPoints.size() >= 2)
                {
                    // triangles share an edge and cannot self-intersect
                    continue;
                }
                else if (sharedPoints.size() == 1)
                {
                    // check if the opposite edge intersects the triangle
                    const point& s = pts[nt[(sharedPoints[0]+1)%3]];
                    const point& e = pts[nt[(sharedPoints[0]+2)%3]];

                    point intersection;
                    const bool lineIntersection =
                        help::triLineIntersection(currTri, s, e, intersection);

                    if (lineIntersection)
                    {
                        intersected[tI] = true;
                        intersected[triJ] = true;
                    }
                }
                else
                {
                    // triangles do not share any vertices
                    const triangle<point, point> neiTri
                    (
                        pts[nt[0]],
                        pts[nt[1]],
                        pts[nt[2]]
                    );

                    const bool intersect =
                        help::doTrianglesIntersect
                        (
                            currTri,
                            neiTri,
                            tol
                        );

                    if (intersect)
                    {
                        intersected[tI] = true;
                        intersected[triJ] = true;

                    }
                }
            }
        }
    }

    forAll(intersected, tI)
    {
        if (intersected[tI])
            badFaces.append(tI);
    }

    return badFaces.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkSelfIntersections
(
    triSurf& surf,
    const word subsetName,
    const scalar tol
)
{
    DynamicList<label> badFacets;

    if (checkSelfIntersections(surf, badFacets, tol))
    {
        label setId = surf.facetSubsetIndex(subsetName);
        if (setId >= 0)
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetName);

        forAll(badFacets, i)
            surf.addFacetToSubset(setId, badFacets[i]);
    }

    return badFacets.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkOverlaps
(
    const triSurf& surf,
    DynamicList<label>& badFaces,
    const scalar tol,
    const scalar angleTol
)
{
    badFaces.clear();

    meshOctree octree(surf);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 30);

    const scalar cosVal = Foam::cos(angleTol*M_PI/180.0);

    const pointField& pts = surf.points();

    boolList intersected(surf.size(), false);

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(intersected, tI)
    {
        const labelledTri& tri = surf[tI];

        const triangle<point, point> currTri
        (
            pts[tri[0]],
            pts[tri[1]],
            pts[tri[2]]
        );

        boundBox bb(pts[tri[0]], pts[tri[0]]);
        for (label i = 1; i < 3; ++i)
        {
            bb.min() = Foam::min(bb.min(), pts[tri[i]]);
            bb.max() = Foam::max(bb.max(), pts[tri[i]]);
        }

        bb.min() -= point(tol, tol, tol);
        bb.max() += point(tol, tol, tol);

        DynamicList<label> leavesInBox;
        octree.findLeavesContainedInBox(bb, leavesInBox);

        forAll(leavesInBox, i)
        {
            const label leafI = leavesInBox[i];

            DynamicList<label> trianglesInBox;

            octree.containedTriangles(leafI, trianglesInBox);

            forAll(trianglesInBox, j)
            {
                const label triJ = trianglesInBox[j];
                const labelledTri& nt = surf[triJ];

                if (tI >= triJ)
                    continue;

                const triangle<point, point> neiTri
                (
                    pts[nt[0]],
                    pts[nt[1]],
                    pts[nt[2]]
                );

                DynamicList<point> commonPolygon;

                const bool intersect =
                    help::doTrianglesOverlap
                    (
                        currTri,
                        neiTri,
                        commonPolygon,
                        tol,
                        cosVal
                    );

                if (intersect)
                {
                    intersected[tI] = true;
                    intersected[triJ] = true;

                }
            }
        }
    }

    forAll(intersected, tI)
    {
        if (intersected[tI])
            badFaces.append(tI);
    }

    return badFaces.size();
}


Foam::label Foam::Module::triSurfaceChecks::checkOverlaps
(
    triSurf& surf,
    const word subsetName,
    const scalar tol,
    const scalar angleTol
)
{
    DynamicList<label> badFaces;

    if (checkOverlaps(surf, badFaces, tol, angleTol))
    {
        label setId = surf.facetSubsetIndex(subsetName);
        if (setId >= 0)
            surf.removeFacetSubset(setId);
        setId = surf.addFacetSubset(subsetName);

        forAll(badFaces, i)
            surf.addFacetToSubset(setId, badFaces[i]);
    }

    return badFaces.size();
}


// ************************************************************************* //
