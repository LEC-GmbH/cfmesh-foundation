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
#include "meshSurfaceOptimizer.H"
#include "plane.H"
#include "Map.H"
#include "surfaceOptimizer.H"
#include "helperFunctions.H"
#include "partTriMeshSimplex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::Module::partTriMesh&
Foam::Module::meshSurfaceOptimizer::triMesh() const
{
    if (!triMeshPtr_)
    {
        calculateTrianglesAndAddressing();
    }

    return *triMeshPtr_;
}


inline void Foam::Module::meshSurfaceOptimizer::updateTriMesh
(
    const DynamicList<label>& selPoints
)
{
    if (!triMeshPtr_)
    {
        FatalErrorInFunction
            << "triMeshPtr_ is not allocated " << abort(FatalError);
    }

    triMeshPtr_->updateVertices(selPoints);
}


inline void Foam::Module::meshSurfaceOptimizer::updateTriMesh()
{
    if (!triMeshPtr_)
    {
        FatalErrorInFunction
            << "triMeshPtr_ is not allocated " << abort(FatalError);
    }

    triMeshPtr_->updateVertices();
}


inline bool Foam::Module::meshSurfaceOptimizer::transformIntoPlane
(
    const label bpI,
    const plane& pl,
    vector& vecX,
    vector& vecY,
    DynamicList<point>& pts,
    DynamicList<triFace>& trias
) const
{
    const partTriMesh& triMesh = this->triMesh();
    const label triPointI = triMesh.meshSurfacePointLabelInTriMesh()[bpI];
    if (triPointI < 0)
    {
        return false;
    }

    partTriMeshSimplex simplex(triMesh, triPointI);

    const DynamicList<point>& sPts = simplex.pts();

    //- create vecX and vecY
    const point& p = simplex.centrePoint();
    pts.setSize(0);
    bool found(false);
    forAll(sPts, pI)
    {
        const point sp = pl.nearestPoint(sPts[pI]);
        const scalar d = mag(sp - p);
        if (d > VSMALL)
        {
            vecX = sp - p;
            vecX /= d;
            vecY = pl.normal() ^ vecX;
            vecY /= mag(vecY);
            found = true;
            break;
        }
    }

    if (!found)
    {
        return false;
    }

    trias = simplex.triangles();

    //- transform the vertices
    pts.setSize(sPts.size());

    forAll(sPts, pI)
    {
        const point pt
        (
            ((sPts[pI] - p) & vecX),
            ((sPts[pI] - p) & vecY),
            0.0
        );

        pts[pI] = pt;
    }

    return found;
}


inline Foam::point Foam::Module::meshSurfaceOptimizer::newPositionLaplacian
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const labelListList2D& pPoints = surfaceEngine_.pointPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    if (vertexType_[bpI] & LOCKED)
    {
        return points[bPoints[bpI]];
    }

    vector newP(vector::zero);
    if (transformIntoPlane)
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if (magSqr(pNormal) < VSMALL)
        {
            return points[bPoints[bpI]];
        }

        plane pl(points[bPoints[bpI]], pNormal);

        DynamicList<point> projectedPoints;
        projectedPoints.setSize(pPoints[bpI].size());
        forAllRow(pPoints, bpI, pI)
        {
            projectedPoints[pI] =
                pl.nearestPoint(points[bPoints[pPoints[bpI][pI]]]);
        }

        forAll(projectedPoints, pI)
        {
            newP += projectedPoints[pI];
        }

        newP /= projectedPoints.size();
    }
    else
    {
        forAllRow(pPoints, bpI, pI)
        {
            newP += points[bPoints[pPoints[bpI][pI]]];
        }

        newP /= pPoints[bpI].size();
    }

    return newP;
}


inline Foam::point Foam::Module::meshSurfaceOptimizer::newPositionLaplacianFC
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const labelListList2D& pointFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    if (vertexType_[bpI] & LOCKED)
    {
        return points[bPoints[bpI]];
    }

    vector newP(vector::zero);

    if (transformIntoPlane)
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if (magSqr(pNormal) < VSMALL)
        {
            return points[bPoints[bpI]];
        }

        plane pl(points[bPoints[bpI]], pNormal);

        DynamicList<point> projectedPoints;
        projectedPoints.setSize(pointFaces[bpI].size());
        forAllRow(pointFaces, bpI, pfI)
        {
            projectedPoints[pfI] =
                pl.nearestPoint(faceCentres[pointFaces[bpI][pfI]]);
        }

        forAll(projectedPoints, pI)
        {
            newP += projectedPoints[pI];
        }

        newP /= projectedPoints.size();
    }
    else
    {
        forAllRow(pointFaces, bpI, pfI)
        {
            newP += faceCentres[pointFaces[bpI][pfI]];
        }

        newP /= pointFaces[bpI].size();
    }

    return newP;
}


inline Foam::point Foam::Module::meshSurfaceOptimizer::newPositionLaplacianWFC
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const labelListList2D& pointFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceAreas = surfaceEngine_.faceNormals();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    if (vertexType_[bpI] & LOCKED)
    {
        return points[bPoints[bpI]];
    }

    vector newP(vector::zero);

    if (transformIntoPlane)
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if (magSqr(pNormal) < VSMALL)
        {
            return points[bPoints[bpI]];
        }

        plane pl(points[bPoints[bpI]], pNormal);

        DynamicList<point> projectedPoints;
        projectedPoints.setSize(pointFaces[bpI].size());
        forAllRow(pointFaces, bpI, pfI)
        {
            projectedPoints[pfI] =
                pl.nearestPoint(faceCentres[pointFaces[bpI][pfI]]);
        }

        scalar sumW(0.0);
        forAll(projectedPoints, pI)
        {
            const label bfI = pointFaces[bpI][pI];
            const scalar w = (faceAreas[bfI] & faceAreas[bfI]);
            newP += w*projectedPoints[pI];
            sumW += w;
        }

        newP /= Foam::max(sumW, VSMALL);
    }
    else
    {
        scalar sumW(0.0);
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces[bpI][pfI];
            const scalar w = (faceAreas[bfI] & faceAreas[bfI]);
            newP += w*faceCentres[pointFaces[bpI][pfI]];
            sumW += w;
        }

        newP /= Foam::max(sumW, VSMALL);
    }

    return newP;
}


inline Foam::point
Foam::Module::meshSurfaceOptimizer::newPositionSurfaceOptimizer
(
    const label bpI,
    const scalar tol
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    if (vertexType_[bpI] & LOCKED)
    {
        return points[bPoints[bpI]];
    }

    //- project vertices onto the plane
    const vector& pNormal = surfaceEngine_.pointNormals()[bpI];
    if (magSqr(pNormal) < VSMALL)
    {
        return points[bPoints[bpI]];
    }

    const plane pl(points[bPoints[bpI]], pNormal);

    DynamicList<point> pts;
    DynamicList<triFace> trias;
    vector vecX, vecY;
    bool success = this->transformIntoPlane(bpI, pl, vecX, vecY, pts, trias);
    if (!success)
    {
        //Warning << "Cannot transform to plane" << endl;
        return points[bPoints[bpI]];
    }

    surfaceOptimizer so(pts, trias);
    const point newPoint = so.optimizePoint(tol);

    const point newP
    (
        points[bPoints[bpI]] +
        vecX*newPoint.x() +
        vecY*newPoint.y()
    );

    if (help::isnan(newP) || help::isinf(newP))
    {
        WarningInFunction
            << "Cannot move point " << bpI << endl;

        return points[bPoints[bpI]];
    }

    return newP;
}


inline Foam::point Foam::Module::meshSurfaceOptimizer::newEdgePositionLaplacian
(
    const label bpI
) const
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeList& edges = surfaceEngine_.edges();
    const labelListList2D& bpEdges = surfaceEngine_.boundaryPointEdges();
    const pointFieldPMG& points = surfaceEngine_.points();

    if (vertexType_[bpI] & LOCKED)
    {
        return points[bPoints[bpI]];
    }

    const labelHashSet& featureEdges = partitionerPtr_->featureEdges();

    DynamicList<label> edgePoints;

    forAllRow(bpEdges, bpI, i)
    {
        const label beI = bpEdges[bpI][i];

        if (featureEdges.found(beI))
        {
            edgePoints.append(edges[beI].otherVertex(bPoints[bpI]));
        }
    }

    if (edgePoints.size() != 2)
    {
        return points[bPoints[bpI]];
    }

    vector pos(vector::zero);
    forAll(edgePoints, epI)
    {
        pos += points[edgePoints[epI]];
    }

    pos /= edgePoints.size();

    return pos;
}


template<class labelListType>
inline void Foam::Module::meshSurfaceOptimizer::lockBoundaryFaces
(
    const labelListType& l
)
{
    lockedSurfaceFaces_ = l;

    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const labelList& bp = surfaceEngine_.bp();

    # pragma omp parallel for schedule(dynamic, 20)
    forAll(lockedSurfaceFaces_, lfI)
    {
        const face& bf = bFaces[lockedSurfaceFaces_[lfI]];

        forAll(bf, pI)
        {
            vertexType_[bp[bf[pI]]] |= LOCKED;
        }
    }
}


template<class labelListType>
inline void Foam::Module::meshSurfaceOptimizer::lockBoundaryPoints
(
    const labelListType& l
)
{
    # pragma omp parallel for schedule(dynamic, 50)
    forAll(l, i)
    {
        const label bpI = l[i];

        vertexType_[bpI] |= LOCKED;
    }
}


inline void Foam::Module::meshSurfaceOptimizer::lockFeatureEdges()
{
    # pragma omp parallel for schedule(dynamic, 50)
    forAll(vertexType_, bpI)
    {
        if (vertexType_[bpI] & (EDGE | CORNER))
        {
            vertexType_[bpI] |= LOCKED;
        }
    }
}


// ************************************************************************* //
