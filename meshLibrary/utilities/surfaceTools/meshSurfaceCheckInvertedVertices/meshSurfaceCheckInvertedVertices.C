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

#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfacePartitioner.H"
#include "boolList.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceCheckInvertedVertices::checkVertices()
{
    const labelList& facePatch = surfacePartitioner_.boundaryFacePatches();
    const meshSurfaceEngine& mse = surfacePartitioner_.surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelList& bp = mse.bp();
    const labelListList2D& pointFaces = mse.pointFaces();
    const labelListList2D& pointInFaces = mse.pointInFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorField& pNormals = mse.pointNormals();
    const vectorField& fCentres = mse.faceCentres();
    const vectorField& fNormals = mse.faceNormals();

    const labelHashSet& corners = surfacePartitioner_.corners();
    const labelHashSet& edgePoints = surfacePartitioner_.edgePoints();

    typedef std::map<label, vector> ltvMap;
    typedef std::map<label, ltvMap> lltvMap;
    lltvMap pointPatchNormal;

    for (const label bpI : corners.toc())
    {
        if (activePointsPtr_ && !activePointsPtr_->operator[](bpI))
            continue;

        ltvMap& patchNormal = pointPatchNormal[bpI];

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces[bpI][pfI];
            const label patchI = facePatch[bfI];

            if (patchNormal.find(patchI) == patchNormal.end())
            {
                patchNormal[patchI] = fNormals[bfI];
            }
            else
            {
                patchNormal[patchI] += fNormals[bfI];
            }
        }
    }

    for (const label bpI : edgePoints.toc())
    {
        if (activePointsPtr_ && !activePointsPtr_->operator[](bpI))
            continue;

        ltvMap& patchNormal = pointPatchNormal[bpI];

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces[bpI][pfI];
            const label patchI = facePatch[bfI];

            if (patchNormal.find(patchI) == patchNormal.end())
            {
                patchNormal[patchI] = fNormals[bfI];
            }
            else
            {
                patchNormal[patchI] += fNormals[bfI];
            }
        }
    }

    forAllIter(lltvMap, pointPatchNormal, it)
    {
        auto& patchNormal = pointPatchNormal[it->first];

        forAllIter(ltvMap, patchNormal, pIt)
        {
            const scalar magv = mag(pIt->second) + VSMALL;

            pIt->second /= magv;
        }
    }

    invertedVertices_.clear();

    # pragma omp parallel for if (pointFaces.size() > 100) \
    schedule(dynamic, 20)
    forAll(pointFaces, bpI)
    {
        if (activePointsPtr_ && !activePointsPtr_->operator[](bpI))
            continue;

        forAllRow(pointFaces, bpI, pfI)
        {
            const label pI = pointInFaces[bpI][pfI];
            const label bfI = pointFaces[bpI][pfI];

            vector pNormal = pNormals[bpI];

            if (pointPatchNormal.find(bpI) != pointPatchNormal.end())
                pNormal = pointPatchNormal[bpI][facePatch[bfI]];

            const face& bf = bFaces[bfI];

            // chech the first triangle (with the next node)
            triangle<point, point> triNext
            (
                points[bf[pI]],
                points[bf.nextLabel(pI)],
                fCentres[bfI]
            );

            vector nNext = triNext.area();
            scalar mNext = mag(nNext);

            // face has zero area
            if (mNext < VSMALL)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }
            else
            {
                nNext /= mNext;
            }

            // collocated points
            if (magSqr(triNext.a() - triNext.b()) < VSMALL)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }
            if (magSqr(triNext.c() - triNext.a()) < VSMALL)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }

            // normal vector is not visible
            if ((nNext & pNormal) < 0.0)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }

            // check the second triangle (with previous node)
            triangle<point, point> triPrev
            (
                points[bf[pI]],
                fCentres[bfI],
                points[bf.prevLabel(pI)]
            );

            vector nPrev = triPrev.area();
            scalar mPrev = mag(nPrev);

            // face has zero area
            if (mPrev < VSMALL)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }
            else
            {
                nPrev /= mPrev;
            }

            // collocated points
            if (magSqr(triPrev.a() - triPrev.b()) < VSMALL)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }
            if (magSqr(triPrev.c() - triPrev.a()) < VSMALL)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }

            // normal vector is not visible
            if ((nPrev & pNormal) < 0.0)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }

            // check whether the normals of both triangles
            // point in the same direction
            if ((nNext & nPrev) < 0.0)
            {
                # pragma omp critical
                invertedVertices_.insert(bf[pI]);

                break;
            }
        }
    }

    // check if there exist concave faces
    # pragma omp parallel for schedule(dynamic, 50)
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        DynamicList<bool> OkPoints;
        if (!help::isFaceConvexAndOk(bf, points, OkPoints))
        {
            forAll(OkPoints, pI)
            {
                if (activePointsPtr_ && !(*activePointsPtr_)[bp[bf[pI]]])
                    continue;

                if (!OkPoints[pI])
                {
                    # pragma omp critical
                    {
                        invertedVertices_.insert(bf[pI]);
                    }
                }
            }
        }
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshSurfaceCheckInvertedVertices::meshSurfaceCheckInvertedVertices
(
    const meshSurfacePartitioner& mpart
)
:
    surfacePartitioner_(mpart),
    activePointsPtr_(nullptr),
    invertedVertices_()
{
    checkVertices();
}


Foam::Module::meshSurfaceCheckInvertedVertices::meshSurfaceCheckInvertedVertices
(
    const meshSurfacePartitioner& mpart,
    const boolList& activePoints
)
:
    surfacePartitioner_(mpart),
    activePointsPtr_(&activePoints),
    invertedVertices_()
{
    checkVertices();
}


// ************************************************************************* //
