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

#include "triSurfaceCurvatureEstimator.H"
#include "matrix3D.H"
#include "quadricFitting.H"
#include "helperFunctions.H"
#include "HashSet.H"
#include "boolList.H"
#include "Map.H"
#include "labelListListOps.H"

#include "OFstream.H"

#include <map>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Module::triSurfaceCurvatureEstimator::calculateEdgeCurvature()
{
    const pointField& points = surface_.points();
    const DynamicList<edge>& edges = surface_.edges();
    const labelListList2D& pointEdges = surface_.pointEdges();
    const labelListList2D& edgeFaces = surface_.edgeFacets();

    edgePointCurvature_.setSize(points.size());
    boolList featureEdge(edges.size());

    # pragma omp parallel
    {
        # pragma omp for schedule(static, 1)
        forAll(edgePointCurvature_, i)
            edgePointCurvature_[i] = 0.0;

        // mark feature edges
        # pragma omp for schedule(static, 1)
        forAll(edgeFaces, eI)
        {
            if (edgeFaces[eI].size() != 2)
            {
                featureEdge[eI] = false;
                continue;
            }

            if (
                 surface_[edgeFaces[eI][0]].region() ==
                 surface_[edgeFaces[eI][1]].region()
            )
            {
                featureEdge[eI] = false;
            }
            else
            {
                featureEdge[eI] = true;
            }
        }

        # pragma omp barrier

        // loop through the points and calculate the curvature for points
        // attached to two feature edges
        # pragma omp for schedule(dynamic, 20)
        forAll(pointEdges, pI)
        {
            DynamicList<label> features;
            forAllRow(pointEdges, pI, peI)
            {
                const label edgeI = pointEdges[pI][peI];
                if (featureEdge[edgeI])
                    features.append(edgeI);
            }

            if (features.size() == 2)
            {
                // calculate the curvature and store it in the map
                vector e1 = edges[features[0]].vec(points);
                const scalar d1 = Foam::max(mag(e1), VSMALL);
                e1 /= d1;
                vector e2 = edges[features[1]].vec(points);
                const scalar d2 = Foam::max(mag(e2), VSMALL);
                e2 /= d2;

                scalar cs = e1 & e2;
                cs = Foam::min(1.0, cs);
                cs = Foam::max(-1.0, cs);

                const scalar curv =
                    Foam::acos(cs)/Foam::max(0.5*(d1 + d2), VSMALL);
                edgePointCurvature_[pI] = Foam::mag(curv);
            }
        }
    }
}


void Foam::Module::triSurfaceCurvatureEstimator::calculateSurfaceCurvatures()
{
    const labelListList2D& pointTriangles = surface_.pointFacets();

    const pointField& points = surface_.points();

    patchPositions_.setSize(surface_.size());
    gaussianCurvature_.setSize(points.size());
    meanCurvature_.setSize(points.size());
    maxCurvature_.setSize(points.size());
    minCurvature_.setSize(points.size());
    maxCurvatureVector_.setSize(points.size());
    minCurvatureVector_.setSize(points.size());

    List<DynamicList<label>> pointPatches(points.size());

    # pragma omp parallel for schedule(dynamic, 40)
    forAll(pointTriangles, pointI)
    {
        Map<DynamicList<label>> regionTriangles;
        Map<labelHashSet> otherLabels;
        Map<vector> normals;

        forAllRow(pointTriangles, pointI, ptI)
        {
            const label triI = pointTriangles[pointI][ptI];
            const label regionI = surface_[triI].region();

            if (!normals.found(regionI))
                normals.insert(regionI, vector::zero);
            if (regionTriangles.found(regionI))
                regionTriangles.insert
                (
                    regionI,
                    DynamicList<label>()
                );
            if (!otherLabels.found(regionI))
                otherLabels.insert(regionI, labelHashSet());

            regionTriangles(regionI).append(triI);

            forAll(surface_[triI], tpI)
            {
                const label pI = surface_[triI][tpI];

                if (pointI == pI)
                    continue;

                otherLabels[regionI].insert(pI);
                normals[regionI] += surface_[triI].area(points);
            }
        }

        forAllConstIter(Map<labelHashSet>, otherLabels, it)
        {
            const labelHashSet& currLabels = it();

            if (otherLabels[it.key()].size() > 5)
                continue;

            labelHashSet additionalPoints;
            for (const label neiPointI : currLabels.toc())
            {
                const DynamicList<label>& pTriangles = pointTriangles[neiPointI];

                forAll(pTriangles, ptI)
                {
                    const labelledTri& nTri = surface_[pTriangles[ptI]];

                    if (nTri.region() != it.key())
                        continue;

                    forAll(nTri, pI)
                        additionalPoints.insert(nTri[pI]);
                }
            }

            for (const label aIter : additionalPoints.toc())
            {
                otherLabels[it.key()].insert(aIter);
            }
        }

        forAllIter(Map<vector>, normals, nit)
        {
            nit() /= mag(nit());
        }

        forAllConstIter(Map<labelHashSet>, otherLabels, it)
        {
            const labelHashSet& labels = it();

            const label regionI = it.key();

            // store the patches
            pointPatches[pointI].append(regionI);

            // find the position of point in the patch triangles
            const DynamicList<label>& rTriangles = regionTriangles[regionI];
            forAll(rTriangles, i)
            {
                const label tI = rTriangles[i];

                label pos(-1);
                forAll(surface_[tI], j)
                    if (surface_[tI][j] == pointI)
                    {
                        pos = j;
                        break;
                    }

                patchPositions_[tI][pos] = gaussianCurvature_[pointI].size();
            }

            // store point coordinates
            DynamicList<point> op;

            for (const label lit : labels.toc())
            {
                op.append(points[lit]);
            }

            // fit the quadric patch to the surface
            quadricFitting qfit(points[pointI], normals[it.key()], op);

            // store curvatures
            gaussianCurvature_[pointI].append(qfit.gaussianCurvature());
            meanCurvature_[pointI].append(qfit.meanCurvature());
            maxCurvature_[pointI].append(qfit.maxCurvature());
            minCurvature_[pointI].append(qfit.minCurvature());
            maxCurvatureVector_[pointI].append(qfit.maxCurvatureVector());
            minCurvatureVector_[pointI].append(qfit.minCurvatureVector());
        }
    }

    // smooth curvatures using weighted Laplace
    List<DynamicList<scalar>> smoothMinCurv(points.size());
    List<DynamicList<scalar>> smoothMaxCurv(points.size());

    # pragma omp parallel
    {
        for (label iteration = 0; iteration < 2; ++iteration)
        {
            # pragma omp for schedule(static, 1)
            forAll(pointTriangles, pointI)
            {
                const DynamicList<label>& pTriangles = pointTriangles[pointI];

                // find neighbouring points for each patch
                Map<DynamicList<label>> patchNeiPoints;
                forAll(pointPatches[pointI], ppI)
                    patchNeiPoints.insert
                    (
                        pointPatches[pointI][ppI],
                        DynamicList<label>()
                    );

                forAll(pTriangles, ptI)
                {
                    const labelledTri& tri = surface_[pTriangles[ptI]];

                    if (!patchNeiPoints.found(tri.region()))
                        patchNeiPoints.insert(tri.region(), DynamicList<label>());

                    forAll(tri, pI)
                    {
                        const label neiPointI = tri[pI];

                        if (neiPointI == pointI)
                            continue;

                        appendUniq(patchNeiPoints[tri.region()], neiPointI);
                    }
                }

                smoothMinCurv[pointI].setSize(minCurvature_[pointI].size());
                smoothMaxCurv[pointI].setSize(maxCurvature_[pointI].size());

                // update min curvature for all point patches
                forAll(minCurvature_[pointI], patchI)
                {
                    const label cPatch = pointPatches[pointI][patchI];

                    scalar minCurv(0.0);
                    scalar maxCurv(0.0);

                    const DynamicList<label>& neiPoints = patchNeiPoints[cPatch];

                    if (neiPoints.size() == 0)
                    {
                        smoothMinCurv[pointI][patchI] =
                            minCurvature_[pointI][patchI];

                        smoothMaxCurv[pointI][patchI] =
                            maxCurvature_[pointI][patchI];
                    }

                    forAll(neiPoints, i)
                    {
                        const label neiPointI = neiPoints[i];

                        const label pos =
                            findIndex(pointPatches[neiPointI], cPatch);

                        minCurv += minCurvature_[neiPointI][pos];
                        maxCurv += maxCurvature_[neiPointI][pos];
                    }

                    minCurv /= neiPoints.size();
                    maxCurv /= neiPoints.size();

                    // store the value
                    smoothMinCurv[pointI][patchI] =
                        0.5*(minCurv + minCurvature_[pointI][patchI]);

                    smoothMaxCurv[pointI][patchI] =
                        0.5*(maxCurv + maxCurvature_[pointI][patchI]);
                }
            }

            # pragma omp barrier

            # pragma omp for schedule(static, 1)
            forAll(minCurvature_, pointI)
            {
                forAll(minCurvature_[pointI], i)
                {
                    minCurvature_[pointI][i] = smoothMinCurv[pointI][i];
                    maxCurvature_[pointI][i] = smoothMaxCurv[pointI][i];
                }
            }
        }

        # pragma omp barrier

        // update Gaussian and mean curvatures
        # pragma omp for schedule(static, 1)
        forAll(minCurvature_, pointI)
        {
            const DynamicList<scalar>& minCurv = minCurvature_[pointI];
            const DynamicList<scalar>& maxCurv = maxCurvature_[pointI];

            forAll(minCurv, i)
            {
                gaussianCurvature_[pointI][i] = minCurv[i]*maxCurv[i];
                meanCurvature_[pointI][i] = 0.5*(minCurv[i] + maxCurv[i]);
            }
        }
    }
}


// ************************************************************************* //
