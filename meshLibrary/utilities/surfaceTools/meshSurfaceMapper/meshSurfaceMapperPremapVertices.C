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
#include "meshOctree.H"
#include "triSurf.H"
#include "labelledPointScalar.H"
#include "helperFunctions.H"
#include "meshSurfaceOptimizer.H"
#include "labelListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceMapper::preMapVertices(const label nIterations)
{
    Info<< "Smoothing mesh surface before mapping." << endl;

    const labelList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const labelListList2D& pointFaces = surfaceEngine_.pointFaces();
    const labelListList2D& pointInFace = surfaceEngine_.pointInFaces();
    const labelList& bp = surfaceEngine_.bp();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

    const triSurf& surf = meshOctree_.surface();

    List<labelledPointScalar> preMapPositions(boundaryPoints.size());
    List<DynamicList<scalar>> faceCentreDistances(bFaces.size());

    # pragma omp parallel for schedule(dynamic, 20)
    forAll(bFaces, bfI)
    {
        const point& c = faceCentres[bfI];
        const face& bf = bFaces[bfI];

        faceCentreDistances[bfI].setSize(bf.size());

        forAll(bf, pI)
        {
            faceCentreDistances[bfI][pI] = magSqr(points[bf[pI]] - c);
        }
    }

    for (label iterI = 0; iterI < nIterations; ++iterI)
    {
        // find patches in the vicinity of a boundary face
        labelListList2D boundaryPointPatches(boundaryPoints.size());
        # pragma omp parallel for schedule(dynamic, 50)
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];
            scalar boxSize(0.0);
            forAll(bf, pI)
            {
                boxSize =
                    Foam::max
                    (
                        boxSize,
                        mag(faceCentres[bfI] - points[bf[pI]])
                    );
            }

            const boundBox bb
            (
                faceCentres[bfI] - vector(boxSize, boxSize, boxSize),
                faceCentres[bfI] + vector(boxSize, boxSize, boxSize)
            );

            DynamicList<label> containedLeaves;
            meshOctree_.findLeavesContainedInBox(bb, containedLeaves);

            DynamicList<label> patches;
            forAll(containedLeaves, clI)
            {
                DynamicList<label> ct;
                meshOctree_.containedTriangles(containedLeaves[clI], ct);

                forAll(ct, i)
                    appendUniq(patches, surf[ct[i]].region());
            }

            scalar metric(VGREAT);
            label bestPatch(-1);
            forAll(patches, ptchI)
            {
                const scalar m = faceMetricInPatch(bfI, patches[ptchI]);

                if (m < metric)
                {
                    metric = m;
                    bestPatch = patches[ptchI];
                }
            }

            forAll(bf, pI)
            {
                # pragma omp critical
                appendUniq(boundaryPointPatches[bp[bf[pI]]], bestPatch);
            }
        }

        // use the shrinking laplace first
        # pragma omp parallel for schedule(dynamic, 40)
        forAll(pointFaces, bpI)
        {
            labelledPointScalar lp(bpI, vector::zero, 0.0);

            const point& p = points[boundaryPoints[bpI]];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces[bpI][pfI];
                const point& fc = faceCentres[pointFaces[bpI][pfI]];
                const label pos = pointInFace[bpI][pfI];
                const scalar w
                (
                    max(magSqr(p - fc) / faceCentreDistances[bfI][pos], SMALL)
                );
                lp.coordinates() += w*faceCentres[bfI];
                lp.scalarValue() += w;
            }

            preMapPositions[bpI] = lp;
        }

        // create the surface modifier and move the surface points
        meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

        # pragma omp parallel for schedule(dynamic, 50)
        forAll(boundaryPoints, bpI)
        {
            labelledPointScalar& lps = preMapPositions[bpI];

            lps.coordinates() /= lps.scalarValue();

            const point& p = points[boundaryPoints[bpI]];

            //label patch, nearestTri;
            point pMap = p;
            scalar dSq;

            if (boundaryPointPatches[bpI].size() == 1)
            {
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    pMap,
                    dSq,
                    nt,
                    boundaryPointPatches[bpI][0],
                    lps.coordinates()
                );
            }
            else
            {
                meshOctree_.findNearestPointToPatches
                (
                    pMap,
                    dSq,
                    lps.coordinates(),
                    boundaryPointPatches[bpI]
                );
            }

            const point newP = p + 0.5*(pMap - p);

            surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);
        }

        // update the surface geometry of the
        surfaceModifier.updateGeometry();

        meshSurfaceOptimizer(surfaceEngine_, meshOctree_).untangleSurface();

        surfaceModifier.updateGeometry();
    }

    Info<< "Finished smoothing mesh surface before mapping." << endl;
}


// ************************************************************************* //
