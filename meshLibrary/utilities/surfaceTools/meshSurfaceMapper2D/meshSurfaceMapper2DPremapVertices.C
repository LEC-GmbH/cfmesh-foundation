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
#include "meshSurfaceMapper2D.H"
#include "polyMeshGen2DEngine.H"
#include "meshOctree.H"
#include "labelledPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceMapper2D::preMapVertices(const label nIterations)
{
    Info<< "Smoothing mesh surface before mapping. Iteration:" << flush;

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bp = surfaceEngine_.bp();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& edgeFaces = surfaceEngine_.edgeFaces();
    const edgeList& edges = surfaceEngine_.edges();

    List<labelledPoint> preMapPositions(activeBoundaryEdges_.size());

    for (label iterI = 0; iterI < nIterations; ++iterI)
    {
        labelLongList parBndEdges;

        // use the shrinking laplace first
        # pragma omp parallel for schedule(dynamic, 40)
        forAll(activeBoundaryEdges_, eI)
        {
            const label beI = activeBoundaryEdges_[eI];

            labelledPoint lp(0, vector::zero);

            if (edgeFaces.sizeOfRow(beI) == 2)
            {
                forAllRow(edgeFaces, beI, efI)
                {
                    ++lp.pointLabel();
                    lp.coordinates() += faceCentres[edgeFaces(beI, efI)];
                }

            }
            else if (edgeFaces.sizeOfRow(beI) == 1)
            {
                ++lp.pointLabel();
                lp.coordinates() += faceCentres[edgeFaces(beI, 0)];

                # pragma omp critical
                parBndEdges.append(eI);
            }

            // store the information
            preMapPositions[eI] = lp;
        }

        // calculate coordinates of points for searching
        # pragma omp parallel for schedule(dynamic, 50)
        forAll(activeBoundaryEdges_, eI)
        {
            labelledPoint& lp = preMapPositions[eI];

            if (lp.pointLabel() == 0)
            {
                Warning << "Surface edge " << activeBoundaryEdges_[eI]
                    << " has no active faces" << endl;
                continue;
            }

            lp.coordinates() /= lp.pointLabel();
        }

        // create the surface modifier and move the surface points
        meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

        // accumulate target positions per boundary point rather than
        // moving directly: adjacent edges share an endpoint, so two
        // threads computing different target positions for the same
        // shared vertex would otherwise race (last-write-wins - unlike a
        // plain index, a position has no natural "smaller wins" tie-break)
        // - average all edges' proposals for a shared vertex instead
        List<labelledPoint> newPositions
        (
            surfaceEngine_.boundaryPoints().size(),
            labelledPoint(0, vector::zero)
        );

        # pragma omp parallel for schedule(dynamic, 50)
        forAll(activeBoundaryEdges_, eI)
        {
            const label beI = activeBoundaryEdges_[eI];
            const edge& e = edges[beI];

            const point& p = points[e.start()];

            label patch, nt;
            point pMap = p;
            scalar dSq;

            meshOctree_.findNearestSurfacePoint
            (
                pMap,
                dSq,
                nt,
                patch,
                preMapPositions[eI].coordinates()
            );

            pMap.z() = p.z();
            point newP = 0.5*(pMap + p);

            const label bpStart = bp[e.start()];
            # pragma omp critical
            {
                ++newPositions[bpStart].pointLabel();
                newPositions[bpStart].coordinates() += newP;
            }

            newP.z() = points[e.end()].z();

            const label bpEnd = bp[e.end()];
            # pragma omp critical
            {
                ++newPositions[bpEnd].pointLabel();
                newPositions[bpEnd].coordinates() += newP;
            }
        }

        # pragma omp parallel for schedule(dynamic, 50)
        forAll(newPositions, bpI)
        {
            labelledPoint& lp = newPositions[bpI];

            if (lp.pointLabel() == 0)
                continue;

            surfaceModifier.moveBoundaryVertexNoUpdate
            (
                bpI,
                lp.coordinates()/lp.pointLabel()
            );
        }

        surfaceModifier.updateGeometry();

        Info<< "." << flush;
    }

    Info<< endl;
}


// ************************************************************************* //
