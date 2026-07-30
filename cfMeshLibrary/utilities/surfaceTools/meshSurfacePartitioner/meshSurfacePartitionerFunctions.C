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

#include "meshSurfacePartitioner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfacePartitioner::calculateCornersEdgesAndAddressing()
{
    const labelList& bPoints = meshSurface_.boundaryPoints();
    const labelList& bp = meshSurface_.bp();
    const edgeList& edges = meshSurface_.edges();
    const labelListList2D& edgeFaces = meshSurface_.edgeFaces();
    const labelListList2D& pointFaces = meshSurface_.pointFaces();

    corners_.clear();
    edgePoints_.clear();

    // count the number of patches
    label nPatches(0);
    # pragma omp parallel
    {
        label localMax(0);

        # pragma omp for
        forAll(facePatch_, bfI)
        {
            localMax = Foam::max(localMax, facePatch_[bfI]);
        }

        # pragma omp critical
        nPatches = Foam::max(localMax, nPatches);
    }
    ++nPatches;

    // set the size and starting creating addressing
    patchPatches_.setSize(nPatches);

    nEdgesAtPoint_.clear();
    nEdgesAtPoint_.setSize(bPoints.size());
    nEdgesAtPoint_ = 0;

    featureEdges_.clear();

    forAll(edgeFaces, edgeI)
    {
        if (edgeFaces[edgeI].size() != 2)
        {
            continue;
        }

        const label patch0 = facePatch_[edgeFaces[edgeI][0]];
        const label patch1 = facePatch_[edgeFaces[edgeI][1]];

        if (patch0 != patch1)
        {
            const edge& e = edges[edgeI];
            ++nEdgesAtPoint_[bp[e.start()]];
            ++nEdgesAtPoint_[bp[e.end()]];

            patchPatches_[patch0].insert(patch1);
            patchPatches_[patch1].insert(patch0);

            featureEdges_.insert(edgeI);
        }
    }

    // mark edges and corners
    forAll(nEdgesAtPoint_, bpI)
    {
        if (nEdgesAtPoint_[bpI] > 2)
        {
            corners_.insert(bpI);
        }
        else if (nEdgesAtPoint_[bpI] == 2)
        {
            edgePoints_.insert(bpI);
        }
    }

    // find patches at a surface points
    pointPatches_.setSize(pointFaces.size());
    forAll(pointFaces, bpI)
    {
        forAllRow(pointFaces, bpI, pfI)
        {
            appendUniq(pointPatches_[bpI], facePatch_[pointFaces[bpI][pfI]]);
        }
    }

    label counter = corners_.size();
    Info<< "Found " << counter
        << " corners at the surface of the volume mesh" << endl;
    counter = edgePoints_.size();
    Info<< "Found " << counter
        << " edge points at the surface of the volume mesh" << endl;
}


// ************************************************************************* //
