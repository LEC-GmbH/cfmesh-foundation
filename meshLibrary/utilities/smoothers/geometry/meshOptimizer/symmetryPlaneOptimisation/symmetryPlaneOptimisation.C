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
#include "symmetryPlaneOptimisation.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"
#include "polyMeshGenChecks.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::symmetryPlaneOptimisation::detectSymmetryPlanes()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    symmetryPlanes_.clear();

    typedef Map<Tuple2<vector, label>> mapType;
    mapType centreSum, normalSum;

    forAll(boundaries, patchI)
    {
        if (boundaries[patchI].patchType() == "symmetryPlane")
        {
            centreSum.insert(patchI, Tuple2<vector, label>(vector::zero, 0));
            Tuple2<vector, label>& cs = centreSum[patchI];

            normalSum.insert(patchI, Tuple2<vector, label>(vector::zero, 0));
            Tuple2<vector, label>& ns = normalSum[patchI];

            const label start = boundaries[patchI].patchStart();
            const label end = start + boundaries[patchI].patchSize();
            for (label faceI = start; faceI < end; ++faceI)
            {
                cs.first() += help::stableFaceCentre(faces[faceI], points);
                ns.first() += help::stableFaceArea(faces[faceI], points);
            }

            cs.second() = ns.second() = boundaries[patchI].patchSize();
        }
    }

    // create planes corresponding to each symmetry plane
    forAllConstIter(mapType, centreSum, it)
    {
        const point c = it().first()/it().second();

        const Tuple2<vector, label>& ns = normalSum[it.key()];
        const vector n = ns.first()/ns.second();

        symmetryPlanes_.insert(it.key(), plane(c, n));
    }
}


bool Foam::Module::symmetryPlaneOptimisation::pointInPlanes
(
    labelListList2D& pointInPlanes
) const
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    pointInPlanes.clear();
    pointInPlanes.setSize(points.size());

    bool foundProblematic(false);

    forAll(boundaries, patchI)
    {
        if (boundaries[patchI].patchType() == "symmetryPlane")
        {
            const label start = boundaries[patchI].patchStart();
            const label end = start + boundaries[patchI].patchSize();
            for (label faceI = start; faceI < end; ++faceI)
            {
                const face& f = faces[faceI];

                const point c = help::stableFaceCentre(f, points);
                scalar maxDist(0.0);
                forAll(f, pI)
                    maxDist = max(maxDist, mag(points[f[pI]] - c));

                forAll(f, pI)
                {
                    Map<plane>::const_iterator it =
                        symmetryPlanes_.find(patchI);
                    if (it != symmetryPlanes_.end())
                    {
                        const scalar d = it().distance(points[f[pI]]);

                        if (d > 0.5*maxDist)
                            foundProblematic = true;
                    }

                    appendUniq(pointInPlanes[f[pI]], patchI);
                }
            }
        }
    }

    return foundProblematic;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::symmetryPlaneOptimisation::symmetryPlaneOptimisation
(
    polyMeshGen& mesh
)
:
    mesh_(mesh)
{
    detectSymmetryPlanes();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::symmetryPlaneOptimisation::optimizeSymmetryPlanes()
{
    pointFieldPMG& points = mesh_.points();

    labelListList2D pointInPlane;
    if (pointInPlanes(pointInPlane))
    {
        Warning
            << "Detected large deviation from some symmetry planes."
            << "Please check your settings and ensure that all your patches"
            << " with type symmetryPlane are planar. If there exists a patch"
            << " that is not planar please use the symmetry type, instead."
            << " Skipping symmetryPlane optimisation..." << endl;

        return;
    }

    forAll(pointInPlane, pointI)
    {
        const label nPlanes = pointInPlane[pointI].size();

        if (nPlanes > 3)
        {
            WarningInFunction
                << "Point " << pointI << " is in more than three symmetry"
                << " planes. Cannot move it" << endl;

            continue;
        }

        point& p = points[pointI];

        if (nPlanes == 1)
        {
            // point is in a plane
            Map<plane>::const_iterator it =
                symmetryPlanes_.find(pointInPlane[pointI][0]);

            p = it().nearestPoint(p);
        }
        else if (nPlanes == 2)
        {
            // point is at the edge between two planes
            const Map<plane>::const_iterator it0 =
                symmetryPlanes_.find(pointInPlane[pointI][0]);
            const Map<plane>::const_iterator it1 =
                symmetryPlanes_.find(pointInPlane[pointI][1]);
            const plane& pl0 = it0();
            const plane& pl1 = it1();

            const plane::ray il = pl0.planeIntersect(pl1);

            vector n = il.dir();
            n /= Foam::max(mag(n), VSMALL);

            const scalar lProj = (p - il.refPoint()) & n;

            p = il.refPoint() + lProj*n;
        }
        else if (nPlanes == 3)
        {
            // points is a corner between three planes
            const Map<plane>::const_iterator it0 =
                symmetryPlanes_.find(pointInPlane[pointI][0]);
            const Map<plane>::const_iterator it1 =
                symmetryPlanes_.find(pointInPlane[pointI][1]);
            const Map<plane>::const_iterator it2 =
                symmetryPlanes_.find(pointInPlane[pointI][2]);
            const plane& pl0 = it0();
            const plane& pl1 = it1();
            const plane& pl2 = it2();

            p = pl0.planePlaneIntersect(pl1, pl2);
        }
    }

    labelHashSet badFaces;
    polyMeshGenChecks::checkFacePyramids(mesh_, false, VSMALL, &badFaces);

    if (badFaces.size())
    {
        WarningInFunction
            << "Bad quality or inverted faces found in the mesh" << endl;

        const label badFacesId = mesh_.addFaceSubset("invalidFaces");
        for (const label facei : badFaces.toc())
        {
            mesh_.addFaceToSubset(badFacesId, facei);
        }
    }
}


// ************************************************************************* //
