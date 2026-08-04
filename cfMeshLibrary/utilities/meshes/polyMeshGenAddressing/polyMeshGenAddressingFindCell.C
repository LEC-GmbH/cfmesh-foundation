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

Description
    Point-in-cell test via tet decomposition (same principle as
    OpenFOAM's own polyMeshTetDecomposition/meshSearch::findCell): each
    cell is decomposed into tets by fanning every face into triangles
    around that face's own centre, then forming a tet from each such
    triangle and the cell centre. A point lies in the cell iff it lies
    in at least one of these tets. This decomposition is exact (the tets
    exactly tile the cell) for any cell that is star-shaped from its own
    centre, which covers the castellated cells this is used on.

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "meshOctree.H"
#include "Map.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

//- six times the signed volume of tet (a, b, c, d)
inline Foam::scalar sixTetVol
(
    const Foam::point& a,
    const Foam::point& b,
    const Foam::point& c,
    const Foam::point& d
)
{
    return ((b - a) ^ (c - a)) & (d - a);
}


//- is p inside (or on the boundary of) tet (a, b, c, d)?
bool pointInTet
(
    const Foam::point& p,
    const Foam::point& a,
    const Foam::point& b,
    const Foam::point& c,
    const Foam::point& d
)
{
    const Foam::scalar v0 = sixTetVol(a, b, c, d);

    if (Foam::mag(v0) < Foam::SMALL)
    {
        // degenerate tet (zero volume) - cannot contain any point
        return false;
    }

    // v1..v4 always sum exactly to v0, for any p - a point is inside the
    // tet iff all four sub-tet volumes have the same sign as v0 (a small
    // tolerance, relative to v0, admits points on a face/edge/vertex)
    const Foam::scalar v1 = sixTetVol(p, b, c, d);
    const Foam::scalar v2 = sixTetVol(a, p, c, d);
    const Foam::scalar v3 = sixTetVol(a, b, p, d);
    const Foam::scalar v4 = sixTetVol(a, b, c, p);

    const Foam::scalar tol = Foam::SMALL*Foam::mag(v0);

    if (v0 > 0)
    {
        return (v1 > -tol) && (v2 > -tol) && (v3 > -tol) && (v4 > -tol);
    }

    return (v1 < tol) && (v2 < tol) && (v3 < tol) && (v4 < tol);
}

} // End anonymous namespace


bool Foam::Module::polyMeshGenAddressing::pointInCell
(
    const label cellI,
    const point& p
) const
{
    const cell& c = mesh_.cells()[cellI];
    const faceListPMG& faces = mesh_.faces();
    const pointFieldPMG& points = mesh_.points();
    const vectorField& fCtrs = faceCentres();
    const point& cc = cellCentres()[cellI];

    forAll(c, fI)
    {
        const face& f = faces[c[fI]];
        const point& fc = fCtrs[c[fI]];

        forAll(f, pI)
        {
            const point& p0 = points[f[pI]];
            const point& p1 = points[f[f.fcIndex(pI)]];

            if (pointInTet(p, cc, fc, p0, p1))
            {
                return true;
            }
        }
    }

    return false;
}


Foam::label Foam::Module::polyMeshGenAddressing::findCell
(
    const point& p
) const
{
    const cellListPMG& cells = mesh_.cells();

    forAll(cells, cellI)
    {
        if (pointInCell(cellI, p))
        {
            return cellI;
        }
    }

    return -1;
}


Foam::labelList Foam::Module::polyMeshGenAddressing::findCells
(
    const pointField& points,
    const meshOctree& octree
) const
{
    // bucket every cell by the octree leaf containing its own centroid -
    // cheap (one octree descent per cell, no point-in-cell testing at
    // all), and narrows each query below from "every cell in the mesh"
    // down to a handful of candidates that actually share a leaf with
    // the query point
    const vectorField& cCtrs = cellCentres();
    const boundBox& rootBox = octree.rootBox();

    Map<DynamicList<label>> leafToCells;
    forAll(cCtrs, cellI)
    {
        if (!rootBox.contains(cCtrs[cellI]))
        {
            // a cell centroid can only fall outside the octree's own
            // root box due to some other, unrelated meshing defect -
            // not this function's concern; skip it, the fallback scan
            // below still covers it if a query point ever needs it
            continue;
        }

        const label leafI = octree.findLeafContainingVertex(cCtrs[cellI]);

        leafToCells(leafI).append(cellI);
    }

    labelList result(points.size(), -1);

    forAll(points, pI)
    {
        const point& p = points[pI];

        if (rootBox.contains(p))
        {
            const label leafI = octree.findLeafContainingVertex(p);

            const Map<DynamicList<label>>::const_iterator iter =
                leafToCells.find(leafI);

            if (iter != leafToCells.end())
            {
                forAll(iter(), i)
                {
                    if (pointInCell(iter()[i], p))
                    {
                        result[pI] = iter()[i];
                        break;
                    }
                }
            }
        }

        if (result[pI] == -1)
        {
            // fast path missed (p outside the octree's root box, its
            // leaf has no cells recorded, or none of that leaf's
            // candidates actually contain p, e.g. a boundary-cut cell's
            // centroid landed in a neighbouring leaf) - fall back to an
            // exhaustive scan so correctness never depends on the
            // narrowing actually succeeding
            result[pI] = findCell(p);
        }
    }

    return result;
}


// ************************************************************************* //
