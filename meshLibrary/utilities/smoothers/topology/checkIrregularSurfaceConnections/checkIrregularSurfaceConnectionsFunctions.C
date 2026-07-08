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

#include "checkIrregularSurfaceConnections.H"
#include "polyMeshGenAddressing.H"
#include "sortEdgesIntoChains.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Module::checkIrregularSurfaceConnections::
checkAndFixCellGroupsAtBndVertices
(
    labelHashSet& badVertices,
    const bool removeConnections
)
{
    Info<< "Checking cells connected to surface vertices" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();

    polyMeshGenModifier meshModifier(mesh_);
    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    const VRWGraph& cellCells = mesh_.addressingData().cellCells();
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();

    // check cells connected to surface nodes
    // cells connected to a vertex must create a single loop when cells
    // are visited over faces from to the other
    label nBadVertices(0);
    DynList<label> frontCells;

    const label size = bPoints.size();
    # pragma omp parallel for private(frontCells) schedule(dynamic, 1000)
    for (label bpI = 0; bpI < size; ++bpI)
    {
        const label pointI = bPoints[bpI];

        Map<label> cellGroup(pointCells.sizeOfRow(pointI));

        label nGroup(0);

        forAllRow(pointCells, pointI, cI)
        {
            const label cellI = pointCells(pointI, cI);

            if (cellGroup.found(cellI))
                continue;

            cellGroup.insert(cellI, nGroup);
            frontCells.clear();
            frontCells.append(cellI);

            while (frontCells.size() != 0)
            {
                const label cLabel = frontCells.remove();

                forAllRow(cellCells, cLabel, nI)
                {
                    const label neiCell = cellCells(cLabel, nI);

                    if (cellGroup.found(neiCell))
                        continue;

                    if (pointCells.found(pointI, neiCell))
                    {
                        cellGroup.insert(neiCell, nGroup);
                        frontCells.append(neiCell);
                    }
                }

            }

            ++nGroup;
        }

        if (nGroup != 1)
        {
            # pragma omp critical
            {
                ++nBadVertices;
                badVertices.insert(pointI);

                if (removeConnections)
                {
                    const label nPoints = points.size();
                    forAllRow(pointCells, pointI, pcI)
                    {
                        const label cellI = pointCells(pointI, pcI);

                        if (cellGroup[cellI] != 0)
                        {
                            const cell& c = cells[cellI];

                            forAll(c, fI)
                            {
                                face& f = faces[c[fI]];

                                const label pos = f.which(pointI);

                                if (pos > -1)
                                    f[pos] = nPoints + cellGroup[cellI] - 1;
                            }
                        }
                    }

                    for (label i = 1; i < nGroup; ++i)
                    {
                        const point p = points[pointI];
                        points.append(p);
                    }
                }
            }
        }
    }

    Info<< "Found " << nBadVertices << " problematic vertices" << endl;
    Info<< "Finished checking cells connected to surface vertices" << endl;

    if (nBadVertices != 0)
    {
        clearMeshEngine();
        mesh_.clearAddressingData();

        return true;
    }

    return false;
}


bool Foam::Module::checkIrregularSurfaceConnections::checkEdgeFaceConnections
(
    labelHashSet& badVertices,
    const bool removeCells
)
{
    Info<< "Checking for non - manifold surface edges" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const VRWGraph& pointEdges = mse.boundaryPointEdges();

    labelHashSet badEdges;

    forAll(edgeFaces, edgeI)
    {
        if (edgeFaces.sizeOfRow(edgeI) > 2)
        {
            badVertices.insert(edges[edgeI].start());
            badVertices.insert(edges[edgeI].end());

            badEdges.insert(edgeI);
        }
    }

    const label nBadEdges = badEdges.size();
    Info<< "Found " << nBadEdges << " non - manifold edges" << endl;
    Info<< "Finished checking for non - manifold surface edges" << endl;

    if (nBadEdges != 0 && removeCells)
    {
        // remove all cells connected to the selected edge
        boolList removeCell(mesh_.cells().size(), false);

        const faceListPMG& faces = mesh_.faces();
        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        forAll(faces, faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                const edge e = f.faceEdge(pI);

                if ((bp[e[0]] != -1) && (bp[e[1]] != -1))
                {
                    const label bpI = bp[e[0]];

                    forAllRow(pointEdges, bpI, peI)
                    {
                        const label edgeI = pointEdges(bpI, peI);
                        if ((edges[edgeI] == e) && badEdges.found(edgeI))
                        {
                            removeCell[owner[faceI]] = true;
                            if (neighbour[faceI] < 0)
                                continue;
                            removeCell[neighbour[faceI]] = true;
                        }
                    }
                }
            }
        }

        polyMeshGenModifier(mesh_).removeCells(removeCell);
        clearMeshEngine();

        return true;
    }

    return false;
}


bool Foam::Module::checkIrregularSurfaceConnections::
checkFaceGroupsAtBndVertices
(
    labelHashSet& badVertices,
    const bool removeCells
)
{
    Info<< "Checking faces connections to surface vertices" << endl;

    labelHashSet invalidVertices(100);

    const meshSurfaceEngine& mse = surfaceEngine();

    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& faceFaces = mse.faceFaces();

    // check number of face groups
    DynList<label> front;
    const label size = pointFaces.size();
    # pragma omp parallel for private(front) schedule(dynamic)
    for (label bpI = 0; bpI < size; ++bpI)
    {
        Map<label> faceGroup(pointFaces.sizeOfRow(bpI));
        label nGroup(0);

        forAllRow(pointFaces, bpI, pfI)
        {
            const label fI = pointFaces(bpI, pfI);

            if (faceGroup.found(fI))
                continue;

            front.clear();
            front.append(fI);
            faceGroup.insert(fI, nGroup);

            while (front.size() != 0)
            {
                const label fLabel = front.remove();

                forAllRow(faceFaces, fLabel, ffI)
                {
                    const label neiFace = faceFaces(fLabel, ffI);

                    if (faceGroup.found(neiFace))
                        continue;

                    if (pointFaces.found(bpI, neiFace))
                    {
                        front.append(neiFace);
                        faceGroup.insert(neiFace, nGroup);
                    }
                }
            }

            ++nGroup;
        }

        if (nGroup != 1)
        {
            const label pointI = bPoints[bpI];
            # pragma omp critical
            {
                badVertices.insert(pointI);
                invalidVertices.insert(pointI);
            }
        }
    }

    const label nBadVertices = invalidVertices.size();
    Info<< "Found " << nBadVertices << " invalid connections" << endl;
    Info<< "Finished checking faces connections to surface vertices" << endl;

    if (nBadVertices != 0 && removeCells)
    {
        const VRWGraph& pointCells = mesh_.addressingData().pointCells();

        boolList removeCell(mesh_.cells().size(), false);

        for (const label verti : invalidVertices)
        {
            forAllRow(pointCells, verti, pcI)
            {
                removeCell[pointCells(verti, pcI)] = true;
            }
        }

        polyMeshGenModifier(mesh_).removeCells(removeCell);
        clearMeshEngine();

        return true;
    }

    return false;
}


// ************************************************************************* //
