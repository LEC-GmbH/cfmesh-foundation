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

#include "topologicalCleaner.H"
#include "polyMeshGenAddressing.H"
#include "DynList.H"
#include "meshSurfaceEngine.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::topologicalCleaner::checkInvalidConnectionsForVerticesCells
(
    labelHashSet* irregularNodesPtr
)
{
    polyMeshGenModifier meshModifier(mesh_);

    label nPoints = mesh_.points().size();
    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();
    const cellListPMG& cells = mesh_.cells();

    const VRWGraph& pointCells = mesh_.addressingData().pointCells();
    const VRWGraph& cellCells = mesh_.addressingData().cellCells();

    meshSurfaceEngine mse(mesh_);
    const labelList& bPoints = mse.boundaryPoints();

    label nInvalidConnections(0);

    forAll(bPoints, bpI)
    {
        const label pointI = bPoints[bpI];

        label material(1);

        labelList materialForCell(pointCells.sizeOfRow(pointI), 0);

        forAllRow(pointCells, pointI, cI)
        if (!materialForCell[cI])
        {
            materialForCell[cI] = material;

            DynList<label> frontCells;
            frontCells.append(cI);

            do
            {
                DynList<label> newFrontCells;

                forAll(frontCells, fcI)
                {
                    const label pointCellI =
                    pointCells(pointI, frontCells[fcI]);

                    forAllRow(cellCells, pointCellI, nI)
                    {
                        forAllRow(pointCells, pointI, pcI)
                        {
                            if (materialForCell[pcI])
                                continue;

                            if (
                                cellCells(pointCellI, nI) ==
                                pointCells(pointI, pcI)
                            )
                            {
                                newFrontCells.append(pcI);
                                materialForCell[pcI] = material;
                                break;
                            }
                        }
                    }
                }

                frontCells = newFrontCells;

            } while (frontCells.size() != 0);

            ++material;
        }

        if (material > 2)
        {
            ++nInvalidConnections;

            if (irregularNodesPtr)
            {
                irregularNodesPtr->insert(pointI);
                continue;
            }

            forAllRow(pointCells, pointI, pcI)
            {
                if (materialForCell[pcI] == 1)
                    continue;

                const cell& c = cells[pointCells(pointI, pcI)];

                forAll(c, fI)
                {
                    face& f = faces[c[fI]];

                    forAll(f, pI)
                        if (f[pI] == pointI)
                            f[pI] = nPoints + materialForCell[pcI] - 1;
                }
            }

            for (label i = 1; i < material; ++i)
            {
                const point p = points[pointI];
                points.append(p);
                ++nPoints;
            }
        }
    }

    Info<< "Found " << nInvalidConnections
        << " invalid cell connections" << endl;

    mesh_.clearAddressingData();

    if (nInvalidConnections != 0)
        meshModifier.removeUnusedVertices();
}


void Foam::Module::topologicalCleaner::checkInvalidConnectionsForVerticesFaces
(
    labelHashSet* /*irregularNodesPtr*/
)
{
    const meshSurfaceEngine mse(mesh_);

    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& faceOwner = mse.faceOwners();

    boolList removeCell(mesh_.cells().size(), false);
    bool changed(false);

    # pragma omp parallel for schedule(static, 1)
    forAll(edgeFaces, edgeI)
        if (edgeFaces.sizeOfRow(edgeI) > 2)
        {
            forAllRow(edgeFaces, edgeI, fI)
                removeCell[faceOwner[edgeFaces(edgeI, fI)]] = true;

            changed = true;
        }

    if (changed)
    {
        polyMeshGenModifier(mesh_).removeCells(removeCell);

        decomposeCell_.setSize(mesh_.cells().size());
        decomposeCell_ = false;
    }
}


// ************************************************************************* //
