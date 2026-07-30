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

\*---------------------------------------------------------------------------*/

#include "voronoiMeshExtractor.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::voronoiMeshExtractor::createPoints()
{
    const DynamicList<point>& tetPoints = tetCreator_.tetPoints();
    const DynamicList<partTet>& tets = tetCreator_.tets();

    pointFieldPMG& points = mesh_.points();
    points.setSize(tets.size());

    # pragma omp parallel for schedule(dynamic, 100)
    forAll(tets, tetI)
    {
        points[tetI] = tets[tetI].centroid(tetPoints);

    }
}

void Foam::Module::voronoiMeshExtractor::createPolyMesh()

{
    const labelListList2D& pointEdges = this->pointEdges();
    const DynamicList<DynamicList<label>>& edgeTets = this->edgeTets();
    const boolList& boundaryEdge = this->boundaryEdge();
    const DynamicList<edge>& edges = this->edges();

    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    // count the number of cells
    label nCells(0);
    labelList cellLabel(pointEdges.size(), -1);

    forAll(pointEdges, pointI)
    {
        bool create(true);
        forAllRow(pointEdges, pointI, eI)
        {
            if (boundaryEdge[pointEdges[pointI][eI]])
            {
                create = false;
                break;
            }
        }

        if (!create || (pointEdges[pointI].size() == 0))
        {
            continue;
        }

        cellLabel[pointI] = nCells;
        ++nCells;
    }

    // count the number of faces
    label nFaces(0);
    labelList faceLabel(edges.size(), -1);

    forAll(boundaryEdge, edgeI)
    {
        if (boundaryEdge[edgeI])
        {
            continue;
        }

        const edge& e = edges[edgeI];
        if (cellLabel[e[0]] < 0 && cellLabel[e[1]] < 0)
        {
            continue;
        }

        faceLabel[edgeI] = nFaces;
        ++nFaces;
    }

    // create faces
    Info<< "Creating faces " << endl;
    faces.setSize(nFaces);
    labelList nFacesInCell(nCells, 0);

    forAll(edges, edgeI)
    {
        if (faceLabel[edgeI] < 0)
        {
            continue;
        }

        const label faceI = faceLabel[edgeI];

        face& f = faces[faceI];
        f.setSize(edgeTets[edgeI].size());

        // fill the faces with the node labels
        forAllRow(edgeTets, edgeI, pI)
        {
            f[pI] = edgeTets[edgeI][pI];
        }

        const edge& e = edges[edgeI];
        const label cOwn = cellLabel[e[0]];
        const label cNei = cellLabel[e[1]];

        if (cOwn < 0 || ((cOwn > cNei) && (cNei != -1)))
        {

            f = f.reverseFace();
        }

        if (cOwn >= 0)
        {
            ++nFacesInCell[cOwn];
        }
        if (cNei >= 0)
        {
            ++nFacesInCell[cNei];
        }
    }

    // create cells

    cells.setSize(nCells);
    forAll(nFacesInCell, cellI)
    {
        cells[cellI].setSize(nFacesInCell[cellI]);
        nFacesInCell[cellI] = 0;
    }

    forAll(edges, edgeI)
    {
        if (faceLabel[edgeI] < 0)
        {
            continue;
        }

        const label faceI = faceLabel[edgeI];
        const edge& e = edges[edgeI];
        const label cOwn = cellLabel[e[0]];
        const label cNei = cellLabel[e[1]];

        if (cOwn >= 0)
        {
            cells[cOwn][nFacesInCell[cOwn]++] = faceI;
        }
        if (cNei >= 0)
        {
            cells[cNei][nFacesInCell[cNei]++] = faceI;
        }
    }

    mesh_.clearAddressingData();
}


// ************************************************************************* //
