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

#include "decomposeCells.H"
#include "helperFunctions.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::decomposeCells::findAddressingForCell
(
    const label cellI,
    DynamicList<label>& vrt,
    DynamicList<edge>& edges,
    labelListList2D& faceEdges,
    DynamicList<DynamicList<label>>& edgeFaces
) const
{
    const cell& c = mesh_.cells()[cellI];

    vrt.clear();
    edges.clear();
    edgeFaces.clear();
    faceEdges.setSize(c.size());

    const faceListPMG& faces = mesh_.faces();
    forAll(faceEdges, feI)
    {
        DynamicList<label>& fEdges = faceEdges[feI];

        fEdges.setSize(faces[c[feI]].size());
        fEdges = -1;
    }

    forAll(c, fI)
    {
        const face& f = faces[c[fI]];

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            bool store(true);
            forAll(vrt, vI)
            {
                if (vrt[vI] == f[eI])
                {
                    store = false;
                    break;
                }
            }
            if (store)
            {
                vrt.append(f[eI]);
            }

            // check if the edge alreready exists
            store = true;

            forAll(edges, eJ)
            {
                if (e == edges[eJ])
                {
                    store = false;
                    faceEdges[fI][eI] = eJ;
                    edgeFaces[eJ].append(fI);
                    break;
                }
            }

            if (store)
            {
                faceEdges[fI][eI] = edges.size();
                DynamicList<label> ef;
                ef.append(fI);
                edgeFaces.append(ef);
                edges.append(e);
            }
        }
    }

    // check if the cell is topologically closed
    forAll(edgeFaces, efI)
    {
        if (edgeFaces[efI].size() != 2)
        {
            forAll(c, fI)
            {
                Info<< "Face " << c[fI] << " is " << faces[c[fI]] << endl;
            }

            Info<< "Edges " << edges << endl;
            Info<< "faceEdges " << faceEdges << endl;
            Info<< "edgeFaces " << edgeFaces << endl;
            mesh_.write();

            FatalErrorInFunction
                << "Cell " << cellI << " is not topologically closed!"
                << abort(FatalError);
        }
    }
}


Foam::label Foam::Module::decomposeCells::findTopVertex
(
    const label cellI,
    const DynamicList<label>& /*vrt*/,
    const DynamicList<edge>& /*edges*/,
    const DynamicList<DynamicList<label>>& /*edgeFaces*/
)
{
    const cell& c = mesh_.cells()[cellI];
    const faceListPMG& faces = mesh_.faces();

    pointFieldPMG& pointsAccess = mesh_.points();

    // there is no vertex in 3 or more patches
    // find boundary faces
    label topVertex(-1);

    const labelList cp = c.labels(faces);
    point p(vector::zero);
    forAll(cp, cpI)
    {
        p += pointsAccess[cp[cpI]];
    }
    p /= cp.size();

    topVertex = pointsAccess.size();
    pointsAccess.append(p);

    return topVertex;
}


void Foam::Module::decomposeCells::decomposeCellIntoPyramids(const label cellI)
{
    const cellListPMG& cells = mesh_.cells();
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();

    const cell& c = cells[cellI];

    // calculate edges, faceEdges and edgeFaces addressing
    DynamicList<label> vrt;
    DynamicList<edge> edges;
    labelListList2D faceEdges;
    faceEdges.setSize(c.size());
    DynamicList<DynamicList<label>> edgeFaces;
    findAddressingForCell(cellI, vrt, edges, faceEdges, edgeFaces);

    // find a vertex which will be the top of the pyramids
    // if there exist a corner vertex which is in 3 or more patches then
    // it is selected as the top vertex
    label topVertex = findTopVertex(cellI, vrt, edges, edgeFaces);

    // start generating pyramids
    forAll(c, fI)
    {
        const face& f = faces[c[fI]];
        faceList cellFaces(f.size() + 1);

        DynamicList<triFace> triFaces;
        triFaces.setSize(f.size());
        forAll(triFaces, pI)
        {
            triFaces[pI][0] = f.nextLabel(pI);
            triFaces[pI][1] = f[pI];
            triFaces[pI][2] = topVertex;
        }

        label cfI(0);
        if (owner[c[fI]] == cellI)
        {
            cellFaces[cfI++] = faces[c[fI]];

            forAll(triFaces, tfI)
            {
                cellFaces[cfI++] = triFaces[tfI];
            }
        }
        else
        {
            cellFaces[cfI++] = faces[c[fI]].reverseFace();

            forAll(triFaces, tfI)
            {
                triFace rf;
                rf[0] = triFaces[tfI][0];
                rf[1] = triFaces[tfI][2];
                rf[2] = triFaces[tfI][1];
                cellFaces[cfI++] = rf;
            }
        }

        facesOfNewCells_.append(cellFaces);
    }
}


// ************************************************************************* //
