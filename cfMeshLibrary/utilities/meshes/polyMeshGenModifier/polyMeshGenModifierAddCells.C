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

#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::polyMeshGenModifier::addCells
(
    const DynamicList<faceList>& cellFaces
)
{
    Info<< "Adding cells to the mesh" << endl;

    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    labelListList2D& pointFaces = this->pointFaces();

    // start adding cells into the mesh
    label nFaces = faces.size();
    label nCells = cells.size();

    forAll(cellFaces, cI)
    {
        const faceList& facesInCell = cellFaces[cI];

        label fI(0);
        cell c(cellFaces[cI].size());

        forAll(facesInCell, faceI)
        {
            const face& f = facesInCell[faceI];

            const label pointI = f[0];

            label fLabel(-1);
            forAllRow(pointFaces, pointI, pfI)
            {
                const label faceI = pointFaces[pointI][pfI];

                if (faces[faceI] == f)
                {
                    fLabel = faceI;
                    break;
                }
            }

            if (fLabel == -1)
            {
                faces.append(f);
                c[fI++] = nFaces;

                forAll(f, pI)
                {
                    pointFaces[f[pI]].append(nFaces);
                }

                ++nFaces;
            }
            else
            {
                c[fI++] = fLabel;
            }
        }

        cells.append(c);
        ++nCells;
    }

    this->clearOut();
    mesh_.clearOut();

    Info<< "Finished adding cells to the mesh" << endl;
}


void Foam::Module::polyMeshGenModifier::addCell(const faceList& cellFaces)
{
    faceListPMG& faces = this->facesAccess();
    cellListPMG& cells = this->cellsAccess();

    label nFaces = faces.size();

    labelListList2D& pointFaces = this->pointFaces();

    cell c(cellFaces.size());
    label fI(0);

    forAll(cellFaces, faceI)
    {
        const face& f = cellFaces[faceI];

        const label pointI = f[0];

        label fLabel(-1);
        forAllRow(pointFaces, pointI, pfI)
        {
            const label faceI = pointFaces[pointI][pfI];

            if (faces[faceI] == f)
            {
                fLabel = faceI;
                break;
            }
        }

        if (fLabel == -1)
        {
            faces.append(f);
            c[fI++] = nFaces;

            forAll(f, pI)
            {
                pointFaces[f[pI]].append(nFaces);
            }

            ++nFaces;
        }
        else
        {
            c[fI++] = fLabel;
        }
    }

    cells.append(c);
    mesh_.clearOut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
