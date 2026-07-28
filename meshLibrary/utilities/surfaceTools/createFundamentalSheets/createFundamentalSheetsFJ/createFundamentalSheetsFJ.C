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

#include "createFundamentalSheetsFJ.H"
#include "demandDrivenData.H"
#include "meshSurfaceEngine.H"
#include "extrudeLayer.H"

#include "addToRunTimeSelectionTable.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
    defineTypeNameAndDebug(createFundamentalSheetsFJ, 0);
    addToRunTimeSelectionTable
    (
        createFundamentalSheets,
        createFundamentalSheetsFJ,
        polyMeshGen
    );
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::createFundamentalSheetsFJ::createInitialSheet()
{
    if (!createWrapperSheet_)
    {
        return;
    }

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    const label start = boundaries[0].patchStart();
    const label end
    (
        boundaries[boundaries.size()-1].patchStart() +
        boundaries[boundaries.size()-1].patchSize()
    );

    faceListPMG::subList bFaces(mesh_.faces(), end - start, start);

    const labelList& owner = mesh_.owner();

    DynamicList<labelPair> extrudeFaces;
    extrudeFaces.setSize(end - start);

    # pragma omp parallel for
    for (label faceI = start; faceI < end; ++faceI)
    {
        extrudeFaces[faceI - start] = labelPair(faceI, owner[faceI]);
    }

    extrudeLayer(mesh_, extrudeFaces);
}


void Foam::Module::createFundamentalSheetsFJ::createSheetsAtFeatureEdges()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(boundaries, patchI)
    {
        const cellListPMG& cells = mesh_.cells();
        boolList patchCell(cells.size(), false);

        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for (label faceI = start; faceI < end; ++faceI)
        {
            patchCell[owner[faceI]] = true;
        }

        DynamicList<labelPair> front;

        for (label faceI = start; faceI < end; ++faceI)
        {
            const cell& c = cells[owner[faceI]];

            forAll(c, fI)
            {
                if (neighbour[c[fI]] < 0)
                {
                    continue;
                }

                label nei = owner[c[fI]];
                if (nei == owner[faceI])
                {
                    nei = neighbour[c[fI]];
                }

                if (!patchCell[nei])
                {
                    front.append(labelPair(c[fI], nei));
                }
            }
        }

        extrudeLayer(mesh_, front);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::createFundamentalSheetsFJ::createFundamentalSheetsFJ
(
    polyMeshGen& mesh,
    const bool createWrapperSheet
)
:
    createFundamentalSheets(mesh, createWrapperSheet)
{
    createInitialSheet();

    createSheetsAtFeatureEdges();
}


// ************************************************************************* //
