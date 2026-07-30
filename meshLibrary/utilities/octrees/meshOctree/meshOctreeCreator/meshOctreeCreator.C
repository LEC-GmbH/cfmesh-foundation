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

#include "meshOctreeCreator.H"
#include "triSurf.H"
#include "boundBox.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshOctreeCreator::meshOctreeCreator
(
    meshOctree& mo,
    const dictionary& dict
)
:
    octree_(mo),
    scalingFactor_(1.0),
    meshDict_(dict),
    haveMeshDict_(&dict != &dictionary::null),
    maxCellSize_(haveMeshDict_ ? readScalar(dict.lookup("maxCellSize")) : 0),
    haveBoundaryCellSize_(dict.found("boundaryCellSize")),
    boundaryCellSize_
    (
        haveBoundaryCellSize_ ? readScalar(dict.lookup("boundaryCellSize")) : 0
    ),
    boundaryCellSizeRefinementThickness_
    (
        haveBoundaryCellSize_
      ? dict.lookupOrDefault<scalar>("boundaryCellSizeRefinementThickness", 0)
      : 0
    ),
    patchCellSizeDict_(dict.subOrEmptyDict("patchCellSize")),
    subsetCellSizeDict_(dict.subOrEmptyDict("subsetCellSize")),
    localRefinementDict_(dict.subOrEmptyDict("localRefinement")),
    objectRefinementsDict_(dict.subOrEmptyDict("objectRefinements")),
    surfaceMeshRefinementDict_(dict.subOrEmptyDict("surfaceMeshRefinement")),
    edgeMeshRefinementDict_(dict.subOrEmptyDict("edgeMeshRefinement")),
    anisotropicSourcesDict_(dict.subOrEmptyDict("anisotropicSources")),
    hexRefinement_(false),
    globalRefLevel_(0),
    surfRefLevel_(mo.surface().size())
{
    if (haveMeshDict_)
    {
        if (maxCellSize_ < 0)
        {
            FatalErrorInFunction
                << "maxCellSize is negative! Cannot generate the mesh!!"
                << exit(FatalError);
        }

        if (haveBoundaryCellSize_)
        {
            if (boundaryCellSize_ < 0)
            {
                WarningInFunction
                    << "Boundary cell size is negative!!" << endl;
            }

            if (boundaryCellSizeRefinementThickness_ < 0)
            {
                WarningInFunction
                    << "Boundary cell size refinement thickness is"
                    << " negative!!" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeCreator::activateHexRefinement()
{
    hexRefinement_ = true;
}


// ************************************************************************* //
