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
    Converts a cfMesh-internal Foam::Module::polyMeshGen into a standalone
    Foam::polyMesh, so that the standard OpenFOAM IOobject/region/instance
    machinery (rather than cfMesh's own hardcoded constant/polyMesh writer)
    determines where the mesh is written from - the same mechanism
    blockMesh itself uses for its own output.

SourceFiles
    polyMeshGenToFoamMesh.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenToFoamMesh_H
#define polyMeshGenToFoamMesh_H

#include "autoPtr.H"
#include "polyMesh.H"

namespace Foam
{

class IOobject;

namespace Module
{

class polyMeshGen;

//- Convert a finished polyMeshGen into a Foam::polyMesh identified by the
//- given IOobject (region/instance/local all taken from io). Point/face/
//- cell subsets are intentionally not transferred - checkMesh already
//- reports bad-quality cells/faces, and OpenFOAM zones are the modern
//- equivalent of subsets.
autoPtr<polyMesh> polyMeshGenToFoamMesh
(
    const polyMeshGen& meshGen,
    const IOobject& io
);

} // End namespace Module
} // End namespace Foam

#endif

// ************************************************************************* //
