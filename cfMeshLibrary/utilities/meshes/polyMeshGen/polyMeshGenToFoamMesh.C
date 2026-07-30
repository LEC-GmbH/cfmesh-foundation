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

#include "polyMeshGenToFoamMesh.H"
#include "polyMeshGen.H"
#include "boundaryPatch.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * //

namespace
{

//- Copy the "used" portion of a cfMesh PMG-style container (its size()
//- may be smaller than the underlying List's allocated capacity)
template<class T, class PMGList>
Foam::List<T> copyUsedElements(const PMGList& l)
{
    return Foam::List<T>(Foam::SubList<T>(l, l.size()));
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMesh> Foam::Module::polyMeshGenToFoamMesh
(
    const polyMeshGen& meshGen,
    const IOobject& io
)
{
    pointField points(copyUsedElements<point>(meshGen.points()));
    faceList faces(copyUsedElements<face>(meshGen.faces()));
    labelList own(meshGen.owner());
    labelList nei(meshGen.neighbour());

    autoPtr<polyMesh> meshPtr
    (
        new polyMesh
        (
            io,
            std::move(points),
            std::move(faces),
            std::move(own),
            std::move(nei)
        )
    );

    polyMesh& mesh = meshPtr();

    const PtrList<boundaryPatch>& boundaries = meshGen.boundaries();

    List<polyPatch*> patches(boundaries.size());

    forAll(boundaries, patchI)
    {
        patches[patchI] = polyPatch::New
        (
            boundaries[patchI].patchType(),
            boundaries[patchI].patchName(),
            boundaries[patchI].patchSize(),
            boundaries[patchI].patchStart(),
            patchI,
            mesh.boundary()
        ).ptr();
    }

    mesh.addPatches(patches);

    return meshPtr;
}

// ************************************************************************* //
