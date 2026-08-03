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

Class
    Foam::Module::renameBoundaryPatches

Description
    Collapses any feature-edge-split boundary patch back to its one
    pre-split name/type (a mesh never exposes a split sub-name like
    "patch1_3" in its written output), fixes empty patches, and
    flattens symmetry planes. If constructed with a meshDict, also
    renames/merges patches according to the optional renameBoundary
    subDict (newPatchNames/defaultName/defaultType, mirroring
    snappyHexMesh's own rename/merge model).

    Callers apply these two constructors at different pipeline points:
    the no-dict overload right after patches are first assigned (so any
    later per-patch-name-keyed config, e.g. boundaryLayers.
    patchBoundaryLayers, still resolves against stable, un-merged
    names), and the meshDict overload only once no further per-patch
    config remains to be resolved - merging patches together any
    earlier would erase the very per-region distinction such config
    depends on.

SourceFiles
    renameBoundaryPatches.C

\*---------------------------------------------------------------------------*/

#ifndef renameBoundaryPatches_H
#define renameBoundaryPatches_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class IOdictionary;

namespace Module
{

/*---------------------------------------------------------------------------*\
                    Class renameBoundaryPatches Declaration
\*---------------------------------------------------------------------------*/

class renameBoundaryPatches
{
    // Private data

        //- mesh
        polyMeshGen& mesh_;

        //- renameBoundary subDict, or dictionary::null if not present;
        //- derived once at construction via subOrEmptyDict()
        const dictionary& renameBoundaryDict_;


    // Private member functions

        //- rename/merge patches per renameBoundaryDict_ and collapse
        //- any remaining feature-edge-split patch back to its
        //- pre-split name/type; replace the existing boundary
        void calculateNewBoundary();

        //- find empty patches in 3D meshes and change them to wall
        void checkEmptyPatches();

        //- check if there exist any symmetry planes
        //- and make sure they get perfectly flat
        void checkSymmetryPlanes();

        //- Disallow default bitwise copy construct
        renameBoundaryPatches(const renameBoundaryPatches&);

        //- Disallow default bitwise assignment
        void operator=(const renameBoundaryPatches&);


public:

    //- Construct from mesh; collapses feature-edge-split patches back
    //- to their pre-split name/type only - no user-configured
    //- rename/merge
    renameBoundaryPatches
    (
        polyMeshGen& mesh,
        const bool allowEmptyPatches = false
    );

    //- Construct from mesh and IOdictionary; additionally applies the
    //- optional renameBoundary rename/merge
    renameBoundaryPatches
    (
        polyMeshGen& mesh,
        const IOdictionary& meshDict,
        const bool allowEmptyPatches = false
    );

    //- Destructor
    ~renameBoundaryPatches() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
