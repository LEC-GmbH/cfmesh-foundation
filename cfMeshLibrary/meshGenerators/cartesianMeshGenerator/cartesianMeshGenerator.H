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
    Foam::Module::cartesianMeshGenerator

Description
    Creates cartesian mesh from the octree

SourceFiles
    cartesianMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianMeshGenerator_H
#define cartesianMeshGenerator_H

#include "polyMeshGen.H"
#include "IOdictionary.H"
#include "meshGenerator.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class Time;

namespace Module
{
class meshOctree;
class triSurf;

/*---------------------------------------------------------------------------*\
                   Class cartesianMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class cartesianMeshGenerator
:
    public meshGenerator
{
    // Private data

        //- reference to Time
        const Time& db_;

        //- pointer to the surface
        const triSurf* surfacePtr_;

        //- pointer to the modified surface
        const triSurf* modSurfacePtr_;

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary& meshDict_;

        //- pointer to the octree
        meshOctree* octreePtr_;

        //- mesh
        polyMeshGen mesh_;


    // Private member functions

        //- create cartesian mesh
        void createCartesianMesh();

        //- prepare mesh surface
        void surfacePreparation();

        //- map mesh to the surface and untangle surface
        void mapMeshToSurface();

        //- capture edges and corners
        void extractPatches();

        //- map points at feature edges and corners
        void mapEdgesAndCorners();

        //- optimise surface mesh
        void optimiseMeshSurface();

        //- add boundary layers
        void generateBoundaryLayers();

        //- refine boundary layers
        void refBoundaryLayers();

        //- mesh optimisation. layersExist selects whether boundary-layer
        //- specific optimisation also runs (see meshGenerator::
        //- optimiseMesh())
        void optimiseFinalMesh(const bool layersExist);

        //- re-project points back on the surface mesh after back-scaling
        void projectSurfaceAfterBackScaling(const bool layersExist);

        //- collapse feature-edge-split boundary patches back to their
        //- pre-split name/type
        void replaceBoundaries();

        //- apply the user's renameBoundary rename/merge, if configured
        void renamePatches();

        //- Disallow default bitwise copy construct
        cartesianMeshGenerator(const cartesianMeshGenerator&);

        //- Disallow default bitwise assignment
        void operator=(const cartesianMeshGenerator&);

public:

    //- Runtime type information
    TypeName("cartesian");

    //- Construct from time, meshDict, and the region/mesh-path this run
    //- should read/write its intermediate mesh state at
    cartesianMeshGenerator
    (
        const Time&,
        IOdictionary& meshDict,
        const word& regionName,
        const fileName& meshPath
    );

    //- Destructor
    ~cartesianMeshGenerator();


    // Member Functions

        //- generate the initial (castellated) mesh from the octree
        void createCastellatedMesh() override;

        //- read an existing mesh from disk instead of generating one
        void readMesh() override;

        //- fit the mesh boundary onto the input surface
        void snap() override;

        //- add boundary layers
        void addLayers() override;

        //- final mesh-quality optimisation
        void optimiseMesh(const bool layersExist) override;

        //- renumber the mesh
        void renumberMesh() override;

        //- convert the generated mesh into a Foam::polyMesh
        autoPtr<polyMesh> generatePolyMesh(const IOobject& io) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
