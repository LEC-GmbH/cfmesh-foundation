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
    Foam::Module::meshGenerator

Description
    Runtime-selectable base class for cfMesh's mesh generators
    (cartesianMeshGenerator, cartesian2DMeshGenerator, tetMeshGenerator,
    voronoiMeshGenerator), selected by mode name (cartesian, cartesian2D,
    tetrahedral, polyhedral).

SourceFiles
    meshGenerator.C
    meshGeneratorNew.C

\*---------------------------------------------------------------------------*/

#ifndef meshGenerator_H
#define meshGenerator_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class Time;
class IOdictionary;
class IOobject;
class polyMesh;

namespace Module
{

/*---------------------------------------------------------------------------*\
                      Class meshGenerator Declaration
\*---------------------------------------------------------------------------*/

class meshGenerator
{
public:

    //- Runtime type information
    TypeName("meshGenerator");


    // Constructors

        //- Construct null
        meshGenerator() = default;

        //- Disallow default bitwise copy construct
        meshGenerator(const meshGenerator&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const meshGenerator&) = delete;


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        meshGenerator,
        dictionary,
        (
            const Time& time,
            IOdictionary& meshDict,
            const word& regionName,
            const fileName& meshPath
        ),
        (time, meshDict, regionName, meshPath)
    );

    //- Select constructed from Time and meshDict, using the given mode name -
    //- regionName/meshPath determine where the generated mesh is read
    //- from/written to at each intermediate stage (the final output
    //- location is set separately, via generatePolyMesh()'s IOobject)
    static autoPtr<meshGenerator> New
    (
        const word& modeName,
        const Time& time,
        IOdictionary& meshDict,
        const word& regionName,
        const fileName& meshPath
    );

    //- Destructor
    virtual ~meshGenerator() = default;


    // Member Functions

        //- generate the initial (castellated) mesh from the octree -
        //- the caller decides whether to call this or readMesh() (see
        //- meshDict's castellatedMesh switch), matching snappyHexMesh's
        //- own castellatedMesh/snap/addLayers stage control
        virtual void createCastellatedMesh() = 0;

        //- read an existing mesh from disk instead of generating one
        virtual void readMesh() = 0;

        //- fit the mesh boundary onto the input surface
        virtual void snap() = 0;

        //- add boundary layers
        virtual void addLayers() = 0;

        //- final mesh-quality optimisation (surface untangling, mesh
        //- smoothing, etc.) - not itself a layer-adding operation, so the
        //- caller runs this even when addLayers is disabled, instead of
        //- silently skipping mesh cleanup along with the layers.
        //- addLayers() also calls this itself as part of its own sequence
        //- when layers are actually added, in the same relative order as
        //- before this method existed. layersExist must reflect whether
        //- boundary layer cells actually exist in the mesh right now -
        //- some of the optimisation performed here assumes they do (e.g.
        //- refining/smoothing the layer itself) and will fatal if asked
        //- to run against a mesh that has none
        virtual void optimiseMesh(const bool layersExist) = 0;

        //- renumber the mesh and apply boundary renaming - always run,
        //- regardless of which stages above were selected
        virtual void finalise() = 0;

        //- convert the generated mesh into a Foam::polyMesh identified by
        //- the given IOobject (region/instance/local all taken from it) -
        //- the caller is responsible for calling write() on the result
        virtual autoPtr<polyMesh> generatePolyMesh(const IOobject& io)
            const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
