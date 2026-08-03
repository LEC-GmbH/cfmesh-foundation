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
        //- smoothing, etc.) - called by snap() right after it inserts the
        //- single geometry-adapted buffer boundary layer, so snap()'s own
        //- output is already good quality (matching snappyHexMesh's own
        //- snap/addLayers split, where a stage's written mesh is always
        //- fully valid); addLayers() only calls this itself if snap()
        //- didn't already (e.g. a staged workflow reading an unoptimised
        //- buffer layer back from an earlier invocation). layersExist
        //- must reflect whether boundary layer cells actually exist in
        //- the mesh right now - some of the optimisation performed here
        //- assumes they do (e.g. refining/smoothing the layer itself) and
        //- will fatal if asked to run against a mesh that has none
        virtual void optimiseMesh(const bool layersExist) = 0;

        //- renumber the mesh - the caller runs this before every write, no
        //- matter which stage it is: every written mesh (e.g. an
        //- intermediate -noOverwrite stage) is a complete, real mesh
        //- someone may inspect or use directly, not just the final one, so
        //- every write should be well-formed, not just the last
        virtual void renumberMesh() = 0;

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
