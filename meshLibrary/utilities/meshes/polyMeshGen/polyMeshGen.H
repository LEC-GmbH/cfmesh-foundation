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
    Foam::Module::polyMeshGen

Description
    Mesh with selections

SourceFiles
    polyMeshGen.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGen_H
#define polyMeshGen_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenCells.H"

namespace Foam
{
namespace Module
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGen
:
    public polyMeshGenCells
{
public:

    // Constructors

        //- Null constructor
        polyMeshGen(const Time& t);

        //- Construct from components without the boundary
        polyMeshGen
        (
            const Time& t,
            const pointField& points,
            const faceList& faces,
            const cellList& cells
        );

        //- Construct from components with the boundary
        polyMeshGen
        (
            const Time& t,
            const pointField& points,
            const faceList& faces,
            const cellList& cells,
            const wordList& patchNames,
            const labelList& patchStart,
            const labelList& nFacesInPatch
        );


    //- Destructor
    ~polyMeshGen() = default;


    // Public member functions

        //- set the region/mesh-path this mesh is read from/written to -
        //- must be called (if at all) before any read()/write()
        void setMeshLocation
        (
            const word& regionName,
            const fileName& meshPath
        );

        //- Read mesh
        void read();

        //- Write mesh
        void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
