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

#include "meshSurfaceMapper2D.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "polyMeshGen2DEngine.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "demandDrivenData.H"
#include "meshOctree.H"
#include "helperFunctions.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceMapper2D::findActiveBoundaryEdges()
{
    const labelListList2D& edgeFaces = surfaceEngine_.edgeFaces();

    const polyMeshGen2DEngine& mesh2DEngine = this->mesh2DEngine();
    const boolList& activeFace = mesh2DEngine.activeFace();

    const label startFace = surfaceEngine_.mesh().boundaries()[0].patchStart();

    activeBoundaryEdges_.clear();

    // check which edges are at the boundary
    forAll(edgeFaces, edgeI)
    {
        if (edgeFaces[edgeI].size() == 2)
        {
            const bool active0 = activeFace[startFace + edgeFaces[edgeI][0]];
            const bool active1 = activeFace[startFace + edgeFaces[edgeI][1]];

            if (active0 && active1)
                activeBoundaryEdges_.append(edgeI);
        }
    }
}


void Foam::Module::meshSurfaceMapper2D::create2DEngine() const
{
    polyMeshGen& mesh = const_cast<polyMeshGen&>(surfaceEngine_.mesh());
    mesh2DEnginePtr_ = new polyMeshGen2DEngine(mesh);
}


void Foam::Module::meshSurfaceMapper2D::createTriSurfacePartitioner() const
{
    surfPartitionerPtr_ = new triSurfacePartitioner(meshOctree_.surface());
}


void Foam::Module::meshSurfaceMapper2D::createMeshSurfacePartitioner() const
{
    meshPartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine_);
}


void Foam::Module::meshSurfaceMapper2D::clearOut()
{
    deleteDemandDrivenData(mesh2DEnginePtr_);
    deleteDemandDrivenData(surfPartitionerPtr_);
    deleteDemandDrivenData(meshPartitionerPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshSurfaceMapper2D::meshSurfaceMapper2D
(
    const meshSurfaceEngine& mse,
    const meshOctree& octree
)
:
    surfaceEngine_(mse),
    meshOctree_(octree),
    mesh2DEnginePtr_(nullptr),
    surfPartitionerPtr_(nullptr),
    meshPartitionerPtr_(nullptr)
{
    findActiveBoundaryEdges();

    createMeshSurfacePartitioner();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::meshSurfaceMapper2D::~meshSurfaceMapper2D()
{
    clearOut();
}


// ************************************************************************* //
