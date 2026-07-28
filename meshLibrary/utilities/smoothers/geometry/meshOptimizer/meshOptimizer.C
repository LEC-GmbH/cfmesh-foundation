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

#include "demandDrivenData.H"
#include "meshOptimizer.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGenChecks.H"
#include "DynamicList.H"

// * * * * * * * * Private member functions * * * * * * * * * * * * * * * * * //

const Foam::Module::meshSurfaceEngine& 
Foam::Module::meshOptimizer::meshSurface() const
{
    if (!msePtr_)
    {
        msePtr_ = new meshSurfaceEngine(mesh_);
    }

    return *msePtr_;
}


void Foam::Module::meshOptimizer::clearSurface()
{
    deleteDemandDrivenData(msePtr_);
}


Foam::label Foam::Module::meshOptimizer::findBadFaces
(
    labelHashSet& badFaces,
    const boolList& changedFace
) const
{
    badFaces.clear();

    polyMeshGenChecks::checkFacePyramids
    (
        mesh_,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkFaceFlatness
    (
        mesh_,
        false,
        0.8,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkCellPartTetrahedra
    (
        mesh_,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    polyMeshGenChecks::checkFaceAreas
    (
        mesh_,
        false,
        VSMALL,
        &badFaces,
        &changedFace
    );

    const label nBadFaces = badFaces.size();

    return nBadFaces;
}


Foam::label Foam::Module::meshOptimizer::findLowQualityFaces
(
    labelHashSet& badFaces,
    const boolList& /*changedFace*/
) const
{
    badFaces.clear();

    polyMeshGenChecks::checkFaceDotProduct
    (
        mesh_,
        false,
        70.0,
        &badFaces
    );

    polyMeshGenChecks::checkFaceSkewness
    (
        mesh_,
        false,
        2.0,
        &badFaces
    );

    const label nBadFaces = badFaces.size();

    return nBadFaces;
}


void Foam::Module::meshOptimizer::calculatePointLocations()
{
    vertexLocation_.setSize(mesh_.points().size());
    vertexLocation_ = INSIDE;

    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bPoints = mse.boundaryPoints();

    // mark boundary vertices
    forAll(bPoints, bpI)
    {
        vertexLocation_[bPoints[bpI]] = BOUNDARY;
    }

    // mark edge vertices
    meshSurfacePartitioner mPart(mse);
    for (const label edgei : mPart.edgePoints().toc())
    {
        vertexLocation_[bPoints[edgei]] = EDGE;
    }

    // mark corner vertices
    for (const label corneri : mPart.corners().toc())
    {
        vertexLocation_[bPoints[corneri]] = CORNER;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshOptimizer::meshOptimizer(polyMeshGen& mesh)
:
    mesh_(mesh),
    vertexLocation_(),
    lockedFaces_(),
    msePtr_(nullptr),
    enforceConstraints_(false),
    badPointsSubsetName_()
{
    calculatePointLocations();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::meshOptimizer::~meshOptimizer()
{
    clearSurface();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOptimizer::enforceConstraints(const word subsetName)
{
    enforceConstraints_ = true;

    badPointsSubsetName_ = subsetName;
}


void Foam::Module::meshOptimizer::lockCellsInSubset(const word& subsetName)
{
    // lock the points in the cell subset with the given name
    label subsetI = mesh_.cellSubsetIndex(subsetName);
    if (subsetI >= 0)
    {
        DynamicList<label> lc;
        mesh_.cellsInSubset(subsetI, lc);

        lockCells(lc);
    }
    else
    {
        Warning << "Subset " << subsetName << " is not a cell subset!"
            << " Cannot lock cells!" << endl;
    }
}


void Foam::Module::meshOptimizer::removeUserConstraints()
{
    lockedFaces_.setSize(0);

    // unlock points
    # pragma omp parallel for schedule(dynamic, 50)
    forAll(vertexLocation_, i)
    {
        if (vertexLocation_[i] & LOCKED)
        {
            vertexLocation_[i] ^= LOCKED;
        }
    }
}


// ************************************************************************* //
