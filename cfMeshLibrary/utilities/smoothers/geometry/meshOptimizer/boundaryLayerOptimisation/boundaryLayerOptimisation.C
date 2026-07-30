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

#include "boundaryLayerOptimisation.H"
#include "meshSurfacePartitioner.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::boundaryLayerOptimisation::boundaryLayerOptimisation
(
    polyMeshGen& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    meshSurfacePtr_(new meshSurfaceEngine(mesh)),
    deleteMeshSurface_(true),
    partitionerPtr_(nullptr),
    maxNumIterations_(5),
    nSmoothNormals_(5),
    relThicknessTol_(0.1),
    featureSizeFactor_(0.3),
    reCalculateNormals_(true)
{
    calculateHairEdges();
    read(dict);
}


Foam::Module::boundaryLayerOptimisation::boundaryLayerOptimisation
(
    polyMeshGen& mesh,
    const meshSurfaceEngine& mse,
    const dictionary& dict
)
:
    mesh_(mesh),
    meshSurfacePtr_(&mse),
    deleteMeshSurface_(false),
    partitionerPtr_(nullptr),
    maxNumIterations_(5),
    nSmoothNormals_(5),
    relThicknessTol_(0.15),
    featureSizeFactor_(0.3),
    reCalculateNormals_(true)
{
    calculateHairEdges();
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::boundaryLayerOptimisation::~boundaryLayerOptimisation()
{
    deleteDemandDrivenData(partitionerPtr_);

    if (deleteMeshSurface_)
    {
        deleteDemandDrivenData(meshSurfacePtr_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::boundaryLayerOptimisation::setMaxNumIterations
(
    const label maxNumIterations
)
{
    maxNumIterations_ = maxNumIterations;
}


void Foam::Module::boundaryLayerOptimisation::setNumNormalsSmoothingIterations
(
    const label nSmoothNormals
)
{
    nSmoothNormals_ = nSmoothNormals;
}


void Foam::Module::boundaryLayerOptimisation::recalculateNormals
(
    const bool shallRecalculate
)
{
    reCalculateNormals_ = shallRecalculate;
}


void Foam::Module::boundaryLayerOptimisation::setRelativeThicknessTolerance
(
    const scalar relThicknessTol
)
{
    relThicknessTol_ = relThicknessTol;
}


void Foam::Module::boundaryLayerOptimisation::setFeatureSizeFactor
(
    const scalar featureSizeFactor
)
{
    featureSizeFactor_ = featureSizeFactor;
}


const Foam::DynamicList<Foam::edge>&
Foam::Module::boundaryLayerOptimisation::hairEdges() const
{
    return hairEdges_;
}


const Foam::Module::labelListList2D&
Foam::Module::boundaryLayerOptimisation::hairEdgesAtBndPoint() const
{
    return hairEdgesAtBndPoint_;
}


const Foam::boolList&
Foam::Module::boundaryLayerOptimisation::isBaseFace() const
{
    return isBndLayerBase_;
}


void Foam::Module::boundaryLayerOptimisation::read(const dictionary& layersDict)
{
    if (!layersDict.lookupOrDefault<bool>("optimiseLayer", true))
    {
        return;
    }

    const dictionary& optParams =
        layersDict.subOrEmptyDict("optimisationParameters");

    recalculateNormals
    (
        optParams.lookupOrDefault("recalculateNormals", reCalculateNormals_)
    );

    const label nSmoothNormals =
        optParams.lookupOrDefault("nSmoothNormals", nSmoothNormals_);

    if (nSmoothNormals < 0)
    {
        FatalErrorInFunction
            << "nSmoothNormals must not be negative!"
            << exit(FatalError);
    }

    setNumNormalsSmoothingIterations(nSmoothNormals);

    const scalar featureSizeFactor =
        optParams.lookupOrDefault("featureSizeFactor", featureSizeFactor_);

    if (featureSizeFactor >= 1.0 || featureSizeFactor < 0.0)
    {
        FatalErrorInFunction
            << "Feature size factor is out"
            << " of a valid range 0 to 1" << exit(FatalError);
    }

    setFeatureSizeFactor(featureSizeFactor);

    const scalar relThicknessTol =
        optParams.lookupOrDefault("relThicknessTol", relThicknessTol_);

    if (relThicknessTol >= 1.0 || relThicknessTol < 0.0)
    {
        FatalErrorInFunction
            << "Relative thickness tolerance is out"
            << " of a valid range 0 to 1" << exit(FatalError);
    }

    setRelativeThicknessTolerance(relThicknessTol);

    const label maxNumIterations =
        optParams.lookupOrDefault("maxNumIterations", maxNumIterations_);

    if (maxNumIterations < 0)
    {
        FatalErrorInFunction
            << "maxNumIterations must not be negative!"
            << exit(FatalError);
    }

    setMaxNumIterations(maxNumIterations);
}


// ************************************************************************* //
