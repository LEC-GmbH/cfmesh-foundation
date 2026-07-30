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

#include "meshOctreeAutomaticRefinement.H"
#include "demandDrivenData.H"
#include "triSurfacePartitioner.H"
#include "triSurfaceCurvatureEstimator.H"
#include "IOdictionary.H"
#include "triSurf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeAutomaticRefinement::createSurfacePartitioner()
const
{
    partitionerPtr_ = new triSurfacePartitioner(octree_.surface());
}


const Foam::Module::triSurfacePartitioner&
Foam::Module::meshOctreeAutomaticRefinement::partitioner() const
{
    if (!partitionerPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate addressing!" << abort(FatalError);
        }

        createSurfacePartitioner();
    }

    return *partitionerPtr_;
}


void Foam::Module::meshOctreeAutomaticRefinement::createCurvatureEstimator()
const
{
    curvaturePtr_ = new triSurfaceCurvatureEstimator(octree_.surface());
}


const Foam::Module::triSurfaceCurvatureEstimator&
Foam::Module::meshOctreeAutomaticRefinement::curvature() const
{
    if (!curvaturePtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate addressing!" << abort(FatalError);
        }

        createCurvatureEstimator();
    }

    return *curvaturePtr_;
}


void Foam::Module::meshOctreeAutomaticRefinement::setMaxRefLevel
(
    const dictionary& dict
)
{
    const boundBox& rootBox = octree_.rootBox();
    const scalar size = rootBox.max().x() - rootBox.min().x();
    maxRefLevel_ = 0;

    if (dict.found("minCellSize"))
    {
        const scalar maxSize(readScalar(dict.lookup("maxCellSize")));
        scalar cs(readScalar(dict.lookup("minCellSize")));

        if (cs < 0)
        {
            FatalErrorInFunction
                << "Minimum cell size for automatic refinement is"
                << " negative!!" << exit(FatalError);
        }

        cs *= (1.0 + SMALL);

        if (cs > maxSize)
        {
            return;
        }

        bool finished;
        do
        {
            finished = false;

            const scalar lSize = size/pow(label(2), label(maxRefLevel_));
            // Or: const scalar lSize = size/(1L << maxRefLevel_);

            if (lSize < cs)
            {
                finished = true;
            }
            else
            {
                ++maxRefLevel_;
            }
        } while (!finished);

        useDATABoxes_ = true;

        Info<< "Requested min cell size corresponds to octree level "
            << label(maxRefLevel_) << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshOctreeAutomaticRefinement::meshOctreeAutomaticRefinement
(
    meshOctree& mo,
    const dictionary& dict,
    bool useDATABoxes
)
:
    octree_(mo),
    useDATABoxes_(useDATABoxes),
    hexRefinement_(false),
    partitionerPtr_(nullptr),
    curvaturePtr_(nullptr),
    maxRefLevel_(0)
{
    if (!useDATABoxes_)
    {
        dict.readIfPresent("keepCellsIntersectingBoundary", useDATABoxes_);
    }

    // calculate maximum allowed refinement level from the minimum cell size
    setMaxRefLevel(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::meshOctreeAutomaticRefinement::~meshOctreeAutomaticRefinement()
{
    deleteDemandDrivenData(curvaturePtr_);
    deleteDemandDrivenData(partitionerPtr_);
}


// ************************************************************************* //
