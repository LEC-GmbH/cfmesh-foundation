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
    Foam::Module::meshOctreeCreator

Description
    Creates octree for mesh generation

SourceFiles
    meshOctreeCreator.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeCreator_H
#define meshOctreeCreator_H

#include "boolList.H"
#include "dictionary.H"
#include "meshOctreeModifier.H"
#include "patchRefinementList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

// Forward declarations
class meshOctreeCube;

/*---------------------------------------------------------------------------*\
                      Class meshOctreeCreator Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeCreator
{
protected:

    // Protected data

        //- Reference to meshOctree
        meshOctree& octree_;

        //- Scaling factor
        scalar scalingFactor_;

        //- Full meshDict, kept only to pass on to meshOctreeAutomaticRefinement
        //- (which needs keepCellsIntersectingBoundary alongside
        //- maxCellSize/minCellSize, so no single subDict covers it) -
        //- every other key this class itself reads is extracted once, in
        //- the constructor, into the narrowly-scoped members below instead
        //- of being re-read ad hoc from here. dictionary::null if no dict
        //- was given.
        const dictionary& meshDict_;

        //- was a real dict given to the constructor (as opposed to the
        //- default dictionary::null)? Needed only for the top-level scalar
        //- reads (maxCellSize) below, which would otherwise FatalError via
        //- .lookup() on a genuinely-absent dict.
        bool haveMeshDict_;

        //- maxCellSize, validated and extracted once at construction
        //- (only meaningful if haveMeshDict_ is true)
        scalar maxCellSize_;

        //- boundaryCellSize/boundaryCellSizeRefinementThickness, if present
        bool haveBoundaryCellSize_;
        scalar boundaryCellSize_;
        scalar boundaryCellSizeRefinementThickness_;

        //- patchCellSize subDict, or dictionary::null if not present
        const dictionary& patchCellSizeDict_;

        //- subsetCellSize subDict, or dictionary::null if not present
        const dictionary& subsetCellSizeDict_;

        //- localRefinement subDict, or dictionary::null if not present
        const dictionary& localRefinementDict_;

        //- objectRefinements subDict, or dictionary::null if not present
        const dictionary& objectRefinementsDict_;

        //- surfaceMeshRefinement subDict, or dictionary::null if not present
        const dictionary& surfaceMeshRefinementDict_;

        //- edgeMeshRefinement subDict, or dictionary::null if not present
        const dictionary& edgeMeshRefinementDict_;

        //- anisotropicSources subDict, or dictionary::null if not present
        const dictionary& anisotropicSourcesDict_;

        //- hex refinement flag
        bool hexRefinement_;


private:

    // Private member functions

        // Octree refinement

        //- refine boxes contained inside the objects for refinement
        void refineBoxesContainedInObjects();

        //- refine boxes intersected by surface meshes
        //- used as refinement sources
        void refineBoxesIntersectingSurfaces();

        //- refine boxes by edge meshes
        //- used as refinement sources
        void refineBoxesIntersectingEdgeMeshes();

        //- refine boxes near DATA boxes to get a nice smooth surface
        void refineBoxesNearDataBoxes(const label nLayers = 1);

        //- refine boxes of the given flag to the given size
        void refineBoxes(const direction refLevel, const direction cubeType);

        // Private octree creation methods
        void refineBoundary();
        void createInsideOutsideInformation();
        void refineInsideAndUnknownBoxes();


        // information about octree refinement
        //- ref level to achieve max cell size
        direction globalRefLevel_;

        //- this list contains ref level for each surface triangle
        List<DynamicList<Tuple2<direction, scalar>>> surfRefLevel_;

        //- set the boundBox such that maxCellSize is achieved
        void setRootCubeSizeAndRefParameters();


    // Private copy constructor

        //- Disallow default bitwise copy construct
        meshOctreeCreator(const meshOctreeCreator&);

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeCreator&);


public:

    // Constructors

        //- Construct from meshOctree and, optionally, meshDict (or its
        //- subDict, e.g. via meshDict.subOrEmptyDict("someKey")).
        //- dictionary::null (the default) means no dict was given.
        meshOctreeCreator
        (
            meshOctree& mo,
            const dictionary& dict = dictionary::null
        );


    //- Destructor
    ~meshOctreeCreator() = default;


    // Member Functions

        //- activate octree refinement which marks all sons of an octree box
        //- for refinement in case a single son is marked for refinement
        //- this type of refinement is necessary for creating hex meshes
        void activateHexRefinement();

        //- create octree boxes
        void createOctreeBoxes();

        //- refine boxes containing surface elements
        //- this is used for proximity searches
        void createOctreeWithRefinedBoundary
        (
            const direction maxLevel,
            const label nTrianglesInLeaf = 15
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
