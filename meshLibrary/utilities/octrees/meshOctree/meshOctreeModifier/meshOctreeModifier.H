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
    Foam::Module::meshOctreeModifier

Description
    Octree for mesh generation

SourceFiles
    meshOctreeModifier.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeModifier_H
#define meshOctreeModifier_H

#include "meshOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

/*---------------------------------------------------------------------------*\
                     Class meshOctreeModifier Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeModifier
{
    // Private data

        //- Reference to the octree
        meshOctree& octree_;


    // Private member functions

        //- correct refinement such that it produces 1-irregular octree
        void ensureCorrectRegularity(labelList& refineBox);

        //- add additional refinement with the following rule. If a son
        //- of a father box is refined, then refine all other sons
        bool ensureCorrectRegularitySons(labelList& refineBox);


    // Private copy constructor

        //- Disallow default bitwise copy construct
        meshOctreeModifier(const meshOctreeModifier&);

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeModifier&);


public:

    //- Construct from octree
    meshOctreeModifier(meshOctree&);

    //- Destructor
    ~meshOctreeModifier() = default;


    // Member Functions

        //- return octree
        inline const meshOctree& octree() const;

        //- return rootBox
        inline boundBox& rootBoxAccess();

        //- return isRootInitialised_
        inline bool& isRootInitialisedAccess();

        //- return search range
        inline scalar& searchRangeAccess();

        //- return initial cube
        inline meshOctreeCube& initialCubeAccess();

        //- return octree slots
        inline List<meshOctreeSlot>& dataSlotsAccess();

        //- return the pointer to the meshOctreeCube at the given position
        inline meshOctreeCube* findCubeForPosition
        (
            const meshOctreeCubeCoordinates&
        ) const;

        //- find leaves contained in a given boundBox
        inline void findLeavesContainedInBox
        (
            const boundBox&,
            DynList<const meshOctreeCube*, 256>&
        ) const;

        //- return leaves
        inline LongList<meshOctreeCube*>& leavesAccess();

        //- create leaves
        inline void createListOfLeaves();

        //- mark additional layers around the leaves selected for refinement
        void markAdditionalLayers
        (
            labelList& refineBox,
            const label nLayers = 1
        ) const;

        //- mark additional layers around the leaves selected for refinement
        void markAdditionalLayersOfFaceNeighbours
        (
            labelList& refineBox,
            const label nLayers = 1
        ) const;

        //- mark additional layers around the leaves selected for refinement
        //- given on a box-by-box basis
        //- returns the number of boxes selected for refinement
        label markAdditionalLayers
        (
            labelList& refineBox,
            labelList& nLayers,
            List<direction>& targetRefLevel
        ) const;

        //- refine leaves marked for refinement
        //- hexRefinement is activated when it is required to refine all
        //- sons of the same father, if a single son gets marked for refinement
        void refineSelectedBoxes
        (
            labelList& refineBox,
            const bool hexRefinement = false
        );

        //- refine selected boxes and the boxes within the given range
        void refineSelectedBoxesAndAdditionalLayers
        (
            labelList& refineBox,
            const scalarList& refThickness
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOctreeModifierI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
