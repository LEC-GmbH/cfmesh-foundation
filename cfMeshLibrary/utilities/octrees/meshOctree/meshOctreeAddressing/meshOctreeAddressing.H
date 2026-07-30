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
    Foam::Module::meshOctreeAddressing

Description
    Automatic octree for mesh generation. This class checks if the current
    octree resolution is enough to resolve all surface features

SourceFiles
    meshOctreeAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeAddressing_H
#define meshOctreeAddressing_H

#include "FixedList.H"
#include "DynamicList.H"
#include "labelListListOps.H"
#include "meshOctreeCubeCoordinates.H"
#include "face.H"
#include "edge.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class IOdictionary;

namespace Module
{
class meshOctree;

/*---------------------------------------------------------------------------*\
                    Class meshOctreeAddressing Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeAddressing
{
    // Private data

        //- reference to the octree
        const meshOctree& octree_;

        //- checkForGluedMesh, extracted once at construction
        bool checkForGluedMesh_;

        //- nonManifoldMeshing, extracted once at construction
        bool nonManifoldMeshing_;

        //- removeCellsIntersectingPatches subDict, or dictionary::null
        //- if not present
        const dictionary& removeCellsIntersectingPatchesDict_;

        //- keepCellsIntersectingPatches subDict, or dictionary::null
        //- if not present
        const dictionary& keepCellsIntersectingPatchesDict_;

        //- use DATA boxes
        bool useDATABoxes_;

        //- number of created octree nodes
        mutable label nNodes_;

        //- coordinates of octree nodes
        mutable pointField* octreePointsPtr_;

        //- node labels
        mutable labelListList2D* nodeLabelsPtr_;

        //- node leaves
        mutable List<FixedList<label, 8>>* nodeLeavesPtr_;

        //- identify which boxes should be used as mesh cells
        mutable List<direction>* boxTypePtr_;

        //- identify created nodes as OUTERNODE or INNERNODE
        mutable List<direction>* nodeTypePtr_;


    // Additional addressing useful for some algorithms

        //- faces of the octree
        mutable labelListList2D* octreeFacesPtr_;
        mutable DynamicList<label>* octreeFacesOwnersPtr_;
        mutable DynamicList<label>* octreeFacesNeighboursPtr_;

        //- octree box-faces addressing
        mutable labelListList2D* leafFacesPtr_;

        //- node-faces addressing
        mutable labelListList2D* nodeFacesPtr_;

        //- leaf-leaves addressing
        mutable labelListList2D* leafLeavesPtr_;

        //- edges of the octree
        mutable DynamicList<edge>* octreeEdgesPtr_;

        //- edges-leaves addressing
        mutable labelListList2D* edgeLeavesPtr_;

        //- leaf-edges addressing
        mutable labelListList2D* leafEdgesPtr_;

        //- node-edges addressing
        mutable labelListList2D* nodeEdgesPtr_;

        //- face-edges addressing
        mutable labelListList2D* faceEdgesPtr_;

        //- edge-faces addressing
        mutable labelListList2D* edgeFacesPtr_;


    // Private member functions

        //- calculate octreePointsPtr_
        void createOctreePoints() const;

        //- calculate nodeLabelsPtr_
        void createNodeLabels() const;

        //- calculate nodeLeavesPtr_
        void createNodeLeaves() const;

        //- assemble boxTypePtr_ list
        void findUsedBoxes() const;

        //- calculate nodeTypePtr_
        void calculateNodeType() const;

        //- calculate faces
        void createOctreeFaces() const;

        //- calculate leaf-faces addressing
        void calculateLeafFaces() const;

        //- calculate node-faces addressing
        void calculateNodeFaces() const;

        //- calculate leaf-leaves addressing
        void calculateLeafLeaves() const;

        //- calculate edges
        void createOctreeEdges() const;

        //- calculate leaf-edges
        void calculateLeafEdges() const;

        //- calculate edge-leaves addressing
        void calculateEdgeLeaves() const;

        //- calculate edge-faces
        void calculateEdgeFaces() const;

        //- Clear allocated data
        void clearOut();
        void clearBoxTypes();
        void clearNodeAddressing();
        void clearOctreeFaces();
        void clearAddressing();

        //- check if distinct parts are glued together
        void checkGluedRegions();


    // Private copy constructor

        //- Disallow default bitwise copy construct
        meshOctreeAddressing(const meshOctreeAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeAddressing&);


public:

    //- Construct from surface and dictionary
    meshOctreeAddressing
    (
        const meshOctree& mo,
        const dictionary& dict,
        bool useDATABoxes = false
    );

    //- Destructor
    ~meshOctreeAddressing();

    // Enumerators

        enum octreeCubeTypes
        {
            NONE = 0,
            MESHCELL = 1,
            BOUNDARY = 2,
            SPLITHEX = 4,
            INNERNODE = 8,
            OUTERNODE = 16,
            BOUNDARYNODE = 32
        };


    // Member Functions

        //- check if there exist any non-signly connected edges and vertices
        //- in the octree mesh. This is primarily needed
        //- for the cartesian template
        void checkAndFixIrregularConnections();

        //- return number of octree nodes
        inline label numberOfNodes() const;

        //- return coordinates of octree vertices
        inline const pointField& octreePoints() const;

        //- return nodeLabels
        inline const labelListList2D& nodeLabels() const;

        //- return nodeLeaves
        inline const List<FixedList<label, 8>>& nodeLeaves() const;

        //- return which octree boxes are used for mesh creation
        inline const List<direction>& boxType() const;

        //- return type of node (INNERNODE,or OUTERNODE)
        inline const List<direction>& nodeType() const;

        //- set box type
        inline void setBoxType(const label boxI, const direction type);

        //- return octree faces, created for MESHCELL boxes
        inline const labelListList2D& octreeFaces() const;

        //- return owners of octree faces
        inline const DynamicList<label>& octreeFaceOwner() const;

        //- return neighbours of octree faces
        inline const DynamicList<label>& octreeFaceNeighbour() const;

        //- return octree box-faces addressing
        inline const labelListList2D& leafFaces() const;

        //- return leaf-leaves addressing
        inline const labelListList2D& leafLeaves() const;

        //- return node-faces addressing
        inline const labelListList2D& nodeFaces() const;

        //- return octree edges, created for MESHCELL boxes
        inline const DynamicList<edge>& octreeEdges() const;

        //- return edge-leaves addressing
        inline const labelListList2D& edgeLeaves() const;

        //- return leaf-edges addressing
        inline const labelListList2D& leafEdges() const;

        void edgeIntersections
        (
            const label eI,
            DynamicList<point>& intersections
        ) const;

        //- return node-edges addressing
        inline const labelListList2D& nodeEdges() const;

        //- return face-edges addressing
        inline const labelListList2D& faceEdges() const;

        //- return edge-faces addressing
        inline const labelListList2D& edgeFaces() const;

        //- return const reference to meshOctree
        inline const meshOctree& octree() const;

        //- find cubes around an edge (cubes must be at the same level)
        void cubesAroundEdge
        (
            const label leafI,
            const direction eI,
            FixedList<label, 4>& edgeCubes
        ) const;

        //- find edge centre if it exists
        label findEdgeCentre(const label leafI, const direction eI) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOctreeAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
