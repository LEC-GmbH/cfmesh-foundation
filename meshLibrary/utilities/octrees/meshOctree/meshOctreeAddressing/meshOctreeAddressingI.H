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

#include "meshOctreeAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::label Foam::Module::meshOctreeAddressing::numberOfNodes() const
{
    if (!nodeLabelsPtr_)
    {
        createNodeLabels();
    }

    return nNodes_;
}


inline const Foam::pointField&
Foam::Module::meshOctreeAddressing::octreePoints() const
{
    if (!octreePointsPtr_)
    {
        createOctreePoints();
    }

    return *octreePointsPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::nodeLabels() const
{
    if (!nodeLabelsPtr_)
    {
        createNodeLabels();
    }

    return *nodeLabelsPtr_;
}


inline const Foam::List<Foam::FixedList<Foam::label, 8>>&
Foam::Module::meshOctreeAddressing::nodeLeaves() const
{
    if (!nodeLeavesPtr_)
    {
        createNodeLeaves();
    }

    return *nodeLeavesPtr_;
}


inline const Foam::List<Foam::direction>&
Foam::Module::meshOctreeAddressing::boxType() const
{
    if (!boxTypePtr_)
    {
        findUsedBoxes();
    }

    return *boxTypePtr_;
}


inline void Foam::Module::meshOctreeAddressing::setBoxType
(
    const label boxI,
    const direction type
)
{
    if (boxTypePtr_)
    {
        List<direction>& boxType = *boxTypePtr_;
        boxType[boxI] |= type;
    }
}


inline const Foam::Module::meshOctree&
Foam::Module::meshOctreeAddressing::octree() const
{
    return octree_;
}


inline const Foam::List<Foam::direction>&
Foam::Module::meshOctreeAddressing::nodeType() const
{
    if (!nodeTypePtr_)
    {
        calculateNodeType();
    }

    return *nodeTypePtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::octreeFaces() const
{
    if (!octreeFacesPtr_)
    {
        createOctreeFaces();
    }

    return *octreeFacesPtr_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::Module::meshOctreeAddressing::octreeFaceOwner() const
{
    if (!octreeFacesOwnersPtr_)
    {
        createOctreeFaces();
    }

    return *octreeFacesOwnersPtr_;
}


inline const Foam::DynamicList<Foam::label>&
Foam::Module::meshOctreeAddressing::octreeFaceNeighbour() const
{
    if (!octreeFacesNeighboursPtr_)
    {
        createOctreeFaces();
    }

    return *octreeFacesNeighboursPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::leafFaces() const
{
    if (!leafFacesPtr_)
    {
        calculateLeafFaces();
    }

    return *leafFacesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::leafLeaves() const
{
    if (!leafLeavesPtr_)
    {
        calculateLeafLeaves();
    }

    return *leafLeavesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::nodeFaces() const
{
    if (!nodeFacesPtr_)
    {
        calculateNodeFaces();
    }

    return *nodeFacesPtr_;
}


inline const Foam::DynamicList<Foam::edge>&
Foam::Module::meshOctreeAddressing::octreeEdges() const
{
    if (!octreeEdgesPtr_)
    {
        createOctreeEdges();
    }

    return *octreeEdgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::edgeLeaves() const
{
    if (!edgeLeavesPtr_)
    {
        calculateEdgeLeaves();
    }

    return *edgeLeavesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::leafEdges() const
{
    if (!leafEdgesPtr_)
    {
        calculateLeafEdges();
    }

    return *leafEdgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::nodeEdges() const
{
    if (!nodeEdgesPtr_)
    {
        createOctreeEdges();
    }

    return *nodeEdgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::faceEdges() const
{
    if (!faceEdgesPtr_)
    {
        createOctreeEdges();
    }

    return *faceEdgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshOctreeAddressing::edgeFaces() const
{
    if (!edgeFacesPtr_)
    {
        calculateEdgeFaces();
    }

    return *edgeFacesPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
