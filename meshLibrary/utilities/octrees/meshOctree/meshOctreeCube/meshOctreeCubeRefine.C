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

#include "meshOctreeCube.H"
#include "labelListListOps.H"
#include "triSurf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeCube::findContainedEdges
(
    const triSurf& surface,
    const boundBox& rootBox
)
{
    const labelListList2D& faceEdges = surface.facetEdges();
    const labelListList2D& edgeFaces = surface.edgeFacets();
    const DynamicList<edge>& edges = surface.edges();
    const pointField& points = surface.points();

    const labelListList2D& containedElements =
        activeSlotPtr_->containedTriangles_;
    labelListList2D& containedEdges = activeSlotPtr_->containedEdges_;

    DynamicList<label> addedEdges;
    labelHashSet addEdge;
    forAllRow(containedElements, containedElementsLabel_, tI)
    {
        const label facetI = containedElements[containedElementsLabel_][tI];

        forAllRow(faceEdges, facetI, feI)
        {
            const label edgeI = faceEdges[facetI][feI];

            if (addEdge.found(edgeI))
                continue;

            if (edgeFaces[edgeI].size() != 2)
                continue;

            if (
                surface[edgeFaces[edgeI][0]].region() !=
                surface[edgeFaces[edgeI][1]].region()
            )
            {
                const edge& edg = edges[edgeI];
                const point& s = points[edg.start()];
                const point& e = points[edg.end()];
                if (intersectsLine(rootBox, s, e))
                {
                    addEdge.insert(edgeI);
                    addedEdges.append(edgeI);
                }
            }
        }
    }

    if (addedEdges.size() != 0)
    {
        const label rowI = activeSlotPtr_->nextContainedEdgesRow_++;
        containedEdgesLabel_ = rowI;
        setRow(containedEdges, rowI, addedEdges);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeCube::refineCube2D
(
    const triSurf& surface,
    const boundBox& rootBox,
    meshOctreeSlot* slotPtr
)
{
    if (!slotPtr)
        slotPtr = activeSlotPtr_;

    // set the cube label to -1
    cubeLabel_ = -1;

    // create subCubes
    FixedList<meshOctreeCube*, 8> subCubes;
    forAll(subCubes, scI)
        subCubes[scI] = nullptr;

    // create new cubes in the Z-order fashion
    for (label scI = 0; scI < 4; ++scI)
    {
        const label cubeI = slotPtr->cubes_.size();
        const meshOctreeCubeCoordinates cc = this->refineForPosition(scI);

        slotPtr->cubes_.append(new meshOctreeCube(cc));

        subCubes[scI] = &slotPtr->cubes_[cubeI];
        subCubes[scI]->activeSlotPtr_ = slotPtr;
        subCubes[scI]->setCubeType(this->cubeType());
    }

    const label subCubesLabel = slotPtr->childCubes_.size();
    slotPtr->childCubes_.append(new FixedList<meshOctreeCube*, 8>(subCubes));
    subCubesPtr_ = &slotPtr->childCubes_[subCubesLabel][0];

    if (hasContainedElements())
    {
        const labelListList2D& containedElements =
            activeSlotPtr_->containedTriangles_;

        // check if the subCube contain the element
        FixedList<DynamicList<label>, 4> elementsInSubCubes;
        forAll(elementsInSubCubes, scI)
            elementsInSubCubes[scI].reserve(512);

        forAllRow(containedElements, containedElementsLabel_, tI)
        {
            const label elI =
                containedElements[containedElementsLabel_][tI];

            bool used(false);
            for (label scI = 0; scI < 4; ++scI)
                if
                (
                    subCubes[scI]->intersectsTriangleExact
                    (
                        surface,
                        rootBox,
                        elI
                    )
                )
                {
                    used = true;
                    elementsInSubCubes[scI].append(elI);
                }

            if (!used)
            {
                Warning << "Triangle " << elI
                    << " is not transferred to the child cubes!" << endl;
            }
        }

        forAll(elementsInSubCubes, scI)
        {
            const DynamicList<label>& elmts = elementsInSubCubes[scI];

            if (elmts.size() != 0)
            {
                const label rowI = slotPtr->nextContainedTrianglesRow_++;
                subCubes[scI]->containedElementsLabel_ = rowI;
                setRow(slotPtr->containedTriangles_, rowI, elmts);
            }
        }

        // find surface edges within the cube
        for (label scI = 0; scI < 4; ++scI)
            if (subCubes[scI]->hasContainedElements())
            {
                subCubes[scI]->findContainedEdges
                (
                    surface,
                    rootBox
                );
            }
            else if (subCubes[scI]->cubeType() & DATA)
            {
                subCubes[scI]->setCubeType(UNKNOWN);
            }
    }

}


void Foam::Module::meshOctreeCube::refineCube
(
    const triSurf& surface,
    const boundBox& rootBox,
    meshOctreeSlot* slotPtr
)
{
    if (!slotPtr)
        slotPtr = activeSlotPtr_;

    // set the cube label to -1
    cubeLabel_ = -1;

    // create subCubes
    FixedList<meshOctreeCube*, 8> subCubes;

    // create new cubes in the Z-order fashion
    for (label scI = 0; scI < 8; ++scI)
    {
        const label cubeI = slotPtr->cubes_.size();
        slotPtr->cubes_.append(new meshOctreeCube(this->refineForPosition(scI)));

        subCubes[scI] = &slotPtr->cubes_[cubeI];
        subCubes[scI]->activeSlotPtr_ = slotPtr;
        subCubes[scI]->setCubeType(this->cubeType());
    }

    const label subCubesLabel = slotPtr->childCubes_.size();
    slotPtr->childCubes_.append(new FixedList<meshOctreeCube*, 8>(subCubes));
    subCubesPtr_ = &slotPtr->childCubes_[subCubesLabel][0];

    if (hasContainedElements())
    {
        const labelListList2D& containedElements =
            activeSlotPtr_->containedTriangles_;

        // check if the subCube contain the element
        FixedList<DynamicList<label>, 8> elementsInSubCubes;
        forAll(elementsInSubCubes, scI)
            elementsInSubCubes[scI].reserve(512);

        forAllRow(containedElements, containedElementsLabel_, tI)
        {
            const label elI =
                containedElements[containedElementsLabel_][tI];

            bool used(false);
            for (label scI = 0; scI < 8; ++scI)
                if (
                    subCubes[scI]->intersectsTriangleExact
                    (
                        surface,
                        rootBox,
                        elI
                    )
                )
                {
                    used = true;
                    elementsInSubCubes[scI].append(elI);
                }

            if (!used)
            {
                Warning << "Triangle " << elI
                    << " is not transferred to the child cubes!" << endl;
            }
        }

        forAll(elementsInSubCubes, scI)
        {
            const DynamicList<label>& elmts = elementsInSubCubes[scI];

            if (elmts.size() != 0)
            {
                const label rowI = slotPtr->nextContainedTrianglesRow_++;
                subCubes[scI]->containedElementsLabel_ = rowI;
                setRow(slotPtr->containedTriangles_, rowI, elmts);
            }
        }

        // find surface edges within the cube
        for (label scI = 0; scI < 8; ++scI)
            if (subCubes[scI]->hasContainedElements())
            {
                subCubes[scI]->findContainedEdges
                (
                    surface,
                    rootBox
                );
            }
            else if (subCubes[scI]->cubeType() & DATA)
            {
                subCubes[scI]->setCubeType(UNKNOWN);
            }
    }

}


// ************************************************************************* //
