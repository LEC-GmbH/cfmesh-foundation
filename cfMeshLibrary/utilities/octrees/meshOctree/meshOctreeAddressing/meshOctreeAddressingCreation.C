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
#include "helperFunctions.H"
#include "labelListListOps.H"
#include "demandDrivenData.H"
#include "meshOctree.H"
#include "DynamicList.H"
#include "triSurf.H"

#include <omp.h>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeAddressing::createOctreePoints() const
{
    const labelListList2D& nodeLabels = this->nodeLabels();
    const boundBox& rootBox = octree_.rootBox();

    octreePointsPtr_ = new pointField(nNodes_);
    pointField& octreePoints = *octreePointsPtr_;

    const label nLeaves = nodeLabels.size();
    # pragma omp parallel for schedule(guided, 100)
    for (label cubeI = 0; cubeI < nLeaves; ++cubeI)
    {
        if (nodeLabels[cubeI].size() == 0)
        {
            continue;
        }

        FixedList<point, 8> vertices;
        const meshOctreeCubeBasic& oc = octree_.returnLeaf(cubeI);
        oc.vertices(rootBox, vertices);

        forAllRow(nodeLabels, cubeI, nI)
        {
            const label nodeI = nodeLabels[cubeI][nI];

            octreePoints[nodeI] = vertices[nI];
        }
    }
}


void Foam::Module::meshOctreeAddressing::createNodeLabels() const
{
    const List<direction>& boxType = this->boxType();

    nodeLabelsPtr_ = new labelListList2D(octree_.numberOfLeaves());
    labelListList2D& nodeLabels = *nodeLabelsPtr_;

    // allocate storage for node labels
    forAll(nodeLabels, leafI)
    {
        if (boxType[leafI])
        {
            nodeLabels[leafI].setSize(8);

            forAllRow(nodeLabels, leafI, i)
            {
                nodeLabels[leafI][i] = -1;
            }
        }
    }

    // start creating node labels
    nNodes_ = 0;
    DynamicList<label> numLocalNodes;
    # pragma omp parallel
    {
        const label nThreads = omp_get_num_threads();
        const label threadI = omp_get_thread_num();

        # pragma omp master
        numLocalNodes.setSize(nThreads);

        # pragma omp barrier

        // count the number of nodes local to each process
        label& nLocalNodes = numLocalNodes[threadI];
        nLocalNodes = 0;

        # pragma omp for schedule(static, 100)
        forAll(nodeLabels, leafI)
        {
            forAllRow(nodeLabels, leafI, nI)
            {
                if (nodeLabels[leafI][nI] != -1)
                {
                    continue;
                }

                FixedList<label, 8> pLeaves;
                octree_.findLeavesForCubeVertex(leafI, nI, pLeaves);

                FixedList<bool, 8> validLeaf(true);
                label minLeaf(leafI);
                forAll(pLeaves, plI)
                {
                    if (pLeaves[plI] > -1)
                    {
                        for (label i = plI + 1; i < 8; ++i)
                        {
                            if (pLeaves[plI] == pLeaves[i])
                            {
                                validLeaf[plI] = false;
                                validLeaf[i] = false;
                            }
                        }

                        if (!boxType[pLeaves[plI]])
                        {
                            validLeaf[plI] = false;
                            pLeaves[plI] = -1;
                        }

                        if (validLeaf[plI])
                        {
                            minLeaf = Foam::min(minLeaf, pLeaves[plI]);
                        }
                    }
                    else
                    {
                        validLeaf[plI] = false;
                    }
                }

                if ((minLeaf == leafI) && validLeaf[7 - nI])
                {
                    forAll(pLeaves, plI)
                    {
                        if (validLeaf[plI])
                        {
                            // set node labels to -2 not to repeat searches
                            nodeLabels[pLeaves[plI]][7 - plI] = -2;
                        }
                    }

                    ++nLocalNodes;
                }
            }
        }

        // set start node for each process
        # pragma omp barrier

        label startNode(0);
        for (label i = 0; i < threadI; ++i)
        {
            startNode += numLocalNodes[i];
        }

        // start creating node labels
        # pragma omp for schedule(static, 100)
        forAll(nodeLabels, leafI)
        {
            forAllRow(nodeLabels, leafI, nI)
            {
                if (nodeLabels[leafI][nI] >= 0)
                {
                    continue;
                }

                FixedList<label, 8> pLeaves;
                octree_.findLeavesForCubeVertex(leafI, nI, pLeaves);

                FixedList<bool, 8> validLeaf(true);
                label minLeaf(leafI);
                forAll(pLeaves, plI)
                {
                    if (pLeaves[plI] > -1)
                    {
                        for (label i = plI + 1; i < 8; ++i)
                        {
                            if (pLeaves[plI] == pLeaves[i])
                            {
                                validLeaf[plI] = false;
                                validLeaf[i] = false;
                            }
                        }

                        if (!boxType[pLeaves[plI]])
                        {
                            validLeaf[plI] = false;
                            pLeaves[plI] = -1;
                        }

                        if (validLeaf[plI])
                        {
                            minLeaf = Foam::min(minLeaf, pLeaves[plI]);
                        }
                    }
                    else
                    {
                        validLeaf[plI] = false;
                    }
                }

                if ((minLeaf == leafI) && validLeaf[7 - nI])
                {
                    forAll(pLeaves, plI)
                    {
                        if (validLeaf[plI])
                        {
                            // store the vertex at the corresponding
                            // location in the cube
                            nodeLabels[pLeaves[plI]][7 - plI] = startNode;
                        }
                    }

                    // store vertex label
                    ++startNode;
                }
            }
        }

        // set the number of nodes
        # pragma omp critical
        {
            nNodes_ = Foam::max(nNodes_, startNode);
        }
    }
}


void Foam::Module::meshOctreeAddressing::createNodeLeaves() const
{
    const List<direction>& boxType = this->boxType();
    const labelListList2D& nodeLabels = this->nodeLabels();

    // allocate nodeLeavesPtr_
    nodeLeavesPtr_ = new List<FixedList<label, 8>>(nNodes_);
    List<FixedList<label, 8>>& nodeLeaves = *nodeLeavesPtr_;

    boolList storedNode(nNodes_, false);
    # pragma omp parallel for schedule(dynamic, 100)
    forAll(nodeLabels, leafI)
    {
        forAllRow(nodeLabels, leafI, nI)
        {
            const label nodeI = nodeLabels[leafI][nI];

            if (storedNode[nodeI])
            {
                continue;
            }

            storedNode[nodeI] = true;

            FixedList<label, 8> pLeaves;
            octree_.findLeavesForCubeVertex(leafI, nI, pLeaves);

            forAll(pLeaves, plI)
            {
                if (pLeaves[plI] < 0)
                {
                    continue;
                }

                if (!boxType[pLeaves[plI]])
                {
                    pLeaves[plI] = -1;
                }
            }

            nodeLeaves[nodeI] = pLeaves;
        }
    }
}


void Foam::Module::meshOctreeAddressing::findUsedBoxes() const
{
    boxTypePtr_ = new List<direction>(octree_.numberOfLeaves(), NONE);
    List<direction>& boxType = *boxTypePtr_;

    # pragma omp parallel for schedule(dynamic, 40)
    forAll(boxType, leafI)
    {
        const meshOctreeCubeBasic& leaf = octree_.returnLeaf(leafI);

        if
        (
            !octree_.hasContainedTriangles(leafI)
         && !octree_.hasContainedEdges(leafI)
         && (leaf.cubeType() & meshOctreeCubeBasic::INSIDE)
        )
            boxType[leafI] |= MESHCELL;
    }

    if (nonManifoldMeshing_)
    {
        # pragma omp parallel for schedule(dynamic, 40)
        forAll(boxType, leafI)
        {
            const meshOctreeCubeBasic& leaf = octree_.returnLeaf(leafI);
            if (leaf.cubeType() & meshOctreeCubeBasic::UNKNOWN)
            {
                boxType[leafI] |= MESHCELL;
            }
        }
    }

    if (useDATABoxes_)
    {
        Info<< "Using DATA boxes" << endl;

        forAll(boxType, leafI)
        {
            if
            (
                octree_.hasContainedTriangles(leafI)
             || octree_.hasContainedEdges(leafI)
            )
            {
                boxType[leafI] |= MESHCELL;
            }
        }

        // do not use boxes intersecting given patches
        if (!removeCellsIntersectingPatchesDict_.empty())
        {
            wordHashSet patchesToRemove;

            const wordList patchNames =
                removeCellsIntersectingPatchesDict_.toc();
            forAll(patchNames, patchI)
            {
                patchesToRemove.insert(patchNames[patchI]);
            }

            const triSurf& ts = octree_.surface();
            boolList removeFacets(ts.size(), false);

            // remove facets in patches
            for (const word& patchName : patchesToRemove.toc())
            {
                const labelList matchedPatches = ts.findPatches(patchName);
                boolList activePatch(ts.patches().size(), false);
                forAll(matchedPatches, ptchI)
                {
                    activePatch[matchedPatches[ptchI]] = true;
                }

                forAll(ts, triI)
                {
                    if (activePatch[ts[triI].region()])
                    {
                        removeFacets[triI] = true;
                    }
                }
            }

            // remove facets in subsets
            for (const word& patchName : patchesToRemove.toc())
            {
                const label subsetID = ts.facetSubsetIndex(patchName);
                if (subsetID >= 0)
                {
                    DynamicList<label> facets;
                    ts.facetsInSubset(subsetID, facets);

                    forAll(facets, i)
                    {
                        removeFacets[facets[i]] = true;
                    }
                }
            }

            // set BOUNDARY flag to boxes intersected by the given facets
            DynamicList<label> containedTriangles;
            forAll(boxType, leafI)
            {
                octree_.containedTriangles(leafI, containedTriangles);

                forAll(containedTriangles, i)
                {
                    if (removeFacets[containedTriangles[i]])
                    {
                        boxType[leafI] = NONE;
                    }
                }
            }
        }
    }
    else if (!keepCellsIntersectingPatchesDict_.empty())
    {
        wordHashSet patchesToKeep;

        const wordList patchNames = keepCellsIntersectingPatchesDict_.toc();

        forAll(patchNames, patchI)
        {
            patchesToKeep.insert(patchNames[patchI]);
        }

        const triSurf& ts = octree_.surface();
        boolList keepFacets(ts.size(), false);

        // keep facets in patches
        for (const word& patchName : patchesToKeep.toc())
        {
            const labelList matchedPatches = ts.findPatches(patchName);
            boolList activePatch(ts.patches().size(), false);
            forAll(matchedPatches, ptchI)
            {
                activePatch[matchedPatches[ptchI]] = true;
            }

            forAll(ts, triI)
            {
                if (activePatch[ts[triI].region()])
                {
                    keepFacets[triI] = true;
                }
            }
        }

        // keep facets in subsets
        for (const word& patchName : patchesToKeep.toc())
        {
            const label subsetID = ts.facetSubsetIndex(patchName);

            if (subsetID >= 0)
            {
                DynamicList<label> facets;
                ts.facetsInSubset(subsetID, facets);

                forAll(facets, i)
                {
                    keepFacets[facets[i]] = true;
                }
            }
        }

        // set MESHCELL flag to boxes intersected by the given facets
        DynamicList<label> containedTriangles;
        forAll(boxType, leafI)
        {
            octree_.containedTriangles(leafI, containedTriangles);

            forAll(containedTriangles, i)
            {
                if (keepFacets[containedTriangles[i]])
                {
                    boxType[leafI] = MESHCELL;
                }
            }
        }
    }

    // set BOUNDARY flag to boxes which do not have a MESHCELL flag
    DynamicList<label> neighs;
    # pragma omp parallel for if (boxType.size() > 1000) \
    private(neighs) schedule(dynamic, 20)
    forAll(boxType, leafI)
    {
        if (boxType[leafI] & MESHCELL)
        {
            for (label i = 0; i < 6; ++i)
            {
                neighs.clear();
                octree_.findNeighboursInDirection(leafI, i, neighs);

                forAll(neighs, neiI)
                {
                    const label neiLabel = neighs[neiI];

                    if (neiLabel < 0)
                    {
                        continue;
                    }

                    if (!(boxType[neiLabel] & MESHCELL))
                    {
                        boxType[neiLabel] = BOUNDARY;
                    }
                }
            }
        }
    }

}


void Foam::Module::meshOctreeAddressing::calculateNodeType() const
{
    const List<FixedList<label, 8>>& nodeLeaves = this->nodeLeaves();

    nodeTypePtr_ = new List<direction>(nNodes_, NONE);
    List<direction>& nodeType = *nodeTypePtr_;

    # pragma omp parallel for schedule(static, 1)
    forAll(nodeLeaves, nodeI)
    {
        for (label nlI = 0; nlI < 8; ++nlI)
        {
            const label leafI = nodeLeaves[nodeI][nlI];

            if (leafI == -1)
            {
                continue;
            }

            const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

            if
            (
                (oc.cubeType() & meshOctreeCubeBasic::OUTSIDE)
             || (oc.cubeType() & meshOctreeCubeBasic::UNKNOWN)
            )
            {
                nodeType[nodeI] |= OUTERNODE;
                break;
            }
            else if
            (
                !octree_.hasContainedTriangles(leafI)
             && (oc.cubeType() & meshOctreeCubeBasic::INSIDE)
            )
            {
                nodeType[nodeI] |= INNERNODE;
                break;
            }
        }
    }
}


void Foam::Module::meshOctreeAddressing::createOctreeFaces() const
{
    octreeFacesPtr_ = new labelListList2D();
    octreeFacesOwnersPtr_ = new DynamicList<label>();
    octreeFacesNeighboursPtr_ = new DynamicList<label>();

    const labelListList2D& nodeLabels = this->nodeLabels();
    const List<direction>& boxType = this->boxType();
    this->nodeLeaves();

    label nFaces(0);
    labelList rowSizes, chunkSizes;

    # pragma omp parallel
    {
        // faces are created and stored into helper arrays, and each thread
        // allocates its own graph for storing faces. The faces are generated
        // by dividing the octree leaves into chunks, and distributing these
        // chunks over the threads. There are four chunks per each thread to
        // improve load balancing. The number of faces generated in each chunk
        // is stored and later in used to store the faces into the octree faces
        // graph in the correct order
        DynamicList<DynamicList<label>> helperFaces;
        DynamicList<label> helperOwner, helperNeighbour;

        const label nThreads = omp_get_num_threads();
        const label threadI = omp_get_thread_num();
        const label nChunks = 4*omp_get_num_threads();
        const label chunkSize = boxType.size() / nChunks + 1;

        # pragma omp master
        {
            chunkSizes.setSize(nChunks);
            chunkSizes = 0;
        }

        # pragma omp barrier

        for
        (
            label chunkI = threadI;
            chunkI < nChunks;
            chunkI += nThreads
        )
        {
            const label start = chunkSize*chunkI;
            const label end = Foam::min(start + chunkSize, boxType.size());

            const label nBefore = helperFaces.size();

            for (label leafI = start; leafI < end; ++leafI)
            {
                const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

                if (boxType[leafI] & MESHCELL)
                {
                    FixedList<label, 12> edgeCentreLabel(-1);
                    for (label i = 0; i < 12; ++i)
                    {
                        edgeCentreLabel[i] = findEdgeCentre(leafI, i);
                    }

                    for (label fI = 0; fI < 6; ++fI)
                    {
                        DynamicList<label> neighbours;
                        octree_.findNeighboursInDirection
                        (
                            leafI,
                            fI,
                            neighbours
                        );

                        if (neighbours.size() != 1)
                        {
                            continue;
                        }

                        const label nei = neighbours[0];

                        // create face
                        DynamicList<label> f;
                        f.reserve(8);
                        for (label pI = 0; pI < 4; ++pI)
                        {
                            const label nI =
                                meshOctreeCubeCoordinates::faceNodes_[fI][pI];
                            const label feI =
                                meshOctreeCubeCoordinates::faceEdges_[fI][pI];

                            f.append(nodeLabels[leafI][nI]);

                            if (edgeCentreLabel[feI] != -1)
                            {
                                f.append(edgeCentreLabel[feI]);
                            }
                        }

                        if (nei < 0)
                        {
                            // face is at the boundary of the octree
                            appendList(helperFaces, f);
                            helperOwner.append(leafI);
                            helperNeighbour.append(-1);
                        }
                        else if (boxType[nei] & MESHCELL)
                        {
                            // face is an internal face
                            if (nei > leafI)
                            {
                                appendList(helperFaces, f);
                                helperOwner.append(leafI);
                                helperNeighbour.append(nei);
                            }
                            else if
                            (
                                octree_.returnLeaf(nei).level() < oc.level()
                            )
                            {
                                // append a reversed face
                                label i(1);
                                for (label j = f.size()-1; j > i; --j)
                                {
                                    const label add = f[j];
                                    f[j] = f[i];
                                    f[i] = add;
                                    ++i;
                                }

                                appendList(helperFaces, f);
                                helperOwner.append(nei);
                                helperNeighbour.append(leafI);
                            }
                        }
                        else if (boxType[nei] & BOUNDARY)
                        {
                            // face is at the boundary of the mesh cells
                            appendList(helperFaces, f);
                            helperOwner.append(leafI);
                            helperNeighbour.append(nei);
                        }
                    }
                }
                else if (boxType[leafI] & BOUNDARY)
                {
                    for (label fI = 0; fI < 6; ++fI)
                    {
                        DynamicList<label> neighbours;
                        octree_.findNeighboursInDirection
                        (
                            leafI,
                            fI,
                            neighbours
                        );

                        if (neighbours.size() != 1)
                        {
                            continue;
                        }
                        const label nei = neighbours[0];
                        if (nei < 0)
                        {
                            continue;
                        }
                        if
                        (
                            (boxType[nei] & MESHCELL)
                         && (octree_.returnLeaf(nei).level() < oc.level())
                        )
                        {
                            // add a boundary face
                            const label* fNodes =
                                meshOctreeCubeCoordinates::faceNodes_[fI];
                            face cf(4);
                            for (label i = 0; i < 4; ++i)
                            {
                                cf[i] = nodeLabels[leafI][fNodes[i]];
                            }

                            appendList(helperFaces, cf.reverseFace());
                            helperOwner.append(nei);
                            helperNeighbour.append(leafI);
                        }
                    }
                }
            }

            // store the size of this chunk
            chunkSizes[chunkI] = helperFaces.size() - nBefore;
        }

        // set the sizes of faces graph
        # pragma omp critical
        nFaces += helperFaces.size();

        # pragma omp barrier

        # pragma omp master
        {
            rowSizes.setSize(nFaces);
            octreeFacesPtr_->setSize(nFaces);
            octreeFacesOwnersPtr_->setSize(nFaces);
            octreeFacesNeighboursPtr_->setSize(nFaces);
        }

        # pragma omp barrier

        // set the size of face graph rows and copy owners and neighbours
        for
        (
            label chunkI = threadI;
            chunkI < nChunks;
            chunkI += nThreads
        )
        {
            label start(0), localStart(0);
            for (label i = 0; i < chunkI; ++i)
            {
                start += chunkSizes[i];
            }
            for (label i = threadI; i < chunkI; i += nThreads)
            {
                localStart += chunkSizes[i];
            }

            for (label faceI = 0; faceI < chunkSizes[chunkI]; ++faceI)
            {
                octreeFacesOwnersPtr_->operator[](start) =
                    helperOwner[localStart];
                octreeFacesNeighboursPtr_->operator[](start) =
                    helperNeighbour[localStart];
                rowSizes[start++] = helperFaces[localStart++].size();
            }
        }

        # pragma omp barrier

        // set the size of octree faces
        # pragma omp master
        Foam::Module::setSizeAndRowSize(*octreeFacesPtr_, rowSizes);

        # pragma omp barrier

        // copy the data into octree faces
        for
        (
            label chunkI = threadI;
            chunkI < nChunks;
            chunkI += nThreads
        )
        {
            label start(0), localStart(0);

            for (label i = 0; i < chunkI; ++i)
            {
                start += chunkSizes[i];
            }
            for (label i = threadI; i < chunkI; i += nThreads)
            {
                localStart += chunkSizes[i];
            }

            for (label faceI = 0; faceI < chunkSizes[chunkI]; ++faceI)
            {
                for (label i = 0; i < helperFaces[localStart].size(); ++i)
                {
                    (*octreeFacesPtr_)[start][i] =
                        helperFaces[localStart][i];
                }

                ++start;
                ++localStart;
            }
        }
    }
}


void Foam::Module::meshOctreeAddressing::calculateLeafFaces() const
{
    const DynamicList<label>& owner = octreeFaceOwner();
    const DynamicList<label>& neighbour = octreeFaceNeighbour();

    leafFacesPtr_ = new labelListList2D(octree_.numberOfLeaves());
    labelListList2D& leafFaces = *leafFacesPtr_;

    labelList nlf(leafFaces.size(), 0);
    forAll(owner, fI)
    {
        ++nlf[owner[fI]];
        if (neighbour[fI] < 0)
        {
            continue;
        }
        ++nlf[neighbour[fI]];
    }

    forAll(nlf, leafI)
    {
        leafFaces[leafI].setSize(nlf[leafI]);
    }
    nlf = 0;

    forAll(owner, fI)
    {
        leafFaces[owner[fI]][nlf[owner[fI]]++] = fI;
        if (neighbour[fI] < 0)
        {
            continue;
        }
        leafFaces[neighbour[fI]][nlf[neighbour[fI]]++] = fI;
    }
}


void Foam::Module::meshOctreeAddressing::calculateNodeFaces() const
{
    const labelListList2D& octreeFaces = this->octreeFaces();
    nodeFacesPtr_ = new labelListList2D(numberOfNodes());
    labelListList2D& nodeFaces = *nodeFacesPtr_;

    Foam::Module::reverseAddressing(nodeFaces, numberOfNodes(), octreeFaces);
}


void Foam::Module::meshOctreeAddressing::calculateLeafLeaves() const
{
    const DynamicList<label>& owner = octreeFaceOwner();
    const DynamicList<label>& neighbour = octreeFaceNeighbour();

    leafLeavesPtr_ = new labelListList2D(octree_.numberOfLeaves());
    labelListList2D& leafLeaves = *leafLeavesPtr_;

    labelList nNei(leafLeaves.size(), 0);
    forAll(owner, faceI)
    {
        if (owner[faceI] < 0)
        {
            continue;
        }
        if (neighbour[faceI] < 0)
        {
            continue;
        }

        ++nNei[owner[faceI]];
        ++nNei[neighbour[faceI]];
    }

    forAll(nNei, leafI)
    {
        leafLeaves[leafI].setSize(nNei[leafI]);
    }

    nNei = 0;

    forAll(owner, faceI)
    {
        if (owner[faceI] < 0)
        {
            continue;
        }
        if (neighbour[faceI] < 0)
        {
            continue;
        }

        leafLeaves[owner[faceI]][nNei[owner[faceI]]++] = neighbour[faceI];
        leafLeaves[neighbour[faceI]][nNei[neighbour[faceI]]++] = owner[faceI];
    }
}


void Foam::Module::meshOctreeAddressing::createOctreeEdges() const
{
    const labelListList2D& faces = this->octreeFaces();

    // allocate memory for edges, face-edges addressing
    // and node-edges addressing
    octreeEdgesPtr_ = new DynamicList<edge>();
    DynamicList<edge>& edges = *octreeEdgesPtr_;
    faceEdgesPtr_ = new labelListList2D(faces.size());
    labelListList2D& faceEdges = *faceEdgesPtr_;
    nodeEdgesPtr_ = new labelListList2D();
    labelListList2D& nodeEdges = *nodeEdgesPtr_;
    nodeEdges.setSize(nNodes_);

    forAll(faces, faceI)
    {
        faceEdges[faceI].setSize(faces[faceI].size());
        forAllRow(faceEdges, faceI, feI)
        {
            faceEdges[faceI][feI] = -1;
        }
    }

    forAll(faces, faceI)
    {
        const label nEdges = faces[faceI].size();

        for (label eI = 0; eI < nEdges; ++eI)
        {
            const edge e
            (
                faces[faceI][eI],
                faces[faceI][(eI + 1)%nEdges]
            );

            label eLabel(-1);
            forAllRow(nodeEdges, e.start(), neI)
            {
                if (edges[nodeEdges[e.start()][neI]] == e )
                {
                    eLabel = nodeEdges[e.start()][neI];
                    break;
                }
            }

            if (eLabel < 0)
            {
                // append new edge
                faceEdges[faceI][eI] = edges.size();
                nodeEdges[e.start()].append(edges.size());
                nodeEdges[e.end()].append(edges.size());

                edges.append(e);
            }
            else
            {
                faceEdges[faceI][eI] = eLabel;
            }
        }
    }
}


void Foam::Module::meshOctreeAddressing::calculateLeafEdges() const
{
    const labelListList2D& edgeLeaves = this->edgeLeaves();

    leafEdgesPtr_ = new labelListList2D();
    labelListList2D& leafEdges = *leafEdgesPtr_;

    Foam::Module::reverseAddressing
    (
        leafEdges,
        octree_.numberOfLeaves(),
        edgeLeaves
    );
}


void Foam::Module::meshOctreeAddressing::calculateEdgeLeaves() const
{
    const labelListList2D& edgeFaces = this->edgeFaces();
    const DynamicList<label>& owner = this->octreeFaceOwner();
    const DynamicList<label>& neighbour = this->octreeFaceNeighbour();

    edgeLeavesPtr_ = new labelListList2D();
    labelListList2D& edgeLeaves = *edgeLeavesPtr_;
    edgeLeaves.setSize(edgeFaces.size());

    forAll(edgeFaces, edgeI)
    {
        forAllRow(edgeFaces, edgeI, efI)
        {
            const label fI = edgeFaces[edgeI][efI];
            const label own = owner[fI];
            const label nei = neighbour[fI];

            Foam::Module::appendUniq(edgeLeaves[edgeI], own);

            if (nei < 0)
            {
                continue;
            }
            Foam::Module::appendUniq(edgeLeaves[edgeI], nei);
        }
    }
}


void Foam::Module::meshOctreeAddressing::calculateEdgeFaces() const
{
    const labelListList2D& faceEdges = this->faceEdges();
    edgeFacesPtr_ = new labelListList2D(octreeEdges().size());
    labelListList2D& edgeFaces = *edgeFacesPtr_;

    Foam::Module::reverseAddressing
    (
        edgeFaces,
        octreeEdges().size(),
        faceEdges
    );
}


// ************************************************************************* //
