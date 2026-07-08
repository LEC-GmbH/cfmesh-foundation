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

#include "meshOctreeModifier.H"
#include "triSurf.H"
#include "HashSet.H"
#include "helperFunctions.H"

#include <omp.h>

#include <set>

#include <sys/stat.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeModifier::markAdditionalLayers
(
    labelList& refineBox,
    const label nLayers
) const
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    DynList<label> neiLeaves;

    for (label i = 1; i<=nLayers; ++i)
    {
        labelLongList activeLeaves;
        forAll(leaves, leafI)
            if (refineBox[leafI] == i)
                activeLeaves.append(leafI);

        # pragma omp parallel for private(neiLeaves) schedule(dynamic, 20)
        forAll(activeLeaves, lI)
        {
            const label leafI = activeLeaves[lI];

            const meshOctreeCubeCoordinates& oc = leaves[leafI]->coordinates();

            neiLeaves.clear();
            octree_.findAllLeafNeighbours(oc, neiLeaves);

            forAll(neiLeaves, posI)
            {
                const label neiLabel = neiLeaves[posI];

                if (neiLabel < 0)
                    continue;

                if (!refineBox[neiLabel])
                    refineBox[neiLabel] = i + 1;
            }
        }

    }
}


void Foam::Module::meshOctreeModifier::markAdditionalLayersOfFaceNeighbours
(
    labelList& refineBox,
    const label nLayers
) const
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    DynList<label> neiLeaves;

    for (label i = 1; i<=nLayers; ++i)
    {
        labelLongList activeLeaves;
        forAll(leaves, leafI)
            if (refineBox[leafI] == i)
                activeLeaves.append(leafI);

        # pragma omp parallel for private(neiLeaves) schedule(dynamic, 20)
        forAll(activeLeaves, lI)
        {
            const label leafI = activeLeaves[lI];

            const meshOctreeCubeCoordinates& oc = leaves[leafI]->coordinates();

            neiLeaves.clear();
            octree_.findNeighboursForLeaf(oc, neiLeaves);

            forAll(neiLeaves, posI)
            {
                const label neiLabel = neiLeaves[posI];

                if (neiLabel < 0)
                    continue;

                if (!refineBox[neiLabel])
                    refineBox[neiLabel] = i + 1;
            }
        }

    }
}


Foam::label Foam::Module::meshOctreeModifier::markAdditionalLayers
(
    labelList& refineBox,
    labelList& nLayers,
    List<direction>& targetRefLevel
) const
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    // sort leaves based on the number of additional layers
    label maxLevel = Foam::max(nLayers);

    List<labelLongList> leavesForLayer(maxLevel + 1);
    forAll(nLayers, leafI)
    {
        if (nLayers[leafI] < 1)
            continue;

        leavesForLayer[nLayers[leafI]].append(leafI);
    }

    // set refinement flag to additional boxes marked separately
    label nMarked(0);

    forAllReverse(leavesForLayer, layerI)
    {
        const labelLongList& activeLeaves = leavesForLayer[layerI];

        if (activeLeaves.size() == 0 )
            continue;

        // find the max required refinement level
        direction maxLevel(0);

        # pragma omp parallel
        {
            direction localMax(0);

            # pragma omp for schedule(dynamic, 50)
            forAll(activeLeaves, i)
                localMax = Foam::max(localMax, targetRefLevel[activeLeaves[i]]);

            # pragma omp critical
            maxLevel = Foam::max(localMax, maxLevel);
        }

        // mark additional boxes for refinement
        for (direction levelI = maxLevel; levelI > 0; --levelI)
        {
            labelList markedBoxes(leaves.size(), 0);

            label counter(0);
            # pragma omp parallel for reduction(+:counter)
            forAll(activeLeaves, lI)
            {
                const label leafI = activeLeaves[lI];

                if (targetRefLevel[leafI] == levelI)
                {
                    markedBoxes[leafI] = 1;
                    ++counter;
                }
            }

            if (counter == 0)
                continue;

            // mark additional cells at this refinement level
            markAdditionalLayersOfFaceNeighbours(markedBoxes, layerI);

            // update the main list
            # pragma omp parallel for schedule(dynamic, 100) \
            reduction(+:nMarked)
            forAll(markedBoxes, leafI)
            {
                if (markedBoxes[leafI] < 2)
                    continue;

                if (leaves[leafI]->level() >= levelI)
                    continue;

                if (!refineBox[leafI])
                {
                    refineBox[leafI] |= 1;
                    ++nMarked;
                }
            }
        }
    }

    return nMarked;
}


void Foam::Module::meshOctreeModifier::refineSelectedBoxes
(
    labelList& refineBox,
    const bool hexRefinement
)
{
    // ensure that refinement will produce 1-irregular octree
    do
    {
        ensureCorrectRegularity(refineBox);
    } while (hexRefinement && ensureCorrectRegularitySons(refineBox));

    const triSurf& surface = octree_.surface();
    const boundBox& rootBox = octree_.rootBox();
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    // this is needed for thread safety
    // such solutions make me a sad bunny :(
    surface.facetEdges();
    surface.edgeFacets();
    surface.edges();

    # pragma omp parallel num_threads(octree_.dataSlots_.size())
    {
        meshOctreeSlot* slotPtr = &octree_.dataSlots_[omp_get_thread_num()];

        if (!octree_.isQuadtree())
        {
            # pragma omp for schedule(dynamic, 100)
            forAll(leaves, leafI)
            {
                if (refineBox[leafI])
                    leaves[leafI]->refineCube(surface, rootBox, slotPtr);
            }
        }
        else
        {
            # pragma omp for schedule(dynamic, 100)
            forAll(leaves, leafI)
            {
                if (refineBox[leafI])
                    leaves[leafI]->refineCube2D(surface, rootBox, slotPtr);
            }
        }
    }

    // update the list of leaves
    createListOfLeaves();

}


void Foam::Module::meshOctreeModifier::refineSelectedBoxesAndAdditionalLayers
(
    labelList& refineBox,
    const scalarList& refThickness
)
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    typedef std::pair<direction, label> mapKey;
    typedef std::map<mapKey, LongList<meshOctreeCube*>> lMap;
    lMap leavesMap;

    // find the maximum refinement level of leaves marked for refinement
    direction maxLevel(0);

    # pragma omp parallel
    {
        direction localMax(0);

        forAll(refineBox, leafI)
        {
            if (refineBox[leafI])
                localMax = Foam::max(localMax, leaves[leafI]->level());
        }

        # pragma omp critical
        maxLevel = Foam::max(maxLevel, localMax);
    }

    // sort leaves based on the current level
    List<labelLongList> leavesForLevel(maxLevel + 1);
    forAll(refineBox, leafI)
    {
        if (!refineBox[leafI])
            continue;

        leavesForLevel[leaves[leafI]->level()].append(leafI);
    }

    // find leaves with the same number of additional layers
    forAllReverse(leavesForLevel, levelI)
    {
        const labelLongList& activeLeaves = leavesForLevel[levelI];

        if (activeLeaves.size() == 0 )
            continue;

        labelLongList nLayersForActive(activeLeaves.size());
        label maxNLayers(0);
        # pragma omp parallel
        {
            label localMaxNLayers(0);

            forAll(activeLeaves, i)
            {
                const label leafI = activeLeaves[i];

                const scalar cs = leaves[leafI]->size(octree_.rootBox());

                nLayersForActive[i] = Foam::max(ceil(refThickness[leafI]/cs), 1);

                localMaxNLayers =
                    Foam::max(localMaxNLayers, nLayersForActive[i]);
            }

            # pragma omp critical
            maxNLayers = Foam::max(maxNLayers, localMaxNLayers);
        }

        for (label layerI = 0; layerI<=maxNLayers; ++layerI)
        {
            mapKey key(levelI, layerI);
            LongList<meshOctreeCube*>& currLeaves = leavesMap[key];
            currLeaves.clear();
        }

        forAll(activeLeaves, i)
        {
            mapKey key(levelI, nLayersForActive[i]);
            leavesMap[key].append(leaves[activeLeaves[i]]);
        }
    }

    // refine leaves
    forAllConstIters(leavesMap, it)
    {
        const label nLayers = it->first.second;
        const direction levelI = it->first.first;

        const LongList<meshOctreeCube*>& selectedLeaves = it->second;

        if (selectedLeaves.size() == 0 )
            continue;

        label nMarked;
        do
        {
            nMarked = 0;

            // mark current leaves for refinement
            labelList markedLeaves(leaves.size(), 0);

            # pragma omp parallel for schedule(dynamic, 50) reduction(+:nMarked)
            forAll(selectedLeaves, i)
            {
                if (!selectedLeaves[i]->isLeaf())
                    continue;

                markedLeaves[selectedLeaves[i]->cubeLabel()] = 1;
                ++nMarked;
            }

            // get out of the do-while loop if there are no selected leaves
            if (nMarked == 0)
                break;

            // mark additional boxes for refinement
            markAdditionalLayers(markedLeaves, nLayers);

            // find the leaves in the additional layers
            labelLongList activeLeaves;
            forAll(markedLeaves, leafI)
            {
                if (markedLeaves[leafI])
                    activeLeaves.append(leafI);
            }

            // check if there exist leaves at lower refinement level
            bool hasLowerLevel(false);

            # pragma omp parallel for schedule(guided)
            forAll(activeLeaves, i)
            {
                const direction level = leaves[activeLeaves[i]]->level();
                if (level < levelI)
                {
                    // found a neighbour at a lower refinement level
                    hasLowerLevel = true;
                }
                else if (level > levelI)
                {
                    // do not allow refinement of leaves at higher
                    // refinement level
                    markedLeaves[activeLeaves[i]] = 0;
                }
            }

            // deselect leaves at the current level
            if (hasLowerLevel)
            {
                # pragma omp parallel for schedule(guided)
                forAll(activeLeaves, i)
                {
                    if (leaves[activeLeaves[i]]->level() == levelI)
                        markedLeaves[activeLeaves[i]] = 0;
                }
            }

            // refine selected octree boxes
            refineSelectedBoxes(markedLeaves);
        } while (nMarked);
    }
}


// ************************************************************************* //
