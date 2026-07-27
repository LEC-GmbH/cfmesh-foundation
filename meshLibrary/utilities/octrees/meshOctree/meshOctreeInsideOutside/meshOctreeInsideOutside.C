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

#include "meshOctreeInsideOutside.H"
#include "triSurf.H"
#include "boundBox.H"
#include "labelPair.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Module::meshOctreeInsideOutside::meshOctreeInsideOutside
(
    meshOctree& octree
)
:
    octreeModifier_(octree),
    cubeGroup_(octree.numberOfLeaves(), -1),
    cubesInGroup_(),
    groupType_(),
    boundaryDATACubes_(),
    hasOutsideNeighbour_(octree.numberOfLeaves(), false),
    neighbouringGroups_()
{
    initialiseBoxes();

    frontalMarking();

    markOutsideCubes();

    reviseDataBoxes();

    markInsideCubes();

    label nInternal(0), nUnknown(0), nData(0), nOutside(0);

    const label nLeaves = octree.numberOfLeaves();
    for (label leafI = 0; leafI < nLeaves; ++leafI)
    {
        const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);

        if (oc.cubeType() & meshOctreeCube::INSIDE)
        {
            ++nInternal;
        }
        else if (oc.cubeType() & meshOctreeCube::UNKNOWN)
        {
            ++nUnknown;
        }
        else if (oc.cubeType() & meshOctreeCube::DATA)
        {
            ++nData;
        }
        else if (oc.cubeType() & meshOctreeCube::OUTSIDE)
        {
            ++nOutside;
        }
    }

    Info<< "Number of internal boxes is " << nInternal << endl;
    Info<< "Number of outside boxes is " << nOutside << endl;
    Info<< "Number of data boxes is " << nData << endl;
    Info<< "Number of unknown boxes is " << nUnknown << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeInsideOutside::initialiseBoxes()
{
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();

    # pragma omp parallel for if (leaves.size() > 1000)
    forAll(leaves, leafI)
    {
        if (leaves[leafI]->hasContainedElements())
        {
            leaves[leafI]->setCubeType(meshOctreeCubeBasic::DATA);
        }
        else
        {
            leaves[leafI]->setCubeType(meshOctreeCubeBasic::UNKNOWN);
        }
    }
}


void Foam::Module::meshOctreeInsideOutside::frontalMarking()
{
    neighbouringGroups_.clear();

    DynamicList<label> frontCubes;
    DynamicList<label> neighbours;

    label nGroup(0), nThreads(1);

    const DynamicList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();

    if (leaves.size() > 1000)
        nThreads = 3*omp_get_num_procs();

    # pragma omp parallel num_threads(nThreads) \
    private(frontCubes, neighbours)
    {
        DynamicList<labelPair> threadCommPairs;

        const label threadI = omp_get_thread_num();

        const label chunkSize = leaves.size() / nThreads + 1;

        const label minLeaf = threadI*chunkSize;

        const label maxLeaf = Foam::min(leaves.size(), minLeaf + chunkSize);

        for (label leafI = minLeaf; leafI < maxLeaf; ++leafI)
        {
            if (leaves[leafI]->hasContainedElements())
                continue;
            if (cubeGroup_[leafI] != -1)
                continue;

            label groupI;
            # pragma omp critical
            groupI = nGroup++;

            direction cType(meshOctreeCubeBasic::UNKNOWN);
            frontCubes.clear();
            frontCubes.append(leafI);
            cubeGroup_[leafI] = groupI;

            DynamicList<label> neiDATACubes;

            while (frontCubes.size())
            {
                const label fLabel = frontCubes.remove();
                octree.findNeighboursForLeaf(fLabel, neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];
                    if ((nei >= minLeaf) && (nei < maxLeaf))
                    {
                        if (cubeGroup_[nei] != -1)
                            continue;

                        if (leaves[nei]->hasContainedElements())
                        {
                            neiDATACubes.append(nei);
                        }
                        else
                        {
                            frontCubes.append(nei);
                            cubeGroup_[nei] = groupI;
                        }
                    }
                    else if (nei == -1)
                    {
                        cType = meshOctreeCubeBasic::OUTSIDE;
                    }
                    else
                    {
                        if (leaves[nei]->hasContainedElements())
                        {
                            neiDATACubes.append(nei);
                        }
                        else
                        {
                            threadCommPairs.append
                            (
                                labelPair(fLabel, nei)
                            );
                        }
                    }
                }
            }

            # pragma omp critical
            {
                if (groupI >= boundaryDATACubes_.size())
                    boundaryDATACubes_.setSize(groupI + 1);

                setRow(boundaryDATACubes_, groupI, neiDATACubes);
                groupType_(groupI) = cType;
            }
        }

        # pragma omp barrier

        // find group to neighbouring groups addressing
        labelListList2D localNeiGroups(nGroup);
        forAll(threadCommPairs, cfI)
        {
            const labelPair& lp = threadCommPairs[cfI];
            const label groupI = cubeGroup_[lp.first()];
            const label neiGroup = cubeGroup_[lp.second()];

            if ((neiGroup >= nGroup) || (groupI >= nGroup))
                FatalError << "neiGroup " << neiGroup
                    << " groupI " << groupI << " are >= than "
                    << "nGroups " << nGroup << abort(FatalError);

            if (neiGroup != -1)
            {
                appendUniq(localNeiGroups[groupI], neiGroup);
                appendUniq(localNeiGroups[neiGroup], groupI);
            }
        }

        # pragma omp critical
        {
            neighbouringGroups_.setSize(nGroup);

            forAll(localNeiGroups, groupI)
            {
                const DynamicList<label>& lGroups = localNeiGroups[groupI];

                appendUniq(neighbouringGroups_[groupI], groupI);

                forAll(lGroups, i)
                {
                    neighbouringGroups_[groupI].append(lGroups[i]);
                }
            }
        }
    }

    // create cubesInGroup_ addressing
    labelList nCubesInGroup(nGroup, 0);
    forAll(cubeGroup_, leafI)
    {
        if (cubeGroup_[leafI] < 0)
            continue;

        ++nCubesInGroup[cubeGroup_[leafI]];
    }

    setSizeAndRowSize(cubesInGroup_, nCubesInGroup);

    forAllReverse(cubeGroup_, leafI)
    {
        const label groupI = cubeGroup_[leafI];

        if (groupI < 0)
            continue;

        cubesInGroup_[groupI][--nCubesInGroup[groupI]] = leafI;
    }

}


void Foam::Module::meshOctreeInsideOutside::markOutsideCubes()
{
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();

    DynamicList<label> neighbours;
    label nChanged;
    bool keepUpdating;

    do
    {
        keepUpdating = false;

        do
        {
            nChanged = 0;

            // make sure that groups created by different threads
            // have the same information
            forAll(neighbouringGroups_, groupI)
            {
                if (groupType_[groupI] & meshOctreeCubeBasic::OUTSIDE)
                {
                    forAllRow(neighbouringGroups_, groupI, i)
                    {
                        const label neiGroup = neighbouringGroups_[groupI][i];
                        if (groupType_[neiGroup] & meshOctreeCube::UNKNOWN)
                        {
                            ++nChanged;
                            groupType_(neiGroup) = meshOctreeCube::OUTSIDE;
                        }
                    }
                }
            }

            if (nChanged != 0)
                keepUpdating = true;

        } while (nChanged != 0);

    } while (keepUpdating);

    // set OUTSIDE type to the cubes in OUTSIDE groups
    forAllConstIter(Map<direction>, groupType_, it)
    {
        if (it.key() < 0)
            continue;

        if (it() & meshOctreeCubeBasic::OUTSIDE)
        {
            const label groupI = it.key();

            // set the cube type to OUTSIDE
            forAllRow(cubesInGroup_, groupI, i)
                leaves[cubesInGroup_[groupI][i]]->setCubeType
                (
                    meshOctreeCube::OUTSIDE
                );

            // set true to the collected DATA boxes
            forAllRow(boundaryDATACubes_, groupI, neiI)
                hasOutsideNeighbour_[boundaryDATACubes_[groupI][neiI]] = true;
        }
    }
}


void Foam::Module::meshOctreeInsideOutside::reviseDataBoxes()
{
    // remove DATA flag from boxes which do not have an OUTSIDE neighbour
    // and are not surrounded with DATA boxes containing different surface
    // triangles in different patches
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();
    const triSurf& surface = octree.surface();
    DynamicList<label> neighbours;

    label nMarked;

    do
    {
        nMarked = 0;

        # pragma omp parallel for if (leaves.size() > 1000) \
        private(neighbours) schedule(dynamic, 20) reduction(+ : nMarked)
        forAll(leaves, leafI)
            if (
                (leaves[leafI]->cubeType() & meshOctreeCube::DATA) &&
                !hasOutsideNeighbour_[leafI]
            )
            {
                meshOctreeCube* oc = leaves[leafI];

                DynamicList<label> patches;
                const labelListList2D& ct =
                    oc->slotPtr()->containedTriangles_;
                const DynamicList<label>& el = ct[oc->containedElements()];
                forAll(el, elI)
                    appendUniq(patches, surface[el[elI]].region());

                if (patches.size() > 1)
                    continue;

                // check if there exist neighbours
                // which have some DATA neighbours
                octree.findAllLeafNeighbours(leafI, neighbours);
                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if (nei < 0)
                        continue;

                    if (hasOutsideNeighbour_[nei])
                    {
                        oc->setCubeType(meshOctreeCube::INSIDE);

                        ++nMarked;
                        break;
                    }
                }
            }

    } while (nMarked != 0);

}


void Foam::Module::meshOctreeInsideOutside::markInsideCubes()
{
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    label nChanged;
    bool keepUpdating;
    DynamicList<label> neighbours;

    // make INSIDE groups for which it is possible
    forAllIter(Map<direction>, groupType_, it)
    {
        const label groupI = it.key();

        if (groupI < 0)
            continue;

        if (it() & meshOctreeCubeBasic::UNKNOWN)
        {
            forAllRow(boundaryDATACubes_, groupI, neiI)
            {
                const label cLabel = boundaryDATACubes_[groupI][neiI];
                if
                (
                    hasOutsideNeighbour_[cLabel]
                 || (leaves[cLabel]->cubeType() & meshOctreeCube::INSIDE)
                )
                {
                    it() = meshOctreeCube::INSIDE;
                    break;
                }
            }
        }
    }

    do
    {
        keepUpdating = false;

        // mark INSIDE groups created by different threads
        do
        {
            nChanged = 0;

            forAll(neighbouringGroups_, groupI)
            {
                if (groupType_[groupI] & meshOctreeCube::INSIDE)
                {
                    forAllRow(neighbouringGroups_, groupI, i)
                    {
                        const label neiGroup = neighbouringGroups_[groupI][i];

                        if (groupType_[neiGroup] & meshOctreeCube::UNKNOWN)
                        {
                            ++nChanged;
                            groupType_(neiGroup) = meshOctreeCube::INSIDE;
                        }
                    }
                }
            }

            if (nChanged != 0)
                keepUpdating = true;

        } while (nChanged != 0);

    } while (keepUpdating);

    // set INSIDE type to the cubes in INSIDE groups
    forAllConstIter(Map<direction>, groupType_, it)
    {
        if (it.key() < 0)
            continue;

        if (it() & meshOctreeCubeBasic::INSIDE)
        {
            const label groupI = it.key();

            // set the cube type to OUTSIDE
            forAllRow(cubesInGroup_, groupI, i)
                leaves[cubesInGroup_[groupI][i]]->setCubeType
                (
                    meshOctreeCube::INSIDE
                );
        }
    }
}


// ************************************************************************* //
