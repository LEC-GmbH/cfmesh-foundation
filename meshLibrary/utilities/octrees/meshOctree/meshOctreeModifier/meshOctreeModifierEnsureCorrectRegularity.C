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
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeModifier::ensureCorrectRegularity
(
    labelList& refineBox
)
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    labelLongList front;
    forAll(refineBox, leafI)
    {
        if (refineBox[leafI])
            front.append(leafI);
    }

    DynList<label> neighbours;

    label nMarked;
    do
    {
        nMarked = 0;

        # pragma omp parallel private(neighbours)
        {
            labelLongList tFront;

            # pragma omp for schedule(dynamic, 50)
            forAll(front, i)
            {
                const label leafI = front[i];
                const meshOctreeCube* oc = leaves[leafI];

                neighbours.clear();
                octree_.findAllLeafNeighbours(*oc, neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if (nei < 0)
                        continue;

                    const meshOctreeCube& noc = *leaves[nei];

                    if (noc.level() >= oc->level())
                        continue;

                    if (!refineBox[nei])
                    {
                        refineBox[nei] = 1;
                        tFront.append(nei);
                    }
                }
            }

            # pragma omp barrier

            # pragma omp single
            front.clear();

            # pragma omp barrier

            # pragma omp critical
            {
                label start = front.size();
                front.setSize(start + tFront.size());

                forAll(tFront, i)
                    front[start++] = tFront[i];
            }
        }

        nMarked = front.size();
    }
    while (nMarked != 0);
}


bool Foam::Module::meshOctreeModifier::ensureCorrectRegularitySons
(
    labelList& refineBox
)
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    label nMarked(0);

    # pragma omp parallel for schedule(dynamic, 100) reduction(+ : nMarked)
    forAll(leaves, leafI)
    {
        if (!refineBox[leafI])
            continue;

        const meshOctreeCubeCoordinates cc = leaves[leafI]->reduceLevelBy(1);

        for (label scI = 0; scI < 8; ++scI)
        {
            const label neiLeaf =
                octree_.findLeafLabelForPosition(cc.refineForPosition(scI));

            if (neiLeaf >= 0 && !refineBox[neiLeaf])
            {
                // mark this leaf for refinement
                ++nMarked;
                refineBox[neiLeaf] = 1;
            }
        }
    }

    if (nMarked)
        return true;

    return false;
}


// ************************************************************************* //
