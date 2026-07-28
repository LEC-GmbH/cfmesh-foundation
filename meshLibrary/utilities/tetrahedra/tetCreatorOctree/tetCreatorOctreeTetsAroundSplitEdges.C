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

#include "tetCreatorOctree.H"
#include "demandDrivenData.H"
#include "meshOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::tetCreatorOctree::createTetsAroundSplitEdges()
{
    Info<< "Creating tets around split edges " << endl;

    const labelList& cubeLabel = *cubeLabelPtr_;
    const meshOctree& octree = octreeCheck_.octree();
    const labelListList2D& nodeLabels = octreeCheck_.nodeLabels();
    const List<FixedList<label, 8>>& pLeaves = octreeCheck_.nodeLeaves();
    const labelListList2D& subNodeLabels = *subNodeLabelsPtr_;
    const labelListList2D& faceCentreLabel = *faceCentreLabelPtr_;

    // find maximum refinement level of octree leaves attached to each vertex
    List<direction> nodeLevel(octreeCheck_.numberOfNodes());

    forAll(pLeaves, nodeI)
    {
        direction level(0);

        for (label plI = 0; plI < 8; ++plI)
        {
            const label leafI = pLeaves[nodeI][plI];

            if (leafI < 0)
            {
                continue;
            }

            level = Foam::max(level, octree.returnLeaf(leafI).level());
        }

        nodeLevel[nodeI] = level;
    }

    // start creating tets around split edges
    label helpNodes[2][8];
    label faceCentres[4];

    forAllReverse(sortedLeaves_, levelI)
    {
        const DynamicList<label>& curLevelLeaves = sortedLeaves_[levelI];

        const direction level = direction(levelI);

        forAll(curLevelLeaves, leafI)
        {
            const label curLabel = curLevelLeaves[leafI];

            if (cubeLabel[curLabel] == -1)
            {
                continue;
            }

            // start checking edges
            for (label eI = 0; eI < 12; ++eI)
            {
                const label start =
                    nodeLabels[curLabel]
                    [
                        meshOctreeCubeCoordinates::edgeNodes_[eI][0]
                    ];
                const label end =
                    nodeLabels[curLabel]
                    [
                        meshOctreeCubeCoordinates::edgeNodes_[eI][1]
                    ];

                if ((nodeLevel[start] == level) && (nodeLevel[end] == level))
                {
                    continue;
                }

                // the edge has at least one vertex at different ref level
                bool create(true);

                FixedList<label, 4> edgeCubes;
                const label fI = 2*(eI/4)+1;

                const label* fNodes =
                    meshOctreeCubeCoordinates::faceNodes_[fI];

                // store octree leaves at this edge
                // they are all adjacent to the start point
                for (label i = 0; i < 4; ++i)
                {
                    edgeCubes[i] = pLeaves[start][fNodes[i]];
                }

                forAll(edgeCubes, i)
                {
                    const label cLabel = edgeCubes[i];

                    if
                    (
                        (cLabel == -1)
                     || (cLabel < curLabel)
                     || (octree.returnLeaf(cLabel).level() != level)
                    )
                    {
                        create = false;
                        break;
                    }

                    for (label j = 0; j < 8; ++j)
                    {
                        if (nodeLabels[cLabel][j] == start)
                        {
                            if (subNodeLabels[cLabel].size() != 0)
                            {
                                helpNodes[0][i] = subNodeLabels[cLabel][j];
                            }
                            else
                            {
                                helpNodes[0][i] = -1;
                            }

                            helpNodes[0][i + 4] = cubeLabel[cLabel];
                        }

                        if (nodeLabels[cLabel][j] == end)
                        {
                            if (subNodeLabels[cLabel].size() != 0)
                            {
                                helpNodes[1][i + 4] = subNodeLabels[cLabel][j];
                            }
                            else
                            {
                                helpNodes[1][i + 4] = -1;
                            }

                            helpNodes[1][i] = cubeLabel[cLabel];
                        }
                    }

                    if (faceCentreLabel[cLabel].size() != 0)
                    {
                        const label helpFace = eI/4;

                        faceCentres[i] =
                            faceCentreLabel
                            [cLabel]
                            [faceCentreHelper_[helpFace][i]];
                    }
                    else
                    {
                        faceCentres[i] = -1;
                    }
                }

                if (!create)
                {
                    continue;
                }

                if (nodeLevel[start] > level)
                {
                    for (label k = 0; k < 4; ++k)
                    {
                        // add 4 tets
                        checkAndAppendTet
                        (
                            partTet
                            (
                                start,
                                helpNodes[0][k],
                                helpNodes[0][(k + 1)%4],
                                tetPoints_.size()
                            )
                        );

                        // 2. tet
                        checkAndAppendTet
                        (
                            partTet
                            (
                                tetPoints_.size(),
                                helpNodes[0][k],
                                helpNodes[0][(k + 1)%4],
                                faceCentres[k]
                            )
                        );

                        // 3. tet
                        checkAndAppendTet
                        (
                            partTet
                            (
                                faceCentres[k],
                                helpNodes[0][k],
                                helpNodes[0][k + 4],
                                tetPoints_.size()
                            )
                        );

                        // 4. tet
                        checkAndAppendTet
                        (
                            partTet
                            (
                                faceCentres[(k + 3)%4],
                                helpNodes[0][k + 4],
                                helpNodes[0][k],
                                tetPoints_.size()
                            )
                        );
                    }
                }
                else
                {
                    for (label k = 0; k < 4; ++k)
                    {
                        checkAndAppendTet
                        (
                            partTet
                            (
                                faceCentres[(k + 3)%4],
                                helpNodes[0][k + 4],
                                start,
                                tetPoints_.size()
                            )
                        );

                        checkAndAppendTet
                        (
                            partTet
                            (
                                faceCentres[k],
                                start,
                                helpNodes[0][k + 4],
                                tetPoints_.size()
                            )
                        );
                    }
                }

                if (nodeLevel[end] > level)
                {
                    for (label k = 0; k < 4; ++k)
                    {
                        // add 4 tets
                        checkAndAppendTet
                        (
                            partTet
                            (
                                tetPoints_.size(),
                                helpNodes[1][k + 4],
                                helpNodes[1][((k + 1)%4)+4],
                                end
                            )
                        );

                        // 2. tet
                        checkAndAppendTet
                        (
                            partTet
                            (
                                faceCentres[k],
                                helpNodes[1][k + 4],
                                helpNodes[1][((k + 1)%4)+4],
                                tetPoints_.size()
                            )
                        );

                        // 3. tet
                        checkAndAppendTet
                        (
                            partTet
                            (
                                tetPoints_.size(),
                                helpNodes[1][k],
                                faceCentres[k],
                                helpNodes[1][k + 4]
                            )
                        );

                        // 4. tet
                        checkAndAppendTet
                        (
                            partTet
                            (
                                faceCentres[(k + 3)%4],
                                helpNodes[1][k],
                                tetPoints_.size(),
                                helpNodes[1][k + 4]
                            )
                        );
                    }
                }
                else
                {
                    for (label k = 0; k < 4; ++k)
                    {
                        checkAndAppendTet
                        (
                            partTet
                            (
                                faceCentres[(k + 3)%4],
                                helpNodes[1][k],
                                tetPoints_.size(),
                                end
                            )
                        );

                        checkAndAppendTet
                        (
                            partTet
                            (
                                helpNodes[1][k],
                                faceCentres[k],
                                tetPoints_.size(),
                                end
                            )
                        );
                    }
                }

                // add the edge centre
                tetPoints_.append(0.5 *(tetPoints_[start] + tetPoints_[end]));
            }
        }
    }
}


// ************************************************************************* //
