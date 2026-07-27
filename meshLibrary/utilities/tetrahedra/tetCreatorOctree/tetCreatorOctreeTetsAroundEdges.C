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
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::tetCreatorOctree::createTetsAroundEdges()
{
    Info<< "Creating tets around edges" << endl;

    const labelList& cubeLabel = *cubeLabelPtr_;
    const labelListList2D& nodeLabels = octreeCheck_.nodeLabels();
    const List<FixedList<label, 8>>& pLeaves = octreeCheck_.nodeLeaves();
    const labelListList2D& subNodeLabels = *subNodeLabelsPtr_;
    const labelListList2D& faceCentreLabel = *faceCentreLabelPtr_;
    const meshOctree& octree = octreeCheck_.octree();
    const FixedList<FixedList<meshOctreeCubeCoordinates, 8>, 8>& vlPos =
        octree.positionsOfLeavesAtNodes();

    // find maximum refinement level of octree leaves attached to each vertex
    List<direction> nodeLevel(octreeCheck_.numberOfNodes());

    forAll(pLeaves, nodeI)
    {
        direction level(0);

        for (label plI = 0; plI < 8; ++plI)
        {
            const label leafI = pLeaves[nodeI][plI];

            if (leafI < 0)
                continue;

            level = Foam::max(level, octree.returnLeaf(leafI).level());
        }

        nodeLevel[nodeI] = level;
    }

    // start creating tets around edges which have both vertices at the same
    // refinement level which is equal to the max refinement level of boxes
    // incident to such edges
    forAllReverse(sortedLeaves_, levelI)
    {
        const DynamicList<label>& curLevelLeaves = sortedLeaves_[levelI];

        forAll(curLevelLeaves, leafI)
        {
            const label curLeaf = curLevelLeaves[leafI];

            if (cubeLabel[curLeaf] == -1)
                continue;

            const meshOctreeCubeCoordinates& cc =
                octree.returnLeaf(curLeaf).coordinates();

            // start checking edges
            for (label eI = 0; eI < 12; ++eI)
            {
                const label startNode =
                    meshOctreeCubeCoordinates::edgeNodes_[eI][0];

                const label start = nodeLabels[curLeaf][startNode];
                const label end =
                    nodeLabels[curLeaf]
                    [
                        meshOctreeCubeCoordinates::edgeNodes_[eI][1]
                    ];

                bool create(true);

                if
                (
                    (nodeLevel[start] == direction(levelI))
                 && (nodeLevel[end] == direction(levelI))
                )
                {
                    // edge has both vertices at the same refinement level
                    // as the current leaf
                    FixedList<label, 4> edgeCubes;
                    const label fI = 2*(eI/4)+1;

                    const label* fNodes =
                        meshOctreeCubeCoordinates::faceNodes_[fI];

                    // store octree leaves at this edge
                    // they are all adjacent to the start point
                    for (label i = 0; i < 4; ++i)
                        edgeCubes[i] = pLeaves[start][fNodes[i]];

                    DynamicList<label> centreNodes;

                    forAll(edgeCubes, i)
                    {
                        const label cLabel = edgeCubes[i];

                        if ((cLabel == -1) || (cubeLabel[cLabel] == -1))
                        {
                            centreNodes.append(-1);
                            continue;
                        }

                        const meshOctreeCubeCoordinates& oc =
                            octree.returnLeaf(cLabel).coordinates();

                        if (oc.level() == direction(levelI))
                        {
                            if (cLabel < curLeaf)
                            {
                                create = false;
                                break;
                            }

                            centreNodes.append(cubeLabel[cLabel]);

                            // adding face centre labels
                            if (faceCentreLabel[cLabel].size() != 0)
                            {
                                const label helpFace = eI/4;

                                const label fcl =
                                    faceCentreLabel
                                    [cLabel]
                                    [faceCentreHelper_[helpFace][i]];

                                if (fcl != -1)
                                    centreNodes.append(fcl);
                            }
                        }
                        else if (oc.level() < direction(levelI))
                        {
                            const meshOctreeCubeCoordinates sc
                            (
                                cc + vlPos[startNode][fNodes[i]]
                            );

                            label pos(-1);

                            for (label j = 0; j < 8; j++)
                            {
                                if (sc == oc.refineForPosition(j))
                                {
                                    pos = j;
                                    break;
                                }
                            }

                            if (pos == -1)
                                FatalErrorInFunction
                                    << "Cannot find cube position"
                                    << abort(FatalError);

                            centreNodes.append(subNodeLabels[cLabel][pos]);
                        }
                    }

                    // create tets around this edge
                    if (create)
                    {
                        const label nCentres = centreNodes.size();

                        forAll(centreNodes, i)
                        {
                            if (centreNodes[i] == -1)
                                continue;
                            if (centreNodes[(i + 1)%nCentres] == -1)
                                continue;

                            partTet tet
                            (
                                centreNodes[i],
                                centreNodes[(i + 1)%nCentres],
                                start,
                                end
                            );

                            tets_.append(tet);
                        }
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
