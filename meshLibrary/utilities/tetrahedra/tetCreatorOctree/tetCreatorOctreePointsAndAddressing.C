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
#include "meshOctree.H"
#include "triSurface.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::tetCreatorOctree::selectElements()
{
    const List<direction>& boxType = octreeCheck_.boxType();
    const meshOctree& octree = octreeCheck_.octree();
    const boundBox& rootBox = octree.rootBox();

    // store nodeLabels first
    const labelListList2D& nodeLabels = octreeCheck_.nodeLabels();
    tetPoints_.setSize(octreeCheck_.numberOfNodes());

    forAll(nodeLabels, leafI)
    {
        if (nodeLabels[leafI].size() != 0)
        {
            const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);
            FixedList<point, 8> lv;
            oc.vertices(rootBox, lv);

            forAll(lv, vI)
            {
                tetPoints_[nodeLabels[leafI][vI]] = lv[vI];
            }
        }
    }

    // create cubeLabel list
    if (!cubeLabelPtr_)
    {
        cubeLabelPtr_ = new labelList();
    }
    labelList& cubeLabel = *cubeLabelPtr_;
    cubeLabel.setSize(octree.numberOfLeaves());
    cubeLabel = -1;

    forAll(boxType, leafI)
    {
        if (boxType[leafI] & meshOctreeAddressing::MESHCELL)
        {
            const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);
            cubeLabel[leafI] = tetPoints_.size();
            tetPoints_.append(oc.centre(rootBox));
        }
    }
}


void Foam::Module::tetCreatorOctree::createPointsAndAddressing()
{
    selectElements();

    const meshOctree& octree = octreeCheck_.octree();
    const boundBox& rootBox = octree.rootBox();
    const labelListList2D& nodeLabels = octreeCheck_.nodeLabels();

    if (!subNodeLabelsPtr_)
    {
        subNodeLabelsPtr_ = new labelListList2D(octree.numberOfLeaves());
    }
    labelListList2D& subNodeLabels = *subNodeLabelsPtr_;

    // store nodeLabels
    direction maxLevel(0);
    forAll(nodeLabels, leafI)
    {
        if (octree.returnLeaf(leafI).level() > maxLevel)
        {
            maxLevel = octree.returnLeaf(leafI).level();
        }
    }

    sortedLeaves_.setSize(maxLevel + 1);
    forAll(sortedLeaves_, levelI)
    {
        sortedLeaves_[levelI].clear();
    }

    forAll(nodeLabels, leafI)
    {
        const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);
        sortedLeaves_[oc.level()].append(leafI);
    }

    // create subNodeLabels
    const List<FixedList<label, 8>>& pointLeaves = octreeCheck_.nodeLeaves();

    forAll(pointLeaves, pointI)
    {
        bool validLeaf[8];
        direction levelI(0);
        for (label i = 0; i < 8; ++i)
        {
            const label pointLeafI = pointLeaves[pointI][i];

            if (pointLeafI == -1)
            {
                validLeaf[i] = false;
            }
            else
            {
                validLeaf[i] = true;
                for (label j = i + 1; j < 8; ++j)
                {
                    if (pointLeafI == pointLeaves[pointI][j])
                    {
                        validLeaf[i] = false;
                        validLeaf[j] = false;
                    }
                }

                const direction level =
                    octree.returnLeaf(pointLeafI).level();

                if (level > levelI)
                {
                    levelI = level;
                }
            }
        }

        for (label plI = 0; plI < 8; ++plI)
        {
            if (validLeaf[plI])
            {
                const label pointLeafI = pointLeaves[pointI][plI];

                const meshOctreeCubeBasic& lc =
                    octree.returnLeaf(pointLeafI);

                if (lc.level() < levelI)
                {
                    if (subNodeLabels[pointLeafI].size() != 8)
                    {
                        subNodeLabels[pointLeafI].setSize(8);
                        forAllRow(subNodeLabels, pointLeafI, k)
                        {
                            subNodeLabels[pointLeafI][k] = -1;
                        }
                    }

                    subNodeLabels[pointLeafI][(7 - plI)] = tetPoints_.size();
                    FixedList<point, 8> lv;
                    lc.vertices(rootBox, lv);

                    tetPoints_.append
                    (
                        0.5*(lv[7 - plI] + lc.centre(rootBox))
                    );
                }
            }
        }
    }

    createFaceCentreLabels();
}


void Foam::Module::tetCreatorOctree::createFaceCentreLabels()
{
    Info<< "Creating face centre labels " << endl;
    const labelList& cubeLabel = *cubeLabelPtr_;
    const labelListList2D& nodeLabels = octreeCheck_.nodeLabels();
    const List<FixedList<label, 8>>& pointLeaves = octreeCheck_.nodeLeaves();
    const meshOctree& octree = octreeCheck_.octree();

    List<direction> nodeLevel(pointLeaves.size(), direction(0));
    forAll(nodeLabels, leafI)
    {
        const direction level = octree.returnLeaf(leafI).level();

        forAllRow(nodeLabels, leafI, nlI)
        {
            const label nLabel = nodeLabels[leafI][nlI];

            if (nodeLevel[nLabel] < level)
            {
                nodeLevel[nLabel] = level;
            }
        }
    }

    if (!faceCentreLabelPtr_)
    {
        faceCentreLabelPtr_ = new labelListList2D(cubeLabel.size());
    }
    labelListList2D& faceCentreLabel = *faceCentreLabelPtr_;

    forAll(cubeLabel, cubeI)
    {
        if (cubeLabel[cubeI] != -1)
        {
            const direction level = octree.returnLeaf(cubeI).level();

            for (label i = 0; i < 6; ++i)
            {
                if (
                    (faceCentreLabel[cubeI].size() != 0) &&
                    (faceCentreLabel[cubeI][i] != -1)
                )
                    continue;

                FixedList<label, 4> faceNodes;
                forAll(faceNodes, fnI)
                {
                    faceNodes[fnI] =
                        meshOctreeCubeCoordinates::faceNodes_[i][fnI];
                }

                label highLevelNode(-1);
                for (label j = 0; j < 4; ++j)
                {
                    if (nodeLevel[nodeLabels[cubeI][faceNodes[j]]] > level)
                    {
                        highLevelNode = j;
                        break;
                    }
                }

                if (highLevelNode == -1)
                {
                    continue;
                }

                DynamicList<label> neighbours;
                octree.findNeighboursInDirection(cubeI, i, neighbours);

                if ((neighbours.size() != 1) || (neighbours[0] == -1))
                {
                    continue;
                }

                if (faceCentreLabel[cubeI].size() == 0)
                {
                    faceCentreLabel[cubeI].setSize(6);
                    forAllRow(faceCentreLabel, cubeI, colI)
                    {
                        faceCentreLabel[cubeI][colI] = -1;
                    }
                }
                const label cNei = neighbours[0];
                if (faceCentreLabel[cNei].size() == 0)
                {
                    faceCentreLabel[cNei].setSize(6);
                    forAllRow(faceCentreLabel, cNei, colI)
                    {
                        faceCentreLabel[cNei][colI] = -1;
                    }
                }

                faceCentreLabel[cubeI][i] = tetPoints_.size();
                switch(i)
                {
                    case 0:
                    {
                        faceCentreLabel[cNei][1] = tetPoints_.size();
                    } break;
                    case 1:
                    {
                        faceCentreLabel[cNei][0] = tetPoints_.size();
                    } break;
                    case 2:
                    {
                        faceCentreLabel[cNei][3] = tetPoints_.size();
                    } break;
                    case 3:
                    {
                        faceCentreLabel[cNei][2] = tetPoints_.size();
                    } break;
                    case 4:
                    {
                        faceCentreLabel[cNei][5] = tetPoints_.size();
                    } break;
                    case 5:
                    {
                        faceCentreLabel[cNei][4] = tetPoints_.size();
                    } break;
                };

                point p(vector::zero);
                for (label j = 0; j < 4; ++j)
                {
                    p += tetPoints_[nodeLabels[cubeI][faceNodes[j]]];
                }
                p /= 4;
                tetPoints_.append(p);
            }
        }
    }
}


// ************************************************************************* //
