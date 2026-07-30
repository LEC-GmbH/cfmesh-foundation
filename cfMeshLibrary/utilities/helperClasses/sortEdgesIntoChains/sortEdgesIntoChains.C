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

#include "sortEdgesIntoChains.H"
#include "helperFunctions.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::sortEdgesIntoChains::createNodeLabels()
{
    label nPoints(0);
    forAll(bEdges_, eI)
    {
        const edge& e = bEdges_[eI];
        if (!newNodeLabel_.found(e.start()))
        {
            newNodeLabel_.insert(e.start(), nPoints++);
        }
        if (!newNodeLabel_.found(e.end()))
        {
            newNodeLabel_.insert(e.end(), nPoints++);
        }
    }

    edgesAtPoint_.setSize(nPoints);
    forAll(bEdges_, eI)
    {
        const edge& e = bEdges_[eI];
        label l = newNodeLabel_[e.start()];
        edgesAtPoint_[l].append(eI);

        l = newNodeLabel_[e.end()];
        edgesAtPoint_[l].append(eI);
    }

    forAll(edgesAtPoint_, pI)
    {
        if (edgesAtPoint_[pI].size() % 2)
        {
            openEdges_ = true;
        }
    }
}


bool Foam::Module::sortEdgesIntoChains::findPointsBelongingToTheChain
(
    const label currPos,
    DynamicList<bool>& chainEdges
) const
{
    chainEdges.setSize(bEdges_.size());
    chainEdges = false;

    if (edgesAtPoint_[currPos].size() != 2)
    {
        return false;
    }

    const label commonVrt =
        bEdges_[edgesAtPoint_[currPos][0]].commonVertex
            (
                bEdges_[edgesAtPoint_[currPos][1]]
            );
    label prevVrt = bEdges_[edgesAtPoint_[currPos][0]].otherVertex(commonVrt);
    label nextVrt = bEdges_[edgesAtPoint_[currPos][1]].otherVertex(commonVrt);
    forAll(edgesAtPoint_[currPos], posI)
    {
        chainEdges[edgesAtPoint_[currPos][posI]] = true;
    }

    bool found;
    do
    {
        found = false;

        const DynamicList<label>& vEdges = edgesAtPoint_[newNodeLabel_[prevVrt]];
        if (vEdges.size() == 2)
        {
            forAll(vEdges, eI)
            {
                if (!chainEdges[vEdges[eI]])
                {
                    found = true;
                    chainEdges[vEdges[eI]] = true;
                    prevVrt = bEdges_[vEdges[eI]].otherVertex(prevVrt);
                }
            }
        }
    } while (found);

    do
    {
        found = false;

        const DynamicList<label>& vEdges = edgesAtPoint_[newNodeLabel_[nextVrt]];
        if (vEdges.size() == 2)
        {
            forAll(vEdges, eI)
            {
                if (!chainEdges[vEdges[eI]])
                {
                    found = true;
                    chainEdges[vEdges[eI]] = true;
                    nextVrt = bEdges_[vEdges[eI]].otherVertex(nextVrt);
                }
            }
        }
    } while (found);

    if
    (
        (edgesAtPoint_[newNodeLabel_[nextVrt]].size() != 2) &&
        (edgesAtPoint_[newNodeLabel_[prevVrt]].size() != 2) &&
        (prevVrt != nextVrt)
    )
    {
        chainEdges = false;
        return false;
    }

    return true;
}


void Foam::Module::sortEdgesIntoChains::shrinkEdges
(
    const DynamicList<bool>& chainEdges
)
{
    forAll(chainEdges, eI)
    {
        if (chainEdges[eI])
        {
            const edge& e = bEdges_[eI];

            DynamicList<label>& atStart =
                edgesAtPoint_[newNodeLabel_[e.start()]];
            atStart.erase(atStart.begin() + findIndex(atStart, eI));

            DynamicList<label>& atEnd =
                edgesAtPoint_[newNodeLabel_[e.end()]];
            atEnd.erase(atEnd.begin() + findIndex(atEnd, eI));
        }
    }
}


void Foam::Module::sortEdgesIntoChains::createChainFromEdges
(
    const DynamicList<bool>& chainEdges
)
{
    label i(0);
    forAll(chainEdges, eI)
    {
        if (chainEdges[eI])
        {
            ++i;
        }
    }


    DynamicList<label> chainPoints;
    chainPoints.setSize(i);
    i = 0;

    forAll(chainEdges, eI)
    {
        if (chainEdges[eI])
        {
            chainPoints[i++] = bEdges_[eI].start();
            chainPoints[i++] = bEdges_[eI].end();

            bool found;
            do
            {
                found = false;
                const DynamicList<label>& pEdges =
                    edgesAtPoint_[newNodeLabel_[chainPoints[i - 1]]];

                forAll(pEdges, peI)
                {
                    if (chainEdges[pEdges[peI]])
                    {
                        const label otherPoint =
                            bEdges_[pEdges[peI]].otherVertex(chainPoints[i - 1]);

                        if (otherPoint == -1)
                        {
                            continue;
                        }
                        if (chainPoints[i - 2] == otherPoint)
                        {
                            continue;
                        }
                        if (chainPoints[0] == otherPoint)
                        {
                            continue;
                        }

                        found = true;
                        chainPoints[i++] = otherPoint;
                    }
                }
            } while (found);

            createdChains_.append(chainPoints);

            break;
        }
    }
}


void Foam::Module::sortEdgesIntoChains::sortEdges()
{
    createNodeLabels();

    if (!openEdges_)
    {
        DynamicList<bool> chainEdges;
        chainEdges.reserve(bEdges_.size());
        forAll(edgesAtPoint_, pI)
        {
            if (findPointsBelongingToTheChain(pI, chainEdges))
            {
                createChainFromEdges(chainEdges);

                shrinkEdges(chainEdges);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Module::sortEdgesIntoChains::sortEdgesIntoChains
(
    const DynamicList<edge>& bEdges
)
:
    bEdges_(bEdges),
    openEdges_(false)
{
    sortEdges();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::DynamicList<Foam::DynamicList<Foam::label>>&
Foam::Module::sortEdgesIntoChains::sortedChains() const
{
    return createdChains_;
}


// ************************************************************************* //
