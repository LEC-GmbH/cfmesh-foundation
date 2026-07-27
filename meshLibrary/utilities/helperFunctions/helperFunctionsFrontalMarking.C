/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | www.cfmesh.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 Creative Fields, Ltd.
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "error.H"
#include "helperFunctionsFrontalMarking.H"
#include "DynamicList.H"
#include "labelPair.H"
#include "HashSet.H"
#include "labelListListOps.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{
namespace Module
{
namespace help
{

class graphNeiOp
{
    // Private data

        // const reference to labelListList2D
        const labelListList2D& neiGroups_;

public:

    // Constructors

        // Construct from labelListList2D
        inline graphNeiOp(const labelListList2D& neiGroups)
        :
            neiGroups_(neiGroups)
        {}


    // Public member functions

        //- The size of the graph
        inline label size() const
        {
            return neiGroups_.size();
        }

        // operator for finding neighbours
        inline void operator()(const label groupI, DynamicList<label>& ng) const
        {
            ng.setSize(neiGroups_[groupI].size());
            forAllRow(neiGroups_, groupI, i)
                ng[i] = neiGroups_[groupI][i];
        }
};


class graphSelectorOp
{
    // Private data

        // const reference to labelListList2D
        const labelListList2D& neiGroups_;


public:

    // Constructors

        // Construct from labelListList2D
        inline graphSelectorOp(const labelListList2D& neiGroups)
        :
            neiGroups_(neiGroups)
        {}


    // Public member functions

        // operator for selecting elements
        inline bool operator()(const label groupI) const
        {
            if ((groupI < 0) || (groupI >= neiGroups_.size()))
                return false;

            return true;
        }
};

} // End namespace help
} // End namespace Module
} // End namespace Foam


template<class labelListType, class neiOp, class filterOp>
void Foam::Module::help::frontalMarking
(
    labelListType& result,
    const label startingIndex,
    const neiOp& neighbourCalculator,
    const filterOp& selector
)
{
    // add the starting element
    result.clear();
    result.append(startingIndex);

    // add the starting element to the front
    DynamicList<label> front;
    front.append(startingIndex);

    // store information which element were already visited
    boolList alreadySelected(neighbourCalculator.size(), false);

    // start with frontal marking
    while (front.size())
    {
        const label eLabel = front.remove();

        // find neighbours of the current element
        DynamicList<label> neighbours;
        neighbourCalculator(eLabel, neighbours);

        forAll(neighbours, neiI)
        {
            const label nei = neighbours[neiI];

            if (nei < 0)
            {
                continue;
            }
            if (alreadySelected[nei])
            {
                continue;
            }

            if (selector(nei))
            {
                alreadySelected[nei] = true;
                front.append(nei);
                result.append(nei);
            }
        }
    }
}


template<class labelListType, class neiOp, class filterOp>
Foam::label Foam::Module::help::groupMarking
(
    labelListType& elementInGroup,
    const neiOp& neighbourCalculator,
    const filterOp& selector
)
{
    label nGroups(0);

    elementInGroup.setSize(neighbourCalculator.size());
    elementInGroup = -1;

    labelListList2D neighbouringGroups;

    label nThreads(1);

    //nThreads = 3*omp_get_num_procs();

    DynamicList<label> nGroupsAtThread(nThreads, 0);

    # pragma omp parallel num_threads(nThreads)
    {
        const label chunkSize =
            Foam::max(1, neighbourCalculator.size() / nThreads);

        const label threadI = omp_get_thread_num();

        DynamicList<labelPair> threadCommPairs;

        const label minEl = threadI*chunkSize;

        label maxEl = minEl + chunkSize;
        if (threadI == (nThreads - 1))
        {
            maxEl = Foam::max(maxEl, neighbourCalculator.size());
        }

        label& groupI = nGroupsAtThread[threadI];
        groupI = 0;

        for (label elI = minEl; elI < maxEl; ++elI)
        {
            if (elementInGroup[elI] != -1)
            {
                continue;
            }
            if (!selector(elI))
            {
                continue;
            }

            elementInGroup[elI] = groupI;
            DynamicList<label> front;
            front.append(elI);

            while (front.size())
            {
                const label eLabel = front.remove();

                DynamicList<label> neighbours;
                neighbourCalculator(eLabel, neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if ((nei < 0) || (elementInGroup[nei] != -1))
                    {
                        continue;
                    }

                    if ((nei < minEl) || (nei >= maxEl))
                    {
                        // this is a communication interface between
                        // two threads
                        threadCommPairs.append(labelPair(elI, nei));
                    }
                    else if (selector(nei))
                    {
                        // this element is part of the same group
                        elementInGroup[nei] = groupI;
                        front.append(nei);
                    }
                }
            }

            ++groupI;
        }

        # pragma omp barrier

        # pragma omp master
        {
            forAll(nGroupsAtThread, i)
            {
                nGroups += nGroupsAtThread[i];
            }
        }

        label startGroup(0);
        for (label i = 0; i < threadI; ++i)
        {
            startGroup += nGroupsAtThread[i];
        }

        for (label elI = minEl; elI < maxEl; ++elI)
        {
            elementInGroup[elI] += startGroup;
        }

        # pragma omp barrier

        // find group to neighbouring groups addressing
        labelListList2D localNeiGroups(nGroups);
        forAll(threadCommPairs, cfI)
        {
            const labelPair& lp = threadCommPairs[cfI];
            const label groupI = elementInGroup[lp.first()];
            const label neiGroup = elementInGroup[lp.second()];

            if ((neiGroup >= nGroups) || (groupI >= nGroups))
            {
                FatalError << "neiGroup " << neiGroup
                    << " groupI " << groupI << " are >= than "
                    << "nGroups " << nGroups << abort(FatalError);
            }

            if (neiGroup != -1)
            {
                appendUniq(localNeiGroups[groupI], neiGroup);
                appendUniq(localNeiGroups[neiGroup], groupI);
            }
        }

        # pragma omp critical
        {
            neighbouringGroups.setSize(nGroups);

            forAll(localNeiGroups, groupI)
            {
                const DynamicList<label>& lGroups = localNeiGroups[groupI];

                appendUniq(neighbouringGroups[groupI], groupI);

                forAll(lGroups, i)
                {
                    appendUniq(neighbouringGroups[groupI], lGroups[i]);
                }
            }
        }
    }

    forAll(neighbouringGroups, i)
    {
        labelList helper(neighbouringGroups[i].size());
        forAllRow(neighbouringGroups, i, j)
        {
            helper[j] = neighbouringGroups[i][j];
        }

        sort(helper);

        setRow(neighbouringGroups, i, helper);
    }

    // start processing connections between the group and merge the connected
    // ones into a new group
    DynamicList<label> globalGroupLabel;
    globalGroupLabel.setSize(nGroups);
    globalGroupLabel = -1;

    // reduce the information about the groups
    label counter(0);

    forAll(neighbouringGroups, groupI)
    {
        if (globalGroupLabel[groupI] != -1)
        {
            continue;
        }

        DynamicList<label> connectedGroups;
        frontalMarking
        (
            connectedGroups,
            groupI,
            graphNeiOp(neighbouringGroups),
            graphSelectorOp(neighbouringGroups)
        );

        forAll(connectedGroups, gI)
        {
            globalGroupLabel[connectedGroups[gI]] = counter;
        }

        ++counter;
    }

    nGroups = counter;

    forAll(neighbouringGroups, groupI)
    {
        if (globalGroupLabel[groupI] != -1)
        {
            continue;
        }

        forAllRow(neighbouringGroups, groupI, ngI)
        {
            globalGroupLabel[neighbouringGroups[groupI][ngI]] = counter;
        }

        ++counter;
    }

    // set the global group label
    # pragma omp parallel for schedule(dynamic, 50)
    forAll(elementInGroup, elI)
    {
        if (elementInGroup[elI] < 0)
        {
            continue;
        }

        elementInGroup[elI] = globalGroupLabel[elementInGroup[elI]];
    }

    return nGroups;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
