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

#include "voronoiMeshExtractor.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::voronoiMeshExtractor::createAddressing() const
{
    if (pointEdgesPtr_ || edgeTetsPtr_ || boundaryEdgePtr_ || edgesPtr_)
    {
        return;
    }

    pointEdgesPtr_ = new labelListList2D(tetCreator_.tetPoints().size());
    labelListList2D& pointEdges = *pointEdgesPtr_;

    edgeTetsPtr_ = new DynamicList<DynamicList<label>>();
    DynamicList<DynamicList<label>>& edgeTets = *edgeTetsPtr_;

    boundaryEdgePtr_ = new boolList();
    boolList& boundaryEdges = *boundaryEdgePtr_;

    edgesPtr_ = new DynamicList<edge>();
    DynamicList<edge>& edges = *edgesPtr_;

    // create edges and edgeTets
    const DynamicList<partTet>& tets = tetCreator_.tets();

    labelListList2D pointTets;
    reverseAddressing(pointTets, tets);

    forAll(tets, tetI)
    {
        const FixedList<edge, 6> tetEdges = tets[tetI].edges();

        forAll(tetEdges, eI)
        {
            const edge& e = tetEdges[eI];
            const label start = e.start();

            const DynamicList<label>& endTets = pointTets[e.end()];

            bool store(true);

            DynamicList<label> eTets;
            forAllRow(pointTets, start, ptI)
            {
                const label tetJ = pointTets[start][ptI];

                if (findIndex(endTets, tetJ) < 0)
                {
                    continue;
                }

                if (tetJ < tetI)
                {
                    store = false;
                    break;
                }

                eTets.append(tetJ);
            }

            if (store)
            {
                appendList(edgeTets, eTets);

                edges.append(e);
            }
        }
    }

    // calculate point-edges addressing
    reverseAddressing(pointEdges, edges);

    // sort edge-tets in circular order
    boundaryEdges.setSize(edgeTets.size());
    boundaryEdges = false;

    # pragma omp parallel for schedule(dynamic, 100)
    forAll(edgeTets, edgeI)
    {
        const edge& e = edges[edgeI];
        DynamicList<label>& eTets = edgeTets[edgeI];

        forAll(eTets, tetI)
        {
            const partTet& pt = tets[eTets[tetI]];

            // find the face shared with the neighbour
            const FixedList<edge, 6> tetEdges = pt.edges();

            label searchPoint(-1);
            forAll(tetEdges, eI)
            {
                if (tetEdges[eI] == e)
                {
                    if (tetEdges[eI].start() == e.start())
                    {
                        searchPoint = pt[sameOrientation_[eI]];
                    }
                    else
                    {
                        searchPoint = pt[oppositeOrientation_[eI]];
                    }

                    break;
                }
            }

            if (searchPoint < 0)
            {
                FatalErrorInFunction
                    << " invalid search point " << abort(FatalError);
            }

            bool found(false);
            forAll(eTets, i)
            {
                if (tetI == i)
                {
                    continue;
                }

                const partTet& ptNei = tets[eTets[i]];

                const label pos = ptNei.whichPosition(searchPoint);

                if (pos < 0)
                {
                    continue;
                }

                if (tetI < eTets.size()-1)
                {
                    const label add = eTets[tetI + 1];
                    eTets[tetI + 1] = eTets[i];
                    eTets[i] = add;
                }
                found = true;
                break;
            }

            if (!found)
            {
                boundaryEdges[edgeI] = true;
            }
        }
    }
}


const Foam::Module::labelListList2D&
Foam::Module::voronoiMeshExtractor::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        createAddressing();
    }

    return *pointEdgesPtr_;
}


const Foam::DynamicList<Foam::edge>&
Foam::Module::voronoiMeshExtractor::edges() const
{
    if (!edgesPtr_)
    {
        createAddressing();
    }

    return *edgesPtr_;
}


const Foam::DynamicList<Foam::DynamicList<Foam::label>>&
Foam::Module::voronoiMeshExtractor::edgeTets() const
{
    if (!edgeTetsPtr_)
    {
        createAddressing();
    }

    return *edgeTetsPtr_;
}


const Foam::boolList& Foam::Module::voronoiMeshExtractor::boundaryEdge() const
{
    if (!boundaryEdgePtr_)
    {
        createAddressing();
    }

    return *boundaryEdgePtr_;
}


// ************************************************************************* //
