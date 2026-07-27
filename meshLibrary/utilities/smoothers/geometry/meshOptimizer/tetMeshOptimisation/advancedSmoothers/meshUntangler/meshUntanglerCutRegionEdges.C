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

#include "demandDrivenData.H"
#include "meshUntangler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshUntangler::cutRegion::findNewEdges()
{
    cEdgesPtr_ = new DynamicList<edge>();
    DynamicList<edge>& cEdges = *cEdgesPtr_;

    const DynamicList<point>& pts = *pointsPtr_;
    const DynamicList<edge>& edges = *edgesPtr_;

    newEdgeLabel_.setSize(edges.size());
    newEdgeLabel_ = -1;

    forAll(edges, eI)
    {
        const edge& e = edges[eI];
        const label start = e.start();
        const label end = e.end();
        const label newStartLabel = newVertexLabel_[start];
        const label newEndLabel = newVertexLabel_[end];

        if ((newStartLabel != -1) && (newEndLabel != -1))
        {
            newEdgeLabel_[eI] = cEdges.size();
            cEdges.append(edge(newStartLabel, newEndLabel));
        }
        else if ((newEndLabel != -1) && (vertexTypes_[end] & KEEP))
        {
            // start edge vertex is not visible, but the other one is
            newEdgeLabel_[eI] = cEdges.size();
            cEdges.append(edge(newEndLabel, cPtsPtr_->size()));

            const scalar t =
                -vertexDistance_[start] /
                (vertexDistance_[end] - vertexDistance_[start]);

            const point newP = (1.0 - t)*pts[start] + t*pts[end];
            cPtsPtr_->append(newP);
        }
        else if ((newStartLabel != -1) && (vertexTypes_[start] & KEEP))
        {
            // end edge vertex is not visible, but the other one is
            newEdgeLabel_[eI] = cEdges.size();
            cEdges.append(edge(newStartLabel, cPtsPtr_->size()));

            const scalar t =
                -vertexDistance_[end] /
                (vertexDistance_[start] - vertexDistance_[end]);

            const point newP = (1.0 - t)*pts[end] + t*pts[start];
            cPtsPtr_->append(newP);
        }
    }
}


// ************************************************************************* //
