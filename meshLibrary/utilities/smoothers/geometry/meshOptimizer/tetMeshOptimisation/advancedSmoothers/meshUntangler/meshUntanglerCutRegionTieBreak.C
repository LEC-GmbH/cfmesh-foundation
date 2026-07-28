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
#include "sortEdgesIntoChains.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshUntangler::cutRegion::tieBreak
(
    const DynamicList<label>& f
)
{
    // There must be only one, singly connected region which is cut off
    // from the feasible region. This may not be the case if one operates in the
    // range of computer tolerances. In order to resolve the tie it is necessary
    // to find a single region which will be cut off from the feasible region.
    // This will be done by finding the node with the most negative distance,
    // and then start marking vertices connected to that vertex via an edge.

    // delete pointer data
    deleteDemandDrivenData(cPtsPtr_);
    deleteDemandDrivenData(cEdgesPtr_);
    deleteDemandDrivenData(cFacesPtr_);

    const DynamicList<edge>& edges = *edgesPtr_;

    DynamicList<edge> faceEdges;
    forAll(f, eI)
    {
        faceEdges.append(edges[f[eI]]);
    }

    const DynamicList<DynamicList<label>> fvertices =
        sortEdgesIntoChains(faceEdges).sortedChains();

    if (fvertices.size() != 1)
    {
        valid_ = false;
        return;
    }

    const DynamicList<label>& fv = fvertices[0];

    DynamicList<label> vertexRegion;
    vertexRegion.setSize(fv.size());
    vertexRegion = 0;

    label region(1);
    forAll(fv, vI)
    {
        if (!vertexTypes_[fv[vI]] && !vertexRegion[vI])
        {
            vertexRegion[vI] = region;

            label fcI = fv.fcIndex(vI);
            label rcI = fv.rcIndex(vI);
            bool found;
            do
            {
                found = false;
                if (!vertexTypes_[fv[fcI]])
                {
                    vertexRegion[fcI] = region;
                    fcI = fv.fcIndex(fcI);
                    found = true;
                }

                if (!vertexTypes_[fv[rcI]])
                {
                    vertexRegion[rcI] = region;
                    rcI = fv.rcIndex(rcI);
                    found = true;
                }
            } while (found);

            ++region;
        }
    }

    if (region > 2)
    {
        // there are more than two regions which need to be cut off
        // there should be only one cut-off region
        // there this region will be determined by the most negative
        // distance from plane
        scalar minDist(VGREAT);
        label minRegion(-1);
        forAll(fv, vI)
        {
            if (vertexRegion[vI] && (vertexDistance_[fv[vI]] < minDist))
            {
                minDist = vertexDistance_[fv[vI]];
                minRegion = vertexRegion[vI];
            }
        }

        forAll(vertexRegion, vI)
        {
            if (vertexRegion[vI] && (vertexRegion[vI] != minRegion))
            {
                vertexTypes_[fv[vI]] |= INPLANE;
            }
        }
    }
    else
    {
        forAll(fv, vI)
        {
            if (
                (vertexTypes_[fv[vI]] & INPLANE) &&
                !vertexRegion[fv.rcIndex(vI)] &&
                !vertexRegion[fv.fcIndex(vI)]
            )
            {
                vertexTypes_[fv[vI]] ^= INPLANE;
        }
                vertexTypes_[fv[vI]] |= KEEP;
            }
    }

    // create new points
    const DynamicList<point>& points = *pointsPtr_;
    cPtsPtr_ = new DynamicList<point>();
    newVertexLabel_ = -1;
    origNumVertices_ = 0;
    forAll(points, pI)
    {
        if (vertexTypes_[pI])
        {
            cPtsPtr_->append(points[pI]);
            newVertexLabel_[pI] = origNumVertices_++;
        }
    }

    // find new edges and continue creating faces
    findNewEdges();
}


// ************************************************************************* //
