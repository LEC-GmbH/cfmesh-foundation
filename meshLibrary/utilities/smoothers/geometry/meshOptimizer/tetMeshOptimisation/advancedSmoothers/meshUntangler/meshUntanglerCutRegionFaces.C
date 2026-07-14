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

void Foam::Module::meshUntangler::cutRegion::findNewFaces()
{

    bool finished;
    do
    {
        finished = true;

        const DynList<DynList<label, 8>, 64>& fcs = *facesPtr_;
        DynList<edge, 128>& cEdges = *cEdgesPtr_;

        cFacesPtr_ = new DynList<DynList<label, 8>, 64>();
        DynList<DynList<label, 8>, 64>& cFaces = *cFacesPtr_;

        DynList<label, 8> faceInPlane;

        DynList<label, 64> pointUsage;
        pointUsage.setSize(cPtsPtr_->size());

        forAll(fcs, fI)
        {
            const DynList<label, 8>& f = fcs[fI];

            pointUsage = 0;

            DynList<label, 8> newFace;

            forAll(f, eI)
            {

                const label edgeLabel = newEdgeLabel_[f[eI]];

                if (edgeLabel != -1)
                {

                    const edge& e = cEdges[edgeLabel];
                    ++pointUsage[e[0]];
                    ++pointUsage[e[1]];
                    newFace.append(edgeLabel);
                }
            }

            if (newFace.size() > 1)
            {
                DynList<label, 4> newEdge;
                forAll(pointUsage, pI)
                    if (pointUsage[pI] == 1)
                        newEdge.append(pI);

                if (newEdge.size() == 2)
                {

                    newFace.append(cEdges.size());
                    cEdges.append(edge(newEdge[0], newEdge[1]));
                }
                else if (newEdge.size() > 2)
                {

                    tieBreak(f);
                    if (!valid_) return;
                    finished = false;
                    break;
                }

                cFaces.append(newFace);
            }
        }

        if (!finished) continue;

        // find edges which form the faceInPlane
        DynList<label, 128> edgeUsage;
        edgeUsage.setSize(cEdges.size());
        edgeUsage = 0;
        forAll(cFaces, fI)
        {
            const DynList<label, 8>& f = cFaces[fI];

            forAll(f, eI)
                ++edgeUsage[f[eI]];
        }

        forAll(edgeUsage, eI)
            if (edgeUsage[eI] == 1)
                faceInPlane.append(eI);

        if (faceInPlane.size() > 2)
        {

            cFaces.append(faceInPlane);
        }

    } while (!finished);
}


// ************************************************************************* //
