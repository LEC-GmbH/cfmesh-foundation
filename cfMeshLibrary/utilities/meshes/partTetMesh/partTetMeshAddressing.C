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
#include "polyMeshGenModifier.H"
#include "partTetMesh.H"
#include "tetMatcher.H"
#include "polyMeshGenAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::partTetMesh::createPointsAndTets
(
    const List<direction>& useCell,
    const boolList& lockedPoints
)
{
    const pointFieldPMG& points = origMesh_.points();
    const faceListPMG& faces = origMesh_.faces();
    const cellListPMG& cells = origMesh_.cells();
    const labelList& owner = origMesh_.owner();
    const labelList& neighbour = origMesh_.neighbour();
    const label nInternalFaces = origMesh_.nInternalFaces();

    // check how many neighbours of a face are marked for smoothing
    labelList usedFace(faces.size(), 0);

    // mark faces
    forAll(faces, faceI)
    {
        if (useCell[owner[faceI]])
        {
            ++usedFace[faceI];
        }

        if (neighbour[faceI] < 0)
        {
            continue;
        }

        if (useCell[neighbour[faceI]])
        {
            ++usedFace[faceI];
        }
    }

    const vectorField& faceCentres = origMesh_.addressingData().faceCentres();
    const vectorField& cellCentres = origMesh_.addressingData().cellCentres();

    DynamicList<label> nodeLabelForPoint(points.size(), -1);
    DynamicList<label> nodeLabelForFace(faces.size(), -1);
    DynamicList<label> nodeLabelForCell(cells.size(), -1);

    points_.clear();
    smoothVertex_.clear();

    // create BOUNDARY points
    // Note: iterate the boundary face range directly rather than via
    // boundaries[patchI].patchStart()/patchSize() - named patches are not
    // guaranteed to tile [nInternalFaces, faces.size()) exactly at every
    // point in the mesh-optimisation pipeline (e.g. mid-untangling, after
    // faces/cells have been added), and a face left outside all patches
    // would silently keep its nodeLabelForFace/nodeLabelForPoint at -1,
    // corrupting later-built tets.
    for (label faceI = nInternalFaces; faceI < faces.size(); ++faceI)
    {
        if (!usedFace[faceI])
        {
            continue;
        }

        const face& f = faces[faceI];

        if (f.size() > 3)
        {
            // create face centre
            nodeLabelForFace[faceI] = points_.size();
            points_.append(faceCentres[faceI]);
            smoothVertex_.append(FACECENTRE);
        }

        // add face corners
        forAll(f, pI)
        {
            const label pointI = f[pI];
            if (nodeLabelForPoint[pointI] == -1)
            {
                nodeLabelForPoint[pointI] = points_.size();
                points_.append(points[pointI]);

                smoothVertex_.append(BOUNDARY);
            }
        }
    }

    // create points for internal faces
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        if (usedFace[faceI])
        {
            const face& f = faces[faceI];

            if (f.size() > 3)
            {
                // create face centre
                nodeLabelForFace[faceI] = points_.size();
                points_.append(faceCentres[faceI]);
                smoothVertex_.append(FACECENTRE);
            }

            // add face corners
            const direction vType = usedFace[faceI]==2?SMOOTH:NONE;
            forAll(f, pI)
            {
                const label pointI = f[pI];
                if (nodeLabelForPoint[pointI] == -1)
                {
                    nodeLabelForPoint[pointI] = points_.size();
                    points_.append(points[pointI]);

                    smoothVertex_.append(vType);
                }
                else if (vType == NONE)
                {
                     smoothVertex_[nodeLabelForPoint[pointI]] = NONE;
                }
            }
        }
    }

    // create tets
    tetMatcher tet;
    forAll(useCell, cI)
    {
        if (useCell[cI])
        {
            const cell& c = cells[cI];

            if (tet.matchShape(false, faces, owner, cI, cells[cI]))
            {
                const labelList& tVrt = tet.vertLabels();

                // add tet
                tets_.append
                (
                    partTet
                    (
                        nodeLabelForPoint[tVrt[0]],
                        nodeLabelForPoint[tVrt[1]],
                        nodeLabelForPoint[tVrt[2]],
                        nodeLabelForPoint[tVrt[3]]
                    )
                );

                continue;
            }

            nodeLabelForCell[cI] = points_.size();
            const label centreLabel = points_.size();
            points_.append(cellCentres[cI]);
            smoothVertex_.append(CELLCENTRE);

            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                if (owner[c[fI]] == cI)
                {
                    if (f.size() == 3)
                    {
                        partTet tet
                        (
                            nodeLabelForPoint[f[0]],
                            nodeLabelForPoint[f[2]],
                            nodeLabelForPoint[f[1]],
                            centreLabel
                        );

                        tets_.append(tet);
                    }
                    else
                    {
                        forAll(f, pI)
                        {
                            partTet tet
                            (
                                nodeLabelForPoint[f[pI]],
                                nodeLabelForPoint[f.prevLabel(pI)],
                                nodeLabelForFace[c[fI]],
                                centreLabel
                            );

                            tets_.append(tet);
                        }
                    }
                }
                else
                {
                    if (f.size() == 3)
                    {
                        partTet tet
                        (
                            nodeLabelForPoint[f[0]],
                            nodeLabelForPoint[f[1]],
                            nodeLabelForPoint[f[2]],
                            centreLabel
                        );

                        tets_.append(tet);
                    }
                    else
                    {
                        forAll(f, pI)
                        {
                            partTet tet
                            (
                                nodeLabelForPoint[f[pI]],
                                nodeLabelForPoint[f.nextLabel(pI)],
                                nodeLabelForFace[c[fI]],
                                centreLabel
                            );

                            tets_.append(tet);
                        }
                    }
                }
            }
        }
    }

    // create node labels in origMesh_
    nodeLabelInOrigMesh_.setSize(points_.size());
    nodeLabelInOrigMesh_ = -1;
    forAll(nodeLabelForPoint, pI)
    {
        if (nodeLabelForPoint[pI] != -1)
        {
            // lock mesh vertices
            if (lockedPoints[pI])
            {
                smoothVertex_[nodeLabelForPoint[pI]] |= LOCKED;
            }

            nodeLabelInOrigMesh_[nodeLabelForPoint[pI]] = pI;
        }
    }

    // create pointTets_
    reverseAddressing(pointTets_, points_.size(), tets_);
}


void Foam::Module::partTetMesh::createSMOOTHPointsOrdering() const
{
    internalPointsOrderPtr_ = new labelListList2D();
    labelListList2D& internalPointsOrder = *internalPointsOrderPtr_;

    internalPointsOrder.setSize(0);
    DynamicList<label> order(points_.size(), -1);
    boolList helper(points_.size());

    bool found;
    do
    {
        found = false;
        helper = false;
        DynamicList<label> selectedPoints;

        forAll(points_, nodeI)
        {
            if (smoothVertex_[nodeI] & SMOOTH)
            {
                if (helper[nodeI])
                {
                    continue;
                }
                if (order[nodeI] != -1)
                {
                    continue;
                }

                // find neighbouring FACECENTRE and CELLCENTRE points
                DynamicList<label> neiCentrePoints, neiSmoothPoints;
                neiCentrePoints.reserve(64); neiSmoothPoints.reserve(64);
                forAllRow(pointTets_, nodeI, ptI)
                {
                    const partTet& tet = tets_[pointTets_[nodeI][ptI]];

                    for (label i = 0; i < 4; ++i)
                    {
                        if (smoothVertex_[tet[i]] & (FACECENTRE + CELLCENTRE))
                        {
                            appendUniq(neiCentrePoints, tet[i]);
                        }
                        else if (smoothVertex_[tet[i]] & SMOOTH)
                        {
                            appendUniq(neiSmoothPoints, tet[i]);
                        }
                    }
                }

                // find neighbouring SMOOTH points
                forAll(neiCentrePoints, ncI)
                {
                    const label centreI = neiCentrePoints[ncI];

                    forAllRow(pointTets_, centreI, ptI)
                    {
                        const partTet& tet = tets_[pointTets_[centreI][ptI]];

                        for (label i = 0; i < 4; ++i)
                        {
                            if (smoothVertex_[tet[i]] & SMOOTH)
                            {
                                appendUniq(neiSmoothPoints, tet[i]);
                            }
                        }
                    }
                }

                // select the point and mark neighbouring SMOOTH points
                selectedPoints.append(nodeI);
                order[nodeI] = internalPointsOrder.size();

                forAll(neiSmoothPoints, i)
                {
                    helper[neiSmoothPoints[i]] = true;
                }
            }
        }

        if (selectedPoints.size() != 0)
        {
            internalPointsOrder.setSize(internalPointsOrder.size() + 1);
            setRow(internalPointsOrder, internalPointsOrder.size() - 1, selectedPoints);
            found = true;
        }

    } while (found);
}


void Foam::Module::partTetMesh::createBOUNDARYPointsOrdering() const
{
    boundaryPointsOrderPtr_ = new labelListList2D();
    labelListList2D& boundaryPointsOrder = *boundaryPointsOrderPtr_;

    boundaryPointsOrder.setSize(0);
    DynamicList<label> order(points_.size(), -1);
    boolList helper(points_.size());

    bool found;
    do
    {
        found = false;
        helper = false;

        DynamicList<label> selectedPoints;
        forAll(points_, nodeI)
        {
            if (smoothVertex_[nodeI] & BOUNDARY)
            {
                if (helper[nodeI])
                {
                    continue;
                }
                if (order[nodeI] != -1)
                {
                    continue;
                }

                // find neighbouring FACECENTRE and CELLCENTRE points
                DynamicList<label> neiCentrePoints, neiSmoothPoints;
                neiCentrePoints.reserve(64); neiSmoothPoints.reserve(64);
                forAllRow(pointTets_, nodeI, ptI)
                {
                    const partTet& tet = tets_[pointTets_[nodeI][ptI]];

                    for (label i = 0; i < 4; ++i)
                    {
                        if (smoothVertex_[tet[i]] & (FACECENTRE + CELLCENTRE))
                        {
                            appendUniq(neiCentrePoints, tet[i]);
                        }
                        else if (smoothVertex_[tet[i]] & BOUNDARY)
                        {
                            appendUniq(neiSmoothPoints, tet[i]);
                        }
                    }
                }

                // find neighbouring BOUNDARY points
                forAll(neiCentrePoints, ncI)
                {
                    const label centreI = neiCentrePoints[ncI];

                    forAllRow(pointTets_, centreI, ptI)
                    {
                        const partTet& tet = tets_[pointTets_[centreI][ptI]];

                        for (label i = 0; i < 4; ++i)
                        {
                            if (smoothVertex_[tet[i]] & BOUNDARY)
                            {
                                appendUniq(neiSmoothPoints, tet[i]);
                            }
                        }
                    }
                }

                // select the point and mark neighbouring  BOUNDARY points
                selectedPoints.append(nodeI);
                order[nodeI] = boundaryPointsOrder.size();

                forAll(neiSmoothPoints, i)
                {
                    helper[neiSmoothPoints[i]] = true;
                }
            }
        }

        if (selectedPoints.size() != 0)
        {
            boundaryPointsOrder.setSize(boundaryPointsOrder.size() + 1);
            setRow(boundaryPointsOrder, boundaryPointsOrder.size() - 1, selectedPoints);
            found = true;
        }

    } while (found);
}


// ************************************************************************* //
