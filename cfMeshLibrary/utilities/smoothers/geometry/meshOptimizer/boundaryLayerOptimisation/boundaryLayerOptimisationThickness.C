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
#include "boundaryLayerOptimisation.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "labelledScalar.H"
#include "polyMeshGenAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::boundaryLayerOptimisation::hairEdgesAtBndFace
(
    const label cellI,
    const label baseFaceI,
    DynamicList<edge>& hairEdges
) const
{
    const faceListPMG& faces = mesh_.faces();

    const cell& c = mesh_.cells()[cellI];

    // check cell topology
    DynamicList<edge> edges;
    edges.reserve(48);
    DynamicList<DynamicList<label>> edgeFaces;
    labelListList2D faceEdges;
    faceEdges.setSize(c.size());
    label baseFace(-1);
    forAll(c, fI)
    {
        if (c[fI] == baseFaceI)
        {
            baseFace = fI;
        }

        const face& f = faces[c[fI]];
        faceEdges[fI].setSize(f.size());

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            label pos = findIndex(edges, e);

            if (pos < 0)
            {
                pos = edges.size();
                edges.append(e);
                edgeFaces.setSize(pos + 1);
            }

            edgeFaces[pos].append(fI);
            faceEdges[fI][eI] = pos;
        }
    }

    const face& bf = faces[c[baseFace]];

    // built up point-by-point via append(), not pre-sized to bf.size() -
    // on an irregular cell (e.g. a boundary-layer transition cell created
    // by allowDiscontinuity) the loop below may stop before reaching
    // every point, and the caller must be able to tell how far it
    // actually got rather than trust the rest to be valid
    hairEdges.clear();
    hairEdges.setCapacity(bf.size());

    forAll(bf, pI)
    {
        const label nextEdge = faceEdges[baseFace][pI];
        const label prevEdge = faceEdges[baseFace][bf.rcIndex(pI)];

        if (edgeFaces[nextEdge].size() != 2 || edgeFaces[prevEdge].size() != 2)
        {
            break;
        }

        // find the face attached to the edge after the current point
        label otherNextFace = edgeFaces[nextEdge][0];
        if (otherNextFace == baseFace)
        {
            otherNextFace = edgeFaces[nextEdge][1];
        }

        // find the face attached to the edge before the current point
        label otherPrevFace = edgeFaces[prevEdge][0];
        if (otherPrevFace == baseFace)
        {
            otherPrevFace = edgeFaces[prevEdge][1];
        }

        // the hair edge at bf[pI] must actually touch bf[pI] - on a
        // simple hex-like cell, otherNextFace and otherPrevFace share
        // exactly one edge (the hair edge itself), but on an irregular
        // cell they can share more than one, and only the one touching
        // bf[pI] (if any) is a valid candidate; the others are
        // unrelated coincidental adjacencies elsewhere in the cell
        label commonEdge(-1);
        forAll(edges, eI)
        {
            if
            (
                (edges[eI][0] == bf[pI] || edges[eI][1] == bf[pI])
             && (findIndex(edgeFaces[eI], otherNextFace) != -1)
             && (findIndex(edgeFaces[eI], otherPrevFace) != -1)
            )
            {
                commonEdge = eI;
                break;
            }
        }

        if (commonEdge == -1)
        {
            break;
        }

        // there exists a common edge which shall be used as a hair
        if (edges[commonEdge].start() == bf[pI])
        {
            hairEdges.append(edges[commonEdge]);
        }
        else
        {
            hairEdges.append(edges[commonEdge].reverseEdge());
        }
    }
}


Foam::scalar Foam::Module::boundaryLayerOptimisation::calculateThickness
(
    const label heI,
    const label heJ
) const
{
    const pointFieldPMG& points = mesh_.points();

    // references to hair edges
    const edge& he = hairEdges_[heI];
    const edge& nhe = hairEdges_[heJ];

    // distance vector between the surface points of hair edges
    const point& sp = points[he[0]];
    const point& ep = points[nhe[0]];
    const vector dv = ep - sp;
    const scalar magDv = mag(dv);

    // calculate layer thickness
    const scalar currThickness = he.mag(points);
    scalar retThickness = currThickness;

    const scalar currNeiThickness = nhe.mag(points);
    scalar suggestedNeiThickness = currNeiThickness;

    // calculate layer height at the current point
    const point npAlpha = help::nearestPointOnTheEdge(sp, ep, points[he[1]]);
    const scalar currHeight = mag(npAlpha - points[he[1]]);
    scalar retHeight = currHeight;
    const scalar cosAlpha = sign((npAlpha - sp) & dv)*mag(npAlpha - sp);
    const scalar alpha =
        Foam::acos
        (
            Foam::max
            (
                -1.0,
                Foam::min(1.0, cosAlpha /Foam::max(currThickness, VSMALL))
            )
        );

    // calculate the height of the layer at the neighbour
    // point
    const point npBeta = help::nearestPointOnTheEdge(ep, sp, points[nhe[1]]);
    const scalar currNeiHeight = mag(npBeta - points[nhe[1]]);
    scalar suggestedNeiHeight = currNeiHeight;
    const scalar cosBeta = sign((npBeta - ep) & -dv)*mag(npBeta - ep);
    const scalar beta =
        Foam::acos
        (
            Foam::max
            (
                -1.0,
                Foam::min(1.0, cosBeta /Foam::max(currNeiThickness, VSMALL))
            )
        );

    // check if the current thickness is Ok for the local curvature
    if ((alpha + beta) < M_PI)
    {
        const scalar gamma = M_PI - (alpha + beta);
        const scalar sinGamma = Foam::max(SMALL, Foam::sin(gamma));
        const scalar sinAlpha = Foam::max(SMALL, Foam::sin(alpha));
        const scalar sinBeta = Foam::max(SMALL, Foam::sin(beta));

        // max allowed thickness and layer height due to curvature
        retThickness =
            Foam::min
            (
                retThickness,
                featureSizeFactor_*magDv*sinBeta/sinGamma
            );
        retHeight *= (retThickness /Foam::max(currThickness, VSMALL));

        // max allowed neighbour hair thickness
        // and layer height due to curvature
        suggestedNeiThickness =
            Foam::min
            (
                suggestedNeiThickness,
                featureSizeFactor_*magDv*sinAlpha/sinGamma
            );
        suggestedNeiHeight *=
            (suggestedNeiThickness /Foam::max(currNeiThickness, VSMALL));
    }

    // check the height variation
    const scalar tanVal =
        (retHeight - suggestedNeiHeight)/Foam::max(magDv, VSMALL);

    if (tanVal > relThicknessTol_)
    {
        retHeight = suggestedNeiHeight + relThicknessTol_*magDv;

        retThickness = (retHeight/currHeight)*currThickness;
    }

    return retThickness;
}


Foam::scalar Foam::Module::boundaryLayerOptimisation::calculateThicknessOverCell
(
    const label heI,
    const label cellI,
    const label baseFaceI
) const
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    const cell& c = mesh_.cells()[cellI];

    const face& bf = faces[baseFaceI];

    const edge& he = hairEdges_[heI];

    const point& sp = points[he[0]];
    const point& ep = points[he[1]];

    scalar maxThickness = he.mag(points);

    // the base face must not contain the hair edge
    // this is the case at exitting layers
    forAll(bf, eI)
    {
        if (bf.faceEdge(eI) == he)
        {
            return maxThickness;
        }
    }

    forAll(c, fI)
    {
        if (c[fI] == baseFaceI)
        {
            continue;
        }

        const face& f = faces[c[fI]];

        if (help::shareAnEdge(bf, f) && (f.which(he.start()) == -1))
        {
            point intersection;

            if (!help::lineFaceIntersection(sp, ep, f, points, intersection))
            {
                continue;
            }

            const scalar maxDist = featureSizeFactor_*mag(intersection - sp);

            maxThickness =
                Foam::min(maxThickness, maxDist);
        }
    }

    return maxThickness;
}


void Foam::Module::boundaryLayerOptimisation::optimiseThicknessVariation
(
    const direction edgeType
)
{
    pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const labelListList2D& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const label start = mesh_.nInternalFaces();
    const labelList& faceOwner = mse.faceOwners();

    vectorField hairDirections(hairEdges_.size());
    scalarField hairLength(hairEdges_.size());

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(hairEdges_, hairEdgeI)
    {
        vector n = hairEdges_[hairEdgeI].vec(points);

        hairLength[hairEdgeI] = (Foam::mag(n) + VSMALL);
        hairDirections[hairEdgeI] = n/hairLength[hairEdgeI];
    }

    // reduce thickness of the layer
    // such that the variation of layer thickness
    // It is an iterative process where the layer is thinned in the regions
    // where the tangent is greater than the tolerance value or the curvature
    // permits thicker boundary layers.
    boolList activeHairEdge(hairEdges_.size(), true);
    bool changed;
    label nIter(0);
    do
    {
        changed = false;

        boolList modifiedHairEdge(hairEdges_.size(), false);
        boolList influencedEdges(hairEdges_.size(), false);

        // check if the hair edge intersects some other face in the cells
        // attached to the hair edge
        # pragma omp parallel for schedule(dynamic, 50)
        forAll(hairEdges_, hairEdgeI)
        {
            if
            (
                (hairEdgeType_[hairEdgeI] & edgeType) &&
                activeHairEdge[hairEdgeI]
            )
            {
                const edge& he = hairEdges_[hairEdgeI];

                const label bpI = bp[he.start()];

                scalar maxThickness = hairLength[hairEdgeI];

                DynamicList<label> influencers;
                influencers.reserve(64);

                forAllRow(pFaces, bpI, pfI)
                {
                    const label bfI = pFaces[bpI][pfI];

                    const face& bf = bFaces[bfI];
                    if (bf.which(he.end()) >= 0)
                    {
                        continue;
                    }

                    const label baseFaceI = start + bfI;
                    const label cOwn = faceOwner[bfI];

                    // check if there exist any self-intersections
                    maxThickness =
                        Foam::min
                        (
                            maxThickness,
                            calculateThicknessOverCell
                            (
                                hairEdgeI,
                                cOwn,
                                baseFaceI
                            )
                        );

                    // check thickness variation over all hair edges -
                    // hairEdges may have fewer entries than bf has
                    // points (e.g. an irregular boundary-layer
                    // transition cell), so iterate it directly rather
                    // than assuming it matches bf point-for-point
                    DynamicList<edge> hairEdges;
                    hairEdgesAtBndFace(faceOwner[bfI], baseFaceI, hairEdges);

                    forAll(hairEdges, pI)
                    {
                        const edge& nhe = hairEdges[pI];

                        const label bpJ = bp[nhe.start()];

                        forAllRow(hairEdgesAtBndPoint_, bpJ, peJ)
                        {
                            const label heJ = hairEdgesAtBndPoint_[bpJ][peJ];

                            if (hairEdgeI == heJ)
                            {
                                continue;
                            }

                            if (nhe == hairEdges_[heJ])
                            {
                                influencers.append(heJ);

                                const scalar edgeThickness =
                                    calculateThickness(hairEdgeI, heJ);

                                maxThickness =
                                    Foam::min(maxThickness, edgeThickness);
                            }
                        }
                    }
                }

                if (hairLength[hairEdgeI] > maxThickness)
                {
                    // make the hair edge shorter
                    hairLength[hairEdgeI] = maxThickness;
                    modifiedHairEdge[hairEdgeI] = true;
                    changed = true;

                    thinnedHairEdge_[hairEdgeI] = true;

                    forAll(influencers, i)
                    {
                        influencedEdges[influencers[i]] = true;
                    }
                }
            }
        }

        if (!changed)
        {
            break;
        }

        // move boundary vertices to the new positions
        # pragma omp parallel for schedule(dynamic, 100)
        forAll(hairEdges_, hairEdgeI)
        {
            if (hairEdgeType_[hairEdgeI] & edgeType)
            {
                const edge& he = hairEdges_[hairEdgeI];
                const vector& hv = hairDirections[hairEdgeI];
                const point& s = points[he.start()];

                points[he.end()] = s + hairLength[hairEdgeI]*hv;
            }
        }

        // mark edges which may be changed
        activeHairEdge.transfer(influencedEdges);
    } while (changed && (++nIter < 1000));
}


// ************************************************************************* //
