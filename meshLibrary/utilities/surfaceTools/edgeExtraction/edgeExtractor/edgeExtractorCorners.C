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

#include "error.H"
#include "polyMeshGenModifier.H"
#include "edgeExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceOptimizer.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"
#include "DynList.H"
#include "labelPair.H"
#include "labelledScalar.H"
#include "labelledPoint.H"
#include "HashSet.H"
#include "triSurfacePartitioner.H"
#include "triSurfaceClassifyEdges.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//


void Foam::Module::edgeExtractor::faceEvaluator::neiFacesOverEdges
(
    const label bfI,
    DynList<label>& neiFaces
) const
{
    const meshSurfaceEngine& mse = extractor_.surfaceEngine();

    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    neiFaces.setSize(faceEdges.sizeOfRow(bfI));

    forAllRow(faceEdges, bfI, feI)
    {
        const label beI = faceEdges(bfI, feI);

        if (edgeFaces.sizeOfRow(beI) == 2)
        {
            neiFaces[feI] = edgeFaces(beI, 0);
            if (neiFaces[feI] == bfI)
                neiFaces[feI] = edgeFaces(beI, 1);
        }
        else
        {
            neiFaces[feI] = -1;
        }
    }
}


void Foam::Module::edgeExtractor::faceEvaluator::neiPatchesOverEdges
(
    const label bfI,
    const labelList& fPatches,
    DynList<label> &neiPatches
) const
{
    const meshSurfaceEngine& mse = extractor_.surfaceEngine();

    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    neiPatches.setSize(faceEdges.sizeOfRow(bfI));

    forAllRow(faceEdges, bfI, feI)
    {
        const label beI = faceEdges(bfI, feI);

        if (edgeFaces.sizeOfRow(beI) == 2)
        {
            label nei = edgeFaces(beI, 0);
            if (nei == bfI)
                nei = edgeFaces(beI, 1);

            neiPatches[feI] = fPatches[nei];
        }
    }
}


Foam::label Foam::Module::edgeExtractor::faceEvaluator::bestPatchTopological
(
    const DynList<label>& neiPatches,
    const label currentPatch
)
{
    // find indices of all neighbour patches
    DynList<label> allNeiPatches;
    forAll(neiPatches, i)
        allNeiPatches.appendUniq(neiPatches[i]);

    if ((allNeiPatches.size() == 1) && (allNeiPatches[0] == currentPatch))
        return currentPatch;

    // counter the number of neighbours in a patch
    Map<label> nNeiInPatch(allNeiPatches.size());
    forAll(allNeiPatches, i)
        nNeiInPatch.insert(allNeiPatches[i], 0);
    forAll(neiPatches, eI)
        ++nNeiInPatch[neiPatches[eI]];

    label newPatch = -1;
    label nNeiEdges(0);
    forAllConstIter(Map<label>, nNeiInPatch, it)
    {
        if (it() > nNeiEdges)
        {
            newPatch = it.key();
            nNeiEdges = it();
        }
        else if
        (
            (it() == nNeiEdges) && (it.key() == currentPatch)
        )
        {
            newPatch = it.key();
        }
    }

    // do not swap if the situation allows for more than one edge
    // shared with faces in other patches than the dominant one
    if (nNeiEdges < (neiPatches.size() - 1))
        return currentPatch;

    // do not swap if the best face is in the current patch
    if ((newPatch < 0) || (newPatch == currentPatch))
        return currentPatch;

    // check whether the edges shared ith the neighbour patch form
    // a singly linked chain
    DynList<bool> sharedEdge;
    sharedEdge.setSize(neiPatches.size());
    sharedEdge = false;

    forAll(neiPatches, eI)
        if (neiPatches[eI] == newPatch)
            sharedEdge[eI] = true;

    if (help::areElementsInChain(sharedEdge))
        return newPatch;

    return currentPatch;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

Foam::Module::edgeExtractor::faceEvaluator::faceEvaluator
(
    const edgeExtractor& ee
)
:
    extractor_(ee),
    newBoundaryPatchesPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

Foam::Module::edgeExtractor::faceEvaluator::~faceEvaluator()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void Foam::Module::edgeExtractor::faceEvaluator::setNewBoundaryPatches
(
    const labelList& newBoudaryPatches
)
{
    newBoundaryPatchesPtr_ = &newBoudaryPatches;
}


void Foam::Module::edgeExtractor::faceEvaluator::neiPatchesOverEdges
(
    const label bfI,
    DynList<label>& neiPatches
) const
{
    neiPatchesOverEdges
    (
        bfI,
        extractor_.facePatch_,
        neiPatches
    );
}


Foam::label Foam::Module::edgeExtractor::faceEvaluator::bestPatchTopological
(
    const label bfI
) const
{
    // get neighbour patches over edges
    DynList<label> neiPatches;
    neiPatchesOverEdges
    (
        bfI,
        extractor_.facePatch_,
        neiPatches
    );

    return bestPatchTopological(neiPatches, extractor_.facePatch_[bfI]);
}


Foam::label
Foam::Module::edgeExtractor::faceEvaluator::bestPatchAfterModification
(
    const label bfI
) const
{
    const label patchI = newBoundaryPatchesPtr_->operator[](bfI);

    if (patchI != extractor_.facePatch_[bfI])
    {
        DynList<label> newNeiPatches, oldNeiPatches;
        neiPatchesOverEdges
        (
            bfI,
            *newBoundaryPatchesPtr_,
            newNeiPatches
        );

        neiPatchesOverEdges
        (
            bfI,
            extractor_.facePatch_,
            oldNeiPatches
        );

        DynList<label> neiFaces;
        neiFacesOverEdges(bfI, neiFaces);

        forAll(neiFaces, eI)
        {
            const label origPatchI = oldNeiPatches[eI];
            const label newPatchI = newNeiPatches[eI];

            if (neiFaces[eI] > bfI)
            {
                if (newPatchI != origPatchI)
                    newNeiPatches[eI] = origPatchI;
            }
        }

        return bestPatchTopological(newNeiPatches, patchI);
    }

    return patchI;
}


void Foam::Module::edgeExtractor::cornerEvaluator::sortedFacesAtPoint
(
    const label bpI,
    DynList<label>& pFaces
) const
{
    const meshSurfaceEngine& mse = extractor_.surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& pointInFace = mse.pointInFaces();

    pFaces = pointFaces[bpI];

    forAll(pFaces, i)
    {
        const face& bf = bFaces[pFaces[i]];
        const label pos = pointFaces.find(bpI, pFaces[i]);

        const edge e = bf.faceEdge(bf.rcIndex(pointInFace(bpI, pos)));

        for (label j = i + 2; j < pFaces.size(); ++j)
        {
            if (bFaces[pFaces[j]].which(e.start()) >= 0)
            {
                label bfJ = pFaces[i + 1];
                pFaces[i + 1] = pFaces[j];
                pFaces[j] = bfJ;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

Foam::Module::edgeExtractor::cornerEvaluator::cornerEvaluator
(
    const edgeExtractor& ee,
    const meshSurfacePartitioner& mPart
)
:
    extractor_(ee),
    partitioner_(mPart)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void Foam::Module::edgeExtractor::cornerEvaluator::improveCorners
(
    labelList& newBoundaryPatches
)
{
    const meshOctree& mo = extractor_.meshOctree_;

    const labelHashSet& corners = partitioner_.corners();

    const meshSurfaceEngine& mse = extractor_.surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const pointFieldPMG& points = mse.points();
    const faceList::subList& bFaces = mse.boundaryFaces();

    for (const label corneri : corners.toc())
    {
        const point& p = points[bPoints[corneri]];

        {
            const labelList& facePatch = extractor_.facePatch_;
            Info<< "Checking corner " << corneri << endl;

            // get faces attached to this point
            // sorted in the counter-clocwise order
            DynList<label> pFaces;
            sortedFacesAtPoint(corneri, pFaces);

            Info<< "Faces at corner " << pFaces << endl;

            // find feature edges attached to the point
            DynList<labelPair> edgePatches;
            DynList<label> edgeIndex;
            forAll(pFaces, i)
            {
                const label patch0 = facePatch[pFaces[i]];
                const label patch1 = facePatch[pFaces[pFaces.fcIndex(i)]];

                if (patch0 != patch1)
                {
                    edgeIndex.append(i);
                    edgePatches.append(labelPair(patch0, patch1));
                }
            }

            Info<< "Edge patches " << edgePatches << endl;
            Info<< "Edge indices " << edgeIndex << endl;

            // find the best aligned feature edge to the feature edge
            DynList<label> bestEdgeIndices(edgePatches.size());
            forAll(edgePatches, i)
            {
                const label patch0 = edgePatches[i].first();
                const label patch1 = edgePatches[i].second();
                DynList<label> patches(2);
                patches[0] = patch0;
                patches[1] = patch1;

                // find the proection of the first point onto the feature edge
                point mps;
                scalar dSqS;
                mo.findNearestPointToPatches(mps, dSqS, p, patches);

                scalar bestAlignment(0.0);
                label bestEdgeIndex(-1);

                forAll(pFaces, pfI)
                {
                    const face& bf = bFaces[pFaces[pfI]];

                    const label pos = bf.which(bPoints[corneri]);

                    const point& pe = points[bf[bf.rcIndex(pos)]];

                    // calculate the vector of current edge
                    vector ev = pe - p;
                    ev /= (mag(ev) + VSMALL);

                    // project the endpoint onto the feature edge
                    point mpe;
                    scalar dSqE;
                    mo.findNearestPointToPatches(mpe, dSqE, pe, patches);

                    vector fv = mpe - mps;
                    fv /= (mag(fv) + VSMALL);

                    // calculate alignment metrix between the current edge
                    // and the feature edge
                    const scalar align = 0.5*(1.0 + (ev & fv));

                    if (align > bestAlignment)
                    {
                        bestAlignment = align;
                        bestEdgeIndex = pfI;
                    }
                }

                bestEdgeIndices[i] = bestEdgeIndex;
            }

            Info<< "\nPatches of edges at corner " << corneri
                << " are " << edgePatches << endl;
            Info<< "Orig edge indices " << edgeIndex << endl;
            Info<< "New edge indices " << bestEdgeIndices << endl;

            DynList<label> newFacePatches(pFaces.size());
            forAll(bestEdgeIndices, i)
            {
                const label patch = edgePatches[i].first();
                label index = bestEdgeIndices[i];

                while (index != bestEdgeIndices[bestEdgeIndices.fcIndex(i)])
                {
                    newFacePatches[index] = patch;
                    index = (index + 1) % pFaces.size();
                }
            }

            bool allPatchesRemain(true);
            forAll(edgePatches, i)
            {
                if (findIndex(newFacePatches, edgePatches[i].first()) == -1)
                {
                    allPatchesRemain = false;
                    break;
                }
            }

            Info<< "New patches at point" << newFacePatches << endl;
            if (allPatchesRemain)
            {
                forAll(pFaces, i)
                    newBoundaryPatches[pFaces[i]] = newFacePatches[i];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

bool Foam::Module::edgeExtractor::findCornerCandidates()
{
    bool changed(false);

    const triSurf& surf = meshOctree_.surface();
    const pointField& sp = surf.points();

    // create the surface partitioner and find the corners on the surface mesh
    triSurfacePartitioner surfPartitioner(surf);

    const labelList& corners = surfPartitioner.corners();

    Map<label> surfacePointToCorner;
    forAll(corners, cornerI)
        surfacePointToCorner.insert(corners[cornerI], cornerI);

    List<labelledScalar> nearestPoint
    (
        corners.size(),
        labelledScalar(0, VGREAT)
    );

    // calculate the search range
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const pointFieldPMG& points = mesh_.points();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const faceList::subList& bFaces = mse.boundaryFaces();

    scalarList searchRange(bPoints.size(), 0.0);
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];
        const point c = bf.centre(static_cast<const pointField&>(points));

        forAll(bf, pI)
        {
            const scalar d = 2.0*Foam::mag(c - points[bf[pI]]);
            const label bpI = bp[bf[pI]];

            searchRange[bpI] = Foam::max(d, searchRange[bpI]);
        }
    }

    DynList<label> containedTriangles;
    # pragma omp parallel for schedule(dynamic, 40) private(containedTriangles)
    forAll(bPoints, bpI)
    {
        const point& p = points[bPoints[bpI]];
        const scalar s = searchRange[bpI];

        boundBox bb(p - vector(s, s, s), p + vector(s, s, s));
        meshOctree_.findTrianglesInBox(bb, containedTriangles);

        // find the nearest corner on the surface mesh
        forAll(containedTriangles, i)
        {
            const labelledTri& tri = surf[containedTriangles[i]];

            forAll(tri, pI)
            {
                const label spI = tri[pI];
                if (!surfacePointToCorner.found(spI))
                    continue;

                const label cornerI = surfacePointToCorner[spI];

                const scalar d = Foam::magSqr(sp[spI] - points[bPoints[bpI]]);

                if (nearestPoint[cornerI].value() > d)
                {
                    // update the nearest point to the found corner
                    # pragma omp critical
                    {
                        nearestPoint[cornerI] = labelledScalar(bpI, d);
                    }
                }
            }
        }
    }

    return changed;
}


bool Foam::Module::edgeExtractor::checkCorners()
{
    bool changed(false);

    const pointFieldPMG& points = mesh_.points();
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    const triSurfacePartitioner& sPartitioner = this->partitioner();
    const List<labelHashSet>& patchPatches = sPartitioner.patchPatches();

    // allocate a copy of boundary patches
    labelList newBoundaryPatches(facePatch_.size());

    label nCorrected;
    Map<label> otherProcNewPatch;

    label nIteration(0);

    do
    {
        nCorrected = 0;
        newBoundaryPatches = facePatch_;

        // update the information which edges are assigned as feature edges
        meshSurfacePartitioner mPart(mse, newBoundaryPatches);

        // find corners in the current constelation
        const labelHashSet& corners = mPart.corners();
        const VRWGraph& pPatches = mPart.pointPatches();

        typedef std::map<label, label> mapType;
        typedef std::map<label, std::pair<point, scalar>> distMapType;

        mapType cornerIndex;
        distMapType nearestToCorner;

        // find nearest corners in the surface mesh
        for (const label bpI : corners.toc())
        {
            const DynList<label> neiPatches = pPatches[bpI];

            const point& p = points[bPoints[bpI]];
            point pMap;
            scalar dSq;
            label nsp;

            if (!meshOctree_.findNearestCorner(pMap, dSq, nsp, p, neiPatches))
            {
                nsp = -1;
                meshOctree_.findNearestPointToPatches(pMap, dSq, p, neiPatches);
            }

            nearestToCorner[bpI] = std::make_pair(pMap, dSq);
            cornerIndex[bpI] = nsp;
        }

        std::map<label, DynList<label>> bpMappedAtSurfaceCorner;
        forAllConstIter(mapType, cornerIndex, mIt)
        {
            if (mIt->second < 0)
            {
                // the corner does not exist in the surface mesh
                // this may be a situation where parts of the surface are in
                // close proximity and are not topologically connected
                Warning << "Should not get in here. Not implemented" << endl;
            }
            else
            {
                bpMappedAtSurfaceCorner[mIt->second].append(mIt->first);
            }
        }

        // for all edge nodes find their nearest counterparts in the surface

        const labelHashSet& featureEdge = mPart.featureEdges();

        distMapType nearestToEdgePoint;
        mapType edgePointIndex;
        std::map<std::pair<label, label>, labelLongList> edgesSharedByPatches;

        for (const label beI : featureEdge.toc())
        {
            // get the patches bounded by this feature edge
            DynList<label> patches(2);

            if (edgeFaces.sizeOfRow(beI) == 2)
            {
                patches[0] = newBoundaryPatches[edgeFaces(beI, 0)];
                patches[1] = newBoundaryPatches[edgeFaces(beI, 1)];
            }
            else if (edgeFaces.sizeOfRow(beI) == 1)
            {
                patches[0] = newBoundaryPatches[edgeFaces(beI, 0)];
                patches[1] = otherProcNewPatch[beI];
            }

            // store the edge into the right group
            std::pair<label, label> patchPair
            (
                Foam::min(patches[0], patches[1]),
                Foam::max(patches[0], patches[1])
            );
            edgesSharedByPatches[patchPair].append(beI);

            // check if some points have already been checked
            const edge& e = edges[beI];
            const label bps = bp[e.start()];
            const label bpe = bp[e.end()];

            DynList<label> checkPoints;
            if
            (
                corners.found(bps) &&
                (nearestToEdgePoint.find(bps) == nearestToEdgePoint.end())
            )
                checkPoints.append(bps);

            if
            (
                corners.found(bpe) &&
                (nearestToEdgePoint.find(bpe) == nearestToEdgePoint.end())
            )
                checkPoints.append(bpe);

            forAll(checkPoints, i)
            {
                const point& p = points[bPoints[checkPoints[i]]];
                point pMap;
                scalar dSq;
                label nse;

                if (!meshOctree_.findNearestEdgePoint(pMap, dSq, nse, p, patches))
                {
                    nse = -1;
                    meshOctree_.findNearestPointToPatches(pMap, dSq, p, patches);
                }

                nearestToEdgePoint[checkPoints[i]] = std::make_pair(pMap, dSq);
                edgePointIndex[checkPoints[i]] = nse;
            }
        }

        labelHashSet invalidFeatureEdges;
        for (std::map<std::pair<label, label>, labelLongList>::const_iterator mIt = edgesSharedByPatches.cbegin(); mIt != edgesSharedByPatches.cend(); ++mIt)
        {
            const bool validConnection =
                patchPatches[mIt->first.first].found(mIt->first.second);

            if (!validConnection)
            {
                // find surface facets in the vicinity of the edge
                // and check whether there exist various disconnected surface
                // parts in the vicinity of this group of edge

                const auto& invalidEdges = mIt->second;

                forAll(invalidEdges, i)
                    invalidFeatureEdges.insert(invalidEdges[i]);
            }
        }

        // check the vicinity of the corner point and check whether
        // it shall be replaced by some other point in the vicinity
        std::map<label, DynList<labelPair>> facesAtCornerNewPatches;
        forAllConstIter(distMapType, nearestToCorner, it)
        {
            const label bpI = it->first;

            // find all points connected to the corner point via a face
            labelHashSet nearPoints;
            std::map<label, DynList<label>> facesContainingPoint;
            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                const face& bf = bFaces[bfI];

                forAll(bf, pI)
                {
                    nearPoints.insert(bf[pI]);
                    facesContainingPoint[bf[pI]].append(bfI);
                }
            }

            // find the nearest point to the location where the corner
            // shall be mapped onto the surface mesh
            label bestPoint(-1);
            scalar minDistSq(VGREAT);
            for (const label pointi : nearPoints.toc())
            {
                const point& p = points[pointi];

                const scalar dSq = magSqr(p - nearestToCorner[bpI].first);

                if (dSq < minDistSq)
                {
                    minDistSq = dSq;
                    bestPoint = pointi;
                }
            }

            if (bestPoint == bPoints[bpI])
            {
                continue;
            }

            // sort faces and edges at the corner in counter-clockwise order
            DynList<label> pFaces, pEdges;

            DynList<label> front;
            front.append(0);

            while (front.size())
            {
                const label fI = front.remove();
                const label bfI = pointFaces(bpI, fI);

                pFaces.append(bfI);

                const face& bf = bFaces[bfI];
                const label pos = bf.which(bPoints[bpI]);

                const label beI = faceEdges(bfI, bf.rcIndex(pos));

                pEdges.append(beI);

                label nei = edgeFaces(beI, 0);
                if (nei == bfI)
                    nei = edgeFaces(beI, 1);

                if ((findIndex(pFaces, nei) != -1) || (nei < 0))
                    continue;

                front.append(pointFaces.find(bpI, nei));
            }

            // find the best fitting edge to move towards the corner
            label bestEdge(-1);
            scalar bestAlignment(0.0);

            vector dv = it->second.first - points[bPoints[bpI]];
            dv /= (mag(dv) + VSMALL);

            if (facesContainingPoint[bestPoint].size() == 1)
            {
                const label bfI = facesContainingPoint[bestPoint][0];

                // only the edges of this face can be candidates
                forAll(pEdges, i)
                {
                    const label beI = pEdges[i];
                    const edge& e = edges[pEdges[i]];

                    if (!(edgeFaces(beI, 0) == bfI || edgeFaces(beI, 1) == bfI))
                        continue;

                    vector ev
                    (
                        points[e.otherVertex(bPoints[bpI])] - points[bPoints[bpI]]
                    );
                    ev /= (mag(ev) + VSMALL);

                    const scalar metric = 0.5*(1.0 + (dv & ev));

                    if (metric > bestAlignment)
                    {
                        bestAlignment = metric;
                        bestEdge = pEdges[i];
                    }
                }
            }
            else
            {
                // the edge containing the best point is the candidate
                forAll(pEdges, i)
                {
                    const edge& e = edges[pEdges[i]];

                    if (e.otherVertex(bestPoint) != -1)
                    {
                        // select the edge which contains best fitting point
                        bestEdge = pEdges[i];
                        break;
                    }
                }
            }

            // find groups of sorted faces at this corner which are bounded by
            // existing feature edges and the new candidate
            DynList<label> faceInGroup(pFaces.size(), -1);
            label groupI(0);
            forAll(pFaces, i)
            {
                if (faceInGroup[i] != -1)
                    continue;

                DynList<label> front;
                front.append(i);
                faceInGroup[i] = groupI;

                while (front.size() != 0)
                {
                    const label currIndex = front.remove();

                    const label nbeI = pEdges[currIndex];
                    const label pbeI = pEdges[pFaces.rcIndex(currIndex)];
                    if ((nbeI != bestEdge) && !featureEdge.found(nbeI))
                    {
                        const label nfI = pFaces.fcIndex(currIndex);

                        if (faceInGroup[nfI] == -1)
                        {
                            faceInGroup[nfI] = groupI;
                            front.append(nfI);
                        }
                    }
                    else if ((pbeI != bestEdge) && !featureEdge.found(pbeI))
                    {
                        const label pfI = pFaces.rcIndex(currIndex);

                        if (faceInGroup[pfI] == -1)
                        {
                            faceInGroup[pfI] = groupI;
                            front.append(pfI);
                        }
                    }
                }

                ++groupI;
            }

            // check which group of faces shall change patch in order to make
            // the best fitting edge a feature edge
            DynList<labelPair> groupPairs, groupPatches;
            labelPair groupsForChanging(-1, -1);
            forAll(pEdges, i)
            {
                const label beI = pEdges[i];

                if (beI == bestEdge)
                {
                    label pos = findIndex(pFaces, edgeFaces(beI, 0));
                    groupsForChanging.first() = faceInGroup[pos];
                    pos = findIndex(pFaces, edgeFaces(beI, 1));
                    groupsForChanging.second() = faceInGroup[pos];
                }
                else if (featureEdge.found(beI))
                {
                    labelPair lp, lpp;
                    label pos = findIndex(pFaces, edgeFaces(beI, 0));
                    lpp.first() = facePatch_[pFaces[pos]];
                    lp.first() = faceInGroup[pos];
                    pos = findIndex(pFaces, edgeFaces(beI, 1));
                    lpp.second() = facePatch_[pFaces[pos]];
                    lp.second() = faceInGroup[pos];

                    groupPairs.append(lp);
                    groupPatches.append(lpp);
                }
            }

            // check which groups shall change their patch
            // desired patches at the best edge are determined by finding
            // the best alignment between the bestEdge and the feature edges
            scalar maxAlignment(0.0);
            label bestGroupsOfPatches(-1);
            forAll(groupPatches, i)
            {
                const labelPair& pp = groupPatches[i];

                const point& bes = points[edges[bestEdge].start()];
                const point& bee = points[edges[bestEdge].end()];

                DynList<label> patches(2);
                patches[0] = pp.first();
                patches[1] = pp.second();

                point mps, mpe;
                scalar dSqS, dSqE;
                label nse;

                meshOctree_.findNearestEdgePoint(mps, dSqS, nse, bes, patches);
                meshOctree_.findNearestEdgePoint(mpe, dSqE, nse, bee, patches);

                const scalar align = 0.5*(1.0 + ((bee - bes) & (mpe - mps)));

                if (align > maxAlignment)
                {
                    maxAlignment = align;
                    bestGroupsOfPatches = i;
                }
            }

            // assign new patches
            FixedList<label, 2> newPatchForGroup(-1);
            DynList<labelPair>& newPatchForFace = facesAtCornerNewPatches[bpI];
            forAll(groupPatches, i)
            {
                if (i == bestGroupsOfPatches)
                    continue;

                const labelPair& gp = groupPairs[i];
                const labelPair& gpp = groupPatches[i];

                if
                (
                    gp.first() == groupsForChanging.first() ||
                    gp.first() == groupsForChanging.second()
                )
                {
                    const label otherPatch = gpp.second();

                    scalar Eold(0.0), Enew(0.0);
                    forAll(faceInGroup, j)
                    {
                        if (faceInGroup[j] != gp.first())
                            continue;

                        const label bfI = pFaces[j];

                        forAllRow(faceEdges, bfI, feI)
                        {
                            const label beI = faceEdges(bfI, feI);

                            const point& ps = points[edges[beI].start()];
                            const point& pe = points[edges[beI].end()];
                            const scalar magE = edges[beI].mag(points) + VSMALL;

                            vector ev = pe - ps;
                            ev /= magE;

                            label nei = edgeFaces(beI, 0);
                            if (nei == bfI)
                                nei = edgeFaces(beI, 1);

                            const label posNei = findIndex(pFaces, nei);
                            if
                            (
                                (posNei < 0) ||
                                (faceInGroup[posNei] != faceInGroup[j])
                            )
                            {
                                if (facePatch_[bfI] != facePatch_[nei])
                                {
                                    DynList<label> patches(2);
                                    patches[0] = facePatch_[bfI];
                                    patches[1] = facePatch_[nei];
                                    // calculate deformation energy
                                    // of the old state
                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mps,
                                        dSqS,
                                        ps,
                                        patches
                                    );
                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mpe,
                                        dSqE,
                                        pe,
                                        patches
                                    );

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Eold +=
                                        1.0/magE*(dSqS + dSqE) + magE*angle;
                                }
                                if (otherPatch != facePatch_[nei])
                                {
                                    // calculate deformation energy
                                    // of the new state
                                    DynList<label> patches(2);
                                    patches[0] = otherPatch;
                                    patches[1] = facePatch_[nei];

                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mps,
                                        dSqS,
                                        ps,
                                        patches
                                    );
                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mpe,
                                        dSqE,
                                        pe,
                                        patches
                                    );

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Enew +=
                                        1.0/magE*(dSqS + dSqE) + magE*angle;
                                }
                            }
                        }
                    }

                    if (Enew <= Eold)
                    {
                        newPatchForGroup[0] = otherPatch;

                        forAll(faceInGroup, j)
                        {
                            if (faceInGroup[j] != gp.first())
                                continue;

                            newPatchForFace.append
                            (
                                labelPair(pFaces[j], otherPatch)
                            );
                        }
                    }
                }
                else if
                (
                    gp.second() == groupsForChanging.first() ||
                    gp.second() == groupsForChanging.second()
                )
                {
                    const label otherPatch = gpp.first();

                    scalar Eold(0.0), Enew(0.0);
                    forAll(faceInGroup, j)
                    {
                        if (faceInGroup[j] != gp.second())
                            continue;

                        const label bfI = pFaces[j];

                        forAllRow(faceEdges, bfI, feI)
                        {
                            const label beI = faceEdges(bfI, feI);

                            const point& ps = points[edges[beI].start()];
                            const point& pe = points[edges[beI].end()];
                            const scalar magE = edges[beI].mag(points) + VSMALL;

                            vector ev = pe - ps;
                            ev /= magE;

                            label nei = edgeFaces(beI, 0);
                            if (nei == bfI)
                                nei = edgeFaces(beI, 1);

                            const label posNei = findIndex(pFaces, nei);
                            if
                            (
                                (posNei < 0) ||
                                (faceInGroup[posNei] != faceInGroup[j])
                            )
                            {
                                if (facePatch_[bfI] != facePatch_[nei])
                                {
                                    DynList<label> patches(2);
                                    patches[0] = facePatch_[bfI];
                                    patches[1] = facePatch_[nei];

                                    // calculate deformation energy
                                    // of the old state
                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mps,
                                        dSqS,
                                        ps,
                                        patches
                                    );
                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mpe,
                                        dSqE,
                                        pe,
                                        patches
                                    );

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Eold +=
                                        1.0/magE*(dSqS + dSqE) + magE*angle;
                                }
                                if (otherPatch != facePatch_[nei])
                                {
                                    // calculate deformation energy
                                    // of the new state
                                    DynList<label> patches(2);
                                    patches[0] = otherPatch;
                                    patches[1] = facePatch_[nei];

                                    point mps, mpe;
                                    scalar dSqS, dSqE;

                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mps,
                                        dSqS,
                                        ps,
                                        patches
                                    );
                                    meshOctree_.findNearestPointToPatches
                                    (
                                        mpe,
                                        dSqE,
                                        pe,
                                        patches
                                    );

                                    vector fv = mpe - mps;
                                    fv /= (mag(fv) + VSMALL);

                                    scalar c = min(fv & ev, 1.0);
                                    c = max(-1.0, c);
                                    const scalar angle = acos(c);

                                    Enew +=
                                        1.0/magE*(dSqS + dSqE) + magE*angle;
                                }
                            }
                        }
                    }

                    if (Enew <= Eold)
                    {
                        newPatchForGroup[1] = otherPatch;

                        forAll(faceInGroup, j)
                        {
                            if (faceInGroup[j] != gp.second())
                                continue;

                            newPatchForFace.append
                            (
                                labelPair(pFaces[j], otherPatch)
                            );
                        }
                    }
                }
            }

        }

        labelHashSet changedPatch;
        for (std::map<label, DynList<labelPair>>::const_iterator it = facesAtCornerNewPatches.cbegin(); it != facesAtCornerNewPatches.cend(); ++it)
        {
            const DynList<labelPair>& lp = it->second;
            forAll(lp, i)
            {
                if (changedPatch.found(lp[i].first()))
                    FatalError << "Face " << lp[i].first()
                               << " is already modified" << abort(FatalError);

                changedPatch.insert(lp[i].first());

                newBoundaryPatches[lp[i].first()] = lp[i].second();
                ++nCorrected;
            }
        }

        if (nCorrected)
        {
            changed = true;
            facePatch_ = newBoundaryPatches;
        }

        if (nIteration++ > 5)
        {
            Info<< "Too many iterations" << endl;
            break;
        }

    } while (nCorrected != 0);

    return changed;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
