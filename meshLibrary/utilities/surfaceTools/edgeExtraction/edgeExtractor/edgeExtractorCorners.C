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
                    ev /= Foam::max(mag(ev), VSMALL);

                    // project the endpoint onto the feature edge
                    point mpe;
                    scalar dSqE;
                    mo.findNearestPointToPatches(mpe, dSqE, pe, patches);

                    vector fv = mpe - mps;
                    fv /= Foam::max(mag(fv), VSMALL);

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
