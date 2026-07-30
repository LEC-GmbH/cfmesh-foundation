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

#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "polyMeshGen2DEngine.H"
#include "meshSurfacePartitioner.H"
#include "detectBoundaryLayers.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Module::refineBoundaryLayers::analyseLayers()
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch = mse.boundaryFacePatches();

    meshSurfacePartitioner mPart(mse);
    detectBoundaryLayers dbl(mPart, is2DMesh_);

    const label nGroups = dbl.nDistinctLayers();
    const labelList& faceInLayer = dbl.faceInLayer();

    // get the hair edges
    splitEdges_ = dbl.hairEdges();

    // create point to split edges addressing
    // Note: use the explicit-nRows overload sized to the CURRENT total
    // mesh point count, not the size-inferring overload (which sizes the
    // graph from the max point label appearing in splitEdges_ alone).
    // refineLayers() can run more than once on the same mesh (e.g. a
    // second "buffer layer" pass in meshOptimizer::optimizeBoundaryLayer),
    // and a later pass's own (smaller, localised) hair-edge set need not
    // reach the highest point labels already present from an earlier
    // pass's generateNewVertices() - undersizing here caused out-of-bounds
    // reads later in refineFace() for faces touching those higher points.
    reverseAddressing(splitEdgesAtPoint_, mesh_.points().size(), splitEdges_);

    // check if the layer is valid
    bool validLayer(true);
    # pragma omp parallel for schedule(dynamic, 40)
    forAll(faceInLayer, bfI)
    {
        if (faceInLayer[bfI] < 0)
        {
            continue;
        }

        const face& bf = bFaces[bfI];

        forAll(bf, pI)
        {
            if (splitEdgesAtPoint_[bf[pI]].size() == 0)
            {
                validLayer = false;
            }
        }
    }

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    // create patch name to index addressing
    HashTable<label> patchNameToIndex;
    forAll(boundaries, patchI)
    {
        patchNameToIndex(boundaries[patchI].patchName()) = patchI;
    }

    // check layer labels over a patch
    layerAtPatch_.setSize(boundaries.size());
    forAll(layerAtPatch_, i)
    {
        layerAtPatch_[i].clear();
    }
    labelListList2D groupsAtPatch(boundaries.size());
    forAll(faceInLayer, bfI)
    {
        appendUniq(groupsAtPatch[facePatch[bfI]], faceInLayer[bfI]);
    }

    // set the information which patches have an extruded layer
    forAll(groupsAtPatch, patchI)
    {
        const DynamicList<label>& layers = groupsAtPatch[patchI];

        forAll(layers, i)
        {
            if (layers[i] < 0)
            {
                layerAtPatch_[patchI].clear();
                break;
            }
            else
            {
                layerAtPatch_[patchI].append(layers[i]);
            }
        }
    }

    // set the information which patches are a single boundary layer face
    patchesInLayer_.setSize(nGroups);
    forAll(layerAtPatch_, patchI)
    {
        const DynamicList<label>& layers = layerAtPatch_[patchI];

        forAll(layers, i)
        {
            patchesInLayer_[layers[i]].append
            (
                boundaries[patchI].patchName()
            );
        }
    }

    // set the number of boundary layers for each patch
    labelList nLayersAtPatch(layerAtPatch_.size(), -1);
    boolList protectedValue(layerAtPatch_.size(), false);

    forAll(patchesInLayer_, layerI)
    {
        const DynamicList<word>& layerPatches = patchesInLayer_[layerI];

        label maxNumLayers(1);
        bool hasLocalValue(false);

        // find the maximum requested number of layers over the layer
        forAll(layerPatches, lpI)
        {
            const word pName = layerPatches[lpI];

            HashTable<label>::const_iterator it =
                numLayersForPatch_.find(pName);

            if (it != numLayersForPatch_.end())
            {
                // check if the layer is interrupted at this patch
                if (discontinuousLayersForPatch_.found(pName))
                {
                    // set the number of layers and lock this location
                    // Note: normalise a disabled (<1) patch to the same
                    // -1 sentinel used everywhere else in this function
                    // (e.g. below and at the non-discontinuous path's
                    // Foam::max(maxNumLayers, ...) which floors at 1) -
                    // without this, a discontinuous patch with nLayers 0
                    // ends up with a literal 0 here, which later gets
                    // multiplied into generateNewCells()'s per-cell
                    // split-count product and zeroes out any cell that
                    // also touches an enabled, multi-layer patch.
                    nLayersAtPatch[patchNameToIndex[pName]] =
                        (it() < 1) ? -1 : it();
                    protectedValue[patchNameToIndex[pName]] = true;
                    hasLocalValue = true;
                }
                else
                {
                    // take the maximum number of layers
                    maxNumLayers = Foam::max(maxNumLayers, it());
                    hasLocalValue = true;
                }
            }
        }

        // apply the global value if no local values exist
        if (!hasLocalValue)
        {
            maxNumLayers = globalNumLayers_;
        }

        // apply the maximum number of ayer of all unprotected patches
        forAll(layerPatches, lpI)
        {
            const label ptchI = patchNameToIndex[layerPatches[lpI]];

            if (!protectedValue[ptchI])
            {
                nLayersAtPatch[ptchI] = maxNumLayers;
            }
        }
    }

    if (is2DMesh_)
    {
        polyMeshGen2DEngine mesh2DEngine(mesh_);
        const boolList& zMinPoint = mesh2DEngine.zMinPoints();
        const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

        const faceList::subList& bFaces = mse.boundaryFaces();

        boolList allZMax(mesh_.boundaries().size(), true);
        boolList allZMin(mesh_.boundaries().size(), true);

        # pragma omp parallel for schedule(dynamic, 50)
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                if (!zMinPoint[bf[pI]])
                {
                    allZMin[facePatch[bfI]] = false;
                }
                if (!zMaxPoint[bf[pI]])
                {
                    allZMax[facePatch[bfI]] = false;
                }
            }
        }

        // mark empty patches as already used
        forAll(allZMin, patchI)
        {
            if (allZMin[patchI] ^ allZMax[patchI])
            {
                nLayersAtPatch[patchI] = -1;
                layerAtPatch_[patchI].clear();
            }
        }
    }

    // set the number of boundary layers which shall be generated above
    // each boundary face
    nLayersAtBndFace_.setSize(facePatch.size());
    nLayersAtBndFace_ = globalNumLayers_;

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(nLayersAtBndFace_, bfI)
    {
        const label patchI = facePatch[bfI];

        if (nLayersAtPatch[patchI] < 0)
        {
            nLayersAtBndFace_[bfI] = 1;
        }
        else
        {
            nLayersAtBndFace_[bfI] = nLayersAtPatch[patchI];

            if (specialMode_)
            {
                ++nLayersAtBndFace_[bfI];
            }
        }
    }

    return validLayer;
}


void Foam::Module::refineBoundaryLayers::generateNewVertices()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelListList2D& pointFaces = mse.pointFaces();
    const labelList& facePatch = mse.boundaryFacePatches();
    const labelList& bp = mse.bp();

    // allocate the data from storing parameters applying to a split edge
    DynamicList<scalar> firstLayerThickness;
    firstLayerThickness.setSize(splitEdges_.size());
    DynamicList<scalar> thicknessRatio;
    thicknessRatio.setSize(splitEdges_.size());
    DynamicList<label> nNodesAtEdge;
    nNodesAtEdge.setSize(splitEdges_.size());
    DynamicList<label> nLayersAtEdge;
    nLayersAtEdge.setSize(splitEdges_.size());

    // count the number of vertices for each split edge
    const label nThreads = 3*omp_get_num_procs();

    # pragma omp parallel num_threads(nThreads)
    {
        // start counting vertices at each thread
        # pragma omp for schedule(static, 1)
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            // get the requested number of boundary layers
            label nLayers(1);
            scalar ratio(globalThicknessRatio_);
            scalar thickness(globalMaxThicknessFirstLayer_);
            bool overridenThickness(false);

            const label bpI = bp[e.start()];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces[bpI][pfI];
                const label pos = help::positionOfEdgeInFace(e, bFaces[bfI]);
                if (pos >= 0)
                {
                    continue;
                }

                const word& patchName =
                    boundaries[facePatch[bfI]].patchName();

                // overrride the global value with the maximum number of layers
                // at this edge
                nLayers = Foam::max(nLayers, nLayersAtBndFace_[bfI]);

                // override with the maximum ratio
                const HashTable<scalar>::const_iterator rIt =
                    thicknessRatioForPatch_.find(patchName);
                if (rIt != thicknessRatioForPatch_.end())
                {
                    ratio = rIt();
                }

                // override with the minimum thickness set for this edge
                const HashTable<scalar>::const_iterator tIt =
                    maxThicknessForPatch_.find(patchName);
                if (tIt != maxThicknessForPatch_.end())
                {
                    if (overridenThickness)
                    {
                        thickness = Foam::min(thickness, tIt());
                    }
                    else
                    {
                        thickness = tIt();
                        overridenThickness = true;
                    }
                }
            }

            // store the information
            firstLayerThickness[seI] = thickness;
            thicknessRatio[seI] = ratio;
            nLayersAtEdge[seI] = nLayers;

            if (!specialMode_)
            {
                nNodesAtEdge[seI] = nLayers + 1;
            }
            else
            {
                nNodesAtEdge[seI] = 3;
            }
        }
    }

    // calculate the number of additional vertices which will be generated
    // on edges of the mesh
    DynamicList<label> numPointsAtThread;
    numPointsAtThread.setSize(nThreads);
    numPointsAtThread = 0;

    # pragma omp parallel for num_threads(nThreads) schedule(static, 1)
    forAll(nNodesAtEdge, seI)
    {
        const label threadI = omp_get_thread_num();

        numPointsAtThread[threadI] += nNodesAtEdge[seI] - 2;
    }

    // allocate the space in a graph storing ids of points on a split edge
    setSizeAndRowSize(newVerticesForSplitEdge_, nNodesAtEdge);

    // calculate the number of points which will be generated
    // on split edges
    label numPoints = points.size();
    forAll(numPointsAtThread, threadI)
    {
        const label nPts = numPointsAtThread[threadI];
        numPointsAtThread[threadI] = numPoints;
        numPoints += nPts;
    }

    points.setSize(numPoints);

    // generate vertices on split edges
    # pragma omp parallel num_threads(nThreads)
    {
        const label threadI = omp_get_thread_num();

        label& nPoints = numPointsAtThread[threadI];

        # pragma omp for schedule(static, 1)
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const vector v = e.vec(points);
            const scalar magv = mag(v);

            const label nLayers = newVerticesForSplitEdge_[seI].size() - 1;

            scalar firstThickness = magv/nLayersAtEdge[seI];
            if (thicknessRatio[seI] >(1. + SMALL))
            {
                firstThickness =
                    magv /
                    (
                        (1 - Foam::pow(thicknessRatio[seI], nLayersAtEdge[seI]))
                        /(1.0 - thicknessRatio[seI])
                    );

            }

            firstThickness =
                Foam::min
                (
                    Foam::max(firstLayerThickness[seI], SMALL),
                    firstThickness
                );

            if (specialMode_)
            {
                scalar t = firstThickness;

                for (label i = 1; i < nLayersAtEdge[seI]-1; ++i)
                {
                    t += firstThickness*Foam::pow(thicknessRatio[seI], i);
                }

                firstThickness = t;
            }

            // generate vertices for this edge
            newVerticesForSplitEdge_[seI][0] = e.start();

            scalar param = firstThickness;
            const vector vec = v /Foam::max(magv, VSMALL);

            for (label pI = 1; pI < nLayers; ++pI)
            {
                // generate the new vertex
                const point newP = points[e.start()] + param*vec;

                param += firstThickness*Foam::pow(thicknessRatio[seI], pI);

                newVerticesForSplitEdge_[seI][pI] = nPoints;
                points[nPoints++] = newP;
            }

            newVerticesForSplitEdge_[seI][nLayers] = e.end();
        }
    }

    if (specialMode_)
    {
        // set the number of layers to 2
        forAll(nLayersAtBndFace_, bfI)
        {
            if (nLayersAtBndFace_[bfI] > 1)
            {
                nLayersAtBndFace_[bfI] = 2;
            }
        }
    }

}


// ************************************************************************* //
