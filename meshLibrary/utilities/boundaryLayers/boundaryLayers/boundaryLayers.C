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

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfacePartitioner.H"
#include "polyMeshGen2DEngine.H"

#include "labelledPoint.H"
#include "EdgeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Module::meshSurfaceEngine&
Foam::Module::boundaryLayers::surfaceEngine() const
{
    if (!msePtr_)
    {
        msePtr_ = new meshSurfaceEngine(mesh_);
    }

    return *msePtr_;
}


const Foam::Module::meshSurfacePartitioner& 
Foam::Module::boundaryLayers::surfacePartitioner() const
{
    if (!meshPartitionerPtr_)
    {
        meshPartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine());
    }

    return *meshPartitionerPtr_;
}


void Foam::Module::boundaryLayers::findPatchesToBeTreatedTogether()
{
    if (geometryAnalysed_)
    {
        return;
    }

    forAll(treatPatchesWithPatch_, patchI)
    {
        treatPatchesWithPatch_[patchI].append(patchI);
    }

    const meshSurfaceEngine& mse = surfaceEngine();

    const pointFieldPMG& points = mesh_.points();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const edgeList& edges = mse.edges();
    const labelListList2D& eFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const labelListList2D& pPatches = mPart.pointPatches();

    // patches must be treated together if there exist a corner where
    // more than three patches meet
    const labelHashSet& corners = mPart.corners();
    for (const label bpI : corners.toc())
    {
        if (mPart.numberOfFeatureEdgesAtPoint(bpI) > 3)
        {
            labelHashSet commonPatches;
            DynamicList<label> allPatches;

            forAllRow(pPatches, bpI, patchI)
            {
                const DynamicList<label>& tpwp =
                    treatPatchesWithPatch_[pPatches[bpI][patchI]];

                forAll(tpwp, pJ)
                {
                    if (commonPatches.found(tpwp[pJ]))
                    {
                        continue;
                    }

                    commonPatches.insert(tpwp[pJ]);
                    allPatches.append(tpwp[pJ]);
                }
            }

            forAllRow(pPatches, bpI, patchI)
            {
                treatPatchesWithPatch_[pPatches[bpI][patchI]] = allPatches;
            }

        }
    }

    // patches must be treated together for concave geometries
    // edgeClassification map counts the number of convex and concave edges
    // for a given patch. The first counts convex edges and the second counts
    // concave ones. If the number of concave edges is of the considerable
    // percentage, it is treated as O-topology
    meshSurfaceCheckInvertedVertices vertexCheck(mse);
    const labelHashSet& invertedVertices = vertexCheck.invertedVertices();

    EdgeMap<Pair<label>> edgeClassification;
    forAll(eFaces, eI)
    {
        if (eFaces[eI].size() != 2)
        {
            continue;
        }

        // check if the any of the face vertices is tangled
        const edge& e = edges[eI];
        if
        (
            !is2DMesh_
         && (invertedVertices.found(e[0]) || invertedVertices.found(e[1]))
        )
        {
            continue;
        }

        const label patch0 = boundaryFacePatches[eFaces[eI][0]];
        const label patch1 = boundaryFacePatches[eFaces[eI][1]];
        if (patch0 != patch1)
        {
            edge pp
            (
                Foam::min(patch0, patch1),
                Foam::max(patch0, patch1)
            );
            if (!edgeClassification.found(pp))
            {
                edgeClassification.insert
                (
                    pp, Pair<label>(0, 0)
                );
            }

            const face& f1 = bFaces[eFaces[eI][0]];
            const face& f2 = bFaces[eFaces[eI][1]];

            if
            (
                !help::isSharedEdgeConvex(points, f1, f2)
             || (help::angleBetweenFaces(points, f1, f2) > 0.75*M_PI)
            )
            {
                ++edgeClassification[pp].second();
            }
            else
            {
                ++edgeClassification[pp].first();
            }
        }
    }

    forAllConstIter(EdgeMap<Pair<label>>, edgeClassification, it)
    {
        const edge& edgePair = it.key();
        const Pair<label>& nConvexAndConcave = it();

        if (nConvexAndConcave.second() != 0)
        {
            // number of concave edges is greater than the number
            // of the convex ones. Treat patches together.
            const label patch0 = edgePair.start();
            const label patch1 = edgePair.end();

            // avoid adding unused patches in case of 2D meshing
            if (treatedPatch_[patch0] || treatedPatch_[patch1])
            {
                continue;
            }

            treatPatchesWithPatch_[patch0].append(patch1);
            treatPatchesWithPatch_[patch1].append(patch0);
        }
    }

    // final adjusting of patches which shall be treated together
    boolList confirmed(treatPatchesWithPatch_.size(), false);
    forAll(treatPatchesWithPatch_, patchI)
    {
        if (treatPatchesWithPatch_[patchI].size() <= 1)
        {
            confirmed[patchI] = true;
            continue;
        }

        if (confirmed[patchI])
        {
            continue;
        }

        labelHashSet commonPatches;
        commonPatches.insert(patchI);

        DynamicList<label> front;
        front.append(patchI);
        confirmed[patchI] = true;

        while (front.size())
        {
            const label fPatch = front.remove();

            forAll(treatPatchesWithPatch_[fPatch], i)
            {
                const label patchJ = treatPatchesWithPatch_[fPatch][i];

                if (confirmed[patchJ])
                {
                    continue;
                }

                front.append(patchJ);
                confirmed[patchJ] = true;
                commonPatches.insert(patchJ);
                forAll(treatPatchesWithPatch_[patchJ], j)
                {
                    commonPatches.insert(treatPatchesWithPatch_[patchJ][j]);
                }
            }
        }

        const labelList commonPatchesSorted(commonPatches.sortedToc());
        forAll(commonPatchesSorted, i)
        {
            const label patchJ = commonPatchesSorted[i];

            treatPatchesWithPatch_[patchJ].clear();
            forAll(commonPatchesSorted, j)
            {
                treatPatchesWithPatch_[patchJ].append(commonPatchesSorted[j]);
            }
        }
    }

    geometryAnalysed_ = true;
}


void Foam::Module::boundaryLayers::addLayerForPatch(const label patchLabel)
{
    if (treatedPatch_[patchLabel])
    {
        return;
    }

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    if (boundaries[patchLabel].patchSize() == 0)
    {
        return;
    }

    boolList treatPatches(boundaries.size(), false);
    if (patchWiseLayers_)
    {
        forAll(treatPatchesWithPatch_[patchLabel], pI)
        {
            treatPatches[treatPatchesWithPatch_[patchLabel][pI]] = true;
        }
    }
    else
    {
        forAll(treatedPatch_, patchI)
        {
            if (!treatedPatch_[patchI])
            {
                treatPatches[patchI] = true;
            }
        }
    }

    newLabelForVertex_.setSize(nPoints_);
    newLabelForVertex_ = -1;
    otherVrts_.clear();
    patchKey_.clear();

    createNewVertices(treatPatches);

    createNewFacesAndCells(treatPatches);

    forAll(treatPatches, patchI)
    {
        if (treatPatches[patchI])
        {
            treatedPatch_[patchI] = true;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::boundaryLayers::boundaryLayers
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    msePtr_(nullptr),
    meshPartitionerPtr_(nullptr),
    patchWiseLayers_(true),
    is2DMesh_(false),
    patchNames_(),
    patchTypes_(),
    treatedPatch_(),
    treatPatchesWithPatch_(),
    newLabelForVertex_(),
    otherVrts_(),
    patchKey_(),
    nPoints_(mesh.points().size()),
    geometryAnalysed_(false)
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    patchNames_.setSize(boundaries.size());
    patchTypes_.setSize(boundaries.size());
    forAll(boundaries, patchI)
    {
        patchNames_[patchI] = boundaries[patchI].patchName();
        patchTypes_[patchI] = boundaries[patchI].patchType();
    }

    treatedPatch_.setSize(boundaries.size());
    treatedPatch_ = false;

    treatPatchesWithPatch_.setSize(boundaries.size());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::boundaryLayers::~boundaryLayers()
{
    clearOut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::boundaryLayers::addLayerForPatch(const word& patchName)
{
    if (!geometryAnalysed_)
    {
        findPatchesToBeTreatedTogether();
    }

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    forAll(boundaries, patchI)
    {
        if (boundaries[patchI].patchName() == patchName)
        {
            addLayerForPatch(patchI);
        }
    }
}


void Foam::Module::boundaryLayers::activate2DMode()
{
    polyMeshGen2DEngine mesh2DEngine(mesh_);
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

    const faceList::subList& bFaces = surfaceEngine().boundaryFaces();
    const labelList& facePatch = surfaceEngine().boundaryFacePatches();

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
            treatedPatch_[patchI] = true;
        }
    }

    forAll(treatPatchesWithPatch_, patchI)
    {
        DynamicList<label>& patches = treatPatchesWithPatch_[patchI];

        for (label i = patches.size()-1; i>=0; --i)
        {
            if (treatedPatch_[patches[i]])
            {
                patches.erase(patches.begin() + i);
            }
        }
    }

    is2DMesh_ = true;
}


void Foam::Module::boundaryLayers::addLayerForAllPatches()
{
    if (!geometryAnalysed_)
    {
        findPatchesToBeTreatedTogether();
    }

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    if (!patchWiseLayers_)
    {
        forAll(boundaries, patchI)
        {
            addLayerForPatch(patchI);
        }
    }
    else
    {
        newLabelForVertex_.setSize(nPoints_);
        newLabelForVertex_ = -1;
        otherVrts_.clear();
        patchKey_.clear();

        // avoid generating bnd layer at empty patches in case of 2D meshing
        label counter(0);
        forAll(treatedPatch_, patchI)
        {
            if (!treatedPatch_[patchI])
            {
                ++counter;
            }
        }

        labelList treatedPatches(counter);
        counter = 0;
        forAll(treatedPatch_, i)
        {
            if (!treatedPatch_[i])
            {
                treatedPatches[counter++] = i;
            }
        }

        // create bnd layer vertices
        createNewVertices(treatedPatches);

        // create bnd layer cells
        createLayerCells(treatedPatches);
    }
}


// ************************************************************************* //
