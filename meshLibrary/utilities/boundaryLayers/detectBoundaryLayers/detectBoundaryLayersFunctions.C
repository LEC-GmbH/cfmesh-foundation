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

#include "detectBoundaryLayers.H"
#include "meshSurfacePartitioner.H"
#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGen2DEngine.H"

#include "labelledScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace bndLayerOps
{

class meshBndLayerNeighbourOperator
{
    const meshSurfaceEngine& mse_;
    const label size_;

public:

    meshBndLayerNeighbourOperator(const meshSurfaceEngine& mse)
    :
        mse_(mse),
        size_(mse.boundaryFaces().size())
    {}


    label size() const
    {
        return size_;
    }

    void operator()(const label bfI, DynamicList<label>& neighbourFaces) const
    {
        neighbourFaces.clear();

        const cellListPMG& cells = mse_.cells();

        const labelList& faceOwner = mse_.faceOwners();
        const label cellI = faceOwner[bfI];
        const cell& c = cells[cellI];

        const labelListList2D& faceEdges = mse_.faceEdges();
        const labelListList2D& edgeFaces = mse_.edgeFaces();

        forAllRow(faceEdges, bfI, feI)
        {
            const label edgeI = faceEdges[bfI][feI];

            if (edgeFaces[edgeI].size() == 2)
            {
                label nei = edgeFaces[edgeI][0];

                if (nei == bfI)
                {
                    nei = edgeFaces[edgeI][1];
                }

                // faces must not be part of the same cell
                if (faceOwner[nei] == cellI)
                {
                    continue;
                }

                // owner cell of the other face must
                // have cellI as its neighbour
                const cell& neiC = cells[faceOwner[nei]];
                bool sharedFace(false);
                forAll(c, fI)
                {
                    forAll(neiC, fJ)
                    {
                        if (c[fI] == neiC[fJ])
                        {
                            sharedFace = true;
                            break;
                        }
                    }

                    if (sharedFace)
                    {
                        break;
                    }
                }

                if (sharedFace)
                {
                    neighbourFaces.append(nei);
                }
            }
        }
    }
};


class meshBndLayerSelectorOperator
{
    const meshSurfaceEngine& mse_;

public:

    meshBndLayerSelectorOperator(const meshSurfaceEngine& mse)
    :
        mse_(mse)
    {}


    bool operator()(const label bfI) const
    {
        const labelList& faceOwner = mse_.faceOwners();
        const polyMeshGen& mesh = mse_.mesh();
        const faceListPMG& faces = mesh.faces();

        const cell& c = mesh.cells()[faceOwner[bfI]];
        const PtrList<boundaryPatch>& boundaries = mesh.boundaries();
        const label start = boundaries[0].patchStart();

        label nBndFaces(0), baseFace(-1), otherBase(-1), nQuads(0);
        forAll(c, fI)
        {
            if (faces[c[fI]].size() == 4)
            {
                ++nQuads;
            }

            if ((c[fI] - start) == bfI)
            {
                baseFace = fI;
                ++nBndFaces;
            }
        }

        if (nQuads == 6)
        {
            // cell is a hex
            return true;
        }

        if ((nQuads + 2) != c.size())
        {
            return false;
        }

        if (nBndFaces != 1)
        {
            return false;
        }

        label nQuadsAttachedToBaseFace(0);
        forAll(c, fI)
        {
            if (fI == baseFace)
            {
                continue;
            }

            const bool sEdge =
                help::shareAnEdge(faces[c[baseFace]], faces[c[fI]]);

            if ((faces[c[fI]].size() == 4) && sEdge)
            {
                ++nQuadsAttachedToBaseFace;
            }
            else if (!sEdge)
            {
                if (otherBase != -1)
                {
                    return false;
                }

                otherBase = fI;
            }
        }

        if
        (
            (nQuads == 6)
         || (
                ((nQuadsAttachedToBaseFace + 2) == c.size()) &&
                otherBase != -1 &&
                !help::shareAnEdge(faces[c[baseFace]], faces[c[otherBase]])
            )
        )
        {
            return true;
        }

        return false;
    }
};


} // End namespace bndLayerOps
} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::detectBoundaryLayers::analyseLayers()
{
    Info<< "Analysing mesh for bnd layer existence" << endl;

    const meshSurfaceEngine& mse = meshSurface_.surfaceEngine();
    const polyMeshGen& mesh = mse.mesh();
    const PtrList<boundaryPatch>& boundaries = mesh.boundaries();

    // allocate data needed in parallel loops
    mse.faceOwners();
    mse.faceEdges();
    mse.edgeFaces();
    mse.edges();
    mse.boundaryPointEdges();

    // find layers in patch
    nFirstLayers_ =
        help::groupMarking
        (
            layerAtBndFace_,
            bndLayerOps::meshBndLayerNeighbourOperator(mse),
            bndLayerOps::meshBndLayerSelectorOperator(mse)
        );

    if (is2DMesh_)
    {
        polyMeshGen2DEngine mesh2DEngine(mse.mesh());
        const boolList& zMinPoint = mesh2DEngine.zMinPoints();
        const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

        const faceList::subList& bFaces = mse.boundaryFaces();

        # pragma omp parallel for schedule(dynamic, 50)
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            bool allZMin(true), allZMax(true);
            forAll(bf, pI)
            {
                if (!zMinPoint[bf[pI]])
                {
                    allZMin = false;
                }
                if (!zMaxPoint[bf[pI]])
                {
                    allZMax = false;
                }
            }

            if (allZMax ^ allZMin)
            {
                layerAtBndFace_[bfI] = -1;
            }
        }
    }

    // check the number of layer found at a patch
    typedef std::map<label, DynamicList<label>> patchToLayerType;
    patchToLayerType patchToLayer;

    const labelList& facePatch = meshSurface_.boundaryFacePatches();

    forAll(facePatch, bfI)
    {
        appendUniq(patchToLayer[facePatch[bfI]], layerAtBndFace_[bfI]);
    }

    // all faces of a patch must be in the same layer
    layerAtPatch_.setSize(boundaries.size());
    forAll(layerAtPatch_, i)
    {
        layerAtPatch_[i].clear();
    }

    forAllConstIter(patchToLayerType, patchToLayer, it)
    {
        const DynamicList<label>& layersAtPatch = it->second;

        forAll(layersAtPatch, i)
        {
            if (layersAtPatch[i] < 0)
            {
                layerAtPatch_[it->first].clear();
                break;
            }
            else
            {
                layerAtPatch_[it->first].append(layersAtPatch[i]);
            }
        }
    }

    // set the layer ID to -1 for all faces where the patch is set to -1
    forAll(facePatch, bfI)
    {
        if (layerAtPatch_[facePatch[bfI]].size() == 0)
        {
            layerAtBndFace_[bfI] = -1;
        }
    }

}


bool Foam::Module::detectBoundaryLayers::findHairsForFace
(
    const label bfI,
    DynamicList<edge>& hairEdges
) const
{
    const meshSurfaceEngine& mse = meshSurface_.surfaceEngine();
    const polyMeshGen& mesh = mse.mesh();

    const label nInternalFaces = mesh.boundaries()[0].patchStart();

    const labelList& faceOwner = mse.faceOwners();

    const faceListPMG& faces = mesh.faces();
    const cell& c = mesh.cells()[faceOwner[bfI]];

    // check cell topology
    DynamicList<edge> edges;
    edges.reserve(48);
    DynamicList<DynamicList<label>> edgeFaces;
    labelListList2D faceEdges;
    faceEdges.setSize(c.size());
    label baseFace(-1);
    forAll(c, fI)
    {
        if (c[fI] - nInternalFaces == bfI)
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

    // check does the base face exist and is the number of faces
    // in the cell corresponding to a prism cell
    if ((baseFace < 0) || ((c.size() - faces[c[baseFace]].size()) != 2))
    {
        return false;
    }

    // check if all faces attached to the base face are quads
    bool isPrism(true);

    const face& bf = faces[c[baseFace]];
    hairEdges.setSize(bf.size());

    forAll(bf, pI)
    {
        const label nextEdge = faceEdges[baseFace][pI];
        const label prevEdge = faceEdges[baseFace][bf.rcIndex(pI)];

        if (edgeFaces[nextEdge].size() != 2 || edgeFaces[prevEdge].size() != 2)
        {
            isPrism = false;
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

        label commonEdge;
        for (commonEdge = 0; commonEdge < edges.size(); ++commonEdge)
        {
            if
            (
                (findIndex(edgeFaces[commonEdge], otherNextFace) != -1)
             && (findIndex(edgeFaces[commonEdge], otherPrevFace) != -1)
            )
            {
                break;
            }
        }


        if (commonEdge == edges.size())
        {
            isPrism = false;
            break;
        }

        // there exists a common edge which shall be used as a hair
        if (edges[commonEdge].start() == bf[pI])
        {
            hairEdges[pI] = edges[commonEdge];
        }
        else
        {
            hairEdges[pI] = edges[commonEdge].reverseEdge();
        }
    }

    return isPrism;
}


void Foam::Module::detectBoundaryLayers::generateHairEdges()
{
    hairEdges_.clear();
    hairEdgesAtBoundaryPoint_.clear();

    const meshSurfaceEngine& mse = meshSurface_.surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    mse.faceOwners();
    const labelListList2D& pFaces = mse.pointFaces();
    const labelList& bp = mse.bp();

    # pragma omp parallel
    {
        // per-thread local buffer; findHairsForFace() is a purely local
        // query (verified: all its intermediate state is declared inside
        // the function, only reads shared mesh data through const refs,
        // writes only to its own out-parameter) so this loop is safe to
        // parallelize as long as the shared hairEdges_ append below is
        // properly guarded
        DynamicList<edge> localHairEdges;

        # pragma omp for schedule(dynamic, 40)
        forAll(bFaces, bfI)
        {
            if (layerAtBndFace_[bfI] < 0)
            {
                continue;
            }

            // find hair edges for this face
            DynamicList<edge> hairEdges;
            if (!findHairsForFace(bfI, hairEdges))
            {
                continue;
            }

            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                // store hair edges in a list
                const edge& he = hairEdges[pI];

                if (he.start() != bf[pI])
                {
                    FatalErrorInFunction
                        << "Wrong starting point" << abort(FatalError);
                }

                localHairEdges.append(he);
            }
        }

        // reserve a disjoint range in the shared list for this thread's
        // edges (mirrors the correct disjoint-range-reservation pattern
        // used elsewhere in this codebase, e.g.
        // polyMeshGenModifierReorderBoundaryFaces.C), then write into it
        // only after every thread has finished reserving
        label startEl;
        # pragma omp critical
        {
            startEl = hairEdges_.size();
            hairEdges_.setSize(startEl + localHairEdges.size());
        }

        # pragma omp barrier

        forAll(localHairEdges, i)
        {
            hairEdges_[startEl + i] = localHairEdges[i];
        }
    }

    // filter out duplicate edges
    labelListList2D pHairEdges;
    reverseAddressing(pHairEdges, hairEdges_);

    boolList duplicateEdge(hairEdges_.size(), false);
    # pragma omp parallel for schedule(dynamic, 50)
    forAll(pHairEdges, pointI)
    {
        forAllRow(pHairEdges, pointI, pheI)
        {
            const label heI = pHairEdges[pointI][pheI];
            const edge& he = hairEdges_[heI];

            for
            (
                label pheJ = pheI + 1;
                pheJ < pHairEdges[pointI].size();
                ++pheJ
            )
            {
                const label heJ = pHairEdges[pointI][pheJ];
                const edge& nhe = hairEdges_[heJ];

                if (he == nhe)
                {
                    duplicateEdge[heJ] = true;
                }
            }
        }
    }

    label counter(0);
    forAll(hairEdges_, heI)
    {
        if (!duplicateEdge[heI])
        {
            if (heI > counter)
            {
                hairEdges_[counter++] = hairEdges_[heI];
            }
            else
            {
                ++counter;
            }
        }
    }

    hairEdges_.setSize(counter);

    // create point to split edges addressing
    hairEdgesAtBoundaryPoint_.setSize(pFaces.size());

    forAll(hairEdges_, heI)
    {
        const edge& he = hairEdges_[heI];
        hairEdgesAtBoundaryPoint_[bp[he.start()]].append(heI);
    }
}


// ************************************************************************* //
