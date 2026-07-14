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
#include "meshSurfaceMapper.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfaceCheckEdgeTypes.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Module::edgeExtractor::calculateValence()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    pointValence_.setSize(mse.boundaryPoints().size());
    pointValence_ = 0;

    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();

    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
            ++pointValence_[bp[bf[pI]]];
    }

}


void Foam::Module::edgeExtractor::calculateSingleCellEdge()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const edgeList& edges = mse.edges();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& faceCells = mse.faceOwners();

    // find the number of boundary faces for each cell in the mesh
    edgeType_.setSize(edgeFaces.size());
    edgeType_ = NONE;

    forAll(edgeFaces, eI)
    {
        if (edgeFaces.sizeOfRow(eI) == 2)
        {
            const label c0 = faceCells[edgeFaces(eI, 0)];
            const label c1 = faceCells[edgeFaces(eI, 1)];

            if (c0 == c1)
                edgeType_[eI] |= SINGLECELLEDGE;
        }
    }

    // calculate the number of cells attache to a boundary edge
    const labelList& bp = mse.bp();
    const cellListPMG& cells = mse.mesh().cells();
    const faceListPMG& faces = mse.faces();

    nCellsAtEdge_.setSize(edgeFaces.size());
    nCellsAtEdge_ = 0;

    # pragma omp parallel for schedule(dynamic, 100)
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];

        DynList<edge> foundEdge;

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, eI)
            {
                const edge e = f.faceEdge(eI);

                const label bps = bp[e.start()];

                if (bps < 0)
                    continue;

                forAllRow(bpEdges, bps, i)
                {
                    const label beI = bpEdges(bps, i);
                    const edge& be = edges[beI];

                    if ((e == be) && findIndex(foundEdge, be) == -1)
                    {
                        foundEdge.append(be);

                        # pragma omp atomic
                        ++nCellsAtEdge_[beI];
                    }
                }
            }
        }
    }
}


void Foam::Module::edgeExtractor::findPatchesNearSurfaceFace()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const triSurf& surface = meshOctree_.surface();

    patchesNearFace_.setSize(bFaces.size());
    labelLongList nPatchesAtFace(bFaces.size());

    # pragma omp parallel
    {
        labelLongList localData;
        DynList<label> nearFacets;

        # pragma omp for schedule(dynamic, 40)
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            const vector c = help::stableFaceCentre(bf, points);

            // find a reasonable searching distance comparable to face size
            scalar d(0.0);
            forAll(bf, pI)
                d = Foam::max(d, Foam::mag(c - points[bf[pI]]));
            d = 2.0*d + VSMALL;

            const boundBox bb(c - vector(d, d, d), c + vector(d, d, d));

            // get the patches near the current boundary face
            meshOctree_.findTrianglesInBox(bb, nearFacets);
            DynList<label> nearPatches;
            forAll(nearFacets, i)
                nearPatches.appendUniq(surface[nearFacets[i]].region());

            localData.append(bfI);
            nPatchesAtFace[bfI] = nearPatches.size();
            forAll(nearPatches, i)
                localData.append(nearPatches[i]);
        }

        # pragma omp barrier

        # pragma omp master
        patchesNearFace_.setSizeAndRowSize(nPatchesAtFace);

        # pragma omp barrier

        // copy the data to the graph
        label counter(0);
        while (counter < localData.size())
        {
            const label edgeI = localData[counter++];

            const label size = nPatchesAtFace[edgeI];

            for (label i = 0; i < size; ++i)
                patchesNearFace_(edgeI, i) = localData[counter++];
        }
    }
}


void Foam::Module::edgeExtractor::findFeatureEdgesNearEdge()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const edgeList& edges = mse.edges();

    featureEdgesNearEdge_.setSize(edges.size());
    labelLongList nFeatureEdgesAtEdge(edges.size());

    # pragma omp parallel
    {
        labelLongList localData;
        DynList<label> nearEdges;

        # pragma omp for schedule(dynamic, 40)
        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];
            const vector c = e.centre(points);
            const scalar d = 1.5*e.mag(points);

            const boundBox bb(c - vector(d, d, d), c + vector(d, d, d));

            // get the edges near the current edge
            meshOctree_.findEdgesInBox(bb, nearEdges);
            forAllReverse(nearEdges, i)
            {
                const label pos = findIndex(nearEdges, nearEdges[i]);

                if (pos < i)
                    nearEdges.removeElement(i);
            }

            localData.append(edgeI);
            nFeatureEdgesAtEdge[edgeI] = nearEdges.size();
            forAll(nearEdges, i)
                localData.append(nearEdges[i]);
        }

        # pragma omp barrier

        # pragma omp master
        featureEdgesNearEdge_.setSizeAndRowSize(nFeatureEdgesAtEdge);

        # pragma omp barrier

        // copy the data to the graph
        label counter(0);
        while (counter < localData.size())
        {
            const label edgeI = localData[counter++];

            const label size = nFeatureEdgesAtEdge[edgeI];

            for (label i = 0; i < size; ++i)
                featureEdgesNearEdge_(edgeI, i) = localData[counter++];
        }
    }
}


void Foam::Module::edgeExtractor::markPatchPoints(boolList& patchPoint)
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& bp = mse.bp();

    patchPoint.setSize(bPoints.size());
    patchPoint = true;

    std::map<label, label> otherProcPatch;

    // set the patchPoint to false for all vertices at feature edges
    # pragma omp parallel for schedule(dynamic, 40)
    forAll(edgeFaces, beI)
    {
        if (edgeFaces.sizeOfRow(beI) == 2)
        {
            // an ordinary edge
            if (facePatch_[edgeFaces(beI, 0)] != facePatch_[edgeFaces(beI, 1)])
            {
                const edge& e = edges[beI];
                patchPoint[bp[e.start()]] = false;
                patchPoint[bp[e.end()]] = false;
            }
        }
        else if (edgeFaces.sizeOfRow(beI) == 1)
        {
            // an edge at a parallel interface
            const label otherPatch = otherProcPatch[beI];

            if (facePatch_[edgeFaces(beI, 0)] != otherPatch)
            {
                const edge& e = edges[beI];
                patchPoint[bp[e.start()]] = false;
                patchPoint[bp[e.end()]] = false;
            }
        }
        else
        {
            // this is a non-manifold edge
            const edge& e = edges[beI];
            patchPoint[bp[e.start()]] = false;
            patchPoint[bp[e.end()]] = false;
        }
    }
}


const Foam::Module::meshSurfaceEngine&
Foam::Module::edgeExtractor::surfaceEngine() const
{
    // the outer unprotected null-check that used to sit here was removed:
    // a thread seeing a non-null pointer there would have no guarantee,
    // under the C++ memory model, of observing a fully-constructed object
    // (this was previously safe only by relying on the first call always
    // happening serially, before any parallel region - true today, but
    // fragile). #pragma omp critical implies a flush on both entry and
    // exit, so always going through it - even for the fast "already
    // constructed" path - makes every call properly synchronized.
    # pragma omp critical
    {
        if (!surfaceEnginePtr_)
        {
            surfaceEnginePtr_ = new meshSurfaceEngine(mesh_);
        }
    }

    return *surfaceEnginePtr_;
}


const Foam::Module::triSurfacePartitioner&
Foam::Module::edgeExtractor::partitioner() const
{
    // see the comment in surfaceEngine() above - same fix, same reasoning
    # pragma omp critical
    {
        if (!surfPartitionerPtr_)
        {
            surfPartitionerPtr_ =
                new triSurfacePartitioner(meshOctree_.surface());
        }
    }

    return *surfPartitionerPtr_;
}


const Foam::Module::triSurfaceClassifyEdges&
Foam::Module::edgeExtractor::edgeClassifier() const
{
    // previously had no synchronization at all (a weaker version of the
    // same lazy-singleton anti-pattern fixed in surfaceEngine()/
    // partitioner() above); added the same critical-guarded check
    # pragma omp critical
    {
        if (!surfEdgeClassificationPtr_)
        {
            surfEdgeClassificationPtr_ =
                new triSurfaceClassifyEdges(meshOctree_);
        }
    }

    return *surfEdgeClassificationPtr_;
}


void Foam::Module::edgeExtractor::findFaceCandidates
(
    labelLongList& faceCandidates,
    const labelList* facePatchPtr,
    const Map<label>* otherFacePatchPtr
) const
{
    faceCandidates.clear();
    if (!facePatchPtr)
        facePatchPtr = &facePatch_;

    const labelList& fPatches = *facePatchPtr;

    bool deleteOtherFacePatchPtr(false);
    if (!otherFacePatchPtr)
    {
        otherFacePatchPtr = new Map<label>();
        deleteOtherFacePatchPtr = true;
    }

    const Map<label>& otherFacePatch = *otherFacePatchPtr;

    const meshSurfaceEngine& mse = surfaceEngine();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    # pragma omp parallel if (faceEdges.size() > 1000)
    {
        labelLongList procCandidates;
        # pragma omp for schedule(dynamic, 40)
        forAll(faceEdges, bfI)
        {
            DynList<label> allNeiPatches;
            forAllRow(faceEdges, bfI, feI)
            {
                const label beI = faceEdges(bfI, feI);

                if (edgeFaces.sizeOfRow(beI) == 2)
                {
                    label fNei = edgeFaces(beI, 0);
                    if (fNei == bfI)
                        fNei = edgeFaces(beI, 1);

                    allNeiPatches.appendUniq(fPatches[fNei]);
                }
                else if (edgeFaces.sizeOfRow(beI) == 1)
                {
                    allNeiPatches.appendUniq(otherFacePatch[beI]);
                }
            }

            if (allNeiPatches.size() > 1)
            {
                // this face is probably near some feature edge
                procCandidates.append(bfI);
            }
        }

        # pragma omp critical
        {
            forAll(procCandidates, i)
                faceCandidates.append(procCandidates[i]);
        }
    }

    if (deleteOtherFacePatchPtr)
        deleteDemandDrivenData(otherFacePatchPtr);
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

Foam::Module::edgeExtractor::edgeExtractor
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    surfaceEnginePtr_(nullptr),
    meshOctree_(octree),
    surfPartitionerPtr_(nullptr),
    surfEdgeClassificationPtr_(nullptr),
    pointValence_(),
    pointPatch_(),
    facePatch_(),
    nCellsAtEdge_(),
    edgeType_(),
    patchesNearFace_(),
    featureEdgesNearEdge_()
{
    calculateValence();

    calculateSingleCellEdge();

    findPatchesNearSurfaceFace();

    findFeatureEdgesNearEdge();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
// Destructor

Foam::Module::edgeExtractor::~edgeExtractor()
{
    deleteDemandDrivenData(surfaceEnginePtr_);
    deleteDemandDrivenData(surfPartitionerPtr_);
    deleteDemandDrivenData(surfEdgeClassificationPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void Foam::Module::edgeExtractor::moveVerticesTowardsDiscontinuities
(
    const label nIterations
)
{
    Info<< "Reducing Hausdorff distance:" << flush;

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointFaces = mse.pointFaces();
    const pointFieldPMG& points = mse.points();
    const faceList::subList& bFaces = mse.boundaryFaces();

    meshSurfaceEngineModifier modifier(mse);

    vectorField faceCentreDisplacement(bFaces.size());
    List<labelledPoint> pointDisplacements(bPoints.size());

    for (label iterI = 0; iterI < nIterations; ++iterI)
    {
        # pragma omp parallel
        {
            // find displacements of face centres
            # pragma omp for schedule(dynamic, 40)
            forAll(bFaces, bfI)
            {
                const vector centre = help::stableFaceCentre(bFaces[bfI], points);

                point newP;
                scalar distSq;
                label patchI, nearestTri;
                meshOctree_.findNearestSurfacePoint
                (
                    newP,
                    distSq,
                    nearestTri,
                    patchI,
                    centre
                );

                faceCentreDisplacement[bfI] = newP - centre;
            }

            // initialise displacements
            # pragma omp for schedule(dynamic, 40)
            forAll(pointDisplacements, bpI)
                pointDisplacements[bpI] = labelledPoint(0, vector::zero);

            # pragma omp barrier

            // calculate displacements of boundary points as the average
            // of face centre displacements
            # pragma omp for schedule(dynamic, 40)
            forAll(pointFaces, bpI)
            {
                forAllRow(pointFaces, bpI, pfI)
                {
                    pointDisplacements[bpI].coordinates() +=
                        faceCentreDisplacement[pointFaces(bpI, pfI)];
                    ++pointDisplacements[bpI].pointLabel();
                }
            }
        }

        # pragma omp parallel for schedule(dynamic, 40)
        forAll(pointDisplacements, bpI)
        {
            const labelledPoint& lp = pointDisplacements[bpI];
            const point mp =
                points[bPoints[bpI]] + lp.coordinates() / lp.pointLabel();

            //Info<< "Original point " << bPoints[bpI] << " has coordinates "
            //        << points[bPoints[bpI]] << endl;
            //Info<< "Displacement vector " << lp.coordinates() / lp.pointLabel() << endl;
            //Info<< "Moved point " << mp << endl;

            point newPoint;
            label patchI, nt;
            scalar distSq;

            meshOctree_.findNearestSurfacePoint(newPoint, distSq, nt, patchI, mp);

            //Info<< "Mapped point " << newPoint << nl << endl;

            modifier.moveBoundaryVertexNoUpdate(bpI, newPoint);
        }

        // update geometry
        modifier.updateGeometry();

        Info<< '.' << flush;
    }

    Info<< endl;
}


void Foam::Module::edgeExtractor::distributeBoundaryFaces()
{
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const pointFieldPMG& points = mse.points();

    // set the size of the facePatch list
    facePatch_.setSize(bFaces.size());

    // check if the mesh already has patches
    if (mesh_.boundaries().size() > 1)
        Warning << "Mesh patches are already assigned!" << endl;

    // set size of patchNames, newBoundaryFaces_ and newBoundaryOwners_
    const triSurf& surface = meshOctree_.surface();
    const label nPatches = surface.patches().size();

    // find patches to which the surface points are mapped to
    pointPatch_.setSize(bPoints.size());

    # pragma omp parallel for schedule(dynamic, 40)
    forAll(bPoints, bpI)
    {
        const point& bp = points[bPoints[bpI]];

        label fPatch, nTri;
        point p;
        scalar distSq;

        meshOctree_.findNearestSurfacePoint(p, distSq, nTri, fPatch, bp);

        if ((fPatch > -1) && (fPatch < nPatches))
        {
            pointPatch_[bpI] = fPatch;
        }
        else
        {
            pointPatch_[bpI] = nPatches;
            FatalErrorInFunction
                << "Cannot distribute a boundary points " << bPoints[bpI]
                << " into any surface patch!. Exiting.." << exit(FatalError);
        }
    }

    // find the patch for face by finding the patch nearest
    // to the face centre
    # pragma omp parallel for schedule(dynamic, 40)
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        const point c = help::stableFaceCentre(bf, points);

        // find the nearest surface patch to face centre
        label fPatch, nTri;
        point p;
        scalar distSq;

        meshOctree_.findNearestSurfacePoint(p, distSq, nTri, fPatch, c);

        if ((fPatch > -1) && (fPatch < nPatches))
        {
            facePatch_[bfI] = fPatch;
        }
        else
        {
            facePatch_[bfI] = nPatches;

            FatalErrorInFunction
                << "Cannot distribute a face " << bFaces[bfI] << " into any "
                << "surface patch!. Exiting.." << exit(FatalError);
        }
    }
}


bool Foam::Module::edgeExtractor::distributeBoundaryFacesNormalAlignment()
{
    bool changed(false);

    const pointFieldPMG& points = mesh_.points();
    const meshSurfaceEngine& mse = this->surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    const triSurf& surf = meshOctree_.surface();
    const pointField& sPoints = surf.points();

    label nCorrected, nIterations(0);
    Map<label> otherProcNewPatch;

    do
    {
        nCorrected = 0;

        // allocate a copy of boundary patches
        labelList newBoundaryPatches(facePatch_);

        // find the faces which have neighbouring faces in other patches
        labelLongList candidates;
        findFaceCandidates(candidates, &facePatch_, &otherProcNewPatch);

        // go through the list of faces and check if they shall remain
        // in the current patch
        # pragma omp parallel for schedule(dynamic, 40) \
        reduction(+ : nCorrected)
        forAll(candidates, i)
        {
            const label bfI = candidates[i];
            const face& bf = bFaces[bfI];

            DynList<label> allNeiPatches;
            DynList<label> neiPatches;
            neiPatches.setSize(faceEdges.sizeOfRow(bfI));

            forAllRow(faceEdges, bfI, eI)
            {
                const label beI = faceEdges(bfI, eI);

                if (edgeFaces.sizeOfRow(beI) == 2)
                {
                    label fNei = edgeFaces(beI, 0);
                    if (fNei == bfI)
                        fNei = edgeFaces(faceEdges(bfI, eI), 1);

                    allNeiPatches.appendUniq(facePatch_[fNei]);
                    neiPatches[eI] = facePatch_[fNei];
                }
                else if (edgeFaces.sizeOfRow(beI) == 1)
                {
                    allNeiPatches.appendUniq(otherProcNewPatch[beI]);
                    neiPatches[eI] = otherProcNewPatch[beI];
                }
            }

            // do not modify faces with all neighbours in the same patch
            if
            (
                (allNeiPatches.size() == 1) &&
                (allNeiPatches[0] == facePatch_[bfI])
            )
                continue;

            // check whether there exist edges which are more suitable for
            // projection onto feature edges than the currently selected ones
            label newPatch(-1);
            DynList<scalar> normalAlignment(allNeiPatches.size());
            DynList<scalar> distanceSq(allNeiPatches.size());
            scalar maxDSq(0.0);
            forAll(allNeiPatches, i)
            {
                point pMap;
                scalar dSq(VGREAT);
                label nearestTriangle;

                point p = help::stableFaceCentre(bf, points);
                meshOctree_.findNearestSurfacePointInRegion
                (
                    pMap,
                    dSq,
                    nearestTriangle,
                    allNeiPatches[i],
                    p
                );

                maxDSq = Foam::max(dSq, maxDSq);

                // calculate normal vectors
                vector tn = surf[nearestTriangle].normal(sPoints);
                vector fn = help::stableFaceNormal(bf, points);

                // calculate alignment
                normalAlignment[i] = mag(tn & fn);
                distanceSq[i] = dSq;
            }

            scalar maxAlignment(0.0);
            forAll(normalAlignment, i)
            {
                const scalar metric
                (
                    sqrt
                    (
                        maxDSq
                       /Foam::max
                        (
                            Foam::max(distanceSq[i], SMALL*maxDSq),
                            VSMALL
                        )
                    )
                   *normalAlignment[i]
                );

                if (metric > maxAlignment)
                {
                    maxAlignment = metric;
                    newPatch = allNeiPatches[i];
                }
            }

            if ((newPatch >= 0) && (newPatch != facePatch_[bfI]))
            {
                newBoundaryPatches[bfI] = newPatch;
                ++nCorrected;
            }
        }


        if (nCorrected)
        {
            changed = true;

            // transfer the new patches back
            facePatch_.transfer(newBoundaryPatches);
        }
    } while ((nCorrected != 0) && (++nIterations < 5));

    return changed;
}


void Foam::Module::edgeExtractor::findEdgeCandidates()
{
    const triSurf& surface = meshOctree_.surface();
    const vectorField& sp = surface.points();
    const VRWGraph& facetEdges = surface.facetEdges();
    const VRWGraph& edgeFacets = surface.edgeFacets();

    const triSurfacePartitioner& partitioner = this->partitioner();

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const VRWGraph& faceEdges = mse.faceEdges();

    Map<label> otherFacePatch;
    labelLongList faceCandidates;
    findFaceCandidates(faceCandidates, &facePatch_, &otherFacePatch);

    # pragma omp parallel for schedule(dynamic, 40) \
    if (faceCandidates.size() > 100)
    forAll(faceCandidates, fcI)
    {
        const label bfI = faceCandidates[fcI];

        forAllRow(faceEdges, bfI, i)
        {
            const label eI = faceEdges(bfI, i);
            edgeType_[eI] |= CANDIDATE;
        }
    }

    // find distances of all vertices supporting CANDIDATE edges
    // from feature edges separating various patches
    const VRWGraph& pEdges = mse.boundaryPointEdges();
    const edgeList& edges = mse.edges();

    List<List<labelledScalar>> featureEdgesNearPoint(bPoints.size());

    DynList<label> containedTriangles;
    # pragma omp parallel for schedule(dynamic, 40) private(containedTriangles)
    forAll(pEdges, bpI)
    {
        // TODO rewrite for execution on distributed machines
        bool check(false);
        forAllRow(pEdges, bpI, peI)
        {
            const label eI = pEdges(bpI, peI);

            if (edgeType_[eI] & CANDIDATE)
            {
                check = true;
                break;
            }
        }

        if (check)
        {
            // check the squared distance from the nearest feature edge
            scalar rSq(0.0);
            forAllRow(pEdges, bpI, peI)
            {
                const label eI = pEdges(bpI, peI);
                const edge& e = edges[eI];
                const scalar dSq = magSqr(points[e.start()] - points[e.end()]);

                rSq = Foam::max(rSq, dSq);
            }

            rSq *= 2.0;
            const scalar r = Foam::sqrt(rSq);

            // create a boundaing box used for searching neighbour edges
            const point& p = points[bPoints[bpI]];
            boundBox bb(p - point(r, r, r), p + point(r, r, r));

            // find the surface triangles in the vicinity of the point
            // check for potential feature edges
            containedTriangles.clear();
            meshOctree_.findTrianglesInBox(bb, containedTriangles);

            DynList<label> featureEdgeCandidates;

            forAll(containedTriangles, ctI)
            {
                const label tI = containedTriangles[ctI];

                forAllRow(facetEdges, tI, feI)
                {
                    const label seI = facetEdges(tI, feI);

                    if (edgeFacets.sizeOfRow(seI) == 2)
                    {
                        const label p0 = surface[edgeFacets(seI, 0)].region();
                        const label p1 = surface[edgeFacets(seI, 1)].region();

                        if (p0 != p1)
                        {
                            featureEdgeCandidates.appendUniq(seI);
                        }
                    }
                    else
                    {
                        featureEdgeCandidates.appendUniq(seI);
                    }
                }
            }

            // check the distance of the vertex from the candidates
            List<labelledScalar>& featureEdgeDistances =
                featureEdgesNearPoint[bpI];
            featureEdgeDistances.setSize(featureEdgeCandidates.size());
            forAll(featureEdgeCandidates, i)
            {
                const label seI = featureEdgeCandidates[i];

                const point s = sp[edges[seI].start()];
                const point e = sp[edges[seI].end()];
                const point np = help::nearestPointOnTheEdgeExact(s, e, p);

                featureEdgeDistances[i] = labelledScalar(seI, magSqr(np - p));
            }

            // find nearest edges
            sort(featureEdgeDistances);
        }
    }

    // start post-processing gathered data
    const labelList& edgeGroup = partitioner.edgeGroups();

    List<List<labelledScalar>> edgeGroupAndWeights(edges.size());

    # pragma omp parallel for schedule(dynamic, 40) \
    if (edges.size() > 1000)
    forAll(edgeType_, edgeI)
    {
        if (edgeType_[edgeI] & CANDIDATE)
        {
            const edge& e = edges[edgeI];

            const List<labelledScalar>& sc =
                featureEdgesNearPoint[bp[e.start()]];
            const List<labelledScalar>& ec =
                featureEdgesNearPoint[bp[e.end()]];

            // find the feature-edge partition for which the sum of
            // node weights is minimal.
            DynList<labelledScalar> weights;
            forAll(sc, i)
            {
                const label sPart = edgeGroup[sc[i].scalarLabel()];

                forAll(ec, j)
                {
                    const label ePart = edgeGroup[ec[j].scalarLabel()];

                    if ((sPart >= 0) && (sPart == ePart))
                    {
                        weights.append
                        (
                            labelledScalar
                            (
                                sPart,
                                sc[i].value() + ec[j].value()
                            )
                        );
                    }
                }
            }

            // store the data
            edgeGroupAndWeights[edgeI].setSize(weights.size());
            forAll(edgeGroupAndWeights[edgeI], epI)
                edgeGroupAndWeights[edgeI][epI] = weights[epI];

            // sort the data according to the weights
            stableSort(edgeGroupAndWeights[edgeI]);
        }
    }

    Info<< "Edge partitions and weights " << edgeGroupAndWeights << endl;
}


void Foam::Module::edgeExtractor::findNeiPatches
(
    const label bfI,
    const Map<label>& otherProcPatch,
    DynList<label>& neiPatches
) const
{
    const meshSurfaceEngine& mse = surfaceEngine();

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

            neiPatches[feI] = facePatch_[nei];
        }
        else if (edgeFaces.sizeOfRow(beI) == 1)
        {
            neiPatches[feI] = otherProcPatch[beI];
        }
    }
}


Foam::scalar Foam::Module::edgeExtractor::calculateAlignmentForEdge
(
    const label beI,
    const label patch0,
    const label patch1
) const
{
    scalar val(0.0);

    DynList<label> patches(2);
    patches[0] = patch0;
    patches[1] = patch1;

    const pointFieldPMG& points = surfaceEnginePtr_->mesh().points();

    const edge& e = surfaceEnginePtr_->edges()[beI];
    const point& ps = points[e.start()];
    const point& pe = points[e.end()];

    vector ev = e.vec(points);
    const scalar magE = mag(ev) + VSMALL;
    ev /= magE;

    point mps, mpe;
    scalar dSqS, dSqE;

    meshOctree_.findNearestPointToPatches(mps, dSqS, ps, patches);
    meshOctree_.findNearestPointToPatches(mpe, dSqE, pe, patches);

    vector fv = mpe - mps;
    fv /= Foam::max(mag(fv), VSMALL);

    val = 0.5*(1.0 + (ev & fv));

    val = min(val, 1.0);
    val = max(val, 0.0);

    return val;
}


Foam::scalar Foam::Module::edgeExtractor::calculateDeformationMetricForEdge
(
    const label beI,
    const label patch0,
    const label patch1
) const
{
    scalar val(0.0);

    DynList<label> patches(2);
    patches[0] = patch0;
    patches[1] = patch1;

    const pointFieldPMG& points = surfaceEnginePtr_->mesh().points();

    const edge& e = surfaceEnginePtr_->edges()[beI];
    const point& ps = points[e.start()];
    const point& pe = points[e.end()];

    vector ev = e.vec(points);
    const scalar magE = mag(ev) + VSMALL;
    ev /= magE;

    point mps, mpe;
    scalar dSqS, dSqE;

    meshOctree_.findNearestPointToPatches(mps, dSqS, ps, patches);
    meshOctree_.findNearestPointToPatches(mpe, dSqE, pe, patches);

    vector fv = mpe - mps;
    fv /= Foam::max(mag(fv), VSMALL);

    scalar c = min(fv & ev, 1.0);
    c = max(-1.0, c);
    const scalar angle = acos(c);

    val = 0.5*(sqrt(dSqS) + sqrt(dSqE)) + magE*angle;

    return val;
}


Foam::scalar Foam::Module::edgeExtractor::calculateDeformationMetricForFace
(
    const label bfI,
    const DynList<label>& neiPatches,
    const label facePatch
) const
{
    scalar Enorm(0.0);

    const VRWGraph& faceEdges = surfaceEnginePtr_->faceEdges();

    if (neiPatches.size() != faceEdges.sizeOfRow(bfI))
    {
        FatalErrorInFunction
            << "Number of neiPatches and faceEdge does not match for face "
            << bfI << abort(FatalError);
    }

    const label patch0 = facePatch == -1 ? facePatch_[bfI] : facePatch;

    forAllRow(faceEdges, bfI, i)
    {
        const label beI = faceEdges(bfI, i);

        if (neiPatches[i] == patch0)
            continue;

        Enorm += calculateDeformationMetricForEdge(beI, patch0, neiPatches[i]);
    }

    return Enorm;
}


bool Foam::Module::edgeExtractor::checkConcaveEdgeCells()
{
    bool changed(false);

    const triSurf& surf = meshOctree_.surface();
    const VRWGraph& edgeFacets = surf.edgeFacets();

    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    const label bndStartFace = mesh_.boundaries()[0].patchStart();

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const labelList& faceCells = mse.faceOwners();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    // analyse the surface mesh and find out which edges are concave or convex
    const triSurfaceClassifyEdges& edgeClassifier = this->edgeClassifier();
    const List<direction>& edgeType = edgeClassifier.edgeTypes();

    // create a copy of facePatch array for local modifications
    labelList newBoundaryPatches(facePatch_);

    // start checking the surface of the mesh
    label nChanged;

    boolList patchPoint(mse.boundaryPoints().size(), false);

    do
    {
        nChanged = 0;

        // check which surface points are surrounded by boundary faces
        // in the same surface patch
        markPatchPoints(patchPoint);

        // check whether exist edges of a single cell which shall be projected
        // onto a concave edge
        # pragma omp parallel for schedule(dynamic, 40) reduction(+ : nChanged)
        forAll(edgeType_, beI)
        {
            if (!(edgeType_[beI] & SINGLECELLEDGE))
                continue;

            // check if all faces are assigned to the same patch
            const label firstPatch = newBoundaryPatches[edgeFaces(beI, 0)];
            const label secondPatch = newBoundaryPatches[edgeFaces(beI, 1)];

            if (firstPatch == secondPatch)
                continue;

            const cell& c = cells[faceCells[edgeFaces(beI, 0)]];

            // find edges within the bounding box determined by the cell
            point pMin(VGREAT, VGREAT, VGREAT);
            point pMax(-VGREAT, -VGREAT, -VGREAT);
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                {
                    pMin = Foam::min(pMin, points[f[pI]]);
                    pMax = Foam::max(pMax, points[f[pI]]);
                }
            }

            const point cc = 0.5*(pMin + pMax);
            const point diff = pMax - pMin;
            boundBox bb(cc - diff, cc + diff);
            DynList<label> containedEdges;
            meshOctree_.findEdgesInBox(bb, containedEdges);

            // check if there exists concave edges boundaing patches
            // assigned to boundary faces of the current cell
            forAll(containedEdges, ceI)
            {
                const label eI = containedEdges[ceI];

                if (edgeFacets.sizeOfRow(eI) != 2)
                    continue;
                if (!(edgeType[eI] & triSurfaceClassifyEdges::FEATUREEDGE))
                    continue;

                if (edgeType[eI] & triSurfaceClassifyEdges::CONCAVEEDGE)
                {
                    const label patch0 = surf[edgeFacets(eI, 0)].region();
                    const label patch1 = surf[edgeFacets(eI, 1)].region();

                    if
                    (
                        ((firstPatch == patch0) && (secondPatch == patch1)) ||
                        ((firstPatch == patch1) && (secondPatch == patch0))
                    )
                    {
                        // guard the whole read-decide-write sequence below:
                        // two different SINGLECELLEDGE edges belonging to
                        // the same cell can reach this block concurrently
                        // on different threads, and the reassignment both
                        // reads and writes newBoundaryPatches[] for that
                        // cell's own boundary faces. OpenMP critical blocks
                        // must be structured (single entry/exit), so a local
                        // flag stands in for the original "break" and the
                        // actual break happens once we've left the block.
                        bool matched = false;

                        # pragma omp critical (checkConcaveEdgeCellsReassign)
                        {
                        DynList<DynList<label>, 2> facesInPatch;
                        facesInPatch.setSize(2);

                        DynList<label, 2> nFacesInPatch;
                        nFacesInPatch.setSize(2);
                        nFacesInPatch = 0;

                        DynList<bool, 2> hasPatchPoints;
                        hasPatchPoints.setSize(2);
                        hasPatchPoints = false;

                        forAll(c, fI)
                        {
                            if (c[fI] < bndStartFace)
                                continue;

                            const label bfI = c[fI] - bndStartFace;
                            const face& bf = bFaces[bfI];

                            if (newBoundaryPatches[bfI] == patch1)
                            {
                                facesInPatch[1].append(bfI);
                                ++nFacesInPatch[1];

                                forAll(bf, pI)
                                {
                                    if (patchPoint[bp[bf[pI]]])
                                        hasPatchPoints[1] = true;
                                }
                            }
                            else if (newBoundaryPatches[bfI] == patch0)
                            {
                                facesInPatch[0].append(bfI);
                                ++nFacesInPatch[0];

                                forAll(bf, pI)
                                {
                                    if (patchPoint[bp[bf[pI]]])
                                        hasPatchPoints[0] = true;
                                }
                            }
                        }

                        if (nFacesInPatch[1] > nFacesInPatch[0])
                        {
                            // there exist more faces in patch 1
                            // assign all boundary faces to the same patch
                            forAll(facesInPatch[0], i)
                                newBoundaryPatches[facesInPatch[0][i]] = patch1;
                            ++nChanged;
                            matched = true;
                        }
                        else if (nFacesInPatch[0] > nFacesInPatch[1])
                        {
                            // there exist more faces in patch 0
                            // assign all boundary faces to the same patch
                            forAll(facesInPatch[1], i)
                                newBoundaryPatches[facesInPatch[1][i]] = patch0;
                            ++nChanged;
                            matched = true;
                        }
                        else
                        {
                            if (hasPatchPoints[0] && !hasPatchPoints[1])
                            {
                                // transfer all faces to patch 1
                                forAll(facesInPatch[0], i)
                                    newBoundaryPatches[facesInPatch[0][i]] =
                                        patch1;
                                ++nChanged;
                                matched = true;
                            }
                            else if (!hasPatchPoints[0] && hasPatchPoints[1])
                            {
                                // transfer all faces to patch 0
                                forAll(facesInPatch[1], i)
                                    newBoundaryPatches[facesInPatch[1][i]] =
                                        patch0;
                                ++nChanged;
                                matched = true;
                            }
                            else
                            {
                                // just transfer all faces to the same patch
                                forAll(facesInPatch[1], i)
                                    newBoundaryPatches[facesInPatch[1][i]] =
                                        patch0;
                                ++nChanged;
                                matched = true;
                            }
                        }
                        } // end omp critical (checkConcaveEdgeCellsReassign)

                        if (matched)
                            break;
                    }
                }
            }
        }


        if (nChanged)
            changed = true;

    } while (nChanged != 0);

    // transfer the information back to facePatch
    facePatch_.transfer(newBoundaryPatches);

    return changed;
}


bool Foam::Module::edgeExtractor::checkFacePatchesTopology()
{
    bool changed(false);

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    label nCorrected;
    Map<label> otherProcNewPatch;

    label nIter(0);
    do
    {
        nCorrected = 0;

        // allocate a copy of boundary patches
        labelList newBoundaryPatches(facePatch_);

        // find the faces which have neighbouring faces in other patches
        labelLongList candidates;
        findFaceCandidates(candidates, &facePatch_, &otherProcNewPatch);

        // go through the list of faces and check if they shall remain
        // in the current patch
        # pragma omp parallel for schedule(dynamic, 40) \
        reduction(+ : nCorrected)
        forAll(candidates, i)
        {
            const label bfI = candidates[i];
            const face& bf = bFaces[bfI];

            // do not change patches of faces where all points are mapped
            // onto the same patch
            bool allInSamePatch(true);
            forAll(bf, pI)
                if (pointPatch_[bp[bf[pI]]] != facePatch_[bfI])
                {
                    allInSamePatch = false;
                    break;
                }

            if (allInSamePatch)
                continue;

            DynList<label> allNeiPatches;
            DynList<label> neiPatches;
            neiPatches.setSize(faceEdges.sizeOfRow(bfI));

            forAllRow(faceEdges, bfI, eI)
            {
                const label beI = faceEdges(bfI, eI);

                if (edgeFaces.sizeOfRow(beI) == 2)
                {
                    label fNei = edgeFaces(beI, 0);
                    if (fNei == bfI)
                        fNei = edgeFaces(faceEdges(bfI, eI), 1);

                    allNeiPatches.appendUniq(facePatch_[fNei]);
                    neiPatches[eI] = facePatch_[fNei];
                }
                else if (edgeFaces.sizeOfRow(beI) == 1)
                {
                    allNeiPatches.appendUniq(otherProcNewPatch[beI]);
                    neiPatches[eI] = otherProcNewPatch[beI];
                }
            }

            // do not modify faces with all neighbours in the same patch
            if
            (
                (allNeiPatches.size() == 1) &&
                (allNeiPatches[0] == facePatch_[bfI])
            )
                continue;

            // check whether there exist edges which are more suitable for
            // projection onto feature edges than the currently selected ones
            label newPatch(-1);

            // check if some faces have to be distributed to another patch
            // in order to reduce the number of feature edges
            Map<label> nNeiInPatch(allNeiPatches.size());
            forAll(allNeiPatches, i)
                nNeiInPatch.insert(allNeiPatches[i], 0);
            forAll(neiPatches, eI)
                ++nNeiInPatch[neiPatches[eI]];

            newPatch = -1;
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
                    (it() == nNeiEdges) && (it.key() == facePatch_[bfI])
                )
                {
                    newPatch = it.key();
                }
            }

            // do not swap in case the
            if ((newPatch < 0) || (newPatch == facePatch_[bfI]))
                continue;

            // check whether the edges shared ith the neighbour patch form
            // a singly linked chain
            DynList<bool> sharedEdge;
            sharedEdge.setSize(bFaces[bfI].size());
            sharedEdge = false;

            forAll(neiPatches, eI)
                if (neiPatches[eI] == newPatch)
                    sharedEdge[eI] = true;

            if (help::areElementsInChain(sharedEdge))
            {
                // change the patch to the newPatch
                ++nCorrected;
                newBoundaryPatches[bfI] = newPatch;
            }
        }

        // evaluate the new situation and ensure that no oscillation occur
        if (nCorrected)
        {
            faceEvaluator faceEvaluator(*this);

            faceEvaluator.setNewBoundaryPatches(newBoundaryPatches);

            # pragma omp parallel for schedule(dynamic, 50)
            forAll(candidates, i)
            {
                const label bfI = candidates[i];

                const label bestPatch =
                    faceEvaluator.bestPatchAfterModification(bfI);

                newBoundaryPatches[bfI] = bestPatch;
            }
        }

        if (nCorrected)
        {
            changed = true;

            // transfer the new patches back
            facePatch_.transfer(newBoundaryPatches);
        }

    } while (nCorrected != 0 && (nIter++< 3));

    return changed;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{
namespace Module
{
namespace featureEdgeHelpers
{

class featureEdgesNeiOp
{
    // Private data

        // reference to meshSurfaceEngine
        const meshSurfaceEngine& mse_;

        // refence to a list holding information which edges are feature edges
        const boolList& isFeatureEdge_;

        // number of feature edges at a surface point
        labelList nFeatureEdgesAtPoint_;


    // Private member functions

        // calculate the number of feature edges connected to a surface vertex
        void calculateNumberOfEdgesAtPoint()
        {
            const labelList& bp = mse_.bp();
            const edgeList& edges = mse_.edges();

            nFeatureEdgesAtPoint_.setSize(mse_.boundaryPoints().size());
            nFeatureEdgesAtPoint_ = 0;

            forAll(isFeatureEdge_, edgeI)
            {
                if (!isFeatureEdge_[edgeI])
                    continue;

                const edge& e = edges[edgeI];
                ++nFeatureEdgesAtPoint_[bp[e.start()]];
                ++nFeatureEdgesAtPoint_[bp[e.end()]];
            }

        }

public:

    featureEdgesNeiOp
    (
        const meshSurfaceEngine& mse,
        const boolList& isFeatureEdge
    )
    :
        mse_(mse),
        isFeatureEdge_(isFeatureEdge),
        nFeatureEdgesAtPoint_()
    {
        calculateNumberOfEdgesAtPoint();
    }

    label size() const
    {
        return isFeatureEdge_.size();
    }

    void operator()(const label beI, DynList<label>& neighbourEdges) const
    {
        neighbourEdges.clear();

        const VRWGraph& bpEdges = mse_.boundaryPointEdges();
        const labelList& bp = mse_.bp();
        const edgeList& edges = mse_.edges();

        const edge& e = edges[beI];

        const label bps = bp[e.start()];
        const label bpe = bp[e.end()];

        if (nFeatureEdgesAtPoint_[bps] == 2)
        {
            forAllRow(bpEdges, bps, peI)
            {
                const label beJ = bpEdges(bps, peI);

                if ((beJ == beI) || !isFeatureEdge_[beJ])
                    continue;

                neighbourEdges.append(beJ);
            }
        }

        if (nFeatureEdgesAtPoint_[bpe] == 2)
        {
            forAllRow(bpEdges, bpe, peI)
            {
                const label beJ = bpEdges(bpe, peI);

                if ((beJ == beI) || !isFeatureEdge_[beJ])
                    continue;

                neighbourEdges.append(beJ);
            }
        }
    }

};


class featureEdgesSelOp
{
    // Private data

        // reference to a list holding information which edge is afeture edge
        const boolList& isFeatureEdge_;

public:

    featureEdgesSelOp(const boolList& isFeatureEdge)
    :
        isFeatureEdge_(isFeatureEdge)
    {}


    bool operator()(const label beI) const
    {
        return isFeatureEdge_[beI];
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace featureEdgeHelpers
} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Module::edgeExtractor::checkFacePatchesGeometry()
{
    bool changed(false);

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();

    // allocate a copy of boundary patches
    labelList newBoundaryPatches(facePatch_.size());

    label nCorrected;
    Map<label> otherProcNewPatch;

    boolList activePoints(bPoints.size(), true);
    labelLongList activePointLabel(bPoints.size());
    forAll(bPoints, bpI)
        activePointLabel[bpI] = bpI;

    label iter(0);

    do
    {
        // create feature edges and corners information
        meshSurfacePartitioner mPart(mse, facePatch_);

        // project vertices onto the surface mesh
        meshSurfaceMapper(mPart, meshOctree_).mapVerticesOntoSurfacePatches
        (
            activePointLabel
        );

        // stop after a certain number of iterations
        if (++iter > 3)
            break;

        // update surface geometry data
        meshSurfaceEngineModifier(mse).updateGeometry();

        // check if there exist any inverted faces
        meshSurfaceCheckInvertedVertices surfCheck(mPart, activePoints);
        const labelHashSet& invertedPoints = surfCheck.invertedVertices();

        if (invertedPoints.size() == 0 )
            return false;

        WarningInFunction
            << "Found " << invertedPoints.size()
            << " points with inverted surface normals. Getting rid of them..."
            << endl;

        // untangle the surface
        activePointLabel.clear();
        activePoints = false;
        for (const label pointi : invertedPoints.toc())
        {
            activePointLabel.append(bp[pointi]);
            activePoints[bp[pointi]] = true;
        }

        // untangle the surface
        meshSurfaceOptimizer mso(mPart, meshOctree_);
        mso.untangleSurface(activePointLabel, 1);

        nCorrected = 0;
        newBoundaryPatches = facePatch_;

        // find the faces which have neighbouring faces in other patches
        labelLongList candidates;
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
                if (invertedPoints.found(bf[pI]))
                {
                    candidates.append(bfI);
                    break;
                }
        }

        # pragma omp parallel for schedule(dynamic, 5) reduction(+ : nCorrected)
        forAll(candidates, i)
        {
            const label bfI = candidates[i];

            DynList<label> neiPatches;
            findNeiPatches(bfI, otherProcNewPatch, neiPatches);

            DynList<label> allNeiPatches;
            forAll(neiPatches, i)
                allNeiPatches.appendUniq(neiPatches[i]);

            // check the deformation energy and find the minimum energy which
            // can be achieved by switching face patch
            scalar minE(VGREAT);
            label minEPatch(-1);
            DynList<scalar> Enorm(allNeiPatches.size());
            forAll(allNeiPatches, i)
            {
                Enorm[i] =
                    calculateDeformationMetricForFace
                    (
                        bfI,
                        neiPatches,
                        allNeiPatches[i]
                    );

                if (Enorm[i] < minE)
                {
                    minE = Enorm[i];
                    minEPatch = allNeiPatches[i];
                }
            }

            if (minEPatch != facePatch_[bfI])
            {
                newBoundaryPatches[bfI] = minEPatch;
                ++nCorrected;
            }
        }

        // check if any faces are re-assigned to some other patch
        if (nCorrected == 0)
            break;

        // update faceEvaluator with information after patches have been
        // altered. It blocks chaning of patches if it causes oscillations
        faceEvaluator facePatchEvaluator(*this);
        facePatchEvaluator.setNewBoundaryPatches(newBoundaryPatches);

        // compare face patches before and after
        // disallow modification which may trigger oscillating behaviour
        labelLongList changedFaces;
        forAll(newBoundaryPatches, bfI)
        {
            if (newBoundaryPatches[bfI] != facePatch_[bfI])
            {
                const label patchI =
                    facePatchEvaluator.bestPatchAfterModification(bfI);
                newBoundaryPatches[bfI] = patchI;

                if (patchI != facePatch_[bfI])
                    changedFaces.append(bfI);
            }
        }

        nCorrected = changedFaces.size();

        if (nCorrected)
        {
            changed = true;
            facePatch_ = newBoundaryPatches;
        }

    } while (nCorrected != 0);

    return changed;
}


void Foam::Module::edgeExtractor::projectDeterminedFeatureVertices()
{
    List<DynList<label, 5>> pointPatches;
    pointPatches.setSize(pointValence_.size());

    const meshSurfaceEngine& mse = surfaceEngine();
    const pointFieldPMG& points = mse.mesh().points();
    const labelList& bPoints = mse.boundaryPoints();
    const labelList& bp = mse.bp();
    const faceList::subList& bFaces = mse.boundaryFaces();
    meshOctree_.surface().pointEdges();

    // calculate patches for each point
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
            pointPatches[bp[bf[pI]]].appendUniq(facePatch_[bfI]);
    }

    meshSurfaceEngineModifier surfMod(mse);

    # pragma omp parallel for schedule(dynamic, 10)
    forAll(pointPatches, bpI)
    {
        if (pointPatches[bpI].size() < 2)
            continue;

        const DynList<label> pPatches = pointPatches[bpI];

        const point& p = points[bPoints[bpI]];

        // find the nearest object on the surface mesh
        point newP;
        scalar dSqExact;
        if (pPatches.size() == 2)
        {
            label nse;
            meshOctree_.findNearestEdgePoint(newP, dSqExact, nse, p, pPatches);
        }
        else
        {
            label nsp;
            meshOctree_.findNearestCorner(newP, dSqExact, nsp, p, pPatches);
        }

        // find the nearest object in an iterative procedure
        point pp(p);
        for (label iterI = 0; iterI < 20; ++iterI)
        {
            point inp(vector::zero);

            forAll(pPatches, i)
            {
                point np;
                scalar dSq;
                label nt;

                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    dSq,
                    nt,
                    pPatches[i],
                    pp
                );

                inp += np;
            }

            inp /= pPatches.size();
            const scalar currDSq = magSqr(inp - pp);
            pp = inp;

            if (currDSq < 1e-2*dSqExact)
                break;
        }

        // check if the exact position of the corner is further away
        // than the iteratively found object
        if (dSqExact > 1.1*magSqr(pp - p))
            newP = pp;

        surfMod.moveBoundaryVertexNoUpdate(bpI, newP);
    }

    surfMod.updateGeometry();
}


bool Foam::Module::edgeExtractor::untangleSurface()
{
    bool changed(false);

    meshSurfaceEngine& mse =
        const_cast<meshSurfaceEngine&>(this->surfaceEngine());
    meshSurfaceOptimizer optimizer(mse, meshOctree_);
    changed = optimizer.untangleSurface();

    return changed;
}


void Foam::Module::edgeExtractor::extractEdges()
{
    distributeBoundaryFaces();

    distributeBoundaryFacesNormalAlignment();

    Info<< "Starting topological adjustment of patches" << endl;
    if (checkFacePatchesTopology())
    {
        Info<< "Finished topological adjustment of patches" << endl;
    }
    else
    {
        Info<< "No topological adjustment was needed" << endl;
    }

    Info<< "Starting geometrical adjustment of patches" << endl;
    if (checkFacePatchesGeometry())
    {
        Info<< "Finished geometrical adjustment of patches" << endl;
    }
    else
    {
        Info<< "No geometrical adjustment was needed" << endl;
    }
}


const Foam::Module::triSurf*
Foam::Module::edgeExtractor::surfaceWithPatches() const
{
    // allocate the memory for the surface mesh
    triSurf* surfPtr = new triSurf();

    // surface of the volume mesh
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const pointFieldPMG& points = mesh_.points();

    // modifier of the new surface mesh
    triSurfModifier surfModifier(*surfPtr);
    surfModifier.patchesAccess() = meshOctree_.surface().patches();
    pointField& sPts = surfModifier.pointsAccess();
    sPts.setSize(mse.boundaryPoints().size());

    // copy points
    forAll(bp, pointI)
    {
        if (bp[pointI] < 0)
            continue;

        sPts[bp[pointI]] = points[pointI];
    }

    // create the triangulation of the volume mesh surface
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        labelledTri tri;
        tri.region() = facePatch_[bfI];
        tri[0] = bp[bf[0]];

        for (label i = bf.size()-2; i > 0; --i)
        {
            tri[1] = bp[bf[i]];
            tri[2] = bp[bf[i + 1]];

            surfPtr->appendTriangle(tri);
        }
    }

    return surfPtr;
}


const Foam::Module::triSurf*
Foam::Module::edgeExtractor::surfaceWithPatches(const label bpI) const
{
    // allocate the memory for the surface mesh
    triSurf* surfPtr = new triSurf();

    // surface of the volume mesh
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pFaces = mse.pointFaces();
    const pointFieldPMG& points = mesh_.points();

    // modifier of the new surface mesh
    triSurfModifier surfModifier(*surfPtr);
    surfModifier.patchesAccess() = meshOctree_.surface().patches();
    pointField& sPts = surfModifier.pointsAccess();

    // create the triangulation of the volume mesh surface
    labelLongList newPointLabel(points.size(), -1);
    label nPoints(0);
    forAllRow(pFaces, bpI, pfI)
    {
        const label bfI = pFaces(bpI, pfI);
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
            if (newPointLabel[bf[pI]] == -1)
                newPointLabel[bf[pI]] = nPoints++;

        labelledTri tri;
        tri.region() = facePatch_[bfI];
        tri[0] = newPointLabel[bf[0]];

        for (label i = bf.size()-2; i > 0; --i)
        {
            tri[1] = newPointLabel[bf[i]];
            tri[2] = newPointLabel[bf[i + 1]];

            surfPtr->appendTriangle(tri);
        }
    }

    // copy points
    sPts.setSize(nPoints);
    forAll(newPointLabel, pointI)
        if (newPointLabel[pointI] != -1)
        {
            sPts[newPointLabel[pointI]] = points[pointI];
        }

    return surfPtr;
}


void Foam::Module::edgeExtractor::updateMeshPatches()
{
    const triSurf& surface = meshOctree_.surface();
    const geometricSurfacePatchList& surfPatches = surface.patches();
    const label nPatches = surfPatches.size();

    const meshSurfaceEngine& mse = this->surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& faceOwner = mse.faceOwners();

    wordList patchNames(nPatches);
    VRWGraph newBoundaryFaces;
    labelLongList newBoundaryOwners(bFaces.size());
    labelLongList newBoundaryPatches(bFaces.size());

    // set patchNames
    forAll(surface.patches(), patchI)
        patchNames[patchI] = surface.patches()[patchI].name();

    // append boundary faces
    forAll(bFaces, bfI)
    {
        newBoundaryFaces.appendList(bFaces[bfI]);
        newBoundaryOwners[bfI] = faceOwner[bfI];
        newBoundaryPatches[bfI] = facePatch_[bfI];
    }

    // replace the boundary with the new patches
    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.replaceBoundary
    (
        patchNames,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );

    // set the new patch types
    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();

    forAll(surfPatches, patchI)
        boundaries[patchI].patchType() = surfPatches[patchI].geometricType();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
