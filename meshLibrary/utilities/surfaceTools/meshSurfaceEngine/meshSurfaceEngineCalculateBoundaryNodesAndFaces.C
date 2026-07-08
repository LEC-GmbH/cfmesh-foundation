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

#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"
#include "boolList.H"
#include "helperFunctions.H"
#include "VRWGraphSMPModifier.H"
#include "labelledPoint.H"
#include "HashSet.H"

#include <map>
#include <set>

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceEngine::calculateBoundaryFaces() const
{
    if (mesh_.boundaries().size() != 0)
    {
        const faceListPMG& faces = mesh_.faces();
        const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

        label nBoundaryFaces(0);
        if (activePatch_ < 0)
        {
            // take all patches
            forAll(boundaries, patchI)
            {
                nBoundaryFaces += boundaries[patchI].patchSize();
            }

            boundaryFacesPtr_ =
                new faceList::subList
                (
                    faces,
                    nBoundaryFaces,
                    boundaries[0].patchStart()
                );
        }
        else if (activePatch_ < boundaries.size())
        {
            nBoundaryFaces = boundaries[activePatch_].patchSize();

            boundaryFacesPtr_ =
                new faceList::subList
                (
                    faces,
                    nBoundaryFaces,
                    boundaries[activePatch_].patchStart()
                );
        }
        else
        {
            FatalErrorInFunction
                << "Cannot select boundary faces. Invalid patch index "
                << activePatch_ << exit(FatalError);
        }

        Info<< "Found " << nBoundaryFaces << " boundary faces " << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Boundary faces are not at the end of the face list!"
            << exit(FatalError);
    }
}


void Foam::Module::meshSurfaceEngine::calculateBoundaryOwners() const
{
    const labelList& owner = mesh_.owner();

    const faceList::subList& boundaryFaces = this->boundaryFaces();

    if (!boundaryFaceOwnersPtr_)
    {
        boundaryFaceOwnersPtr_ = new labelList(boundaryFaces.size());
    }

    labelList& owners = *boundaryFaceOwnersPtr_;

    const label start = mesh_.boundaries()[0].patchStart();

    # pragma omp parallel for schedule(static, 1)
    forAll(boundaryFaces, fI)
    owners[fI] = owner[start + fI];
}


void Foam::Module::meshSurfaceEngine::calculateBoundaryNodes() const
{
    // mark boundary points
    label pointI(0);
    if (!bppPtr_)
    {
        bppPtr_ = new labelList(mesh_.points().size(), -1);
    }
    labelList& bp = *bppPtr_;

    const faceList::subList& boundaryFaces = this->boundaryFaces();

    boolList isBndPoint(bp.size(), false);

    const label nThreads = 3*omp_get_num_procs();
    # pragma omp parallel for num_threads(nThreads) schedule(static, 1)
    forAll(boundaryFaces, bfI)
    {
        const face& bf = boundaryFaces[bfI];

        forAll(bf, pI)
        {
            isBndPoint[bf[pI]] = true;
        }
    }

    forAll(isBndPoint, pI)
    {
        if (isBndPoint[pI])
        {
            bp[pI] = pointI++;
        }
    }


    if (!boundaryPointsPtr_)
    {
        boundaryPointsPtr_ = new labelList();
    }
    labelList& boundaryPoints = *boundaryPointsPtr_;
    boundaryPoints.setSize(pointI);

    // fill the boundaryPoints list
    # pragma omp parallel for num_threads(nThreads) schedule(static, 1)
    forAll(bp, bpI)
    {
        if (bp[bpI] != -1)
        {
            boundaryPoints[bp[bpI]] = bpI;
        }
    }
}


void Foam::Module::meshSurfaceEngine::calculateBoundaryFacePatches() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    boundaryFacePatchPtr_ = new labelList(bFaces.size());
    labelList& facePatch = *boundaryFacePatchPtr_;

    label faceI(0);
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        const label nFaces = boundaries[patchI].patchSize();
        for (label patchFaceI = 0; patchFaceI < nFaces; ++patchFaceI)
        {
            facePatch[faceI] = patchI;
            ++faceI;
        }
    }
}


void Foam::Module::meshSurfaceEngine::calculatePointFaces() const
{
    // fill pointFacesAddr
    if (!pointFacesPtr_)
    {
        pointFacesPtr_ = new VRWGraph();
    }
    VRWGraph& pointFacesAddr = *pointFacesPtr_;

    if (!pointInFacePtr_)
    {
        pointInFacePtr_ = new VRWGraph();
    }
    VRWGraph& pointInFaceAddr = *pointInFacePtr_;

    const labelList& bPoints = this->boundaryPoints();
    const faceList::subList& bFaces = this->boundaryFaces();

    // create boundary points
    const labelList& bp = this->bp();

    labelLongList npf;

    label nThreads = 3*omp_get_num_procs();
    if (bPoints.size() < 1000)
    {
        nThreads = 1;
    }

    label minRow(INT_MAX), maxRow(0);
    List<List<LongList<labelPair>>> dataForOtherThreads(nThreads);

    # pragma omp parallel num_threads(nThreads)
    {
        const label threadI = omp_get_thread_num();

        List<LongList<labelPair>>& dot = dataForOtherThreads[threadI];
        dot.setSize(nThreads);

        // find min and max entry in the graph
        // they are used for assigning ranges of values local for each process
        label localMinRow(minRow), localMaxRow(0);
        # pragma omp for schedule(static)
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];
            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];
                localMaxRow = Foam::max(localMaxRow, bpI);
                localMinRow = Foam::min(localMinRow, bpI);
            }
        }

        ++localMaxRow;

        # pragma omp critical
        {
            minRow = Foam::min(minRow, localMinRow);
            minRow = Foam::max(minRow, 0);
            maxRow = Foam::max(maxRow, localMaxRow);

            npf.setSize(maxRow);
        }

        # pragma omp barrier

        // initialise appearances
        # pragma omp for schedule(static)
        for (label i = 0; i < maxRow; ++i)
        {
            npf[i] = 0;
        }
        # pragma omp barrier

        const label range = (maxRow - minRow) / nThreads + 1;
        const label localMin = minRow + threadI*range;
        const label localMax = Foam::min(localMin + range, maxRow);

        // find the number of appearances of each element in the original graph
        # pragma omp for schedule(static)
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];

                const label threadNo = (bpI - minRow) / range;

                if (threadNo == threadI)
                {
                    ++npf[bpI];
                }
                else
                {
                    dot[threadNo].append(labelPair(bpI, bfI));
                }
            }
        }

        # pragma omp barrier

        // count the appearances which are not local to the processor
        for (label i = 0; i < nThreads; ++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];

            forAll(data, j)
            {
                ++npf[data[j].first()];
            }
        }

        # pragma omp barrier

        // allocate graph
        # pragma omp master
        {
            VRWGraphSMPModifier(pointFacesAddr).setSizeAndRowSize(npf);
            VRWGraphSMPModifier(pointInFaceAddr).setSizeAndRowSize(npf);
        }

        # pragma omp barrier

        for (label i = localMin; i < localMax; ++i)
        {
            npf[i] = 0;
        }

        // start filling reverse addressing graph
        // update data from processors with smaller labels
        for (label i = 0; i < threadI; ++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];

            forAll(data, j)
            {
                const label bpI = data[j].first();
                const label bfI = data[j].second();

                pointFacesAddr(bpI, npf[bpI]) = bfI;
                pointInFaceAddr(bpI, npf[bpI]) =
                    bFaces[bfI].which(bPoints[bpI]);

                ++npf[bpI];
            }
        }

        // update data local to the processor
        # pragma omp for schedule(static)
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];

                if ((bpI >= localMin) && (bpI < localMax))
                {
                    pointInFaceAddr(bpI, npf[bpI]) = pI;
                    pointFacesAddr(bpI, npf[bpI]++) = bfI;
                }
            }
        }

        // update data from the processors with higher labels
        for (label i = threadI + 1; i < nThreads; ++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];

            forAll(data, j)
            {
                const label bpI = data[j].first();
                const label bfI = data[j].second();

                pointFacesAddr(bpI, npf[bpI]) = bfI;
                pointInFaceAddr(bpI, npf[bpI]) =
                    bFaces[bfI].which(bPoints[bpI]);

                ++npf[bpI];
            }
        }
    }

    pointFacesAddr.setSize(bPoints.size());
    pointInFaceAddr.setSize(bPoints.size());
}


void Foam::Module::meshSurfaceEngine::calculatePointPatches() const
{
    if (!pointPatchesPtr_)
    {
        pointPatchesPtr_ = new VRWGraph();
    }
    VRWGraph& pPatches = *pointPatchesPtr_;

    const labelList& facePatch = boundaryFacePatches();
    const VRWGraph& pFaces = pointFaces();

    const label nThreads = 3*omp_get_num_procs();

    labelList npPatches(pFaces.size());

    # pragma omp parallel num_threads(nThreads)
    {
        # pragma omp for schedule(static)
        forAll(npPatches, i)
        {
            npPatches[i] = 0;
        }

        # pragma omp for schedule(static)
        forAll(pFaces, bpI)
        {
            DynList<label> pf;
            forAllRow(pFaces, bpI, pfI)
            {
                pf.appendUniq(facePatch[pFaces(bpI, pfI)]);
            }

            npPatches[bpI] = pf.size();
        }

        # pragma omp barrier

        # pragma omp master
        VRWGraphSMPModifier(pPatches).setSizeAndRowSize(npPatches);

        # pragma omp barrier

        # pragma omp for schedule(static)
        forAll(pFaces, bpI)
        {
            DynList<label> pf;
            forAllRow(pFaces, bpI, pfI)
            {
                pf.appendUniq(facePatch[pFaces(bpI, pfI)]);
            }

            pPatches.setRow(bpI, pf);
        }
    }
}


void Foam::Module::meshSurfaceEngine::calculatePointPoints() const
{
    if (!pointPointsPtr_)
    {
        pointPointsPtr_ = new VRWGraph();
    }

    VRWGraph& pointPoints = *pointPointsPtr_;

    const labelList& boundaryPoints = this->boundaryPoints();
    const faceList::subList& bFaces = this->boundaryFaces();
    const VRWGraph& pFaces = this->pointFaces();
    const labelList& bp = this->bp();

    const label nThreads = 3*omp_get_num_procs();

    labelList npp(boundaryPoints.size());

    # pragma omp parallel num_threads(nThreads)
    {
        # pragma omp for schedule(static)
        forAll(npp, i)
        {
            npp[i] = 0;
        }

        # pragma omp for schedule(static)
        forAll(pFaces, bpI)
        {
            DynList<label> pPoints;

            forAllRow(pFaces, bpI, pfI)
            {
                const face& bf = bFaces[pFaces(bpI, pfI)];

                const label pos = bf.which(boundaryPoints[bpI]);

                pPoints.appendUniq(bp[bf.nextLabel(pos)]);
                pPoints.appendUniq(bp[bf.prevLabel(pos)]);
            }

            npp[bpI] = pPoints.size();
        }

        # pragma omp barrier

        # pragma omp master
        VRWGraphSMPModifier(pointPoints).setSizeAndRowSize(npp);

        # pragma omp barrier

        # pragma omp for schedule(static)
        forAll(pFaces, bpI)
        {
            DynList<label> pPoints;

            forAllRow(pFaces, bpI, pfI)
            {
                const face& bf = bFaces[pFaces(bpI, pfI)];

                const label pos = bf.which(boundaryPoints[bpI]);

                pPoints.appendUniq(bp[bf.nextLabel(pos)]);
                pPoints.appendUniq(bp[bf.prevLabel(pos)]);
            }

            pointPoints.setRow(bpI, pPoints);
        }
    }
}


void Foam::Module::meshSurfaceEngine::calculatePointNormals() const
{
    const VRWGraph& pFaces = pointFaces();
    const vectorField& fNormals = faceNormals();

    pointNormalsPtr_ = new vectorField(pFaces.size());

    # pragma omp parallel for if (pFaces.size() > 1000) schedule(dynamic, 50)
    forAll(pFaces, pI)
    {
        vector normal(vector::zero);

        forAllRow(pFaces, pI, pfI)
        {
            normal += fNormals[pFaces(pI, pfI)];
        }

        const scalar d = mag(normal);
        if (d > VSMALL)
        {
            normal /= d;
        }
        else
        {
            normal = vector::zero;
        }

        (*pointNormalsPtr_)[pI] = normal;
    }
}


void Foam::Module::meshSurfaceEngine::calculateFaceNormals() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    const pointFieldPMG& points = mesh_.points();

    faceNormalsPtr_ = new vectorField(bFaces.size());

    # pragma omp parallel for if (bFaces.size() > 1000)
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        faceNormalsPtr_->operator[](bfI) = bf.areaNormal(points);
    }
}


void Foam::Module::meshSurfaceEngine::calculateFaceCentres() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    const pointFieldPMG& points = mesh_.points();

    faceCentresPtr_ = new vectorField(bFaces.size());

    # pragma omp parallel for if (bFaces.size() > 1000)
    forAll(bFaces, bfI)
    {
        faceCentresPtr_->operator[](bfI) = bFaces[bfI].centre(points);
    }
}



void Foam::Module::meshSurfaceEngine::calculateEdgesAndAddressing() const
{
    const VRWGraph& pFaces = pointFaces();
    const faceList::subList& bFaces = boundaryFaces();
    const labelList& bPoints = boundaryPoints();
    const labelList& bp = this->bp();

    edgesPtr_ = new edgeList();
    edgeList& edges = *edgesPtr_;

    bpEdgesPtr_ = new VRWGraph();
    VRWGraph& bpEdges = *bpEdgesPtr_;

    label nThreads = 3*omp_get_num_procs();
    if (pFaces.size() < 1000)
    {
        nThreads = 1;
    }

    labelList nEdgesForThread(nThreads);

    label edgeI(0);

    # pragma omp parallel num_threads(nThreads)
    {
        edgeLongList edgesHelper;

        # pragma omp for schedule(static)
        forAll(pFaces, bpI)
        {
            std::set<std::pair<label, label>> edgesAtPoint;

            forAllRow(pFaces, bpI, pfI)
            {
                const label bfI = pFaces(bpI, pfI);
                const face& bf = bFaces[bfI];

                const label pos = bf.which(bPoints[bpI]);

                if (bp[bf.nextLabel(pos)] >= bpI)
                {
                    edgesAtPoint.insert
                    (
                        std::make_pair(bf[pos], bf.nextLabel(pos))
                    );
                }
                if (bp[bf.prevLabel(pos)] >= bpI)
                {
                    edgesAtPoint.insert
                    (
                        std::make_pair(bf[pos], bf.prevLabel(pos))
                    );
                }
            }

            forAllConstIters(edgesAtPoint, it)
            {
                edgesHelper.append(edge(it->first, it->second));
            }
        }

        // this enables other threads to see the number of edges
        // generated by each thread
        const label threadI = omp_get_thread_num();
        nEdgesForThread[threadI] = edgesHelper.size();

        # pragma omp critical
        edgeI += edgesHelper.size();

        # pragma omp barrier

        # pragma omp master
        edgesPtr_->setSize(edgeI);

        # pragma omp barrier

        // find the starting position of the edges generated by this thread
        // in the global list of edges
        label localStart(0);
        for (label i = 0; i < threadI; ++i)
        {
            localStart += nEdgesForThread[i];
        }

        // store edges into the global list
        forAll(edgesHelper, i)
        {
            edgesPtr_->operator[](localStart + i) = edgesHelper[i];
        }
    }

    // set the bpEdges
    VRWGraphSMPModifier(bpEdges).reverseAddressing(bp, edges);
    bpEdges.setSize(pFaces.size());
}


void Foam::Module::meshSurfaceEngine::calculateFaceEdgesAddressing() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    const labelList& bp = this->bp();
    const edgeList& edges = this->edges();
    const VRWGraph& bpEdges = this->boundaryPointEdges();

    faceEdgesPtr_ = new VRWGraph(bFaces.size());
    VRWGraph& faceEdges = *faceEdgesPtr_;

    labelList nfe(bFaces.size());

    const label nThreads = 3*omp_get_num_procs();

    # pragma omp parallel num_threads(nThreads)
    {
        # pragma omp for schedule(static)
        forAll(bFaces, bfI)
        {
            nfe[bfI] = bFaces[bfI].size();
        }

        # pragma omp barrier

        # pragma omp master
        {

        VRWGraphSMPModifier(faceEdges).setSizeAndRowSize(nfe);

        }

        # pragma omp barrier

        # pragma omp for schedule(dynamic, 100)
        forAll(faceEdges, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, eI)
            {
                const edge e = bf.faceEdge(eI);

                const label bps = bp[e.start()];

                forAllRow(bpEdges, bps, peI)
                {
                    const label beI = bpEdges(bps, peI);
                    const edge& ee = edges[beI];

                    if (e == ee)
                    {
                        faceEdges(bfI, eI) = beI;
                        break;
                    }
                }
            }
        }
    }
}


void Foam::Module::meshSurfaceEngine::calculateEdgeFacesAddressing() const
{
    const faceList::subList& bFaces = this->boundaryFaces();
    const VRWGraph& pointFaces = this->pointFaces();
    const edgeList& edges = this->edges();
    const labelList& bp = this->bp();

    edgeFacesPtr_ = new VRWGraph();
    VRWGraph& edgeFaces = *edgeFacesPtr_;

    labelList nef(edges.size());

    const label nThreads = 3*omp_get_num_procs();

    # pragma omp parallel num_threads(nThreads)
    {
        # pragma omp for schedule(static)
        forAll(nef, edgeI)
        {
            nef[edgeI] = 0;
        }

        # pragma omp for schedule(static)
        forAll(edges, edgeI)
        {
            const edge& ee = edges[edgeI];
            const label bpI = bp[ee.start()];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);

                const face& bf = bFaces[bfI];

                forAll(bf, eI)
                {
                    if (bf.faceEdge(eI) == ee)
                    {
                        ++nef[edgeI];
                        break;
                    }
                }
            }
        }

        # pragma omp barrier

        # pragma omp master
        VRWGraphSMPModifier(edgeFaces).setSizeAndRowSize(nef);

        # pragma omp barrier

        # pragma omp for schedule(static)
        forAll(edges, edgeI)
        {
            const edge& ee = edges[edgeI];
            const label bpI = bp[ee.start()];

            // find boundary faces attached to this edge
            DynList<label> eFaces;
            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);

                const face& bf = bFaces[bfI];

                forAll(bf, eI)
                {
                    if (bf.faceEdge(eI) == ee)
                    {
                        eFaces.append(bfI);
                        break;
                    }
                }
            }

            // the face that owns the edge shall be the first one in the list
            // TODO: find out whether this will be necessary
            if (eFaces.size() == 2)
            {
                const face& bf = bFaces[eFaces[1]];

                const label pos = bf.which(ee.start());

                if (bf.nextLabel(pos) == ee.end())
                {
                    // this face shall be the first one in the list
                    const label helper = eFaces[0];
                    eFaces[0] = eFaces[1];
                    eFaces[1] = helper;
                }
            }

            edgeFaces.setRow(edgeI, eFaces);
        }
    }
}


void Foam::Module::meshSurfaceEngine::calculateEdgePatchesAddressing() const
{
    edgePatchesPtr_ = new VRWGraph();
    VRWGraph& edgePatches = *edgePatchesPtr_;

    const VRWGraph& edgeFaces = this->edgeFaces();
    const labelList& facePatch = this->boundaryFacePatches();

    edgePatches.setSize(edgeFaces.size());

    forAll(edgeFaces, eI)
    {
        DynList<label> ePatches;

        forAllRow(edgeFaces, eI, i)
        {
            const label patchI = facePatch[edgeFaces(eI, i)];

            ePatches.appendUniq(patchI);
        }

        edgePatches.setRow(eI, ePatches);
    }
}


void Foam::Module::meshSurfaceEngine::calculateFaceFacesAddressing() const
{
    const VRWGraph& edgeFaces = this->edgeFaces();

    const faceList::subList& bFaces = boundaryFaces();
    faceFacesPtr_ = new VRWGraph(bFaces.size());
    VRWGraph& faceFaces = *faceFacesPtr_;

    forAll(bFaces, bfI)
    {
        faceFaces.setRowSize(bfI, bFaces[bfI].size());
    }

    labelList nAppearances(bFaces.size(), 0);

    forAll(edgeFaces, efI)
    {
        if (edgeFaces.sizeOfRow(efI) == 2)
        {
            const label f0 = edgeFaces(efI, 0);
            const label f1 = edgeFaces(efI, 1);

            faceFaces(f0, nAppearances[f0]++) = f1;
            faceFaces(f1, nAppearances[f1]++) = f0;
        }
    }
}


// ************************************************************************* //
