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
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshOctree.H"
#include "triangle.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceMapper2D.H"
#include "polyMeshGen2DEngine.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGenChecks.H"
#include "labelledPoint.H"
#include "FIFOStack.H"

#include <map>
#include <stdexcept>

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::Module::meshSurfaceOptimizer::findInvertedVertices
(
    boolList& smoothVertex,
    const label nAdditionalLayers
) const
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const labelListList2D& pPoints = surfaceEngine_.pointPoints();

    if (smoothVertex.size() != bPoints.size())
    {
        smoothVertex.setSize(bPoints.size());
        smoothVertex = true;
    }

    label nInvertedTria(0);

    // check the vertices at the surface
    // mark the ones where the mesh is tangled
    meshSurfaceCheckInvertedVertices vrtCheck(*partitionerPtr_, smoothVertex);
    const labelHashSet& inverted = vrtCheck.invertedVertices();

    smoothVertex = false;
    forAll(bPoints, bpI)
    {
        if (inverted.found(bPoints[bpI]))
        {
            ++nInvertedTria;
            smoothVertex[bpI] = true;
        }
    }

    Info<< "Number of inverted boundary faces is " << nInvertedTria << endl;

    if (nInvertedTria == 0)
        return 0;

    // add additional layers around inverted points
    for (label i = 0; i < nAdditionalLayers; ++i)
    {
        boolList originallySelected = smoothVertex;
        forAll(smoothVertex, bpI)
            if (originallySelected[bpI])
                forAllRow(pPoints, bpI, ppI)
                    smoothVertex[pPoints[bpI][ppI]] = true;

    }

    return nInvertedTria;
}


void Foam::Module::meshSurfaceOptimizer::smoothEdgePoints
(
    const DynamicList<label>& edgePoints
)
{
    List<DynamicList<labelledPoint>> newPositions(1);
    newPositions.setSize(omp_get_num_procs());

    // smooth edge vertices
    # pragma omp parallel num_threads(newPositions.size())
    {
        DynamicList<labelledPoint>& newPos =
            newPositions[omp_get_thread_num()];

        # pragma omp for schedule(dynamic, 40)
        forAll(edgePoints, i)
        {
            const label bpI = edgePoints[i];

            if (vertexType_[bpI] & LOCKED)
                continue;

            newPos.append(labelledPoint(bpI, newEdgePositionLaplacian(bpI)));
        }
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    forAll(newPositions, threadI)
    {
        const DynamicList<labelledPoint>& newPos = newPositions[threadI];

        forAll(newPos, i)
            surfaceModifier.moveBoundaryVertexNoUpdate
            (
                newPos[i].pointLabel(),
                newPos[i].coordinates()
            );
    }

    surfaceModifier.updateGeometry(edgePoints);
}


void Foam::Module::meshSurfaceOptimizer::smoothLaplacianFC
(
    const DynamicList<label>& selectedPoints,
    const bool transform
)
{
    List<DynamicList<labelledPoint>> newPositions(1);
    newPositions.setSize(omp_get_num_procs());

    # pragma omp parallel num_threads(newPositions.size())
    {
        DynamicList<labelledPoint>& newPos =
            newPositions[omp_get_thread_num()];

        # pragma omp for schedule(dynamic, 40)
        forAll(selectedPoints, i)
        {
            const label bpI = selectedPoints[i];

            if (vertexType_[bpI] & LOCKED)
                continue;

            newPos.append
            (
                labelledPoint(bpI, newPositionLaplacianFC(bpI, transform))
            );
        }
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    # pragma omp parallel num_threads(newPositions.size())
    {
        const label threadI = omp_get_thread_num();

        const DynamicList<labelledPoint>& newPos = newPositions[threadI];

        forAll(newPos, i)
            surfaceModifier.moveBoundaryVertexNoUpdate
            (
                newPos[i].pointLabel(),
                newPos[i].coordinates()
            );
    }

    surfaceModifier.updateGeometry(selectedPoints);
}


void Foam::Module::meshSurfaceOptimizer::smoothSurfaceOptimizer
(
    const DynamicList<label>& selectedPoints
)
{
    // create partTriMesh is it is not yet present
    this->triMesh();

    // update coordinates of the triangulation
    updateTriMesh(selectedPoints);

    pointField newPositions(selectedPoints.size());

    # pragma omp parallel for schedule(dynamic, 20)
    forAll(selectedPoints, i)
    {
        const label bpI = selectedPoints[i];

        newPositions[i] = newPositionSurfaceOptimizer(bpI);
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    # pragma omp parallel for schedule(dynamic, 100)
    forAll(newPositions, i)
    {
        const label bpI = selectedPoints[i];

        surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newPositions[i]);
    }

    // update geometry addressing for moved points
    surfaceModifier.updateGeometry(selectedPoints);
}


bool Foam::Module::meshSurfaceOptimizer::untangleSurface
(
    const DynamicList<label>& selectedBoundaryPoints,
    const label nAdditionalLayers
)
{
    Info<< "Starting untangling the surface of the volume mesh" << endl;

    bool changed(false);

    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    boolList smoothVertex(bPoints.size(), false);
    # pragma omp parallel for schedule(dynamic, 50)
    forAll(selectedBoundaryPoints, i)
    {
        if (vertexType_[selectedBoundaryPoints[i]] & LOCKED)
            continue;

        smoothVertex[selectedBoundaryPoints[i]] = true;
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    meshSurfaceMapper* mapperPtr = nullptr;
    if (octreePtr_)
        mapperPtr = new meshSurfaceMapper(*partitionerPtr_, *octreePtr_);

    bool remapVertex(true);
    label nInvertedTria;
    label nGlobalIter(0);

    DynamicList<label> movedPoints;
    DynamicList<label> movedEdgePoints;

    label minNumInverted(bPoints.size());
    FIFOStack<label> nInvertedHistory;
    pointField minInvertedPoints(bPoints.size());

    do
    {
        label nIter(0), nAfterRefresh(0);

        do
        {
            nInvertedTria =
                findInvertedVertices(smoothVertex, nAdditionalLayers);

            if (nInvertedTria == 0)
            {
                break;
            }
            else if (enforceConstraints_ && !remapVertex)
            {
                polyMeshGen& mesh =
                    const_cast<polyMeshGen&>(surfaceEngine_.mesh());

                const label subsetId =
                    mesh.addPointSubset(badPointsSubsetName_);

                forAll(smoothVertex, bpI)
                    if (smoothVertex[bpI])
                        mesh.addPointToSubset(subsetId, bPoints[bpI]);

                WarningInFunction
                    << "Writing mesh with " << badPointsSubsetName_
                    << " subset. These points cannot be untangled"
                    << " without sacrificing geometry constraints. Exitting.."
                    << endl;

                throw std::logic_error
                (
                    "bool meshSurfaceOptimizer::untangleSurface"
                    "(const DynamicList<label>&, const label)"
                    "Cannot untangle mesh!!"
                );
            }

            // find the min number of inverted points and
            // add the last number to the stack
            if (nInvertedTria < minNumInverted)
            {
                minNumInverted = nInvertedTria;
                nAfterRefresh = 0;

                # pragma omp parallel for schedule(dynamic, 100)
                forAll(bPoints, bpI)
                    minInvertedPoints[bpI] = points[bPoints[bpI]];
            }

            // count the number of iteration after the last minimum occurence
            ++nAfterRefresh;

            // update the stack
            nInvertedHistory.push(nInvertedTria);
            if (nInvertedHistory.size() > 2)
                nInvertedHistory.pop();

            // check if the number of inverted points reduces
            bool minimumInStack(false);
            forAllConstIter(FIFOStack<label>, nInvertedHistory, it)
            {
                if (it() == minNumInverted)
                {
                    minimumInStack = true;
                }
            }

            // stop if the procedure does not minimise
            // the number of inverted points
            if (!minimumInStack || (nAfterRefresh > 2))
                break;

            // find points which will be handled by the smoothers
            changed = true;

            movedPoints.clear();
            movedEdgePoints.clear();

            forAll(bPoints, bpI)
            {
                if (!smoothVertex[bpI])
                    continue;

                if (vertexType_[bpI] & PARTITION)
                {
                    movedPoints.append(bpI);
                }
                else if (vertexType_[bpI] & EDGE)
                {
                    movedEdgePoints.append(bpI);
                }
            }

            // smooth edge vertices
            smoothEdgePoints(movedEdgePoints);
            if (remapVertex && mapperPtr)
                mapperPtr->mapEdgeNodes(movedEdgePoints);
            surfaceModifier.updateGeometry(movedEdgePoints);

            // use laplacian smoothing
            smoothLaplacianFC(movedPoints);
            surfaceModifier.updateGeometry(movedPoints);

            // use surface optimizer
            smoothSurfaceOptimizer(movedPoints);

            if (remapVertex && mapperPtr)
                mapperPtr->mapVerticesOntoSurface(movedPoints);

            // update normals and other geometric data
            surfaceModifier.updateGeometry(movedPoints);

        } while (nInvertedTria && (++nIter < 20));

        if (nInvertedTria > 0)
        {
            // use the combination with the minimum number of inverted points
            meshSurfaceEngineModifier sMod(surfaceEngine_);
            forAll(minInvertedPoints, bpI)
                sMod.moveBoundaryVertexNoUpdate(bpI, minInvertedPoints[bpI]);

            sMod.updateGeometry();
        }

        if (nInvertedTria)
        {
            Info<< "Smoothing remaining inverted vertices " << endl;

            movedPoints.clear();
            forAll(smoothVertex, bpI)
                if (smoothVertex[bpI])
                {
                    movedPoints.append(bpI);
                }

            smoothLaplacianFC(movedPoints, false);

            if (remapVertex && mapperPtr)
                mapperPtr->mapVerticesOntoSurface(movedPoints);

            // update normals and other geometric data
            surfaceModifier.updateGeometry(movedPoints);

            if (nGlobalIter > 5)
                remapVertex = false;
        }

    } while (nInvertedTria && (++nGlobalIter < 10));

    deleteDemandDrivenData(mapperPtr);

    if (nInvertedTria != 0)
    {
        // the procedure has given up without success
        // there exist some remaining inverted faces in the mesh
        polyMeshGen& mesh =
            const_cast<polyMeshGen&>(surfaceEngine_.mesh());

        label subsetId = mesh.pointSubsetIndex(badPointsSubsetName_);
        if (subsetId >= 0)
            mesh.removePointSubset(subsetId);
        subsetId = mesh.addPointSubset(badPointsSubsetName_);

        forAll(smoothVertex, bpI)
            if (smoothVertex[bpI])
                mesh.addPointToSubset(subsetId, bPoints[bpI]);
    }

    Info<< "Finished untangling the surface of the volume mesh" << endl;

    return changed;
}


bool Foam::Module::meshSurfaceOptimizer::untangleSurface
(
    const label nAdditionalLayers
)
{
    DynamicList<label> selectedPts;
    selectedPts.setSize(surfaceEngine_.boundaryPoints().size());
    forAll(selectedPts, i)
        selectedPts[i] = i;

    return untangleSurface(selectedPts, nAdditionalLayers);
}


void Foam::Module::meshSurfaceOptimizer::optimizeSurface
(
    const label nIterations
)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    // needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    meshSurfaceMapper* mapperPtr = nullptr;
    if (octreePtr_)
        mapperPtr = new meshSurfaceMapper(*partitionerPtr_, *octreePtr_);

    DynamicList<label> edgePoints, partitionPoints;
    forAll(bPoints, bpI)
    {
        if (vertexType_[bpI] & LOCKED)
            continue;

        if (vertexType_[bpI] & EDGE)
        {
            edgePoints.append(bpI);
        }
        else if (vertexType_[bpI] & PARTITION)
        {
            partitionPoints.append(bpI);
        }
    }

    // optimize edge vertices
    Info<< "Optimizing edges. Iteration:" << flush;
    for (label i = 0; i < nIterations; ++i)
    {
        Info<< "." << flush;

        meshSurfaceEngineModifier bMod(surfaceEngine_);

        smoothEdgePoints(edgePoints);

        // project vertices back onto the boundary
        if (mapperPtr)
            mapperPtr->mapEdgeNodes(edgePoints);

        // update the geometry information
        bMod.updateGeometry(edgePoints);
    }
    Info<< endl;

    // delete the mapper
    deleteDemandDrivenData(mapperPtr);

    // optimize positions of surface vertices which are not on surface edges
    Info<< "Optimizing surface vertices. Iteration:";
    for (label i = 0; i < nIterations; ++i)
    {
        smoothLaplacianFC(partitionPoints, true);

        smoothSurfaceOptimizer(partitionPoints);

        Info<< "." << flush;
    }

    Info<< endl;

    untangleSurface(0);
}


void Foam::Module::meshSurfaceOptimizer::optimizeSurface2D
(
    const label nIterations
)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeList& edges = surfaceEngine_.edges();
    const labelList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine
    (
        const_cast<polyMeshGen&>(surfaceEngine_.mesh())
    );
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();

    // needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    DynamicList<label> movedPoints, activeEdges, updatePoints;
    forAll(edges, beI)
    {
        const edge& e = edges[beI];

        if (zMinPoint[e.start()] ^ zMinPoint[e.end()])
        {
            label bpI = bp[e.start()];
            if (!zMinPoint[e.start()])
                bpI = bp[e.end()];

            if (vertexType_[bpI] & EDGE)
            {
                activeEdges.append(beI);

                updatePoints.append(bp[e.start()]);
                updatePoints.append(bp[e.end()]);

                movedPoints.append(bpI);
            }
        }
    }

    meshSurfaceMapper2D* mapperPtr = nullptr;
    if (octreePtr_)
        mapperPtr = new meshSurfaceMapper2D(surfaceEngine_, *octreePtr_);

    // optimize edge vertices
    meshSurfaceEngineModifier bMod(surfaceEngine_);

    Info<< "Optimizing edges. Iteration:" << flush;
    for (label i = 0; i < nIterations; ++i)
    {
        Info<< "." << flush;

        smoothEdgePoints(movedPoints);

        // move points with maximum z coordinate
        mesh2DEngine.correctPoints();

        // map boundary edges to the surface
        mapperPtr->mapVerticesOntoSurfacePatches(activeEdges);

        // update normal, centres, etc, after the surface has been modified
        bMod.updateGeometry(updatePoints);
    }
    Info<< endl;

    // optimize Pts of surface vertices which are not on surface edges
    movedPoints.clear();
    forAll(bPoints, bpI)
        if (zMinPoint[bPoints[bpI]] && (vertexType_[bpI] & PARTITION))
        {
            movedPoints.append(bpI);
        }
    Info<< "Optimizing surface vertices. Iteration:";
    for (label i = 0; i < nIterations; ++i)
    {
        Info<< "." << flush;

        smoothSurfaceOptimizer(movedPoints);

        // move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        // update geometrical data due to movement of vertices
        bMod.updateGeometry();
    }

    Info<< endl;

    deleteDemandDrivenData(mapperPtr);
}


void Foam::Module::meshSurfaceOptimizer::untangleSurface2D()
{
    const polyMeshGen& mesh = surfaceEngine_.mesh();
    const faceListPMG& faces = mesh.faces();
    const labelListList2D& pointFaces = mesh.addressingData().pointFaces();

    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const labelList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine(const_cast<polyMeshGen&>(mesh));
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& activeFace = mesh2DEngine.activeFace();

    // needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    boolList activeBoundaryPoint(bPoints.size());
    boolList changedFace(activeFace.size(), true);

    label iterationI(0);
    do
    {
        labelHashSet badFaces;
        const label nBadFaces =
            polyMeshGenChecks::findBadFaces
            (
                mesh,
                badFaces,
                false,
                &changedFace
            );

        Info<< "Iteration " << iterationI
             << ". Number of bad faces " << nBadFaces << endl;

        if (nBadFaces == 0)
            break;

        // update active points and faces affected by the movement
        // of active points
        activeBoundaryPoint = false;
        changedFace = false;
        for (const label facei : badFaces.toc())
        {
            const face& f = faces[facei];

            forAll(f, pI)
            {
                if (zMinPoint[f[pI]])
                {
                    activeBoundaryPoint[bp[f[pI]]] = true;

                    forAllRow(pointFaces, f[pI], pfI)
                        changedFace[pointFaces[f[pI]][pfI]] = true;
                }
            }
        }

        // apply smoothing to the activated points
        meshSurfaceEngineModifier bMod(surfaceEngine_);

        DynamicList<label> movedPts, edgePts;
        forAll(bPoints, bpI)
        {
            if (!activeBoundaryPoint[bpI])
                continue;

            if (vertexType_[bpI] & EDGE)
            {
                edgePts.append(bpI);
            }
            else if (vertexType_[bpI] & PARTITION)
            {
                movedPts.append(bpI);
            }
        }

        for (label i = 0; i < 5; ++i)
        {
            smoothEdgePoints(edgePts);

            bMod.updateGeometry(edgePts);

            smoothSurfaceOptimizer(movedPts);

            bMod.updateGeometry(movedPts);
        }

        // move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        // update geometrical data due to movement of vertices
        bMod.updateGeometry();

        // update cell centres and face centres
        const_cast<polyMeshGenAddressing&>
        (
            mesh.addressingData()
        ).updateGeometry(changedFace);

    } while (++iterationI < 20);

    // delete invalid data
    mesh.clearAddressingData();
}


// ************************************************************************* //
