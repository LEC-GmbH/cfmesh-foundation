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

#include <stdexcept>

#include "demandDrivenData.H"
#include "meshOptimizer.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGenChecks.H"
#include "partTetMesh.H"
#include "HashSet.H"

#include "tetMeshOptimisation.H"
#include "boundaryLayerOptimisation.H"
#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOptimizer::untangleMeshFV
(
    const label maxNumGlobalIterations,
    const label maxNumIterations,
    const label maxNumSurfaceIterations,
    const bool relaxedCheck
)
{
    Info<< "Starting untangling the mesh" << endl;

    label nBadFaces, nGlobalIter(0), nIter;

    const faceListPMG& faces = mesh_.faces();

    boolList changedFace(faces.size(), true);

    // check if any points in the tet mesh shall not move
    DynamicList<label> lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
        {
            lockedPoints.append(pointI);
        }
    }

    labelHashSet badFaces;

    do
    {
        nIter = 0;

        label minNumBadFaces(10*faces.size()), minIter(-1);
        do
        {
            if (!relaxedCheck)
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFaces
                    (
                        mesh_,
                        badFaces,
                        false,
                        &changedFace,
                        qualitySettingsDict_
                    );
            }
            else
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFacesRelaxed
                    (
                        mesh_,
                        badFaces,
                        false,
                        &changedFace,
                        qualitySettingsDict_
                    );
            }

            Info<< "Iteration " << nIter
                << ". Number of bad faces is " << nBadFaces << endl;

            // perform optimisation
            if (nBadFaces == 0)
            {
                break;
            }

            if (nBadFaces < minNumBadFaces)
            {
                minNumBadFaces = nBadFaces;
                minIter = nIter;
            }

            // create a tet mesh from the mesh and the labels of bad faces
            partTetMesh tetMesh
            (
                mesh_,
                lockedPoints,
                badFaces,
                (nGlobalIter/2) + 1
            );

            // construct tetMeshOptimisation and improve positions of
            // points in the tet mesh
            tetMeshOptimisation tmo(tetMesh);

            tmo.optimiseUsingKnuppMetric();

            tmo.optimiseUsingMeshUntangler();

            tmo.optimiseUsingVolumeOptimizer();

            // update points in the mesh from the coordinates in the tet mesh
            tetMesh.updateOrigMesh(&changedFace);

        } while ((nIter < minIter + 5) && (++nIter < maxNumIterations));

        if ((nBadFaces == 0) || (++nGlobalIter >= maxNumGlobalIterations))
        {
            break;
        }

        // move boundary vertices
        nIter = 0;

        while (nIter++< maxNumSurfaceIterations)
        {
            if (!relaxedCheck)
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFaces
                    (
                        mesh_,
                        badFaces,
                        false,
                        &changedFace,
                        qualitySettingsDict_
                    );
            }
            else
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFacesRelaxed
                    (
                        mesh_,
                        badFaces,
                        false,
                        &changedFace,
                        qualitySettingsDict_
                    );
            }

            Info<< "Iteration " << nIter
                << ". Number of bad faces is " << nBadFaces << endl;

            // perform optimisation
            if (nBadFaces == 0)
            {
                break;
            }
            else if (enforceConstraints_)
            {
                const label subsetId =
                    mesh_.addPointSubset(badPointsSubsetName_);

                for (const label facei : badFaces.toc())
                {
                    const face& f = faces[facei];
                    forAll(f, pI)
                    {
                        mesh_.addPointToSubset(subsetId, f[pI]);
                    }
                }

                WarningInFunction
                    << "Writing mesh with " << badPointsSubsetName_
                    << " subset. These points cannot be untangled"
                    << " without sacrificing geometry constraints. Exitting.."
                    << endl;

                throw std::logic_error
                (
                    "void meshOptimizer::untangleMeshFV()"
                    "Cannot untangle mesh!!"
                );
            }

            // create tethrahedral mesh from the cells which shall be smoothed
            partTetMesh tetMesh(mesh_, lockedPoints, badFaces, 0);

            // contruct tetMeshOptimisation
            tetMeshOptimisation tmo(tetMesh);

            if (nGlobalIter < 2)
            {
                // the point stays in the plane determined by the point normal
                tmo.optimiseBoundaryVolumeOptimizer(true);
            }
            else if (nGlobalIter < 5)
            {
                // move points without any constraints on the movement
                tmo.optimiseBoundarySurfaceLaplace();
            }
            else
            {
                // move boundary points without any constraints
                tmo.optimiseBoundaryVolumeOptimizer(false);
            }

            tetMesh.updateOrigMesh(&changedFace);

        }

    } while (nBadFaces);

    if (nBadFaces != 0)
    {
        label subsetId = mesh_.faceSubsetIndex("badFaces");
        if (subsetId >= 0)
        {
            mesh_.removeFaceSubset(subsetId);
        }
        subsetId = mesh_.addFaceSubset("badFaces");

        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        const label badCellsId =
                mesh_.addCellSubset("badCells");

        for (const label facei : badFaces.toc())
        {
            mesh_.addFaceToSubset(subsetId, facei);
            mesh_.addCellToSubset(badCellsId, owner[facei]);
            if (neighbour[facei] < 0)
            {
                continue;
            }
            mesh_.addCellToSubset(badCellsId, neighbour[facei]);
        }
    }

    Info<< "Finished untangling the mesh" << endl;
}


void Foam::Module::meshOptimizer::optimizeBoundaryLayer
(
    const bool addBufferLayer
)
{
    if (!meshDict_.empty())
    {
        const dictionary& meshDict = meshDict_;

        bool smoothLayer(false);

        if (meshDict.found("boundaryLayers"))
        {
            const dictionary& layersDict = meshDict.subDict("boundaryLayers");

            layersDict.readIfPresent("optimiseLayer", smoothLayer);
        }

        if (!smoothLayer)
        {
            return;
        }

        if (addBufferLayer)
        {
            // create a buffer layer which will not be modified by the smoother
            refineBoundaryLayers refLayers(mesh_);

            refineBoundaryLayers::readSettings(meshDict, refLayers);

            refLayers.activateSpecialMode();

            refLayers.refineLayers();

            clearSurface();
            calculatePointLocations();
        }

        Info<< "Starting optimising boundary layer" << endl;

        const meshSurfaceEngine& mse = meshSurface();
        const labelList& faceOwner = mse.faceOwners();

        boundaryLayerOptimisation optimiser(mesh_, mse);

        boundaryLayerOptimisation::readSettings(meshDict, optimiser);

        optimiser.optimiseLayer();

        // check if the bnd layer is tangled somewhere
        DynamicList<label> bndLayerCells;
        const boolList& baseFace = optimiser.isBaseFace();

        forAll(baseFace, bfI)
        {
            if (baseFace[bfI])
            {
                bndLayerCells.append(faceOwner[bfI]);
            }
        }

        clearSurface();
        mesh_.clearAddressingData();

        // lock boundary layer points, faces and cells
        lockCells(bndLayerCells);

        // optimize mesh quality
        optimizeMeshFV(5, 1, 50, 0);

        // untangle remaining faces and lock the boundary layer cells
        untangleMeshFV(2, 50, 0);

        // unlock bnd layer points
        removeUserConstraints();

        Info<< "Finished optimising boundary layer" << endl;
    }
}


void Foam::Module::meshOptimizer::untangleBoundaryLayer()
{
    bool untangleLayer(true);
    if (!meshDict_.empty())
    {
        const dictionary& meshDict = meshDict_;

        if (meshDict.found("boundaryLayers") )
        {
            const dictionary& subdict = meshDict.subDict("boundaryLayers");

            subdict.readIfPresent("untangleLayers", untangleLayer);
        }
    }

    if (!untangleLayer)
    {
        labelHashSet badFaces;
        polyMeshGenChecks::checkFacePyramids(mesh_, false, VSMALL, &badFaces);

        const label nInvalidFaces = badFaces.size();

        if (nInvalidFaces != 0)
        {
            const labelList& owner = mesh_.owner();
            const labelList& neighbour = mesh_.neighbour();

            const label badBlCellsId =
                mesh_.addCellSubset("invalidBoundaryLayerCells");

            for (const label facei : badFaces.toc())
            {
                mesh_.addCellToSubset(badBlCellsId, owner[facei]);

                if (neighbour[facei] < 0)
                {
                    continue;
                }

                mesh_.addCellToSubset(badBlCellsId, neighbour[facei]);
            }

            throw std::logic_error
            (
                "void meshOptimizer::untangleBoundaryLayer()"
                "Found invalid faces in the boundary layer."
                " Cannot untangle mesh!!"
            );
        }
    }
    else
    {
        optimizeLowQualityFaces();
        removeUserConstraints();
        untangleMeshFV(2, 50, 1, true);
    }
}


void Foam::Module::meshOptimizer::optimizeLowQualityFaces
(
    const label maxNumIterations
)
{
    label nBadFaces, nIter(0);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    // check if any points in the tet mesh shall not move
    DynamicList<label> lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
        {
            lockedPoints.append(pointI);
        }
    }

    do
    {
        labelHashSet lowQualityFaces;
        nBadFaces =
            polyMeshGenChecks::findLowQualityFaces
            (
                mesh_,
                lowQualityFaces,
                false,
                &changedFace,
                qualitySettingsDict_
            );

        changedFace = false;
        for (const label facei : lowQualityFaces.toc())
        {
            changedFace[facei] = true;
        }

        Info<< "Iteration " << nIter
            << ". Number of bad faces is " << nBadFaces << endl;

        // perform optimisation
        if (nBadFaces == 0)
        {
            break;
        }

        partTetMesh tetMesh(mesh_, lockedPoints, lowQualityFaces, 2);

        // construct tetMeshOptimisation and improve positions
        // of points in the tet mesh
        tetMeshOptimisation tmo(tetMesh);

        tmo.optimiseUsingVolumeOptimizer();

        // update points in the mesh from the new coordinates in the tet mesh
        tetMesh.updateOrigMesh(&changedFace);

    } while (++nIter < maxNumIterations);
}


void Foam::Module::meshOptimizer::optimizeMeshNearBoundaries
(
    const label maxNumIterations,
    const label numLayersOfCells
)
{
    label nIter(0);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    // check if any points in the tet mesh shall not move
    DynamicList<label> lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
        {
            lockedPoints.append(pointI);
        }
    }

    partTetMesh tetMesh(mesh_, lockedPoints, numLayersOfCells);
    tetMeshOptimisation tmo(tetMesh);
    Info<< "Iteration:" << flush;
    do
    {
        tmo.optimiseUsingVolumeOptimizer(1);

        tetMesh.updateOrigMesh(&changedFace);

        Info<< "." << flush;

    } while (++nIter < maxNumIterations);

    Info<< endl;
}


void Foam::Module::meshOptimizer::optimizeMeshFV
(
    const label numLaplaceIterations,
    const label maxNumGlobalIterations,
    const label maxNumIterations,
    const label maxNumSurfaceIterations
)
{
    Info<< "Starting smoothing the mesh" << endl;

    laplaceSmoother lps(mesh_, vertexLocation_);
    lps.optimizeLaplacianPC(numLaplaceIterations);

    untangleMeshFV
    (
        maxNumGlobalIterations,
        maxNumIterations,
        maxNumSurfaceIterations
    );

    Info<< "Finished smoothing the mesh" << endl;
}


void Foam::Module::meshOptimizer::optimizeMeshFVBestQuality
(
    const label maxNumIterations,
    const scalar threshold
)
{
    label nBadFaces, nIter(0);
    label minIter(-1);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    // check if any points in the tet mesh shall not move
    DynamicList<label> lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
        {
            lockedPoints.append(pointI);
        }
    }

    label minNumBadFaces(10*faces.size());
    do
    {
        labelHashSet lowQualityFaces;
        nBadFaces =
            polyMeshGenChecks::findWorstQualityFaces
            (
                mesh_,
                lowQualityFaces,
                false,
                &changedFace,
                threshold
            );

        changedFace = false;
        for (const label facei : lowQualityFaces.toc())
        {
            changedFace[facei] = true;
        }

        Info<< "Iteration " << nIter
            << ". Number of worst quality faces is " << nBadFaces << endl;

        // perform optimisation
        if (nBadFaces == 0)
        {
            break;
        }

        if (nBadFaces < minNumBadFaces)
        {
            minNumBadFaces = nBadFaces;

            // update the iteration number when the minimum is achieved
            minIter = nIter;
        }

        partTetMesh tetMesh(mesh_, lockedPoints, lowQualityFaces, 2);

        // construct tetMeshOptimisation and improve positions
        // of points in the tet mesh
        tetMeshOptimisation tmo(tetMesh);

        tmo.optimiseUsingVolumeOptimizer(20);

        // update points in the mesh from the new coordinates in the tet mesh
        tetMesh.updateOrigMesh(&changedFace);

    } while ((nIter < minIter + 5) && (++nIter < maxNumIterations));
}


// ************************************************************************* //
