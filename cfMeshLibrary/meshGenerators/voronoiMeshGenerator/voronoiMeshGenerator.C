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

#include "voronoiMeshGenerator.H"
#include "triSurf.H"
#include "demandDrivenData.H"
#include "meshOctreeCreator.H"
#include "voronoiMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "surfaceMorpherCells.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "triSurfacePatchManipulator.H"
#include "surfaceMeshGeometryModification.H"
#include "polyMeshGenGeometryModification.H"
#include "edgeExtractor.H"
#include "meshSurfaceEdgeExtractorFUN.H"
#include "polyMeshGenToFoamMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
    defineTypeNameAndDebug(voronoiMeshGenerator, 0);
    addToRunTimeSelectionTable(meshGenerator, voronoiMeshGenerator, dictionary);
}
}

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void Foam::Module::voronoiMeshGenerator::createVoronoiMesh()
{
    // create voronoi mesh from octree and Delaunay tets
    voronoiMeshExtractor vme(*octreePtr_, meshDict_, mesh_);

    vme.createMesh();

}


void Foam::Module::voronoiMeshGenerator::surfacePreparation()
{
    // removes unnecessary cells and morphs the boundary
    // such that there exists only one boundary face per cell

    surfaceMorpherCells* cmPtr = new surfaceMorpherCells(mesh_);
    cmPtr->morphMesh();
    deleteDemandDrivenData(cmPtr);

}


void Foam::Module::voronoiMeshGenerator::mapMeshToSurface()
{
    // calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

    // map mesh surface on the geometry surface
    meshSurfaceMapper mapper(*msePtr, *octreePtr_);
    mapper.preMapVertices();
    mapper.mapVerticesOntoSurface();

    // untangle surface faces
    meshSurfaceOptimizer(*msePtr, *octreePtr_).untangleSurface();

    deleteDemandDrivenData(msePtr);
}


void Foam::Module::voronoiMeshGenerator::extractPatches()
{
    edgeExtractor extractor(mesh_, *octreePtr_);

    Info<< "Extracting edges" << endl;
    extractor.extractEdges();

    extractor.updateMeshPatches();
}


void Foam::Module::voronoiMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractorFUN(mesh_, *octreePtr_, false);

}


void Foam::Module::voronoiMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer surfOptimiser(mse, *octreePtr_);
    surfOptimiser.optimizeSurface();
    surfOptimiser.untangleSurface();

}


void Foam::Module::voronoiMeshGenerator::generateBoundaryLayers()
{
    if (meshDict_.found("boundaryLayers"))
    {
        boundaryLayers bl(mesh_);

        const dictionary& bndLayers = meshDict_.subDict("boundaryLayers");

        label nLayers;
        if (bndLayers.readIfPresent("nLayers", nLayers))
        {
            if (nLayers > 0)
            {
                bl.addLayerForAllPatches();
            }
        }
        else if (bndLayers.found("patchBoundaryLayers"))
        {
            const dictionary& patchLayers =
                bndLayers.subDict("patchBoundaryLayers");

            for (const word& patchName : patchLayers.toc())
            {
                bl.addLayerForPatch(patchName);
            }
        }
    }

}


void Foam::Module::voronoiMeshGenerator::refBoundaryLayers()
{
    if (meshDict_.isDict("boundaryLayers"))
    {
        refineBoundaryLayers refLayers
        (
            mesh_,
            meshDict_.subDict("boundaryLayers")
        );

        refLayers.refineLayers();

        DynamicList<label> pointsInLayer;
        refLayers.pointsInBndLayer(pointsInLayer);

        meshOptimizer mOpt(mesh_, meshDict_);
        mOpt.lockPoints(pointsInLayer);
        mOpt.untangleBoundaryLayer();
    }
}


void Foam::Module::voronoiMeshGenerator::optimiseFinalMesh
(
    const bool layersExist
)
{
    // untangle the surface if needed
    const bool enforceConstraints =
        meshDict_.lookupOrDefault<bool>("enforceGeometryConstraints", false);

    if (true)
    {
        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer surfOpt(mse, *octreePtr_);

        if (enforceConstraints)
        {
            surfOpt.enforceConstraints();
        }

        surfOpt.optimizeSurface();
    }

    deleteDemandDrivenData(octreePtr_);

    // final optimisation
    meshOptimizer optimizer(mesh_, meshDict_);
    if (enforceConstraints)
    {
        optimizer.enforceConstraints();
    }
    optimizer.optimizeMeshFV();

    optimizer.optimizeLowQualityFaces();
    if (layersExist)
    {
        optimizer.optimizeBoundaryLayer(false);
    }
    optimizer.untangleMeshFV();

    mesh_.clearAddressingData();

    if (modSurfacePtr_)
    {
        polyMeshGenGeometryModification meshMod
        (
            mesh_,
            meshDict_.subOrEmptyDict("anisotropicSources")
        );

        // revert the mesh into the original space
        meshMod.revertGeometryModification();

        // delete modified surface mesh
        deleteDemandDrivenData(modSurfacePtr_);
    }

}


void Foam::Module::voronoiMeshGenerator::projectSurfaceAfterBackScaling
(
    const bool layersExist
)
{
    if (!meshDict_.found("anisotropicSources"))
    {
        return;
    }

    deleteDemandDrivenData(octreePtr_);
    octreePtr_ = new meshOctree(*surfacePtr_);

    meshOctreeCreator
    (
        *octreePtr_,
        meshDict_
    ).createOctreeWithRefinedBoundary(20, 30);

    // calculate mesh surface
    meshSurfaceEngine mse(mesh_);

    // pre-map mesh surface
    meshSurfaceMapper mapper(mse, *octreePtr_);

    // map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    optimiseFinalMesh(layersExist);
}


void Foam::Module::voronoiMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_);
}


void Foam::Module::voronoiMeshGenerator::renamePatches()
{
    // the user's renameBoundary rename/merge, if configured, runs here
    // (end of addLayers()) rather than in replaceBoundaries() (end of
    // snap()) because boundaryLayers.patchBoundaryLayers - consumed by
    // refBoundaryLayers(), which addLayers() has just called - resolves
    // per-patch layer settings by patch name; merging patches together
    // any earlier would erase the very per-region distinction that
    // config depends on
    renameBoundaryPatches rbp(mesh_, meshDict_);
}


void Foam::Module::voronoiMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();
}


void Foam::Module::voronoiMeshGenerator::createCastellatedMesh()
{
    createVoronoiMesh();

    surfacePreparation();
}


void Foam::Module::voronoiMeshGenerator::readMesh()
{
    mesh_.read();
}


void Foam::Module::voronoiMeshGenerator::snap()
{
    mapMeshToSurface();

    // this is where boundary faces actually get assigned to their real
    // (possibly feature-edge-split) patches - createCastellatedMesh()'s
    // own output has no such assignment yet (every boundary face is
    // still in a single default patch), so collapsing split patches
    // back to their pre-split names has to happen after this, not
    // inside createCastellatedMesh()
    extractPatches();

    mapEdgesAndCorners();

    optimiseMeshSurface();

    generateBoundaryLayers();

    // optimise the freshly-inserted buffer layer now, so snap()'s own
    // output is already good quality - matching snappyHexMesh's own
    // snap/addLayers split, where a stage's written mesh is always
    // fully valid, never "needs a later stage to become good quality"
    optimiseMesh(true);

    replaceBoundaries();
}


void Foam::Module::voronoiMeshGenerator::addLayers()
{
    // assumes snap() already ran (this invocation or an earlier one) and
    // left an optimised buffer layer in place - matching snappyHexMesh's
    // own addLayers stage, which never compensates for a skipped/deferred
    // snap either; running addLayers without snap having run is on the
    // caller, not something to special-case here
    refBoundaryLayers();

    renamePatches();
}


void Foam::Module::voronoiMeshGenerator::optimiseMesh(const bool layersExist)
{
    optimiseFinalMesh(layersExist);

    projectSurfaceAfterBackScaling(layersExist);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::voronoiMeshGenerator::voronoiMeshGenerator
(
    const Time& time,
    IOdictionary& meshDict,
    const word& regionName,
    const fileName& meshPath
)
:
    runTime_(time),
    surfacePtr_(nullptr),
    modSurfacePtr_(nullptr),
    octreePtr_(nullptr),
    pointRegionsPtr_(nullptr),
    meshDict_(meshDict),
    mesh_(time)
{
    mesh_.setMeshLocation(regionName, meshPath);

    const fileName surfaceFile(meshDict_.lookup("surfaceFile"));

    surfacePtr_ = new triSurf(runTime_.path()/surfaceFile);

    if (surfacePtr_->featureEdges().size() != 0)
    {
        // create surface patches based on the feature edges
        // and update the meshDict based on the given data
        triSurfacePatchManipulator manipulator(*surfacePtr_);

        const triSurf* surfaceWithPatches =
            manipulator.surfaceWithPatches();

        // delete the old surface and assign the new one
        deleteDemandDrivenData(surfacePtr_);
        surfacePtr_ = surfaceWithPatches;
    }

    if (meshDict_.found("anisotropicSources"))
    {
        surfaceMeshGeometryModification surfMod
        (
            *surfacePtr_,
            meshDict_.subOrEmptyDict("anisotropicSources")
        );

        modSurfacePtr_ = surfMod.modifyGeometry();

        octreePtr_ = new meshOctree(*modSurfacePtr_);
    }
    else
    {
        octreePtr_ = new meshOctree(*surfacePtr_);
    }

    meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::voronoiMeshGenerator::~voronoiMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(modSurfacePtr_);
    deleteDemandDrivenData(octreePtr_);
    deleteDemandDrivenData(pointRegionsPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMesh>
Foam::Module::voronoiMeshGenerator::generatePolyMesh
(
    const IOobject& io
) const
{
    return polyMeshGenToFoamMesh(mesh_, io);
}


// ************************************************************************* //
