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

#include "cartesianMeshGenerator.H"
#include "triSurf.H"
#include "triSurfacePatchManipulator.H"
#include "demandDrivenData.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "edgeExtractor.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "checkBoundaryFacesSharingTwoEdges.H"
#include "polyMeshGenGeometryModification.H"
#include "surfaceMeshGeometryModification.H"
#include "coordinateModifier.H"
#include "polyMeshGenToFoamMesh.H"
#include "addToRunTimeSelectionTable.H"

//#define DEBUG

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
    defineTypeNameAndDebug(cartesianMeshGenerator, 0);
    addToRunTimeSelectionTable(meshGenerator, cartesianMeshGenerator, dictionary);
}
}

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void Foam::Module::cartesianMeshGenerator::createCartesianMesh()
{
    // create polyMesh from octree boxes
    cartesianMeshExtractor cme(*octreePtr_, meshDict_, mesh_);

    if
    (
        meshDict_.lookupOrDefault<bool>
        (
            "decomposePolyhedraIntoTetsAndPyrs",
            false
        )
    )
    {
        cme.decomposeSplitHexes();
    }

    cme.createMesh();
}


void Foam::Module::cartesianMeshGenerator::surfacePreparation()
{
    // removes unnecessary cells and morph the boundary
    // such that there is only one boundary face per cell
    // It also checks topology of cells after morphing is performed

    // which disconnected region to keep if the surface encloses more
    // than one - the region containing insidePoint if given (resolved
    // into the same modified space anisotropicSources may already have
    // put the mesh into), otherwise the region with the most cells
    point insidePoint;
    const bool haveInsidePoint =
        meshDict_.readIfPresent("insidePoint", insidePoint);

    if (haveInsidePoint)
    {
        const dictionary& anisotropicSourcesDict =
            meshDict_.subOrEmptyDict("anisotropicSources");

        if (!anisotropicSourcesDict.empty())
        {
            insidePoint =
                coordinateModifier(anisotropicSourcesDict)
               .modifiedPoint(insidePoint);
        }
    }

    bool changed;

    do
    {
        changed = false;

        checkIrregularSurfaceConnections checkConnections(mesh_);
        if (checkConnections.checkAndFixIrregularConnections())
        {
            changed = true;
        }

        if (checkNonMappableCellConnections(mesh_).removeCells())
        {
            changed = true;
        }

        if
        (
            checkCellConnectionsOverFaces(mesh_).checkCellGroups
            (
                haveInsidePoint,
                insidePoint
            )
        )
        {
            changed = true;
        }
    } while (changed);

    checkBoundaryFacesSharingTwoEdges(mesh_).improveTopology();
}


void Foam::Module::cartesianMeshGenerator::mapMeshToSurface()
{
    // calculate mesh surface
    meshSurfaceEngine mse(mesh_);

    // pre-map mesh surface
    meshSurfaceMapper mapper(mse, *octreePtr_);
    mapper.preMapVertices();

    // map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    // untangle surface faces
    meshSurfaceOptimizer(mse, *octreePtr_).untangleSurface();
}


void Foam::Module::cartesianMeshGenerator::extractPatches()
{
    edgeExtractor extractor(mesh_, *octreePtr_);

    Info<< "Extracting edges" << endl;
    extractor.extractEdges();

    extractor.updateMeshPatches();
}


void Foam::Module::cartesianMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);
}


void Foam::Module::cartesianMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();
}


void Foam::Module::cartesianMeshGenerator::generateBoundaryLayers()
{
    // insert a single geometry-adapted (non-constant width) buffer layer
    // for every patch - unconditional, regardless of whether the user
    // configured "real" (multi-cell, thicknessRatio-graded) boundary
    // layers at all; refBoundaryLayers() only subdivides this same base
    // layer further when boundaryLayers is actually configured. Runs at
    // the end of snap() (not the start of addLayers()) so a snap-only
    // run - or the snap half of a staged snap/addLayers split - gets
    // this quality benefit too, not just full/addLayers-enabled runs
    boundaryLayers bl(mesh_);
    bl.addLayerForAllPatches();
}


void Foam::Module::cartesianMeshGenerator::refBoundaryLayers()
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


void Foam::Module::cartesianMeshGenerator::optimiseFinalMesh
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
        optimizer.optimizeBoundaryLayer(modSurfacePtr_ == nullptr);
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


void Foam::Module::cartesianMeshGenerator::projectSurfaceAfterBackScaling
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


void Foam::Module::cartesianMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_);
}


void Foam::Module::cartesianMeshGenerator::renamePatches()
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


void Foam::Module::cartesianMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();
}


void Foam::Module::cartesianMeshGenerator::createCastellatedMesh()
{
    createCartesianMesh();

    surfacePreparation();
}


void Foam::Module::cartesianMeshGenerator::readMesh()
{
    mesh_.read();
}


void Foam::Module::cartesianMeshGenerator::snap()
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


void Foam::Module::cartesianMeshGenerator::addLayers()
{
    // assumes snap() already ran (this invocation or an earlier one) and
    // left an optimised buffer layer in place - matching snappyHexMesh's
    // own addLayers stage, which never compensates for a skipped/deferred
    // snap either; running addLayers without snap having run is on the
    // caller, not something to special-case here
    refBoundaryLayers();

    renamePatches();
}


void Foam::Module::cartesianMeshGenerator::optimiseMesh(const bool layersExist)
{
    optimiseFinalMesh(layersExist);

    projectSurfaceAfterBackScaling(layersExist);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::cartesianMeshGenerator::cartesianMeshGenerator
(
    const Time& time,
    IOdictionary& meshDict,
    const word& regionName,
    const fileName& meshPath
)
:
    db_(time),
    surfacePtr_(nullptr),
    modSurfacePtr_(nullptr),
    meshDict_(meshDict),
    octreePtr_(nullptr),
    mesh_(time)
{
    mesh_.setMeshLocation(regionName, meshPath);

    fileName surfaceFile(meshDict_.lookup("surfaceFile"));

    surfacePtr_ = new triSurf(db_.path()/surfaceFile);

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

Foam::Module::cartesianMeshGenerator::~cartesianMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(modSurfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMesh>
Foam::Module::cartesianMeshGenerator::generatePolyMesh
(
    const IOobject& io
) const
{
    return polyMeshGenToFoamMesh(mesh_, io);
}


// ************************************************************************* //
