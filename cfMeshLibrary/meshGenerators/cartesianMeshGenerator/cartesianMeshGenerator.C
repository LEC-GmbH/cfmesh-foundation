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

        if (checkCellConnectionsOverFaces(mesh_).checkCellGroups())
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
    // add boundary layers
    boundaryLayers bl(mesh_);
    bl.addLayerForAllPatches();
}


void Foam::Module::cartesianMeshGenerator::refBoundaryLayers()
{
    if (meshDict_.isDict("boundaryLayers"))
    {
        refineBoundaryLayers refLayers(mesh_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

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

    extractPatches();

    mapEdgesAndCorners();

    optimiseMeshSurface();
}


void Foam::Module::cartesianMeshGenerator::addLayers()
{
    generateBoundaryLayers();

    optimiseMesh(true);

    refBoundaryLayers();
}


void Foam::Module::cartesianMeshGenerator::optimiseMesh(const bool layersExist)
{
    optimiseFinalMesh(layersExist);

    projectSurfaceAfterBackScaling(layersExist);
}


void Foam::Module::cartesianMeshGenerator::finalise()
{
    renumberMesh();

    replaceBoundaries();
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
