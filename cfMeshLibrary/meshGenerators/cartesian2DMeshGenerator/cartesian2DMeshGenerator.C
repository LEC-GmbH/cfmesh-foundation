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

#include "cartesian2DMeshGenerator.H"
#include "triSurface2DCheck.H"
#include "polyMeshGen2DEngine.H"
#include "triSurf.H"
#include "triSurfacePatchManipulator.H"
#include "triSurfaceCleanupDuplicateTriangles.H"
#include "demandDrivenData.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper2D.H"
#include "meshSurfaceEdgeExtractor2D.H"
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

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
    defineTypeNameAndDebug(cartesian2DMeshGenerator, 0);
    addToRunTimeSelectionTable
    (
        meshGenerator,
        cartesian2DMeshGenerator,
        dictionary
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::Module::cartesian2DMeshGenerator::createCartesianMesh()
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


void Foam::Module::cartesian2DMeshGenerator::surfacePreparation()
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


void Foam::Module::cartesian2DMeshGenerator::mapMeshToSurface()
{
    // calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(mesh_);

    // pre-map mesh surface
    meshSurfaceMapper2D mapper(*msePtr, *octreePtr_);

    mapper.adjustZCoordinates();

    mapper.preMapVertices();

    // map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    deleteDemandDrivenData(msePtr);
}


void Foam::Module::cartesian2DMeshGenerator::extractPatches()
{
    meshSurfaceEdgeExtractor2D(mesh_, *octreePtr_).distributeBoundaryFaces();
}


void Foam::Module::cartesian2DMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractor2D(mesh_, *octreePtr_).remapBoundaryPoints();
}


void Foam::Module::cartesian2DMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer optimizer(mse, *octreePtr_);
    optimizer.optimizeSurface2D();
    optimizer.untangleSurface2D(meshDict_.subOrEmptyDict("meshQualitySettings"));
}


void Foam::Module::cartesian2DMeshGenerator::generateBoundaryLayers()
{
    boundaryLayers bl(mesh_);

    bl.activate2DMode();

    bl.addLayerForAllPatches();

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

        // delete the octree
        deleteDemandDrivenData(octreePtr_);

        // contruct a new octree from the input surface
        octreePtr_ = new meshOctree(*surfacePtr_, true);
        meshOctreeCreator(*octreePtr_).createOctreeWithRefinedBoundary(20);

        mapEdgesAndCorners();

        optimiseMeshSurface();
    }
}


void Foam::Module::cartesian2DMeshGenerator::refBoundaryLayers()
{
    if (meshDict_.isDict("boundaryLayers"))
    {
        refineBoundaryLayers refLayers(mesh_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

        refLayers.activate2DMode();

        refLayers.refineLayers();

        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer optimizer(mse, *octreePtr_);

        optimizer.untangleSurface2D
        (
            meshDict_.subOrEmptyDict("meshQualitySettings")
        );
    }
}


void Foam::Module::cartesian2DMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_, true);
}


void Foam::Module::cartesian2DMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();
}


void Foam::Module::cartesian2DMeshGenerator::createCastellatedMesh()
{
    createCartesianMesh();

    surfacePreparation();
}


void Foam::Module::cartesian2DMeshGenerator::readMesh()
{
    mesh_.read();
}


void Foam::Module::cartesian2DMeshGenerator::snap()
{
    mapMeshToSurface();

    extractPatches();

    mapEdgesAndCorners();

    optimiseMeshSurface();
}


void Foam::Module::cartesian2DMeshGenerator::addLayers()
{
    generateBoundaryLayers();

    optimiseMesh(true);

    refBoundaryLayers();
}


void Foam::Module::cartesian2DMeshGenerator::optimiseMesh
(
    const bool /*layersExist*/
)
{
    // no boundary-layer-specific optimisation in this generator
    optimiseMeshSurface();
}


void Foam::Module::cartesian2DMeshGenerator::finalise()
{
    renumberMesh();

    replaceBoundaries();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::cartesian2DMeshGenerator::cartesian2DMeshGenerator
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

    if (true)
    {
        triSurface2DCheck surfCheck(*surfacePtr_);
        if (!surfCheck.is2DSurface())
        {
            surfCheck.createSubsets();

            Info<< "Writting surface with subsets to file "
                << "badSurfaceWithSubsets.fms" << endl;

            surfacePtr_->writeSurface("badSurfaceWithSubsets.fms");
        }
    }

    if (surfacePtr_->featureEdges().size() != 0)
    {
        // get rid of duplicate triangles as they cause strange problems
        triSurfaceCleanupDuplicateTriangles
        (
            const_cast<triSurf&>(*surfacePtr_)
        );

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

        octreePtr_ = new meshOctree(*modSurfacePtr_, true);
    }
    else
    {
        octreePtr_ = new meshOctree(*surfacePtr_, true);
    }

    meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::cartesian2DMeshGenerator::~cartesian2DMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(modSurfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyMesh>
Foam::Module::cartesian2DMeshGenerator::generatePolyMesh
(
    const IOobject& io
) const
{
    return polyMeshGenToFoamMesh(mesh_, io);
}


// ************************************************************************* //
