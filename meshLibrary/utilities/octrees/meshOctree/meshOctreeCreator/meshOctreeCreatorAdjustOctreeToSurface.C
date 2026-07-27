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

#include "meshOctreeCreator.H"
#include "triSurf.H"
#include "edgeMesh.H"
#include "boundBox.H"
#include "demandDrivenData.H"
#include "objectRefinementList.H"
#include "meshOctreeModifier.H"
#include "helperFunctions.H"
#include "HashSet.H"

#include "coordinateModifier.H"
#include "surfaceMeshGeometryModification.H"
#include "edgeMeshGeometryModification.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeCreator::refineBoundary()
{
    meshOctreeModifier octreeModifier(octree_);
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    // refine DATA boxes to the given level
    Info<< "Refining boundary boxes to the given size" << endl;

    bool changed;
    do
    {
        changed = false;

        labelList refineCubes(leaves.size(), 0);
        scalarList rThickness(leaves.size(), 0.0);
        bool useNLayers(false);

        // select boxes which need to be refined
        # pragma omp parallel for schedule(dynamic, 100)
        forAll(leaves, leafI)
        {
            const meshOctreeCube& oc = *leaves[leafI];

            if (oc.hasContainedElements())
            {
                const label elRowI = oc.containedElements();
                const labelListList2D& containedTriangles =
                    oc.slotPtr()->containedTriangles_;

                bool refine(false);
                direction minRequestedLevel(255);
                scalar maxThicknessForLevel(0.0);
                forAllRow(containedTriangles, elRowI, tI)
                {
                    const label triI = containedTriangles[elRowI][tI];
                    DynamicList<std::pair<direction, scalar>>& refRequests =
                        surfRefLevel_[triI];

                    forAllReverse(refRequests, i)
                    {
                        const std::pair<direction, scalar>& rp = refRequests[i];
                        if (rp.first <= oc.level())
                        {
                            continue;
                        }

                        if (rp.first < minRequestedLevel)
                        {
                            refine = true;
                            minRequestedLevel = rp.first;
                            maxThicknessForLevel = rp.second;
                        }
                        else if
                        (
                            (rp.first == minRequestedLevel) &&
                            (rp.second > maxThicknessForLevel)
                        )
                        {
                            maxThicknessForLevel = rp.second;
                        }
                    }
                }

                if (refine)
                {
                    refineCubes[leafI] = 1;

                    if (maxThicknessForLevel > VSMALL)
                    {
                        rThickness[leafI] = maxThicknessForLevel;
                        useNLayers = true;
                    }

                    changed = true;
                }
            }
        }

        // mark additional boxes for refinement to achieve
        // correct refinement distance
        if (useNLayers && changed)
        {
            octreeModifier.refineSelectedBoxesAndAdditionalLayers
            (
                refineCubes,
                rThickness
            );
        }
        else if (changed)
        {
            // refine boxes
            octreeModifier.refineSelectedBoxes(refineCubes, hexRefinement_);
        }

    } while (changed);

    Info<< "Finished refining boundary boxes" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshOctreeCreator::refineBoxesContainedInObjects()
{
    if (!meshDictPtr_ || !meshDictPtr_->found("objectRefinements"))
    {
        return;
    }

    Info<< "Refining boxes inside objects" << endl;
    objectRefinementList refObjects;

    // Read objects
    if (meshDictPtr_->isDict("objectRefinements"))
    {
        const dictionary& dict = meshDictPtr_->subDict("objectRefinements");
        const wordList objectNames = dict.toc();

        refObjects.setSize(objectNames.size());

        forAll(refObjects, objectI)
        {
            const entry& objectEntry =
                dict.lookupEntry(objectNames[objectI], false, false);

            refObjects.set
            (
                objectI,
                objectRefinement::New
                (
                    objectEntry.keyword(),
                    objectEntry.dict()
                )
            );
        }
    }
    else
    {
        Istream& is = meshDictPtr_->lookup("objectRefinements");

        PtrList<entry> objectEntries(is);
        refObjects.setSize(objectEntries.size());

        forAll(refObjects, objectI)
        {
            refObjects.set
            (
                objectI,
                objectRefinement::New
                (
                    objectEntries[objectI].keyword(),
                    objectEntries[objectI].dict()
                )
            );
        }

        objectEntries.clear();
    }

    coordinateModifier* cModPtr = nullptr;

    if (meshDictPtr_->found("anisotropicSources"))
    {
        cModPtr =
            new coordinateModifier
            (
                meshDictPtr_->subDict("anisotropicSources")
            );
    }

    // calculate refinement levels
    const scalar s(readScalar(meshDictPtr_->lookup("maxCellSize")));

    List<direction> refLevels(refObjects.size(), globalRefLevel_);

    forAll(refObjects, oI)
    {
        // calculate additional refinement levels from cell size
        refObjects[oI].calculateAdditionalRefLevels(s);

        refLevels[oI] += refObjects[oI].additionalRefinementLevels();
    }

    // read refinement thickness
    scalarList refThickness(refObjects.size(), 0.0);

    forAll(refThickness, oI)
    {
        refThickness[oI] = refObjects[oI].refinementThickness();
    }

    forAll(refLevels, i)
    {
        Info<< "Ref level for object " << refObjects[i].name()
            << " is " << label(refLevels[i])
            << " refinement thickness " << refThickness[i] << endl;
    }

    // start refining boxes inside the objects
    const boundBox& rootBox = octree_.rootBox();
    meshOctreeModifier octreeModifier(octree_);
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    bool changed;
    do
    {
        changed = false;

        labelList refineCubes(leaves.size(), 0);
        scalarList rThickness(leaves.size(), 0.0);
        bool useNLayers(false);

        // select boxes which need to be refined
        # pragma omp parallel for if (leaves.size() > 1000) \
        schedule(dynamic, 20)
        forAll(leaves, leafI)
        {
            const meshOctreeCube& oc = *leaves[leafI];

            if (oc.cubeType() & meshOctreeCubeBasic::OUTSIDE)
                continue;

            boundBox bb;
            oc.cubeBox(rootBox, bb.min(), bb.max());

            if (cModPtr)
            {
                // transform bounding boxes into non-modified coordinates
                bb.min() = cModPtr->backwardModifiedPoint(bb.min());
                bb.max() = cModPtr->backwardModifiedPoint(bb.max());
            }

            bool refine(false);
            forAll(refObjects, oI)
            {
                if (refObjects[oI].intersectsObject(bb))
                {
                    if (oc.level() < refLevels[oI])
                        refine = true;

                    if (refThickness[oI] > VSMALL)
                    {
                        rThickness[leafI] =
                            Foam::max(rThickness[leafI], refThickness[oI]);

                        useNLayers = true;
                    }
                }
            }

            if (refine)
            {
                refineCubes[leafI] = 1;
                changed = true;
            }
        }

        // mark additional boxes for refinement to achieve
        // correct refinement distance
        if (useNLayers && changed)
        {
            octreeModifier.refineSelectedBoxesAndAdditionalLayers
            (
                refineCubes,
                rThickness
            );
        }
        else if (changed)
        {
            // refine boxes
            octreeModifier.refineSelectedBoxes(refineCubes, hexRefinement_);
        }

    } while (changed);

    // delete coordinate modifier if it exists
    deleteDemandDrivenData(cModPtr);

    Info<< "Finished refinement of boxes inside objects" << endl;

    // set up inside-outside information
    createInsideOutsideInformation();
}


void Foam::Module::meshOctreeCreator::refineBoxesIntersectingSurfaces()
{
    if (!meshDictPtr_ || !meshDictPtr_->found("surfaceMeshRefinement"))
    {
        return;
    }

    Info<< "Refining boxes intersecting surface meshes" << endl;

    label nMarked;

    // read surface meshes and calculate the refinement level for each
    // surface mesh
    const dictionary& surfDict = meshDictPtr_->subDict("surfaceMeshRefinement");
    const wordList surfaces = surfDict.toc();
    PtrList<const triSurf> surfaceMeshesPtr(surfaces.size());
    List<direction> refLevels(surfaces.size(), globalRefLevel_);
    scalarList refThickness(surfaces.size(), 0.0);

    // load surface meshes into memory
    forAll(surfaces, surfI)
    {
        const dictionary& dict = surfDict.subDict(surfaces[surfI]);

        const fileName fName(dict.lookup("surfaceFile"));

        if (meshDictPtr_->found("anisotropicSources"))
        {
            triSurf origSurf(fName);
            surfaceMeshGeometryModification surfMod(origSurf, *meshDictPtr_);

            surfaceMeshesPtr.set(surfI, surfMod.modifyGeometry());
        }
        else
        {
            surfaceMeshesPtr.set(surfI, new triSurf(fName));
        }

        direction addLevel(0);
        if (dict.found("cellSize"))
        {
            scalar s(readScalar(meshDictPtr_->lookup("maxCellSize")));

            const scalar cs = readScalar(dict.lookup("cellSize"));

            do
            {
                nMarked = 0;
                if (cs <= s*(1.+SMALL))
                {
                    ++nMarked;
                    ++addLevel;
                }

                s /= 2.0;

            } while (nMarked != 0);
        }
        else
        {
            label addLevelTmp(addLevel);
            dict.readIfPresent("additionalRefinementLevels", addLevelTmp);
            addLevel = direction(addLevelTmp);
        }

        dict.readIfPresent("refinementThickness", refThickness[surfI]);

        // set the refinement level for the current surface
        refLevels[surfI] += addLevel;
    }

    // start refining boxes intersecting triangles in each refinement surface
    const boundBox& rootBox = octree_.rootBox();
    const vector tol = SMALL*rootBox.span();
    meshOctreeModifier octreeModifier(octree_);
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();
    DynamicList<label> leavesInBox;

    bool changed;
    do
    {
        changed = false;

        labelList refineCubes(leaves.size(), 0);
        labelList nLayers(leaves.size(), 0);
        scalarField rThickness(leaves.size(), 0.0);
        bool useNLayers(false);

        // select boxes that need to be refined
        forAll(surfaceMeshesPtr, surfI)
        {
            const triSurf& surf = surfaceMeshesPtr[surfI];
            const pointField& points = surf.points();

            const direction surfLevel = refLevels[surfI];
            const scalar sThickness = refThickness[surfI];

            # pragma omp parallel for \
            reduction(+ : nMarked) schedule(dynamic, 10) private(leavesInBox)
            forAll(surf, triI)
            {
                // find the bounding box of the current triangle
                const labelledTri& tri = surf[triI];

                boundBox triBB(points[tri[0]], points[tri[0]]);
                for (label pI = 1; pI < 3; ++pI)
                {
                    triBB.min() = Foam::min(triBB.min(), points[tri[pI]]);
                    triBB.max() = Foam::max(triBB.max(), points[tri[pI]]);
                }

                triBB.min() -= tol;
                triBB.max() += tol;

                // find octree leaves inside the bounding box
                leavesInBox.clear();
                octree_.findLeavesContainedInBox(triBB, leavesInBox);

                // check which of the leaves are intersected by the triangle
                forAll(leavesInBox, i)
                {
                    const label leafI = leavesInBox[i];

                    const meshOctreeCube& oc = *leaves[leafI];

                    if (oc.intersectsTriangleExact(surf, rootBox, triI))
                    {

                        if (oc.level() < surfLevel)
                        {
                            // mark boxes for refinement
                            changed = true;
                            refineCubes[leafI] = 1;
                        }

                        if (sThickness > VSMALL)
                        {
                            useNLayers = true;

                            const scalar cs = oc.size(rootBox);
                            const label numLayers =
                                ceil(sThickness/cs);

                            nLayers[leafI] =
                                Foam::max
                                (
                                    nLayers[leafI],
                                    Foam::max(numLayers, 1)
                                );

                            rThickness[leafI] =
                                max(rThickness[leafI], sThickness);
                        }
                    }
                }
            }
        }

        // mark additional boxes for refinement to achieve
        // correct refinement distance

        if (useNLayers && changed)
        {
            octreeModifier.refineSelectedBoxesAndAdditionalLayers
            (
                refineCubes,
                rThickness
            );
        }
        else if (changed)
        {
            // refine boxes
            octreeModifier.refineSelectedBoxes(refineCubes, hexRefinement_);
        }

    } while (changed);

    Info<< "Finished refinement of boxes intersecting surface meshes" << endl;
}


void Foam::Module::meshOctreeCreator::refineBoxesIntersectingEdgeMeshes()
{
    if (!meshDictPtr_ || !meshDictPtr_->found("edgeMeshRefinement"))
    {
        return;
    }

    Info<< "Refining boxes intersecting edge meshes" << endl;

    // read edge meshes and calculate the refinement level for each
    // edge mesh
    const dictionary& edgeDict = meshDictPtr_->subDict("edgeMeshRefinement");
    const wordList edgeMeshNames = edgeDict.toc();

    PtrList<const edgeMesh> edgeMeshesPtr(edgeMeshNames.size());
    List<direction> refLevels(edgeMeshNames.size(), globalRefLevel_);
    scalarList refThickness(edgeMeshNames.size(), 0.0);

    // load edge meshes into memory
    forAll(edgeMeshNames, emI)
    {
        const dictionary& dict = edgeDict.subDict(edgeMeshNames[emI]);

        const fileName fName(dict.lookup("edgeFile"));

        if (meshDictPtr_->found("anisotropicSources"))
        {
            edgeMesh origEdgeMesh(fName);
            edgeMeshGeometryModification edgeMod(origEdgeMesh, *meshDictPtr_);

            edgeMeshesPtr.set(emI, edgeMod.modifyGeometry());
        }
        else
        {
            edgeMeshesPtr.set(emI, new edgeMesh(fName));
        }

        label nMarked;
        direction addLevel(0);
        if (dict.found("cellSize"))
        {
            scalar s(readScalar(meshDictPtr_->lookup("maxCellSize")));

            const scalar cs = readScalar(dict.lookup("cellSize"));

            do
            {
                nMarked = 0;
                if (cs <= s*(1.+SMALL))
                {
                    ++nMarked;
                    ++addLevel;
                }

                s /= 2.0;

            } while (nMarked != 0);
        }
        else
        {
            label addLevelTmp(addLevel);
            dict.readIfPresent("additionalRefinementLevels", addLevelTmp);
            addLevel = direction(addLevelTmp);
        }

        dict.readIfPresent("refinementThickness", refThickness[emI]);

        // set the refinement level for the current mesh
        refLevels[emI] += addLevel;
    }

    // start refining boxes intersecting triangles in each refinement surface
    const boundBox& rootBox = octree_.rootBox();
    const vector tol = SMALL*rootBox.span();
    meshOctreeModifier octreeModifier(octree_);
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();
    DynamicList<label> leavesInBox, intersectedLeaves;

    bool changed;
    do
    {
        changed = false;

        labelList refineCubes(leaves.size(), 0);
        scalarList rThickness(leaves.size(), 0.0);
        bool useNLayers(false);

        // select boxes which need to be refined
        forAll(edgeMeshNames, emI)
        {
            const edgeMesh& edgeMesh = edgeMeshesPtr[emI];
            const pointField& points = edgeMesh.points();
            const edgeList& edges = edgeMesh.edges();

            # pragma omp parallel for schedule(dynamic, 10) \
            private(leavesInBox, intersectedLeaves)
            forAll(edges, edgeI)
            {
                // find the bounding box of the current triangle
                const edge& e = edges[edgeI];

                const point& sp = points[e.start()];
                const point& ep = points[e.end()];

                boundBox edgeBB(min(sp, ep), max(sp, ep));

                edgeBB.min() -= tol;
                edgeBB.max() += tol;

                // find octree leaves inside the bounding box
                leavesInBox.clear();
                octree_.findLeavesContainedInBox(edgeBB, leavesInBox);

                // check which of the leaves are intersected by the triangle
                intersectedLeaves.clear();
                forAll(leavesInBox, i)
                {
                    const label leafI = leavesInBox[i];

                    const meshOctreeCube& oc = *leaves[leafI];

                    if (oc.intersectsLine(rootBox, sp, ep))
                    {
                        intersectedLeaves.append(leafI);

                        if (oc.level() < refLevels[emI])
                        {

                            refineCubes[leafI] = 1;
                            changed = true;
                        }
                    }
                }

                if (refThickness[emI] > VSMALL)
                {
                    useNLayers = true;

                    forAll(intersectedLeaves, i)
                    {
                        const label leafI = intersectedLeaves[i];

                        rThickness[leafI] =
                            Foam::max(rThickness[leafI], refThickness[emI]);
                    }
                }
            }
        }

        // mark additional boxes for refinement to achieve
        // correct refinement distance
        if (useNLayers && changed)
        {
            octreeModifier.refineSelectedBoxesAndAdditionalLayers
            (
                refineCubes,
                rThickness
            );
        }
        else if (changed)
        {
            // refine boxes
            octreeModifier.refineSelectedBoxes(refineCubes, hexRefinement_);
        }

    } while (changed);

    Info<< "Finished refinement of boxes intersecting edge meshes" << endl;
}


void Foam::Module::meshOctreeCreator::refineBoxesNearDataBoxes
(
    const label nLayers
)
{
    const FixedList<meshOctreeCubeCoordinates, 26>& rp =
        octree_.regularityPositions();

    Info<< "Refining boxes near DATA boxes" << endl;

    // ensure one to one matching with unknown boxes
    meshOctreeModifier octreeModifier(octree_);
    const DynamicList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    labelList refineBox(leaves.size(), 0);

    # pragma omp parallel for if (leaves.size() > 1000) \
    schedule(dynamic, 20)
    forAll(leaves, leafI)
    {
        if (leaves[leafI]->hasContainedElements())
        {
            const meshOctreeCube& oc = *leaves[leafI];

            for (label k = 0; k < 26; ++k)
            {
                const label neiLabel =
                    octree_.findLeafLabelForPosition
                    (
                        oc.coordinates() + rp[k]
                    );

                if (neiLabel == -1)
                {
                    continue;
                }

                const meshOctreeCube* nei = leaves[neiLabel];
                if (nei->level() == oc.level())
                {
                    continue;
                }
                if (nei->cubeType() & meshOctreeCubeBasic::OUTSIDE)
                {
                    continue;
                }

                refineBox[nei->cubeLabel()] = 1;
            }
        }
    }

    for (label i = 1; i < nLayers; i++)
    {
        # pragma omp parallel for if (leaves.size() > 1000) \
        schedule(dynamic, 20)
        forAll(leaves, leafI)
        {
            if (refineBox[leafI] == i)
            {
                const meshOctreeCube& oc = *leaves[leafI];

                for (label k = 0; k < 26; ++k)
                {
                    const label neiLabel =
                        octree_.findLeafLabelForPosition
                        (
                            oc.coordinates() + rp[k]
                        );

                    if (neiLabel == -1)
                    {
                        continue;
                    }

                    const meshOctreeCube* nei = leaves[neiLabel];
                    if (nei->level() == oc.level())
                    {
                        continue;
                    }
                    if (nei->cubeType() & meshOctreeCubeBasic::OUTSIDE)
                    {
                        continue;
                    }

                    refineBox[nei->cubeLabel()] = i + 1;
                }
            }
        }

    }

    // refine cubes
    octreeModifier.refineSelectedBoxes(refineBox, hexRefinement_);

    // set up inside-outside information
    createInsideOutsideInformation();

    Info<< "Finished refining boxes near DATA boxes" << endl;
}


void Foam::Module::meshOctreeCreator::refineBoxes
(
    const direction refLevel,
    const direction cubeType
)
{
    label nRefined;
    meshOctreeModifier octreeMod(octree_);

    do
    {
        nRefined = 0;

        const DynamicList<meshOctreeCube*>& leaves = octreeMod.leavesAccess();

        labelList refineCubes(leaves.size(), direction(0));

        # pragma omp parallel for if (leaves.size() > 1000) \
        reduction(+ : nRefined) schedule(dynamic, 20)
        forAll(leaves, leafI)
        {
            const meshOctreeCube& oc = *leaves[leafI];

            if ((oc.cubeType() & cubeType) && (oc.level() < refLevel))
            {
                refineCubes[leafI] = 1;
                ++nRefined;
            }
        }

        // refine selected cubes
        octreeMod.refineSelectedBoxes(refineCubes, hexRefinement_);

    } while (nRefined != 0);
}


// ************************************************************************* //
