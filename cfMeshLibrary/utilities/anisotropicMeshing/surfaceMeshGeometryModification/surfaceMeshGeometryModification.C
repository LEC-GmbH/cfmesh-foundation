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

#include "surfaceMeshGeometryModification.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "triSurf.H"

// * * * * * * * * * * * * * * Private member functions* * * * * * * * * * * //

void Foam::Module::surfaceMeshGeometryModification::checkModification
(
    const dictionary& anisotropicSourcesDict
)
{
    // active whenever the caller passed a real dict at all (i.e. the
    // anisotropicSources key was present in meshDict), matching this
    // class's historical behaviour - even an empty anisotropicSources {}
    // block must activate (constructing a no-op coordinateModifier),
    // NOT just a non-empty one: modifyGeometry() unconditionally returns
    // a real triSurf* only when active, and every caller dereferences
    // that pointer without a null check
    if (&anisotropicSourcesDict != &dictionary::null)
    {
        modificationActive_ = true;

        coordinateModifierPtr_ = new coordinateModifier(anisotropicSourcesDict);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::surfaceMeshGeometryModification::surfaceMeshGeometryModification
(
    const triSurf& surf,
    const dictionary& anisotropicSourcesDict
)
:
    surf_(surf),
    coordinateModifierPtr_(nullptr),
    modificationActive_(false)
{
    checkModification(anisotropicSourcesDict);
}


Foam::Module::surfaceMeshGeometryModification::
~surfaceMeshGeometryModification()
{
    deleteDemandDrivenData(coordinateModifierPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Module::surfaceMeshGeometryModification::activeModification() const
{
    return modificationActive_;
}


const Foam::Module::triSurf*
Foam::Module::surfaceMeshGeometryModification::modifyGeometry() const
{
    if (!modificationActive_)
    {
        WarningInFunction
            << "Modification is not active" << endl;

        return nullptr;
    }

    const pointField& pts = surf_.points();

    pointField newPts(pts.size());

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(pts, pointI)
    {
        newPts[pointI] = coordinateModifierPtr_->modifiedPoint(pts[pointI]);
    }

    triSurf* newSurf =
        new triSurf
        (
            surf_.facets(),
            surf_.patches(),
            surf_.featureEdges(),
            newPts
        );

    // copy subsets
    DynamicList<label> sIds;

    // copy facet subsets
    surf_.facetSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addFacetSubset(surf_.facetSubsetName(sIds[i]));

        DynamicList<label> facetsInSubset;
        surf_.facetsInSubset(sIds[i], facetsInSubset);

        forAll(facetsInSubset, fI)
        {
            newSurf->addFacetToSubset(newId, facetsInSubset[fI]);
        }
    }

    // copy point subsets
    surf_.pointSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addPointSubset(surf_.pointSubsetName(sIds[i]));

        DynamicList<label> pointsInSubset;
        surf_.pointsInSubset(sIds[i], pointsInSubset);

        forAll(pointsInSubset, pI)
        {
            newSurf->addPointToSubset(newId, pointsInSubset[pI]);
        }
    }

    // copy subsets of feature edges
    surf_.edgeSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addEdgeSubset(surf_.edgeSubsetName(sIds[i]));

        DynamicList<label> edgesInSubset;
        surf_.edgesInSubset(sIds[i], edgesInSubset);

        forAll(edgesInSubset, eI)
        {
            newSurf->addEdgeToSubset(newId, edgesInSubset[eI]);
        }
    }

    // carry the pre-split -> post-split patch name/type map across, so
    // downstream patch-name resolution (findPatches()) keeps working on
    // the coordinate-modified copy exactly as it did on surf_ - without
    // this, any surface using anisotropicSources silently loses the map
    // the moment it's coordinate-modified
    newSurf->setOrigPatchNamesMap(surf_.origPatchNamesMap());
    newSurf->setOrigPatchTypesMap(surf_.origPatchTypesMap());

    return newSurf;
}


const Foam::Module::triSurf*
Foam::Module::surfaceMeshGeometryModification::revertGeometryModification()
const
{
    if (!modificationActive_)
    {
        WarningInFunction
            << "Modification is not active" << endl;

        return nullptr;
    }

    const pointField& pts = surf_.points();

    pointField newPts(pts.size());

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(pts, pointI)
    {
        newPts[pointI] =
            coordinateModifierPtr_->backwardModifiedPoint(pts[pointI]);
    }

    triSurf* newSurf =
        new triSurf
        (
            surf_.facets(),
            surf_.patches(),
            surf_.featureEdges(),
            newPts
        );

    // copy subsets
    DynamicList<label> sIds;

    // copy facet subsets
    surf_.facetSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addFacetSubset(surf_.facetSubsetName(sIds[i]));

        DynamicList<label> facetsInSubset;
        surf_.facetsInSubset(sIds[i], facetsInSubset);

        forAll(facetsInSubset, fI)
        {
            newSurf->addFacetToSubset(newId, facetsInSubset[fI]);
        }
    }

    // copy point subsets
    surf_.pointSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addPointSubset(surf_.pointSubsetName(sIds[i]));

        DynamicList<label> pointsInSubset;
        surf_.pointsInSubset(sIds[i], pointsInSubset);

        forAll(pointsInSubset, pI)
        {
            newSurf->addPointToSubset(newId, pointsInSubset[pI]);
        }
    }

    // copy subsets of feature edges
    surf_.edgeSubsetIndices(sIds);
    forAll(sIds, i)
    {
        const label newId =
            newSurf->addEdgeSubset(surf_.edgeSubsetName(sIds[i]));

        DynamicList<label> edgesInSubset;
        surf_.edgesInSubset(sIds[i], edgesInSubset);

        forAll(edgesInSubset, eI)
        {
            newSurf->addEdgeToSubset(newId, edgesInSubset[eI]);
        }
    }

    // carry the pre-split -> post-split patch name/type map across, so
    // downstream patch-name resolution (findPatches()) keeps working on
    // the reverted copy exactly as it did on surf_
    newSurf->setOrigPatchNamesMap(surf_.origPatchNamesMap());
    newSurf->setOrigPatchTypesMap(surf_.origPatchTypesMap());

    return newSurf;
}


// ************************************************************************* //
