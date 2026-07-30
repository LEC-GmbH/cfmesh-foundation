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

#include "polyMeshGenGeometryModification.H"
#include "dictionary.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::Module::polyMeshGenGeometryModification::checkModification
(
    const dictionary& anisotropicSourcesDict
)
{
    // active whenever the caller passed a real dict at all (i.e. the
    // anisotropicSources key was present in meshDict), matching this
    // class's historical behaviour - even an empty anisotropicSources {}
    // block must activate (constructing a no-op coordinateModifier), NOT
    // just a non-empty one: modifyGeometry()/revertGeometryModification()
    // are only ever called by callers who already checked
    // activeModification() themselves for this class, but the sibling
    // surface-/edge-mesh variants are dereferenced unconditionally by
    // their callers, so keep the semantics identical across all three
    if (&anisotropicSourcesDict != &dictionary::null)
    {
        modificationActive_ = true;

        coordinateModifierPtr_ = new coordinateModifier(anisotropicSourcesDict);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::polyMeshGenGeometryModification::polyMeshGenGeometryModification
(
    polyMeshGen& mesh,
    const dictionary& anisotropicSourcesDict
)
:
    mesh_(mesh),
    coordinateModifierPtr_(nullptr),
    modificationActive_(false)
{
    checkModification(anisotropicSourcesDict);
}


Foam::Module::polyMeshGenGeometryModification::
~polyMeshGenGeometryModification()
{
    deleteDemandDrivenData(coordinateModifierPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Module::polyMeshGenGeometryModification::activeModification() const
{
    return modificationActive_;
}


void Foam::Module::polyMeshGenGeometryModification::modifyGeometry()
{
    if (!modificationActive_)
    {
        WarningInFunction
            << "Modification is not active" << endl;

        return;
    }

    pointFieldPMG& pts = mesh_.points();

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(pts, pointI)
    {
        pts[pointI] = coordinateModifierPtr_->modifiedPoint(pts[pointI]);
    }
}


void Foam::Module::polyMeshGenGeometryModification::revertGeometryModification()
{
    if (!modificationActive_)
    {
        WarningInFunction
            << "Modification is not active" << endl;

        return;
    }

    pointFieldPMG& pts = mesh_.points();

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(pts, pointI)
    {
        pts[pointI] =
            coordinateModifierPtr_->backwardModifiedPoint(pts[pointI]);
    }
}


// ************************************************************************* //
