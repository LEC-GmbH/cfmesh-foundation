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

#include "edgeMeshGeometryModification.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "edgeMesh.H"

// * * * * * * * * * * * * * * Private member functions* * * * * * * * * * * //

void Foam::Module::edgeMeshGeometryModification::checkModification
(
    const dictionary& anisotropicSourcesDict
)
{
    // active whenever the caller passed a real dict at all (i.e. the
    // anisotropicSources key was present in meshDict), matching this
    // class's historical behaviour - even an empty anisotropicSources {}
    // block must activate (constructing a no-op coordinateModifier),
    // NOT just a non-empty one: modifyGeometry() unconditionally returns
    // a real edgeMesh* only when active, and every caller dereferences
    // that pointer without a null check
    if (&anisotropicSourcesDict != &dictionary::null)
    {
        modificationActive_ = true;

        coordinateModifierPtr_ = new coordinateModifier(anisotropicSourcesDict);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::edgeMeshGeometryModification::edgeMeshGeometryModification
(
    const edgeMesh& em,
    const dictionary& anisotropicSourcesDict
)
:
    edgeMesh_(em),
    coordinateModifierPtr_(nullptr),
    modificationActive_(false)
{
    checkModification(anisotropicSourcesDict);
}


Foam::Module::edgeMeshGeometryModification::~edgeMeshGeometryModification()
{
    deleteDemandDrivenData(coordinateModifierPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Module::edgeMeshGeometryModification::activeModification() const
{
    return modificationActive_;
}


const Foam::edgeMesh*
Foam::Module::edgeMeshGeometryModification::modifyGeometry() const
{
    if (!modificationActive_)
    {
        WarningInFunction
            << "Modification is not active" << endl;

        return nullptr;
    }

    const pointField& pts = edgeMesh_.points();

    pointField newPts(pts.size());

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(pts, pointI)
    {
        newPts[pointI] = coordinateModifierPtr_->modifiedPoint(pts[pointI]);
    }

    const edgeMesh* newEdgeMesh = new edgeMesh(newPts, edgeMesh_.edges());

    return newEdgeMesh;
}


const Foam::edgeMesh*
Foam::Module::edgeMeshGeometryModification::revertGeometryModification() const
{
    if (!modificationActive_)
    {
        WarningInFunction
            << "Modification is not active" << endl;

        return nullptr;
    }

    const pointField& pts = edgeMesh_.points();

    pointField newPts(pts.size());

    # pragma omp parallel for schedule(dynamic, 50)
    forAll(pts, pointI)
    {
        newPts[pointI] =
            coordinateModifierPtr_->backwardModifiedPoint(pts[pointI]);
    }

    const edgeMesh* newEdgeMeshPtr = new edgeMesh(newPts, edgeMesh_.edges());

    return newEdgeMeshPtr;
}


// ************************************************************************* //
