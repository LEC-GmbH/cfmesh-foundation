/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | www.cfmesh.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 LEC GmbH
-------------------------------------------------------------------------------
Author
    Clemens Goessnitzer <Clemens.Goessnitzer@lec.tugraz.at>

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

#include "surfaceRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "orientedSurface.H"

#include <omp.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
    defineTypeNameAndDebug(surfaceRefinement, 0);
    addToRunTimeSelectionTable
    (
        objectRefinement,
        surfaceRefinement,
        dictionary
    );
}
}


// * * * * * * * * * * * *  Private Member Functions   * * * * * * * * * * * //

Foam::label Foam::Module::surfaceRefinement::maxThreads() const
{
    return omp_get_max_threads();
}


Foam::label Foam::Module::surfaceRefinement::threadIdx() const
{
    return omp_get_thread_num();
}


void Foam::Module::surfaceRefinement::initialise()
{
    forAll(search_, procI)
    {
        surface_[procI].reset(new triSurface(fName_));
        orientedSurface::orient
        (
            surface_[procI](),
            boundBox(surface_[procI]().points()).min() - vector(1, 1, 1),
            true
        );
        search_[procI].reset(new triSurfaceSearch(surface_[procI]()));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::surfaceRefinement::surfaceRefinement()
:
    objectRefinement(),
    distance_(-1),
    surface_(maxThreads()),
    search_(maxThreads())
{}


Foam::Module::surfaceRefinement::surfaceRefinement
(
    const word& name,
    const scalar cellSize,
    const direction additionalRefLevels,
    const fileName fName
)
:
    objectRefinement(),
    fName_(fName),
    distance_(-1),
    surface_(maxThreads()),
    search_(maxThreads())
{
    setName(name);
    setCellSize(cellSize);
    setAdditionalRefinementLevels(additionalRefLevels);
    initialise();
}


Foam::Module::surfaceRefinement::surfaceRefinement
(
    const word& name,
    const dictionary& dict
)
:
    objectRefinement(name, dict),
    distance_(-1),
    surface_(maxThreads()),
    search_(maxThreads())
{
    this->operator=(dict);

    initialise();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Module::surfaceRefinement::intersectsObject(const boundBox& bb) const
{
    boolList isInside(bb.points()().size());

    if (distance_ > 0)
    {
        const vector span(distance_, distance_, distance_);
        pointIndexHit inter = search_[threadIdx()]().nearest
        (
            bb.midpoint(),
            span
        );

        if (inter.hit() && mag(inter.hitPoint() - bb.midpoint()) < distance_)
        {
            return true;
        }
    }

    isInside = search_[threadIdx()]().calcInside(bb.points());
    forAll(isInside, i)
    {
        if (isInside[i])
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary
Foam::Module::surfaceRefinement::dict(bool /*ignoreType*/) const
{
    dictionary dict;

    if (additionalRefinementLevels() == 0 && cellSize() >= 0.0)
    {
        dict.add("cellSize", cellSize());
    }
    else
    {
        dict.add("additionalRefinementLevels", additionalRefinementLevels());
    }

    dict.add("type", type());

    dict.add("file", fName_);

    if (distance_ > 0)
    {
        dict.add("distance", distance_);
    }

    return dict;
}


void Foam::Module::surfaceRefinement::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " file:   " << fName_;
}


void Foam::Module::surfaceRefinement::writeDict
(
    Ostream& os,
    bool subDict
) const
{
    if (subDict)
    {
        os << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    if (additionalRefinementLevels() == 0 && cellSize() >= 0.0)
    {
        os.writeKeyword("cellSize") << cellSize() << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeKeyword("additionalRefinementLevels")
            << additionalRefinementLevels()
            << token::END_STATEMENT << nl;
    }

    // only write type for derived types
    if (type() != typeName_())
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    }

    os.writeKeyword("file") << fName_ << token::END_STATEMENT << nl;

    if (distance_ > 0)
    {
        os.writeKeyword("distance") << distance_ << token::END_STATEMENT << nl;
    }

    if (subDict)
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}


void Foam::Module::surfaceRefinement::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        d.found(typeName_())
      ? d.subDict(typeName_())
      : d
    );

    fName_ = dict.lookup<fileName>("file");
    distance_ = dict.lookupOrDefault<scalar>("distance", -1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::Module::surfaceRefinement::operator<<
(
    Ostream& os
) const
{
    os << "name " << name() << nl;
    os << "cellSize " << cellSize() << nl;
    os << "additionalRefinementLevels " << additionalRefinementLevels() << endl;

    write(os);
    return os;
}


// ************************************************************************* //
