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

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::Module::labelListList2D&
Foam::Module::triSurfAddressing::pointFacets() const
{
    if (!pointFacetsPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate pointFacets" << abort(FatalError);
        }

        calculatePointFacets();
    }

    return *pointFacetsPtr_;
}


inline const Foam::DynamicList<Foam::edge>&
Foam::Module::triSurfAddressing::edges() const
{
    if (!edgesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate edges" << abort(FatalError);
        }

        calculateEdges();
    }

    return *edgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::triSurfAddressing::facetEdges() const
{
    if (!facetEdgesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate facetEdges" << abort(FatalError);
        }

        calculateFacetEdges();
    }

    return *facetEdgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::triSurfAddressing::edgeFacets() const
{
    if (!edgeFacetsPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate edgeFacets" << abort(FatalError);
        }

        calculateEdgeFacets();
    }

    return *edgeFacetsPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::triSurfAddressing::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate pointEdges" << abort(FatalError);
        }

        calculatePointEdges();
    }

    return *pointEdgesPtr_;
}


inline const Foam::vectorField&
Foam::Module::triSurfAddressing::pointNormals() const
{
    if (!pointNormalsPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate pointNormals" << abort(FatalError);
        }

        calculatePointNormals();
    }

    return *pointNormalsPtr_;
}


inline const Foam::vectorField&
Foam::Module::triSurfAddressing::facetNormals() const
{
    if (!facetNormalsPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate facetNormals" << abort(FatalError);
        }

        calculateFacetNormals();
    }

    return *facetNormalsPtr_;
}


inline const Foam::vectorField&
Foam::Module::triSurfAddressing::facetCentres() const
{
    if (!facetCentresPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Cannot calculate facetCentres" << abort(FatalError);
        }

        calculateFacetCentres();
    }

    return *facetCentresPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
