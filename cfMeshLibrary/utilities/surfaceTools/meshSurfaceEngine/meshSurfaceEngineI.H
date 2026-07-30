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

#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::Module::polyMeshGen&
Foam::Module::meshSurfaceEngine::mesh() const
{
    return mesh_;
}


inline const Foam::Module::pointFieldPMG&
Foam::Module::meshSurfaceEngine::points() const
{
    return mesh_.points();
}


inline const Foam::Module::faceListPMG&
Foam::Module::meshSurfaceEngine::faces() const
{
    return mesh_.faces();
}


inline const Foam::Module::cellListPMG&
Foam::Module::meshSurfaceEngine::cells() const
{
    return mesh_.cells();
}


inline const Foam::labelList&
Foam::Module::meshSurfaceEngine::bp() const
{
    if (!bppPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateBoundaryFaces();
        calculateBoundaryNodes();
    }

    return *bppPtr_;
}


inline const Foam::labelList&
Foam::Module::meshSurfaceEngine::boundaryPoints() const
{
    if (!boundaryPointsPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateBoundaryNodes();
    }

    return *boundaryPointsPtr_;
}


inline const Foam::faceList::subList&
Foam::Module::meshSurfaceEngine::boundaryFaces() const
{
    if (!boundaryFacesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateBoundaryFaces();
    }

    return *boundaryFacesPtr_;
}


inline const Foam::labelList&
Foam::Module::meshSurfaceEngine::boundaryFacePatches() const
{
    if (!boundaryFacePatchPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateBoundaryFacePatches();
    }

    return *boundaryFacePatchPtr_;
}


inline const Foam::labelList&
Foam::Module::meshSurfaceEngine::faceOwners() const
{
    if (!boundaryFaceOwnersPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateBoundaryOwners();
    }

    return *boundaryFaceOwnersPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshSurfaceEngine::pointFaces() const
{
    if (!pointFacesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculatePointFaces();
    }

    return *pointFacesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshSurfaceEngine::pointInFaces() const
{
    if (!pointInFacePtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculatePointFaces();
    }

    return *pointInFacePtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshSurfaceEngine::pointPoints() const
{
    if (!pointPointsPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculatePointPoints();
    }

    return *pointPointsPtr_;
}


inline const Foam::vectorField&
Foam::Module::meshSurfaceEngine::pointNormals() const
{
    if (!pointNormalsPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculatePointNormals();
    }

    return *pointNormalsPtr_;
}


inline const Foam::vectorField&
Foam::Module::meshSurfaceEngine::faceNormals() const
{
    if (!faceNormalsPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateFaceNormals();
    }

    return *faceNormalsPtr_;
}


inline const Foam::vectorField&
Foam::Module::meshSurfaceEngine::faceCentres() const
{
    if (!faceCentresPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateFaceCentres();
    }

    return *faceCentresPtr_;
}


inline const Foam::edgeList&
Foam::Module::meshSurfaceEngine::edges() const
{
    if (!edgesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateEdgesAndAddressing();
    }

    return *edgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshSurfaceEngine::boundaryPointEdges() const
{
    if (!bpEdgesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateEdgesAndAddressing();
    }

    return *bpEdgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshSurfaceEngine::edgeFaces() const
{
    if (!edgeFacesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateEdgeFacesAddressing();
    }

    return *edgeFacesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshSurfaceEngine::faceEdges() const
{
    if (!faceEdgesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateFaceEdgesAddressing();
    }

    return *faceEdgesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshSurfaceEngine::edgePatches() const
{
    if (!edgePatchesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateEdgePatchesAddressing();
    }

    return *edgePatchesPtr_;
}


inline const Foam::Module::labelListList2D&
Foam::Module::meshSurfaceEngine::faceFaces() const
{
    if (!faceFacesPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateFaceFacesAddressing();
    }

    return *faceFacesPtr_;
}

// ************************************************************************* //
