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

#include "meshSurfaceEngineModifier.H"
#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"

#include "labelledPoint.H"
#include "helperFunctionsGeometryQueries.H"

using namespace Foam::Module::help;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::meshSurfaceEngineModifier::meshSurfaceEngineModifier
(
    meshSurfaceEngine& surfaceEngine
)
:
    surfaceEngine_(surfaceEngine)
{}


Foam::Module::meshSurfaceEngineModifier::meshSurfaceEngineModifier
(
    const meshSurfaceEngine& surfaceEngine
)
:
    surfaceEngine_(const_cast<meshSurfaceEngine&>(surfaceEngine))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::meshSurfaceEngineModifier::moveBoundaryVertexNoUpdate
(
    const label bpI,
    const point& newP
)
{
    surfaceEngine_.mesh_.points()[surfaceEngine_.boundaryPoints()[bpI]] = newP;
}


void Foam::Module::meshSurfaceEngineModifier::moveBoundaryVertex
(
    const label bpI,
    const point& newP
)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    pointFieldPMG& points = surfaceEngine_.mesh_.points();
    points[bPoints[bpI]] = newP;

    if (surfaceEngine_.faceCentresPtr_)
    {
        vectorField& faceCentres = *surfaceEngine_.faceCentresPtr_;
        const labelListList2D& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces[bpI][pfI];

            faceCentres[bfI] = stableFaceCentre(bFaces[bfI], points);
        }
    }

    if (surfaceEngine_.faceNormalsPtr_)
    {
        vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;
        const labelListList2D& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces[bpI][pfI];

            faceNormals[bfI] = stableFaceArea(bFaces[bfI], points);
        }
    }

    if (surfaceEngine_.pointNormalsPtr_)
    {
        const vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;
        const labelListList2D& pFaces = surfaceEngine_.pointFaces();
        const labelListList2D& pPoints = surfaceEngine_.pointPoints();

        vectorField& pn = *surfaceEngine_.pointNormalsPtr_;
        vector n(vector::zero);
        forAllRow(pFaces, bpI, pfI)
        {
            n += faceNormals[pFaces[bpI][pfI]];
        }

        const scalar l = mag(n);
        if (l > VSMALL)
        {
            n /= l;
        }
        else
        {
            n = vector::zero;
        }

        pn[bpI] = n;

        // change normal of vertices connected to bpI
        forAllRow(pPoints, bpI, ppI)
        {
            const label bpJ = pPoints[bpI][ppI];
            n = vector::zero;
            forAllRow(pFaces, bpJ, pfI)
            {
                n += faceNormals[pFaces[bpJ][pfI]];
            }

            const scalar d = mag(n);
            if (d > VSMALL)
            {
                n /= d;
            }
            else
            {
                n = vector::zero;
            }

            pn[bpJ] = n;
        }
    }
}


void Foam::Module::meshSurfaceEngineModifier::updateGeometry
(
    const DynamicList<label>& updateBndNodes
)
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const labelListList2D& pFaces = surfaceEngine_.pointFaces();
    const labelList& bp = surfaceEngine_.bp();

    boolList updateFaces(bFaces.size(), false);
    # pragma omp parallel for if (updateBndNodes.size() > 1000)
    forAll(updateBndNodes, i)
    {
        const label bpI = updateBndNodes[i];
        forAllRow(pFaces, bpI, j)
        {
            updateFaces[pFaces[bpI][j]] = true;
        }
    }

    if (surfaceEngine_.faceCentresPtr_)
    {
        vectorField& faceCentres = *surfaceEngine_.faceCentresPtr_;

        # pragma omp parallel for if (updateFaces.size() > 1000) \
        schedule(dynamic, 100)
        forAll(updateFaces, bfI)
        {
            if (updateFaces[bfI])
            {
                faceCentres[bfI] = stableFaceCentre(bFaces[bfI], points);
            }
        }
    }

    if (surfaceEngine_.faceNormalsPtr_)
    {
        vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;

        # pragma omp parallel for if (updateFaces.size() > 1000) \
        schedule(dynamic, 100)
        forAll(updateFaces, bfI)
        {
            if (updateFaces[bfI])
            {
                faceNormals[bfI] = stableFaceArea(bFaces[bfI], points);
            }
        }
    }

    if (surfaceEngine_.pointNormalsPtr_)
    {
        const vectorField& faceNormals = surfaceEngine_.faceNormals();

        boolList updateBndPoint(pFaces.size(), false);
        # pragma omp parallel for schedule(dynamic, 50)
        forAll(updateBndNodes, i)
        {
            const label bpI = updateBndNodes[i];

            forAllRow(pFaces, bpI, pfI)
            {
                const face& bf = bFaces[pFaces[bpI][pfI]];

                forAll(bf, pI)
                {
                    updateBndPoint[bp[bf[pI]]] = true;
                }
            }
        }

        vectorField& pn = *surfaceEngine_.pointNormalsPtr_;
        # pragma omp parallel for schedule(dynamic, 100)
        forAll(updateBndPoint, bpI)
        {
            if (!updateBndPoint[bpI])
            {
                continue;
            }

            vector n(vector::zero);
            forAllRow(pFaces, bpI, pfI)
            {
                n += faceNormals[pFaces[bpI][pfI]];
            }

            const scalar l = mag(n);
            if (l > VSMALL)
            {
                n /= l;
            }
            else
            {
                n = vector::zero;
            }

            pn[bpI] = n;
        }
    }
}


void Foam::Module::meshSurfaceEngineModifier::updateGeometry()
{
    DynamicList<label> updateBndNodes;
    updateBndNodes.setSize(surfaceEngine_.boundaryPoints().size());

    # pragma omp parallel for if (updateBndNodes.size() > 10000)
    forAll(updateBndNodes, bpI)
    {
        updateBndNodes[bpI] = bpI;
    }

    updateGeometry(updateBndNodes);
}


// ************************************************************************* //
