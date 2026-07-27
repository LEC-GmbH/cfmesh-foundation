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

#include "polyMeshGenAddressing.H"

#include <omp.h>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Module::polyMeshGenAddressing::calcPointPoints() const
{
    if (ppPtr_)
    {
        FatalErrorInFunction
            << "pointPoints already calculated"
            << abort(FatalError);
    }
    else
    {
        ppPtr_ = new labelListList2D();
        labelListList2D& pp = *ppPtr_;

        const faceListPMG& faces = mesh_.faces();
        const labelListList2D& pointFaces = this->pointFaces();

        labelList nPoints(pointFaces.size());

        const label nThreads = 3*omp_get_num_procs();
        # pragma omp parallel num_threads(nThreads) if (nPoints.size() > 10000)
        {
            # pragma omp for schedule(static)
            forAll(nPoints, i)
                nPoints[i] = 0;

            # pragma omp for schedule(static)
            forAll(pointFaces, pointI)
            {
                DynamicList<label> helper; helper.reserve(32);

                forAllRow(pointFaces, pointI, pfI)
                {
                    const face& f = faces[pointFaces[pointI][pfI]];

                    const label pos = f.which(pointI);
                    appendUniq(helper, f.prevLabel(pos));
                    appendUniq(helper, f.nextLabel(pos));
                }

                nPoints[pointI] = helper.size();
            }

            # pragma omp barrier

            # pragma omp master
            setSizeAndRowSize(pp, nPoints);

            # pragma omp barrier

            # pragma omp for schedule(static)
            forAll(pointFaces, pointI)
            {
                DynamicList<label> helper; helper.reserve(32);

                forAllRow(pointFaces, pointI, pfI)
                {
                    const face& f = faces[pointFaces[pointI][pfI]];

                    const label pos = f.which(pointI);
                    const label pLabel = f.prevLabel(pos);
                    const label nLabel = f.nextLabel(pos);

                    appendUniq(helper, nLabel);
                    appendUniq(helper, pLabel);
                }

                setRow(pp, pointI, helper);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Module::labelListList2D&
Foam::Module::polyMeshGenAddressing::pointPoints() const
{
    if (!ppPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calcPointPoints();
    }

    return *ppPtr_;
}


// ************************************************************************* //
