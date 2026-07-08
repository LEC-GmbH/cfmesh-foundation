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
#include "DynList.H"
#include "VRWGraphSMPModifier.H"

#include <omp.h>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Module::polyMeshGenAddressing::calcCellEdges() const
{
    if (cePtr_)
    {
        FatalErrorInFunction
            << "cellEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        const cellListPMG& cells = mesh_.cells();
        const VRWGraph& fe = faceEdges();

        cePtr_ = new VRWGraph();
        VRWGraph& cellEdgeAddr = *cePtr_;

        labelList nEdges(cells.size());

        const label nThreads = 3*omp_get_num_procs();

        # pragma omp parallel num_threads(nThreads) if (cells.size() > 10000)
        {
            # pragma omp for schedule(static)
            forAll(nEdges, i)
            {
                nEdges[i] = 0;
            }

            # pragma omp for schedule(static)
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];

                DynList<label, 32> cEdges;

                forAll(c, fI)
                {
                    const label faceI = c[fI];

                    forAllRow(fe, faceI, eI)
                    {
                        cEdges.appendUniq(fe(faceI, eI));
                    }
                }

                nEdges[cellI] = cEdges.size();
            }

            # pragma omp barrier

            # pragma omp master
            VRWGraphSMPModifier(cellEdgeAddr).setSizeAndRowSize(nEdges);

            # pragma omp barrier

            # pragma omp for schedule(static)
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];

                DynList<label, 32> cEdges;

                forAll(c, fI)
                {
                    const label faceI = c[fI];

                    forAllRow(fe, faceI, eI)
                    {
                        cEdges.appendUniq(fe(faceI, eI));
                    }
                }

                cellEdgeAddr.setRow(cellI, cEdges);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Module::VRWGraph&
Foam::Module::polyMeshGenAddressing::cellEdges() const
{
    if (!cePtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calcCellEdges();
    }

    return *cePtr_;
}


// ************************************************************************* //
