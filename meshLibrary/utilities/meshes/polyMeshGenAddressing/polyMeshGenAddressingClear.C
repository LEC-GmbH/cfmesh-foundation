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
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Module::polyMeshGenAddressing::clearGeom()
{
    if (debug)
    {
        Pout<< "polyMeshGenAddressing::clearGeom() : "
            << "clearing geometric data"
            << endl;
    }

    deleteDemandDrivenData(cellCentresPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(cellVolumesPtr_);
    deleteDemandDrivenData(faceAreasPtr_);
}


void Foam::Module::polyMeshGenAddressing::clearAddressing()
{
    if (debug)
    {
        Pout<< "polyMeshGenAddressing::clearAddressing() : "
            << "clearing topology"
            << endl;
    }

    clearOutEdges();

    deleteDemandDrivenData(ccPtr_);
    deleteDemandDrivenData(pcPtr_);

    deleteDemandDrivenData(efPtr_);
    deleteDemandDrivenData(pfPtr_);

    deleteDemandDrivenData(cePtr_);
    deleteDemandDrivenData(fePtr_);
    deleteDemandDrivenData(pePtr_);
    deleteDemandDrivenData(ppPtr_);
    deleteDemandDrivenData(cpPtr_);
}


void Foam::Module::polyMeshGenAddressing::clearOut()
{
    clearGeom();
    clearAddressing();
}


void Foam::Module::polyMeshGenAddressing::clearAll()
{
    clearGeom();
    clearAddressing();
}


// ************************************************************************* //
