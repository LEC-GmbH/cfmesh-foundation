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

Description
    Translates FOAM mesh to AVL's FPMA format

\*---------------------------------------------------------------------------*/

#include "polyMeshGen.H"
#include "meshSurfaceEngine.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "fileName.H"

#include "fpmaMesh.H"
#include "writeMeshFPMA.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::writeMeshFPMA(const polyMeshGen& mesh, const word& fName)
{
    const Time& time = mesh.returnTime();

    const word postProcDir = "FPMA";

    fileName postProcPath = time.path()/postProcDir;

    if (!Foam::isDir(postProcPath))
    {
        mkDir(postProcPath);
    }

    // Open the Case file
    const fileName fpmaFileName = fName + ".fpma";

    Info<< "Writing mesh into " << fpmaFileName << endl;

/*    OFstream fpmaGeometryFile
    (
        postProcPath/fpmaFileName,
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );
*/

    OFstream fpmaGeometryFile(postProcPath/fpmaFileName);

    // Construct the FIRE mesh
    fpmaMesh Mesh(mesh);
    Mesh.write(fpmaGeometryFile);
}


// ************************************************************************* //
