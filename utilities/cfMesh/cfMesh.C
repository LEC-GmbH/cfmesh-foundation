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

Application
    cfMesh

Description
    Takes a triangulated surface and generates a mesh, using the mesh
    generation mode selected by -mode or meshDict's meshMode entry
    (cartesian, cartesian2D, tetrahedral, polyhedral).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "systemDict.H"
#include "polyMesh.H"
#include "meshGenerator.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Write the mesh as it stands after a stage, matching snappyHexMesh's own
// writeMesh()-after-every-enabled-phase behaviour: by default (overwrite)
// every stage lands in the same constant/ location; with -noOverwrite, each
// stage advances runTime and lands in its own time directory instead,
// mirroring snappyHexMesh.C's mesh.time()++. Always renumbers first - every
// written mesh is a complete, real mesh someone may inspect or use
// directly (e.g. an intermediate -noOverwrite stage), not just the final
// one, so every write should be well-formed, not just the last.
void writeStageMesh
(
    const string& stageName,
    Module::meshGenerator& generator,
    const word& regionName,
    const fileName& meshPath,
    const bool overwrite,
    Time& runTime
)
{
    generator.renumberMesh();

    if (!overwrite)
    {
        runTime++;
    }

    IOobject meshIO
    (
        regionName,
        overwrite ? runTime.constant() : runTime.name(),
        meshPath,
        runTime
    );

    Info<< "Writing " << stageName << " mesh to " << meshIO.instance()
        << endl;

    autoPtr<polyMesh> meshPtr(generator.generatePolyMesh(meshIO));
    meshPtr->write();
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "(cfmesh)\n"
        "Takes a triangulated surface and generates a mesh, using the mesh"
        " generation mode selected by -mode or meshDict's meshMode entry"
    );

    argList::addOption
    (
        "mode",
        "mode",
        "mesh generation mode (cartesian | cartesian2D | tetrahedral |"
        " polyhedral); overrides meshDict's meshMode entry if both are given"
    );

    #include "addDictOption.H"
    #include "addMeshOption.H"
    #include "addRegionOption.H"
    #include "addNoOverwriteOption.H"

    // None of cfMesh's mesh generation modes support parallel execution
    argList::noParallel();

    #include "setRootCaseNoFunctionObjects.H"
    #include "setMeshPath.H"
    #include "createTimeNoFunctionObjects.H"
    #include "setRegionName.H"
    #include "setNoOverwrite.H"

    typeIOobject<IOdictionary> meshDictIO
    (
        systemDictIO("meshDict", args, runTime, regionName, meshPath)
    );

    if (!meshDictIO.headerOk())
    {
        FatalErrorInFunction
            << "Cannot find file " << meshDictIO.relativeObjectPath()
            << nl
            << exit(FatalError);
    }

    IOdictionary meshDict(meshDictIO);

    word modeName;
    if (!args.optionReadIfPresent("mode", modeName))
    {
        meshDict.readIfPresent("meshMode", modeName);
    }

    autoPtr<Module::meshGenerator> generatorPtr
    (
        Module::meshGenerator::New
        (
            modeName,
            runTime,
            meshDict,
            regionName,
            meshPath
        )
    );

    // Stage control mirroring snappyHexMeshDict's castellatedMesh/snap/
    // addLayers switches exactly, including the number of stages (3, not
    // 4 - there is no separate "finalise" stage; the implicit collapse of
    // feature-edge-split patches back to their pre-split names happens
    // inside createCastellatedMesh() itself, and renumbering happens on
    // every write, inside writeStageMesh() itself, same as snappyHexMesh
    // has no separate post-addLayers stage either) - no cfMesh-side
    // bookkeeping of what already ran; a disabled castellatedMesh just
    // means reading back whatever mesh currently exists instead of
    // generating one.
    if (meshDict.lookupOrDefault<bool>("castellatedMesh", true))
    {
        generatorPtr->createCastellatedMesh();
        writeStageMesh
        (
            "castellated",
            generatorPtr(),
            regionName,
            meshPath,
            overwrite,
            runTime
        );
    }
    else
    {
        generatorPtr->readMesh();
    }

    if (meshDict.lookupOrDefault<bool>("snap", true))
    {
        generatorPtr->snap();
        writeStageMesh
        (
            "snapped",
            generatorPtr(),
            regionName,
            meshPath,
            overwrite,
            runTime
        );
    }

    if (meshDict.lookupOrDefault<bool>("addLayers", true))
    {
        generatorPtr->addLayers();
        writeStageMesh
        (
            "layer",
            generatorPtr(),
            regionName,
            meshPath,
            overwrite,
            runTime
        );
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n"
        << "ClockTime = " << runTime.elapsedClockTime() << " s" << endl;

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
