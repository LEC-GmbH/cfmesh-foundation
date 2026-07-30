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
// mirroring snappyHexMesh.C's mesh.time()++.
void writeStageMesh
(
    const string& stageName,
    const Module::meshGenerator& generator,
    const word& regionName,
    const fileName& meshPath,
    const bool overwrite,
    Time& runTime
)
{
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
    // addLayers switches exactly - no cfMesh-side bookkeeping of what
    // already ran; a disabled castellatedMesh just means reading back
    // whatever mesh currently exists instead of generating one.
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

    const bool snapEnabled = meshDict.lookupOrDefault<bool>("snap", true);
    if (snapEnabled)
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
    else if
    (
        snapEnabled
     && meshDict.lookupOrDefault<bool>("optimiseMesh", true)
    )
    {
        // addLayers() also runs final mesh-quality optimisation as part
        // of its own sequence - run it here too, so disabling addLayers
        // doesn't silently skip mesh cleanup along with the layers.
        // Only meaningful once the surface has actually been snapped -
        // if snap was disabled too, there is nothing yet worth
        // optimising against (e.g. a deliberate castellatedMesh-only
        // stop, matching workflowControls.stopAfter's old
        // templateGeneration checkpoint).
        // optimiseMesh defaults to true (matching addLayers/snap/
        // castellatedMesh's own always-on-by-default convention), but
        // can be explicitly disabled for a staged workflow that intends
        // to run addLayers in a later invocation - addLayers() runs this
        // same optimisation itself once layers actually exist, so
        // running it prematurely here would perturb the mesh addLayers
        // is about to work with, changing the eventual result
        generatorPtr->optimiseMesh(false);
        writeStageMesh
        (
            "optimised",
            generatorPtr(),
            regionName,
            meshPath,
            overwrite,
            runTime
        );
    }

    if (meshDict.lookupOrDefault<bool>("finalise", true))
    {
        // renumberMesh() is harmless regardless of when it runs, but
        // replaceBoundaries() (patch renaming, per meshDict's
        // renameBoundary) is not: refineBoundaryLayers/boundaryLayers'
        // own patchBoundaryLayers matching (called from addLayers(),
        // possibly in a *later* invocation of a staged workflow) matches
        // against whatever patch names currently exist on the mesh - if
        // renaming already happened here, those patch-name regexes see
        // the post-rename names instead of the original ones they were
        // written against, and silently match fewer (or no) patches.
        // Defer finalise (set finalise false) in any invocation that
        // isn't the one actually running addLayers, in a staged
        // workflow that splits snap and addLayers into separate
        // invocations
        generatorPtr->finalise();
        writeStageMesh
        (
            "final",
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
