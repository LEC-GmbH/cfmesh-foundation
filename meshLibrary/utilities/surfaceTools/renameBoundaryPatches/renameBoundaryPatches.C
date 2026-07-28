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

#include "renameBoundaryPatches.H"
#include "symmetryPlaneOptimisation.H"
#include "demandDrivenData.H"
#include "IOdictionary.H"
#include "DynamicList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::renameBoundaryPatches::calculateNewBoundary()
{
    Info<< "Renaming boundary patches" << endl;

    const dictionary& dict = meshDict_.subDict("renameBoundary");

    HashTable<label> patchToLabel;
    forAll(mesh_.boundaries(), patchI)
    {
        patchToLabel.insert
        (
            mesh_.boundaries()[patchI].patchName(),
            patchI
        );
    }

    labelList patchToNew(mesh_.boundaries().size(), -1);

    // Map every boundary face to the patch it currently belongs to.
    // Note: named patches are not guaranteed to tile
    // [nIntFaces, faces.size()) exactly at every point in the mesh
    // pipeline (e.g. after mesh-optimisation faces/cells have been added
    // without the patch list being kept in sync) - a face left outside
    // all patches (-1) is routed to the default patch below, the same
    // way a face from an administratively-removed patch is handled.
    const label nIntFaces = mesh_.nInternalFaces();
    const faceListPMG& meshFaces = mesh_.faces();

    labelList facePatch(meshFaces.size() - nIntFaces, -1);
    bool hasOrphanFaces(false);
    forAll(mesh_.boundaries(), patchI)
    {
        const boundaryPatch& wp = mesh_.boundaries()[patchI];
        const label start = wp.patchStart();
        const label end = start + wp.patchSize();

        for (label faceI = start; faceI < end; ++faceI)
        {
            facePatch[faceI - nIntFaces] = patchI;
        }
    }
    forAll(facePatch, i)
    {
        if (facePatch[i] == -1)
        {
            hasOrphanFaces = true;
            break;
        }
    }

    // +1: room for a default patch beyond one-slot-per-original-patch,
    // needed when orphaned faces force a default patch into existence
    // even though every original patch already got its own new slot
    wordList newPatchNames(patchToNew.size() + 1);
    wordList newPatchTypes(patchToNew.size() + 1);
    HashTable<label> newNameToPos;
    label newPatchI(0);

    // read new patch names and types
    if (dict.found("newPatchNames"))
    {
        PtrList<entry> patchesToRename;

        if (dict.isDict("newPatchNames"))
        {
            const dictionary& newPatchNames = dict.subDict("newPatchNames");
            const wordList keys = newPatchNames.toc();

            patchesToRename.setSize(keys.size());

            forAll(keys, patchI)
            {
                patchesToRename.set
                (
                    patchI,
                    newPatchNames.lookupEntry
                    (
                        keys[patchI],
                        false,
                        false
                    ).clone()
                );
            }
        }
        else
        {
            PtrList<entry> copyPatchesToRename(dict.lookup("newPatchNames"));
            patchesToRename.transfer(copyPatchesToRename);
        }

        forAll(patchesToRename, patchI)
        {
            const word patchName = patchesToRename[patchI].keyword();

            const labelList matchedPatches = mesh_.findPatches(patchName);

            if (matchedPatches.empty())
            {
                Warning << "No matches for " << patchName << " found!!" << endl;
                continue;
            }

            if (!patchesToRename[patchI].isDict())
            {
                Warning << "Cannot rename patch " << patchName << endl;
                Warning << "This is due to incorrect settings! Exitting."
                    << endl;
                return;
            }

            const dictionary pDict = patchesToRename[patchI].dict();

            forAll(matchedPatches, matchI)
            {
                word newName(mesh_.getPatchName(matchedPatches[matchI]));

                pDict.readIfPresent("newName", newName);

                if (newNameToPos.found(newName))
                {
                    // patch with the same name already exists
                    patchToNew[matchedPatches[matchI]] = newNameToPos[newName];
                    continue;
                }

                // add a new patch
                newNameToPos.insert(newName, newPatchI);
                newPatchNames[newPatchI] = newName;
                newPatchTypes[newPatchI] = "wall";

                pDict.readIfPresent("type", newPatchTypes[newPatchI]);

                patchToNew[matchedPatches[matchI]] = newPatchI;
                ++newPatchI;
            }
        }
    }

    word defaultName("walls");
    dict.readIfPresent("defaultName", defaultName);

    word defaultType("wall");
    dict.readIfPresent("defaultType", defaultType);

    bool addDefaultPatch(hasOrphanFaces);
    forAll(patchToNew, patchI)
    {
        if (patchToNew[patchI] == -1)
        {
            addDefaultPatch = true;
            break;
        }
    }

    if (dict.found("defaultName"))
    {
        if (addDefaultPatch)
        {
            newNameToPos.insert(defaultName, newPatchI);
            newPatchNames[newPatchI] = defaultName;
            newPatchTypes[newPatchI] = defaultType;
            ++newPatchI;
        }
    }
    else
    {
        forAll(patchToNew, patchI)
        {
            if (patchToNew[patchI] != -1) continue;

            patchToNew[patchI] = newPatchI;
            newPatchNames[newPatchI] = mesh_.boundaries()[patchI].patchName();
            newPatchTypes[newPatchI] = mesh_.boundaries()[patchI].patchType();
            ++newPatchI;
        }

        if (hasOrphanFaces)
        {
            newNameToPos.insert(defaultName, newPatchI);
            newPatchNames[newPatchI] = defaultName;
            newPatchTypes[newPatchI] = defaultType;
            ++newPatchI;
        }
    }

    if (newPatchI == 0)
    {
        return;
    }

    newPatchNames.setSize(newPatchI);
    newPatchTypes.setSize(newPatchI);

    // start creating new boundary
    DynamicList<DynamicList<label>> newBoundaryFaces;
    DynamicList<label> newBoundaryOwners;
    DynamicList<label> newBoundaryPatches;

    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();

    for (label faceI = nIntFaces; faceI < faces.size(); ++faceI)
    {
        // a face left behind by cell removal without being compacted out
        // of the mesh has no valid owner - skip it, following the same
        // "owner[faceI] < 0 => logically deleted" convention already used
        // elsewhere in this codebase (e.g. extrudeLayer.C)
        if (owner[faceI] < 0)
        {
            continue;
        }

        const label patchI = facePatch[faceI - nIntFaces];

        appendList(newBoundaryFaces, faces[faceI]);
        newBoundaryOwners.append(owner[faceI]);

        if ((patchI == -1) || (patchToNew[patchI] == -1))
        {
            // no patch, or this patch is moved to the default patch
            newBoundaryPatches.append(newPatchI - 1);
        }
        else
        {
            // this patch is renamed
            newBoundaryPatches.append(patchToNew[patchI]);
        }
    }

    // execute the modifier
    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.replaceBoundary
    (
        newPatchNames,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );
    forAll(meshModifier.boundariesAccess(), patchI)
    {
        meshModifier.boundariesAccess()[patchI].patchType() =
            newPatchTypes[patchI];
    }

    Info<< "Finished renaming boundary patches" << endl;
}


void Foam::Module::renameBoundaryPatches::checkEmptyPatches()
{
    polyMeshGenModifier meshModifier(mesh_);

    forAll(mesh_.boundaries(), patchI)
    {
        boundaryPatch& patch = meshModifier.boundariesAccess()[patchI];
        if (patch.patchType() == "empty")
        {
            patch.patchType() = "wall";
        }
    }
}


void Foam::Module::renameBoundaryPatches::checkSymmetryPlanes()
{
    symmetryPlaneOptimisation symmSmother(mesh_);

    symmSmother.optimizeSymmetryPlanes();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::renameBoundaryPatches::renameBoundaryPatches
(
    polyMeshGen& mesh,
    const IOdictionary& meshDict,
    const bool allowEmptyPatches
)
:
    mesh_(mesh),
    meshDict_(meshDict)
{
    if (meshDict.found("renameBoundary"))
    {
        calculateNewBoundary();
    }

    if (!allowEmptyPatches)
    {
        checkEmptyPatches();
    }

    checkSymmetryPlanes();
}


// ************************************************************************* //
