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

#include "polyMeshGenFaces.H"
#include "faceIOList.H"
#include "IOPtrList.H"
#include "IOobjectList.H"
#include "faceSet.H"
#include "demandDrivenData.H"
#include "stringListOps.H"
#include "DynamicList.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::Module::polyMeshGenFaces::clearOut() const
{
    deleteDemandDrivenData(ownerPtr_);
    deleteDemandDrivenData(neighbourPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Module::polyMeshGenFaces::polyMeshGenFaces(const Time& runTime)
:
    polyMeshGenPoints(runTime),
    faces_
    (
        IOobject
        (
            "faces",
            runTime.constant(),
            "polyMesh",
            runTime
        ),
        0
    ),
    boundaries_(),
    faceSubsets_(),
    nIntFaces_(0),
    ownerPtr_(nullptr),
    neighbourPtr_(nullptr)
{}


Foam::Module::polyMeshGenFaces::polyMeshGenFaces
(
    const Time& runTime,
    const pointField& points,
    const faceList& faces
)
:
    polyMeshGenPoints(runTime, points),
    faces_
    (
        IOobject
        (
            "faces",
            runTime.constant(),
            "polyMesh",
            runTime
        ),
        faces
    ),
    boundaries_(),
    faceSubsets_(),
    nIntFaces_(0),
    ownerPtr_(nullptr),
    neighbourPtr_(nullptr)
{}


Foam::Module::polyMeshGenFaces::polyMeshGenFaces
(
    const Time& runTime,
    const pointField& points,
    const faceList& faces,
    const wordList& patchNames,
    const labelList& patchStart,
    const labelList& nFacesInPatch
)
:
    polyMeshGenPoints(runTime, points),
    faces_
    (
        IOobject
        (
            "faces",
            runTime.constant(),
            "polyMesh",
            runTime
        ),
        faces
    ),
    boundaries_(),
    faceSubsets_(),
    nIntFaces_(0),
    ownerPtr_(nullptr),
    neighbourPtr_(nullptr)
{
    boundaries_.setSize(patchNames.size());
    forAll(patchNames, patchI)
    {
        boundaries_.set
        (
            patchI,
            new boundaryPatch
            (
                patchNames[patchI],
                "patch",
                nFacesInPatch[patchI],
                patchStart[patchI]
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Module::polyMeshGenFaces::~polyMeshGenFaces()
{
    clearOut();
}


void Foam::Module::polyMeshGenFaces::setMeshLocation
(
    const word& regionName,
    const fileName& meshPath
)
{
    polyMeshGenPoints::setMeshLocation(regionName, meshPath);

    faces_.instance() = resolvedInstance();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::Module::polyMeshGenFaces::faceIsInPatch
(
    const label faceLabel
) const
{
    const label i = boundaries_.size() - 1;
    if (faceLabel >= (boundaries_[i].patchStart()+boundaries_[i].patchSize()))
    {
        return -1;
    }

    forAllReverse(boundaries_, patchI)
    {
        if (faceLabel >= boundaries_[patchI].patchStart())
        {
            return patchI;
        }
    }

    return -1;
}


Foam::wordList Foam::Module::polyMeshGenFaces::patchNames() const
{
    wordList t(boundaries_.size());

    forAll(boundaries_, patchI)
    {
        t[patchI] = boundaries_[patchI].patchName();
    }

    return t;
}


Foam::label Foam::Module::polyMeshGenFaces::getPatchID
(
    const word& patchName
) const
{
    forAll(boundaries_, patchI)
    {
        if (boundaries_.set(patchI))
        {
            if (boundaries_[patchI].patchName() == patchName)
            {
                return patchI;
            }
        }
    }

    // If the code gets here, it implies that the patch was not found.
    // return a -1 in this case
    return -1;
}


Foam::word Foam::Module::polyMeshGenFaces::getPatchName
(
    const label patchID
) const
{
    if ((patchID < 0) || (patchID >= boundaries_.size()))
    {
         FatalErrorInFunction
             << "invalid patch ID supplied"
             << abort(FatalError);
    }

    return boundaries_[patchID].patchName();
}


Foam::labelList Foam::Module::polyMeshGenFaces::findPatches
(
    const word& patchName
) const
{
    if (origPatchNames_.found(patchName))
    {
        const wordList& newNames = origPatchNames_[patchName];

        DynamicList<label> patchIDs;

        forAll(newNames, i)
        {
            patchIDs.append(findMatchingStrings(regExp(newNames[i]), patchNames()));
        }

        return patchIDs;
    }

    const labelList patchIDs =
        findMatchingStrings(regExp(patchName), patchNames());

    if (patchIDs.empty())
    {
        WarningInFunction
            << "Cannot find any patch names matching " << patchName << endl;
    }

    return patchIDs;
}


Foam::label Foam::Module::polyMeshGenFaces::addFaceSubset(const word& setName)
{
    label id = faceSubsetIndex(setName);
    if (id >= 0)
    {
        Warning << "Face subset " << setName << " already exists!" << endl;
        return id;
    }

    id = 0;
    forAllConstIter(Map<meshSubset>, faceSubsets_, it)
    {
        id = Foam::max(id, it.key() + 1);
    }

    faceSubsets_.insert
    (
        id,
        meshSubset(setName, meshSubset::FACESUBSET)
    );

    return id;
}


void Foam::Module::polyMeshGenFaces::removeFaceSubset(const label setI)
{
    if (!faceSubsets_.found(setI))
    {
        return;
    }

    faceSubsets_.erase(setI);
}


Foam::word Foam::Module::polyMeshGenFaces::faceSubsetName
(
    const label setI
) const
{
    Map<meshSubset>::const_iterator it =
        faceSubsets_.find(setI);
    if (it == faceSubsets_.end())
    {
        Warning << "Subset " << setI << " is not a face subset" << endl;
        return word();
    }

    return it().name();
}


Foam::label Foam::Module::polyMeshGenFaces::faceSubsetIndex
(
    const word& setName
) const
{
    forAllConstIter(Map<meshSubset>, faceSubsets_, it)
    {
        if (it().name() == setName)
        {
            return it.key();
        }
    }

    return -1;
}


void Foam::Module::polyMeshGenFaces::read()
{
    polyMeshGenPoints::read();

    // faceCompactIOList (not the plain faceIOList) since with binary
    // writeFormat, OpenFOAM-dev's own faces writer emits the real
    // "faceCompactList" compact encoding (a flat start-offsets labelList
    // + flat elems labelList, not a nested list-of-lists) - only
    // faceCompactIOList's reader recognises and decodes that header;
    // faceIOList only understands the plain "faceList" encoding (which
    // ASCII always uses, since CompactIOList itself falls back to plain
    // encoding for ASCII) and fatals on a genuinely compact binary file
    faceCompactIOList fcs
    (
        IOobject
        (
            "faces",
            resolvedInstance(),
            "polyMesh",
            runTime_,
            IOobject::MUST_READ
        )
    );
    faces_ = fcs;

    deleteDemandDrivenData(ownerPtr_);
    deleteDemandDrivenData(neighbourPtr_);

    ownerPtr_ =
        new labelIOList
        (
            IOobject
            (
                "owner",
                resolvedInstance(),
                "polyMesh",
                runTime_,
                IOobject::MUST_READ
            )
        );

    neighbourPtr_ =
        new labelIOList
        (
            IOobject
            (
                "neighbour",
                resolvedInstance(),
                "polyMesh",
                runTime_,
                IOobject::MUST_READ
            )
        );

    if (neighbourPtr_->size() != ownerPtr_->size())
    {
        neighbourPtr_->setSize(ownerPtr_->size(), -1);
    }

    // read boundary information
    IOPtrList<boundaryPatchBase> patches
    (
        IOobject
        (
            "boundary",
            resolvedInstance(),
            "polyMesh",
            runTime_,
            IOobject::MUST_READ
        )
    );

    boundaries_.setSize(patches.size());
    forAll(patches, patchI)
    {
        boundaries_.set
        (
            patchI,
            new boundaryPatch
            (
                patches[patchI].patchName(),
                patches[patchI].patchType(),
                patches[patchI].patchSize(),
                patches[patchI].patchStart()
            )
        );
    }

    nIntFaces_ = boundaries_[0].patchStart();

    // read face subsets
    IOobjectList allSets
    (
        runTime_,
        resolvedInstance(),
        "polyMesh/sets"
    );

    wordList setNames = allSets.names("faceSet");
    forAll(setNames, setI)
    {
        IOobject* obj = allSets.lookup(setNames[setI]);

        faceSet fSet(*obj);
        const labelList content = fSet.toc();
        const label id = addFaceSubset(setNames[setI]);

        faceSubsets_[id].updateSubset(content);
    }
}


void Foam::Module::polyMeshGenFaces::write() const
{
    polyMeshGenPoints::write();

    faces_.write();

    if (!ownerPtr_ || !neighbourPtr_)
    {
        calculateOwnersAndNeighbours();
    }
    ownerPtr_->write();
    neighbourPtr_->write();

    // write boundary data
    PtrList<boundaryPatchBase> ptchs(boundaries_.size());

    forAll(boundaries_, patchI)
    {
        dictionary dict;
        dict.add("type", boundaries_[patchI].patchType());
        dict.add("nFaces", boundaries_[patchI].patchSize());
        dict.add("startFace", boundaries_[patchI].patchStart());
        ptchs.set
        (
            patchI,
            boundaryPatchBase::New
            (
                boundaries_[patchI].patchName(),
                dict
            )
        );
    }

    IOPtrList<boundaryPatchBase> patches
    (
        IOobject
        (
            "boundary",
            resolvedInstance(),
            "polyMesh",
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ptchs
    );

    patches.write();

    // write face subsets
    forAllConstIter(Map<meshSubset>, faceSubsets_, setIt)
    {
        faceSet set
        (
            IOobject
            (
                setIt().name(),
                resolvedInstance(),
                "polyMesh/sets",
                runTime_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );

        DynamicList<label> containedElements;
        setIt().containedElements(containedElements);

        forAll(containedElements, i)
        {
            set.insert(containedElements[i]);
        }
        set.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
