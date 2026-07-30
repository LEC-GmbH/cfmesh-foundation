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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::Module::faceListPMG&
Foam::Module::polyMeshGenFaces::faces() const
{
    return faces_;
}


inline Foam::label Foam::Module::polyMeshGenFaces::nInternalFaces() const
{
    if (!(ownerPtr_ && neighbourPtr_))
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateOwnersAndNeighbours();
    }

    return nIntFaces_;
}


inline const Foam::labelList& Foam::Module::polyMeshGenFaces::owner() const
{
    if (!ownerPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateOwnersAndNeighbours();
    }

    return *ownerPtr_;
}


inline const Foam::labelList& Foam::Module::polyMeshGenFaces::neighbour() const
{
    if (!neighbourPtr_)
    {
        if (omp_in_parallel())
        {
            FatalErrorInFunction
                << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        }

        calculateOwnersAndNeighbours();
    }

    return *neighbourPtr_;
}


inline const Foam::PtrList<Foam::Module::boundaryPatch>&
Foam::Module::polyMeshGenFaces::boundaries() const
{
    return boundaries_;
}


inline void Foam::Module::polyMeshGenFaces::addFaceToSubset
(
    const label setID,
    const label faceI
)
{
    Map<meshSubset>::iterator it = faceSubsets_.find(setID);
    if (it == faceSubsets_.end())
    {
        return;
    }

    it().addElement(faceI);
}


inline void Foam::Module::polyMeshGenFaces::removeFaceFromSubset
(
    const label setI,
    const label faceI)
{
    Map<meshSubset>::iterator it = faceSubsets_.find(setI);
    if (it == faceSubsets_.end())
    {
        return;
    }

    it().removeElement(faceI);
}


inline void Foam::Module::polyMeshGenFaces::faceInSubsets
(
    const label faceI,
    DynamicList<label>& faceSubsets
) const
{
    faceSubsets.clear();

    forAllConstIter(Map<meshSubset>, faceSubsets_, it)
    {
        if (it().found(faceI))
        {
            faceSubsets.append(it.key());
        }
    }
}


inline void Foam::Module::polyMeshGenFaces::faceSubsetIndices
(
    DynamicList<label>& indices
) const
{
    indices.clear();

    forAllConstIter(Map<meshSubset>, faceSubsets_, it)
    {
        indices.append(it.key());
    }
}


template<class ListType>
inline void Foam::Module::polyMeshGenFaces::facesInSubset
(
    const label setI,
    ListType& faceLabels
) const
{
    faceLabels.clear();

    Map<meshSubset>::const_iterator it = faceSubsets_.find(setI);
    if (it == faceSubsets_.end())
    {
        return;
    }

    it().containedElements(faceLabels);
}


template<class ListType>
inline void Foam::Module::polyMeshGenFaces::updateFaceSubsets
(
    const ListType& newFaceLabels
)
{
    forAllIter(Map<meshSubset>, faceSubsets_, it)
    {
        it().updateSubset(newFaceLabels);
    }
}


inline void Foam::Module::polyMeshGenFaces::updateFaceSubsets
(
    const labelListList2D& newFacesForFace
)
{
    forAllIter(Map<meshSubset>, faceSubsets_, it)
    {
        it().updateSubset(newFacesForFace);
    }
}


inline const Foam::HashTable<Foam::wordList>&
Foam::Module::polyMeshGenFaces::origPatchNamesMap() const
{
    return origPatchNames_;
}


inline void Foam::Module::polyMeshGenFaces::setOrigPatchNamesMap
(
    const HashTable<wordList>& origPatchNames
)
{
    origPatchNames_ = origPatchNames;
}


inline const Foam::HashTable<Foam::word>&
Foam::Module::polyMeshGenFaces::origPatchTypesMap() const
{
    return origPatchTypes_;
}


inline void Foam::Module::polyMeshGenFaces::setOrigPatchTypesMap
(
    const HashTable<word>& origPatchTypes
)
{
    origPatchTypes_ = origPatchTypes;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
