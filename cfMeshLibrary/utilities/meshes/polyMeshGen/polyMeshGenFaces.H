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

Class
    Foam::Module::polyMeshGenFaces

Description
    Mesh with geometry and faces. It can be used for volume and surface meshes

SourceFiles
    polyMeshGenFaces.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenFaces_H
#define polyMeshGenFaces_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenPoints.H"
#include "faceListPMG.H"
#include "labelIOList.H"
#include "boundaryPatch.H"
#include "DynamicList.H"
#include "HashTable.H"

namespace Foam
{
namespace Module
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenFaces
:
    public polyMeshGenPoints
{
protected:

    // Protected Data

        //- list of faces
        faceListPMG faces_;

        //- boundary data
        PtrList<boundaryPatch> boundaries_;

        //- face subsets
        Map<meshSubset> faceSubsets_;

        //- map from a patch's name before triSurfacePatchManipulator
        //- split it by feature edges to the literal names of the
        //- patches it was split into; empty if this mesh's surface was
        //- not produced by such a split. Consulted by findPatches() so
        //- patch-name-keyed meshDict entries (boundaryLayers.
        //- patchBoundaryLayers, ...) written against the pre-split name
        //- keep resolving correctly. Copied from the originating
        //- triSurf's own map (see triSurfFacets::origPatchNamesMap())
        //- when boundary patches are created during mesh extraction.
        HashTable<wordList> origPatchNames_;

        //- geometricType of each pre-split patch, keyed by its
        //- (pre-split) name; empty if this mesh's surface was not
        //- produced by a split. Consulted by renameBoundaryPatches to
        //- collapse every split patch back to its original name/type.
        //- Copied from the originating triSurf's own map (see
        //- triSurfFacets::origPatchTypesMap()) when boundary patches are
        //- created during mesh extraction.
        HashTable<word> origPatchTypes_;


    // Addressing data

        //- number of internal faces, owner and neighbour
        mutable label nIntFaces_;
        mutable labelIOList* ownerPtr_;
        mutable labelIOList* neighbourPtr_;


    // Private Member Functions

        //- calculate owner and neighbour addressing
        virtual void calculateOwnersAndNeighbours() const = 0;

        //- clear all pointer data
        void clearOut() const;


        //- No copy construct
        polyMeshGenFaces(const polyMeshGenFaces&) = delete;

        //- No copy assignment
        void operator=(const polyMeshGenFaces&) = delete;


public:

    friend class polyMeshGenModifier;

    // Constructors

        //- Null constructor
        polyMeshGenFaces(const Time&);

        //- Construct from components without the boundary
        polyMeshGenFaces
        (
            const Time&,
            const pointField& points,
            const faceList& faces
        );

        //- Construct from components with the boundary
        polyMeshGenFaces
        (
            const Time&,
            const pointField& points,
            const faceList& faces,
            const wordList& patchNames,
            const labelList& patchStart,
            const labelList& nFacesInPatch
        );


    //- Destructor, calls clearOut()
    virtual ~polyMeshGenFaces();


    // Member functions

        //- set the region/mesh-path this mesh is read from/written to -
        //- must be called (if at all) before any read()/write(), since it
        //- also fixes up faces_'s already-constructed IOobject instance
        void setMeshLocation
        (
            const word& regionName,
            const fileName& meshPath
        );

        //- access to faces
        inline const faceListPMG& faces() const;

        //- return number of internal faces
        inline label nInternalFaces() const;

        //- owner and neighbour cells for faces
        inline const labelList& owner() const;
        inline const labelList& neighbour() const;


    // Boundary data

        //- ordinary boundaries
        inline const PtrList<boundaryPatch>& boundaries() const;

        //- return patch label for the given face label
        label faceIsInPatch(const label faceLabel) const;

        //- return list of patches in the boundary
        wordList patchNames() const;

        //- return the index of a patch given its name
        label getPatchID(const word& patchName) const;

        //- return the name of a patch given its ID
        word getPatchName(const label patchID) const;

        //- return a list of patch indices corresponding to the given
        // name, expanding regular expressions and pre-split patch names
        // (see origPatchNamesMap())
        labelList findPatches(const word& patchName) const;

        //- access to the pre-split -> post-split patch name map
        //- (see origPatchNames_)
        inline const HashTable<wordList>& origPatchNamesMap() const;

        //- set the pre-split -> post-split patch name map
        inline void setOrigPatchNamesMap(const HashTable<wordList>&);

        //- access to the pre-split patch type map (see origPatchTypes_)
        inline const HashTable<word>& origPatchTypesMap() const;

        //- set the pre-split patch type map
        inline void setOrigPatchTypesMap(const HashTable<word>&);


    // Subsets

        label addFaceSubset(const word&);
        void removeFaceSubset(const label);
        word faceSubsetName(const label) const;
        label faceSubsetIndex(const word&) const;
        inline void addFaceToSubset(const label, const label);
        inline void removeFaceFromSubset(const label, const label);
        inline void faceInSubsets(const label, DynamicList<label>&) const;
        inline void faceSubsetIndices(DynamicList<label>&) const;
        template<class ListType>
        inline void facesInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updateFaceSubsets(const ListType&);
        inline void updateFaceSubsets(const labelListList2D&);

        //- Read mesh
        void read();

        //- Write mesh
        void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenFacesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
