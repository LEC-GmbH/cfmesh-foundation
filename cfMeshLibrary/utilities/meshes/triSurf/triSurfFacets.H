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
    Foam::Module::triSurfFacets

Description
    Facets for the triangulated surface

SourceFiles
    triSurfFacets.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfFacets_H
#define triSurfFacets_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfAddressing.H"
#include "meshSubset.H"
#include "geometricSurfacePatchList.H"
#include "labelledTri.H"
#include "DynamicList.H"
#include "Map.H"
#include "HashTable.H"

namespace Foam
{
namespace Module
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfFacets
{
protected:

    // Protected data

        //- list of triangles
        DynamicList<labelledTri> triangles_;

        //- list of boundary patches and their properties
        geometricSurfacePatchList patches_;

        //- map of point subsets
        Map<meshSubset> facetSubsets_;

        //- map from a patch's name before triSurfacePatchManipulator
        //- split it by feature edges to the literal names of the
        //- patches it was split into; empty if this surface was not
        //- produced by such a split. Consulted by findPatches() so
        //- patch-name-keyed meshDict entries (patchCellSize,
        //- localRefinement, keepCellsIntersectingPatches, ...) written
        //- against the pre-split name keep resolving correctly.
        HashTable<wordList> origPatchNames_;

        //- geometricType of each pre-split patch, keyed by its
        //- (pre-split) name; empty if this surface was not produced by
        //- a split. Consulted by renameBoundaryPatches to collapse
        //- every split patch back to its original name/type.
        HashTable<word> origPatchTypes_;


    // Disallow bitwise assignment

        void operator=(const triSurfFacets&);

        triSurfFacets(const triSurfFacets&);


public:

    // Constructors

        //- Null constructor
        triSurfFacets();

        //- Construct from components without the boundary
        triSurfFacets(const DynamicList<labelledTri>& triangles);

        //- Construct from components
        triSurfFacets
        (
            const DynamicList<labelledTri>& triangles,
            const geometricSurfacePatchList& patches
        );


    //- Destructor
    ~triSurfFacets() = default;


    // Member functions

        //- return the number of triangles
        inline label size() const;

        //- access to facets
        inline const DynamicList<labelledTri>& facets() const;

        //- access to patches
        inline const geometricSurfacePatchList& patches() const;

        //- return list of patches in the boundary
        wordList patchNames() const;

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

        //- append a triangle to the end of the list
        inline void appendTriangle(const labelledTri& tria);

        //- point subsets
        label addFacetSubset(const word&);
        void removeFacetSubset(const label);
        word facetSubsetName(const label) const;
        label facetSubsetIndex(const word&) const;
        inline void addFacetToSubset(const label, const label);
        inline void removeFacetFromSubset(const label, const label);
        inline void facetInSubsets(const label, DynamicList<label>&) const;
        inline void facetSubsetIndices(DynamicList<label>&) const;
        template<class ListType>
        inline void facetsInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updateFacetsSubsets(const ListType&);


    // Operators

        //- access to a triangle
        inline const labelledTri& operator[](const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triSurfFacetsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
