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
    Foam::Module::refineBoundaryLayers

Description
    Refine existing boundary layers

SourceFiles
    refineBoundaryLayers.C
    refineBoundaryLayersFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef refineBoundaryLayers_H
#define refineBoundaryLayers_H

#include "polyMeshGenModifier.H"
#include "meshSurfaceEngine.H"
#include "DynamicList.H"
#include "labelPair.H"
#include "labelListListOps.H"
#include "HashSet.H"
#include "HashTable.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                    Class refineBoundaryLayers Declaration
\*---------------------------------------------------------------------------*/

class refineBoundaryLayers
{
    // Private data

        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- pointer to mesh surface engine
        mutable meshSurfaceEngine* msePtr_;

        //- global number of boundary layers
        label globalNumLayers_;

        //- global thickness ratio
        scalar globalThicknessRatio_;

        //- global maximum thickness of the first layer
        scalar globalMaxThicknessFirstLayer_;

        //- number of boundary layers for user-selected patches
        HashTable<label> numLayersForPatch_;

        //- local thickness ratio for selected patches
        HashTable<scalar> thicknessRatioForPatch_;

        //- local maximum layer thickness for selected patches
        HashTable<scalar> maxThicknessForPatch_;

        //- allow discontinuous layers for patch
        wordHashSet discontinuousLayersForPatch_;

        //- the name of a subset which shall contain the indices
        //- of the cells in the boundary layer
        word cellSubsetName_;

        //- check whether the refinement is already executed
        bool done_;

        //- a flag whether a 2D mesh generation is active or not
        bool is2DMesh_;

        //- shall the layer be refined into two layers with the first one
        //- matching the thickness of the n-1 ones
        bool specialMode_;

        //- information about existing boundary layers at patches
        //- only available layers
        labelListList2D layerAtPatch_;

        //- which patches are part of a single layer
        List<DynamicList<word>> patchesInLayer_;

        //- a containing the number of layers which shall be generated above
        //- a boundary face
        labelList nLayersAtBndFace_;

        //- a list of edges which shall be refined
        DynamicList<edge> splitEdges_;

        //- split edges at point
        labelListList2D splitEdgesAtPoint_;

        //- new vertices for on edges which shall be refined
        labelListList2D newVerticesForSplitEdge_;

        //- a graph containing information which new faces were generated
        //- from an existing face
        labelListList2D facesFromFace_;

        //- a graph containing faces after layer refinement
        //  DynamicList<DynamicList<label>>, not labelListList2D: grown one
        //  new row at a time up to mesh-face scale in generateNewFaces()/
        //  generateNewCells() - needs O(1) amortized outer-array growth,
        //  not List's O(N) per setSize() call (see labelListListOps.H's
        //  appendList overloads).
        DynamicList<DynamicList<label>> newFaces_;


    // Private member functions

        //- Return reference to meshSurfaceEngine
        const meshSurfaceEngine& surfaceEngine() const;

        //- analyse layers to check their topology
        bool analyseLayers();

        //- generate new points on edges, faces and in cells
        void generateNewVertices();

        //- point generated at layerI on split edge seI - fatal if seI
        //- was not generated with enough layers for this (happens when
        //- patches with very different nLayers meet without enough
        //- separation between them - e.g. no real feature edge - for
        //- cfMesh to represent both independently)
        label vertexOnSplitEdge(const label seI, const label layerI) const;

        //- refine a given face and return the new faces
        //- generates new points at cross-split faces
        void refineFace
        (
            const face& f,
            const FixedList<label, 2>& nLayersInDirection,
            DynamicList<DynamicList<label>>& newFaces
        );

        //- generate a matrix of points generated by splitting a face
        //- and return them in the local i, j system of the face
        void sortFacePoints
        (
            const label faceI,
            labelListList2D& facePoints,
            const label transpose = false
        ) const;

        //- generate a matrix of faces generated by splitting a face
        //- and return them in the local i, j, system of the face
        void sortFaceFaces
        (
            const label faceI,
            labelListList2D& faceFaces,
            const label transpose = false
        ) const;

        //- map split edges onto a cell
        void generateNewFaces();

        //- generate new cells for a prism with one boundary face
        void generateNewCellsPrism
        (
            const label cellI,
            List<DynamicList<DynamicList<label>>>& cellsFromCell
        ) const;

        //- a helper function which stores faces generated from
        //- an existing face into new cells
        void storeFacesIntoCells
        (
            const label faceI,
            const bool reverseOrientation,
            const label normalDirection,
            const bool maxCoordinate,
            const label nLayersI,
            const label nLayersJ,
            const label nLayersK,
            List<DynamicList<DynamicList<label>>>& cellsFromCell
        ) const;

        //- generate new cells and add them to the mesh
        void generateNewCells();


    // Nested classes

        class refineEdgeHexCell
        {
            // Private data
                //- label of cell
                const label cellI_;

                //- number of cells in local direction i
                label nLayersI_;

                //- number of cells in local direction j
                label nLayersJ_;

                //- container for new cells
                List<DynamicList<DynamicList<label>>> cellsFromCell_;

                //- const reference to the boundary layer class
                const refineBoundaryLayers& bndLayers_;

                //- faces sorted into directions of a hex shape
                FixedList<label, 6> faceInDirection_;

                //- information about orientation of faces
                //- false means the orientation as expected
                //- true means wrong orientation
                FixedList<bool, 6> faceOrientation_;

                //- points on cross-split faces
                FixedList<labelListList2D, 2> cellPoints_;


            // Private member functions

                //- populate faceInDirection_ and wrongFaceOrientation_
                void determineFacesInDirections();

                //- populate new cells with new faces generated from already
                //- existing faces
                void populateExistingFaces();

                //- generate new internal faces and tore them to new cells
                void generateMissingFaces();

        public:

            //- construct from cell label and the refineBoundaryLayers
            refineEdgeHexCell
            (
                const label cellI,
                const refineBoundaryLayers& ref
            );


            // Public member functions

                inline const List<DynamicList<DynamicList<label>>>&
                newCells() const
                {
                    return cellsFromCell_;
                }
        };

        class refineCornerHexCell
        {
            // Private data

                //- label of cell
                const label cellI_;

                //- number of cells in local direction i
                label nLayersI_;

                //- number of cells in local direction j
                label nLayersJ_;

                //- number of cells in local direction k
                label nLayersK_;

                //- split edge in directions
                FixedList<label, 3> splitEdgeInDirection_;

                //- container for new cells
                List<DynamicList<DynamicList<label>>> cellsFromCell_;

                //- reference to the boundary layer class
                const refineBoundaryLayers& bndLayers_;

                //- faces sorted into directions of a hex shape
                FixedList<label, 6> faceInDirection_;

                //- information about orientation of faces
                //- false means the orientation as expected
                //- true means wrong orientation
                FixedList<bool, 6> faceOrientation_;

                //- points on cross-split faces
                FixedList<labelListList2D, 6> facePoints_;

                //- points inside the cell
                List<List<DynamicList<label>>> cellPoints_;


            // Private member functions

                //- populate faceInDirection_ and wrongFaceOrientation_
                void determineFacesInDirections();

                //- populate new cells with new faces generated from already
                //- existing faces
                void populateExistingFaces();

                //- generate missing points inside the cell
                void generateNewPoints();

                //- generate new internal faces and tore them to new cells
                void generateMissingFaces();


        public:

            //- construct from cell label and the refineBoundaryLayers
            refineCornerHexCell
            (
                const label cellI,
                const refineBoundaryLayers& ref
            );


            // Public member functions

                inline const List<DynamicList<DynamicList<label>>>&
                newCells() const
                {
                    return cellsFromCell_;
                }
        };


    // Private member functions

        //- Disallow bitwise copy construct
        refineBoundaryLayers(const refineBoundaryLayers&);

        //- Disallow bitwise assignment
        void operator=(const refineBoundaryLayers&);


public:

    //- Construct from mesh reference, optionally reading settings from
    //- the given boundaryLayers subDict (dictionary::null for none)
    refineBoundaryLayers
    (
        polyMeshGen& mesh,
        const dictionary& dict = dictionary::null
    );

    //- Destructor
     ~refineBoundaryLayers();


    // Public member functions

        // Settings

        //- read the settings from the boundaryLayers subDict (or an
        //- empty dictionary if the caller has none)
        void read(const dictionary&);

        //- set no refinement flag
        void avoidRefinement();

        //- activate 2D layer refinement
        void activate2DMode();

        //- set the global number of boundary layers
        void setGlobalNumberOfLayers(const label nLayers);

        //- set the global thickness ratio (default is 1)
        void setGlobalThicknessRatio(const scalar thicknessRatio);

        //- set the maximum thickness of the first boundary layer
        void setGlobalMaxThicknessOfFirstLayer(const scalar maxThickness);

        //- set the number of layers for a patch
        //- the settings override the global settings
        void setNumberOfLayersForPatch
        (
            const word& patchName,
            const label nLayers
        );

        //- set the thickness ratio for a patch
        //- it overrides the global settings
        void setThicknessRatioForPatch
        (
            const word& patchName,
            const scalar thicknessRatio
        );

        //- set the maximum thickness of the first layer for a patch
        void setMaxThicknessOfFirstLayerForPatch
        (
            const word& patchName,
            const scalar maxThickness
        );

        //- set whether the settings for a given patch are valid for the
        //- patch only, or whether they extend over th whole sheet
        //- the selected patch belongs to
        //- the default behaviour is to apply the patch settings to the whole
        //- sheet
        void setInteruptForPatch(const word& patchName);

        //- this flag refinement of a layer into a single one matching
        //- the thickness of the n layers. This functionality is intended to
        //- create additional layer that shall be improved by the boundary
        //- layer smoother in order to leave the inside points intact
        void activateSpecialMode();


    // Worker

        //- performs refinement based on the given settings
        void refineLayers();


    // Post information

        //- provide the list of points in the boundary layer
        void pointsInBndLayer(DynamicList<label>&);

        //- create a subset containing points in the boundary layer
        void pointsInBndLayer(const word subsetName="pointsInBndLayer");
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
