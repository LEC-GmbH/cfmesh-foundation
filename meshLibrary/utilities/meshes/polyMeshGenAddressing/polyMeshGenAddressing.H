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
    Foam::Module::polyMeshGenAddressing

Description
    Cell-face mesh analysis engine

SourceFiles
    polyMeshGenAddressing.C
    polyMeshGenAddressingClear.C
    polyMeshGenAddressingCellCells.C
    polyMeshGenAddressingEdgeCells.C
    polyMeshGenAddressingPointCells.C
    polyMeshGenAddressingCells.C
    polyMeshGenAddressingEdgeFaces.C
    polyMeshGenAddressingPointFaces.C
    polyMeshGenAddressingCellEdges.C
    polyMeshGenAddressingPointEdges.C
    polyMeshGenAddressingPointPoints.C
    polyMeshGenAddressingEdges.C
    polyMeshGenAddressingCellCentresAndVols.C
    polyMeshGenAddressingFaceCentresAndAreas.C
    polyMeshGenAddressingUpdateGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenAddressing_H
#define polyMeshGenAddressing_H

#include "edgeList.H"
#include "polyMeshGenCells.H"
#include "boolList.H"
#include "HashSet.H"
#include "DynList.H"
#include "Map.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{

/*---------------------------------------------------------------------------*\
                    Class polyMeshGenAddressing Declaration
\*---------------------------------------------------------------------------*/

class polyMeshGenAddressing
{
    // Mesh data

        //- reference to the mesh
        const polyMeshGenCells& mesh_;


    // Shapes

        //- Edges
        mutable edgeList* edgesPtr_;


    // Connectivity

        //- Cell-cells
        mutable VRWGraph* ccPtr_;

        //- Point-cells
        mutable VRWGraph* pcPtr_;

        //- Edge-faces
        mutable VRWGraph* efPtr_;

        //- Point-faces
        mutable VRWGraph* pfPtr_;

        //- Cell-edges
        mutable VRWGraph* cePtr_;

        //- Face-edges
        mutable VRWGraph* fePtr_;

        //- Point-edges
        mutable VRWGraph* pePtr_;

        //- Point-points
        mutable VRWGraph* ppPtr_;

        //- Cell-points
        mutable VRWGraph* cpPtr_;


    // Geometric data

        //- Cell centres
        mutable vectorField* cellCentresPtr_;

        //- Face centres
        mutable vectorField* faceCentresPtr_;

        //- Cell volumes
        mutable scalarField* cellVolumesPtr_;

        //- Face areas
        mutable vectorField* faceAreasPtr_;


    // Private member functions

        //- Disallow construct as copy
        polyMeshGenAddressing(const polyMeshGenAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const polyMeshGenAddressing&);


    // Private member functions to calculate demand driven data

    // Topological calculations

        //- Calculate cell-cell addressing
        void calcCellCells() const;

        //- Calculate point-cell addressing
        void calcPointCells() const;

        //- Calculate edge-face addressing
        void calcEdgeFaces() const;

        //- Calculate face-edge addressing
        void calcFaceEdges() const;

        //- Calculate point-face addressing
        void calcPointFaces() const;

        //- Calculate edge list
        void calcCellEdges() const;

        //- Calculate point-edge addressing
        void calcPointEdges() const;

        //- Calculate point-point addressing
        void calcPointPoints() const;

        //- Calculate cell-point addressing
        void calcCellPoints() const;

        //- During edge calculation, a larger set of data is assembled.
        // Create and destroy as a set, using clearOutEdges()
        void calcEdges() const;
        void clearOutEdges();


    // Geometrical calculations

        //- Calculate face centres and areas
        void calcFaceCentresAndAreas() const;
        void makeFaceCentresAndAreas
        (
            const pointFieldPMG& p,
            vectorField& fCtrs,
            vectorField& fAreas
        ) const;

        //- Calculate cell centres and volumes
        void calcCellCentresAndVols() const;
        void makeCellCentresAndVols
        (
            const vectorField& fCtrs,
            const vectorField& fAreas,
            vectorField& cellCtrs,
            scalarField& cellVols
        ) const;

        //- Calculate edge vectors
        void calcEdgeVectors() const;


    // Disallow bitwise construct

        //- Default constructor
        polyMeshGenAddressing();


public:

    // Static data

        ClassName("polyMeshGenAddressing");

    //- Construct from components
    polyMeshGenAddressing(const polyMeshGenCells& mesh);


    //- Destructor
    ~polyMeshGenAddressing();


    // Member Functions

        // Primitive mesh data

            //- Return mesh edges
            const edgeList& edges() const;


        // Return mesh connectivity

            const VRWGraph& cellCells() const;
            // faceCells given as faceOwner and faceNeighbour
            const VRWGraph& pointCells() const;

            // faceFaces considered unnecessary
            const VRWGraph& edgeFaces() const;
            const VRWGraph& pointFaces() const;

            const VRWGraph& cellEdges() const;
            const VRWGraph& faceEdges() const;
            // edgeEdges considered unnecessary
            const VRWGraph& pointEdges() const;
            const VRWGraph& pointPoints() const;
            const VRWGraph& cellPoints() const;


        // Geometric data (raw!)

            const vectorField& cellCentres() const;
            const vectorField& faceCentres() const;
            const scalarField& cellVolumes() const;
            const vectorField& faceAreas() const;


    //  Storage management

        //- Clear geometry
        void clearGeom();

        //- Clear topological data
        void clearAddressing();

        //- Clear all geometry and addressing unnecessary for CFD
        void clearOut();

        //- Clear everything primitive, geometry and addressing
        void clearAll();

        //- Update geometry data
        void updateGeometry(const boolList& changedFace);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
