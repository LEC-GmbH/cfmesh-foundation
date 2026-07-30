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
    Foam::Module::partTetMesh

Description
    Mesh smoothing without any topological changes

SourceFiles
    partTetMesh.C

\*---------------------------------------------------------------------------*/

#ifndef partTetMesh_H
#define partTetMesh_H

#include "boolList.H"
#include "DynamicList.H"
#include "labelListListOps.H"
#include "partTet.H"
#include "HashSet.H"
#include "labelledPoint.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
// Forward declarations
class polyMeshGen;

/*---------------------------------------------------------------------------*\
                         Class partTetMesh Declaration
\*---------------------------------------------------------------------------*/

class partTetMesh
{
    // Private data

        //- reference to the original mesh
        polyMeshGen& origMesh_;

        //- points in the tet mesh
        DynamicList<point> points_;

        //- tetrahedra making the mesh
        DynamicList<partTet> tets_;

        //- label of node in the polyMeshGen
        DynamicList<label> nodeLabelInOrigMesh_;

        //- shall a node be used for smoothing or not
        DynamicList<direction> smoothVertex_;

        //- addressing data
        labelListList2D pointTets_;

        //- internal point ordering for parallel runs
        mutable labelListList2D* internalPointsOrderPtr_;

        //- boundary point ordering for parallel runs
        mutable labelListList2D* boundaryPointsOrderPtr_;


    // Private member functions

        //- create points and tets
        void createPointsAndTets
        (
            const List<direction>& usedCells,
            const boolList& lockedPoints
        );

        //- create ordering of internal points for parallel execution
        void createSMOOTHPointsOrdering() const;

        //- create order of boundary points for parallel execution
        void createBOUNDARYPointsOrdering() const;


public:

    // Constructors

        //- construct from polyMeshGen and locked points
        partTetMesh(polyMeshGen& mesh, const DynamicList<label>& lockedPoints);

        //- construct from polyMeshGen, locked points and the number of layers
        //- from the boundary
        partTetMesh
        (
            polyMeshGen& mesh,
            const DynamicList<label>& lockedPoints,
            const direction nLayers = 2
        );

        //- construct from polyMeshGen, bad faces and the number
        //- of additional layers
        partTetMesh
        (
            polyMeshGen& mesh,
            const DynamicList<label>& lockedPoints,
            labelHashSet& badFaces,
            const direction additionalLayers = 0
        );


    // Enumerators

        enum vertexTypes
        {
            NONE = 0,
            SMOOTH = 1,
            FACECENTRE = 2,
            CELLCENTRE = 4,
            BOUNDARY = 16,
            LOCKED = 32
        };


    //- Destructor
    ~partTetMesh();


    // Member functions

        //- access to points, tets and other data
        inline const DynamicList<point>& points() const
        {
            return points_;
        }

        inline const DynamicList<partTet>& tets() const
        {
            return tets_;
        }

        inline const labelListList2D& pointTets() const
        {
            return pointTets_;
        }

        inline const DynamicList<direction>& smoothVertex() const
        {
            return smoothVertex_;
        }

        inline const DynamicList<label>& nodeLabelInOrigMesh() const
        {
            return nodeLabelInOrigMesh_;
        }


    // Access to point ordering for Gauss-Seidel type of iteration

        //- return vertex ordering which can be executed in parallel
        //- the points in each row can be treated in parallel
        //- make sure that points added to different rows of this graph are
        //- not treated concurrently
        const labelListList2D& internalPointOrdering() const;

        //- return vertex ordering for boundary points
        const labelListList2D& boundaryPointOrdering() const;


    // Modifiers

        //- move the vertex to a new position
        void updateVertex(const label pointI, const point& newP);

        //- move vertices to their new positions
        //- intended for SMP parallelisation
        void updateVerticesSMP(const List<DynamicList<labelledPoint>>&);

        //- updates the vertices of the original polyMeshGen
        void updateOrigMesh(boolList* changedFacePtr = nullptr);

        //- creates polyMeshGen from this partTetMesh
        void createPolyMesh(polyMeshGen& pmg) const;

        //- add the LOCKED flag to mesh points which shall not move
        template<class labelListType>
        void lockPoints(const labelListType& l)
        {
            forAll(l, pI)
            {
                smoothVertex_[l[pI]] |= LOCKED;
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
