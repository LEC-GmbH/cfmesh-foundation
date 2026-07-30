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
    Foam::Module::meshSurfaceEngine

Description
    Calculates surface of the mesh

SourceFiles
    meshSurfaceEngine.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEngine_H
#define meshSurfaceEngine_H

#include "polyMeshGenModifier.H"
#include "SubList.H"
#include "boolList.H"
#include "labelListListOps.H"
#include "triangle.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Module
{
// Forward declarations
class meshSurfaceEngineModifier;

/*---------------------------------------------------------------------------*\
                      Class meshSurfaceEngine Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEngine
{
    // Private data

        //- reference to the mesh
        polyMeshGen& mesh_;

        //- number of active patch
        const label activePatch_;

        //- boundary points
        mutable labelList* boundaryPointsPtr_;

        //- boundary faces
        mutable faceList::subList* boundaryFacesPtr_;

        //- patches boundary faces are in
        mutable labelList* boundaryFacePatchPtr_;

        //- face owners
        mutable labelList* boundaryFaceOwnersPtr_;

        //- point faces addressing
        mutable labelListList2D* pointFacesPtr_;
        mutable labelListList2D* pointInFacePtr_;

        //- pointBoundaryPoint addressing
        mutable labelList* bppPtr_;

        //- point points addressing
        mutable labelListList2D* pointPointsPtr_;

        //- edges
        mutable edgeList* edgesPtr_;

        //- boundary point-edges addressing
        mutable labelListList2D* bpEdgesPtr_;

        //- edge faces addressing
        mutable labelListList2D* edgeFacesPtr_;

        //- face edges addressing
        mutable labelListList2D* faceEdgesPtr_;

        //- edge-patches addressing
        mutable labelListList2D* edgePatchesPtr_;

        //- face-faces addressing
        mutable labelListList2D* faceFacesPtr_;

        //- point normals
        mutable vectorField* pointNormalsPtr_;

        //- face normals
        mutable vectorField* faceNormalsPtr_;

        //- face centres
        mutable vectorField* faceCentresPtr_;


    // Private member functions

        //- calculate boundary nodes, faces and addressing
        void calculateBoundaryFaces() const;
        void calculateBoundaryOwners() const;
        void calculateBoundaryNodes() const;

        void calculateBoundaryFacePatches() const;

        void calculatePointFaces() const;

        void calculatePointPoints() const;

        void calculatePointNormals() const;
        void calculateFaceNormals() const;
        void calculateFaceCentres() const;

        void calculateEdgesAndAddressing() const;
        void calculateFaceEdgesAddressing() const;
        void calculateEdgeFacesAddressing() const;
        void calculateEdgePatchesAddressing() const;
        void calculateFaceFacesAddressing() const;

        void clearOut();

        //- Disallow default bitwise copy construct
        meshSurfaceEngine(const meshSurfaceEngine&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEngine&);


public:

        friend class meshSurfaceEngineModifier;


    // Constructors

        //- Construct from mesh
        meshSurfaceEngine(polyMeshGen& mesh);

        //- Construct from mesh and patch
        meshSurfaceEngine(polyMeshGen &mesh, const label patchI);

    //- Destructor
    ~meshSurfaceEngine();


    // Member Functions

        inline const polyMeshGen& mesh() const;
        inline const pointFieldPMG& points() const;
        inline const faceListPMG& faces() const;
        inline const cellListPMG& cells() const;

        inline const labelList& boundaryPoints() const;

        inline const faceList::subList& boundaryFaces() const;

        //- patch label for each boundary face
        inline const labelList& boundaryFacePatches() const;

        inline const labelList& faceOwners() const;

        inline const labelListList2D& pointFaces() const;
        inline const labelListList2D& pointInFaces() const;

        inline const labelListList2D& pointPoints() const;

        inline const vectorField& pointNormals() const;

        inline const vectorField& faceNormals() const;

        inline const vectorField& faceCentres() const;

        inline const labelList& bp() const;

        inline const edgeList& edges() const;

        inline const labelListList2D& boundaryPointEdges() const;

        inline const labelListList2D& edgeFaces() const;

        inline const labelListList2D& faceEdges() const;

        inline const labelListList2D& edgePatches() const;

        inline const labelListList2D& faceFaces() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshSurfaceEngineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
