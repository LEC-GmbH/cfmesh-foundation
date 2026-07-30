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
    Foam::Module::polyMeshGenModifier

Description
    Modifier for polyMeshGen

SourceFiles
    polyMeshGenModifier.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenModifier_H
#define polyMeshGenModifier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGen.H"
#include "boolList.H"
#include "labelListListOps.H"
#include "demandDrivenData.H"

namespace Foam
{
namespace Module
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenModifier
{
    // Private data

        //- reference to the mesh
        polyMeshGen& mesh_;

        //- helper data for adding cells
        labelListList2D* pointFacesPtr_;


    // Private member functions

        //- calculate and return point faces
        inline void calculatePointFaces()
        {
            const faceListPMG& faces = mesh_.faces();

            pointFacesPtr_ = new labelListList2D();
            labelListList2D& pointFaces = *pointFacesPtr_;

            reverseAddressing(pointFaces, mesh_.points().size(), faces);
        }


protected:

    // Protected member functions

        labelListList2D& pointFaces()
        {
            if (!pointFacesPtr_)
            {
                calculatePointFaces();
            }

            return *pointFacesPtr_;
        }


public:

    //- Construct from the reference to the mesh
    polyMeshGenModifier(polyMeshGen& mesh)
    :
        mesh_(mesh),
        pointFacesPtr_(nullptr)
    {
        //mesh_.clearOut();
    }


    //- Destructor
    ~polyMeshGenModifier()
    {
        this->clearOut();
    }


    // Member functions

        //- access to mesh points
        inline pointFieldPMG& pointsAccess()
        {
            return mesh_.points_;
        }

        //- access to mesh faces
        inline faceListPMG& facesAccess()
        {
            return mesh_.faces_;
        }

        //- access to cells
        inline cellListPMG& cellsAccess()
        {
            return mesh_.cells_;
        }

        //- access to boundary data
        inline PtrList<boundaryPatch>& boundariesAccess()
        {
            return mesh_.boundaries_;
        }

        //- functions which change the mesh
        //- reorder boundary faces
        void reorderBoundaryFaces();

        //- remove unused vertices
        void removeUnusedVertices();

        //- remove faces
        void removeFaces(const boolList& removeFace);

        //- remove cells
        void removeCells
        (
            const boolList& removeCell
        );

        //- add cells (vertices must be added)
        void addCells(const DynamicList<faceList>& cellFaces);
        void addCell(const faceList& cellFaces);

        //- replace the boundary with new boundary faces
        void replaceBoundary
        (
            const wordList& patchNames,
            const DynamicList<DynamicList<label>>& boundaryFaces,
            const DynamicList<label>& faceOwners,
            const DynamicList<label>& facePatches
        );

        //- reorder the cells and faces to reduce the matrix bandwidth
        void renumberMesh();

        //- clear out unnecessary data (pointFacesPtr_);
        inline void clearOut()
        {
            deleteDemandDrivenData(pointFacesPtr_);
        }

        //- clear out all allocated data
        inline void clearAll()
        {
            clearOut();
            mesh_.clearOut();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
