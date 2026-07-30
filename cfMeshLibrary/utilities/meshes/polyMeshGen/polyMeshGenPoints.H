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
    Foam::Module::polyMeshGenPoints

Description
    Mesh with geometry information

SourceFiles
    polyMeshGenPoints.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenPoints_H
#define polyMeshGenPoints_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Time.H"
#include "meshSubset.H"
#include "pointFieldPMG.H"
#include "DynamicList.H"
#include "Map.H"
#include "polyMesh.H"

namespace Foam
{
namespace Module
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenPoints
{
protected:

    // Protected data

        //- reference to the Time registry
        const Time& runTime_;

        //- region name and mesh path this mesh is read from/written to
        //- (defaults reproduce the historical constant/polyMesh location)
        word regionName_;
        fileName meshPath_;

        //- list of vertices
        pointFieldPMG points_;

        //- map of point subsets
        Map<meshSubset> pointSubsets_;


    // Protected member functions

        //- instance folding in region + mesh path - local() stays the fixed
        //- "polyMesh" literal everywhere, since local() is immutable once
        //- an IOobject/regIOobject is constructed, unlike instance()
        inline fileName resolvedInstance() const
        {
            return
                runTime_.constant()
               /meshPath_
               /(regionName_ == polyMesh::defaultRegion ? word::null : regionName_);
        }


    // Disallow bitwise assignment

        void operator=(const polyMeshGenPoints&);

        polyMeshGenPoints(const polyMeshGenPoints&);


public:

    friend class polyMeshGenModifier;


    // Constructors

        //- Null constructor
        polyMeshGenPoints(const Time&);

        //- Construct from components without the boundary
        polyMeshGenPoints
        (
            const Time&,
            const pointField& points
        );


    //- Destructor
    ~polyMeshGenPoints() = default;


    // Member functions

        //- access to Time
        inline const Time& returnTime() const;

        //- set the region/mesh-path this mesh is read from/written to -
        //- must be called (if at all) before any read()/write(), since it
        //- also fixes up points_'s already-constructed IOobject instance
        void setMeshLocation
        (
            const word& regionName,
            const fileName& meshPath
        );

        //- access to points
        inline const pointFieldPMG& points() const;

        //- non-const access to points
        inline pointFieldPMG& points();

        //- append a vertex to the end of the list
        inline void appendVertex(const point& p);

        //- point subsets
        label addPointSubset(const word&);
        void removePointSubset(const label);
        word pointSubsetName(const label) const;
        label pointSubsetIndex(const word&) const;
        inline void addPointToSubset(const label, const label);
        inline void removePointFromSubset(const label, const label);
        inline void pointInSubsets(const label, DynamicList<label>&) const;
        inline void pointSubsetIndices(DynamicList<label>&) const;
        template<class ListType>
        inline void pointsInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updatePointSubsets(const ListType&);

        //- Read mesh
        void read();

        //- Write mesh
        void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Module
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenPointsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
