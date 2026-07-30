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

#include "Map.H"
#include "partTetMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Module::partTetMeshSimplex::partTetMeshSimplex
(
    const partTetMesh& tm,
    const label pI
)
{
    const DynamicList<point>& points = tm.points();
    const DynamicList<partTet>& tets = tm.tets();
    const labelListList2D& pt = tm.pointTets();

    tets_.setSize(pt[pI].size());
    label counter(0);

    Map<label> addr(2*pt[pI].size());
    forAllRow(pt, pI, tetI)
    {
        const partTet& tet = tets[pt[pI][tetI]];
        for (label i = 0; i < 4; ++i)
        {
            const label tpI = tet[i];
            if (!addr.found(tpI))
            {
                addr.insert(tpI, counter);
                pts_.append(points[tpI]);
                ++counter;
            }
        }

        const label pos = tet.whichPosition(pI);
        switch(pos)
        {
            case 0:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet.b()],
                        addr[tet.d()],
                        addr[tet.c()],
                        addr[tet.a()]
                    );
            } break;
            case 1:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet.a()],
                        addr[tet.c()],
                        addr[tet.d()],
                        addr[tet.b()]
                    );
            } break;
            case 2:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet.a()],
                        addr[tet.d()],
                        addr[tet.b()],
                        addr[tet.c()]
                    );
            } break;
            case 3:
            {
                tets_[tetI] =
                    partTet
                    (
                        addr[tet.a()],
                        addr[tet.b()],
                        addr[tet.c()],
                        addr[tet.d()]
                    );
            } break;
            default:
            {
                FatalErrorInFunction
                    << "Point " << pI << " is not present in tet" << tet
                    << abort(FatalError);
            }
        }
    }
}


Foam::Module::partTetMeshSimplex::partTetMeshSimplex
(
    const DynamicList<point>& pts,
    const DynamicList<partTet>& tets,
    const label pointI
)
:
    pts_(pts)
{
    tets_.setSize(tets.size());

    forAll(tets, tetI)
    {
        const partTet& tet = tets[tetI];

        const label pos = tet.whichPosition(pointI);

        switch(pos)
        {
            case 0:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.b(),
                        tet.d(),
                        tet.c(),
                        tet.a()
                    );
            } break;
            case 1:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.c(),
                        tet.d(),
                        tet.b()
                    );
            } break;
            case 2:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.d(),
                        tet.b(),
                        tet.c()
                    );
            } break;
            case 3:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.b(),
                        tet.c(),
                        tet.d()
                    );
            } break;
            default:
            {
                FatalErrorInFunction
                    << "Point " << pointI << " is not present in tet" << tet
                    << abort(FatalError);
            }
        }
    }
}


// ************************************************************************* //
