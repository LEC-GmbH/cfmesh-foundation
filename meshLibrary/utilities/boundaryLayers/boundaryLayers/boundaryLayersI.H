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

Foam::label Foam::Module::boundaryLayers::findNewNodeLabel
(
    const label pointI,
    const label pKey
) const
{
    const Map<EdgeMap<label>>::const_iterator it =
        otherVrts_.find(pointI);

    if (it != otherVrts_.end())
    {
        const auto& m = it();

        if (m.size() == 2)
        {
            forAllConstIter(EdgeMap<label>, m, mit)
            {
                if (mit.key().start() != pKey)
                {
                    return mit();
                }
            }
        }
        else
        {
            forAllConstIter(EdgeMap<label>, m, mit)
            {
                if (mit.key().start() == pKey)
                {
                    continue;
                }
                if (mit.key().end() == pKey)
                {
                    continue;
                }
                if (mit.key().start() == mit.key().end())
                {
                    continue;
                }

                return mit();
            }
        }
    }

    return newLabelForVertex_[pointI];
}


inline void Foam::Module::boundaryLayers::createNewCellFromEdge
(
    const edge& e,
    const label pKeyI,
    const label pKeyJ,
    faceList& cellFaces
) const
{
    cellFaces.setSize(6);
    forAll(cellFaces, i)
    {
        cellFaces[i].setSize(4);
    }

    const Map<EdgeMap<label>>::const_iterator
        itS = otherVrts_.find(e.start());
    const Map<EdgeMap<label>>::const_iterator
        itE = otherVrts_.find(e.end());
    const EdgeMap<label>& ms = itS();
    const EdgeMap<label>& me = itE();

    label p0s(-1), p1s(-1), ns(-1), p0e(-1), p1e(-1), ne(-1);

    if (ms.size() == 2)
    {
        p0s = ms.find(edge(pKeyI, pKeyI))();
        p1s = ms.find(edge(pKeyJ, pKeyJ))();
        ns = newLabelForVertex_[e.start()];
    }
    else
    {
        forAllConstIter(EdgeMap<label>, ms, it)
        {
            if
            (
                (it.key().start() != pKeyJ) && (it.key().end() != pKeyJ)
             && (it.key().start() != it.key().end())
            )
            {
                p0s = it();
            }
            else if
            (
                (it.key().start() != pKeyI) && (it.key().end() != pKeyI)
             && (it.key().start() != it.key().end())
            )
            {
                p1s = it();
            }
            else if
            (
                (it.key().start() == it.key().end())
             && (it.key().start() != pKeyI) && (it.key().start() != pKeyJ)
            )
            {
                ns = it();
            }
        }
    }
    if (me.size() == 2)
    {
        p0e = me.find(edge(pKeyI, pKeyI))();
        p1e = me.find(edge(pKeyJ, pKeyJ))();
        ne = newLabelForVertex_[e.end()];
    }
    else
    {
        forAllConstIter(EdgeMap<label>, me, it)
        {
            if
            (
                (it.key().start() != pKeyJ) && (it.key().end() != pKeyJ)
             && (it.key().start() != it.key().end())
            )
            {
                p0e = it();
            }
            else if
            (
                (it.key().start() != pKeyI) && (it.key().end() != pKeyI)
             && (it.key().start() != it.key().end())
            )
            {
                p1e = it();
            }
            else if
            (
                (it.key().start() == it.key().end())
             && (it.key().start() != pKeyI) && (it.key().start() != pKeyJ)
            )
            {
                ne = it();
            }
        }
    }

    //- F0
    cellFaces[0][0] = e.end();
    cellFaces[0][1] = e.start();
    cellFaces[0][2] = p1s;
    cellFaces[0][3] = p1e;

    //- F1
    cellFaces[1][0] = p0e;
    cellFaces[1][1] = ne;
    cellFaces[1][2] = ns;
    cellFaces[1][3] = p0s;

    //- F2
    cellFaces[2][0] = e.start();
    cellFaces[2][1] = e.end();
    cellFaces[2][2] = p0e;
    cellFaces[2][3] = p0s;

    //- F3
    cellFaces[3][0] = p1s;
    cellFaces[3][1] = ns;
    cellFaces[3][2] = ne;
    cellFaces[3][3] = p1e;

    //- F4
    cellFaces[4][0] = e.end();
    cellFaces[4][1] = p1e;
    cellFaces[4][2] = ne;
    cellFaces[4][3] = p0e;

    //- F5
    cellFaces[5][0] = e.start();
    cellFaces[5][1] = p0s;
    cellFaces[5][2] = ns;
    cellFaces[5][3] = p1s;
}


inline void Foam::Module::boundaryLayers::createNewCellFromNode
(
    const label pointI,
    const DynamicList<label>& pKeys,
    faceList& cellFaces
) const
{
    cellFaces.setSize(6);
    forAll(cellFaces, i)
    {
        cellFaces[i].setSize(4);
    }

    const Map<EdgeMap<label>>::const_iterator
        itP = otherVrts_.find(pointI);
    const EdgeMap<label>& m = itP();

    //- create labels before creating cells
    const label n = newLabelForVertex_[pointI];
    const label p00 =
        m.find(edge(pKeys[0], pKeys[0]))();
    const label p11 =
        m.find(edge(pKeys[1], pKeys[1]))();
    const label p22 =
        m.find(edge(pKeys[2], pKeys[2]))();

    edge pr;
    pr.start() = pKeys[0];
    pr.end() = pKeys[1];
    if (!m.found(pr))
    {
        pr.start() = pKeys[1];
        pr.end() = pKeys[0];
    }
    const label p01 = m.find(pr)();

    pr.start() = pKeys[0];
    pr.end() = pKeys[2];
    if (!m.found(pr))
    {
        pr.start() = pKeys[2];
        pr.end() = pKeys[0];
    }
    const label p02 = m.find(pr)();

    pr.start() = pKeys[1];
    pr.end() = pKeys[2];
    if (!m.found(pr))
    {
        pr.start() = pKeys[2];
        pr.end() = pKeys[1];
    }
    const label p12 = m.find(pr)();

    //- create the cell and append it
    //- F0
    cellFaces[0][0] = pointI;
    cellFaces[0][1] = p02;
    cellFaces[0][2] = p00;
    cellFaces[0][3] = p01;
    //- F1
    cellFaces[1][0] = p12;
    cellFaces[1][1] = p11;
    cellFaces[1][2] = n;
    cellFaces[1][3] = p22;

    //- F2
    cellFaces[2][0] = pointI;
    cellFaces[2][1] = p01;
    cellFaces[2][2] = p11;
    cellFaces[2][3] = p12;
    //- F3
    cellFaces[3][0] = p02;
    cellFaces[3][1] = p22;
    cellFaces[3][2] = n;
    cellFaces[3][3] = p00;
    //- F4
    cellFaces[4][0] = pointI;
    cellFaces[4][1] = p12;
    cellFaces[4][2] = p22;
    cellFaces[4][3] = p02;
    //- F5
    cellFaces[5][0] = p01;
    cellFaces[5][1] = p00;
    cellFaces[5][2] = n;
    cellFaces[5][3] = p11;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
