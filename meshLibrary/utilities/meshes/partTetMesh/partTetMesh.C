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

#include "demandDrivenData.H"
#include "partTetMesh.H"
#include "polyMeshGenModifier.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"

#include <map>

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Module::partTetMesh::partTetMesh
(
    polyMeshGen& mesh,
    const DynamicList<label>& lockedPoints
)
:
    origMesh_(mesh),
    points_(),
    tets_(),
    nodeLabelInOrigMesh_(),
    smoothVertex_(),
    pointTets_(),
    internalPointsOrderPtr_(nullptr),
    boundaryPointsOrderPtr_(nullptr)
{
    List<direction> useCell(mesh.cells().size(), direction(1));

    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
    {
        lockedPoint[lockedPoints[i]] = true;
    }

    createPointsAndTets(useCell, lockedPoint);
}


Foam::Module::partTetMesh::partTetMesh
(
    polyMeshGen& mesh,
    const DynamicList<label>& lockedPoints,
    const direction nLayers
)
:
    origMesh_(mesh),
    points_(),
    tets_(),
    nodeLabelInOrigMesh_(),
    smoothVertex_(),
    pointTets_(),
    internalPointsOrderPtr_(nullptr),
    boundaryPointsOrderPtr_(nullptr)
{
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const labelList& owner = mesh.owner();
    const PtrList<boundaryPatch>& boundaries = mesh.boundaries();
    const labelListList2D& pointCells = mesh.addressingData().pointCells();

    List<direction> useCell(cells.size(), direction(0));

    // select cells containing at least one vertex of the bad faces
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label size = boundaries[patchI].patchSize();

        for (label fI = 0; fI < size; ++fI)
        {
            useCell[owner[start + fI]] = 1;
        }
    }

    // add additional layer of cells
    for (direction layerI = 1; layerI<(nLayers + 1); ++layerI)
    {
        forAll(useCell, cI)
        {
            if (useCell[cI] == layerI)
            {
                const cell& c = cells[cI];

                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, pI)
                    {
                        forAllRow(pointCells, f[pI], pcI)
                        {
                            const label cLabel = pointCells[f[pI]][pcI];
                            if (!useCell[cLabel])
                            {
                                useCell[cLabel] = layerI + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
    {
        lockedPoint[lockedPoints[i]] = true;
    }

    createPointsAndTets(useCell, lockedPoint);
}


Foam::Module::partTetMesh::partTetMesh
(
    polyMeshGen& mesh,
    const DynamicList<label>& lockedPoints,
    labelHashSet& badFaces,
    const direction additionalLayers
)
:
    origMesh_(mesh),
    points_(),
    tets_(),
    nodeLabelInOrigMesh_(),
    smoothVertex_(),
    pointTets_(),
    internalPointsOrderPtr_(nullptr),
    boundaryPointsOrderPtr_(nullptr)
{
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const labelListList2D& pointCells = mesh.addressingData().pointCells();

    List<direction> useCell(cells.size(), direction(0));

    // select cells containing at least one vertex of the bad faces
    forAll(faces, faceI)
    {
        if (badFaces.found(faceI))
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                forAllRow(pointCells, f[pI], pcI)
                {
                    useCell[pointCells[f[pI]][pcI]] = 1;
                }
            }
        }
    }

    // add additional layer of cells
    for (direction layerI = 1; layerI<(additionalLayers + 1); ++layerI)
    {
        forAll(useCell, cI)
        {
            if (useCell[cI] == layerI)
            {
                const cell& c = cells[cI];

                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, pI)
                    {
                        forAllRow(pointCells, f[pI], pcI)
                        {
                            const label cLabel = pointCells[f[pI]][pcI];
                            if (!useCell[cLabel])
                            {
                                useCell[cLabel] = layerI + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
    {
        lockedPoint[lockedPoints[i]] = true;
    }

    createPointsAndTets(useCell, lockedPoint);
}


Foam::Module::partTetMesh::~partTetMesh()
{
    deleteDemandDrivenData(internalPointsOrderPtr_);
    deleteDemandDrivenData(boundaryPointsOrderPtr_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Module::labelListList2D&
Foam::Module::partTetMesh::internalPointOrdering() const
{
    if (omp_in_parallel())
    {
        FatalErrorInFunction
            << "Calculating addressing inside a parallel region."
            << " This is not thread safe" << exit(FatalError);
    }

    if (!internalPointsOrderPtr_)
    {
        createSMOOTHPointsOrdering();
    }

    return *internalPointsOrderPtr_;
}


const Foam::Module::labelListList2D&
Foam::Module::partTetMesh::boundaryPointOrdering() const
{
    if (omp_in_parallel())
    {
        FatalErrorInFunction
            << "Calculating addressing inside a parallel region."
            << " This is not thread safe" << exit(FatalError);
    }

    if (!boundaryPointsOrderPtr_)
    {
        createBOUNDARYPointsOrdering();
    }

    return *boundaryPointsOrderPtr_;
}


void Foam::Module::partTetMesh::updateVertex
(
    const label pointI,
    const point& newP
)
{
    points_[pointI] = newP;

    if (smoothVertex_[pointI] & (FACECENTRE + CELLCENTRE))
    {
        Warning << "Smoothing auxiliary vertex."
            << " This has no effect on the original mesh" << endl;
        return;
    }

    // find face centres attached
    DynamicList<label> helper;
    helper.reserve(64);
    forAllRow(pointTets_, pointI, ptI)
    {
        const label centreI = tets_[pointTets_[pointI][ptI]][2];
        if (smoothVertex_[centreI] & FACECENTRE)
        {
            appendUniq(helper, centreI);
        }
    }

    // update coordinates of FACECENTRE vertices
    forAll(helper, i)
    {
        const label centreI = helper[i];

        point centre(vector::zero);
        scalar faceArea(0.0);
        forAllRow(pointTets_, centreI, ptI)
        {
            const partTet& tet = tets_[pointTets_[centreI][ptI]];
            point c(vector::zero);
            for (label i = 0; i < 3; ++i)
            {
                c += points_[tet[i]];
            }
            c /= 3;
            const scalar area = Foam::mag(tet.Sd(points_)) + VSMALL;

            centre += c*area;
            faceArea += area;
        }

        points_[centreI] = centre/faceArea;
    }

    // find cell centres attached
    helper.clear();
    forAllRow(pointTets_, pointI, ptI)
    {
        const label centreI = tets_[pointTets_[pointI][ptI]][3];
        if (smoothVertex_[centreI] & CELLCENTRE)
        {
            appendUniq(helper, centreI);
        }
    }

    // update coordinates of CELLCENTRE vertices
    forAll(helper, i)
    {
        const label centreI = helper[i];

        point centre(vector::zero);
        scalar cellVol(0.0);
        forAllRow(pointTets_, centreI, ptI)
        {
            const partTet& tet = tets_[pointTets_[centreI][ptI]];
            const point c = tet.centroid(points_);
            const scalar vol = Foam::mag(tet.mag(points_)) + VSMALL;

            centre += c*vol;
            cellVol += vol;
        }

        points_[centreI] = centre/cellVol;
    }
}


void Foam::Module::partTetMesh::updateVerticesSMP
(
    const List<DynamicList<labelledPoint>>& np
)
{
    List<direction> updateType(points_.size(), direction(0));

    # pragma omp parallel num_threads(np.size())
    {
        const DynamicList<labelledPoint>& newPoints = np[omp_get_thread_num()];

        forAll(newPoints, i)
        {
            const labelledPoint& lp = newPoints[i];
            const label pointI = lp.pointLabel();

            points_[pointI] = lp.coordinates();

            forAllRow(pointTets_, pointI, ptI)
            {
                const partTet& pt = tets_[pointTets_[pointI][ptI]];

                if (smoothVertex_[pt[3]] & CELLCENTRE)
                {
                    # pragma omp atomic update
                    updateType[pt[3]] |= CELLCENTRE;
                }
                if (smoothVertex_[pt[2]] & FACECENTRE)
                {
                    # pragma omp atomic update
                    updateType[pt[2]] |= FACECENTRE;
                }
            }
        }

        # pragma omp barrier

        # pragma omp flush(updateType)

        // update coordinates of CELLCENTRE and FACECENTRE vertices
        # pragma omp for schedule(dynamic, 20)
        forAll(updateType, pI)
        {
            if (updateType[pI] & CELLCENTRE)
            {
                point centre(vector::zero);
                scalar cellVol(0.0);
                forAllRow(pointTets_, pI, ptI)
                {
                    const partTet& tet = tets_[pointTets_[pI][ptI]];
                    const point c = tet.centroid(points_);
                    const scalar vol = Foam::mag(tet.mag(points_)) + VSMALL;

                    centre += c*vol;
                    cellVol += vol;
                }

                points_[pI] = centre/cellVol;
            }
            else if (updateType[pI] & FACECENTRE)
            {
                point centre(vector::zero);
                scalar faceArea(0.0);
                forAllRow(pointTets_, pI, ptI)
                {
                    const partTet& tet = tets_[pointTets_[pI][ptI]];
                    point c(vector::zero);
                    for (label i = 0; i < 3; ++i)
                    {
                        c += points_[tet[i]];
                    }
                    c /= 3;
                    const scalar area = Foam::mag(tet.Sd(points_)) + VSMALL;

                    centre += c*area;
                    faceArea += area;
                }

                points_[pI] = centre/faceArea;
            }
        }
    }
}


void Foam::Module::partTetMesh::updateOrigMesh(boolList* changedFacePtr)
{
    pointFieldPMG& pts = origMesh_.points();

    boolList changedNode(pts.size(), false);

    # pragma omp parallel for if (pts.size() > 1000) \
    schedule(guided, 10)
    forAll(nodeLabelInOrigMesh_, pI)
    {
        if (nodeLabelInOrigMesh_[pI] != -1)
        {
            changedNode[nodeLabelInOrigMesh_[pI]] = true;
            pts[nodeLabelInOrigMesh_[pI]] = points_[pI];
        }
    }

    if (changedFacePtr)
    {
        boolList& chF = *changedFacePtr;
        chF = false;

        const cellListPMG& cells = origMesh_.cells();
        const labelListList2D& pointCells = origMesh_.addressingData().pointCells();

        # pragma omp parallel for if (pointCells.size() > 100) \
        schedule(dynamic, 20)
        forAll(pointCells, pointI)
        {
            if (changedNode[pointI])
            {
                forAllRow(pointCells, pointI, pcI)
                {
                    const cell& c = cells[pointCells[pointI][pcI]];

                    forAll(c, fI)
                    {
                        chF[c[fI]] = true;
                    }
                }
            }
        }

        // update geometry information
        const_cast<polyMeshGenAddressing&>
        (
            origMesh_.addressingData()
        ).updateGeometry(chF);
    }
    else
    {
        const_cast<polyMeshGenAddressing&>
        (
            origMesh_.addressingData()
        ).clearGeom();
    }
}


void Foam::Module::partTetMesh::createPolyMesh(polyMeshGen& pmg) const
{
    if (omp_in_parallel())
    {
        FatalErrorInFunction
            << "Building cellFaces inside a parallel region."
            << " This is not thread safe" << exit(FatalError);
    }

    polyMeshGenModifier meshModifier(pmg);

    pointFieldPMG& pAccess = meshModifier.pointsAccess();
    pAccess.setSize(points_.size());
    forAll(points_, pI)
    {
        pAccess[pI] = points_[pI];
    }

    DynamicList<faceList> cellFaces;

    forAll(tets_, tetI)
    {
        const partTet& tet = tets_[tetI];

        faceList tetFaces(4);

        // face 0
        tetFaces[0].setSize(3);
        tetFaces[0][0] = tet[0];
        tetFaces[0][1] = tet[2];
        tetFaces[0][2] = tet[1];

        // face 1
        tetFaces[1].setSize(3);
        tetFaces[1][0] = tet[0];
        tetFaces[1][1] = tet[1];
        tetFaces[1][2] = tet[3];

        // face 2
        tetFaces[2].setSize(3);
        tetFaces[2][0] = tet[0];
        tetFaces[2][1] = tet[3];
        tetFaces[2][2] = tet[2];

        // face 3
        tetFaces[3].setSize(3);
        tetFaces[3][0] = tet[1];
        tetFaces[3][1] = tet[2];
        tetFaces[3][2] = tet[3];

        cellFaces.append(tetFaces);
    }

    meshModifier.addCells(cellFaces);
    meshModifier.reorderBoundaryFaces();

    // store points into subsets
    const label bndPointID = pmg.addPointSubset("boundaryPoints");
    const label smoothPointID = pmg.addPointSubset("smoothPoints");
    const label faceCentreID = pmg.addPointSubset("faceCentres");
    const label cellCentreID = pmg.addPointSubset("cellCentres");

    forAll(smoothVertex_, pointI)
    {
        if (smoothVertex_[pointI] & SMOOTH)
        {
            pmg.addPointToSubset(smoothPointID, pointI);
        }
        if (smoothVertex_[pointI] & BOUNDARY)
        {
            pmg.addPointToSubset(bndPointID, pointI);
        }
        if (smoothVertex_[pointI] & FACECENTRE)
        {
            pmg.addPointToSubset(faceCentreID, pointI);
        }
        if (smoothVertex_[pointI] & CELLCENTRE)
        {
            pmg.addPointToSubset(cellCentreID, pointI);
        }
    }

    const labelListList2D& internalOrdering = internalPointOrdering();
    const labelListList2D& boundaryOrdering = boundaryPointOrdering();

    forAll(internalOrdering, i)
    {
        const word name = "smoothPoints_" + Foam::name(i);
        const label orderID = pmg.addPointSubset(name);

        forAllRow(internalOrdering, i, nI)
        {
            pmg.addPointToSubset(orderID, internalOrdering[i][nI]);
        }
    }

    forAll(boundaryOrdering, i)
    {
        const word name = "boundaryPoints_" + Foam::name(i);
        const label orderID = pmg.addPointSubset(name);

        forAllRow(boundaryOrdering, i, nI)
        {
            pmg.addPointToSubset(orderID, boundaryOrdering[i][nI]);
        }
    }
}


// ************************************************************************* //
