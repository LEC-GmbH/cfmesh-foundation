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

#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::refineBoundaryLayers::generateNewCellsPrism
(
    const label cellI,
    List<DynamicList<DynamicList<label>>>& cellsFromCell
) const
{
    cellsFromCell.clear();

    const cell& c = mesh_.cells()[cellI];
    const labelList& owner = mesh_.owner();

    const label startBoundary = mesh_.boundaries()[0].patchStart();

    // find the number of lyers for this cell
    label nLayers(1), baseFace(-1);
    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if ((bfI < 0) || (bfI >= nLayersAtBndFace_.size()))
        {
            continue;
        }

        if (nLayersAtBndFace_[bfI] < 2)
        {
            continue;
        }

        nLayers = nLayersAtBndFace_[bfI];
        baseFace = fI;
    }

    // set the number of layers
    cellsFromCell.setSize(nLayers);

    // distribute existing faces into new cells
    label otherBaseFace(-1);
    forAll(c, fI)
    {
        if (fI == baseFace)
        {
            const label faceI = facesFromFace_[c[fI]][0];
            DynamicList<label> f;
            f = newFaces_[faceI];
            cellsFromCell[nLayers - 1].append(f);
        }
        else if (facesFromFace_[c[fI]].size() == 1)
        {
            const label faceI = facesFromFace_[c[fI]][0];
            otherBaseFace = fI;
            DynamicList<label> f;
            f = newFaces_[faceI];
            cellsFromCell[0].append(f);
        }
        else
        {
            forAllRow(facesFromFace_, c[fI], cfI)
            {
                const label nfI = facesFromFace_[c[fI]][cfI];

                DynamicList<label> cf;
                cf = newFaces_[nfI];

                if (owner[c[fI]] != cellI)
                {
                    cf = help::reverseFace(cf);
                }

                cellsFromCell[Foam::max(nLayers - 1 - cfI, 0)].append(cf);
            }
        }
    }

    // generate missing faces
    const faceListPMG& faces = mesh_.faces();
    const face& bf = faces[c[baseFace]];
    const face& obf = faces[c[otherBaseFace]];
    for (label layerI = 1; layerI < nLayers; ++layerI)
    {
        // create new face from points at the same height
        DynamicList<label> cf;
        cf.reserve(8);
        forAll(bf, pI)
        {
            const label pointI = bf[pI];

            label seI(-1);
            if (splitEdgesAtPoint_[pointI].size() == 1)
            {
                seI = splitEdgesAtPoint_[pointI][0];
            }
            else
            {
                forAllRow(splitEdgesAtPoint_, pointI, sepI)
                {
                    const label seJ = splitEdgesAtPoint_[pointI][sepI];
                    const edge& se = splitEdges_[seJ];

                    if (obf.which(se.end()) >= 0 || obf.which(se.start()) >= 0)
                    {
                        seI = seJ;
                        break;
                    }
                }
            }

            cf.append(vertexOnSplitEdge(seI, layerI));
        }

        // add faces to cells
        cellsFromCell[nLayers - layerI].append(cf);
        cellsFromCell[nLayers - 1 - layerI].append(cf);
    }
}


void Foam::Module::refineBoundaryLayers::storeFacesIntoCells
(
    const label faceI,
    const bool reverseOrientation,
    const label normalDirection,
    const bool maxCoordinate,
    const label nLayersI,
    const label nLayersJ,
    const label nLayersK,
    List<DynamicList<DynamicList<label>>>& cellsFromCell
) const
{
    labelListList2D faceFaces;
    sortFaceFaces(faceI, faceFaces, reverseOrientation);

    const label maxI = nLayersI - 1;
    const label maxJ = nLayersJ - 1;
    const label maxK = nLayersK - 1;

    label i(-1), j(-1), k(-1);

    forAll(faceFaces, nI)
    {
        forAll(faceFaces[nI], nJ)
        {
            const label nfI = faceFaces[nI][nJ];

            if (normalDirection == 0)
            {
                // k is const
                i = Foam::min(nI, maxI);
                j = Foam::min(nJ, maxJ);
                k = maxCoordinate ? maxK : 0;
            }
            else if (normalDirection == 1)
            {
                // j is const
                i = Foam::min(nJ, maxI);
                j = maxCoordinate ? maxJ : 0;
                k = Foam::min(nI, maxK);
            }
            else if (normalDirection == 2)
            {
                // i is const
                i = maxCoordinate ? maxI : 0;
                j = Foam::min(nI, maxJ);
                k = Foam::min(nJ, maxK);
            }

            // store the face into a new cell
            const label cI
            (
                j*nLayersI +
                k*nLayersI*nLayersJ +
                i
            );

            cellsFromCell[cI].append(newFaces_[nfI]);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::refineBoundaryLayers::refineEdgeHexCell::
determineFacesInDirections()
{
    const labelList& nLayersAtBndFace = bndLayers_.nLayersAtBndFace_;
    const polyMeshGen& mesh = bndLayers_.mesh_;
    const faceListPMG& faces = mesh.faces();
    const cell& c = mesh.cells()[cellI_];

    const PtrList<boundaryPatch>& bnd = mesh.boundaries();
    const label startBoundary = bnd[0].patchStart();

    // find the number of layers for this cell
    FixedList<label, 2> layersInDirection(-1), dirFace;
    label currDir(0);

    FixedList<bool, 6> determinedFace(false);

    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if ((bfI < 0) || (bfI >= nLayersAtBndFace.size()))
        {
            continue;
        }

        if (nLayersAtBndFace[bfI] < 2)
        {
            continue;
        }

        layersInDirection[currDir] = nLayersAtBndFace[bfI];
        dirFace[currDir] = fI;
        ++currDir;
    }

    // set the number of newly create cells
    nLayersI_ = layersInDirection[0];
    nLayersJ_ = layersInDirection[1];
    cellsFromCell_.setSize(nLayersI_*nLayersJ_);

    // find the shared edge between the boundary faces
    const edge commonEdge =
        help::sharedEdge(faces[c[dirFace[0]]], faces[c[dirFace[1]]]);

    // faces at i = const in the local coordinate system
    faceInDirection_[4] = dirFace[0];
    determinedFace[dirFace[0]] = true;
    forAll(c, fI)
    {
        if (determinedFace[fI])
        {
            continue;
        }

        if (!help::shareAnEdge(faces[c[dirFace[0]]], faces[c[fI]]))
        {
            faceInDirection_[5] = fI;
            determinedFace[fI] = true;
            break;
        }
    }

    // faces k = const in the local coordinate system
    faceInDirection_[2] = dirFace[1];
    determinedFace[dirFace[1]] = true;
    forAll(c, fI)
    {
        if (determinedFace[fI])
        {
            continue;
        }

        if (!help::shareAnEdge(faces[c[dirFace[1]]], faces[c[fI]]))
        {
            faceInDirection_[3] = fI;
            determinedFace[fI] = true;
            break;
        }
    }

    // find the face attached to the starting point of the edge and
    // the face attached to the end point of the edge
    forAll(c, fI)
    {
        if (determinedFace[fI])
        {
            continue;
        }

        if
        (
            (faces[c[fI]].which(commonEdge.start()) >= 0) &&
            (help::positionOfEdgeInFace(commonEdge, faces[c[fI]]) < 0)
        )
        {
            faceInDirection_[0] = fI;
        }

        if
        (
            (faces[c[fI]].which(commonEdge.end()) >= 0) &&
            (help::positionOfEdgeInFace(commonEdge, faces[c[fI]]) < 0)
        )
        {
            faceInDirection_[1] = fI;
        }
    }

    // check the orientation of faces
    const labelList& owner = mesh.owner();

    // checking face at direction k = 0
    faceOrientation_[0] = owner[c[faceInDirection_[0]]] == cellI_?true:false;

    // checking face in direction k = 1
    faceOrientation_[1] = owner[c[faceInDirection_[1]]] == cellI_?false:true;

    // set orientation flag for face in direction j = 0
    faceOrientation_[2] = true;

    // checking face in direction j = nLayersJ_
    faceOrientation_[3] = owner[c[faceInDirection_[3]]] == cellI_?false:true;

    // set orientation flag for face in direction i = 0
    faceOrientation_[4] = true;

    // checking face in direction i = nLayersI_
    faceOrientation_[5] = owner[c[faceInDirection_[5]]] == cellI_?false:true;

}


void Foam::Module::refineBoundaryLayers::refineEdgeHexCell::
populateExistingFaces()
{
    const cell& c = bndLayers_.mesh_.cells()[cellI_];
    const labelListList2D& facesFromFace = bndLayers_.facesFromFace_;
    const DynamicList<DynamicList<label>>& newFaces = bndLayers_.newFaces_;

    cellsFromCell_.setSize(nLayersI_*nLayersJ_);
    forAll(cellsFromCell_, cI)
    {
        cellsFromCell_[cI].clear();
    }

    // store new faces at k = 0
    bndLayers_.storeFacesIntoCells
    (
        c[faceInDirection_[0]], faceOrientation_[0],
        0, 0,
        nLayersI_, nLayersJ_, 1,
        cellsFromCell_
    );

    // store new faces at k = 1
    bndLayers_.storeFacesIntoCells
    (
        c[faceInDirection_[1]], faceOrientation_[1],
        0, 1,
        nLayersI_, nLayersJ_, 1,
        cellsFromCell_
    );

    // store new faces at j = 0
    forAllRow(facesFromFace, c[faceInDirection_[2]], i)
    {
        const label faceI = facesFromFace[c[faceInDirection_[2]]][i];
        cellsFromCell_[i].append(newFaces[faceI]);
    }

    // store faces at j = nLayersJ
    const label maxJ = nLayersJ_ - 1;
    forAllRow(facesFromFace, c[faceInDirection_[3]], i)
    {
        const label faceI = facesFromFace[c[faceInDirection_[3]]][i];
        cellsFromCell_[i + maxJ*nLayersI_].append(newFaces[faceI]);
    }

    // store new faces at i = 0
    forAllRow(facesFromFace, c[faceInDirection_[4]], j)
    {
        const label faceI = facesFromFace[c[faceInDirection_[4]]][j];
        cellsFromCell_[j*nLayersI_].append(newFaces[faceI]);
    }

    // store new faces at i = nLayersI
    const label maxI = nLayersI_ - 1;
    forAllRow(facesFromFace, c[faceInDirection_[5]], j)
    {
        const label faceI = facesFromFace[c[faceInDirection_[5]]][j];
        cellsFromCell_[j*nLayersI_ + maxI].append(newFaces[faceI]);
    }
}


void Foam::Module::refineBoundaryLayers::refineEdgeHexCell::
generateMissingFaces()
{
    const cell& c = bndLayers_.mesh_.cells()[cellI_];

    // fill up the matrix of points for this cell
    // the matrix is used for generation of new cells
    FixedList<labelListList2D, 2> cellPoints;

    // fill in the data for a cross-split faces
    bndLayers_.sortFacePoints
    (
        c[faceInDirection_[0]],
        cellPoints[0],
        faceOrientation_[0]
    );
    bndLayers_.sortFacePoints
    (
        c[faceInDirection_[1]],
        cellPoints[1],
        faceOrientation_[1]
    );

    // generate new internal faces for this cell
    // generate faces with normal in the i direction
    const label maxI = nLayersI_ - 1;
    const label maxJ = nLayersJ_ - 1;

    for (label i = 1; i < nLayersI_; ++i)
    {
        for (label j = 0; j < nLayersJ_; ++j)
        {
            const label own = j*nLayersI_ + i - 1;
            const label nei = own + 1;

            if (j < maxJ)
            {
                // generate a quad face
                DynamicList<label> mf;
                mf.setSize(4);

                // populate the points form cellPoints
                mf[0] = cellPoints[0][i][j];
                mf[1] = cellPoints[0][i][j + 1];
                mf[2] = cellPoints[1][i][j + 1];
                mf[3] = cellPoints[1][i][j];

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
            else
            {
                DynamicList<label> mf;
                for (label index = j; index < cellPoints[0][i].size(); ++index)
                {
                    mf.append(cellPoints[0][i][index]);
                }
                for (label index = cellPoints[1][i].size()-1; index>=j; --index)
                {
                    mf.append(cellPoints[1][i][index]);
                }

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            };
        }
    }

    // generate faces with the normal in j direction
    for (label i = 0; i < nLayersI_; ++i)
    {
        for (label j = 1; j < nLayersJ_; ++j)
        {
            const label nei = j*nLayersI_ + i;
            const label own = (j - 1)*nLayersI_ + i;

            if (i < maxI)
            {
                // generate a quad face
                DynamicList<label> mf;
                mf.setSize(4);

                // populate the points form cellPoints
                mf[0] = cellPoints[0][i][j];
                mf[1] = cellPoints[1][i][j];
                mf[2] = cellPoints[1][i + 1][j];
                mf[3] = cellPoints[0][i + 1][j];

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
            else
            {
                DynamicList<label> mf;
                for (label index = i; index < cellPoints[1].size(); ++index)
                {
                    mf.append(cellPoints[1][index][j]);
                }
                for (label index = cellPoints[0].size()-1; index>=i; --index)
                {
                    mf.append(cellPoints[0][index][j]);
                }

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            };
        }
    }
}


Foam::Module::refineBoundaryLayers::refineEdgeHexCell::refineEdgeHexCell
(
    const label cellI,
    const refineBoundaryLayers& ref
)
:
    cellI_(cellI),
    nLayersI_(),
    nLayersJ_(),
    bndLayers_(ref)
{
    determineFacesInDirections();

    populateExistingFaces();

    generateMissingFaces();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::refineBoundaryLayers::refineCornerHexCell::
determineFacesInDirections()
{
    const polyMeshGen& mesh = bndLayers_.mesh_;
    const cell& c = mesh.cells()[cellI_];
    const faceListPMG& faces = mesh.faces();
    const labelList& nLayersAtBndFace = bndLayers_.nLayersAtBndFace_;

    const label startBoundary = mesh.boundaries()[0].patchStart();

    // find the number of layers for this cell
    FixedList<label, 3> layersInDirection(-1), dirFace;
    FixedList<bool, 6> usedDirection(false);
    label currDir(0);

    forAll(c, fI)
    {
        const label bfI = c[fI] - startBoundary;

        if ((bfI < 0) || (bfI >= nLayersAtBndFace.size()))
        {
            continue;
        }

        if (nLayersAtBndFace[bfI] < 2)
        {
            continue;
        }

        usedDirection[fI] = true;
        layersInDirection[currDir] = nLayersAtBndFace[bfI];
        dirFace[currDir] = fI;
        ++currDir;
    }

    // find a common point for all three boundary faces
    FixedList<DynamicList<label>, 3> bndFaces;
    forAll(dirFace, i)
    {
        bndFaces[i] = faces[c[dirFace[i]]];
    }

    const label commonPoint = help::sharedVertex(bndFaces);

    // find the position of the common point in each boundary face
    const DynamicList<edge>& splitEdges = bndLayers_.splitEdges_;
    const labelListList2D& splitEdgesAtPoint = bndLayers_.splitEdgesAtPoint_;

    const face& baseFace = faces[c[dirFace[0]]];
    const label posInBndFace = baseFace.which(commonPoint);

    // find split edges starting at the commonPoints
    forAllRow(splitEdgesAtPoint, commonPoint, i)
    {
        const edge& se = splitEdges[splitEdgesAtPoint[commonPoint][i]];

        if (se == baseFace.faceEdge(posInBndFace))
        {
            // this edge is in j direction
            splitEdgeInDirection_[1] = splitEdgesAtPoint[commonPoint][i];
        }
        else if (se == baseFace.faceEdge(baseFace.rcIndex(posInBndFace)))
        {
            // this edge is in i diretion
            splitEdgeInDirection_[0] = splitEdgesAtPoint[commonPoint][i];
        }
        else if (splitEdgesAtPoint[commonPoint].size() == 3)
        {
            // this point is in k direction
            splitEdgeInDirection_[2] = splitEdgesAtPoint[commonPoint][i];
        }
        else
        {
            // this situation is not allowed
            FatalErrorInFunction
                << "Cannot refine layer for cell " << cellI_
                << abort(FatalError);
        }
    }

    // find the direction of other boundary faces
    // in the local coordinate system
    FixedList<label, 3> permutation;
    permutation[0] = 0;

    label helper = help::positionOfEdgeInFace
    (
        baseFace.faceEdge(baseFace.rcIndex(posInBndFace)),
        faces[c[dirFace[1]]]
    );

    if (helper >= 0)
    {
        permutation[1] = 1;
        permutation[2] = 2;
    }
    else
    {
        permutation[1] = 2;
        permutation[2] = 1;
    }

    // find the number of layers and a split in each direction
    nLayersI_ = layersInDirection[permutation[2]];
    nLayersJ_ = layersInDirection[permutation[1]];
    nLayersK_ = layersInDirection[permutation[0]];

    // determine the directions of cell faces
    // store boundary faces first. Their normals point in the wrong direction
    // face at k = 0
    faceInDirection_[0] = dirFace[permutation[0]];
    faceOrientation_[0] = true;
    // face at j = 0
    faceInDirection_[2] = dirFace[permutation[1]];
    faceOrientation_[2] = true;
    // face at i = 0
    faceInDirection_[4] = dirFace[permutation[2]];
    faceOrientation_[4] = true;

    // find directions of other faces and their orientation
    const labelList& owner = mesh.owner();
    forAll(c, fI)
    {
        if (usedDirection[fI])
        {
            continue;
        }

        const bool orientation = owner[c[fI]] == cellI_ ? false : true;

        if (!help::shareAnEdge(faces[c[fI]], faces[c[faceInDirection_[0]]]))
        {
            // face at k = nLayersK_
            faceInDirection_[1] = fI;
            faceOrientation_[1] = orientation;
        }
        else if
        (
            !help::shareAnEdge(faces[c[fI]], faces[c[faceInDirection_[2]]])
        )
        {
            // face at j = nLayersJ_
            faceInDirection_[3] = fI;
            faceOrientation_[3] = orientation;
        }
        else if
        (
            !help::shareAnEdge(faces[c[fI]], faces[c[faceInDirection_[4]]])
        )
        {
            // face at i = nLayersI_
            faceInDirection_[5] = fI;
            faceOrientation_[5] = orientation;
        }
    }
}


void Foam::Module::refineBoundaryLayers::refineCornerHexCell::
populateExistingFaces()
{
    const cell& c = bndLayers_.mesh_.cells()[cellI_];

    // set the number of cells
    cellsFromCell_.setSize(nLayersI_*nLayersJ_*nLayersK_);
    forAll(cellsFromCell_, i)
    {
        cellsFromCell_[i].clear();
    }

    // add new faces from existing faces into new cells
    forAll(faceInDirection_, dirI)
    {
        bndLayers_.storeFacesIntoCells
        (
            c[faceInDirection_[dirI]], faceOrientation_[dirI],
            dirI/2, dirI % 2,
            nLayersI_, nLayersJ_, nLayersK_,
            cellsFromCell_
        );
    }

}


void Foam::Module::refineBoundaryLayers::refineCornerHexCell::
generateNewPoints()
{
    const cell& c = bndLayers_.mesh_.cells()[cellI_];

    // allocate space for points generated inside the cell
    cellPoints_.setSize(nLayersI_ + 1);
    forAll(cellPoints_, i)
    {
        cellPoints_[i].setSize(nLayersJ_ + 1);

        forAll(cellPoints_[i], j)
        {
            cellPoints_[i][j].setSize(nLayersK_ + 1);
            cellPoints_[i][j] = -1;
        }
    }

    // collect information about points generated on faces of the cell
    forAll(faceInDirection_, dirI)
    {
        bndLayers_.sortFacePoints
        (
            c[faceInDirection_[dirI]],
            facePoints_[dirI],
            faceOrientation_[dirI]
        );
    }

    // fill in cellPoints at the boundary
    forAll(cellPoints_, i)
    {
        forAll(cellPoints_[i], j)
        {
            cellPoints_[i][j][0] = facePoints_[0][i][j];
            cellPoints_[i][j][nLayersK_] = facePoints_[1][i][j];
        }
    }

    forAll(cellPoints_, i)
    {
        forAll(cellPoints_[i][0], k)
        {
            cellPoints_[i][0][k] = facePoints_[2][k][i];
            cellPoints_[i][nLayersJ_][k] = facePoints_[3][k][i];
        }
    }

    forAll(cellPoints_[0], j)
    {
        forAll(cellPoints_[0][j], k)
        {
            cellPoints_[0][j][k] = facePoints_[4][j][k];
            cellPoints_[nLayersI_][j][k] = facePoints_[5][j][k];
        }
    }

    // useful data for generating missing points
    const DynamicList<edge>& splitEdges = bndLayers_.splitEdges_;
    const edge& seDirI = splitEdges[splitEdgeInDirection_[0]];
    const edge& seDirJ = splitEdges[splitEdgeInDirection_[1]];
    const edge& seDirK = splitEdges[splitEdgeInDirection_[2]];

    // const references to vertices of the cell ordered in a local
    // i, j, k coordinate system
    pointFieldPMG& points = bndLayers_.mesh_.points();
    const point v000 = points[seDirI.start()];
    const point v100 = points[seDirI.end()];
    const point v110 = points[facePoints_[0].last().last()];
    const point v010 = points[seDirJ.end()];
    const point v001 = points[seDirK.end()];
    const point v101 = points[facePoints_[1].last()[0]];
    const point v111 = points[facePoints_[1].last().last()];
    const point v011 = points[facePoints_[1][0].last()];

    for (label i = 1; i < nLayersI_; ++i)
    {
        const scalar u
        (
            Foam::mag
            (
                points[bndLayers_.vertexOnSplitEdge(splitEdgeInDirection_[0], i)] -
                points[seDirI.start()]
            ) /
            seDirI.mag(points)
        );

        for (label j = 1; j < nLayersJ_; ++j)
        {
            const scalar v
            (
                Foam::mag
                (
                    points[bndLayers_.vertexOnSplitEdge(splitEdgeInDirection_[1], j)] -
                    points[seDirJ.start()]
                ) /
                seDirJ.mag(points)
            );

            for (label k = 1; k < nLayersK_; ++k)
            {
                const scalar w
                (
                    Foam::mag
                    (
                        points[bndLayers_.vertexOnSplitEdge(splitEdgeInDirection_[2], k)] -
                        points[seDirK.start()]
                    ) /
                    seDirK.mag(points)
                );

                // calculate coordinates of the new vertex
                const point newP =
                    (1.0 - u)*(1.0 - v)*(1.0 - w)*v000 +
                    u*(1.0 - v)*(1.0 - w)*v100 +
                    u*v*(1.0 - w)*v110 +
                    (1.0 - u)*v*(1.0 - w)*v010 +
                    (1.0 - u)*(1.0 - v)*w*v001 +
                    u*(1.0 - v)*w*v101 +
                    u*v*w*v111 +
                    (1.0 - u)*v*w*v011;

                // add the point to the mesh
                cellPoints_[i][j][k] = points.size();
                points.append(newP);
            }
        }
    }

}


void Foam::Module::refineBoundaryLayers::refineCornerHexCell::
generateMissingFaces()
{
    // generate face in direction i
    for (label i = 1; i < nLayersI_; ++i)
    {
        // generate quad faces
        for (label j = 0; j < nLayersJ_; ++j)
        {
            for (label k = 0; k < nLayersK_; ++k)
            {
                // skip generating last face because it might not be a quad
                if ((j == (nLayersJ_ - 1)) && (k == (nLayersK_ - 1)))
                {
                    continue;
                }

                const label own
                (
                    k*nLayersI_*nLayersJ_ +
                    j*nLayersI_ +
                    i - 1
                );
                const label nei = own + 1;

                DynamicList<label> mf;
                mf.setSize(4);

                mf[0] = cellPoints_[i][j][k];
                mf[1] = cellPoints_[i][j + 1][k];
                mf[2] = cellPoints_[i][j + 1][k + 1];
                mf[3] = cellPoints_[i][j][k + 1];

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
        }

        // generate faces which might not be a quads
        DynamicList<label> mf;

        mf.append(cellPoints_[i][nLayersJ_ - 1][nLayersK_ - 1]);

        // this face might not be a quad
        // add points fom the last face in direction j
        const labelListList2D& f3 = facePoints_[3];
        for (label index = nLayersK_ - 1; index < f3.size()-1; ++index)
        {
            mf.append(f3[index][i]);
        }

        // add points from the last face in direction k
        const labelListList2D& f1 = facePoints_[1];
        for (label index = f1[i].size()-1; index>=nLayersJ_ - 1; --index)
        {
            mf.append(f1[i][index]);
        }

        const label own
        (
            (nLayersK_ - 1)*nLayersI_*nLayersJ_ +
            (nLayersJ_ - 1)*nLayersI_ +
            i - 1
        );

        const label nei = own + 1;

        cellsFromCell_[own].append(mf);
        cellsFromCell_[nei].append(help::reverseFace(mf));
    }

    // generate faces in direction j
    for (label j = 1; j < nLayersJ_; ++j)
    {
        // generate quad faces
        for (label i = 0; i < nLayersI_; ++i)
        {
            for (label k = 0; k < nLayersK_; ++k)
            {
                // skip generating late face because it might not be a quad
                if ((i == (nLayersI_ - 1)) && (k == (nLayersK_ - 1)))
                {
                    continue;
                }

                const label own
                (
                    k*nLayersI_*nLayersJ_ +
                    (j - 1)*nLayersI_ +
                    i
                );

                const label nei
                (
                    k*nLayersI_*nLayersJ_ +
                    j*nLayersI_ +
                    i
                );

                DynamicList<label> mf;
                mf.setSize(4);

                mf[0] = cellPoints_[i][j][k];
                mf[1] = cellPoints_[i][j][k + 1];
                mf[2] = cellPoints_[i + 1][j][k + 1];
                mf[3] = cellPoints_[i + 1][j][k];

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
        }

        // generate a face which might not be a quad
        DynamicList<label> mf;

        mf.append(cellPoints_[nLayersI_ - 1][j][nLayersK_ - 1]);

        // add points from the last face in direction k
        const labelListList2D& fp1 = facePoints_[1];
        for (label index = nLayersI_ - 1; index < fp1.size()-1; ++index)
        {
            mf.append(fp1[index][j]);
        }

        // add points from the last face in direction i
        const labelListList2D& fp5 = facePoints_[5];
        for (label index = fp5[j].size()-1; index>=nLayersK_ - 1; --index)
        {
            mf.append(fp5[j][index]);
        }

        const label own
        (
            (nLayersK_ - 1)*nLayersI_*nLayersJ_ +
            (j - 1)*nLayersI_ +
            (nLayersI_ - 1)
        );

        const label nei
        (
            (nLayersK_ - 1)*nLayersI_*nLayersJ_ +
            j*nLayersI_ +
            (nLayersI_ - 1)
        );

        cellsFromCell_[own].append(mf);
        cellsFromCell_[nei].append(help::reverseFace(mf));
    }

    // generate faces in direction k
    for (label k = 1; k < nLayersK_; ++k)
    {
        // generate quad faces
        for (label i = 0; i < nLayersI_; ++i)
        {
            for (label j = 0; j < nLayersJ_; ++j)
            {
                // skip the last face because it might not be a quad
                if ((i == (nLayersI_ - 1)) && (j == (nLayersJ_ - 1)))
                {
                    continue;
                }

                const label own
                (
                    (k - 1)*nLayersI_*nLayersJ_ +
                    j*nLayersI_ +
                    i
                );

                const label nei
                (
                    k*nLayersI_*nLayersJ_ +
                    j*nLayersI_ +
                    i
                );

                DynamicList<label> mf;
                mf.setSize(4);

                mf[0] = cellPoints_[i][j][k];
                mf[1] = cellPoints_[i + 1][j][k];
                mf[2] = cellPoints_[i + 1][j + 1][k];
                mf[3] = cellPoints_[i][j + 1][k];

                cellsFromCell_[own].append(mf);
                cellsFromCell_[nei].append(help::reverseFace(mf));
            }
        }

        // generate a face which might not be a quad
        DynamicList<label> mf;

        mf.append(cellPoints_[nLayersI_ - 1][nLayersJ_ - 1][k]);

        // this face might not be a quad
        // add points from the last face in direction i
        const labelListList2D& fp5 = facePoints_[5];
        for (label index = nLayersJ_ - 1; index < fp5.size()-1; ++index)
        {
            mf.append(fp5[index][k]);
        }

        // add points from the last face in direction j
        const labelListList2D& fp3 = facePoints_[3];
        for (label index = fp3[k].size()-1; index>=nLayersI_ - 1; --index)
        {
            mf.append(fp3[k][index]);
        }

        const label own
        (
            (k - 1)*nLayersI_*nLayersJ_ +
            (nLayersJ_ - 1)*nLayersI_ +
            (nLayersI_ - 1)
        );

        const label nei
        (
            k*nLayersI_*nLayersJ_ +
            (nLayersJ_ - 1)*nLayersI_ +
            (nLayersI_ - 1)
        );

        cellsFromCell_[own].append(mf);
        cellsFromCell_[nei].append(help::reverseFace(mf));
    }

}


Foam::Module::refineBoundaryLayers::refineCornerHexCell::refineCornerHexCell
(
    const label cellI,
    const refineBoundaryLayers& ref
)
:
    cellI_(cellI),
    nLayersI_(),
    nLayersJ_(),
    nLayersK_(),
    bndLayers_(ref)
{
    determineFacesInDirections();

    populateExistingFaces();

    generateNewPoints();

    generateMissingFaces();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::refineBoundaryLayers::generateNewCells()
{
    labelList nCellsFromCell(mesh_.cells().size(), 1);
    labelList refType(mesh_.cells().size(), 0);

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& faceOwners = mse.faceOwners();

    // calculate the number new cells generated from a cell
    forAll(faceOwners, bfI)
    {
        const label cellI = faceOwners[bfI];

        nCellsFromCell[cellI] *= nLayersAtBndFace_[bfI];

        if (nLayersAtBndFace_[bfI] > 1)
        {
            ++refType[cellI];
        }
    }

    // add cells which shall be refined in a subset
    if (cellSubsetName_ != "")
    {
        label subsetI = mesh_.cellSubsetIndex(cellSubsetName_);
        if (subsetI >= 0)
        {
            Warning << "The subset with name " << cellSubsetName_
                    << " already exists. Skipping generation of a new subset"
                    << endl;
        }

        subsetI = mesh_.addCellSubset(cellSubsetName_);

        forAll(nCellsFromCell, cI)
        {
            if (nCellsFromCell[cI] > 1)
            {
                mesh_.addCellToSubset(subsetI, cI);
            }
        }
    }

    // check the number of cells which will be generated
    label nNewCells(0);
    forAll(nCellsFromCell, cellI)
    {
        nNewCells += (nCellsFromCell[cellI] - 1);
    }

    const label totalNumNewCells = nNewCells;
    Info<< "Number of newly generated cells " << totalNumNewCells << endl;

    // create mesh modifier
    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();

    const label numFacesBefore = newFaces_.size();

    // set the number of cells to the new value
    cellListPMG& cells = meshModifier.cellsAccess();
    label nCells = cells.size();
    cells.setSize(nCells + nNewCells);

    // start creating new cells
    // store the information which new cells were generated from
    // an existing cell
    labelListList2D newCellsFromCell(refType.size());

    labelListList2D pointNewFaces;
    reverseAddressing(pointNewFaces, newFaces_);

    forAll(nCellsFromCell, cellI)
    {
        if (refType[cellI] == 0)
        {
            // this cell is not refined
            // update face labels
            newCellsFromCell[cellI].append(cellI);

            cell& c = cells[cellI];

            // copy the new faces of this cell
            DynamicList<label> newC;
            newC.reserve(64);
            forAll(c, fI)
            {
                forAllRow(facesFromFace_, c[fI], cfI)
                {
                    newC.append(facesFromFace_[c[fI]][cfI]);
                }
            }

            // update the cell
            c.setSize(newC.size());
            forAll(c, fI)
            {
                c[fI] = newC[fI];
            }
        }
        else if (refType[cellI] == 1)
        {
            // generate new cells from this prism refined in one direction
            List<DynamicList<DynamicList<label>>> cellsFromCell;
            generateNewCellsPrism(cellI, cellsFromCell);

            forAll(cellsFromCell, cI)
            {
                const DynamicList<DynamicList<label>>& nc = cellsFromCell[cI];

                const label newCellI = cI == 0?cellI:nCells++;

                newCellsFromCell[cellI].append(newCellI);

                cell& c = cells[newCellI];
                c.setSize(nc.size());

                // find face labels for this cell
                forAll(nc, fI)
                {
                    const DynamicList<label>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces[nf[0]][pfI];

                        if (help::areFacesEqual(nf, newFaces_[nfI]))
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if (faceLabel < 0)
                    {
                        forAll(nf, pI)
                        {
                            pointNewFaces[nf[pI]].append(newFaces_.size());
                        }
                        c[fI] = newFaces_.size();
                        appendList(newFaces_, nf);
                    }
                }
            }
        }
        else if (refType[cellI] == 2)
        {
            // generate new cell from a hex cell where two layers intersect
            // generate mostly hex cells;
            refineEdgeHexCell refEdgeHex(cellI, *this);
            const List<DynamicList<DynamicList<label>>>& cellsFromCell =
                refEdgeHex.newCells();

            forAll(cellsFromCell, cI)
            {
                const DynamicList<DynamicList<label>>& nc = cellsFromCell[cI];

                const label newCellI = cI == 0?cellI:nCells++;

                newCellsFromCell[cellI].append(newCellI);

                cell& c = cells[newCellI];
                c.setSize(nc.size());

                // find face labels for this cell
                forAll(nc, fI)
                {
                    const DynamicList<label>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces[nf[0]][pfI];

                        if (help::areFacesEqual(nf, newFaces_[nfI]))
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if (faceLabel < 0)
                    {
                        forAll(nf, pI)
                        {
                            pointNewFaces[nf[pI]].append(newFaces_.size());
                        }
                        c[fI] = newFaces_.size();
                        appendList(newFaces_, nf);
                    }
                }
            }
        }
        else if (refType[cellI] == 3)
        {
            // generate new cells from a hex at a corner where three
            // layers intersect
            // generate mostly hex cells
            refineCornerHexCell refCell(cellI, *this);
            const List<DynamicList<DynamicList<label>>>& cellsFromCell =
                refCell.newCells();

            // new points have been generated
            pointNewFaces.setSize(mesh_.points().size());

            // recognise face cells in the graph of new faces
            forAll(cellsFromCell, cI)
            {
                const DynamicList<DynamicList<label>>& nc = cellsFromCell[cI];

                const label newCellI = cI == 0?cellI:nCells++;

                newCellsFromCell[cellI].append(newCellI);

                cell& c = cells[newCellI];
                c.setSize(nc.size());

                // find face labels for this cell
                forAll(nc, fI)
                {
                    const DynamicList<label>& nf = nc[fI];

                    label faceLabel(-1);
                    forAllRow(pointNewFaces, nf[0], pfI)
                    {
                        const label nfI = pointNewFaces[nf[0]][pfI];

                        if (help::areFacesEqual(nf, newFaces_[nfI]))
                        {
                            c[fI] = nfI;
                            faceLabel = nfI;
                            break;
                        }
                    }

                    if (faceLabel < 0)
                    {
                        forAll(nf, pI)
                        {
                            pointNewFaces[nf[pI]].append(newFaces_.size());
                        }
                        c[fI] = newFaces_.size();
                        appendList(newFaces_, nf);
                    }
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot refine boundary layer for cell "
                << cellI << abort(FatalError);
        }
    }

    // check the orientation of new faces, because owner and neighbour cells
    // may require a face to be flipped
    const label nOrigInternalFaces = mesh_.nInternalFaces();

    # pragma omp parallel
    {
        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        # pragma omp for schedule(dynamic, 40)
        for (label fI = 0; fI < nOrigInternalFaces; ++fI)
        {
            const label own = owner[fI];
            const label nei = neighbour[fI];

            if (facesFromFace_[fI].size() == 1)
            {
                continue;
            }

            forAllRow(facesFromFace_, fI, cfI)
            {
                const label nfI = facesFromFace_[fI][cfI];

                // find the new owner and neighbour cells of the new face
                label newOwner(-1), newNeighbour(-1);
                forAllRow(newCellsFromCell, own, cI)
                {
                    const cell& cOwn = cells[newCellsFromCell[own][cI]];

                    const label pos = help::positionInList(nfI, cOwn);

                    if (pos >= 0)
                    {
                        newOwner = newCellsFromCell[own][cI];
                        break;
                    }
                }

                forAllRow(newCellsFromCell, nei, cI)
                {
                    const cell& cNei = cells[newCellsFromCell[nei][cI]];

                    const label pos = help::positionInList(nfI, cNei);

                    if (pos >= 0)
                    {
                        newNeighbour = newCellsFromCell[nei][cI];
                        break;
                    }
                }

                if (newOwner > newNeighbour)
                {
                    DynamicList<label> rf;
                    rf.setSize(newFaces_[nfI].size());

                    forAll(rf, i)
                    {
                        rf[i] = newFaces_[nfI][i];
                    }

                    rf = help::reverseFace(rf);

                    setRow(newFaces_, nfI, rf);
                }
            }
        }
    }

    // update cell sets
    mesh_.updateCellSubsets(newCellsFromCell);
    newCellsFromCell.setSize(0);

    // point-faces addressing is not needed any more
    pointNewFaces.setSize(0);

    // copy newFaces to the mesh

    // store internal faces originating from existing faces
    DynamicList<label> newFaceLabel;
    newFaceLabel.setSize(newFaces_.size());
    faces.setSize(newFaces_.size());

    label currFace = 0;
    for (label faceI = 0; faceI < nOrigInternalFaces; ++faceI)
    {
        forAllRow(facesFromFace_, faceI, ffI)
        {
            face& f = faces[currFace];
            newFaceLabel[currFace] = currFace;
            ++currFace;

            const label newFaceI = facesFromFace_[faceI][ffI];

            f.setSize(newFaces_[newFaceI].size());

            forAll(f, pI)
            {
                f[pI] = newFaces_[newFaceI][pI];
            }
        }
    }

    // store newly-generated internal faces

    for (label faceI = numFacesBefore; faceI < newFaces_.size(); ++faceI)
    {
        newFaceLabel[faceI] = currFace;
        face& f = faces[currFace];
        ++currFace;

        f.setSize(newFaces_[faceI].size());

        forAll(f, pI)
        {
            f[pI] = newFaces_[faceI][pI];
        }
    }

    // store new boundary faces

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label size = boundaries[patchI].patchSize();

        const label newStart = currFace;
        label nNewFacesInPatch(0);
        for (label fI = 0; fI < size; ++fI)
        {
            const label faceI = start + fI;

            forAllRow(facesFromFace_, faceI, nfI)
            {
                face& f = faces[currFace];

                // update the new label
                const label origFaceI = facesFromFace_[faceI][nfI];
                newFaceLabel[origFaceI] = currFace;
                facesFromFace_[faceI][nfI] = currFace;
                ++currFace;

                // copy the face into the mesh
                f.setSize(newFaces_[origFaceI].size());
                forAll(f, pI)
                {
                    f[pI] = newFaces_[origFaceI][pI];
                }

                ++nNewFacesInPatch;
            }
        }

        // update patch
        boundaries[patchI].patchStart() = newStart;
        boundaries[patchI].patchSize() = nNewFacesInPatch;
    }

    // update face subsets
    mesh_.updateFaceSubsets(facesFromFace_);
    facesFromFace_.setSize(0);
    newFaces_.setSize(0);

    // update cells to match the faces

    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
        {
            c[fI] = newFaceLabel[c[fI]];
        }
    }

    // delete all adressing which is no longer up-to-date
    meshModifier.clearAll();
    deleteDemandDrivenData(msePtr_);

    Info<< "Finished generating new cells " << endl;
}


// ************************************************************************* //
