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
#include "demandDrivenData.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::refineBoundaryLayers::refineFace
(
    const face& f,
    const FixedList<label, 2>& nLayersInDirection,
    DynamicList<DynamicList<label>>& newFaces
)
{
    // this face must be a quad
    if (f.size() != 4)
    {
        WarningInFunction
            << "Face " << f << " is not a quad" << endl;
        return;
    }

    // direction 0 represents edges 0 and 2
    // direction 1 represents edges 1 and 3
    if ((nLayersInDirection[0] <= 1) && (nLayersInDirection[1] <= 1))
    {
        // this face may comprise of some split edges
        DynamicList<label> newF;
        newF.reserve(64);
        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            // add the current point label
            newF.append(f[eI]);

            // check if a split edge matches this face edge
            forAllRow(splitEdgesAtPoint_, f[eI], peI)
            {
                const label seI = splitEdgesAtPoint_[f[eI]][peI];
                const edge& se = splitEdges_[seI];

                if (e == se)
                {
                    // check the orientation and add new vertices created
                    // on this edge
                    const label s = newVerticesForSplitEdge_[seI].size() - 1;
                    if (e.start() == se.start())
                    {
                        for (label pI = 1; pI < s; ++pI)
                        {
                            newF.append(newVerticesForSplitEdge_[seI][pI]);
                        }
                    }
                    else
                    {
                        for (label pI = s - 1; pI > 0; --pI)
                        {
                            newF.append(newVerticesForSplitEdge_[seI][pI]);
                        }
                    }
                }
            }
        }

        newFaces.setSize(1);
        newFaces[0] = newF;
        return;
    }

    // check which face edge is a direction 0 and which one is a direction 1
    label dir0(-1), dir1(-1);
    labelPair dir0Edges(-1, -1), dir1Edges(-1, -1);
    forAll(f, eI)
    {
        const edge e = f.faceEdge(eI);

        label ses(-1), see(-1);
        bool start(false), end(false);
        forAllRow(splitEdgesAtPoint_, e.start(), i)
        {
            const edge& se = splitEdges_[splitEdgesAtPoint_[e.start()][i]];

            if ((se.start() == e.start()) && (se.end() == f.prevLabel(eI)))
            {
                ses = splitEdgesAtPoint_[e.start()][i];
                start = true;
                break;
            }
        }

        forAllRow(splitEdgesAtPoint_, e.end(), i)
        {
            const edge& se = splitEdges_[splitEdgesAtPoint_[e.end()][i]];

            if ((se.start() == e.end()) && (se.end() == f[(eI + 2)%4]))
            {
                see = splitEdgesAtPoint_[e.end()][i];
                end = true;
                break;
            }
        }

        if (start && end)
        {
            if (dir0 == -1)
            {
                dir0 = eI;
                dir0Edges = labelPair(ses, see);
            }
            else if (dir1 == -1)
            {
                dir1 = eI;
                dir1Edges = labelPair(ses, see);
            }
            else
            {
                FatalErrorInFunction
                    << "More than two split directions for a face"
                    << abort(FatalError);
            }
        }
    }

    if ((dir0 < 0) && (dir1 < 0))
    {
        Pout << "Refining face " << f << endl;
        forAll(f, pI)
        {
            if (splitEdgesAtPoint_.size() >= f[pI])
            {
                Pout << "Split edges at point " << f[pI]
                     << " are " << splitEdgesAtPoint_[f[pI]] << endl;
            }
        }
        Pout << "Splits in direction " << nLayersInDirection << endl;
        Pout << "Here " << endl;
        Pout << "Dir0 " << dir0 << endl;
        Pout << "dir0Edges " << dir0Edges << endl;
        Pout << "Dir1 " << dir1 << endl;
        Pout << "dir1Edges " << dir1Edges << endl;

        FatalErrorInFunction
            << "Cannot find split edges for a face" << abort(FatalError);
    }

    // in case of only one refinement direction, it must direction 0
    if ((dir1 != -1) && (dir0 == -1))
    {
        dir0 = dir1;
        dir0Edges = dir1Edges;
        dir1 = -1;
    }
    else if ((dir0 != -1) && (dir1 != -1) && (dir1 != f.fcIndex(dir0)))
    {
        // alternate value to preserve correct face orientation
        const label add = dir0;
        dir0 = dir1;
        dir1 = add;

        const labelPair lpAdd = dir0Edges;
        dir0Edges = dir1Edges;
        dir1Edges = lpAdd;
    }

    // permutate the number of refinements in each direction
    const label nLayersDir0 = dir0>=0?nLayersInDirection[dir0%2]:1;
    const label nLayersDir1 = dir1>=0?nLayersInDirection[dir1%2]:1;

    // map the face onto a matrix for easier orientation
    labelListList2D facePoints;
    facePoints.setSize(nLayersDir0 + 1);
    forAll(facePoints, i)
    {
        facePoints[i].setSize(nLayersDir1 + 1);
        facePoints[i] = -1;
    }

    // add points in the matrix
    for (label i = 0; i < nLayersDir0; ++i)
    {
        facePoints[i][0] = newVerticesForSplitEdge_[dir0Edges.second()][i];
        facePoints[i][nLayersDir1] =
            newVerticesForSplitEdge_[dir0Edges.first()][i];
    }
    facePoints[nLayersDir0][0] = splitEdges_[dir0Edges.second()].end();
    facePoints[nLayersDir0][nLayersDir1] = splitEdges_[dir0Edges.first()].end();

    for (label i = 1; i < nLayersDir1; ++i)
    {
        facePoints[0][i] = newVerticesForSplitEdge_[dir1Edges.first()][i];
        facePoints[nLayersDir0][i] =
            newVerticesForSplitEdge_[dir1Edges.second()][i];
    }

    // create missing vertices if there are any
    pointFieldPMG& points = mesh_.points();
    const point v00 = points[facePoints[0][0]];
    const point v10 = points[facePoints[nLayersDir0][0]];
    const point v01 = points[facePoints[0][nLayersDir1]];
    const point v11 = points[facePoints[nLayersDir0][nLayersDir1]];

    forAll(facePoints, i)
    {
        forAll(facePoints[i], j)
        {
            if (facePoints[i][j] < 0)
            {

                const scalar u
                (
                    Foam::mag(points[facePoints[i][0]] - v00) /
                    splitEdges_[dir0Edges.second()].mag(points)
                );

                const scalar v
                (
                    Foam::mag(points[facePoints[0][j]] - v00) /
                    splitEdges_[dir1Edges.first()].mag(points)
                );

                // calculate the coordinates of the missing point via
                // transfinite interpolation
                const point newP
                (
                    (1.0 - u)*(1.0 - v)*v00 +
                    u*(1.0 - v)*v10 +
                    u*v*v11 +
                    (1.0 - u)*v*v01
                );

                // add the vertex to the mesh
                facePoints[i][j] = points.size();
                points.append(newP);
            }
        }
    }

    // Finally, create the faces
    for (label j = 0; j < nLayersDir1; ++j)
    {
        for (label i = 0; i < nLayersDir0; ++i)
        {
            // create quad face
            DynamicList<label> f;
            f.reserve(4);

            f.append(facePoints[i][j]);

            if ((i == (nLayersDir0 - 1)) && (j == 0))
            {

                // add additional points on edge
                const label eLabel = dir0Edges.second();
                const label size =
                    newVerticesForSplitEdge_[eLabel].size() - 1;

                for (label index = i + 1; index < size; ++index)
                {
                    f.append(newVerticesForSplitEdge_[eLabel][index]);
                }
            }

            f.append(facePoints[i + 1][j]);

            if
            (
                (dir1 != -1) &&
                (i == (nLayersDir0 - 1)) &&
                (j == (nLayersDir1 - 1))
            )
            {

                // add additional points on edge
                const label eLabel = dir1Edges.second();
                const label size =
                    newVerticesForSplitEdge_[eLabel].size() - 1;

                for (label index = j + 1; index < size; ++index)
                {
                    f.append(newVerticesForSplitEdge_[eLabel][index]);
                }
            }

            f.append(facePoints[i + 1][j + 1]);

            if ((i == (nLayersDir0 - 1)) && (j == (nLayersDir1 - 1)))
            {

                const label eLabel = dir0Edges.first();
                const label size =
                    newVerticesForSplitEdge_[eLabel].size() - 2;
                for (label index = size; index > i; --index)
                {
                    f.append(newVerticesForSplitEdge_[eLabel][index]);
                }
            }

            f.append(facePoints[i][j + 1]);

            if ((dir1 != -1) && (i == 0) && (j == (nLayersDir1 - 1)))
            {

                const label eLabel = dir1Edges.first();
                const label size =
                    newVerticesForSplitEdge_[eLabel].size() - 2;
                for (label index = size; index > j; --index)
                {
                    f.append(newVerticesForSplitEdge_[eLabel][index]);
                }
            }

            newFaces.append(f);
        }
    }

}


void Foam::Module::refineBoundaryLayers::sortFacePoints
(
    const label faceI,
    labelListList2D& facePoints,
    const label transpose
) const
{
    const faceListPMG& faces = mesh_.faces();
    const face& f = faces[faceI];

    // start the procedure by finding the number of splits in
    // both i and j direction
    label numSplitsI(1);

    const label pos = f.which(newFaces_[facesFromFace_[faceI][0]][0]);

    forAllRow(facesFromFace_, faceI, i)
    {
        const label nfI = facesFromFace_[faceI][i];

        if ((numSplitsI == 1) && found(newFaces_, nfI, f.nextLabel(pos)))
        {
            numSplitsI = i + 1;
            break;
        }
    }

    const label numSplitsJ = (facesFromFace_[faceI].size() / numSplitsI);

    facePoints.setSize(numSplitsI + 1);
    forAll(facePoints, i)
    {
        facePoints[i].setSize(numSplitsJ + 1);
    }

    // start filling in the matrix
    forAllRow(facesFromFace_, faceI, fI)
    {
        const label nfI = facesFromFace_[faceI][fI];

        const label i = fI % numSplitsI;
        const label j = fI/numSplitsI;

        if (newFaces_[nfI].size() == 4)
        {
            facePoints[i][j] = newFaces_[nfI][0];
            facePoints[i + 1][j] = newFaces_[nfI][1];
            facePoints[i + 1][j + 1] = newFaces_[nfI][2];
            facePoints[i][j + 1] = newFaces_[nfI][3];
        }
        else
        {
            if (j == 0)
            {
                forAllRow(newFaces_, nfI, pI)
                {
                    if (f.which(newFaces_[nfI][pI]) >= 0)
                    {
                        facePoints[i + 1][0] = newFaces_[nfI][pI];
                        break;
                    }
                }
            }
            else if (i == 0)
            {
                forAllRow(newFaces_, nfI, pI)
                {
                    if (f.which(newFaces_[nfI][pI]) >= 0)
                    {
                        facePoints[0][j + 1] = newFaces_[nfI][pI];
                        break;
                    }
                }
            }
            else
            {
                forAllRow(newFaces_, nfI, pI)
                {
                    if (f.which(newFaces_[nfI][pI]) >= 0)
                    {
                        facePoints[i + 1][j + 1] = newFaces_[nfI][pI];
                        break;
                    }
                }
            }
        }
    }

    if (transpose)
    {
        labelListList2D transposedFacePoints;
        transposedFacePoints.setSize(facePoints[0].size());
        forAll(transposedFacePoints, j)
        {
            transposedFacePoints[j].setSize(facePoints.size());
        }

        forAll(facePoints, i)
        {
            forAll(facePoints[i], j)
            {
                transposedFacePoints[j][i] = facePoints[i][j];
            }
        }

        facePoints = transposedFacePoints;

    }
}


void Foam::Module::refineBoundaryLayers::sortFaceFaces
(
    const label faceI,
    labelListList2D& faceFaces,
    const label transpose
) const
{
    const faceListPMG& faces = mesh_.faces();
    const face& f = faces[faceI];

    // start the procedure by finding the number of splits in
    // both i and j direction
    label numSplitsI(1);

    const label pos = f.which(newFaces_[facesFromFace_[faceI][0]][0]);

    forAllRow(facesFromFace_, faceI, i)
    {
        const label nfI = facesFromFace_[faceI][i];

        if ((numSplitsI == 1) && found(newFaces_, nfI, f.nextLabel(pos)))
        {
            numSplitsI = i + 1;
            break;
        }
    }

    label numSplitsJ = (facesFromFace_[faceI].size() / numSplitsI);

    faceFaces.setSize(numSplitsI);
    forAll(faceFaces, i)
    {
        faceFaces[i].setSize(numSplitsJ);
    }

    // start filling in the matrix
    forAllRow(facesFromFace_, faceI, fI)
    {
        const label nfI = facesFromFace_[faceI][fI];

        const label i = fI % numSplitsI;
        const label j = fI/numSplitsI;

        faceFaces[i][j] = nfI;
    }

    if (transpose)
    {
        labelListList2D transposedFaceFaces;
        transposedFaceFaces.setSize(faceFaces[0].size());
        forAll(transposedFaceFaces, j)
        {
            transposedFaceFaces[j].setSize(faceFaces.size());
        }

        forAll(faceFaces, i)
        {
            forAll(faceFaces[i], j)
            {
                transposedFaceFaces[j][i] = faceFaces[i][j];
            }
        }


        faceFaces = transposedFaceFaces;

    }
}


void Foam::Module::refineBoundaryLayers::generateNewFaces()
{
    // generate new boundary and inter-processor faces
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatches = mse.boundaryFacePatches();
    const edgeList& edges = mse.edges();
    const labelList& bp = mse.bp();
    const labelListList2D& bfEdges = mse.faceEdges();
    const labelListList2D& bpEdges = mse.boundaryPointEdges();
    const labelListList2D& beFaces = mse.edgeFaces();

    // mesh data
    const label nInternalFaces = mesh_.nInternalFaces();
    const faceListPMG& faces = mesh_.faces();

    // container for faces
    facesFromFace_.setSize(faces.size());
    newFaces_.clear();

    // split internal faces
    for (label faceI = 0; faceI < nInternalFaces; ++faceI)
    {
        const face& f = faces[faceI];

        // only quad faces can be split
        if (f.size() != 4)
        {
            facesFromFace_[faceI].append(newFaces_.size());
            appendList(newFaces_, f);
            continue;
        }

        // check if there exist an edge of the face at the boundary
        FixedList<label, 2> nRefinementInDirection(1);

        forAll(f, eI)
        {
            const edge fe = f.faceEdge(eI);

            const label bps = bp[fe.start()];

            if (bps < 0)
            {
                continue;
            }

            forAllRow(bpEdges, bps, bpsI)
            {
                const label beI = bpEdges[bps][bpsI];

                if (edges[beI] == fe)
                {
                    // this edge is attached to the boundary
                    // get the number of layers for neighbouring cells
                    const label nSplits0 = nLayersAtBndFace_[beFaces[beI][0]];
                    const label nSplits1 = nLayersAtBndFace_[beFaces[beI][1]];

                    // set the number of layers for the given direction
                    const label dir = eI % 2;
                    nRefinementInDirection[dir] =
                        Foam::max
                        (
                            nRefinementInDirection[dir],
                            Foam::max(nSplits0, nSplits1)
                        );
                }
            }
        }

        // refine the face
        DynamicList<DynamicList<label>> newFacesForFace;
        refineFace(f, nRefinementInDirection, newFacesForFace);

        // store decomposed faces
        forAll(newFacesForFace, fI)
        {
            facesFromFace_[faceI].append(newFaces_.size());
            appendList(newFaces_, newFacesForFace[fI]);
        }

    }

    // refine boundary faces where needed
    // it is required in locations where two or three layers intersect
    const label startingBoundaryFace = mesh_.boundaries()[0].patchStart();
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];
        const label faceI = startingBoundaryFace + bfI;

        // only quad faces can be split
        if (bf.size() != 4)
        {
            facesFromFace_[faceI].append(newFaces_.size());
            appendList(newFaces_, bf);
            continue;
        }

        // check whether this face shall be refined and in which directions
        FixedList<label, 2> nRefinementInDirection(1);

        forAll(bf, eI)
        {
            const label beI = bfEdges[bfI][eI];

            if (beFaces[beI].size() != 2)
            {
                continue;
            }

            // get the neighbour face over the edge
            label neiFace = beFaces[beI][0];

            if (neiFace == bfI)
            {
                neiFace = beFaces[beI][1];
            }

            // faces cannot be in the same layer
            const DynamicList<label>& neiLayers =
                layerAtPatch_[facePatches[neiFace]];

            if (neiLayers.size() == 0)
            {
                continue;
            }

            const DynamicList<label>& currLayers = layerAtPatch_[facePatches[bfI]];

            bool foundSame(false);

            forAll(currLayers, i)
            {
                if (findIndex(neiLayers, currLayers[i]) != -1)
                {
                    foundSame = true;
                    break;
                }
            }

            if (foundSame || (neiLayers.size() == 0))
            {
                continue;
            }

            // set the refinement direction for this face
            nRefinementInDirection[eI%2] = nLayersAtBndFace_[neiFace];
        }

        // refine the face
        DynamicList<DynamicList<label>> newFacesForFace;
        refineFace(bf, nRefinementInDirection, newFacesForFace);

        // store the refined faces
        forAll(newFacesForFace, fI)
        {
            facesFromFace_[faceI].append(newFaces_.size());
            appendList(newFaces_, newFacesForFace[fI]);
        }

    }

    Info<< "Finished refining boundary - layer faces " << endl;
}


// ************************************************************************* //
