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

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"
#include "DynamicList.H"

#include "HashSet.H"

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::boundaryLayers::createLayerCells
(
    const labelList& patchLabels
)
{
    if (omp_in_parallel())
    {
        FatalErrorInFunction
            << "Building cellsToAdd inside a parallel region."
            << " This is not thread safe" << exit(FatalError);
    }

    Info<< "Starting creating layer cells" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const edgeList& edges = mse.edges();
    const labelListList2D& faceEdges = mse.faceEdges();
    const labelListList2D& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const labelList& faceOwners = mse.faceOwners();
    const labelList& bp = mse.bp();
    const labelListList2D& pointFaces = mse.pointFaces();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const labelListList2D& pointPatches = mPart.pointPatches();

    // mark patches which will be extruded into layer cells
    boolList treatPatches(mesh_.boundaries().size(), false);
    forAll(patchLabels, patchI)
    {
        const label pLabel = patchLabels[patchI];
        forAll(treatPatchesWithPatch_[pLabel], i)
        {
            treatPatches[treatPatchesWithPatch_[pLabel][i]] = true;
        }
    }

    // create lists for new boundary faces
    DynamicList<DynamicList<label>> newBoundaryFaces;
    DynamicList<label> newBoundaryOwners;
    DynamicList<label> newBoundaryPatches;

    // create storage for new cells
    DynamicList<faceList> cellsToAdd;

    // create layer cells and store boundary faces
    const label nOldCells = mesh_.cells().size();
    forAll(bFaces, bfI)
    {
        if (treatPatches[boundaryFacePatches[bfI]])
        {
            const face& f = bFaces[bfI];
            const label pKey = patchKey_[boundaryFacePatches[bfI]];

            DynamicList<face> cellFaces;

            DynamicList<label> newF;

            // store the current boundary face
            newF.clear();
            newF.append(f[0]);
            for (label pI = f.size()-1; pI > 0; --pI)
            {
                newF.append(f[pI]);
            }

            {
                face nf(newF.size());
                forAll(newF, i)
                {
                    nf[i] = newF[i];
                }
                cellFaces.append(nf);
            }

            // create parallel face
            forAll(f, pI)
            {
                newF[pI] = findNewNodeLabel(f[pI], pKey);
            }

            {
                face nf(newF.size());
                forAll(newF, i)
                {
                    nf[i] = newF[i];
                }
                cellFaces.append(nf);
            }

            appendList(newBoundaryFaces,newF);
            newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);

            // create quad faces
            newF.setSize(4);
            forAll(f, pI)
            {
                newF[0] = f[pI];
                newF[1] = f.nextLabel(pI);
                newF[2] = findNewNodeLabel(newF[1], pKey);
                newF[3] = findNewNodeLabel(f[pI], pKey);

                {
                    face nf(newF.size());
                    forAll(newF, i)
                    {
                        nf[i] = newF[i];
                    }
                    cellFaces.append(nf);
                }

                // check if the face is at the boundary
                // of the treated partitions
                const label edgeI = faceEdges[bfI][pI];
                if (edgeFaces[edgeI].size() == 2)
                {
                    label neiFace = edgeFaces[edgeI][0];
                    if (neiFace == bfI)
                    {
                        neiFace = edgeFaces[edgeI][1];
                    }

                    if (!treatPatches[boundaryFacePatches[neiFace]])
                    {
                        appendList(newBoundaryFaces,newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(boundaryFacePatches[neiFace]);
                    }
                }
            }

            cellsToAdd.append(cellFaces);
        }
        else
        {
            appendList(newBoundaryFaces,bFaces[bfI]);
            newBoundaryOwners.append(faceOwners[bfI]);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);
        }
    }

    // create cells at edges
    forAll(edgeFaces, edgeI)
    {
        // do not consider edges with no faces attached to it
        if (edgeFaces[edgeI].size() == 0)
        {
            continue;
        }

        // check if the edge is a feature edge
        const label patchI = boundaryFacePatches[edgeFaces[edgeI][0]];

        const label patchJ = boundaryFacePatches[edgeFaces[edgeI][1]];

        if (patchI == patchJ)
        {
            continue;
        }

        // check if the faces attached to the edge have different keys
        const label pKeyI = patchKey_[patchI];
        const label pKeyJ = patchKey_[patchJ];

        if (pKeyI < 0 || pKeyJ < 0)
        {
            continue;
        }

        if (pKeyI == pKeyJ)
        {
            continue;
        }

        const edge& e = edges[edgeI];
        if (!otherVrts_.found(e.start()))
        {
            continue;
        }
        if (!otherVrts_.found(e.end()))
        {
            continue;
        }

        // generate faces of the bnd layer cell
        faceList cellFaces;
        createNewCellFromEdge(e, pKeyI, pKeyJ, cellFaces);

        // store boundary faces
        appendList(newBoundaryFaces,cellFaces[1]);
        newBoundaryOwners.append(nOldCells + cellsToAdd.size());
        newBoundaryPatches.append(patchJ);

        appendList(newBoundaryFaces,cellFaces[3]);
        newBoundaryOwners.append(nOldCells + cellsToAdd.size());
        newBoundaryPatches.append(patchI);

        // check if face 5 is a boundary face
        const label bps = bp[e.start()];
        label unusedPatch(-1);
        forAllRow(pointPatches, bps, i)
        {
            const label ptchI = pointPatches[bps][i];

            if (ptchI == patchI)
            {
                continue;
            }
            if (ptchI == patchJ)
            {
                continue;
            }
            if (unusedPatch != -1)
            {
                unusedPatch = -1;
                break;
            }

            unusedPatch = ptchI;
        }

        if (unusedPatch != -1 && treatedPatch_[unusedPatch])
        {
            // add a face in the empty patch in case of 2D layer generation
            appendList(newBoundaryFaces,cellFaces[5]);
            newBoundaryOwners.append(nOldCells + cellsToAdd.size());
            newBoundaryPatches.append(unusedPatch);
        }

        // check if face 4 is a boundary face
        const label bpe = bp[e.end()];
        unusedPatch = -1;
        forAllRow(pointPatches, bpe, i)
        {
            const label ptchI = pointPatches[bpe][i];

            if (ptchI == patchI)
            {
                continue;
            }
            if (ptchI == patchJ)
            {
                continue;
            }
            if (unusedPatch != -1)
            {
                unusedPatch = -1;
                break;
            }

            unusedPatch = ptchI;
        }

        if (unusedPatch != -1 && treatedPatch_[unusedPatch])
        {
            // add a face in the empty patch in case of 2D layer generation
            appendList(newBoundaryFaces,cellFaces[4]);
            newBoundaryOwners.append(nOldCells + cellsToAdd.size());
            newBoundaryPatches.append(unusedPatch);
        }

        // append cell to the queue
        cellsToAdd.append(cellFaces);
    }

    // create cells for corner nodes
    typedef Map<DynamicList<label>> ppType;

    ppType nodePatches;

    forAllConstIter(Map<EdgeMap<label>>, otherVrts_, iter)
    {
        if (iter().size() == 2)
        {
            continue;
        }

        const label bpI = bp[iter.key()];

        // ensure correct orientation
        DynamicList<label> pFaces;
        pFaces.setSize(pointFaces[bpI].size());
        forAll(pFaces, fI)
        {
            pFaces[fI] = pointFaces[bpI][fI];
        }

        for (label i = 0; i < pFaces.size(); ++i)
        {
            const face& bf = bFaces[pFaces[i]];
            const edge e = bf.faceEdge(bf.which(iter.key()));

            for (label j = i + 1; j < pFaces.size(); ++j)
            {
                const face& obf = bFaces[pFaces[j]];
                if
                (
                    (obf.which(e.start()) >= 0)
                 && (obf.which(e.end()) >= 0)
                )
                {
                    const label add = pFaces[i + 1];
                    pFaces[i + 1] = pFaces[j];
                    pFaces[j] = add;
                    break;
                }
            }
        }

        DynamicList<label> patchIDs;
        patchIDs.reserve(3);
        forAll(pFaces, patchI)
        {
            appendUniq(patchIDs, boundaryFacePatches[pFaces[patchI]]);
        }

        nodePatches.insert(iter.key(), patchIDs);
    }

    // create layer cells for corner nodes
    const labelList nodePatchesKeys(nodePatches.sortedToc());
    forAll(nodePatchesKeys, keyI)
    {
        const label nodeI = nodePatchesKeys[keyI];
        const DynamicList<label>& patchIDs = nodePatches[nodeI];
        DynamicList<label> pKeys;
        pKeys.reserve(3);
        forAll(patchIDs, patchI)
        {
            const label pKey = patchKey_[patchIDs[patchI]];

            if (pKey < 0)
            {
                continue;
            }

            appendUniq(pKeys, pKey);
        }

        if (pKeys.size() != 3)
        {
            continue;
        }

        faceList cellFaces;
        createNewCellFromNode(nodeI, pKeys, cellFaces);

        // store boundary faces
        appendList(newBoundaryFaces,cellFaces[1]);
        newBoundaryOwners.append(nOldCells + cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[0]);

        appendList(newBoundaryFaces,cellFaces[3]);
        newBoundaryOwners.append(nOldCells + cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[1]);

        appendList(newBoundaryFaces,cellFaces[5]);
        newBoundaryOwners.append(nOldCells + cellsToAdd.size());
        newBoundaryPatches.append(patchIDs[2]);

        // append cell to the queue
        cellsToAdd.append(cellFaces);
    }

    // create mesh modifier
    polyMeshGenModifier meshModifier(mesh_);

    meshModifier.addCells(cellsToAdd);

    cellsToAdd.clear();
    meshModifier.reorderBoundaryFaces();
    meshModifier.replaceBoundary
    (
        patchNames_,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );

    PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();
    forAll(boundaries, patchI)
    {
        boundaries[patchI].patchType() = patchTypes_[patchI];
    }

    // delete meshSurfaceEngine
    this->clearOut();

    Info<< "Finished creating layer cells" << endl;
}


// ************************************************************************* //
