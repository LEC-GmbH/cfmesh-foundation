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

#include <map>

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::boundaryLayers::findPatchVertices
(
    const boolList& treatPatches,
    List<direction>& pVertices
) const
{
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    pVertices.setSize(pPatches.size());
    pVertices = NONE;

    # pragma omp parallel for if (pPatches.size() > 1000) \
    schedule(dynamic, Foam::max(10, pPatches.size()/(2*omp_get_num_threads())))
    forAll(pPatches, bpI)
    {
        bool hasTreated(false);
        bool hasNotTreated(false);

        forAllRow(pPatches, bpI, patchI)
        {
            const label patch = pPatches(bpI, patchI);
            if (treatPatches[patch])
            {
                hasTreated = true;
            }
            else
            {
                hasNotTreated = true;
            }
        }

        if (hasTreated)
        {
            pVertices[bpI] |= PATCHNODE;

            if (hasNotTreated)
            {
                pVertices[bpI] |= EDGENODE;
            }
        }
    }

}


Foam::point Foam::Module::boundaryLayers::createNewVertex
(
    const label bpI,
    const boolList& treatPatches,
    const List<direction>& patchVertex
) const
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorField& pNormals = mse.pointNormals();
    const VRWGraph& pFaces = mse.pointFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const VRWGraph& pointPoints = mse.pointPoints();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    const pointFieldPMG& points = mesh_.points();

    vector normal(vector::zero);
    scalar dist(VGREAT);
    const point& p = points[bPoints[bpI]];
    if (patchVertex[bpI] & EDGENODE)
    {
        DynList<label> otherPatches;
        forAllRow(pPatches, bpI, patchI)
        {
            if (!treatPatches[pPatches(bpI, patchI)])
            {
                otherPatches.appendUniq
                (
                    pPatches(bpI, patchI)
                );
            }
        }

        if (otherPatches.size() == 1)
        {
            // vertex is on an edge
            vector v(vector::zero);

            forAllRow(pFaces, bpI, pfI)
            {
                const face& f = bFaces[pFaces(bpI, pfI)];
                const label patchLabel =
                    boundaryFacePatches[pFaces(bpI, pfI)];

                if (treatPatches[patchLabel])
                {
                    normal += f.areaNormal(points);
                }
                else
                {
                    v += f.areaNormal(points);
                }
            }

            v.normalise();

            normal -= (normal & v)*v;

            normal.normalise();

            forAllRow(pointPoints, bpI, ppI)
            {
                if (patchVertex[pointPoints(bpI, ppI)])
                {
                    continue;
                }

                const vector vec = points[bPoints[pointPoints(bpI, ppI)]] - p;
                const scalar prod = 0.5*mag(vec & normal);

                if (prod < dist)
                {
                    dist = prod;
                }
            }
        }
        else if (otherPatches.size() == 2)
        {
            label otherVertex(-1);
            forAllRow(pointPoints, bpI, ppI)
            {
                const label bpJ = pointPoints(bpI, ppI);

                bool found(true);
                forAll(otherPatches, opI)
                {
                    if (!pPatches.found(bpJ, otherPatches[opI]))
                    {
                        found = false;
                        break;
                    }
                }

                if (found)
                {
                    otherVertex = bpJ;
                    break;
                }
            }

            if (otherVertex == -1)
            {
                FatalErrorInFunction
                    << "Cannot find moving vertex!" << exit(FatalError);
            }

            // normal vector is co-linear with that edge
            normal = p - points[bPoints[otherVertex]];
            dist = 0.5*mag(normal) + VSMALL;

            normal /= 2.0*dist;
        }
        else
        {
            FatalErrorInFunction
                << "There are more than 3 patches meeting at this vertex!"
                << pPatches[bpI] << abort(FatalError);
        }

        // limit distances
        forAllRow(pFaces, bpI, pfI)
        {
            const label faceLabel = pFaces(bpI, pfI);
            if (otherPatches.found(boundaryFacePatches[faceLabel]))
            {
                const face& f = bFaces[faceLabel];
                const label pos = f.which(bPoints[bpI]);

                if (pos != -1)
                {
                    const point& ep1 = points[f.prevLabel(pos)];
                    const point& ep2 = points[f.nextLabel(pos)];

                    const scalar dst =
                        help::distanceOfPointFromTheEdge(ep1, ep2, p);

                    if (dst < dist)
                    {
                        dist = 0.9*dst;
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << "Face does not contains this vertex!"
                        << abort(FatalError);
                }
            }
        }
    }
    else
    {
        normal = pNormals[bpI];

        forAllRow(pointPoints, bpI, ppI)
        {
            const scalar d =
                0.5*mag
                (
                    points[bPoints[pointPoints(bpI, ppI)]] -
                    p
                );

            if (d < dist)
            {
                dist = d;
            }
        }
    }

    // create new vertex
    dist = Foam::max(dist, VSMALL);

    const point newP = p - dist*normal;

    if (help::isnan(newP) || help::isinf(newP))
    {
        return p;
    }

    return newP;
}


void Foam::Module::boundaryLayers::createNewVertices
(
    const boolList& treatPatches
)
{
    Info<< "Creating vertices for layer cells" << endl;

    List<direction> patchVertex;
    findPatchVertices(treatPatches, patchVertex);

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();

    pointFieldPMG& points = mesh_.points();

    label nExtrudedVertices(0);
    forAll(patchVertex, bpI)
    {
        if (patchVertex[bpI])
        {
            ++nExtrudedVertices;
        }
    }

    points.setSize(points.size() + nExtrudedVertices);

    forAll(bPoints, bpI)
    {
        if (patchVertex[bpI])
        {
            points[nPoints_] = createNewVertex(bpI, treatPatches, patchVertex);
            newLabelForVertex_[bPoints[bpI]] = nPoints_;
            ++nPoints_;
        }
    }

    // swap coordinates of new and old points
    forAll(bPoints, bpI)
    {
        const label pLabel = newLabelForVertex_[bPoints[bpI]];
        if (pLabel != -1)
        {
            const point p = points[pLabel];
            points[pLabel] = points[bPoints[bpI]];
            points[bPoints[bpI]] = p;
        }
    }

    if (nPoints_ != points.size())
    {
        FatalErrorInFunction
            << "Number of vertices " << nPoints_
            << " does not match the list size "
            << abort(FatalError);
    }

    Info<< "Finished creating layer vertices" << endl;
}


void Foam::Module::boundaryLayers::createNewVertices
(
    const labelList& patchLabels
)
{
    otherVrts_.clear();

    patchKey_.setSize(mesh_.boundaries().size());
    patchKey_ = -1;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();

    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    // the following is needed for parallel runs
    // it is ugly, but must stay for now :(
    mse.boundaryFaces();
    mse.pointNormals();
    mse.pointFaces();
    mse.pointPoints();

    pointFieldPMG& points = mesh_.points();
    boolList treatPatches(mesh_.boundaries().size());
    List<direction> patchVertex(bPoints.size());

    // make sure than the points are never re-allocated during the process
    points.reserve(points.size() + 2*bPoints.size());

    // generate new layer vertices for each patch
    forAll(patchLabels, patchI)
    {
        const label pLabel = patchLabels[patchI];
        treatPatches = false;

        bool treat(true);
        forAll(treatPatchesWithPatch_[pLabel], pI)
        {
            const label otherPatch = treatPatchesWithPatch_[pLabel][pI];
            treatPatches[otherPatch] = true;

            if (patchKey_[otherPatch] == -1)
            {
                patchKey_[otherPatch] = patchI;
            }
            else
            {
                treat = false;
            }
        }

        if (!treat)
        {
            continue;
        }

        const label pKey = patchKey_[pLabel];

        // classify vertices belonging to this patch
        findPatchVertices(treatPatches, patchVertex);

        // create indices and allocate maps for new points
        labelLongList patchPoints;
        forAll(bPoints, bpI)
        {
            if (!patchVertex[bpI])
            {
                continue;
            }

            patchPoints.append(bpI);
            const label pointI = bPoints[bpI];

            if (patchVertex[bpI] & EDGENODE)
            {
                if (otherVrts_.find(pointI) == otherVrts_.end())
                {
                    std::map<std::pair<label, label>, label> m;

                    otherVrts_.insert(std::make_pair(pointI, m));
                }

                std::pair<label, label> pr(pKey, pKey);
                otherVrts_[pointI].insert(std::make_pair(pr, nPoints_++));
            }
            else
            {
                // this the only new point
                newLabelForVertex_[pointI] = nPoints_++;
            }
        }

        // set the size of points
        points.setSize(nPoints_);

        // calculate coordinates of new points
        # pragma omp parallel for schedule(dynamic, 50)
        forAll(patchPoints, i)
        {
            const label bpI = patchPoints[i];

            const label pointI = bPoints[bpI];

            // create new point
            const point p = createNewVertex(bpI, treatPatches, patchVertex);

            if (patchVertex[bpI] & EDGENODE)
            {
                // set the new point or an edge point
                if (otherVrts_.find(pointI) == otherVrts_.end())
                {
                    std::map<std::pair<label, label>, label> m;

                    otherVrts_.insert(std::make_pair(pointI, m));
                }

                std::pair<label, label> pr(pKey, pKey);
                const label npI = otherVrts_[pointI][pr];
                points[npI] = p;
            }
            else
            {
                // set the new point
                points[newLabelForVertex_[pointI]] = p;
            }
        }
    }

    // create missing vertices for edge and corner vertices
    // they should be stored in the otherNodes map
    forAll(bPoints, bpI)
    {
        const label pointI = bPoints[bpI];

        if (otherVrts_.find(pointI) == otherVrts_.end())
        {
            continue;
        }

        const point& p = points[pointI];

        std::map<std::pair<label, label>, label>& m = otherVrts_[pointI];
        DynList<label> usedPatches;
        DynList<label> newNodeLabel;
        DynList<vector> newPatchPenetrationVector;
        forAllRow(pPatches, bpI, patchI)
        {
            const label pKey = patchKey_[pPatches(bpI, patchI)];
            const std::pair<label, label> pr(pKey, pKey);
            const std::map<std::pair<label, label>, label>::const_iterator it =
                m.find(pr);
            if ((it != m.end()) && !usedPatches.found(pKey))
            {
                usedPatches.append(pKey);
                newNodeLabel.append(it->second);
                newPatchPenetrationVector.append(points[it->second] - p);
            }
        }

        if (newNodeLabel.size() == 1)
        {
            // only one patch is treated
            newLabelForVertex_[pointI] = newNodeLabel[0];
            otherVrts_.erase(pointI);
        }
        else if (newNodeLabel.size() == 2)
        {
            // point is located at an extrusion edge
            // create the new position for the existing point
            point newP(p);
            newP += newPatchPenetrationVector[0];
            newP += newPatchPenetrationVector[1];

            if (!help::isnan(newP) && !help::isinf(newP))
            {
                points.append(newP);
            }
            else
            {
                points.append(p);
            }
            newLabelForVertex_[pointI] = nPoints_;
            ++nPoints_;
        }
        else if (newNodeLabel.size() == 3)
        {
            // point is located at an extrusion corner
            // create 3 points and the new position for the existing point
            point newP(p);
            for (label i = 0; i < 3; ++i)
            {
                newP += newPatchPenetrationVector[i];
                for (label j = i + 1; j < 3; ++j)
                {
                    const point np =
                        p + newPatchPenetrationVector[i] +
                        newPatchPenetrationVector[j];

                    if (!help::isnan(np) && !help::isinf(np))
                    {
                        points.append(np);
                    }
                    else
                    {
                        points.append(p);
                    }

                    m.insert
                    (
                        std::make_pair
                        (
                            std::make_pair(usedPatches[i], usedPatches[j]),
                            nPoints_
                        )
                    );
                    ++nPoints_;
                }
            }

            // create new position for the existing point
            if (!help::isnan(newP) && !help::isinf(newP))
            {
                points.append(newP);
            }
            else
            {
                points.append(p);
            }

            newLabelForVertex_[pointI] = nPoints_;
            ++nPoints_;
        }
        else
        {
            FatalErrorInFunction
                << "Boundary node " << bpI << " is not at an edge!"
                << abort(FatalError);
        }
    }

    // swap coordinates of new and old points
    # pragma omp parallel for if (bPoints.size() > 1000) \
    schedule(dynamic, 100)
    forAll(bPoints, bpI)
    {
        const label pLabel = newLabelForVertex_[bPoints[bpI]];

        if (pLabel != -1)
        {
            const point p = points[pLabel];
            points[pLabel] = points[bPoints[bpI]];
            points[bPoints[bpI]] = p;
        }
    }
}


// ************************************************************************* //
