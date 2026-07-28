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
#include "partTriMesh.H"
#include "meshSurfacePartitioner.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"

#include <map>

#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Module::partTriMesh::partTriMesh(const meshSurfacePartitioner& mPart)
:
    mPart_(mPart),
    surf_(),
    pointLabelInMeshSurface_(),
    meshSurfacePointLabelInTriMesh_(),
    pointType_()
{
    const meshSurfaceEngine& meshSurface = mPart.surfaceEngine();
    List<direction> useFace(meshSurface.boundaryFaces().size(), direction(1));

    createPointsAndTrias(useFace);
}


Foam::Module::partTriMesh::partTriMesh
(
    const meshSurfacePartitioner& mPart,
    const labelHashSet& invertedPoints,
    const label additionalLayers
)
:
    mPart_(mPart),
    surf_(),
    pointLabelInMeshSurface_(),
    meshSurfacePointLabelInTriMesh_(),
    pointType_()
{
    const meshSurfaceEngine& meshSurface = mPart.surfaceEngine();
    const labelListList2D& pointFaces = meshSurface.pointFaces();
    const faceList::subList& bFaces = meshSurface.boundaryFaces();
    const labelList& bPoints = meshSurface.boundaryPoints();
    const labelList& bp = meshSurface.bp();

    List<direction> useFace(bFaces.size(), direction(0));

    // select cells containing at least one vertex of the bad faces
    forAll(pointFaces, bpI)
    {
        if (invertedPoints.found(bPoints[bpI]))
        {
            forAllRow(pointFaces, bpI, pfI)
            {
                useFace[pointFaces[bpI][pfI]] = 1;
            }
        }
    }

    // add additional layer of cells
    for
    (
        direction layerI = 1;
        layerI<(additionalLayers + direction(1));
        ++layerI
    )
    {
        forAll(useFace, bfI)
        {
            if (useFace[bfI] == layerI)
            {
                const face& bf = bFaces[bfI];

                forAll(bf, pI)
                {
                    const label bpI = bp[bf[pI]];

                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label fLabel = pointFaces[bpI][pfI];

                        if (!useFace[fLabel])
                        {
                            useFace[fLabel] = layerI + 1;
                        }
                    }
                }
            }
        }
    }

    createPointsAndTrias(useFace);
}


Foam::Module::partTriMesh::~partTriMesh()
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::Module::partTriMesh::updateVertex(const label pI, const point& newP)
{
    triSurfModifier sMod(surf_);
    pointField& pts = sMod.pointsAccess();
    const labelListList2D& pointFacets = surf_.pointFacets();

    pts[pI] = newP;

    if (pointType_[pI] & FACECENTRE)
    {
        Warning << "Smoothing auxiliary vertex."
            << " This has no effect on the original mesh" << endl;
        return;
    }

    // find face centres attached
    DynamicList<label> helper;
    forAllRow(pointFacets, pI, ptI)
    {
        const label centreI = surf_[pointFacets[pI][ptI]][2];
        if (pointType_[centreI] & FACECENTRE)
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
        forAllRow(pointFacets, centreI, ptI)
        {
            const labelledTri& tri = surf_[pointFacets[centreI][ptI]];
            point c(vector::zero);
            for (label i = 0; i < 3; ++i)
            {
                c += pts[tri[i]];
            }

            c /= 3;
            const scalar area = tri.mag(pts) + VSMALL;

            centre += c*area;
            faceArea += area;
        }

        pts[centreI] = centre/faceArea;
    }
}


void Foam::Module::partTriMesh::updateVerticesSMP
(
    const List<DynamicList<labelledPoint>>& np
)
{
    triSurfModifier sMod(surf_);
    pointField& pts = sMod.pointsAccess();
    const labelListList2D& pointFacets = surf_.pointFacets();

    List<direction> updateType(pts.size(), direction(0));

    # pragma omp parallel num_threads(np.size())
    {
        const DynamicList<labelledPoint>& newPoints = np[omp_get_thread_num()];

        forAll(newPoints, i)
        {
            const labelledPoint& lp = newPoints[i];
            const label pointI = lp.pointLabel();

            pts[pointI] = lp.coordinates();
            updateType[pointI] |= SMOOTH;

            forAllRow(pointFacets, pointI, ptI)
            {
                const labelledTri& tri = surf_[pointFacets[pointI][ptI]];

                if (pointType_[tri[2]] & FACECENTRE)
                {
                    updateType[tri[2]] |= FACECENTRE;
                }
            }
        }
    }


    // update coordinates of FACECENTRE vertices
    # pragma omp parallel for schedule(dynamic, 20)
    forAll(updateType, pI)
    {
        if (updateType[pI] & FACECENTRE)
        {
            point centre(vector::zero);
            scalar faceArea(0.0);
            forAllRow(pointFacets, pI, ptI)
            {
                const labelledTri& tri = surf_[pointFacets[pI][ptI]];
                point c(vector::zero);
                for (label i = 0; i < 3; ++i)
                {
                    c += pts[tri[i]];
                }
                c /= 3;
                const scalar area = tri.mag(pts) + VSMALL;

                centre += c*area;
                faceArea += area;
            }

            pts[pI] = centre/faceArea;
        }
    }

}


void Foam::Module::partTriMesh::updateVertices()
{
    const meshSurfaceEngine& mse = mPart_.surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    DynamicList<label> movedPoints;
    movedPoints.setSize(bPoints.size());
    forAll(movedPoints, i)
    {
        movedPoints[i] = i;
    }

    updateVertices(movedPoints);
}


void Foam::Module::partTriMesh::updateVertices(const DynamicList<label>& movedPoints)
{
    const meshSurfaceEngine& mse = mPart_.surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelList& bPoints = mse.boundaryPoints();

    triSurfModifier sMod(surf_);
    pointField& pts = sMod.pointsAccess();
    const labelListList2D& pointFacets = surf_.pointFacets();

    List<direction> updateType(pts.size(), direction(0));

    // update coordinates of vertices which exist in the surface
    // of the volume mesh
    # pragma omp parallel for schedule(dynamic, 50)
    forAll(movedPoints, i)
    {
        const label bpI = movedPoints[i];
        const label pointI = bPoints[bpI];
        const label triPointI = meshSurfacePointLabelInTriMesh_[bpI];

        if (triPointI < 0)
        {
            continue;
        }

        pts[triPointI] = points[pointI];
        updateType[triPointI] |= SMOOTH;

        forAllRow(pointFacets, triPointI, ptI)
        {
            const labelledTri& tri = surf_[pointFacets[triPointI][ptI]];

            if (pointType_[tri[2]] & FACECENTRE)
            {
                updateType[tri[2]] |= FACECENTRE;
            }
        }
    }


    // update coordinates of FACECENTRE vertices
    # pragma omp parallel for schedule(dynamic, 20)
    forAll(updateType, pI)
    {
        if (updateType[pI] & FACECENTRE)
        {
            point centre(vector::zero);
            scalar faceArea(0.0);
            forAllRow(pointFacets, pI, ptI)
            {
                const labelledTri& tri = surf_[pointFacets[pI][ptI]];
                point c(vector::zero);
                for (label i = 0; i < 3; ++i)
                {
                    c += pts[tri[i]];
                }

                c /= 3;
                const scalar area = tri.mag(pts) + VSMALL;

                centre += c*area;
                faceArea += area;
            }

            pts[pI] = centre/faceArea;
        }
    }

}


// ************************************************************************* //
