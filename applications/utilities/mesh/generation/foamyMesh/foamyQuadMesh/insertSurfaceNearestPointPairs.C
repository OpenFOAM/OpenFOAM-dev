/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "CV2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::CV2D::dualCellSurfaceIntersection
(
    const Triangulation::Finite_vertices_iterator& vit
) const
{
    Triangulation::Edge_circulator ecStart = incident_edges(vit);
    Triangulation::Edge_circulator ec = ecStart;

    do
    {
        if (!is_infinite(ec))
        {
            point e0 = toPoint3D(circumcenter(ec->first));

            // If edge end is outside bounding box then edge cuts boundary
            if (!qSurf_.globalBounds().contains(e0))
            {
                return true;
            }

            point e1 = toPoint3D(circumcenter(ec->first->neighbor(ec->second)));

            // If other edge end is outside bounding box then edge cuts boundary
            if (!qSurf_.globalBounds().contains(e1))
            {
                return true;
            }

            if (magSqr(e1 - e0) > meshControls().minEdgeLen2())
            {
                if (qSurf_.findSurfaceAnyIntersection(e0, e1))
                {
                    return true;
                }
            }
        }

    } while (++ec != ecStart);

    return false;
}


void Foam::CV2D::insertPointPairs
(
    const DynamicList<point2D>& nearSurfacePoints,
    const DynamicList<point2D>& surfacePoints,
    const DynamicList<label>& surfaceTris,
    const DynamicList<label>& surfaceHits,
    const fileName fName
)
{
    if (meshControls().mirrorPoints())
    {
        forAll(surfacePoints, ppi)
        {
            insertMirrorPoint
            (
                nearSurfacePoints[ppi],
                surfacePoints[ppi]
            );
        }
    }
    else
    {
        forAll(surfacePoints, ppi)
        {
            pointIndexHit pHit
            (
                true,
                toPoint3D(surfacePoints[ppi]),
                surfaceTris[ppi]
            );

            vectorField norm(1);
            qSurf_.geometry()[surfaceHits[ppi]].getNormal
            (
                List<pointIndexHit>(1, pHit),
                norm
            );

            insertPointPair
            (
                meshControls().ppDist(),
                surfacePoints[ppi],
                toPoint2D(norm[0])
            );
        }
    }

    Info<< surfacePoints.size() << " point-pairs inserted" << endl;

    if (meshControls().objOutput())
    {
        OFstream str(fName);
        label vertI = 0;

        forAll(surfacePoints, ppi)
        {
            meshTools::writeOBJ(str, toPoint3D(surfacePoints[ppi]));
            vertI++;
        }

        Info<< "insertPointPairs: Written " << surfacePoints.size()
            << " inserted point-pair locations to file "
            << str.name() << endl;
    }
}


void Foam::CV2D::insertSurfaceNearestPointPairs()
{
    Info<< "insertSurfaceNearestPointPairs: ";

    label nSurfacePointsEst =
        min
        (
            number_of_vertices(),
            size_t(10*sqrt(scalar(number_of_vertices())))
        );

    DynamicList<point2D> nearSurfacePoints(nSurfacePointsEst);
    DynamicList<point2D> surfacePoints(nSurfacePointsEst);
    DynamicList<label> surfaceTris(nSurfacePointsEst);
    DynamicList<label> surfaceHits(nSurfacePointsEst);

    // Local references to surface mesh addressing
//    const pointField& localPoints = qSurf_.localPoints();
//    const labelListList& edgeFaces = qSurf_.edgeFaces();
//    const vectorField& faceNormals = qSurf_.faceNormals();
//    const labelListList& faceEdges = qSurf_.faceEdges();

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalPoint())
        {
            point2DFromPoint vert(toPoint2D(vit->point()));

            pointIndexHit pHit;
            label hitSurface = -1;

            qSurf_.findSurfaceNearest
            (
                toPoint3D(vert),
                4*meshControls().minCellSize2(),
                pHit,
                hitSurface
            );

            if (pHit.hit())
            {
                vit->setNearBoundary();

                // Reference to the nearest triangle
//                const labelledTri& f = qSurf_[hitSurface];

//                // Find where point is on triangle.
//                // Note tolerance needed is relative one
//                // (used in comparing normalized [0..1] triangle coordinates).
//                label nearType, nearLabel;
//                triPointRef
//                (
//                    localPoints[f[0]],
//                    localPoints[f[1]],
//                    localPoints[f[2]]
//                ).classify(pHit.hitPoint(), nearType, nearLabel);

//                // If point is on a edge check if it is an internal feature

//                bool internalFeatureEdge = false;

//                if (nearType == triPointRef::EDGE)
//                {
//                    label edgeI = faceEdges[hitSurface][nearLabel];
//                    const labelList& eFaces = edgeFaces[edgeI];

//                    if
//                    (
//                        eFaces.size() == 2
//                     && (faceNormals[eFaces[0]] & faceNormals[eFaces[1]])
//                       < -0.2
//                    )
//                    {
//                        internalFeatureEdge = true;
//                    }
//                }

                if (dualCellSurfaceIntersection(vit)) //&& !internalFeatureEdge)
                {
                    nearSurfacePoints.append(vert);
                    surfacePoints.append(toPoint2D(pHit.hitPoint()));
                    surfaceTris.append(pHit.index());
                    surfaceHits.append(hitSurface);
                }
            }
        }
    }

    insertPointPairs
    (
        nearSurfacePoints,
        surfacePoints,
        surfaceTris,
        surfaceHits,
        "surfaceNearestIntersections.obj"
    );
}


// ************************************************************************* //
