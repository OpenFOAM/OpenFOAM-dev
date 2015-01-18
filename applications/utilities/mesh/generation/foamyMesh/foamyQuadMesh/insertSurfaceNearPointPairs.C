/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

\*----------------------------------------------------------------------------*/

#include "CV2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV2D::insertSurfaceNearPointPairs()
{
    Info<< "insertSurfaceNearPointPairs: ";

    label nNearPoints = 0;

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        eit++
    )
    {
        Vertex_handle v0h = eit->first->vertex(cw(eit->second));
        Vertex_handle v1h = eit->first->vertex(ccw(eit->second));

        if (v0h->ppMaster() && v1h->ppMaster())
        {
            point2DFromPoint v0(toPoint2D(v0h->point()));
            point2DFromPoint v1(toPoint2D(v1h->point()));

            // Check that the two triangle vertices are further apart than the
            // minimum cell size
            if (magSqr(v1 - v0) > meshControls().minCellSize2())
            {
                point2D e0(toPoint2D(circumcenter(eit->first)));

                point2D e1
                (
                    toPoint2D(circumcenter(eit->first->neighbor(eit->second)))
                );

                // Calculate the length^2 of the edge normal to the surface
                scalar edgeLen2 = magSqr(e0 - e1);

                if (edgeLen2 < meshControls().minNearPointDist2())
                {
                    pointIndexHit pHit;
                    label hitSurface = -1;

                    qSurf_.findSurfaceNearest
                    (
                        toPoint3D(e0),
                        meshControls().minEdgeLen2(),
                        pHit,
                        hitSurface
                    );

                    if (pHit.hit())
                    {
                        vectorField norm(1);
                        qSurf_.geometry()[hitSurface].getNormal
                        (
                            List<pointIndexHit>(1, pHit),
                            norm
                        );

                        insertPointPair
                        (
                            meshControls().ppDist(),
                            toPoint2D(pHit.hitPoint()),
                            toPoint2D(norm[0])
                        );

                        nNearPoints++;

                        // Correct the edge iterator for the change in the
                        // number of edges following the point-pair insertion
                        eit = Finite_edges_iterator
                        (
                            finite_edges_end().base(),
                            eit.predicate(),
                            eit.base()
                        );
                    }
                }
            }
        }
    }

    Info<< nNearPoints << " point-pairs inserted" << endl;
}


// ************************************************************************* //
