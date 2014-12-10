/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::searchableSurfacesQueries

Description
    A collection of tools for searchableSurfaces.

SourceFiles
    searchableSurfacesQueries.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfacesQueries_H
#define searchableSurfacesQueries_H

#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class plane;

/*---------------------------------------------------------------------------*\
                           Class searchableSurfacesQueries Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfacesQueries
{
    // Private data

    // Private Member Functions

        //- Temporary wrapper around findNearest. Used in facesIntersection only
        static pointIndexHit tempFindNearest
        (
            const searchableSurface&,
            const point& pt,
            const scalar initDistSqr
        );

        //- Calculate sum of distances to nearest point on surfaces. Is used
        //  in minimisation to find intersection. Returns sum of (square of)
        //  distances to the surfaces.
        static scalar sumDistSqr
        (
            const PtrList<searchableSurface>&,
            const labelList& surfacesToTest,
            const scalar initialDistSqr,    // search box
            const point& pt
        );

        //- Takes the tet (points p) and reflects the point with the
        //  highest value around the centre (pSum). Checks if it gets closer
        //  and updates p, y if so.
        static scalar tryMorphTet
        (
            const PtrList<searchableSurface>&,
            const labelList& surfacesToTest,
            const scalar initialDistSqr,
            List<vector>& p,
            List<scalar>& y,
            vector& pSum,
            const label ihi,
            const scalar fac
        );

        //- Downhill simplex method: find the point with min cumulative
        //  distance to all surfaces. Does so by morphing a tet (points p).
        //  Returns the point on the 0th surface or hit if not reached within
        //  maxIters iterations.
        static bool morphTet
        (
            const PtrList<searchableSurface>&,
            const labelList& surfacesToTest,
            const scalar initialDistSqr,
            const scalar convergenceDistSqr,
            const label maxIter,
            List<vector>& p,
            List<scalar>& y
        );

        //static void findAllIntersections
        //(
        //    const searchableSurface& s,
        //    const pointField& start,
        //    const pointField& end,
        //    const vectorField& smallVec,
        //    List<List<pointIndexHit> >&
        //);

        static void mergeHits
        (
            const point& start,

            const label surfI,
            const List<pointIndexHit>& surfHits,

            labelList& allSurfaces,
            List<pointIndexHit>& allInfo,
            scalarList& allDistSqr
        );

public:

    // Declare name of the class and its debug switch
    ClassName("searchableSurfacesQueries");


        // Multiple point queries.

            //- Find any intersection. Return hit point information and
            //  index in surfacesToTest. If multiple surfaces hit the first
            //  surface is returned, not necessarily the nearest (to start).
            static void findAnyIntersection
            (
                const PtrList<searchableSurface>&,
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            );

            //- Find all intersections in order from start to end. Returns for
            //  every hit the index in surfacesToTest and the hit info.
            static void findAllIntersections
            (
                const PtrList<searchableSurface>&,
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,
                labelListList& surfaces,
                List<List<pointIndexHit> >& surfaceHits
            );

            //Find intersections of edge nearest to both endpoints.
            static void findNearestIntersection
            (
                const PtrList<searchableSurface>& allSurfaces,
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,
                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& surface2,
                List<pointIndexHit>& hit2
            );

            //- Find nearest. Return -1 (and a miss()) or surface and nearest
            //  point.
            static void findNearest
            (
                const PtrList<searchableSurface>&,
                const labelList& surfacesToTest,
                const pointField&,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            );

            //- Find nearest points to a specific region of the surface
            static void findNearest
            (
                const PtrList<searchableSurface>& allSurfaces,
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                const labelList& regionIndices,
                labelList& nearestSurfaces,
                List<pointIndexHit>& nearestInfo
            );

            //- Find signed distance to nearest surface. Outside is positive.
            //  illegalHandling: how to handle non-inside or outside
            //      OUTSIDE : treat as outside
            //      INSIDE  : treat as inside
            //      UNKNOWN : throw fatal error
            static void signedDistance
            (
                const PtrList<searchableSurface>& allSurfaces,
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                const volumeType illegalHandling,
                labelList& nearestSurfaces,
                scalarField& distance
            );

            //- Find the boundBox of the selected surfaces
            static boundBox bounds
            (
                const PtrList<searchableSurface>& allSurfaces,
                const labelList& surfacesToTest
            );

        // Single point queries

            //- Calculate point which is on a set of surfaces. WIP.
            static pointIndexHit facesIntersection
            (
                const PtrList<searchableSurface>& allSurfaces,
                const labelList& surfacesToTest,
                const scalar initDistSqr,
                const scalar convergenceDistSqr,
                const point& start
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
