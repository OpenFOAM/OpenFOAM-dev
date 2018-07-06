/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
class pointConstraint;

/*---------------------------------------------------------------------------*\
                           Class searchableSurfacesQueries Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfacesQueries
{
    // Private Member Functions

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
                List<List<pointIndexHit>>& surfaceHits
            );

            //- Find intersections of edge nearest to both endpoints.
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

            //- Find nearest points that are on all supplied surfaces
            //  (nearest point if single surface; nearest intersection by
            //   steepst descent if on multiple surfaces). Returns current
            //   best guess). Wip.
            static void findNearest
            (
                const PtrList<searchableSurface>& allSurfaces,
                const labelList& surfacesToTest,
                const pointField& start,
                const scalarField& distSqr,
                pointField& near,
                List<pointConstraint>& constraint,
                const label nIter = 20
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
