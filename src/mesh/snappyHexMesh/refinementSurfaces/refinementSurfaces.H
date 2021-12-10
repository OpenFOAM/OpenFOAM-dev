/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::refinementSurfaces

Description
    Container for data on surfaces used for surface-driven refinement.
    Contains all the data about the level of refinement needed per
    surface.

SourceFiles
    refinementSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef refinementSurfaces_H
#define refinementSurfaces_H

#include "triSurfaceGeoMesh.H"
#include "triSurfaceFields.H"
#include "vectorList.H"
#include "pointIndexHit.H"
#include "surfaceZonesInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;
class refinementRegions;
class triSurfaceMesh;

typedef List<point> pointList;

/*---------------------------------------------------------------------------*\
                     Class refinementSurfaces Declaration
\*---------------------------------------------------------------------------*/

class refinementSurfaces
{
    // Private Data

        //- Reference to all geometry.
        const searchableSurfaces& allGeometry_;

        //- Indices of surfaces that are refinement ones
        labelList surfaces_;

        //- Surface name (word)
        wordList names_;

        //- List of surface zone (face and cell zone) information
        PtrList<surfaceZonesInfo> surfZones_;

        //- From local region number to global region number
        labelList regionOffset_;

        //- From global region number to refinement level
        labelList minLevel_;

        //- From global region number to refinement level
        labelList maxLevel_;

        //- From global region number to small-gap level
        labelList gapLevel_;

        //- From global region number to perpendicular angle
        scalarField perpendicularAngle_;

        //- From global region number to patchType
        PtrList<dictionary> patchInfo_;


public:

    // Constructors

        //- Construct from surfaces and dictionary
        refinementSurfaces
        (
            const searchableSurfaces& allGeometry,
            const dictionary&,
            const label gapLevelIncrement
        );

        //- Construct from components
        refinementSurfaces
        (
            const searchableSurfaces& allGeometry,
            const labelList& surfaces,
            const wordList& names,
            const PtrList<surfaceZonesInfo>& surfZones,
            const labelList& regionOffset,
            const labelList& minLevel,
            const labelList& maxLevel,
            const labelList& gapLevel,
            const scalarField& perpendicularAngle,
            PtrList<dictionary>& patchInfo
        );

        //- Disallow default bitwise copy construction
        refinementSurfaces(const refinementSurfaces&) = delete;


    // Member Functions

        // Access

            const searchableSurfaces& geometry() const
            {
                return allGeometry_;
            }

            const labelList& surfaces() const
            {
                return surfaces_;
            }

            //- Names of surfaces
            const wordList& names() const
            {
                return names_;
            }

            const PtrList<surfaceZonesInfo>& surfZones() const
            {
                return surfZones_;
            }

            //- From local region number to global region number
            const labelList& regionOffset() const
            {
                return regionOffset_;
            }

            //- From global region number to refinement level
            const labelList& minLevel() const
            {
                return minLevel_;
            }

            //- From global region number to refinement level
            const labelList& maxLevel() const
            {
                return maxLevel_;
            }

            //- From global region number to small gap refinement level
            const labelList& gapLevel() const
            {
                return gapLevel_;
            }

            //- From global region number to perpendicular angle
            const scalarField& perpendicularAngle() const
            {
                return perpendicularAngle_;
            }

            //- From global region number to patch type
            const PtrList<dictionary>& patchInfo() const
            {
                return patchInfo_;
            }


        // Helper

            //- From surface and region on surface to global region
            label globalRegion(const label surfi, const label regioni) const
            {
                return regionOffset_[surfi]+regioni;
            }

            //- Min level for surface and region on surface
            label minLevel(const label surfi, const label regioni) const
            {
                return minLevel_[globalRegion(surfi, regioni)];
            }

            //- Max level for surface and region on surface
            label maxLevel(const label surfi, const label regioni) const
            {
                return maxLevel_[globalRegion(surfi, regioni)];
            }

            label nRegions() const
            {
                return minLevel_.size();
            }

            //- Calculate the refinement level for every element
            //  of the searchablesurface
            void setMinLevelFields
            (
                const refinementRegions& shells,
                const scalar level0EdgeLength,
                const bool extendedRefinementSpan
            );


        // Searching

            //- Find intersection of edge. Return -1 or first surface
            //  with higher (than currentLevel) minlevel.
            //  Return surface number and level.
            void findHigherIntersection
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                labelList& surfaces,
                labelList& surfaceLevel
            ) const;

            //- Find all intersections of edge. Unsorted order.
            void findAllHigherIntersections
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level
                const labelList& globalRegionLevel, // level per surfregion

                List<vectorList>& surfaceNormal,
                labelListList& surfaceLevel
            ) const;

            //- Find all intersections of edge. Unsorted order.
            void findAllHigherIntersections
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level
                const labelList& globalRegionLevel, // level per surfregion

                List<pointList>& surfaceLocation,
                List<vectorList>& surfaceNormal,
                labelListList& surfaceLevel
            ) const;

            //- Find intersection nearest to the endpoints. surface1,2 are
            //  not indices into surfacesToTest but refinement surface indices.
            //  Returns surface, region on surface (so not global surface)
            //  and position on surface.
            void findNearestIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2
            ) const;

            //- findNearestIntersection but also get normals
            void findNearestIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                vectorField& normal1,

                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2,
                vectorField& normal2
            ) const;

            //- Used for debugging only: find intersection of edge.
            void findAnyIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find nearest point on surfaces.
            void findNearest
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find nearest point on surfaces. Return surface and region on
            //  surface (so not global surface)
            void findNearestRegion
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& hitSurface,
                labelList& hitRegion
            ) const;

            //- Find nearest point on surfaces. Return surface, region and
            //  normal on surface (so not global surface)
            void findNearestRegion
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& hitSurface,
                List<pointIndexHit>& hitinfo,
                labelList& hitRegion,
                vectorField& hitNormal
            ) const;

            //- Detect if a point is 'inside' (closed) surfaces.
            //  Returns -1 if not, returns first surface it is.
            void findInside
            (
                const labelList& surfacesToTest,
                const pointField& pt,
                labelList& insideSurfaces
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const refinementSurfaces&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
