/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::searchableSurfaces

Description
    Container for searchableSurfaces.

SourceFiles
    searchableSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaces_H
#define searchableSurfaces_H

#include "searchableSurface.H"
#include "labelPair.H"
#include "writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class searchableSurfaces Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaces
:
    public PtrList<searchableSurface>
{
    // Private data

        //- Surface names
        wordList names_;

        //- Region names per surface
        List<wordList> regionNames_;

        ////- From global region name to surface and region on surface
        // HashTable<labelPair> regionNames_;

        //- Indices of all surfaces. Precalculated and stored.
        labelList allSurfaces_;


    // Private Member Functions

        //- Is edge on face
        static bool connected
        (
            const triSurface& s,
            const label edgeI,
            const pointIndexHit& hit
        );


        //- Disallow default bitwise copy construct
        searchableSurfaces(const searchableSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const searchableSurfaces&);


public:

    ClassName("searchableSurfaces");


    // Constructors

        //- Construct with length specified. Fill later.
        explicit searchableSurfaces(const label);


        ////- Construct from list of dictionaries
        // searchableSurfaces(const IOobject&, const PtrList<dictionary>&);

        //- Construct from dictionary and whether to construct names always
        //  as surfaceName "_" regionName (singleRegionName false) or
        //  for single region surfaces as surfaceName only (singleRegionName
        //  true)
        searchableSurfaces
        (
            const IOobject&,
            const dictionary&,
            const bool singleRegionName
        );


    // Member Functions

        const wordList& names() const
        {
            return names_;
        }

        wordList& names()
        {
            return names_;
        }

        const List<wordList>& regionNames() const
        {
            return regionNames_;
        }

        List<wordList>& regionNames()
        {
            return regionNames_;
        }


        ////- If adding surfaces 'by hand'
        // HashTable<labelPair>& regionNames()
        //{
        //    return regionNames_;
        //}
        ////- Get surface and region for a name
        // const labelPair& surfaceRegion(const word& globalRegion) const
        //{
        //    return regionNames_[globalRegion];
        //}

        //- Find index of surface. Return -1 if not found.
        label findSurfaceID(const word& name) const;

        label findSurfaceRegionID
        (
            const word& surfaceName,
            const word& regionName
        ) const;

        // Multiple point queries.

            //- Find any intersection. Return hit point information and
            //  surface number. If multiple surfaces hit the first surface
            //  is returned, not necessarily the nearest (to start).
            void findAnyIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find all intersections in order from start to end. Returns for
            //  every hit the surface and the hit info.
            void findAllIntersections
            (
                const pointField& start,
                const pointField& end,
                labelListList& surfaces,
                List<List<pointIndexHit>>&
            ) const;

            // Find intersections of edge nearest to both endpoints.
            void findNearestIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& surface2,
                List<pointIndexHit>& hit2
            ) const;

            //- Find nearest. Return -1 (and a miss()) or surface and nearest
            //  point.
            void findNearest
            (
                const pointField&,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            void findNearest
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                const labelList& regionIndices,
                labelList& nearestSurfaces,
                List<pointIndexHit>& nearestInfo
            ) const;

            //- Calculate bounding box
            boundBox bounds() const;

        // Checking

            //- Are all surfaces closed and manifold
            bool checkClosed(const bool report) const;

            //- Are all (triangulated) surfaces consistent normal orientation
            bool checkNormalOrientation(const bool report) const;

            //- Are all bounding boxes of similar size
            bool checkSizes(const scalar maxRatio, const bool report) const;

            //- Do surfaces self-intersect or intersect others
            bool checkIntersection
            (
                const scalar tol,
                const autoPtr<writer<scalar>>&,
                const bool report
            ) const;

            //- Check triangle quality
            bool checkQuality
            (
                const scalar minQuality,
                const bool report
            ) const;

            //- All topological checks. Return number of failed checks
            label checkTopology(const bool report) const;

            //- All geometric checks. Return number of failed checks
            label checkGeometry
            (
                const scalar maxRatio,
                const scalar tolerance,
                const autoPtr<writer<scalar>>& setWriter,
                const scalar minQuality,
                const bool report
            ) const;

            //- Write some stats
            void writeStats(const List<wordList>&, Ostream&) const;


    // Member Operators

        //- Return const and non-const reference to searchableSurface by index.
        using PtrList<searchableSurface>::operator[];

        //- Return const reference to searchableSurface by name.
        const searchableSurface& operator[](const word&) const;

        //- Return reference to searchableSurface by name.
        searchableSurface& operator[](const word&);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
