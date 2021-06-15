/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::polygonTriangulate

Description
    Triangulation of three-dimensional polygons

SourceFiles
    polygonTriangulate.C

\*---------------------------------------------------------------------------*/

#ifndef polygonTriangulate_H
#define polygonTriangulate_H

#include "DynamicList.H"
#include "HashSet.H"
#include "labelPair.H"
#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class polygonTriangulate Declaration
\*---------------------------------------------------------------------------*/

class polygonTriangulate
{
    // Private Data

        //- Indirect addressing into point field to form the current as yet
        //  un-triangulated polygon
        DynamicList<label> pointis_;

        //- Edge indices of the un-triangulated polygon
        DynamicList<label> edges_;

        //- The interior angles of the vertices of the un-triangulated polygon
        DynamicList<scalar> angle_;

        //- Flags to denote whether points of the un-triangulated polygon are
        //  "ears"; i.e., that the triangle formed by adding a diagonal between
        //  two adjacent points contains no other point
        DynamicList<bool> ear_;

        //- Set of tri-faces used to mark previously considered configurations
        //  in optimisation and prevent unnecessary repetition
        HashSet<triFace, Hash<triFace>> triPointsSet_;

        //- Point array. Used if the input points cannot be SubList-d.
        DynamicList<point> points_;

        //- Tri-points
        DynamicList<triFace> triPoints_;

        //- Tri-edge associations
        DynamicList<FixedList<label, 3>> triEdges_;

        //- Edge-tri associations
        DynamicList<labelPair> edgeTris_;


    // Private Static Member Functions

        // Renumber a polygon label to global indexing
        static inline label renumberPolyToGlobal
        (
            const label triPoly,
            const UList<label>& polyToGlobal
        );

        // Renumber a container of polygon labels to global indexing
        template<class Type>
        static inline Type renumberPolyToGlobal
        (
            const Type& triPoly,
            const UList<label>& polyToGlobal
        );

        //- Return the area of a triangle projected in a normal direction
        template<class PointField>
        static scalar area
        (
            const triFace& triPoints,
            const PointField& points,
            const vector& normal
        );

        //- Return the quality of a triangle projected in a normal direction
        template<class PointField>
        static scalar quality
        (
            const triFace& triPoints,
            const PointField& points,
            const vector& normal
        );

        //- Return whether or not two edges intersect in the plane defined by a
        //  normal
        template<class PointField>
        static bool intersection
        (
            const edge& edgePointsA,
            const edge& edgePointsB,
            const PointField& points,
            const vector& normal
        );

        //- Return the number of intersections of an edge within a polygon
        template<class PointField>
        static label nIntersections
        (
            const edge& edgePoints,
            const PointField& points,
            const vector& normal
        );

        //- Return the interior angle of the point in a polygon in the plane
        //  defined by a normal
        template<class PointField>
        static scalar angle
        (
            const label pointi,
            const PointField& points,
            const vector& normal
        );

        //- Return whether a point in a polygon is an "ear"; i.e., that the
        //  triangle formed by adding a diagonal between the two adjacent
        //  points contains no other point
        template<class PointField>
        static bool ear
        (
            const label pointi,
            const PointField& points,
            const vector& normal
        );


    // Private Member Functions

        //- Optimise the triangulation quality by flipping edges
        template<class PointField>
        void optimiseTriangulation
        (
            const label trii,
            const PointField& points,
            const vector& normal,
            UList<triFace>& triPoints,
            UList<FixedList<label, 3>>& triEdges,
            UList<labelPair>& edgeTris
        );

        //- Perform triangulation of a polygon by ear clipping, assuming that
        //  the polygon is simple/non-intersecting
        void simpleTriangulate
        (
            const UList<point>& points,
            const vector& normal,
            UList<triFace>& triPoints,
            UList<FixedList<label, 3>>& triEdges,
            UList<labelPair>& edgeTris,
            const bool optimal
        );

        //- Perform triangulation of a polygon, first by adding a spanning
        //  triangle between a point and an edge, then recursing into the
        //  sub-polygons on either side
        void partitionTriangulate
        (
            const UList<point>& points,
            const vector& normal,
            const label spanPointi,
            const label spanEdgei,
            UList<triFace>& triPoints,
            UList<FixedList<label, 3>>& triEdges,
            UList<labelPair>& edgeTris,
            const bool simple,
            const bool optimal
        );

        //- Perform triangulation of the given polygon, detecting self
        //  intersections and resolving them by adding spanning triangles prior
        //  to simple triangulation
        void complexTriangulate
        (
            const UList<point>& points,
            const vector& normal,
            UList<triFace>& triPoints,
            UList<FixedList<label, 3>>& triEdges,
            UList<labelPair>& edgeTris,
            const bool optimal
        );

        //- Perform triangulation of the given polygon
        void triangulate
        (
            const UList<point>& points,
            const vector& normal,
            UList<triFace>& triPoints,
            UList<FixedList<label, 3>>& triEdges,
            UList<labelPair>& edgeTris,
            const bool simple = false,
            const bool optimal = true
        );


public:

    // Static Member Functions

        //- Generate a random polygon for testing purposes
        static List<point> randomPolygon
        (
            Random& rndGen,
            const label n,
            const scalar error
        );



    // Constructors

        //- Construct null
        polygonTriangulate();

        //- Disallow default bitwise copy construction
        polygonTriangulate(const polygonTriangulate&) = delete;


    //- Destructor
    ~polygonTriangulate();


    // Member Functions

        // Access

            //- Get the triangles' points
            inline const UList<triFace>& triPoints() const;

            //- Get the triangles' renumbered points
            inline List<triFace> triPoints
            (
                const UList<label>& polyPoints
            ) const;

            //- Get a triangle's renumbered points
            inline triFace triPoints
            (
                const label trii,
                const UList<label>& polyPoints
            ) const;

            //- Get the triangles' edges
            inline const UList<FixedList<label, 3>>& triEdges() const;

            //- Get the triangles' renumbered edges
            inline List<FixedList<label, 3>> triEdges
            (
                const UList<label>& polyEdges
            ) const;

            //- Get a triangle's renumbered edges
            inline FixedList<label, 3> triEdges
            (
                const label trii,
                const UList<label>& polyEdges
            ) const;


        //- Perform triangulation of the given polygon
        const UList<triFace>& triangulate
        (
            const UList<point>& points,
            const vector& normal,
            const bool simple = false,
            const bool optimal = true
        );

        //- Perform triangulation of the given polygon
        const UList<triFace>& triangulate
        (
            const UList<point>& points,
            const bool simple = false,
            const bool optimal = true
        );

        //- Perform triangulation of the given polygon
        inline const UList<triFace>& triangulate
        (
            const UIndirectList<point>& points,
            const vector& normal,
            const bool simple = false,
            const bool optimal = true
        );

        //- Perform triangulation of the given polygon
        inline const UList<triFace>& triangulate
        (
            const UIndirectList<point>& points,
            const bool simple = false,
            const bool optimal = true
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polygonTriangulate&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polygonTriangulateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
