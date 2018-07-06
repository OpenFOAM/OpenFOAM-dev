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
    Foam::line

Description
    A line primitive.

SourceFiles
    lineI.H

\*---------------------------------------------------------------------------*/

#ifndef line_H
#define line_H

#include "vector.H"
#include "PointHit.H"
#include "point2D.H"
#include "FixedList.H"
#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class Istream;
class Ostream;


// Forward declaration of friend functions and operators

template<class Point, class PointRef> class line;

template<class Point, class PointRef>
inline Istream& operator>>(Istream&, line<Point, PointRef>&);

template<class Point, class PointRef>
inline Ostream& operator<<(Ostream&, const line<Point, PointRef>&);


/*---------------------------------------------------------------------------*\
                           Class line Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class line
{
    // Private data

        PointRef a_, b_;


public:

    // Constructors

        //- Construct from two points
        inline line(const Point& start, const Point& end);

        //- Construct from two points in the list of points
        //  The indices could be from edge etc.
        inline line
        (
            const UList<Point>&,
            const FixedList<label, 2>& indices
        );

        //- Construct from Istream
        inline line(Istream&);


    // Member functions

        // Access

            //- Return first vertex
            inline PointRef start() const;

            //- Return second vertex
            inline PointRef end() const;


        // Properties

            //- Return centre (centroid)
            inline Point centre() const;

            //- Return scalar magnitude
            inline scalar mag() const;

            //- Return start-end vector
            inline Point vec() const;

            //- Return nearest distance to line from a given point
            //  If the nearest point is on the line, return a hit
            PointHit<Point> nearestDist(const Point& p) const;

            //- Return nearest distance from line to line. Returns distance
            //  and sets both points (one on *this, one on the provided
            //  linePointRef.
            scalar nearestDist
            (
                const line<Point, const Point&>& edge,
                Point& thisPoint,
                Point& edgePoint
            ) const;

            bool insideBoundBox(const Point&) const;


    // Ostream operator

        friend Istream& operator>> <Point, PointRef>
        (
            Istream&,
            line&
        );

        friend Ostream& operator<< <Point, PointRef>
        (
            Ostream&,
            const line&
        );
};


//- 2D specialisation
template<>
scalar line<point2D, const point2D&>::nearestDist
(
    const line<point2D, const point2D&>& edge,
    point2D& thisPoint,
    point2D& edgePoint
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
