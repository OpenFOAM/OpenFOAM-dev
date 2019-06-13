/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::boundBox

Description
    A bounding box defined in terms of the points at its extremities.

\*---------------------------------------------------------------------------*/

#ifndef boundBox_H
#define boundBox_H

#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class boundBox;
template<class T> class tmp;

bool operator==(const boundBox&, const boundBox&);
bool operator!=(const boundBox&, const boundBox&);

Istream& operator>>(Istream&, boundBox&);
Ostream& operator<<(Ostream&, const boundBox&);


/*---------------------------------------------------------------------------*\
                          Class boundBox Declaration
\*---------------------------------------------------------------------------*/

class boundBox
{
    // Private Data

        //- Minimum and maximum describing the bounding box
        point min_, max_;

    // Private Member Functions

        //- Calculate the bounding box from the given points.
        //  Does parallel communication (doReduce = true)
        void calculate(const UList<point>&, const bool doReduce = true);

public:

    // Static Data Members

        //- The great value used for greatBox and invertedBox
        static const scalar great;

        //- A very large boundBox: min/max == -/+ vGreat
        static const boundBox greatBox;

        //- A very large inverted boundBox: min/max == +/- vGreat
        static const boundBox invertedBox;


    // Constructors

        //- Construct null, setting points to zero
        inline boundBox();

        //- Construct from components
        inline boundBox(const point& min, const point& max);

        //- Construct as the bounding box of the given points
        //  Does parallel communication (doReduce = true)
        boundBox(const UList<point>&, const bool doReduce = true);

        //- Construct as the bounding box of the given temporary pointField.
        //  Does parallel communication (doReduce = true)
        boundBox(const tmp<pointField>&, const bool doReduce = true);

        //- Construct bounding box as subset of the pointField.
        //  The indices could be from cell/face etc.
        //  Does parallel communication (doReduce = true)
        boundBox
        (
            const UList<point>&,
            const labelUList& indices,
            const bool doReduce = true
        );

        //- Construct bounding box as subset of the pointField.
        //  The indices could be from edge/triFace etc.
        //  Does parallel communication (doReduce = true)
        template<unsigned Size>
        boundBox
        (
            const UList<point>&,
            const FixedList<label, Size>& indices,
            const bool doReduce = true
        );

        //- Construct from Istream
        inline boundBox(Istream&);


    // Member Functions

        // Access

            //- Minimum describing the bounding box
            inline const point& min() const;

            //- Maximum describing the bounding box
            inline const point& max() const;

            //- Minimum describing the bounding box, non-const access
            inline point& min();

            //- Maximum describing the bounding box, non-const access
            inline point& max();

            //- The midpoint of the bounding box
            inline point midpoint() const;

            //- The bounding box span (from minimum to maximum)
            inline vector span() const;

            //- The magnitude of the bounding box span
            inline scalar mag() const;

            //- The volume of the bound box
            inline scalar volume() const;

            //- Smallest length/height/width dimension
            inline scalar minDim() const;

            //- Largest length/height/width dimension
            inline scalar maxDim() const;

            //- Average length/height/width dimension
            inline scalar avgDim() const;

            //- Return corner points in an order corresponding to a 'hex' cell
            tmp<pointField> points() const;

            //- Return faces with correct point order
            static faceList faces();


        // Manipulate

            //- Inflate box by factor*mag(span) in all dimensions
            void inflate(const scalar s);


        // Query

            //- Overlaps/touches boundingBox?
            inline bool overlaps(const boundBox&) const;

            //- Overlaps boundingSphere (centre + sqr(radius))?
            inline bool overlaps(const point&, const scalar radiusSqr) const;

            //- Contains point? (inside or on edge)
            inline bool contains(const point&) const;

            //- Fully contains other boundingBox?
            inline bool contains(const boundBox&) const;

            //- Contains point? (inside only)
            inline bool containsInside(const point&) const;

            //- Contains all of the points? (inside or on edge)
            bool contains(const UList<point>&) const;

            //- Contains all of the points? (inside or on edge)
            bool contains
            (
                const UList<point>&,
                const labelUList& indices
            ) const;

            //- Contains all of the points? (inside or on edge)
            template<unsigned Size>
            bool contains
            (
                const UList<point>&,
                const FixedList<label, Size>& indices
            ) const;


            //- Contains any of the points? (inside or on edge)
            bool containsAny(const UList<point>&) const;

            //- Contains any of the points? (inside or on edge)
            bool containsAny
            (
                const UList<point>&,
                const labelUList& indices
            ) const;

            //- Contains any of the points? (inside or on edge)
            template<unsigned Size>
            bool containsAny
            (
                const UList<point>&,
                const FixedList<label, Size>& indices
            ) const;

            //- Return the nearest point on the boundBox to the supplied point.
            //  If point is inside the boundBox then the point is returned
            //  unchanged.
            point nearest(const point&) const;


    // Friend Operators

        inline friend bool operator==(const boundBox&, const boundBox&);
        inline friend bool operator!=(const boundBox&, const boundBox&);


    // IOstream operator

        friend Istream& operator>>(Istream&, boundBox&);
        friend Ostream& operator<<(Ostream&, const boundBox&);
};


//- Data associated with boundBox type are contiguous
template<>
inline bool contiguous<boundBox>() {return contiguous<point>();}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "boundBoxI.H"

#ifdef NoRepository
    #include "boundBoxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
