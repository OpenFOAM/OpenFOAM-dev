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
    Foam::labelRange

Description
    A label range specifier.

SourceFiles
    labelRange.C

\*---------------------------------------------------------------------------*/

#ifndef labelRange_H
#define labelRange_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class labelRange;
Istream& operator>>(Istream&, labelRange&);
Ostream& operator<<(Ostream&, const labelRange&);

/*---------------------------------------------------------------------------*\
                         Class labelRange Declaration
\*---------------------------------------------------------------------------*/

class labelRange
{
    // Private Data

        label start_;
        label size_;

public:

        static int debug;


    // Public classes

        //- Less function class for sorting labelRange
        class less
        {
        public:

            bool operator()(const labelRange& a, const labelRange& b)
            {
                return a.first() < b.first();
            }
        };


    // Constructors

        //- Construct an empty range
        inline labelRange();

        //- Construct a range
        //  A negative size is automatically changed to zero.
        inline labelRange(const label start, const label size);

        //- Construct from Istream.
        labelRange(Istream&);


    // Member Functions

        //- Reset to zero size
        inline void clear();

        //- Is the range empty?
        inline bool empty() const;

        //- Return the effective size of the range
        inline label size() const;

        //- The (inclusive) lower value of the range
        inline label first() const;

        //- The (inclusive) upper value of the range
        inline label last() const;

        //- Return true if the value is within the range
        inline bool contains(const label) const;

        //- Return true if the ranges intersect
        //  Optional test for ranges that also just touch each other
        bool intersects(const labelRange&, const bool touches = false) const;

        //- Return a joined range, squashing any gaps in between
        //  A prior intersects() check can be used to avoid squashing gaps.
        labelRange join(const labelRange&) const;


    // Member Operators

        //- Return element in range, no bounds checking
        inline label operator[](const label) const;

        //- Comparison function for sorting, compares the start
        inline bool operator<(const labelRange& rhs) const;

        //- Join ranges, squashing any gaps in between
        //  A prior intersects() check can be used to avoid squashing gaps.
        labelRange& operator+=(const labelRange&);

        inline bool operator==(const labelRange&) const;
        inline bool operator!=(const labelRange&) const;


    // STL iterator

        //- An STL const_iterator
        class const_iterator
        {
            // Private Data

                //- Reference to the range for which this is an iterator
                const labelRange& range_;

                //- Current index
                label index_;

        public:

            // Constructors

                //- Construct null - equivalent to an 'end' position
                inline const_iterator();

                //- Construct from range, moving to its 'begin' position
                inline explicit const_iterator(const labelRange&);


            // Member Operators

                inline bool operator==(const const_iterator&) const;

                inline bool operator!=(const const_iterator&) const;

                inline label operator*();
                inline label operator()();

                inline const_iterator& operator++();
                inline const_iterator operator++(int);
        };


        //- const_iterator set to the beginning of the range
        inline const_iterator cbegin() const;

        //- const_iterator set to beyond the end of the range
        inline const const_iterator& cend() const;

        //- const_iterator set to the beginning of the range
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the range
        inline const const_iterator& end() const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, labelRange&);
        friend Ostream& operator<<(Ostream&, const labelRange&);


private:

        //- const labelRange held by endIter_
        static const labelRange endLabelRange_;

        //- const_iterator returned by end(), cend()
        static const const_iterator endIter_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelRangeI.H"

#endif

// ************************************************************************* //
