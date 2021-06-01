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
    Foam::edge

Description
    An edge is a list of two point labels. The functionality it provides
    supports the discretisation on a 2-D flat mesh.

SourceFiles
    edgeI.H

\*---------------------------------------------------------------------------*/

#ifndef edge_H
#define edge_H

#include "FixedList.H"
#include "pointField.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class edge;
inline bool operator==(const edge& a, const edge& b);
inline bool operator!=(const edge& a, const edge& b);


/*---------------------------------------------------------------------------*\
                            Class edge Declaration
\*---------------------------------------------------------------------------*/

class edge
:
    public FixedList<label, 2>
{

public:

    // Static Data Members

        static const char* const typeName;


    // Constructors

        //- Null constructor for lists
        inline edge();

        //- Construct from components
        inline edge(const label a, const label b);

        //- Construct from FixedList
        inline edge(const FixedList<label, 2>&);

        //- Construct from Istream
        inline edge(Istream&);


    // Member Functions

        //- Return start vertex label
        inline label start() const;

        //- Return start vertex label
        inline label& start();

        //- Return end vertex label
        inline label end() const;

        //- Return end vertex label
        inline label& end();

        //- Return true if connected to given edge
        inline bool connected(const edge& a) const;

        //- Return common vertex
        inline label commonVertex(const edge& a) const;

        //- Given one vertex, return the other
        inline label otherVertex(const label a) const;

        //- Flip the edge in-place.
        inline void flip();

        //- Return reverse edge
        inline edge reverseEdge() const;

        //- Return centre (centroid)
        inline point centre(const pointField&) const;

        //- Return the vector (end - start)
        inline vector vec(const pointField&) const;

        //- Return scalar magnitude
        inline scalar mag(const pointField&) const;

        //- Return edge line
        inline linePointRef line(const pointField&) const;

        //- Compare edges
        //  Returns:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        static inline int compare(const edge&, const edge&);


    // Friend Operators

        friend bool operator==(const edge& a, const edge& b);
        friend bool operator!=(const edge& a, const edge& b);
};


//- Hash specialisation for hashing edges - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<edge>::operator()(const edge& e, unsigned seed) const
{
    unsigned val = seed;

    if (e[0] < e[1])
    {
        val = Hash<label>()(e[0], val);
        val = Hash<label>()(e[1], val);
    }
    else
    {
        val = Hash<label>()(e[1], val);
        val = Hash<label>()(e[0], val);
    }

    return val;
}


//- Hash specialisation for hashing edges - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<edge>::operator()(const edge& e) const
{
    return Hash<edge>()(e, 0);
}


template<>
inline bool contiguous<edge>()  {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
