/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::Tuple2

Description
    A 2-tuple for storing two objects of different types.

See also
    Foam::Pair for storing two objects of identical types.

\*---------------------------------------------------------------------------*/

#ifndef Tuple2_H
#define Tuple2_H

#include "Istream.H"
#include "Hash.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type1, class Type2>
class Tuple2;

template<class Type1, class Type2>
inline Istream& operator>>(Istream&, Tuple2<Type1, Type2>&);

template<class Type1, class Type2>
inline Ostream& operator<<(Ostream&, const Tuple2<Type1, Type2>&);


/*---------------------------------------------------------------------------*\
                           class Tuple2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type1, class Type2>
class Tuple2
{
    // Private Data

        Type1 f_;
        Type2 s_;


public:

    //- Hashing function class
    template<class HashT1=Hash<Type1>, class HashT2=Hash<Type2>>
    class Hash
    {
    public:
        Hash()
        {}

        inline unsigned operator()
        (
            const Tuple2<Type1, Type2>&,
            unsigned seed = 0
        ) const;
    };


    // Static Data Members

        static const char* const typeName;


    // Constructors

        //- Null constructor for lists
        inline Tuple2()
        {}

        //- Construct from components
        inline Tuple2(const Type1& f, const Type2& s)
        :
            f_(f),
            s_(s)
        {}

        //- Construct from Istream
        inline Tuple2(Istream& is)
        {
            is >> *this;
        }


    // Member Functions

        //- Return first
        inline const Type1& first() const
        {
            return f_;
        }

        //- Return first
        inline Type1& first()
        {
            return f_;
        }

        //- Return second
        inline const Type2& second() const
        {
            return s_;
        }

        //- Return second
        inline Type2& second()
        {
            return s_;
        }


    // IOstream Operators

        //- Read Tuple2 from Istream, discarding contents of existing Tuple2.
        friend Istream& operator>> <Type1, Type2>
        (
            Istream& is,
            Tuple2<Type1, Type2>& t2
        );

        // Write Tuple2 to Ostream.
        friend Ostream& operator<< <Type1, Type2>
        (
            Ostream& os,
            const Tuple2<Type1, Type2>& t2
        );
};


template<class Type1, class Type2>
template<class HashT1, class HashT2>
inline unsigned Tuple2<Type1, Type2>::Hash<HashT1, HashT2>::operator()
(
    const Tuple2<Type1, Type2>& t,
    unsigned seed
) const
{
    // Hash incrementally
    unsigned val = seed;
    val = HashT1()(t.first(), val);
    val = HashT2()(t.second(), val);
    return val;
}


//- Return reverse of a tuple2
template<class Type1, class Type2>
inline Tuple2<Type2, Type1> reverse(const Tuple2<Type1, Type2>& t)
{
    return Tuple2<Type2, Type1>(t.second(), t.first());
}


template<class Type1, class Type2>
inline bool operator==
(
    const Tuple2<Type1, Type2>& a,
    const Tuple2<Type1, Type2>& b
)
{
    return (a.first() == b.first() && a.second() == b.second());
}


template<class Type1, class Type2>
inline bool operator!=
(
    const Tuple2<Type1, Type2>& a,
    const Tuple2<Type1, Type2>& b
)
{
    return !(a == b);
}


template<class Type1, class Type2>
inline Istream& operator>>(Istream& is, Tuple2<Type1, Type2>& t2)
{
    is.readBegin("Tuple2");
    is >> t2.f_ >> t2.s_;
    is.readEnd("Tuple2");

    // Check state of Istream
    is.check("operator>>(Istream&, Tuple2<Type1, Type2>&)");

    return is;
}


template<class Type1, class Type2>
inline Ostream& operator<<(Ostream& os, const Tuple2<Type1, Type2>& t2)
{
    os  << token::BEGIN_LIST
        << t2.f_ << token::SPACE << t2.s_
        << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
