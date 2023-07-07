/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::Tuple3

Description
    A 3-tuple for storing three objects of different types.

\*---------------------------------------------------------------------------*/

#ifndef Tuple3_H
#define Tuple3_H

#include "Istream.H"
#include "Hash.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type1, class Type2, class Type3>
class Tuple3;

template<class Type1, class Type2, class Type3>
inline Istream& operator>>(Istream&, Tuple3<Type1, Type2, Type3>&);

template<class Type1, class Type2, class Type3>
inline Ostream& operator<<(Ostream&, const Tuple3<Type1, Type2, Type3>&);


/*---------------------------------------------------------------------------*\
                           Class Tuple3 Declaration
\*---------------------------------------------------------------------------*/

template<class Type1, class Type2, class Type3>
class Tuple3
{
    // Private Data

        Type1 f_;
        Type2 s_;
        Type3 t_;


public:

    //- Hashing function class
    template
    <
        class HashT1=Hash<Type1>,
        class HashT2=Hash<Type2>,
        class HashT3=Hash<Type3>
    >
    class Hash
    {
    public:
        Hash()
        {}

        inline unsigned operator()
        (
            const Tuple3<Type1, Type2, Type3>&,
            unsigned seed = 0
        ) const;
    };


    // Static Data Members

        static const char* const typeName;


    // Constructors

        //- Null constructor for lists
        inline Tuple3()
        {}

        //- Construct from components
        inline Tuple3(const Type1& f, const Type2& s, const Type3& t)
        :
            f_(f),
            s_(s),
            t_(t)
        {}

        //- Construct from Istream
        inline Tuple3(Istream& is)
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

        //- Return third
        inline const Type3& third() const
        {
            return t_;
        }

        //- Return third
        inline Type3& third()
        {
            return t_;
        }

    // IOstream Operators

        //- Read Tuple3 from Istream, discarding contents of existing Tuple3.
        friend Istream& operator>> <Type1, Type2, Type3>
        (
            Istream& is,
            Tuple3<Type1, Type2, Type3>& t3
        );

        // Write Tuple3 to Ostream.
        friend Ostream& operator<< <Type1, Type2, Type3>
        (
            Ostream& os,
            const Tuple3<Type1, Type2, Type3>& t3
        );
};


template<class Type1, class Type2, class Type3>
template<class HashT1, class HashT2, class HashT3>
inline unsigned
Tuple3<Type1, Type2, Type3>::Hash<HashT1, HashT2, HashT3>::operator()
(
    const Tuple3<Type1, Type2, Type3>& t,
    unsigned seed
) const
{
    // Hash incrementally
    unsigned val = seed;
    val = HashT1()(t.first(), val);
    val = HashT2()(t.second(), val);
    val = HashT3()(t.third(), val);
    return val;
}


//- Return reverse of a tuple3
template<class Type1, class Type2, class Type3>
inline Tuple3<Type3, Type2, Type1> reverse
(
    const Tuple3<Type1, Type2, Type3>& t
)
{
    return Tuple3<Type3, Type2, Type1>(t.third(), t.second(), t.first());
}


template<class Type1, class Type2, class Type3>
inline bool operator==
(
    const Tuple3<Type1, Type2, Type3>& a,
    const Tuple3<Type1, Type2, Type3>& b
)
{
    return
    (
        a.first() == b.first()
     && a.second() == b.second()
     && a.third() == b.third()
    );
}


template<class Type1, class Type2, class Type3>
inline bool operator!=
(
    const Tuple3<Type1, Type2, Type3>& a,
    const Tuple3<Type1, Type2, Type3>& b
)
{
    return !(a == b);
}


template<class Type1, class Type2, class Type3>
inline Istream& operator>>(Istream& is, Tuple3<Type1, Type2, Type3>& t3)
{
    is.readBegin("Tuple3");
    is >> t3.f_ >> t3.s_ >> t3.t_;
    is.readEnd("Tuple3");

    // Check state of Istream
    is.check("operator>>(Istream&, Tuple3<Type1, Type2, Type3>&)");

    return is;
}


template<class Type1, class Type2, class Type3>
inline Ostream& operator<<(Ostream& os, const Tuple3<Type1, Type2, Type3>& t2)
{
    os  << token::BEGIN_LIST
        << t2.f_ << token::SPACE << t2.s_ << token::SPACE << t2.t_
        << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
