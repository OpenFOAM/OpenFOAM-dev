/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::Pair

Description
    An ordered pair of two objects of type \<Type\> with first() and second()
    elements.

See also
    Foam::Tuple2 for storing two objects of dissimilar types.

\*---------------------------------------------------------------------------*/

#ifndef Pair_H
#define Pair_H

#include "Hash.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class Pair;

template<class Type>
void writeEntry(Ostream& os, const Pair<Type>&);

template<class Type>
inline Istream& operator>>(Istream&, Pair<Type>&);

template<class Type>
inline Ostream& operator<<(Ostream&, const Pair<Type>&);


/*---------------------------------------------------------------------------*\
                            Class Pair Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Pair
{
    // Private Data

        //- First object
        Type f_;

        //- Second object
        Type s_;


public:

    // Public Sub-Classes

        //- Hashing function class
        template<class HashT=Hash<Type>>
        class Hash
        {
        public:

            // Constructors

                //- Null constructor
                inline Hash();


            // Member Operators

                //- Hash a pair
                inline unsigned operator()
                (
                    const Pair<Type>& p,
                    unsigned seed = 0
                ) const;
        };


    // Static Member Functions

        //- Return a null pair
        static inline const Pair<Type>& null();

        //- Return the size
        static inline label size();

        //- Check index i is within valid range (0 ... 1)
        static inline void checkIndex(const label i);

        //- Compare two pairs. Return 0 if they are different, +1 if they are
        //  identical, and -1 if the elements are the same but reversed.
        static inline int compare(const Pair<Type>& a, const Pair<Type>& b);


    // Constructors

        //- Null constructor
        inline Pair();

        //- Construct from components
        inline Pair(const Type& f, const Type& s);

        //- Construct from Istream
        inline Pair(Istream& is);


    // Member Functions

        //- Return first
        inline const Type& first() const;

        //- Return first
        inline Type& first();

        //- Return second
        inline const Type& second() const;

        //- Return second
        inline Type& second();

        //- Return other
        inline const Type& other(const Type& a) const;


    // Member Operators

        //- Access an element by index
        inline Type& operator[](const label i);

        //- Const-access an element by index
        inline const Type& operator[](const label i) const;


    // STL type definitions

        //- Type of values the pair contains
        typedef Type value_type;

        //- Type that can be used for storing into
        //  pair::value_type objects
        typedef Type& reference;

        //- Type that can be used for storing into
        //  constant pair::value_type objects
        typedef const Type& const_reference;
};


// Global Functions

//- Reverse the elements in a pair
template<class Type>
inline Pair<Type> reverse(const Pair<Type>& p);

//- Equality comparison
template<class Type>
inline bool operator==(const Pair<Type>& a, const Pair<Type>& b);

//- Inequality comparison
template<class Type>
inline bool operator!=(const Pair<Type>& a, const Pair<Type>& b);

//- Compare lexographic order
template<class Type>
inline bool operator<(const Pair<Type>& a, const Pair<Type>& b);

//- Compare lexographic order
template<class Type>
inline bool operator<=(const Pair<Type>& a, const Pair<Type>& b);

//- Compare lexographic order
template<class Type>
inline bool operator>(const Pair<Type>& a, const Pair<Type>& b);

//- Compare lexographic order
template<class Type>
inline bool operator>=(const Pair<Type>& a, const Pair<Type>& b);

//- Read a from a stream
template<class Type>
inline Istream& operator>>(Istream& is, Pair<Type>& p);

//- Write to a stream
template<class Type>
inline Ostream& operator<<(Ostream& os, const Pair<Type>& p);

//- Write as a dictionary entry
template<class Type>
inline void writeEntry(Ostream& os, const Pair<Type>& p);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
