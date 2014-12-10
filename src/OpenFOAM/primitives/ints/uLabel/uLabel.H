/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Typedef
    Foam::uLabel

Description
    A uLabel is an unsigned label.

SeeAlso
    label.H

\*---------------------------------------------------------------------------*/

#ifndef uLabel_H
#define uLabel_H

#include <climits>
#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#if FOAM_LABEL64
#    define FOAM_ULABEL_MAX 18000000000000000000u
#else
#    define FOAM_ULABEL_MAX 4000000000u
#endif


#if UINT_MAX > FOAM_ULABEL_MAX

// Define uLabel as an unsigned int

# undef  FOAM_ULABEL_MAX
# define FOAM_ULABEL_MAX UINT_MAX

# include "uint.H"

namespace Foam
{
    typedef unsigned int uLabel;

    static const uLabel uLabelMin = 0;
    static const uLabel uLabelMax = UINT_MAX;

    inline uLabel readULabel(Istream& is)
    {
        return readUint(is);
    }

} // End namespace Foam


#elif ULONG_MAX > FOAM_ULABEL_MAX

// Define uLabel as an unsigned long

# undef  FOAM_ULABEL_MAX
# define FOAM_ULABEL_MAX ULONG_MAX

# include "uint.H"
# include "ulong.H"

namespace Foam
{
    typedef unsigned long uLabel;

    static const uLabel uLabelMin = 0;
    static const uLabel uLabelMax = ULONG_MAX;

    inline uLabel readULabel(Istream& is)
    {
        return readUlong(is);
    }

} // End namespace Foam


#elif ULLONG_MAX > FOAM_ULABEL_MAX

// Define uLabel as an unsigned long long

# undef  FOAM_ULABEL_MAX

# error "Not implemented yet"

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pTraits.H"
#include "direction.H"

namespace Foam
{

//- template specialization for pTraits<uLabel>
template<>
class pTraits<uLabel>
{
    uLabel p_;

public:

    //- Component type
    typedef uLabel cmptType;

    // Member constants

        enum
        {
            dim = 3,         //!< Dimensionality of space
            rank = 0,        //!< Rank of uLabel is 0
            nComponents = 1  //!< Number of components in uLabel is 1
        };


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const uLabel zero;
        static const uLabel one;
        static const uLabel max;
        static const uLabel min;
        static const uLabel rootMax;
        static const uLabel rootMin;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const uLabel&);

        //- Construct from Istream
        pTraits(Istream&);


    // Member Functions

        //- Access to the uLabel value
        operator uLabel() const
        {
            return p_;
        }

        //- Access to the uLabel value
        operator uLabel&()
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline uLabel& setComponent(uLabel& l, const direction)
{
    return l;
}

inline uLabel component(const uLabel l, const direction)
{
    return l;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
