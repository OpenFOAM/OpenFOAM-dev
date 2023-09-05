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

\*---------------------------------------------------------------------------*/

#include "scalarAndError.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Arithmetic operations between two scalar-and-errors
#define SCALAR_AND_ERROR_BINARY_OPERATION(Op)                                  \
                                                                               \
    inline scalarAndError operator Op                                          \
    (                                                                          \
        const scalarAndError& a,                                               \
        const scalarAndError& b                                                \
    )                                                                          \
    {                                                                          \
        scalarAndError c(a);                                                   \
        c Op##= b;                                                             \
        return c;                                                              \
    }
SCALAR_AND_ERROR_BINARY_OPERATION(+);
SCALAR_AND_ERROR_BINARY_OPERATION(-);
SCALAR_AND_ERROR_BINARY_OPERATION(*);
SCALAR_AND_ERROR_BINARY_OPERATION(/);


//- Equality compare two scalar-and-error-s
inline bool operator==(const scalarAndError& a, const scalarAndError& b)
{
    return a.value == b.value && a.error == b.error;
}


//- Inequality compare two scalar-and-error-s
inline bool operator!=(const scalarAndError& a, const scalarAndError& b)
{
    return !(a == b);
}


//- Comparisons between two scalar-and-errors
#define SCALAR_AND_ERROR_BINARY_COMPARE(Op)                                    \
                                                                               \
    inline bool operator Op                                                    \
    (                                                                          \
        const scalarAndError& a,                                               \
        const scalarAndError& b                                                \
    )                                                                          \
    {                                                                          \
        return a.value Op b.value;                                             \
    }
SCALAR_AND_ERROR_BINARY_COMPARE(<);
SCALAR_AND_ERROR_BINARY_COMPARE(<=);
SCALAR_AND_ERROR_BINARY_COMPARE(>=);
SCALAR_AND_ERROR_BINARY_COMPARE(>);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the magnitude square of a scalar-and-error
inline scalarAndError magSqr(const scalarAndError& sae)
{
    return sae*sae;
}


//- Return the square root of a scalar-and-error
inline scalarAndError sqrt(const scalarAndError& sae)
{
    const scalar sqrtV = ::sqrt(sae.value);
    return scalarAndError(sqrtV, 0.5*sqrtV*sae.error);
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

//- Write a scalar-and-error
inline Ostream& operator<<(Ostream& os, const scalarAndError& sae)
{
    return os << sae.value << "(+-" << sae.error << ')';
}


//- Read a scalar-and-error
inline Istream& operator>>(Istream& is, scalarAndError& sae)
{
    auto readPunctuationToken = [&](const token::punctuationToken p)
    {
        token t(is);

        if (t != p)
        {
            is.setBad();

            FatalIOErrorInFunction(is)
                << "Expected a '" << p
                << "' while reading scalarAndError"
                << ", found " << t.info()
                << exit(FatalIOError);
        }
    };

    is >> sae.value;
    readPunctuationToken(token::BEGIN_LIST);
    readPunctuationToken(token::ADD);
    readPunctuationToken(token::SUBTRACT);
    is >> sae.error;
    readPunctuationToken(token::END_LIST);

    return is;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
