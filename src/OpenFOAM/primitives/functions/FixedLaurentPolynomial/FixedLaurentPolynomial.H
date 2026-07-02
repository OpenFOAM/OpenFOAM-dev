/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::FixedLaurentPolynomial

Description
    Fixed-size Laurent polynomial. A polynomial, but also including negative
    coefficients. Evaluation is done using Horner's method to minimise the
    number of operations.

SourceFiles
    FixedLaurentPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef FixedLaurentPolynomial_H
#define FixedLaurentPolynomial_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class FixedLaurentPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class T, int8_t P0, int8_t PN>
class FixedLaurentPolynomial
:
    public VectorSpace<FixedLaurentPolynomial<T, P0, PN>, T, PN - P0>
{
public:

    // Constructors

        //- Default construct
        FixedLaurentPolynomial();

        //- Construct zero
        FixedLaurentPolynomial(const zero);

        //- Construct from C-array of coefficients
        FixedLaurentPolynomial(const T (&coeffs)[PN - P0]);

        //- Construct from a list of coefficients
        FixedLaurentPolynomial(const UList<T>& coeffs);

        //- Construct from a stream
        FixedLaurentPolynomial(Istream& is);

        //- Construct from a dictionary
        FixedLaurentPolynomial
        (
            const char* symbols,
            const Function1s::unitSets& units,
            const dictionary& dict
        );


    // Member Functions

        //- Return the value at a coordinate
        T value(const scalar x) const;

        //- Return the value of the derivative at a coordinate
        T derivative(const scalar x) const;

        //- Return the value of the (indefinite) integral at a coordinate
        T integral(const scalar x) const;

        //- Return the value of the (definite) integral between coordinates
        T integral(const scalar x0, const scalar x1) const;

        //- Return info proxy
        inline InfoProxy<FixedLaurentPolynomial<T, P0, PN>> info() const
        {
            return *this;
        }

        //- Return the polynomial multiplied by the coordinate. Just cast to a
        //  type with the powers increased by one.
        const FixedLaurentPolynomial<T, P0 + 1, PN + 1>& timesX() const
        {
            return
                reinterpret_cast
                <
                    const FixedLaurentPolynomial<T, P0 + 1, PN + 1>&
                >(*this);
        }

        //- Return the polynomial divided by the coordinate. Just cast to a
        //  type with the powers decreased by one.
        const FixedLaurentPolynomial<T, P0 - 1, PN - 1>& byX() const
        {
            return
                reinterpret_cast
                <
                    const FixedLaurentPolynomial<T, P0 - 1, PN - 1>&
                >(*this);
        }

        //- Write to stream
        void write
        (
            const char* symbols,
            const Function1s::unitSets& units,
            Ostream& os
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Convert to standard units after reading
template<class T, int8_t P0, int8_t PN>
void convert
(
    FixedLaurentPolynomial<T, P0, PN>& poly,
    const Function1s::unitSets& units,
    const unitSet::makeStandardOp&
);

//- Convert to user units before writing
template<class T, int8_t P0, int8_t PN>
void convert
(
    FixedLaurentPolynomial<T, P0, PN>& poly,
    const Function1s::unitSets& units,
    const unitSet::makeUserOp&
);

//- Write info to stream
template<class T, int8_t P0, int8_t PN>
Ostream& operator<<
(
    Ostream&,
    const InfoProxy<FixedLaurentPolynomial<T, P0, PN>>&
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FixedLaurentPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
