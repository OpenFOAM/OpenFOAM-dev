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
    Foam::FixedPolynomial

Description
    Fixed-size polynomial. Evaluation is done using Horner's method to minimise
    the number of operations.

SourceFiles
    FixedPolynomial.C

\*---------------------------------------------------------------------------*/

#ifndef FixedPolynomial_H
#define FixedPolynomial_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class FixedPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class T, uint8_t PN>
class FixedPolynomial
:
    public VectorSpace<FixedPolynomial<T, PN>, T, PN>
{
public:

    // Constructors

        //- Default construct
        FixedPolynomial();

        //- Construct zero
        FixedPolynomial(const zero);

        //- Construct from C-array of coefficients
        FixedPolynomial(const T (&coeffs)[PN]);

        //- Construct from a list of coefficients
        FixedPolynomial(const UList<T>& coeffs);

        //- Construct from a stream
        FixedPolynomial(Istream& is);

        //- Construct from a dictionary
        FixedPolynomial
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

        //- Return the polynomial of the (definite) integral from a coordinate
        FixedPolynomial<T, PN + 1> integralPolynomial(const scalar x) const;

        //- Return info proxy
        inline InfoProxy<FixedPolynomial<T, PN>> info() const
        {
            return *this;
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
template<class T, uint8_t PN>
void convert
(
    FixedPolynomial<T, PN>& poly,
    const Function1s::unitSets& units,
    const unitSet::makeStandardOp&
);

//- Convert to user units before writing
template<class T, uint8_t PN>
void convert
(
    FixedPolynomial<T, PN>& poly,
    const Function1s::unitSets& units,
    const unitSet::makeUserOp&
);

//- Write info to stream
template<class T, uint8_t PN>
Ostream& operator<<
(
    Ostream&,
    const InfoProxy<FixedPolynomial<T, PN>>&
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FixedPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
