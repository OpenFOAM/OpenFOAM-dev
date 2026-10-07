/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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
    Foam::specieExponent

Description
    Type for exponents of species in reaction. Can take either integer or
    scalar value, and powers are optimised in the integer case.

SourceFiles
    specieExponentI.H

\*---------------------------------------------------------------------------*/

#ifndef specieExponent_H
#define specieExponent_H

#include "int32.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class specieExponent;

inline scalar pow(scalar c, const specieExponent& e);

inline scalar dpowdc(const scalar c, const specieExponent& e);

inline specieExponent operator+
(
    const specieExponent& a,
    const specieExponent& b
);

inline specieExponent operator-
(
    const specieExponent& a,
    const specieExponent& b
);

inline specieExponent operator-(const specieExponent& a, const int32_t b);

inline Ostream& operator<<(Ostream&, const specieExponent&);


/*---------------------------------------------------------------------------*\
                          Class specieExponent Declaration
\*---------------------------------------------------------------------------*/

class specieExponent
{
    // Private Data

        //- Marker value to indicate that there is no possible integer
        //  representation of this exponent
        static const int32_t noIntegerExponent_ = INT32_MAX;

        //- Integer exponent
        int32_t integerExponent_;

        //- Scalar exponent
        scalar scalarExponent_;


    // Private Member Functions

        //- Return whether or not this exponent has an integer representation
        inline bool hasIntegerExponent() const;


public:

    // Constructors

        //- Construct null
        inline specieExponent();

        //- Construct from integer
        inline specieExponent(const int32_t integerExponent);

        //- Construct from scalar
        inline specieExponent(const scalar scalarExponent);


    // Member Operators

        //- Cast to scalar
        inline operator scalar() const;

        //- Assign to integer
        inline specieExponent& operator=(const int32_t integerExponent);

        //- Assign to scalar
        inline specieExponent& operator=(const scalar scalarExponent);

        //- Negate a specie exponent
        inline specieExponent operator-() const;


    // Friend Functions

        //- Compute the power of a number to a specie exponent
        inline friend scalar pow(const scalar c, const specieExponent& e);

        //- Compute the derivative of the power of a number to a specie exponent
        inline friend scalar dpowdc(const scalar c, const specieExponent& e);


    // Friend Operators

        //- Sum two specie exponents
        inline friend specieExponent operator+
        (
            const specieExponent& a,
            const specieExponent& b
        );

        //- Subtract two specie exponents
        inline friend specieExponent operator-
        (
            const specieExponent& a,
            const specieExponent& b
        );

        //- Subtract integer from specie exponent
        inline friend specieExponent operator-
        (
            const specieExponent& a,
            const int32_t b
        );


    // IOstream Operators

        //- Write to output stream
        inline friend Ostream& operator<<(Ostream& os, const specieExponent& e);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specieExponentI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
