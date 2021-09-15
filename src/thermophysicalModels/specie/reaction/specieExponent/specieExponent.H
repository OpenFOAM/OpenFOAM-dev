/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

SourceFiles
    specieExponentI.H
    specieExponent.C
    specieExponentIO.C

\*---------------------------------------------------------------------------*/

#ifndef specieExponent_H
#define specieExponent_H

#include "label.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class specieExponent;
scalar pow(scalar x, const specieExponent& e);
specieExponent operator+(const specieExponent& a, const specieExponent& b);
specieExponent operator-(const specieExponent& a, const specieExponent& b);
Ostream& operator<<(Ostream&, const specieExponent&);


/*---------------------------------------------------------------------------*\
                          Class specieExponent Declaration
\*---------------------------------------------------------------------------*/

class specieExponent
{
    // Private Data

        //- Marker value to indicate that there is no possible integer
        //  representation of this exponent
        static const label noIntegerExponent_ = labelMax;

        //- Integer exponent
        label integerExponent_;

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
        inline specieExponent(const label integerExponent);

        //- Construct from scalar
        inline specieExponent(const scalar scalarExponent);


    // Member Operators

        //- Cast to scalar
        inline operator scalar() const;

        //- Assign to integer
        inline specieExponent& operator=(const label integerExponent);

        //- Assign to scalar
        inline specieExponent& operator=(const scalar scalarExponent);

        //- Negate a specie exponent
        inline specieExponent operator-() const;


    // Friend Functions

        //- Compute the power of a number to a specie exponent
        inline friend scalar pow(const scalar x, const specieExponent& e);


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
