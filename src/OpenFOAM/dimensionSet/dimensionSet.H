/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::dimensionSet

Description
    Dimension set for the base types.

    This type may be used to implement rigorous dimension checking
    for algebraic manipulation.

SourceFiles
    dimensionSet.C
    dimensionSetIO.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionSet_H
#define dimensionSet_H

#include "bool.H"
#include "className.H"
#include "NamedEnum.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class dimensionSet;
class unitConversion;
template<class Type> class dimensioned;
typedef dimensioned<scalar> dimensionedScalar;

// Friend Functions

dimensionSet max(const dimensionSet&, const dimensionSet&);
dimensionSet min(const dimensionSet&, const dimensionSet&);
dimensionSet cmptMultiply(const dimensionSet&, const dimensionSet&);
dimensionSet cmptDivide(const dimensionSet&, const dimensionSet&);
dimensionSet cmptMag(const dimensionSet&);

dimensionSet pow(const dimensionSet&, const scalar);
dimensionSet pow(const dimensionSet&, const dimensionedScalar&);
dimensionSet pow(const dimensionedScalar&, const dimensionSet&);

dimensionSet sqr(const dimensionSet&);
dimensionSet pow3(const dimensionSet&);
dimensionSet pow4(const dimensionSet&);
dimensionSet pow5(const dimensionSet&);
dimensionSet pow6(const dimensionSet&);
dimensionSet pow025(const dimensionSet&);

dimensionSet sqrt(const dimensionSet&);
dimensionSet cbrt(const dimensionSet&);
dimensionSet magSqr(const dimensionSet&);
dimensionSet mag(const dimensionSet&);
dimensionSet sign(const dimensionSet&);
dimensionSet pos(const dimensionSet&);
dimensionSet pos0(const dimensionSet&);
dimensionSet neg(const dimensionSet&);
dimensionSet neg0(const dimensionSet&);
dimensionSet posPart(const dimensionSet&);
dimensionSet negPart(const dimensionSet&);
dimensionSet inv(const dimensionSet&);

// Function to check the argument is dimensionless
//  for transcendental functions
dimensionSet trans(const dimensionSet&);

dimensionSet atan2(const dimensionSet&, const dimensionSet&);

// Return the argument; transformations do not change the dimensions
dimensionSet transform(const dimensionSet&);

dimensionSet normalised(const dimensionSet&);
dimensionSet perpendicular(const dimensionSet&);

// Friend operators

dimensionSet operator-(const dimensionSet&);
dimensionSet operator+(const dimensionSet&, const dimensionSet&);
dimensionSet operator-(const dimensionSet&, const dimensionSet&);
dimensionSet operator*(const dimensionSet&, const dimensionSet&);
dimensionSet operator/(const dimensionSet&, const dimensionSet&);
dimensionSet operator&(const dimensionSet&, const dimensionSet&);
dimensionSet operator^(const dimensionSet&, const dimensionSet&);
dimensionSet operator&&(const dimensionSet&, const dimensionSet&);

// IOstream Operators

Istream& operator>>(Istream&, dimensionSet&);
Ostream& operator<<(Ostream&, const dimensionSet&);
Ostream& operator<<(Ostream&, const InfoProxy<dimensionSet>&);


/*---------------------------------------------------------------------------*\
                        Class dimensionSet Declaration
\*---------------------------------------------------------------------------*/

class dimensionSet
{

public:

    // Member constants

        //- Define an enumeration for the number of dimensions
        enum
        {
            nDimensions = 7    // Number of dimensions in SI is 7
        };

        //- Define an enumeration for the names of the dimension exponents
        enum dimensionType
        {
            MASS,               // kilogram   kg
            LENGTH,             // metre      m
            TIME,               // second     s
            TEMPERATURE,        // Kelvin     K
            MOLES,              // kilomole   kmol
            CURRENT,            // Ampere     A
            LUMINOUS_INTENSITY  // Candela    Cd
        };

        //- Names of the dimensions
        static const NamedEnum<dimensionType, 7> dimensionTypeNames_;


    // Static Data Members

        //- A small exponent with which to perform inexact comparisons
        static const scalar smallExponent;


private:

    // Private Data

        //- Array of dimension exponents
        scalar exponents_[nDimensions];


    // Private Member Functions

        //- Reset exponents to nearest integer if close to it. Used to
        //  handle reading with insufficient precision.
        void round(const scalar tol);


public:

    // Declare name of the class and its debug switch
    ClassName("dimensionSet");


    // Constructors

        //- Construct given individual dimension exponents for all
        //  seven dimensions
        dimensionSet
        (
            const scalar mass,
            const scalar length,
            const scalar time,
            const scalar temperature,
            const scalar moles,
            const scalar current,
            const scalar luminousIntensity
        );

        //- Construct given individual dimension exponents for first
        //  five dimensions
        dimensionSet
        (
            const scalar mass,
            const scalar length,
            const scalar time,
            const scalar temperature,
            const scalar moles
        );

        //- Copy constructor
        dimensionSet(const dimensionSet& ds);

        //- Construct and return a clone
        autoPtr<dimensionSet> clone() const
        {
            return autoPtr<dimensionSet>(new dimensionSet(*this));
        }

        //- Construct from Istream
        dimensionSet(Istream&);


    // Member Functions

        //- Return true if it is dimensionless
        bool dimensionless() const;

        void reset(const dimensionSet&);


    // I/O

        //- Read without the square brackets
        Istream& readNoBeginOrEnd(Istream& is);

        //- Read
        Istream& read(Istream& is);

        //- Write without the square brackets
        Ostream& writeNoBeginOrEnd(Ostream& os) const;

        //- Write
        Ostream& write(Ostream& os) const;

        //- Write info without the square brackets
        Ostream& writeInfoNoBeginOrEnd(Ostream& os) const;

        //- Return info proxy
        inline InfoProxy<dimensionSet> info() const
        {
            return *this;
        }


    // Member Operators

        scalar operator[](const dimensionType) const;
        scalar& operator[](const dimensionType);

        scalar operator[](const label) const;
        scalar& operator[](const label);

        bool operator==(const dimensionSet&) const;
        bool operator!=(const dimensionSet&) const;

        bool operator=(const dimensionSet&) const;

        bool operator+=(const dimensionSet&) const;
        bool operator-=(const dimensionSet&) const;
        bool operator*=(const dimensionSet&);
        bool operator/=(const dimensionSet&);


    // Friend Functions

        friend dimensionSet max(const dimensionSet&, const dimensionSet&);
        friend dimensionSet min(const dimensionSet&, const dimensionSet&);
        friend dimensionSet cmptMultiply
        (
            const dimensionSet&,
            const dimensionSet&
        );
        friend dimensionSet cmptDivide
        (
            const dimensionSet&,
            const dimensionSet&
        );
        friend dimensionSet cmptMag(const dimensionSet&);

        friend dimensionSet pow(const dimensionSet&, const scalar);
        friend dimensionSet pow(const dimensionSet&, const dimensionedScalar&);
        friend dimensionSet pow(const dimensionedScalar&, const dimensionSet&);

        friend dimensionSet sqr(const dimensionSet&);
        friend dimensionSet pow3(const dimensionSet&);
        friend dimensionSet pow4(const dimensionSet&);
        friend dimensionSet pow5(const dimensionSet&);
        friend dimensionSet pow6(const dimensionSet&);
        friend dimensionSet pow025(const dimensionSet&);

        friend dimensionSet sqrt(const dimensionSet&);
        friend dimensionSet magSqr(const dimensionSet&);
        friend dimensionSet mag(const dimensionSet&);
        friend dimensionSet sign(const dimensionSet&);
        friend dimensionSet pos0(const dimensionSet&);
        friend dimensionSet neg(const dimensionSet&);
        friend dimensionSet inv(const dimensionSet&);

        //- Function to check the argument is dimensionless
        //  for transcendental functions
        friend dimensionSet trans(const dimensionSet&);

        friend dimensionSet atan2(const dimensionSet&, const dimensionSet&);

        friend dimensionSet transform(const dimensionSet&);

        friend dimensionSet normalised(const dimensionSet&);
        friend dimensionSet perpendicular(const dimensionSet&);


    // Friend Operators

        friend dimensionSet operator-(const dimensionSet&);

        friend dimensionSet operator+
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator-
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator*
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator/
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator&
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator^
        (
            const dimensionSet&,
            const dimensionSet&
        );

        friend dimensionSet operator&&
        (
            const dimensionSet&,
            const dimensionSet&
        );


    // IOstream Operators

        friend Istream& operator>>(Istream&, dimensionSet&);
        friend Ostream& operator<<(Ostream&, const dimensionSet&);
        friend Ostream& operator<<(Ostream&, const InfoProxy<dimensionSet>&);
};


void writeEntry(Ostream& os, const dimensionSet& value);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dimensionSets.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
