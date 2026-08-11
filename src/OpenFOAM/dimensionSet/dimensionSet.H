/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
#include "fieldTypes.H"
#include "className.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class dimensionSet;

template<class Type> class dimensioned;
typedef dimensioned<scalar> dimensionedScalar;

// Friend Functions

dimensionSet pow(const dimensionSet& ds, const dimensionedScalar& dS);
dimensionSet pow(const dimensionedScalar& dS, const dimensionSet& ds);

// Friend operators

dimensionSet operator*(const dimensionSet&, const dimensionSet&);
dimensionSet operator/(const dimensionSet&, const dimensionSet&);

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
        static const NamedEnum<dimensionType, 7>& dimensionTypeNames_;


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

        //- Construct a named dimensionSet given individual dimension exponents
        //  for all seven dimensions.
        //  The dimensionSet is added to the dimensions::table
        dimensionSet
        (
            const word& name,
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

        //- Named copy constructor
        dimensionSet(const word& name, const dimensionSet& ds);

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

        inline scalar operator[](const dimensionType) const;
        inline scalar& operator[](const dimensionType);

        inline scalar operator[](const label) const;
        inline scalar& operator[](const label);

        bool operator==(const dimensionSet&) const;
        inline bool operator!=(const dimensionSet&) const;

        inline bool operator=(const dimensionSet&) const;

        inline bool operator+=(const dimensionSet&) const;
        inline bool operator-=(const dimensionSet&) const;
        inline bool operator*=(const dimensionSet&);
        inline bool operator/=(const dimensionSet&);

        template<class ... Args>
        using typeForArgs = typename typeOfNcmpts<sizeof ... (Args)>::type;

        //- Return a dimensionedType with these dimensions
        //  and rank corresponding to the number of arguments
        template<class ... Args>
        dimensioned<typeForArgs<Args ...>> operator()
        (
            const Args& ... args
        ) const
        {
            return dimensioned<typeForArgs<Args ...>>
            (
                *this,
                typeForArgs<Args ...>(args ...)
            );
        }

        //- Return a named dimensionedType with these dimensions
        //  and rank corresponding to the number of arguments
        template<class ... Args>
        dimensioned<typeForArgs<Args ...>> operator()
        (
            const word& name,
            const Args& ... args
        ) const
        {
            return dimensioned<typeForArgs<Args ...>>
            (
                name,
                *this,
                typeForArgs<Args ...>(args ...)
            );
        }


    // Friend Functions

        friend dimensionSet pow(const dimensionSet&, const dimensionedScalar&);
        friend dimensionSet pow(const dimensionedScalar&, const dimensionSet&);


    // Friend Operators

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


    // IOstream Operators

        friend Istream& operator>>(Istream&, dimensionSet&);
        friend Ostream& operator<<(Ostream&, const dimensionSet&);
        friend Ostream& operator<<(Ostream&, const InfoProxy<dimensionSet>&);
};


void writeEntry(Ostream& os, const dimensionSet& value);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dimensions.H"
#include "dimensionSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
