/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::unitConversion

Description
    Unit conversion structure. Contains the associated dimensions and the
    multiplier with which to convert values.

SourceFiles
    unitConversion.C
    unitConversionIO.C
    unitConversionI.H

\*---------------------------------------------------------------------------*/

#ifndef unitConversion_H
#define unitConversion_H

#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;
template<class Type> class List;
template<class Type> class Field;

// Forward declaration of friend functions and operators
class unitConversion;
unitConversion pow(const unitConversion&, const scalar);
const unitConversion& operator+(const unitConversion&, const unitConversion&);
unitConversion operator*(const unitConversion&, const unitConversion&);
unitConversion operator/(const unitConversion&, const unitConversion&);
Istream& operator>>(Istream&, unitConversion&);
Ostream& operator<<(Ostream&, const unitConversion&);
Ostream& operator<<(Ostream&, const InfoProxy<unitConversion>&);

/*---------------------------------------------------------------------------*\
                          Class unitConversion Declaration
\*---------------------------------------------------------------------------*/

class unitConversion
{
public:

    // Member constants

        //- Define an enumeration for the number of dimensionless units
        enum
        {
            nDimlessUnits = 2
        };

        //- Define an enumeration for the names of the dimensionless unit
        //  exponents
        enum dimlessUnitType
        {
            FRACTION,           // fraction   %
            ANGLE               // angle      rad, rot, deg
        };

        //- Names of the dimensionless units
        static const NamedEnum<dimlessUnitType, 2> dimlessUnitTypeNames_;


    // Static Data Members

        //- A small exponent with which to perform inexact comparisons
        static const scalar smallExponent;


private:

    // Private Data

        //- The dimensions
        dimensionSet dimensions_;

        //- Array of dimensionless unit exponents
        scalar exponents_[nDimlessUnits];

        //- The conversion multiplier with which to multiply quantities
        //  in this unit in order to recover the quantity in standard units
        scalar multiplier_;


    // Private Member Functions

        //- Compare two unit conversions and return if they are the same.
        //  An additional control determines whether the multiplier is compared
        //  as well as the dimensions and dimensionless units.
        static bool compare
        (
            const unitConversion&,
            const unitConversion&,
            const bool compareMultiplier
        );


public:

    // Static Data Members

        //- Run time type information
        ClassName("unitConversion");


    // Constructors

        //- Construct from components
        unitConversion
        (
            const dimensionSet&,
            const scalar fraction,
            const scalar angle,
            const scalar multiplier
        );

        //- Construct from a dimension set. No dimensionless units. Unity
        //  multiplier.
        unitConversion(const dimensionSet&);

        //- Copy constructor
        unitConversion(const unitConversion&) = default;

        //- Move constructor
        unitConversion(unitConversion&&) = default;

        //- Construct from stream
        unitConversion(Istream& is);


    // Member Functions

        //- Access the dimensions
        inline const dimensionSet& dimensions() const;

        //- Convert a value to standard units
        template<class T>
        T toStandard(const T&) const;

        //- Convert a pair of values to standard units
        template<class T>
        Pair<T> toStandard(const Pair<T>&) const;

        //- Convert a list of values to standard units
        template<class T>
        List<T> toStandard(const List<T>&) const;

        //- Convert a field of values to standard units
        template<class T>
        tmp<Field<T>> toStandard(const Field<T>&) const;

        //- Convert a tmp field of values to standard units
        template<class T>
        tmp<Field<T>> toStandard(const tmp<Field<T>>&) const;

        //- Convert a value to standard units
        template<class T>
        void makeStandard(T&) const;

        //- Convert a pair of values to standard units
        template<class T>
        void makeStandard(Pair<T>&) const;

        //- Convert a list of values to standard units
        template<class T>
        void makeStandard(List<T>&) const;

        //- Convert a value to user units
        template<class T>
        T toUser(const T&) const;

        //- Convert a pair of values to user units
        template<class T>
        Pair<T> toUser(const Pair<T>&) const;

        //- Convert a list of values to user units
        template<class T>
        List<T> toUser(const List<T>&) const;

        //- Convert a field of values to user units
        template<class T>
        tmp<Field<T>> toUser(const Field<T>&) const;

        //- Convert a tmp field of values to user units
        template<class T>
        tmp<Field<T>> toUser(const tmp<Field<T>>&) const;

        //- Return whether this is the "any" unit. I.e., the case where
        //  dimensions and dimensionless units are not checked, and any
        //  conversion is permitted.
        inline bool any() const;

        //- Return whether this is the "none" unit. I.e., the case where unit
        //  conversions are prohibited.
        inline bool none() const;

        //- Return whether this unit is standard. I.e., is its multiplier one?
        inline bool standard() const;

        //- Reset the unit conversion
        void reset(const unitConversion&);

        //- Update
        void read(const word& keyword, const dictionary&);

        //- Update
        void read(Istream& is);

        //- Update
        void read(const word& keyword, const dictionary&, Istream& is);

        //- Update if found in the dictionary
        bool readIfPresent(const word& keyword, const dictionary&);

        //- Update if found on the stream
        bool readIfPresent(Istream& is);

        //- Update if found on the dictionary stream
        bool readIfPresent(const word& keyword, const dictionary&, Istream& is);

        //- Return info proxy
        inline InfoProxy<unitConversion> info() const
        {
            return *this;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const unitConversion&) = delete;

        //- Disallow default bitwise move assignment
        void operator=(const unitConversion&&) = delete;


    // Friend Functions

        //- Raise to a power
        friend unitConversion pow(const unitConversion&, const scalar);


    // Friend Operators

        //- Combine
        friend const unitConversion& operator+
        (
            const unitConversion&,
            const unitConversion&
        );

        //- Multiply
        friend unitConversion operator*
        (
            const unitConversion&,
            const unitConversion&
        );

        //- Divide
        friend unitConversion operator/
        (
            const unitConversion&,
            const unitConversion&
        );


    // IOstream Operators

        //- Read from stream
        friend Istream& operator>>(Istream&, unitConversion&);

        //- Write to stream
        friend Ostream& operator<<(Ostream&, const unitConversion&);

        //- Write info to stream
        friend Ostream& operator<<(Ostream&, const InfoProxy<unitConversion>&);
};


// Global Functions

//- Write a type with a given unit conversion
template<class Type>
void writeEntry(Ostream& os, const unitConversion&, const Type& t);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "unitConversionI.H"

#ifdef NoRepository
    #include "unitConversionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
