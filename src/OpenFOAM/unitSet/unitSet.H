/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2026 OpenFOAM Foundation
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
    Foam::unitSet

Description
    Unit conversion structure. Contains the associated dimensions and the
    multiplier with which to convert values.

SourceFiles
    unitSet.C
    unitSetIO.C
    unitSetI.H

\*---------------------------------------------------------------------------*/

#ifndef unitSet_H
#define unitSet_H

#include "dimensionSet.H"
#include "nil.H"
#include "scalable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class unitSet;
unitSet pow(const unitSet&, const scalar);
const unitSet& operator+(const unitSet&, const unitSet&);
unitSet operator*(const unitSet&, const unitSet&);
unitSet operator/(const unitSet&, const unitSet&);
Istream& operator>>(Istream&, unitSet&);
Ostream& operator<<(Ostream&, const unitSet&);
Ostream& operator<<(Ostream&, const InfoProxy<unitSet>&);

/*---------------------------------------------------------------------------*\
                          Class unitSet Declaration
\*---------------------------------------------------------------------------*/

class unitSet
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
        static const NamedEnum<dimlessUnitType, 2>& dimlessUnitTypeNames_;


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
            const unitSet&,
            const unitSet&,
            const bool compareMultiplier
        );


public:

    // Static Data Members

        //- Run time type information
        ClassName("unitSet");


    // Static Member Functions

        //- Return a new "unitless" unit set
        inline static unitSet newUnitless();

        //- Return a new "any" unit set
        inline static unitSet newAny();

        //- Return a new "none" unit set
        inline static unitSet newNone();


    // Constructors

        //- Construct from components
        unitSet
        (
            const dimensionSet&,
            const scalar fraction,
            const scalar angle,
            const scalar multiplier
        );

        //- Construct from a dimension set. No dimensionless units. Unity
        //  multiplier.
        unitSet(const dimensionSet&);

        //- Copy constructor
        unitSet(const unitSet&) = default;

        //- Move constructor
        unitSet(unitSet&&) = default;

        //- Construct from stream
        unitSet(Istream& is);


    // Member Functions

        //- Access the dimensions
        inline const dimensionSet& dimensions() const;

        //- Return true if it is dimensionless
        inline bool dimensionless() const;

        //- Convert a value to standard units
        template<class T>
        T toStandard(const T&) const;

        //- Convert a value to standard units
        template<class T>
        void makeStandard(T&) const;

        //- Convert a value to user units
        template<class T>
        T toUser(const T&) const;

        //- Convert a value to user units
        template<class T>
        void makeUser(T&) const;

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
        void reset(const unitSet&);

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
        inline InfoProxy<unitSet> info() const
        {
            return *this;
        }


    // Member Operators

        //- Access the dimension exponents
        inline scalar operator[](const dimensionSet::dimensionType type) const;

        //- Access the dimensionless unit exponents
        inline scalar operator[](const dimlessUnitType type) const;

        //- Disallow default bitwise assignment
        void operator=(const unitSet&) = delete;

        //- Disallow default bitwise move assignment
        void operator=(const unitSet&&) = delete;


    // Friend Functions

        //- Raise to a power
        friend unitSet pow(const unitSet&, const scalar);


    // Friend Operators

        //- Combine
        friend const unitSet& operator+
        (
            const unitSet&,
            const unitSet&
        );

        //- Multiply
        friend unitSet operator*
        (
            const unitSet&,
            const unitSet&
        );

        //- Divide
        friend unitSet operator/
        (
            const unitSet&,
            const unitSet&
        );


    // IOstream Operators

        //- Read from stream
        friend Istream& operator>>(Istream&, unitSet&);

        //- Write to stream
        friend Ostream& operator<<(Ostream&, const unitSet&);

        //- Write info to stream
        friend Ostream& operator<<(Ostream&, const InfoProxy<unitSet>&);


    // Public Types

        //- Functor to convert to standard
        struct makeStandardOp
        {
            template<class Type>
            inline void operator()(Type& t, const unitSet& units) const
            {
                units.makeStandard(t);
            }
        };

        //- Functor to convert to user
        struct makeUserOp
        {
            template<class Type>
            inline void operator()(Type& t, const unitSet& units) const
            {
                units.makeUser(t);
            }
        };
};


// Global Functions

//- Convert a primitive type
template<class Type, class Convert>
enableIfScalarCmptType<Type> convert
(
    Type& t,
    const unitSet& units,
    const Convert&
);

//- Get the units type for a given type. This returns the unit set for
//  everything except label, for which it returns nil. This can be used to
//  modify the type of the units passed to dictionary, such that a single
//  template instantiation can be used for label as well as the field types.
//  The label specialisation checks the unit set provided to make sure it is
//  both dimensionless and standard.
template<class Type>
struct typeUnitsType { typedef unitSet type; };

template<>
struct typeUnitsType<label> { typedef nil type; };

template<class Type>
const typename typeUnitsType<Type>::type& typeUnits
(
    const unitSet&
);

template<class Type>
typename typeUnitsType<Type>::type typeUnits
(
    const dimensionSet&
);

template<>
inline const typename typeUnitsType<label>::type& typeUnits<label>
(
    const unitSet&
);

template<>
inline typename typeUnitsType<label>::type typeUnits<label>
(
    const dimensionSet&
);

//- Read a type which supports unit conversion
template<class Type>
Type readAndConvert(Istream&, const unitSet&);

//- Read a type which supports unit conversion
template<class Type>
Type readAndConvert(Istream&, const dimensionSet&);

//- Read a type which supports unit conversion
template<class Type, class Units>
Type readAndConvert(Istream&, const Units&);

//- Read a type which does not support unit conversion
template<class Type>
Type readAndConvert(Istream&, const nil&);

//- Read a type which may or may not support unit conversion
template<class Type>
enableIfScalable<Type, Type> readAndMaybeConvert(Istream&);

//- Read a type which may or may not support unit conversion
template<class Type>
enableIfNotScalable<Type, Type> readAndMaybeConvert(Istream&);

//- Read a unit set
template<>
inline unitSet readAndMaybeConvert<unitSet>(Istream&);

//- Read a dimension set
template<>
inline dimensionSet readAndMaybeConvert<dimensionSet>(Istream&);

//- Write a type with a given unit conversion
template<class Type>
void writeEntry(Ostream&, const unitSet&, const Type&);

//- Write a type without unit conversion
template<class Type>
void writeEntry(Ostream&, const nil&, const Type&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "unitSetI.H"

#ifdef NoRepository
    #include "unitSetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
