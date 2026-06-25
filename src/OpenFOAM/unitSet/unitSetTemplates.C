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

\*---------------------------------------------------------------------------*/

#include "scalable.H"
#include "unitSet.H"
#include "typeName.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
T Foam::unitSet::toStandard(const T& t) const
{
    return standard() ? t : t*multiplier_;
}


template<class T>
void Foam::unitSet::makeStandard(T& t) const
{
    if (!standard())
    {
        t *= multiplier_;
    }
}


template<class T>
T Foam::unitSet::toUser(const T& t) const
{
    return standard() ? t : t/multiplier_;
}


template<class T>
void Foam::unitSet::makeUser(T& t) const
{
    if (!standard())
    {
        t *= 1/multiplier_;
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, class Convert>
Foam::enableIfScalarCmptType<Type> Foam::convert
(
    Type& t,
    const unitSet& units,
    const Convert& convert
)
{
    convert(t, units);
}


template<class Type>
const typename Foam::typeUnitsType<Type>::type&
Foam::typeUnits(const unitSet& units)
{
    return units;
}


template<class Type>
typename Foam::typeUnitsType<Type>::type
Foam::typeUnits(const dimensionSet& dimensions)
{
    return dimensions;
}


template<>
inline const typename Foam::typeUnitsType<Foam::label>::type&
Foam::typeUnits<Foam::label>(const unitSet& units)
{
    typeUnits<label>(units.dimensions());

    if (!units.standard())
    {
        FatalErrorInFunction
            << "Unit conversions are not supported for "
            << pTraits<label>::typeName << "s"
            << exit(FatalError);
    }

    static const Foam::nil nil;

    return nil;
}


template<>
inline typename Foam::typeUnitsType<Foam::label>::type
Foam::typeUnits<Foam::label>(const dimensionSet& dimensions)
{
    if (!dimensions.dimensionless())
    {
        FatalErrorInFunction
            << pTraits<label>::typeName << "s must be dimensionless"
            << exit(FatalError);
    }

    return nil();
}


template<class Type>
Type Foam::readAndConvert(Istream& is, const unitSet& defaultUnits)
{
    // Read the units if they are before the value
    unitSet units(defaultUnits);
    const bool haveUnits = units.readIfPresent(is);

    // Read the value
    Type value = Foam::read<Type>(is);

    // Read the units if they are after the value
    if (!haveUnits && !is.eof())
    {
        units.readIfPresent(is);
    }

    // Modify the value by the unit conversion
    convert(value, units, unitSet::makeStandardOp());

    return value;
}


template<class Type>
Type Foam::readAndConvert(Istream& is, const dimensionSet& dimensions)
{
    return readAndConvert<Type>(is, unitSet(dimensions));
}


template<class Type, class Units>
Type Foam::readAndConvert(Istream& is, const Units& defaultUnits)
{
    // Read the value
    Type value = Foam::read<Type>(is);

    // If there is more then read units and use them to convert the value
    if (!is.eof())
    {
        Units units(defaultUnits);
        units.read(is);
        convert(value, units, unitSet::makeStandardOp());
    }
    // Otherwise convert the value using the default units
    else
    {
        convert(value, defaultUnits, unitSet::makeStandardOp());
    }

    return value;
}


template<class Type>
Type Foam::readAndConvert(Istream& is, const nil&)
{
    auto error = [&]()
    {
        FatalIOErrorInFunction(is)
            << "Unit conversion is not supported for entries of type "
            << typeName<Type>() << abort(FatalIOError);
    };

    unitSet units(unitSet::newAny());

    if (units.readIfPresent(is)) error();

    const Type value = Foam::read<Type>(is);

    if (!is.eof() && units.readIfPresent(is)) error();

    return value;
}


template<class Type>
Foam::enableIfScalable<Type, Type> Foam::readAndMaybeConvert(Istream& is)
{
    return readAndConvert<Type>(is, unitSet::newAny());
}


template<class Type>
Foam::enableIfNotScalable<Type, Type> Foam::readAndMaybeConvert(Istream& is)
{
    return readAndConvert<Type>(is, nil());
}


template<>
inline Foam::unitSet Foam::readAndMaybeConvert<Foam::unitSet>
(
    Istream& is
)
{
    return unitSet(is);
}


template<>
inline Foam::dimensionSet Foam::readAndMaybeConvert<Foam::dimensionSet>
(
    Istream& is
)
{
    return dimensionSet(is);
}


template<class Type>
void Foam::writeEntry(Ostream& os, const unitSet& defaultUnits, const Type& t)
{
    if (defaultUnits.standard())
    {
        writeEntry(os, t);
    }
    else
    {
        Type tUser(t);
        convert(tUser, defaultUnits, unitSet::makeUserOp());
        return writeEntry(os, tUser);
    }
}


template<class Type>
void Foam::writeEntry(Ostream& os, const nil&, const Type& t)
{
    writeEntry(os, t);
}


// ************************************************************************* //
