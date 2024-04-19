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

\*---------------------------------------------------------------------------*/

#include "Product2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, Foam::direction rank>
Foam::Function2s::ProductFunction1s<Type, rank>::ProductFunction1s
(
    const unitConversions& units,
    const dictionary& dict,
    const Pair<Tuple2<word, label>>& typeAndRanks
)
:
    ProductFunction1s<Type, rank - 1>(units, dict, typeAndRanks)
{
    forAll(fs, i)
    {
        if (typeAndRanks[i].second() == rank)
        {
            fs[i] =
                function1Type::New
                (
                    valueName(i, typeAndRanks[i]),
                    {i ? units.y : units.x, unitAny},
                    dict
                );
        }
    }
}


template<class Type>
Foam::Function2s::ProductFunction1s<Type, 0>::ProductFunction1s
(
    const unitConversions& units,
    const dictionary& dict,
    const Pair<Tuple2<word, label>>& typeAndRanks
)
{
    forAll(fs, i)
    {
        if (typeAndRanks[i].second() == 0)
        {
            fs[i] =
                Function1<scalar>::New
                (
                    valueName(i, typeAndRanks[i]),
                    {i ? units.y : units.x, unitAny},
                    dict
                );
        }
    }
}


template<class Type, Foam::direction rank>
Foam::Function2s::ProductFunction1s<Type, rank>::ProductFunction1s
(
    const ProductFunction1s<Type, rank>& p2f1s
)
:
    ProductFunction1s<Type, rank - 1>(p2f1s),
    fs
    (
        autoPtr<function1Type>(p2f1s.fs.first(), false),
        autoPtr<function1Type>(p2f1s.fs.second(), false)
    )
{}


template<class Type>
Foam::Function2s::ProductFunction1s<Type, 0>::ProductFunction1s
(
    const ProductFunction1s<Type, 0>& p2f1s
)
:
    fs
    (
        autoPtr<Function1<scalar>>(p2f1s.fs.first(), false),
        autoPtr<Function1<scalar>>(p2f1s.fs.second(), false)
    )
{}


template<class Type>
Foam::Function2s::Product<Type>::Product
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction2<Type, Product<Type>>(name),
    fs_(units, dict, lookupValueTypeAndRanks<Type>(dict))
{}


template<class Type>
Foam::Function2s::Product<Type>::Product(const Product<Type>& se)
:
    FieldFunction2<Type, Product<Type>>(se),
    fs_(se.fs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Product<Type>::~Product()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, Foam::direction rank>
void Foam::Function2s::ProductFunction1s<Type, rank>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    ProductFunction1s<Type, rank - 1>::write(os, units);

    forAll(fs, i)
    {
        if (fs[i].valid())
        {
            writeEntry(os, {i ? units.y : units.x, unitAny}, fs[i]());
        }
    }
}


template<class Type>
void Foam::Function2s::ProductFunction1s<Type, 0>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    forAll(fs, i)
    {
        if (fs[i].valid())
        {
            writeEntry(os, {i ? units.y : units.x, unitAny}, fs[i]());
        }
    }
}


template<class Type>
void Foam::Function2s::Product<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    fs_.write(os, units);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, class ValueType>
void Foam::Function2s::lookupValueTypeAndRank
(
    const dictionary& dict,
    const direction argument,
    Tuple2<word, label>& typeAndRank,
    label& found
)
{
    if (dict.found(valueName<ValueType>(argument)))
    {
        if (found != -1)
        {
            FatalIOErrorInFunction(dict)
                << "Multiple " << valueName(argument) << " and/or "
                << valueName(argument, "Type") << "-s specified"
                << exit(FatalIOError);
        }

        typeAndRank =
            Tuple2<word, label>
            (
                pTraits<ValueType>::typeName,
                pTraits<ValueType>::rank
            );

        found = ProductValueTypeIsValid<Type, ValueType>::value;
    }
}


template<class Type>
Foam::Tuple2<Foam::word, Foam::label> Foam::Function2s::lookupValueTypeAndRank
(
    const dictionary& dict,
    const direction argument
)
{
    Tuple2<word, label> typeAndRank(word::null, -1);
    label found = dict.found(valueName(argument)) ? 1 : -1;

    #define LOOKUP_VALUE_TYPE_AND_RANK(ValueType, nullArg)                     \
        lookupValueTypeAndRank<Type, ValueType>                                \
        (                                                                      \
            dict,                                                              \
            argument,                                                          \
            typeAndRank,                                                       \
            found                                                              \
        );
    FOR_ALL_FIELD_TYPES(LOOKUP_VALUE_TYPE_AND_RANK);
    #undef LOOKUP_VALUE_TYPE_AND_RANK

    if (found == -1)
    {
        FatalIOErrorInFunction(dict)
            << "Function " << valueName(argument)
            << " undefined in dictionary " << dict.name()
            << exit(FatalIOError);
    }

    if (found == 0)
    {
        FatalIOErrorInFunction(dict)
            << "Function " << valueName(argument, typeAndRank)
            << " returns a type that cannot be used to produce a product"
            << " of type " << pTraits<Type>::typeName
            << exit(FatalIOError);
    }

    return typeAndRank;
}


template<class Type>
Foam::Pair<Foam::Tuple2<Foam::word, Foam::label>>
Foam::Function2s::lookupValueTypeAndRanks(const dictionary& dict)
{
    Pair<Tuple2<word, label>> typeAndRanks
    (
        lookupValueTypeAndRank<Type>(dict, 0),
        lookupValueTypeAndRank<Type>(dict, 1)
    );

    // If this is a non-scalar type then at least one of the value entries must
    // have specified the type
    if
    (
        pTraits<Type>::rank > 0
     && typeAndRanks.first().second() == -1
     && typeAndRanks.second().second() == -1
    )
    {
        FatalIOErrorInFunction(dict)
            << "One of the functions " << valueName(0) << " and "
            << valueName(1) << " needs to specify the return type, e.g., as "
            << valueName<Type>(0) << exit(FatalIOError);
    }

    // If both types are specified then the sum of their ranks must equal the
    // rank of the function
    if
    (
        typeAndRanks.first().second() != -1
     && typeAndRanks.second().second() != -1
     && typeAndRanks.first().second()
      + typeAndRanks.second().second()
     != pTraits<Type>::rank
    )
    {
        FatalIOErrorInFunction(dict)
            << "The functions " << valueName(0, typeAndRanks.first())
            << " and " << valueName(1, typeAndRanks.second()) << " return "
            << "types for which the product is not of type "
            << pTraits<Type>::typeName << exit(FatalIOError);
    }

    // If this is a scalar type, then neither entry needs to specify the type.
    // They both must be scalars.
    if
    (
        pTraits<Type>::rank == 0
     && typeAndRanks.first().second() == -1
     && typeAndRanks.second().second() == -1
    )
    {
        typeAndRanks.first().second() = 0;
        typeAndRanks.second().second() = 0;
    }

    // Determine remaining unspecified ranks
    forAll(typeAndRanks, i)
    {
        if (typeAndRanks[i].second() == -1)
        {
            typeAndRanks[i].second() =
                pTraits<Type>::rank - typeAndRanks[!i].second();
        }
    }

    return typeAndRanks;
}


// ************************************************************************* //
