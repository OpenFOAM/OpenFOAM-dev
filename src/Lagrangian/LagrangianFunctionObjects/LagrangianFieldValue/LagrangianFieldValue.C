/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianFieldValue.H"
#include "LagrangianFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
struct ValueLocation
{
    Type value;
    FixedList<label, pTraits<Type>::nComponents> proci;
    FixedList<label, pTraits<Type>::nComponents> elementi;
    FixedList<point, pTraits<Type>::nComponents> position;
};

template<class Type>
Ostream& operator<<(Ostream& os, const ValueLocation<Type>& vl)
{
    return os
        << vl.value << token::SPACE
        << vl.proci << token::SPACE
        << vl.elementi << token::SPACE
        << vl.position;
}

template<class Type>
Istream& operator>>(Istream& is, ValueLocation<Type>& vl)
{
    return is
        >> vl.value
        >> vl.proci
        >> vl.elementi
        >> vl.position;
}

#define DefineContiguousValueLocationType(Type, nullArg)                       \
    template<>                                                                 \
    inline bool contiguous<ValueLocation<Type>>() { return true; }
DefineContiguousValueLocationType(bool, )
DefineContiguousValueLocationType(label, )
FOR_ALL_FIELD_TYPES(DefineContiguousValueLocationType)
#undef DefineContiguousValueLocationType

}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(LagrangianFieldValue, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        LagrangianFieldValue,
        dictionary
    );
}
}


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::LagrangianFieldValue::operationType,
    6
>::names[] = {"sum", "average", "min", "max", "minMag", "maxMag"};

const Foam::NamedEnum
<
    Foam::functionObjects::LagrangianFieldValue::operationType,
    6
> Foam::functionObjects::LagrangianFieldValue::operationTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::LagrangianFieldValue::readCoeffs
(
    const dictionary& dict
)
{
    // Read the fields
    const bool haveFields = dict.found("fields");
    const bool haveField = dict.found("field");
    if (haveFields == haveField)
    {
        FatalIOErrorInFunction(dict)
            << "keywords fields and field both "
            << (haveFields ? "" : "un") << "defined in "
            << "dictionary " << dict.name()
            << exit(FatalIOError);
    }
    else if (haveFields)
    {
        dict.lookup("fields") >> fields_;
    }
    else if (haveField)
    {
        fields_.resize(1);
        dict.lookup("field") >> fields_.first();
    }

    // Read the weight fields
    const bool haveWeightFields = dict.found("weightFields");
    const bool haveWeightField = dict.found("weightField");
    if (haveWeightFields && haveWeightField)
    {
        FatalIOErrorInFunction(dict)
            << "keywords weightFields and weightField both "
            << "defined in dictionary " << dict.name()
            << exit(FatalIOError);
    }
    else if (haveWeightFields)
    {
        dict.lookup("weightFields") >> weightFields_;
    }
    else if (haveWeightField)
    {
        weightFields_.resize(1);
        dict.lookup("weightField") >> weightFields_.first();
    }
    else
    {
        // No weights
        weightFields_.clear();
    }

    // Whether or not we are writing the location
    dict.readIfPresent<Switch>("writeLocation", writeLocation_);

    resetName(typeName);
}


template<class Type>
void Foam::functionObjects::LagrangianFieldValue::writeName
(
    const word& name
)
{
    static const direction nComponents = pTraits<Type>::nComponents;

    const auto& componentNames = pTraits<Type>::componentNames;

    if (Pstream::master())
    {
        for (direction d = 0; d < nComponents; ++ d)
        {
            writeTabbed
            (
                file(),
                word(operationTypeNames_[operation_])
              + "("
              + name
              + (word(componentNames[d]).empty() ? "" : "_")
              + word(componentNames[d])
              + ")"
            );
        }
    }
}


template<class Type, class LocationType>
void Foam::functionObjects::LagrangianFieldValue::writeLocationName
(
    const word& name,
    const word& locationName
)
{
    static const direction nComponents = pTraits<Type>::nComponents;
    static const direction lnComponents = pTraits<LocationType>::nComponents;

    const auto& componentNames = pTraits<Type>::componentNames;
    const auto& lcomponentNames = pTraits<LocationType>::componentNames;

    if (Pstream::master())
    {
        for (direction d = 0; d < nComponents; ++ d)
        {
            for (direction ld = 0; ld < lnComponents; ++ ld)
            {
                writeTabbed
                (
                    file(),
                    word(operationTypeNames_[operation_])
                  + "("
                  + name
                  + (word(componentNames[d]).empty() ? "" : "_")
                  + word(componentNames[d])
                  + ")"
                  + ":"
                  + locationName
                  + (word(lcomponentNames[ld]).empty() ? "" : "_")
                  + word(lcomponentNames[ld])
                );
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::LagrangianFieldValue::writeNameAndLocationNames
(
    const word& name
)
{
    writeName<Type>(name);

    if (writeLocation_)
    {
        if (Pstream::parRun())
        {
            writeLocationName<Type, label>(name, "processor");
        }
        writeLocationName<Type, label>(name, "element");
        writeLocationName<Type, point>(name, "position");
    }
}


template<class Type>
void Foam::functionObjects::LagrangianFieldValue::writeValue
(
    const Type& value
)
{
    if (Pstream::master())
    {
        for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
        {
            file() << valueWidth(1) << component(value, d);
        }
    }
}


template<class Type, class LocationType>
void Foam::functionObjects::LagrangianFieldValue::writeLocationValue
(
    const FixedList<LocationType, pTraits<Type>::nComponents>& value
)
{
    if (Pstream::master())
    {
        for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
        {
            writeValue(value[d]);
        }
    }
}


template<class Type, class Op>
void Foam::functionObjects::LagrangianFieldValue::writeValueAndLocationValues
(
    const tmp<Field<Type>>& tField,
    const scalar emptyValue,
    const Op& op
)
{
    const Field<Type>& field = tField();

    ValueLocation<Type> result;

    for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
    {
        setComponent(result.value, d) = emptyValue;
        result.proci[d] = -1;
        result.elementi[d] = -1;
        result.position[d] = point::uniform(NaN);
    }

    if (field.size())
    {
        FixedList<label, pTraits<Type>::nComponents>& ei = result.elementi;

        ei = 0;

        for (label i = 1; i < field.size(); ++ i)
        {
            for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
            {
                if (op(component(field[i], d), component(field[ei[d]], d)))
                {
                    ei[d] = i;
                }
            }
        }

        for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
        {
            setComponent(result.value, d) = component(field[ei[d]], d);
            result.proci[d] = Pstream::parRun() ? Pstream::myProcNo() : -1;
            result.elementi[d] = ei[d];
            result.position[d] = mesh().position(ei[d]);
        }
    }

    reduce
    (
        result,
        [&op](const ValueLocation<Type>& a, const ValueLocation<Type>& b)
        {
            ValueLocation<Type> result;

            for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
            {
                const ValueLocation<Type>& r =
                    op(component(a.value, d), component(b.value, d)) ? a : b;

                setComponent(result.value, d) = component(r.value, d);
                result.proci[d] = r.proci[d];
                result.elementi[d] = r.elementi[d];
                result.position[d] = r.position[d];
            }

            return result;
        }
    );

    writeValue<Type>(result.value);

    if (writeLocation_)
    {
        if (Pstream::parRun())
        {
            writeLocationValue<Type, label>(result.proci);
        }
        writeLocationValue<Type, label>(result.elementi);
        writeLocationValue<Type, point>(result.position);
    }
}


template<template<class> class GeoField>
bool Foam::functionObjects::LagrangianFieldValue::multiplyWeight
(
    const word& weightFieldName,
    scalarField& weight
) const
{
    if (!mesh().foundObject<GeoField<scalar>>(weightFieldName)) return false;

    const GeoField<scalar>& w =
        mesh().lookupObject<GeoField<scalar>>(weightFieldName);

    weight *= w;

    return true;
}


template<template<class> class GeoField, class Type>
bool Foam::functionObjects::LagrangianFieldValue::writeFieldName
(
    const word& fieldName
)
{
    if (!mesh().foundObject<GeoField<Type>>(fieldName)) return false;

    switch (operation_)
    {
        case operationType::sum:
        case operationType::average:
            writeName<Type>(fieldName);
            break;
        case operationType::min:
        case operationType::max:
            writeNameAndLocationNames<Type>(fieldName);
            break;
        case operationType::minMag:
        case operationType::maxMag:
            writeNameAndLocationNames<scalar>(fieldName);
            break;
    }

    return true;
}


template<template<class> class GeoField, class Type>
bool Foam::functionObjects::LagrangianFieldValue::writeFieldValue
(
    const scalarField& weight,
    const word& fieldName
)
{
    if (!mesh().foundObject<GeoField<Type>>(fieldName)) return false;

    const typename GeoField<Type>::FieldType& field =
        mesh().lookupObject<GeoField<Type>>(fieldName).primitiveField();

    switch (operation_)
    {
        case operationType::sum:
            writeValue(gSum(weight*field));
            break;
        case operationType::average:
            writeValue(gSum(weight*field)/gSum(weight));
            break;
        case operationType::min:
            writeValueAndLocationValues
            (
                weight*field,
                vGreat,
                lessOp<scalar>()
            );
            break;
        case operationType::max:
            writeValueAndLocationValues
            (
                weight*field,
                -vGreat,
                greaterOp<scalar>()
            );
            break;
        case operationType::minMag:
            writeValueAndLocationValues
            (
                mag(weight*field),
                vGreat,
                lessOp<scalar>()
            );
            break;
        case operationType::maxMag:
            writeValueAndLocationValues
            (
                mag(weight*field),
                -vGreat,
                greaterOp<scalar>()
            );
            break;
    }

    return true;
}


void Foam::functionObjects::LagrangianFieldValue::writeFileHeader(const label)
{
    writeHeader(file(), "Lagrangian Sum");
    writeCommented(file(), "Time");

    forAll(fields_, fieldi)
    {
        #define WRITE_FIELD_NAME(Type, GeoField) \
            && !writeFieldName<GeoField, Type>(fields_[fieldi])

        if
        (
            true
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_NAME, LagrangianField)
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_NAME, LagrangianDynamicField)
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_NAME, LagrangianInternalField)
        )
        {
            cannotFindObject(fields_[fieldi]);
        }

        #undef WRITE_COLUMN_HEADER
    }

    file().endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::LagrangianFieldValue::LagrangianFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    LagrangianMeshFunctionObject(name, runTime, dict),
    logFiles(mesh(), name),
    fields_(),
    weightFields_(),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    writeLocation_(false)
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::LagrangianFieldValue::~LagrangianFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::LagrangianFieldValue::read(const dictionary& dict)
{
    if (LagrangianMeshFunctionObject::read(dict))
    {
        readCoeffs(dict);
        return true;
    }
    else
    {
        return false;
    }
}


Foam::wordList Foam::functionObjects::LagrangianFieldValue::fields() const
{
    wordList result(fields_);
    result.append(weightFields_);
    return result;
}


bool Foam::functionObjects::LagrangianFieldValue::execute()
{
    return true;
}


bool Foam::functionObjects::LagrangianFieldValue::write()
{
    logFiles::write();

    // Filter out operations that don't make sense if there are no elements
    if (returnReduce(mesh().size(), sumOp<label>()) == 0)
    {
        switch (operation_)
        {
            case operationType::sum:
                break;
            case operationType::average:
            case operationType::min:
            case operationType::max:
            case operationType::minMag:
            case operationType::maxMag:
                return true;
        }
    }

    if (Pstream::master())
    {
        writeTime(file());
    }

    // Construct the weights
    scalarField weight(mesh().size(), 1);
    forAll(weightFields_, weightFieldi)
    {
        const word& weightFieldName = weightFields_[weightFieldi];

        if
        (
            !multiplyWeight<LagrangianField>(weightFieldName, weight)
         && !multiplyWeight<LagrangianDynamicField>(weightFieldName, weight)
         && !multiplyWeight<LagrangianInternalField>(weightFieldName, weight)
        )
        {
            FatalErrorInFunction
                << "Weight field " << weightFieldName << " was not found"
                << exit(FatalError);
        }
    }

    // Write the field values
    forAll(fields_, fieldi)
    {
        const word& fieldName = fields_[fieldi];

        #define WRITE_FIELD_VALUE(Type, GeoField) \
            && !writeFieldValue<GeoField, Type>(weight, fieldName)

        if
        (
            true
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_VALUE, LagrangianField)
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_VALUE, LagrangianDynamicField)
            FOR_ALL_FIELD_TYPES(WRITE_FIELD_VALUE, LagrangianInternalField)
        )
        {
            cannotFindObject(fields_[fieldi]);
        }

        #undef WRITE_COLUMN_VALUE
    }

    if (Pstream::master())
    {
        file().endl();
    }

    return true;
}


// ************************************************************************* //
