/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "volFieldValue.H"
#include "fvMesh.H"
#include "volFields.H"
#include "writeFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(volFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, volFieldValue, dictionary);
    addToRunTimeSelectionTable(functionObject, volFieldValue, dictionary);
}
}
}

template<>
const char*
Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volFieldValue::operationType,
    11
>::names[] =
{
    "none",
    "sum",
    "sumMag",
    "average",
    "volAverage",
    "volIntegrate",
    "min",
    "max",
    "minMag",
    "maxMag",
    "CoV"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volFieldValue::operationType,
    11
> Foam::functionObjects::fieldValues::volFieldValue::operationTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::volFieldValue::initialise
(
    const dictionary& dict
)
{
    dict.readIfPresent<Switch>("writeLocation", writeLocation_);

    if (dict.readIfPresent("weightFields", weightFieldNames_))
    {
        Info<< name() << " " << operationTypeNames_[operation_]
            << " weight fields " << weightFieldNames_;
    }
    else if (dict.found("weightField"))
    {
        weightFieldNames_.setSize(1);
        dict.lookup("weightField") >> weightFieldNames_[0];

        Info<< name() << " " << operationTypeNames_[operation_]
            << " weight field " << weightFieldNames_[0];
    }

    if (dict.readIfPresent("scaleFactor", scaleFactor_))
    {
        Info<< "    scale factor = " << scaleFactor_ << nl;
    }

    Info<< nl << endl;
}


template<class Type>
void Foam::functionObjects::fieldValues::volFieldValue::
writeFileHeaderLocation()
{
    switch (operation_)
    {
        case operationType::minMag:
        case operationType::maxMag:
            file() << tab << "location" << tab << "cell";
            if (Pstream::parRun())
            {
                file() << tab << "processor";
            }
            break;
        default:
            break;
    }
}


template<>
void Foam::functionObjects::fieldValues::volFieldValue::
writeFileHeaderLocation<Foam::scalar>()
{
    switch (operation_)
    {
        case operationType::min:
        case operationType::max:
        case operationType::minMag:
        case operationType::maxMag:
            file() << tab << "location" << tab << "cell";
            if (Pstream::parRun())
            {
                file() << tab << "processor";
            }
            break;
        default:
            break;
    }
}


void Foam::functionObjects::fieldValues::volFieldValue::writeFileHeader
(
    const label i
)
{
    fvCellSet::writeFileHeader(*this, file());

    writeCommented(file(), "Time");

    forAll(fields_, fieldi)
    {
        file() << tab << operationTypeNames_[operation_] << "(";

        forAll(weightFieldNames_, i)
        {
            file() << weightFieldNames_[i] << ',';
        }

        file() << fields_[fieldi] << ")";

        if (writeLocation_)
        {
            #define writeFileHeaderLocationFieldType(fieldType, none)          \
                if (validField<fieldType>(fields_[fieldi]))                    \
                {                                                              \
                    writeFileHeaderLocation<fieldType>();                      \
                }
            FOR_ALL_FIELD_TYPES(writeFileHeaderLocationFieldType)
            #undef writeHeaderLocationFieldType
        }
    }

    file() << endl;
}


bool Foam::functionObjects::fieldValues::volFieldValue::processValues
(
    const Field<scalar>& values,
    const scalarField& weights,
    const scalarField& V,
    Result<scalar>& result
) const
{
    switch (operation_)
    {
        case operationType::min:
        {
            compareScalars(values, vGreat, result, lessOp<scalar>());
            return true;
        }
        case operationType::minMag:
        {
            compareScalars(mag(values), vGreat, result, lessOp<scalar>());
            return true;
        }
        case operationType::max:
        {
            compareScalars(values, -vGreat, result, greaterOp<scalar>());
            return true;
        }
        case operationType::maxMag:
        {
            compareScalars(mag(values), -vGreat, result, greaterOp<scalar>());
            return true;
        }
        default:
        {
            // Fall through to same-type operations
            return processValuesTypeType(values, weights, V, result);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    fvCellSet(fieldValue::mesh_, dict),
    writeLocation_(false),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    scaleFactor_(1)
{
    read(dict);
}


Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    fvCellSet(fieldValue::mesh_, dict),
    writeLocation_(false),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    scaleFactor_(1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::~volFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::volFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);

    // No additional info to read
    initialise(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::volFieldValue::write()
{
    fieldValue::write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    // Construct the weight field and the volumes
    scalarField weights(nCells(), 1);
    forAll(weightFieldNames_, i)
    {
        weights *= getFieldValues<scalar>(weightFieldNames_[i]);
    }
    const scalarField V(filterField(fieldValue::mesh_.V()));

    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        bool ok = false;

        #define writeValuesFieldType(fieldType, none)                          \
            ok = ok || writeValues<fieldType>(fieldName, weights, V);
        FOR_ALL_FIELD_TYPES(writeValuesFieldType)
        #undef writeValuesFieldType

        if (!ok)
        {
            cannotFindObject(fieldName);
        }
    }

    if (operation_ != operationType::none && Pstream::master())
    {
        file() << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
