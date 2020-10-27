/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    13
>::names[] =
{
    "none",
    "sum",
    "weightedSum",
    "sumMag",
    "average",
    "weightedAverage",
    "volAverage",
    "weightedVolAverage",
    "volIntegrate",
    "weightedVolIntegrate",
    "min",
    "max",
    "CoV"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::volFieldValue::operationType,
    13
> Foam::functionObjects::fieldValues::volFieldValue::operationTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::volFieldValue::initialise
(
    const dictionary& dict
)
{
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

    Info<< nl << endl;
}


void Foam::functionObjects::fieldValues::volFieldValue::writeFileHeader
(
    const label i
)
{
    volRegion::writeFileHeader(*this, file());

    writeCommented(file(), "Time");

    forAll(fields_, fieldi)
    {
        file() << tab << operationTypeNames_[operation_] << "(";

        forAll(weightFieldNames_, i)
        {
            file() << weightFieldNames_[i] << ',';
        }

        file() << fields_[fieldi] << ")";
    }

    file() << endl;
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
    volRegion(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.read(dict.lookup("operation")))
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
    volRegion(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.read(dict.lookup("operation")))
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

    forAll(fields_, i)
    {
        const word& fieldName = fields_[i];
        bool processed = false;

        #define processType(fieldType, none)                                   \
            processed = processed || writeValues<fieldType>(fieldName);
        FOR_ALL_FIELD_TYPES(processType)

        if (!processed)
        {
            cannotFindObject(fieldName);
        }
    }

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
