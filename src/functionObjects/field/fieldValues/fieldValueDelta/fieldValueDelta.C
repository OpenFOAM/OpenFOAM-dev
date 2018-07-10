/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "fieldValueDelta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(fieldValueDelta, 0);
    addToRunTimeSelectionTable(functionObject, fieldValueDelta, dictionary);
}
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::fieldValueDelta::operationType,
    5
>::names[] =
{
    "add",
    "subtract",
    "min",
    "max",
    "average"
};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldValues::fieldValueDelta::operationType,
    5
> Foam::functionObjects::fieldValues::fieldValueDelta::operationTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::fieldValueDelta::writeFileHeader
(
    const label i
)
{
    const wordList& fields1 = region1Ptr_->fields();
    const wordList& fields2 = region2Ptr_->fields();

    DynamicList<word> commonFields(fields1.size());
    forAll(fields1, fieldi)
    {
        label index = findIndex(fields2, fields1[fieldi]);
        if (index != -1)
        {
            commonFields.append(fields1[fieldi]);
        }
    }

    Ostream& os = file();

    writeHeaderValue(os, "Source1", region1Ptr_->name());
    writeHeaderValue(os, "Source2", region2Ptr_->name());
    writeHeaderValue(os, "Operation", operationTypeNames_[operation_]);
    writeCommented(os, "Time");

    forAll(commonFields, i)
    {
        os  << tab << commonFields[i];
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::fieldValueDelta::fieldValueDelta
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    operation_(opSubtract),
    region1Ptr_(nullptr),
    region2Ptr_(nullptr)
{
    read(dict);
    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::fieldValueDelta::~fieldValueDelta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::fieldValueDelta::read
(
    const dictionary& dict
)
{
    regionFunctionObject::read(dict);

    region1Ptr_.reset
    (
        fieldValue::New
        (
            name() + ".region1",
            obr_,
            dict.subDict("region1"),
            false
        ).ptr()
    );
    region2Ptr_.reset
    (
        fieldValue::New
        (
            name() + ".region2",
            obr_,
            dict.subDict("region2"),
            false
        ).ptr()
    );

    operation_ = operationTypeNames_.read(dict.lookup("operation"));

    return true;
}


bool Foam::functionObjects::fieldValues::fieldValueDelta::write()
{
    logFiles::write();

    region1Ptr_->write();
    region2Ptr_->write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    Log << type() << " " << name() << " write:" << endl;

    bool found = false;
    processFields<scalar>(found);
    processFields<vector>(found);
    processFields<sphericalTensor>(found);
    processFields<symmTensor>(found);
    processFields<tensor>(found);

    if (Pstream::master())
    {
        file()<< endl;
    }

    if (!found)
    {
        Log << "    none" << endl;
    }
    else
    {
        Log << endl;
    }

    return true;
}


bool Foam::functionObjects::fieldValues::fieldValueDelta::execute()
{
    return true;
}


// ************************************************************************* //
