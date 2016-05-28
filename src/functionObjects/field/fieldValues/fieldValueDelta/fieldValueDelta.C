/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    const wordList& fields1 = source1Ptr_->fields();
    const wordList& fields2 = source2Ptr_->fields();

    DynamicList<word> commonFields(fields1.size());
    forAll(fields1, i)
    {
        label index = findIndex(fields2, fields1[i]);
        if (index != -1)
        {
            commonFields.append(fields1[i]);
        }
    }

    Ostream& os = file();

    writeHeaderValue(os, "Source1", source1Ptr_->name());
    writeHeaderValue(os, "Source2", source2Ptr_->name());
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
    writeFiles(name, runTime, dict, name),
    operation_(opSubtract),
    source1Ptr_(NULL),
    source2Ptr_(NULL)
{
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

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
    writeFiles::read(dict);

    source1Ptr_.reset
    (
        fieldValue::New
        (
            name() + ".source1",
            obr_,
            dict.subDict("source1"),
            false
        ).ptr()
    );
    source2Ptr_.reset
    (
        fieldValue::New
        (
            name() + ".source2",
            obr_,
            dict.subDict("source2"),
            false
        ).ptr()
    );

    operation_ = operationTypeNames_.read(dict.lookup("operation"));

    return true;
}


bool Foam::functionObjects::fieldValues::fieldValueDelta::write
(
    const bool postProcess
)
{
    writeFiles::write();

    source1Ptr_->write();
    source2Ptr_->write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    Log << type() << " " << name() << " output:" << endl;

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


bool Foam::functionObjects::fieldValues::fieldValueDelta::execute
(
    const bool postProcess
)
{
    return true;
}


// ************************************************************************* //
