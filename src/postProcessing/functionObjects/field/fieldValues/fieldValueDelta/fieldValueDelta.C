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
#include "ListOps.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(fieldValueDelta, 0);
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
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFiles(obr, name, typeName),
    name_(name),
    obr_(obr),
    loadFromFiles_(loadFromFiles),
    log_(true),
    operation_(opSubtract),
    source1Ptr_(NULL),
    source2Ptr_(NULL)
{
    if (!isA<fvMesh>(obr))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::fieldValueDelta::~fieldValueDelta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::fieldValues::fieldValueDelta::read
(
    const dictionary& dict
)
{
    log_ = dict.lookupOrDefault<Switch>("log", true);
    source1Ptr_.reset
    (
        fieldValue::New
        (
            name_ + ".source1",
            obr_,
            dict.subDict("source1"),
            loadFromFiles_,
            false
        ).ptr()
    );
    source2Ptr_.reset
    (
        fieldValue::New
        (
            name_ + ".source2",
            obr_,
            dict.subDict("source2"),
            loadFromFiles_,
            false
        ).ptr()
    );

    operation_ = operationTypeNames_.read(dict.lookup("operation"));
}


void Foam::functionObjects::fieldValues::fieldValueDelta::write()
{
    functionObjectFiles::write();

    source1Ptr_->write();
    source2Ptr_->write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    if (log_) Info<< type() << " " << name_ << " output:" << endl;

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

    if (log_)
    {
        if (!found)
        {
            Info<< "    none" << endl;
        }
        else
        {
            Info<< endl;
        }
    }
}


void Foam::functionObjects::fieldValues::fieldValueDelta::execute()
{}


void Foam::functionObjects::fieldValues::fieldValueDelta::end()
{}


void Foam::functionObjects::fieldValues::fieldValueDelta::timeSet()
{}


void Foam::functionObjects::fieldValues::fieldValueDelta::updateMesh
(
    const mapPolyMesh&
)
{}


void Foam::functionObjects::fieldValues::fieldValueDelta::movePoints
(
    const polyMesh&
)
{}


// ************************************************************************* //
