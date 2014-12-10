/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "addSubtract.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(addSubtract, 0);
        addToRunTimeSelectionTable(calcType, addSubtract, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::calcTypes::addSubtract::writeAddSubtractFields
(
    const Time& runTime,
    const fvMesh& mesh,
    const IOobject& baseFieldHeader
)
{
    bool processed = false;

    IOobject addSubtractFieldHeader
    (
        addSubtractFieldName_,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (addSubtractFieldHeader.headerOk())
    {
        writeAddSubtractField<scalar>
        (
            baseFieldHeader,
            addSubtractFieldHeader,
            mesh,
            processed
        );
        writeAddSubtractField<vector>
        (
            baseFieldHeader,
            addSubtractFieldHeader,
            mesh,
            processed
        );
        writeAddSubtractField<sphericalTensor>
        (
            baseFieldHeader,
            addSubtractFieldHeader,
            mesh,
            processed
        );
        writeAddSubtractField<symmTensor>
        (
            baseFieldHeader,
            addSubtractFieldHeader,
            mesh,
            processed
        );
        writeAddSubtractField<tensor>
        (
            baseFieldHeader,
            addSubtractFieldHeader,
            mesh,
            processed
        );

        if (!processed)
        {
            FatalError
                << "Unable to process " << baseFieldName_
                << " + " << addSubtractFieldName_ << nl
                << "No call to addSubtract for fields of type "
                << baseFieldHeader.headerClassName() << " + "
                << addSubtractFieldHeader.headerClassName() << nl << nl
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn("calcTypes::addSubtract::writeAddSubtractFields()")
            << "Unable to read addSubtract field: " << addSubtractFieldName_
            << nl << exit(FatalError);
    }
}


void Foam::calcTypes::addSubtract::writeAddSubtractValues
(
    const Time& runTime,
    const fvMesh& mesh,
    const IOobject& baseFieldHeader
)
{
    bool processed = false;

    writeAddSubtractValue<scalar>
    (
        baseFieldHeader,
        addSubtractValueStr_,
        mesh,
        processed
    );
    writeAddSubtractValue<vector>
    (
        baseFieldHeader,
        addSubtractValueStr_,
        mesh,
        processed
    );
    writeAddSubtractValue<sphericalTensor>
    (
        baseFieldHeader,
        addSubtractValueStr_,
        mesh,
        processed
    );
    writeAddSubtractValue<symmTensor>
    (
        baseFieldHeader,
        addSubtractValueStr_,
        mesh,
        processed
    );
    writeAddSubtractValue<tensor>
    (
        baseFieldHeader,
        addSubtractValueStr_,
        mesh,
        processed
    );

    if (!processed)
    {
        FatalErrorIn("calcTypes::addSubtract::writeAddSubtractValue()")
            << "Unable to process " << baseFieldName_
            << " + " << addSubtractValueStr_ << nl
            << "No call to addSubtract for fields of type "
            << baseFieldHeader.headerClassName() << nl << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::addSubtract::addSubtract()
:
    calcType(),
    baseFieldName_(""),
    calcType_(FIELD),
    addSubtractFieldName_(""),
    addSubtractValueStr_(""),
    resultName_(""),
    calcMode_(ADD)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::addSubtract::~addSubtract()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::addSubtract::init()
{
    argList::validArgs.append("add");
    argList::validArgs.append("baseField");
    argList::validArgs.append("calcMode");
    argList::validOptions.insert("field", "fieldName");
    argList::validOptions.insert("value", "valueString");
    argList::validOptions.insert("resultName", "fieldName");
}


void Foam::calcTypes::addSubtract::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    baseFieldName_ = args[2];
    const word calcModeName = args[3];

    if (calcModeName == "add")
    {
        calcMode_ = ADD;
    }
    else if (calcModeName == "subtract")
    {
        calcMode_ = SUBTRACT;
    }
    else
    {
        FatalErrorIn("calcTypes::addSubtract::preCalc")
            << "Invalid calcMode: " << calcModeName << nl
            << "    Valid calcModes are add and subtract" << nl
            << exit(FatalError);
    }

    if (args.optionReadIfPresent("field", addSubtractFieldName_))
    {
        calcType_ = FIELD;
    }
    else if (args.optionReadIfPresent("value", addSubtractValueStr_))
    {
        calcType_ = VALUE;
    }
    else
    {
        FatalErrorIn("calcTypes::addSubtract::preCalc")
            << "addSubtract requires either -field or -value option"
            << nl << exit(FatalError);
    }

    args.optionReadIfPresent("resultName", resultName_);
}


void Foam::calcTypes::addSubtract::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    IOobject baseFieldHeader
    (
        baseFieldName_,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (baseFieldHeader.headerOk())
    {
        switch (calcType_)
        {
            case FIELD:
            {
                writeAddSubtractFields(runTime, mesh, baseFieldHeader);
                break;
            }
            case VALUE:
            {
                writeAddSubtractValues(runTime, mesh, baseFieldHeader);
                break;
            }
            default:
            {
                FatalErrorIn("calcTypes::addSubtract::calc")
                    << "unknown calcType " << calcType_ << nl
                    << abort(FatalError);
            }
        }
    }
    else
    {
        FatalErrorIn("calcTypes::addSubtract::calc")
            << "Unable to read base field: " << baseFieldName_
            << nl << exit(FatalError);
    }
}


// ************************************************************************* //
