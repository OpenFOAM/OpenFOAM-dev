/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "codedFunctionObject.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(codedFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        codedFunctionObject,
        dictionary
    );
}


const Foam::wordList Foam::codedFunctionObject::codeKeys
{
    "codeData",
    "codeEnd",
    "codeExecute",
    "codeInclude",
    "codeRead",
    "codeFields",
    "codeWrite",
    "localCode"
};

const Foam::wordList Foam::codedFunctionObject::codeDictVars
{
    word::null,
    word::null,
    word::null,
    word::null,
    "dict",
    word::null,
    word::null,
    word::null,
};

const Foam::word Foam::codedFunctionObject::codeOptions
(
    "codedFunctionObjectOptions"
);

const Foam::wordList Foam::codedFunctionObject::compileFiles
{
    "codedFunctionObjectTemplate.C"
};

const Foam::wordList Foam::codedFunctionObject::copyFiles
{
    "codedFunctionObjectTemplate.H"
};


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::codedFunctionObject::updateLibrary(const dictionary& dict)
{
    redirectFunctionObjectPtr_.clear();

    codedBase::updateLibrary(dict);

    dictionary constructDict(dict);
    constructDict.set("type", codeName());

    redirectFunctionObjectPtr_ = functionObject::New
    (
        codeName(),
        time_,
        constructDict
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedFunctionObject::codedFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject(name, time, dict),
    codedBase
    (
        name,
        dict,
        codeKeys,
        codeDictVars,
        codeOptions,
        compileFiles,
        copyFiles
    )
{
    setFilterVariable("typeName", codeName());

    // Make verbose if debugging
    setFilterVariable("verbose", Foam::name(bool(debug)));

    updateLibrary(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::codedFunctionObject::~codedFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::functionObject& Foam::codedFunctionObject::redirectFunctionObject() const
{
    if (!redirectFunctionObjectPtr_.valid())
    {
        FatalErrorInFunction
            << "redirectFunctionObject not allocated" << exit(FatalError);
    }

    return redirectFunctionObjectPtr_();
}


Foam::wordList Foam::codedFunctionObject::fields() const
{
    return redirectFunctionObject().fields();
}


bool Foam::codedFunctionObject::execute()
{
    return redirectFunctionObject().execute();
}


bool Foam::codedFunctionObject::write()
{
    return redirectFunctionObject().write();
}


bool Foam::codedFunctionObject::end()
{
    return redirectFunctionObject().end();
}


bool Foam::codedFunctionObject::read(const dictionary& dict)
{
    if (functionObject::read(dict))
    {
        codedBase::read(dict);
        updateLibrary(dict);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
