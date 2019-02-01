/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
#include "SHA1Digest.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * Protected Static Data Members * * * * * * * * * * * //

const Foam::wordList Foam::codedFunctionObject::codeKeys_ =
    {
        "codeData",
        "codeEnd",
        "codeExecute",
        "codeInclude",
        "codeRead",
        "codeWrite",
        "localCode"
    };


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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::codedFunctionObject::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    dynCode.setFilterVariable("typeName", name_);

    // Compile filtered C template
    dynCode.addCompileFile("functionObjectTemplate.C");

    // Copy filtered H template
    dynCode.addCopyFile("functionObjectTemplate.H");

    // Debugging: make BC verbose
    // dynCode.setFilterVariable("verbose", "true");
    // Info<<"compile " << name_ << " sha1: "
    //     << context.sha1() << endl;

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
      + "    -lOpenFOAM \\\n"
      + "    -lfiniteVolume \\\n"
      + "    -lmeshTools \\\n"
      + context.libs()
    );
}


Foam::dlLibraryTable& Foam::codedFunctionObject::libs() const
{
    return const_cast<Time&>(time_).libs();
}


Foam::string Foam::codedFunctionObject::description() const
{
    return "functionObject " + name();
}


void Foam::codedFunctionObject::clearRedirect() const
{
    redirectFunctionObjectPtr_.clear();
}


const Foam::dictionary& Foam::codedFunctionObject::codeDict() const
{
    return dict_;
}


const Foam::wordList& Foam::codedFunctionObject::codeKeys() const
{
    return codeKeys_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedFunctionObject::codedFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject(name),
    codedBase(),
    time_(time),
    dict_(dict)
{
    read(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::codedFunctionObject::~codedFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::functionObject& Foam::codedFunctionObject::redirectFunctionObject() const
{
    if (!redirectFunctionObjectPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", name_);

        redirectFunctionObjectPtr_ = functionObject::New
        (
            name_,
            time_,
            constructDict
        );
    }
    return redirectFunctionObjectPtr_();
}


bool Foam::codedFunctionObject::execute()
{
    updateLibrary(name_);
    return redirectFunctionObject().execute();
}


bool Foam::codedFunctionObject::write()
{
    updateLibrary(name_);
    return redirectFunctionObject().write();
}


bool Foam::codedFunctionObject::end()
{
    updateLibrary(name_);
    return redirectFunctionObject().end();
}


bool Foam::codedFunctionObject::read(const dictionary& dict)
{
    // The name keyword is "name". "redirectType" is also maintained here
    // for backwards compatibility, but "name" is taken in preference and
    // is printed in the error message if neither keyword is present.
    name_ = word::null;
    name_ = dict.lookupOrDefault("redirectType", name_);
    name_ = dict.lookupOrDefault("name", name_);
    if (name_ == word::null)
    {
        dict.lookup("name"); // <-- generate error message with "name" in it
    }

    updateLibrary(name_);
    return redirectFunctionObject().read(dict);
}


// ************************************************************************* //
