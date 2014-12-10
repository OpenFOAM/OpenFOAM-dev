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

#include "codedFunctionObject.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "SHA1Digest.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::codedFunctionObject::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", redirectType_);
    dynCode.setFilterVariable("codeRead", codeRead_);
    dynCode.setFilterVariable("codeExecute", codeExecute_);
    dynCode.setFilterVariable("codeEnd", codeEnd_);
    dynCode.setFilterVariable("codeData", codeData_);
    dynCode.setFilterVariable("codeTimeSet", codeTimeSet_);
    //dynCode.setFilterVariable("codeWrite", codeWrite_);

    // compile filtered C template
    dynCode.addCompileFile("functionObjectTemplate.C");
    dynCode.addCompileFile("FilterFunctionObjectTemplate.C");

    // copy filtered H template
    dynCode.addCopyFile("FilterFunctionObjectTemplate.H");
    dynCode.addCopyFile("functionObjectTemplate.H");
    dynCode.addCopyFile("IOfunctionObjectTemplate.H");

    // debugging: make BC verbose
    //         dynCode.setFilterVariable("verbose", "true");
    //         Info<<"compile " << redirectType_ << " sha1: "
    //             << context.sha1() << endl;

    // define Make/options
    dynCode.setMakeOptions
        (
            "EXE_INC = -g \\\n"
            "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
            + context.options()
            + "\n\nLIB_LIBS = \\\n"
            + "    -lOpenFOAM \\\n"
            + "    -lfiniteVolume \\\n"
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedFunctionObject::codedFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict,
    bool readNow
)
:
    functionObject(name),
    codedBase(),
    time_(time),
    dict_(dict)
{
    if (readNow)
    {
        read(dict_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::codedFunctionObject::~codedFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::functionObject&
Foam::codedFunctionObject::redirectFunctionObject() const
{
    if (!redirectFunctionObjectPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", redirectType_);

        redirectFunctionObjectPtr_ = functionObject::New
        (
            redirectType_,
            time_,
            constructDict
        );
    }
    return redirectFunctionObjectPtr_();
}


bool Foam::codedFunctionObject::start()
{
    updateLibrary(redirectType_);
    return redirectFunctionObject().start();
}


bool Foam::codedFunctionObject::execute(const bool forceWrite)
{
    updateLibrary(redirectType_);
    return redirectFunctionObject().execute(forceWrite);
}


bool Foam::codedFunctionObject::end()
{
    updateLibrary(redirectType_);
    return redirectFunctionObject().end();
}


bool Foam::codedFunctionObject::timeSet()
{
    updateLibrary(redirectType_);
    return redirectFunctionObject().timeSet();
}


bool Foam::codedFunctionObject::read(const dictionary& dict)
{
    dict.lookup("redirectType") >> redirectType_;

    const entry* dataPtr = dict.lookupEntryPtr
    (
        "codeData",
        false,
        false
    );
    if (dataPtr)
    {
        codeData_ = stringOps::trim(dataPtr->stream());
        stringOps::inplaceExpand(codeData_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeData_,
            dataPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* readPtr = dict.lookupEntryPtr
    (
        "codeRead",
        false,
        false
    );
    if (readPtr)
    {
        codeRead_ = stringOps::trim(readPtr->stream());
        stringOps::inplaceExpand(codeRead_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeRead_,
            readPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* execPtr = dict.lookupEntryPtr
    (
        "codeExecute",
        false,
        false
    );
    if (execPtr)
    {
        codeExecute_ = stringOps::trim(execPtr->stream());
        stringOps::inplaceExpand(codeExecute_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeExecute_,
            execPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* endPtr = dict.lookupEntryPtr
    (
        "codeEnd",
        false,
        false
    );
    if (endPtr)
    {
        codeEnd_ = stringOps::trim(endPtr->stream());
        stringOps::inplaceExpand(codeEnd_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeEnd_,
            endPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* timeSetPtr = dict.lookupEntryPtr
    (
        "codeTimeSet",
        false,
        false
    );
    if (timeSetPtr)
    {
        codeTimeSet_ = stringOps::trim(timeSetPtr->stream());
        stringOps::inplaceExpand(codeTimeSet_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeTimeSet_,
            timeSetPtr->startLineNumber(),
            dict.name()
        );
    }

    updateLibrary(redirectType_);
    return redirectFunctionObject().read(dict);
}


void Foam::codedFunctionObject::updateMesh(const mapPolyMesh&)
{}


void Foam::codedFunctionObject::movePoints(const polyMesh&)
{}


// ************************************************************************* //
