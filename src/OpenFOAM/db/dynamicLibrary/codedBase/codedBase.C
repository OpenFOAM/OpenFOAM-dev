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

#include "codedBase.H"
#include "dlLibraryTable.H"
#include "regIOobject.H"
#include "OSspecific.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(codedBase, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::codedBase::filterCodeName(const word& name)
{
    word result(name);

    if (!isalpha(result[0]))
    {
        FatalErrorInFunction
            << "Cannot construct code name from function name \"" << name
            << "\" as the first character is not alphabetic"
            << exit(FatalError);
    }

    for (word::size_type i = 1; i < name.size(); ++ i)
    {
        const bool valid = isalnum(result[i]) || result[i] == '_';

        if (!valid)
        {
            result[i] = '_';
        }
    }

    return result;
}


void Foam::codedBase::checkLibrary
(
    const dictionary& dict,
    const fileName& libPath,
    void* lib,
    const string& uniqueFuncName,
    const bool load
) const
{
    // Provision for manual execution of code before loading/unloading
    if (dlSymFound(lib, uniqueFuncName))
    {
        const loaderFunctionType function =
            reinterpret_cast<loaderFunctionType>
            (
                dlSym(lib, uniqueFuncName)
            );

        if (function)
        {
            (*function)(load);
            return;
        }
    }

    FatalIOErrorInFunction(dict)
        << "Failed looking up symbol " << uniqueFuncName << nl
        << "from " << libPath << exit(FatalIOError);
}


void* Foam::codedBase::loadLibrary
(
    const dictionary& dict,
    const fileName& libPath,
    const string& uniqueFuncName
) const
{
    void* lib = 0;

    if (!libPath.empty())
    {
        if (libs.open(libPath, false))
        {
            lib = libs.findLibrary(libPath);

            // Verify the loaded version
            if (lib)
            {
                checkLibrary(dict, libPath, lib, uniqueFuncName, true);
            }
        }
    }

    return lib;
}


void Foam::codedBase::unloadLibrary
(
    const dictionary& dict,
    const fileName& libPath,
    const string& uniqueFuncName
) const
{
    if (libPath.empty())
    {
        return;
    }

    void* lib = libs.findLibrary(libPath);

    if (!lib)
    {
        return;
    }

    checkLibrary(dict, libPath, lib, uniqueFuncName, false);

    if (!libs.close(libPath, false))
    {
        FatalIOErrorInFunction(dict)
            << "Failed unloading library " << libPath
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::verbatimString Foam::codedBase::expandCodeString
(
    const dictionary& dict,
    const word& codeKey,
    const word& codeDictVar
) const
{
    verbatimString codeString;

    if (dict.found(codeKey))
    {
        codeString = dict.lookupOrDefault<verbatimString>
        (
            codeKey,
            verbatimString::null
        );

        stringOps::inplaceExpandCodeString
        (
            codeString,
            dict,
            codeDictVar
        );
    }

    return codeString;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedBase::codedBase
(
    const word& name,
    const dictionary& dict,
    const wordList& codeKeys,
    const wordList& codeDictVars,
    const word& codeOptionsFileName,
    const wordList& compileFiles,
    const wordList& copyFiles
)
:
    dynamicCode
    (
        dict,
        filterCodeName(name),
        filterCodeName(name),
        codeKeys,
        codeDictVars,
        codeOptionsFileName,
        compileFiles,
        copyFiles
    )
{}


Foam::codedBase::codedBase
(
    const dictionary& dict,
    const wordList& codeKeys,
    const wordList& codeDictVars,
    const word& codeOptionsFileName,
    const wordList& compileFiles,
    const wordList& copyFiles
)
:
    codedBase
    (
        dict.lookup("name"),
        dict,
        codeKeys,
        codeDictVars,
        codeOptionsFileName,
        compileFiles,
        copyFiles
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::codedBase::~codedBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::codedBase::updateLibrary(const dictionary& dict) const
{
    const fileName libPath = this->libPath();

    // The correct library was already loaded => we are done
    if (libs.findLibrary(libPath))
    {
        return false;
    }

    Info<< "Using dynamicCode for " << type() << " " << codeName()
        << " at line " << dict.startLineNumber()
        << " in " << dict.name() << endl;

    // May need to unload old library
    unloadLibrary
    (
        dict,
        oldLibPath_,
        dynamicCode::libraryBaseName(oldLibPath_)
    );

    // Try loading an existing library (avoid compilation when possible)
    if (!loadLibrary(dict, libPath, codeSha1Name()))
    {
        createLibrary(dict);

        if (!loadLibrary(dict, libPath, codeSha1Name()))
        {
            FatalIOErrorInFunction(dict)
                << "Failed to load " << libPath << exit(FatalIOError);
        }
    }

    // Retain for future reference
    oldLibPath_ = libPath;

    return true;
}


// ************************************************************************* //
