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


void* Foam::codedBase::loadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& dict
) const
{
    void* lib = 0;

    // Avoid compilation by loading an existing library
    if (!libPath.empty())
    {
        if (libs.open(libPath, false))
        {
            lib = libs.findLibrary(libPath);

            // Verify the loaded version and unload if needed
            if (lib)
            {
                // Provision for manual execution of code after loading
                if (dlSymFound(lib, globalFuncName))
                {
                    const loaderFunctionType function =
                        reinterpret_cast<loaderFunctionType>
                        (
                            dlSym(lib, globalFuncName)
                        );

                    if (function)
                    {
                        (*function)(true);    // Force load
                    }
                    else
                    {
                        FatalIOErrorInFunction
                        (
                            dict
                        )   << "Failed looking up symbol " << globalFuncName
                            << nl << "from " << libPath << exit(FatalIOError);
                    }
                }
                else
                {
                    FatalIOErrorInFunction
                    (
                        dict
                    )   << "Failed looking up symbol " << globalFuncName << nl
                        << "from " << libPath << exit(FatalIOError);

                    lib = 0;
                    if (!libs.close(libPath, false))
                    {
                        FatalIOErrorInFunction
                        (
                            dict
                        )   << "Failed unloading library "
                            << libPath
                            << exit(FatalIOError);
                    }
                }
            }
        }
    }

    return lib;
}


void Foam::codedBase::unloadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& dict
) const
{
    void* lib = 0;

    if (libPath.empty())
    {
        return;
    }

    lib = libs.findLibrary(libPath);

    if (!lib)
    {
        return;
    }

    // Provision for manual execution of code before unloading
    if (dlSymFound(lib, globalFuncName))
    {
        const loaderFunctionType function =
            reinterpret_cast<loaderFunctionType>
            (
                dlSym(lib, globalFuncName)
            );

        if (function)
        {
            (*function)(false);    // Force unload
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Failed looking up symbol " << globalFuncName << nl
                << "from " << libPath << exit(FatalIOError);
        }
    }

    if (!libs.close(libPath, false))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Failed unloading library " << libPath
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::verbatimString Foam::codedBase::expandCodeString
(
    const word& codeKey,
    const word& codeDictVar,
    const dictionary& dict
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
        oldLibPath_,
        dynamicCode::libraryBaseName(oldLibPath_),
        dict
    );

    // Try loading an existing library (avoid compilation when possible)
    if (!loadLibrary(libPath, codeSha1Name(), dict))
    {
        createLibrary(dict);

        if (!loadLibrary(libPath, codeSha1Name(), dict))
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
