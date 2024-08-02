/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "dynamicCode.H"
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

Foam::word Foam::codedBase::codeName(const word& name)
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
    const dictionary& contextDict
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
                            contextDict
                        )   << "Failed looking up symbol " << globalFuncName
                            << nl << "from " << libPath << exit(FatalIOError);
                    }
                }
                else
                {
                    FatalIOErrorInFunction
                    (
                        contextDict
                    )   << "Failed looking up symbol " << globalFuncName << nl
                        << "from " << libPath << exit(FatalIOError);

                    lib = 0;
                    if (!libs.close(libPath, false))
                    {
                        FatalIOErrorInFunction
                        (
                            contextDict
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
    const dictionary& contextDict
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
                contextDict
            )   << "Failed looking up symbol " << globalFuncName << nl
                << "from " << libPath << exit(FatalIOError);
        }
    }

    if (!libs.close(libPath, false))
    {
        FatalIOErrorInFunction
        (
            contextDict
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


void Foam::codedBase::createLibrary
(
    const dictionary& dict,
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    const bool create =
        Pstream::master()
     || (regIOobject::fileModificationSkew <= 0);   // Not NFS

    if (create)
    {
        // Write files for new library
        if (!dynCode.upToDate(context))
        {
            // Filter with this context
            dynCode.reset(context);

            this->prepare(dynCode, context);

            if (!dynCode.copyOrCreateFiles(true))
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Failed writing files for" << nl
                    << dynCode.libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        if (!dynCode.wmakeLibso())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Failed wmake " << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }
    }


    // All processes must wait for compile to finish
    if (regIOobject::fileModificationSkew > 0)
    {
        //- Since the library has only been compiled on the master the
        //  other nodes need to pick this library up through NFS
        //  We do this by just polling a few times using the
        //  fileModificationSkew.

        const fileName libPath = dynCode.libPath();

        off_t mySize = fileSize(libPath);
        off_t masterSize = mySize;
        Pstream::scatter(masterSize);

        if (debug)
        {
            Pout<< endl<< "on processor " << Pstream::myProcNo()
                << " have masterSize:" << masterSize
                << " and localSize:" << mySize
                << endl;
        }

        if (mySize < masterSize)
        {
            if (debug)
            {
                Pout<< "Local file " << libPath
                    << " not of same size (" << mySize
                    << ") as master ("
                    << masterSize << "). Waiting for "
                    << regIOobject::fileModificationSkew
                    << " seconds." << endl;
            }
            sleep(regIOobject::fileModificationSkew);

            // Recheck local size
            mySize = Foam::fileSize(libPath);

            if (mySize < masterSize)
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Cannot read (NFS mounted) library " << nl
                    << libPath << nl
                    << "on processor " << Pstream::myProcNo()
                    << " detected size " << mySize
                    << " whereas master size is " << masterSize
                    << " bytes." << nl
                    << "If your case is not NFS mounted"
                    << " (so distributed) set fileModificationSkew"
                    << " to 0"
                    << exit(FatalIOError);
            }
        }

        if (debug)
        {
            Pout<< endl<< "on processor " << Pstream::myProcNo()
                << " after waiting: have masterSize:" << masterSize
                << " and localSize:" << mySize
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedBase::codedBase
(
    const word& name,
    const dictionary& dict,
    const wordList& codeKeys,
    const wordList& codeDictVars
)
:
    codeName_(codeName(name)),
    codeContext_(dict, codeKeys, codeDictVars)
{}


Foam::codedBase::codedBase
(
    const dictionary& dict,
    const wordList& codeKeys,
    const wordList& codeDictVars
)
:
    codedBase(codeName(dict.lookup("name")), dict, codeKeys, codeDictVars)
{}


Foam::codedBase::codedBase(const codedBase& cb)
:
    codeName_(cb.codeName_),
    codeContext_(cb.codeContext_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::codedBase::~codedBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::codedBase::codeName() const
{
    return codeName_;
}


Foam::string Foam::codedBase::description() const
{
    return this->type() + " " + codeName();
}


Foam::word Foam::codedBase::codeTemplateC(const word& baseTypeName) const
{
    return baseTypeName + "Template.C";
}


Foam::word Foam::codedBase::codeTemplateH(const word& baseTypeName) const
{
    return baseTypeName + "Template.H";
}


bool Foam::codedBase::updateLibrary(const dictionary& dict) const
{
    const word& name = codeName();

    dynamicCode::checkSecurity
    (
        "codedBase::updateLibrary()",
        dict
    );

    // codeName: name + _<sha1>
    // codeDir : name
    dynamicCode dynCode
    (
        name + codeContext_.sha1().str(true),
        name
    );
    const fileName libPath = dynCode.libPath();


    // The correct library was already loaded => we are done
    if (libs.findLibrary(libPath))
    {
        return false;
    }

    Info<< "Using dynamicCode for " << this->description().c_str()
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
    if (!loadLibrary(libPath, dynCode.codeName(), dict))
    {
        createLibrary(dict, dynCode, codeContext_);

        if (!loadLibrary(libPath, dynCode.codeName(), dict))
        {
            FatalIOErrorInFunction(dict)
                << "Failed to load " << libPath << exit(FatalIOError);
        }
    }

    // Retain for future reference
    oldLibPath_ = libPath;

    return true;
}


void Foam::codedBase::read(const dictionary& dict)
{
    codeContext_.read(dict);
}


void Foam::codedBase::write(Ostream& os) const
{
    if (codeName().size())
    {
        writeEntry(os, "name", codeName());
    }

    codeContext_.write(os);
}


// ************************************************************************* //
