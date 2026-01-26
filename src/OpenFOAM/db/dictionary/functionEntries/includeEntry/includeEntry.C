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

#include "includeEntry.H"
#include "stringOps.H"
#include "IOobject.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::log(false);

namespace Foam
{
namespace functionEntries
{
    defineFunctionTypeNameAndDebug(includeEntry, 0);
    addToRunTimeSelectionTable(functionEntry, includeEntry, dictionary);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::List<Foam::Tuple3<Foam::word, Foam::string, Foam::label>>
Foam::functionEntries::includeEntry::insertNamedArgs
(
    dictionary& contextDict,
    const Tuple2<string, label>& fNameArgs
) const
{
    List<Tuple3<word, string, label>> namedArgs;

    // If the next character is a '(' process the arguments
    if (fNameArgs.first().size())
    {
        // Parse the argument string
        word funcType;
        List<Tuple2<wordRe, label>> args;
        dictArgList(fNameArgs, funcType, args, namedArgs);

        // Add the named arguments as entries into the contextDict
        // temporarily renaming any existing entries with the same name
        forAll(namedArgs, i)
        {
            const Pair<word> dAk(dictAndKeyword(namedArgs[i].first()));
            dictionary& subDict(contextDict.scopedDict(dAk.first()));

            // Rename the original entry adding a '_'
            if (subDict.found(dAk.second()))
            {
                keyType tmpName(dAk.second());
                tmpName += '_';
                subDict.changeKeyword(dAk.second(), tmpName);
            }

            // Add the temporary argument entry
            IStringStream entryStream
            (
                dAk.second()
              + ' '
              + expandArg
                (
                    namedArgs[i].second(),
                    contextDict,
                    namedArgs[i].third()
                )
              + ';'
            );
            subDict.set(entry::New(entryStream).ptr());
        }
    }

    return namedArgs;
}


void Foam::functionEntries::includeEntry::removeInsertNamedArgs
(
    dictionary& contextDict,
    const List<Tuple3<word, string, label>>& namedArgs
) const
{
    forAll(namedArgs, i)
    {
        // Remove the temporary argument entry
        contextDict.remove(namedArgs[i].first());

        // Reinstate the original entry
        const Pair<word> dAk(dictAndKeyword(namedArgs[i].first()));
        dictionary& subDict(contextDict.scopedDict(dAk.first()));
        keyType tmpName(dAk.second());
        tmpName += '_';
        subDict.changeKeyword(tmpName, dAk.second());
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEntry::includeFileName
(
    const Istream& is,
    const fileName& f,
    const dictionary& dict
)
{
    fileName fName(f);

    // Substitute dictionary and environment variables. Allow empty
    // substitutions.
    stringOps::inplaceExpandEntry(fName, dict, true, true);

    if (fName.empty() || fName.isAbsolute())
    {
        return fName;
    }
    else
    {
        // relative name
        return fileName(is.name()).path()/fName;
    }
}


Foam::fileName Foam::functionEntries::includeEntry::includeFileName
(
    const fileName& dir,
    const fileName& f,
    const dictionary& dict
)
{
    fileName fName(f);

    // Substitute dictionary and environment variables. Allow empty
    // substitutions.
    stringOps::inplaceExpandEntry(fName, dict, true, true);

    if (fName.empty() || fName.isAbsolute())
    {
        return fName;
    }
    else
    {
        // relative name
        return dir/fName;
    }
}


bool Foam::functionEntries::includeEntry::execute
(
    dictionary& contextDict,
    Istream& is,
    const bool ifPresent
)
{
    const fileName fName
    (
        includeFileName(is.name().path(), this->fName(), contextDict)
    );

    if (!fileHandler().exists(fName))
    {
        if (ifPresent)
        {
            return true;
        }
        else
        {
            FatalIOErrorInFunction
            (
                is
            )   << "Cannot find include file " << fName
                << " while reading dictionary " << contextDict.name()
                << exit(FatalIOError);
        }
    }

    // Cache the optional named arguments
    // temporarily inserted into contextDict
    List<Tuple3<word, string, label>> namedArgs
    (
        insertNamedArgs(contextDict, args())
    );

    autoPtr<ISstream> ifsPtr
    (
        fileHandler().NewIFstream(fName, is.format(), is.version())
    );
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }

        // Cache the FoamFile entry if present
        dictionary foamFileDict;
        if (contextDict.found(IOobject::foamFile))
        {
            foamFileDict = contextDict.subDict(IOobject::foamFile);
        }

        // Read and clear the FoamFile entry
        contextDict.read(ifs);

        // Reinstate original FoamFile entry
        if (foamFileDict.size() != 0)
        {
            dictionary contextDictTmp(contextDict);
            contextDict.clear();
            contextDict.add(IOobject::foamFile, foamFileDict);
            contextDict += contextDictTmp;
        }
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open include file "
            << (ifs.name().size() ? ifs.name() : this->fName())
            << " while reading dictionary " << contextDict.name()
            << exit(FatalIOError);
    }

    // Remove named argument entries from contextDict
    // renaming any existing entries which had the same name
    removeInsertNamedArgs(contextDict, namedArgs);

    return true;
}


bool Foam::functionEntries::includeEntry::execute
(
    const dictionary& contextDict,
    primitiveEntry& contextEntry,
    Istream& is,
    const bool ifPresent
)
{
    const includeEntry ie(is.lineNumber(), contextDict, is);

    const fileName fName
    (
        includeFileName(is.name().path(), ie.fName(), contextDict)
    );

    if (!fileHandler().exists(fName))
    {
        if (ifPresent)
        {
            return true;
        }
        else
        {
            FatalIOErrorInFunction
            (
                is
            )   << "Cannot find include file " << fName
                << " while reading dictionary " << contextDict.name()
                << exit(FatalIOError);
        }
    }

    // Cache the optional named arguments
    // temporarily inserted into contextDict
    List<Tuple3<word, string, label>> namedArgs
    (
        ie.insertNamedArgs(const_cast<dictionary&>(contextDict), ie.args())
    );

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }

        contextEntry.read(contextDict, ifs);
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open include file "
            << (ifs.name().size() ? ifs.name() : ie.fName())
            << " while reading dictionary " << contextDict.name()
            << exit(FatalIOError);
    }

    // Remove named argument entries from contextDict
    // renaming any existing entries which had the same name
    ie.removeInsertNamedArgs(const_cast<dictionary&>(contextDict), namedArgs);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionEntries::includeEntry::includeEntry
(
    const functionName& functionType,
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    functionEntry
    (
        functionType,
        lineNumber,
        parentDict,
        is,
        readFileNameArgList(functionType, is)
    )
{}


Foam::functionEntries::includeEntry::includeEntry
(
    const label lineNumber,
    const dictionary& parentDict,
    Istream& is
)
:
    includeEntry(typeName, lineNumber, parentDict, is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::execute
(
    dictionary& contextDict,
    Istream& is
)
{
    return execute(contextDict, is, false);
}


bool Foam::functionEntries::includeEntry::execute
(
    const dictionary& contextDict,
    primitiveEntry& contextEntry,
    Istream& is
)
{
    return execute(contextDict, contextEntry, is, false);
}


// ************************************************************************* //
