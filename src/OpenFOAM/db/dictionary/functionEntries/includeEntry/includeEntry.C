/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "dictionary.H"
#include "IFstream.H"
#include "addToMemberFunctionSelectionTable.H"
#include "stringOps.H"
#include "IOobject.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::includeEntry::typeName
(
    Foam::functionEntries::includeEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include includeEntry
int Foam::functionEntries::includeEntry::debug(0);

bool Foam::functionEntries::includeEntry::log(false);


namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        primitiveEntryIstream
    );
}
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEntry::includeFileName
(
    Istream& is,
    const dictionary& dict
)
{
    fileName fName(is);
    // Substitute dictionary and environment variables. Allow empty
    // substitutions.
    stringOps::inplaceExpand(fName, dict, true, true);

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
    stringOps::inplaceExpand(fName, dict, true, true);

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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const fileName rawFName(is);
    const fileName fName
    (
        includeFileName(is.name().path(), rawFName, parentDict)
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
        if (parentDict.found(IOobject::foamFile))
        {
            foamFileDict = parentDict.subDict(IOobject::foamFile);
        }

        // Read and clear the FoamFile entry
        parentDict.read(ifs);

        // Reinstate original FoamFile entry
        if (foamFileDict.size() != 0)
        {
            dictionary parentDictTmp(parentDict);
            parentDict.clear();
            parentDict.add(IOobject::foamFile, foamFileDict);
            parentDict += parentDictTmp;
        }

        return true;
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open include file "
            << (ifs.name().size() ? ifs.name() : rawFName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


bool Foam::functionEntries::includeEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const fileName rawFName(is);
    const fileName fName
    (
        includeFileName(is.name().path(), rawFName, parentDict)
    );

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }
        entry.read(parentDict, ifs);
        return true;
    }
    else
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open include file "
            << (ifs.name().size() ? ifs.name() : rawFName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


// ************************************************************************* //
