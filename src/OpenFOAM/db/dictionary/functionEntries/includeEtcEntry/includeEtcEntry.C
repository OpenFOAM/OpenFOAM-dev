/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "includeEtcEntry.H"
#include "etcFiles.H"
#include "stringOps.H"
#include "addToMemberFunctionSelectionTable.H"
#include "IOobject.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::includeEtcEntry::typeName
(
    Foam::functionEntries::includeEtcEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include includeEtcEntry
int Foam::functionEntries::includeEtcEntry::debug(0);

bool Foam::functionEntries::includeEtcEntry::log(false);


namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEtcEntry,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEtcEntry,
        execute,
        primitiveEntryIstream
    );
}
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEtcEntry::includeEtcFileName
(
    const fileName& f,
    const dictionary& dict
)
{
    fileName fName(f);

    // Substitute dictionary and environment variables.
    // Allow empty substitutions.
    stringOps::inplaceExpand(fName, dict, true, true);

    if (fName.empty() || fName.isAbsolute())
    {
        return fName;
    }
    else
    {
        // Search the etc directories for the file
        return findEtcFile(fName);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEtcEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const fileName rawFName(is);
    const fileName fName
    (
        includeEtcFileName(rawFName, parentDict)
    );

    // IFstream ifs(fName);
    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEtcEntry::log)
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
        )   << "Cannot open etc file "
            << (ifs.name().size() ? ifs.name() : rawFName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


bool Foam::functionEntries::includeEtcEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const fileName rawFName(is);
    const fileName fName
    (
        includeEtcFileName(rawFName, parentDict)
    );

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEtcEntry::log)
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
        )   << "Cannot open etc file "
            << (ifs.name().size() ? ifs.name() : rawFName)
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}


// ************************************************************************* //
