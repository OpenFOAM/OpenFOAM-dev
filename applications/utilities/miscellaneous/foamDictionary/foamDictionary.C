/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

Application
    foamDictionary

Description
    Interrogates and manipulates dictionaries.

Usage
    \b foamDictionary [OPTION] dictionary

      - \par -entry \<name\>
        Selects an entry

      - \par -keywords \<name\>
        Prints the keywords (of the selected entry or of the top level if
        no entry was selected

      - \par -add \<value\>
        Adds the entry (should not exist yet)

      - \par -set \<value\>
        Adds or replaces the entry

      - \par -merge \<value\>
        Merges the entry

      - \par -dict
        Set, add or merge entry from a dictionary

      - \par -remove
        Remove the selected entry

      - \par -diff \<dictionary\>
        Write differences with respect to the specified dictionary
        (or sub entry if -entry specified)

      - \par -expand
        Read the specified dictionary file, expand the macros etc. and write
        the resulting dictionary to standard output.

      - \par -includes
        List the \c #include and \c #includeIfPresent files to standard output

      - \par -disableFunctionEntries
        Do not expand macros or directives (#include etc)

    Example usage:
      - Change simulation to run for one timestep only:
        \verbatim
          foamDictionary system/controlDict -entry stopAt -set writeNow
        \endverbatim

      - Change solver:
        \verbatim
           foamDictionary system/fvSolution -entry solvers.p.solver -set PCG
        \endverbatim

      - Print bc type:
        \verbatim
           foamDictionary 0/U -entry boundaryField.movingWall.type
        \endverbatim

      - Change bc parameter:
        \verbatim
           foamDictionary 0/U -entry boundaryField.movingWall.value \
             -set "uniform (2 0 0)"
        \endverbatim

      - Change whole bc type:
        \verbatim
          foamDictionary 0/U -entry boundaryField.movingWall \
            -set "{type uniformFixedValue; uniformValue (2 0 0);}"
        \endverbatim

      - Write the differences with respect to a template dictionary:
        \verbatim
          foamDictionary 0/U -diff $FOAM_ETC/templates/closedVolume/0/U
        \endverbatim

      - Write the differences in boundaryField with respect to a
        template dictionary:
        \verbatim
          foamDictionary 0/U -diff $FOAM_ETC/templates/closedVolume/0/U \
            -entry boundaryField
        \endverbatim

      - Change patch type:
        \verbatim
          foamDictionary constant/polyMesh/boundary \
            -entry entry0.fixedWalls.type -set patch
        \endverbatim
        This uses special parsing of Lists which stores these in the
        dictionary with keyword 'entryDDD' where DDD is the position
        in the dictionary (after ignoring the FoamFile entry).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "Pair.H"
#include "IFstream.H"
#include "OFstream.H"
#include "includeEntry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Read dictionary from file and return
//  Sets steam to binary mode if specified in the optional header
IOstream::streamFormat readDict(dictionary& dict, const fileName& dictFileName)
{
    IOstream::streamFormat dictFormat = IOstream::ASCII;

    IFstream dictFile(dictFileName);
    if (!dictFile().good())
    {
        FatalErrorInFunction
            << "Cannot open file " << dictFileName
            << exit(FatalError, 1);
    }

    // Read the first entry from the dictionary
    autoPtr<entry> firstEntry(entry::New(dictFile()));

    // If the first entry is the "FoamFile" header dictionary
    // read and set the stream format
    if (firstEntry->isDict() && firstEntry->keyword() == "FoamFile")
    {
        dictFormat = IOstream::formatEnum(firstEntry->dict().lookup("format"));
        dictFile().format(dictFormat);
    }

    // Add the first entry to the dictionary
    dict.add(firstEntry);

    // Read and add the rest of the dictionary entries
    // preserving the "FoamFile" header dictionary if present
    dict.read(dictFile(), true);

    return dictFormat;
}


//- Converts old scope syntax to new syntax
word scope(const fileName& entryName)
{
    if (entryName.find(':') != string::npos)
    {
        wordList entryNames(entryName.components(':'));

        word entry(entryNames[0]);
        for (label i = 1; i < entryNames.size(); i++)
        {
            entry += word('.') + entryNames[i];
        }
        return entry;
    }
    else
    {
        return entryName;
    }
}


//- Extracts dict name and keyword
Pair<word> dictAndKeyword(const word& scopedName)
{
    string::size_type i = scopedName.find_last_of(".");
    if (i != string::npos)
    {
        return Pair<word>
        (
            scopedName.substr(0, i),
            scopedName.substr(i+1, string::npos)
        );
    }
    else
    {
        return Pair<word>("", scopedName);
    }
}


const dictionary& lookupScopedDict
(
    const dictionary& dict,
    const word& subDictName
)
{
    if (subDictName == "")
    {
        return dict;
    }
    else
    {
        const entry* entPtr = dict.lookupScopedEntryPtr
        (
            subDictName,
            false,
            false
        );
        if (!entPtr || !entPtr->isDict())
        {
            FatalIOErrorInFunction(dict)
                << "keyword " << subDictName
                << " is undefined in dictionary "
                << dict.name() << " or is not a dictionary"
                << endl
                << "Valid keywords are " << dict.keys()
                << exit(FatalIOError);
        }
        return entPtr->dict();
    }
}


void remove(dictionary& dict, const dictionary& removeDict)
{
    forAllConstIter(dictionary, removeDict, iter)
    {
        const entry* entPtr = dict.lookupEntryPtr
        (
            iter().keyword(),
            false,
            false
        );

        if (entPtr)
        {
            if (entPtr->isDict())
            {
                if (iter().isDict())
                {
                    remove
                    (
                        const_cast<dictionary&>(entPtr->dict()),
                        iter().dict()
                    );

                    // Check if dictionary is empty
                    if (!entPtr->dict().size())
                    {
                        dict.remove(iter().keyword());
                    }
                }
            }
            else if (!iter().isDict())
            {
                if (*entPtr == iter())
                {
                    dict.remove(iter().keyword());
                }
            }
        }
    }
}


int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    writeInfoHeader = false;

    argList::addNote("manipulates dictionaries");
    argList::validArgs.append("dictionary file");
    argList::addBoolOption("keywords", "list keywords");
    argList::addOption("entry", "name", "report/select the named entry");
    argList::addBoolOption
    (
        "value",
        "Print entry value"
    );
    argList::addOption
    (
        "set",
        "value",
        "Set entry value or add new entry"
    );
    argList::addOption
    (
        "add",
        "value",
        "Add a new entry"
    );
    argList::addOption
    (
        "merge",
        "value",
        "Merge entry"
    );
    argList::addBoolOption
    (
        "dict",
        "Set, add or merge entry from a dictionary."
    );
    argList::addBoolOption
    (
        "remove",
        "Remove the entry."
    );
    argList::addOption
    (
        "diff",
        "dict",
        "Write differences with respect to the specified dictionary"
    );
    argList::addBoolOption
    (
        "includes",
        "List the #include/#includeIfPresent files to standard output"
    );
    argList::addBoolOption
    (
        "expand",
        "Read the specified dictionary file, expand the macros etc. and write "
        "the resulting dictionary to standard output"
    );
    argList::addBoolOption
    (
        "disableFunctionEntries",
        "Disable expansion of dictionary directives - #include, #codeStream etc"
    );

    argList args(argc, argv);

    const bool listIncludes = args.optionFound("includes");

    if (listIncludes)
    {
        Foam::functionEntries::includeEntry::log = true;
    }

    if (args.optionFound("disableFunctionEntries"))
    {
        entry::disableFunctionEntries = true;
    }


    const fileName dictFileName(args[1]);
    dictionary dict;
    IOstream::streamFormat dictFormat = readDict(dict, dictFileName);

    bool changed = false;

    if (listIncludes)
    {
        return 0;
    }
    else if (args.optionFound("expand"))
    {
        IOobject::writeBanner(Info)
            <<"//\n// " << dictFileName << "\n//\n";
        dict.write(Info, false);
        IOobject::writeDivider(Info);

        return 0;
    }


    // Second dictionary for -diff
    fileName diffFileName;
    dictionary diffDict;

    if (args.optionReadIfPresent("diff", diffFileName))
    {
        readDict(diffDict, diffFileName);
    }


    word entryName;
    if (args.optionReadIfPresent("entry", entryName))
    {
        word scopedName(scope(entryName));

        string newValue;
        if
        (
            args.optionReadIfPresent("set", newValue)
         || args.optionReadIfPresent("add", newValue)
         || args.optionReadIfPresent("merge", newValue)
        )
        {
            const bool overwrite = args.optionFound("set");
            const bool merge = args.optionFound("merge");

            Pair<word> dAk(dictAndKeyword(scopedName));
            const dictionary& d(lookupScopedDict(dict, dAk.first()));

            entry* ePtr = nullptr;

            if (args.optionFound("dict"))
            {
                const fileName fromDictFileName(newValue);
                dictionary fromDict;
                readDict(fromDict, fromDictFileName);

                const entry* fePtr
                (
                    fromDict.lookupScopedEntryPtr
                    (
                        scopedName,
                        false,
                        true            // Support wildcards
                    )
                );

                if (!fePtr)
                {
                    FatalErrorInFunction
                        << "Cannot find entry " << entryName
                        << " in file " << fromDictFileName
                        << exit(FatalError, 1);
                }

                ePtr = fePtr->clone().ptr();
            }
            else
            {
                IStringStream str(string(dAk.second()) + ' ' + newValue + ';');
                ePtr = entry::New(str).ptr();
            }

            if (overwrite)
            {
                const_cast<dictionary&>(d).set(ePtr);
            }
            else
            {
                const_cast<dictionary&>(d).add(ePtr, merge);
            }
            changed = true;

            // Print the changed entry
            // const entry* entPtr = dict.lookupScopedEntryPtr
            // (
            //     scopedName,
            //     false,
            //     true            // Support wildcards
            // );
            // if (entPtr)
            // {
            //     Info<< *entPtr;
            // }
        }
        else if (args.optionFound("remove"))
        {
            // Extract dictionary name and keyword
            Pair<word> dAk(dictAndKeyword(scopedName));

            const dictionary& d(lookupScopedDict(dict, dAk.first()));
            const_cast<dictionary&>(d).remove(dAk.second());
            changed = true;
        }
        else
        {
            // Optionally remove a second dictionary
            if (args.optionFound("diff"))
            {
                Pair<word> dAk(dictAndKeyword(scopedName));

                const dictionary& d(lookupScopedDict(dict, dAk.first()));
                const dictionary& d2(lookupScopedDict(diffDict, dAk.first()));

                const entry* ePtr =
                    d.lookupEntryPtr(dAk.second(), false, true);
                const entry* e2Ptr =
                    d2.lookupEntryPtr(dAk.second(), false, true);

                if (ePtr && e2Ptr)
                {
                    if (*ePtr == *e2Ptr)
                    {
                        const_cast<dictionary&>(d).remove(dAk.second());
                    }
                    else if (ePtr->isDict() && e2Ptr->isDict())
                    {
                        remove
                        (
                            const_cast<dictionary&>(ePtr->dict()),
                            e2Ptr->dict()
                        );
                    }
                }
            }


            const entry* entPtr = dict.lookupScopedEntryPtr
            (
                scopedName,
                false,
                true            // Support wildcards
            );

            if (entPtr)
            {
                if (args.optionFound("keywords"))
                {
                    const dictionary& dict = entPtr->dict();
                    forAllConstIter(dictionary, dict, iter)
                    {
                        Info<< iter().keyword() << endl;
                    }
                }
                else
                {
                    if (args.optionFound("value"))
                    {
                        if (entPtr->isStream())
                        {
                            const tokenList& tokens = entPtr->stream();
                            forAll(tokens, i)
                            {
                                Info<< tokens[i];
                                if (i < tokens.size() - 1)
                                {
                                    Info<< token::SPACE;
                                }
                            }
                            Info<< endl;
                        }
                        else if (entPtr->isDict())
                        {
                            Info<< entPtr->dict();
                        }
                    }
                    else
                    {
                        Info<< *entPtr;
                    }
                }
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Cannot find entry " << entryName
                    << exit(FatalIOError, 2);
            }
        }
    }
    else if (args.optionFound("keywords"))
    {
        forAllConstIter(dictionary, dict, iter)
        {
            Info<< iter().keyword() << endl;
        }
    }
    else if (args.optionFound("diff"))
    {
        remove(dict, diffDict);
        dict.write(Info, false);
    }
    else
    {
        dict.write(Info, false);
    }

    if (changed)
    {
        OFstream os(dictFileName, dictFormat);
        IOobject::writeBanner(os);
        dict.write(os, false);
        IOobject::writeEndDivider(os);
    }

    return 0;
}


// ************************************************************************* //
