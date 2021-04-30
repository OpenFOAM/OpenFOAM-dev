/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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

    Supports parallel operation for decomposed dictionary files associated with
    a case.  These may be mesh or field files or any other decomposed
    dictionaries.

Usage
    \b foamDictionary [OPTION] dictionary
      - \par -parallel
        Specify case as a parallel job

      - \par -doc
        Display the documentation in browser

      - \par -srcDoc
        Display the source documentation in browser

      - \par -help
        Print the usage

      - \par -entry \<name\>
        Selects an entry

      - \par -keywords \<name\>
        Prints the keywords (of the selected entry or of the top level if
        no entry was selected

      - \par -add \<value\>
        Adds the entry (should not exist yet)

      - \par -set \<value\>
        Adds or replaces the entry selected by \c -entry

      - \par -set \<substitutions\>
        Applies the list of substitutions

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

    Example usage:
      - Change simulation to run for one timestep only:
        \verbatim
          foamDictionary system/controlDict -entry stopAt -set writeNow
        \endverbatim

      - Change solver:
        \verbatim
           foamDictionary system/fvSolution -entry solvers/p/solver -set PCG
        \endverbatim

      - Print bc type:
        \verbatim
           foamDictionary 0/U -entry boundaryField/movingWall/type
        \endverbatim

      - Change bc parameter:
        \verbatim
           foamDictionary 0/U -entry boundaryField/movingWall/value \
             -set "uniform (2 0 0)"
        \endverbatim

      - Change bc parameter in parallel:
        \verbatim
           mpirun -np 4 foamDictionary 0.5/U \
             -entry boundaryField/movingWall/value \
             -set "uniform (2 0 0)" -parallel
        \endverbatim

      - Change whole bc type:
        \verbatim
          foamDictionary 0/U -entry boundaryField/movingWall \
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
            -entry entry0/fixedWalls/type -set patch
        \endverbatim
        This uses special parsing of Lists which stores these in the
        dictionary with keyword 'entryDDD' where DDD is the position
        in the dictionary (after ignoring the FoamFile entry).

      - Substitute multiple entries:
        \verbatim
          foamDictionary system/controlDict -set "startTime=2000, endTime=3000"
        \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "localIOdictionary.H"
#include "Pair.H"
#include "IFstream.H"
#include "OFstream.H"
#include "includeEntry.H"
#include "inputSyntaxEntry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Read dictionary from file and return
//  Sets steam to binary mode if specified in the optional header
IOstream::streamFormat readDict(dictionary& dict, const fileName& dictFileName)
{
    IOstream::streamFormat dictFormat = IOstream::ASCII;

    // Read the first entry and if it is FoamFile set the file format
    {
        IFstream dictFile(dictFileName);
        if (!dictFile().good())
        {
            FatalErrorInFunction
                << "Cannot open file " << dictFileName
                << exit(FatalError, 1);
        }

        // Check if the first token in the file is "FoamFile"
        // to avoid problems if the first entry is a variable or function
        token firstToken;
        dictFile.read(firstToken);
        if (firstToken.isWord() && firstToken.wordToken() == IOobject::foamFile)
        {
            dictFile.putBack(firstToken);

            // Read the first entry from the dictionary
            autoPtr<entry> firstEntry(entry::New(dictFile()));

            // If the first entry is the "FoamFile" header
            // read and set the stream format
            if
            (
                firstEntry->isDict()
             && firstEntry->keyword() == IOobject::foamFile
            )
            {
                dictFormat = IOstream::formatEnum
                (
                    firstEntry->dict().lookup("format")
                );
            }
        }
    }

    IFstream dictFile(dictFileName, dictFormat);

    // Read and add the rest of the dictionary entries
    // preserving the IOobject::foamFile header dictionary if present
    dict.read(dictFile(), true);

    return dictFormat;
}


//- Convert keyword syntax to "dot" if the dictionary is "dot" syntax
word dotToSlash(const fileName& entryName)
{
    if
    (
        functionEntries::inputSyntaxEntry::dot()
     && entryName.find('/') != string::npos
    )
    {
        wordList entryNames(entryName.components('/'));

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


void remove(dictionary& dict, const dictionary& removeDict)
{
    forAllConstIter(dictionary, removeDict, iter)
    {
        entry* entPtr = dict.lookupEntryPtr
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
                    remove(entPtr->dict(), iter().dict());

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


void substitute(dictionary& dict, string substitutions)
{
    wordReList args;
    List<Tuple2<word, string>> namedArgs;
    dictArgList(substitutions, args, namedArgs);

    forAll(namedArgs, i)
    {
        const Pair<word> dAk(dictAndKeyword(dotToSlash(namedArgs[i].first())));
        dictionary& subDict(dict.scopedDict(dAk.first()));
        IStringStream entryStream
        (
            dAk.second() + ' ' + namedArgs[i].second() + ';'
        );
        subDict.set(entry::New(entryStream).ptr());
    }
}


int main(int argc, char *argv[])
{
    argList::removeOption("case");
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
        "Set entry value, add new entry or apply list of substitutions"
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
    argList::addOption
    (
        "writePrecision",
        "label",
        "Write with the specified precision"
    );

    argList args(argc, argv);

    const bool listIncludes = args.optionFound("includes");

    if (listIncludes)
    {
        functionEntries::includeEntry::log = true;
    }

    // Do not expand functionEntries except during dictionary expansion
    // with the -expand option
    if (!args.optionFound("expand"))
    {
        entry::disableFunctionEntries = true;
    }

    // Set write precision
    if (args.optionFound("writePrecision"))
    {
        const label writePrecision = args.optionRead<label>("writePrecision");
        IOstream::defaultPrecision(writePrecision);
        Sout.precision(writePrecision);
        Pout.precision(IOstream::defaultPrecision());
    }

    const fileName dictPath(args[1]);

    Time* runTimePtr = nullptr;
    localIOdictionary* localDictPtr = nullptr;

    dictionary* dictPtr = nullptr;
    IOstream::streamFormat dictFormat = IOstream::ASCII;

    // When running in parallel read the dictionary as a case localIOdictionary
    // supporting file handlers
    if (Pstream::parRun())
    {
        if (!args.checkRootCase())
        {
            FatalError.exit();
        }

        runTimePtr = new Time(Time::controlDictName, args);

        const wordList dictPathComponents(dictPath.components());

        if (dictPathComponents.size() == 1)
        {
            FatalErrorInFunction
                << "File name " << dictPath
                << " does not contain an instance path needed in parallel"
                << exit(FatalError, 1);
        }

        const word instance = dictPathComponents[0];
        const fileName dictFileName
        (
            SubList<word>(dictPathComponents, dictPathComponents.size() - 1, 1)
        );

        scalar time;
        if (readScalar(instance.c_str(), time))
        {
            runTimePtr->setTime(time, 0);
        }

        localDictPtr = new localIOdictionary
        (
            IOobject
            (
                dictFileName,
                instance,
                *runTimePtr,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
    }
    else
    {
        dictPtr = new dictionary(dictPath);
        dictFormat = readDict(*dictPtr, dictPath);
    }

    dictionary& dict = localDictPtr ? *localDictPtr : *dictPtr;

    bool changed = false;

    if (listIncludes)
    {
        return 0;
    }
    else if (args.optionFound("expand") && !args.optionFound("entry"))
    {
        IOobject::writeBanner(Info)
            <<"//\n// " << dictPath << "\n//\n";

        // Change the format to ASCII
        if (dict.found(IOobject::foamFile))
        {
            dict.subDict(IOobject::foamFile).add
            (
                "format",
                IOstream::ASCII,
                true
            );
        }

        dict.dictionary::write(Info, false);
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
        const word scopedName(dotToSlash(entryName));

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

            const Pair<word> dAk(dictAndKeyword(scopedName));
            dictionary& subDict(dict.scopedDict(dAk.first()));

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
                Info << "New entry " << *ePtr << endl;
                subDict.set(ePtr);
            }
            else
            {
                subDict.add(ePtr, merge);
            }
            changed = true;
        }
        else if (args.optionFound("remove"))
        {
            // Extract dictionary name and keyword
            const Pair<word> dAk(dictAndKeyword(scopedName));

            dictionary& subDict(dict.scopedDict(dAk.first()));
            subDict.remove(dAk.second());
            changed = true;
        }
        else
        {
            // Optionally remove a second dictionary
            if (args.optionFound("diff"))
            {
                const Pair<word> dAk(dictAndKeyword(scopedName));

                dictionary& subDict(dict.scopedDict(dAk.first()));
                const dictionary& subDict2(diffDict.scopedDict(dAk.first()));

                entry* ePtr =
                    subDict.lookupEntryPtr(dAk.second(), false, true);
                const entry* e2Ptr =
                    subDict2.lookupEntryPtr(dAk.second(), false, true);

                if (ePtr && e2Ptr)
                {
                    if (*ePtr == *e2Ptr)
                    {
                        subDict.remove(dAk.second());
                    }
                    else if (ePtr->isDict() && e2Ptr->isDict())
                    {
                        remove(ePtr->dict(), e2Ptr->dict());
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
    else if (args.optionFound("set"))
    {
        const string substitutions(args.optionRead<string>("set"));
        substitute(dict, substitutions);
        changed = true;
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
        dict.dictionary::write(Info, false);
    }
    else
    {
        dict.dictionary::write(Info, false);
    }

    if (changed)
    {
        if (localDictPtr)
        {
            localDictPtr->regIOobject::write();
        }
        else if (dictPtr)
        {
            OFstream os(dictPath, dictFormat);
            IOobject::writeBanner(os);
            if (dictPtr->found(IOobject::foamFile))
            {
                os << IOobject::foamFile;
                dictPtr->subDict(IOobject::foamFile).write(os);
                dictPtr->remove(IOobject::foamFile);
                IOobject::writeDivider(os) << nl;
            }
            dictPtr->write(os, false);
            IOobject::writeEndDivider(os);
        }
    }

    delete dictPtr;
    delete localDictPtr;
    delete runTimePtr;

    return 0;
}


// ************************************************************************* //
