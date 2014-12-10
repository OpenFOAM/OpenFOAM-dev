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

Application
    foamInfoExec

Description
    Interrogates a case and prints information to stdout.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dictionary.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "interrogates a case and prints information to stdout"
    );

    argList::noBanner();
    argList::addBoolOption("times", "list available times");
    argList::addBoolOption("latestTime", "list last time");
    argList::addBoolOption
    (
        "keywords",
        "report keywords for the specified dictionary"
    );
    #include "addDictOption.H"
    argList::addOption
    (
        "entry",
        "name",
        "report the named entry for the specified dictionary"
    );

    #include "setRootCase.H"

    if (args.optionFound("times"))
    {
        instantList times
        (
            Foam::Time::findTimes(args.rootPath()/args.caseName())
        );

        forAll(times, i)
        {
            Info<< times[i].name() << endl;
        }
    }
    else if (args.optionFound("latestTime"))
    {
        instantList times
        (
            Foam::Time::findTimes(args.rootPath()/args.caseName())
        );

        Info<< times.last().name() << endl;
    }

    if (args.optionFound("dict"))
    {
        fileName dictPath = args["dict"];
        const fileName dictFileName
        (
            dictPath.isAbsolute()
          ? dictPath
          : args.rootPath()/args.caseName()/args["dict"]
        );

        IFstream dictFile(dictFileName);

        if (dictFile.good())
        {
            dictionary dict(dictFile);

            if (args.optionFound("entry"))
            {
                fileName entryName(args.option("entry"));

                const entry* entPtr = NULL;

                if (entryName.find('.') != string::npos)
                {
                    // New syntax
                    entPtr = dict.lookupScopedEntryPtr
                    (
                        entryName,
                        false,
                        true            // wildcards
                    );
                }
                else
                {
                    // Old syntax
                    wordList entryNames(entryName.components(':'));
                    if (dict.found(entryNames[0]))
                    {
                        entPtr = &dict.lookupEntry
                        (
                            entryNames[0],
                            false,
                            true            // wildcards
                        );

                        for (int i=1; i<entryNames.size(); ++i)
                        {
                            if (entPtr->dict().found(entryNames[i]))
                            {
                                entPtr = &entPtr->dict().lookupEntry
                                (
                                    entryNames[i],
                                    false,
                                    true    // wildcards
                                );
                            }
                            else
                            {
                                FatalErrorIn(args.executable())
                                    << "Cannot find sub-entry " << entryNames[i]
                                    << " in entry " << args["entry"]
                                    << " in dictionary " << dictFileName;
                                FatalError.exit(3);
                            }
                        }
                    }
                }


                if (entPtr)
                {
                    if (args.optionFound("keywords"))
                    {
                        /*
                        if (ent[1] != token::BEGIN_BLOCK)
                        {
                            FatalErrorIn(args.executable())
                                << "Cannot find entry "
                                << args["entry"]
                                << " in dictionary " << dictFileName
                                << " is not a sub-dictionary";
                            FatalError.exit(4);
                        }
                        */

                        const dictionary& dict = entPtr->dict();
                        forAllConstIter(dictionary, dict, iter)
                        {
                            Info<< iter().keyword() << endl;
                        }
                    }
                    else
                    {
                        Info<< *entPtr << endl;
                    }
                }
                else
                {
                    FatalErrorIn(args.executable())
                        << "Cannot find entry "
                        << entryName
                        << " in dictionary " << dictFileName;
                    FatalError.exit(2);
                }
            }
            else if (args.optionFound("keywords"))
            {
                forAllConstIter(dictionary, dict, iter)
                {
                    Info<< iter().keyword() << endl;
                }
            }
            else
            {
                Info<< dict;
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Cannot open file " << dictFileName;
            FatalError.exit(1);
        }
    }

    return 0;
}


// ************************************************************************* //
