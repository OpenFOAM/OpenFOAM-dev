/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    foamList

Description
    Print the table of contents of selectable switches, classes etc. in the
    OpenFOAM libraries

    \par Command-line options
    \param -debug \n
        Print the DebugSwitches, InfoSwitches and OptimisationSwitches
        \param -unset \n
            print switches declared in libraries but not set in etc/controlDict
        \param -redundant \n
            print switches not declared in libraries but set in etc/controlDict

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "simpleObjectRegistry.H"
#include "IFstream.H"
#include "IOobject.H"
#include "HashSet.H"
#include "etcFiles.H"

using namespace Foam;

void listDebug(const argList& args)
{
    // Switches declared in libraries
    wordList libDebug(debug::debugObjects().sortedToc());
    wordList libInfo(debug::infoObjects().sortedToc());
    wordList libOpt(debug::optimisationObjects().sortedToc());

    if (args.optionFound("redundant") || args.optionFound("unset"))
    {
        fileNameList controlDictFiles = findEtcFiles("controlDict", true);
        dictionary controlDict;
        forAllReverse(controlDictFiles, cdfi)
        {
            controlDict.merge(dictionary(IFstream(controlDictFiles[cdfi])()));
        }

        wordHashSet controlDictDebug
        (
            controlDict.subDict("DebugSwitches").toc()
        );

        wordHashSet controlDictInfo
        (
            controlDict.subDict("InfoSwitches").toc()
        );

        wordHashSet controlDictOpt
        (
            controlDict.subDict("OptimisationSwitches").toc()
        );


        wordHashSet hashset;
        wordList listing;


        // List redundant switches
        if (args.optionFound("redundant"))
        {
            IOobject::writeDivider(Info);

            hashset = controlDictDebug;
            hashset -= libDebug;
            listing = hashset.toc();
            sort(listing);
            Info<< "Redundant DebugSwitches: " << listing << endl;

            hashset = controlDictInfo;
            hashset -= libInfo;
            listing = hashset.toc();
            sort(listing);
            Info<< "Redundant InfoSwitches: " << listing << endl;

            hashset = controlDictOpt;
            hashset -= libOpt;
            listing = hashset.toc();
            sort(listing);
            Info<< "Redundant OptimisationSwitches: " << listing << endl;
        }

        // List unset switches
        if (args.optionFound("unset"))
        {
            IOobject::writeDivider(Info);

            hashset = libDebug;
            hashset -= controlDictDebug;

            listing = hashset.toc();
            sort(listing);
            Info<< "Unset DebugSwitches: " << listing << endl;

            hashset = libInfo;
            hashset -= controlDictInfo;
            listing = hashset.toc();
            sort(listing);
            Info<< "Unset InfoSwitches: " << listing << endl;

            hashset = libOpt;
            hashset -= controlDictOpt;
            listing = hashset.toc();
            sort(listing);
            Info<< "Unset OptimisationSwitches: " << listing << endl;
        }
    }
    else
    {
        IOobject::writeDivider(Info);

        sort(libDebug);
        Info<< "DebugSwitches: " << libDebug << endl;

        sort(libInfo);
        Info<< "InfoSwitches: " << libInfo << endl;

        sort(libOpt);
        Info<< "OptimisationSwitches: " << libOpt << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "debug",
        "switches declared in libraries but not set in etc/controlDict"
    );
    argList::addBoolOption
    (
        "unset",
        "switches declared in libraries but not set in etc/controlDict"
    );
    argList::addBoolOption
    (
        "redundant",
        "switches not declared in libraries but set in etc/controlDict"
    );

    argList args(argc, argv);

    if (!args.options().size())
    {
        args.printUsage();
    }
    else if (args.optionFound("debug"))
    {
        listDebug(args);
    }

    Info<< "done" << endl;

    return 0;
}


// ************************************************************************* //
