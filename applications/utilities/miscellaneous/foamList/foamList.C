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
    \param -switches \n
        Print the DebugSwitches, InfoSwitches and OptimisationSwitches
    \param -registeredSwitches \n
        Print the registered DebugSwitches, InfoSwitches and
        OptimisationSwitches supporting run-time modification
    \param -unset \n
        print switches declared in libraries but not set in etc/controlDict

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "simpleObjectRegistry.H"
#include "IFstream.H"
#include "IOobject.H"
#include "HashSet.H"
#include "etcFiles.H"

using namespace Foam;

void listSwitches
(
    const wordList& debugSwitches,
    const wordList& infoSwitches,
    const wordList& optSwitches,
    const bool unset
)
{
    if (unset)
    {
        fileNameList controlDictFiles = findEtcFiles("controlDict", true);
        dictionary controlDict;
        forAllReverse(controlDictFiles, cdfi)
        {
            controlDict.merge(dictionary(IFstream(controlDictFiles[cdfi])()));
        }

        wordHashSet controlDictDebug
        (
            controlDict.subDict("DebugSwitches").sortedToc()
        );

        wordHashSet controlDictInfo
        (
            controlDict.subDict("InfoSwitches").sortedToc()
        );

        wordHashSet controlDictOpt
        (
            controlDict.subDict("OptimisationSwitches").sortedToc()
        );


        IOobject::writeDivider(Info);

        wordHashSet hashset;
        hashset = debugSwitches;
        hashset -= controlDictDebug;
        Info<< "Unset DebugSwitches: " << hashset.sortedToc() << endl;

        hashset = infoSwitches;
        hashset -= controlDictInfo;
        Info<< "Unset InfoSwitches: " << hashset.sortedToc() << endl;

        hashset = optSwitches;
        hashset -= controlDictOpt;
        Info<< "Unset OptimisationSwitches: " << hashset.sortedToc() << endl;
    }
    else
    {
        IOobject::writeDivider(Info);
        Info<< "DebugSwitches: " << debugSwitches << endl;
        Info<< "InfoSwitches: " << infoSwitches << endl;
        Info<< "OptimisationSwitches: " << optSwitches << endl;
    }
}


void listSwitches(const argList& args)
{
    if (args.optionFound("registeredSwitches"))
    {
        listSwitches
        (
            debug::debugObjects().sortedToc(),
            debug::infoObjects().sortedToc(),
            debug::optimisationObjects().sortedToc(),
            args.optionFound("unset")
        );
    }
    else
    {
        listSwitches
        (
            debug::debugSwitches().sortedToc(),
            debug::infoSwitches().sortedToc(),
            debug::optimisationSwitches().sortedToc(),
            args.optionFound("unset")
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "switches",
        "Switches declared in libraries but not set in etc/controlDict"
    );
    argList::addBoolOption
    (
        "registeredSwitches",
        "Switches registered for run-time modification"
    );
    argList::addBoolOption
    (
        "unset",
        "Switches declared in libraries but not set in etc/controlDict"
    );

    argList args(argc, argv);

    if (!args.options().size())
    {
        args.printUsage();
    }
    else if
    (
        args.optionFound("switches")
     || args.optionFound("registeredSwitches")
    )
    {
        listSwitches(args);
    }

    Info<< "done" << endl;

    return 0;
}


// ************************************************************************* //
