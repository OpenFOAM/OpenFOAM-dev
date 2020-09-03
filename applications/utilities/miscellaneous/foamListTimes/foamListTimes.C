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

Application
    foamListTimes

Description
    List times using timeSelector.

Usage
    \b foamListTimes [OPTION]

    Options:
      - \par -rm
        Remove selected time directories

      - \par -processor
        List times from processor0/ directory

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    entry::disableFunctionEntries = true;
    writeInfoHeader = false;

    argList::addNote("List times using timeSelector");
    timeSelector::addOptions(true, true);
    argList::noParallel();
    argList::addBoolOption
    (
        "processor",
        "list times from processor0/ directory"
    );
    argList::addBoolOption
    (
        "rm",
        "remove selected time directories"
    );
    argList::addBoolOption
    (
        "withFunctionObjects",
        "execute functionObjects"
    );
    argList::addBoolOption
    (
        "withFunctionEntries",
        "execute functionEntries"
    );

    #include "setRootCase.H"

    entry::disableFunctionEntries = !args.optionFound("withFunctionEntries");

    label nProcs = 0;

    // Create the processor databases
    PtrList<Time> databases(1);

    if (args.optionFound("processor"))
    {
        // Determine the processor count
        nProcs = fileHandler().nProcs(args.path());

        if (!nProcs)
        {
            FatalErrorInFunction
                << "No processor* directories found"
                << exit(FatalError);
        }

        // Create the processor databases
        databases.setSize(nProcs);

        forAll(databases, proci)
        {
            databases.set
            (
                proci,
                new Time
                (
                    Time::controlDictName,
                    args.rootPath(),
                    args.caseName()/fileName(word("processor") + name(proci))
                )
            );
        }
    }
    else
    {
        databases.set
        (
            0,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()
            )
        );
    }

    // Use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    if (args.optionFound("rm"))
    {
        if (args.optionFound("processor"))
        {
            for (label proci=0; proci<nProcs; proci++)
            {
                const fileName procPath
                (
                    args.path()/(word("processor") + name(proci))
                );

                forAll(timeDirs, timeI)
                {
                    const fileName procTimePath
                    (
                        fileHandler().filePath(procPath/timeDirs[timeI].name())
                    );

                    if (isDir(procTimePath))
                    {
                        rmDir(procTimePath);
                    }
                }
            }
        }
        else
        {
            forAll(timeDirs, timeI)
            {
                const fileName procTimePath
                (
                    fileHandler().filePath(args.path()/timeDirs[timeI].name())
                );

                rmDir(procTimePath);
            }
        }
    }
    else
    {
        forAll(timeDirs, timeI)
        {
            Info<< timeDirs[timeI].name() << endl;
        }
    }

    return 0;
}


// ************************************************************************* //
