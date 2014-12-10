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
    foamUpgradeFvSolution

Description
    Simple tool to upgrade the syntax of system/fvSolution.solvers.

Usage
    foamUpgradeFvSolution [OPTION]

    \param -test \n
    Suppress writing the updated fvSolution file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "solution.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "upgrade the syntax of system/fvSolution::solvers"
    );

    argList::noParallel();
    argList::addBoolOption
    (
        "test",
        "suppress writing the updated system/fvSolution file"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary solutionDict
    (
        IOobject
        (
            "fvSolution",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    label nChanged = 0;
    entry* e = solutionDict.lookupEntryPtr("solvers", false, false);
    if (e && e->isDict())
    {
        nChanged = solution::upgradeSolverDict(e->dict(), true);
    }

    Info<< nChanged << " solver settings changed" << nl << endl;
    if (nChanged)
    {
        if (args.optionFound("test"))
        {
            Info<< "-test option: no changes made" << nl << endl;
        }
        else
        {
            if (mvBak(solutionDict.objectPath(), "old"))
            {
                Info<< "Backup to    "
                    << (solutionDict.objectPath() + ".old") << nl;
            }

            solutionDict.writeObject
            (
                IOstream::ASCII,
                IOstream::currentVersion,
                IOstream::UNCOMPRESSED
            );

            Info<< "Write  to    "
                << solutionDict.objectPath() << nl << endl;
        }
    }

    return 0;
}


// ************************************************************************* //
