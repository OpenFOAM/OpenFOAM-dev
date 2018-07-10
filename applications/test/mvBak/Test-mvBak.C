/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    writeInfoHeader = false;

    argList::noParallel();
    argList::validArgs.insert("file .. fileN");

    argList::removeOption("case");
    argList::addOption("ext", "bak");

    argList args(argc, argv, false, true);

    if (args.size() <= 1)
    {
        args.printUsage();
    }

    label ok = 0;

    for (label argI=1; argI < args.size(); ++argI)
    {
        const string& srcFile = args[argI];

        if (args.optionFound("ext"))
        {
            if (mvBak(srcFile, args["ext"]))
            {
                ok++;
            }
        }
        else
        {
            if (mvBak(srcFile))
            {
                ok++;
            }
        }
    }

    Info<< "mvBak called for " << args.size()-1
        << " files (moved " << ok << ")\n" << endl;

    return 0;
}


// ************************************************************************* //
