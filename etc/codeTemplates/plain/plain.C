/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    NAME

Description
    Example plain (non-CFD) application, which can be compiled and tested with
    the following commands, which should print "42" followed by "true".

        NAME 42 guess && echo true
        NAME -multiply 4 3 number

\*---------------------------------------------------------------------------*/

#include "argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Minimal default options and do not print header when running
    #include "removeCaseOptions.H"
    writeInfoHeader = false;

    // Character width of the options listed by "-help"
    argList::usageMin = 24;

    // Example "-multiply <factor>" command line option
    argList::addOption
    (
        "multiply",
        "factor",
        "multiply the scalar by the factor"
    );

    // Example mandatory arguments
    argList::validArgs.append("scalar");
    argList::validArgs.append("name");

    // Example description for "-help"
    argList::addNote
    (
        "Prints '<name> = <scalar>*<factor>', where <factor> is optional\n"
        "Returns true if the word == 'guess'"
    );

    argList args(argc, argv);

    const scalar m = args.optionLookupOrDefault<scalar>("multiply", 1.0);
    const scalar s = args.argRead<scalar>(1);
    const word w(args.argRead<word>(2));

    Info<< w << " = " << s*m << endl;

    return w == "guess" ? 0 : 1 ;
}

// ************************************************************************* //
