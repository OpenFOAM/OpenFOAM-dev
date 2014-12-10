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
    chemkinToFoam

Description
    Converts CHEMKINIII thermodynamics and reaction data files into
    OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "chemkinReader.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("CHEMKINFile");
    argList::validArgs.append("CHEMKINThermodynamicsFile");
    argList::validArgs.append("FOAMChemistryFile");
    argList::validArgs.append("FOAMThermodynamicsFile");

    argList::addBoolOption
    (
        "newFormat",
        "read Chemkin thermo file in new format"
    );

    argList args(argc, argv);

    bool newFormat = args.optionFound("newFormat");

    speciesTable species;

    chemkinReader cr(args[1], species, args[2], newFormat);

    OFstream reactionsFile(args[3]);
    reactionsFile
        << "species" << cr.species() << token::END_STATEMENT << nl << nl;

    cr.reactions().write(reactionsFile);


    OFstream thermoFile(args[4]);
    cr.speciesThermo().write(thermoFile);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
