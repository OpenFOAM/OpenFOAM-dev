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

Application
    chemkinToFoam

Description
    Converts CHEMKINIII thermodynamics and reaction data files into
    OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "chemkinReader.H"
#include "OFstream.H"
#include "OStringStream.H"
#include "IStringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    // Increase the precision of the output for JANAF coefficients
    Ostream::defaultPrecision(10);

    argList::validArgs.append("CHEMKIN file");
    argList::validArgs.append("CHEMKIN thermodynamics file");
    argList::validArgs.append("CHEMKIN transport file");
    argList::validArgs.append("OpenFOAM chemistry file");
    argList::validArgs.append("OpenFOAM thermodynamics file");

    argList::addBoolOption
    (
        "newFormat",
        "read Chemkin thermo file in new format"
    );

    argList args(argc, argv);

    bool newFormat = args.optionFound("newFormat");

    speciesTable species;

    chemkinReader cr(species, args[1], args[3], args[2], newFormat);

    OFstream reactionsFile(args[4]);
    reactionsFile.writeKeyword("elements")
        << cr.elementNames() << token::END_STATEMENT << nl << nl;
    reactionsFile.writeKeyword("species")
        << cr.species() << token::END_STATEMENT << nl << nl;
    cr.reactions().write(reactionsFile);

    // Temporary hack to splice the specie composition data into the thermo file
    // pending complete integration into the thermodynamics structure

    OStringStream os;
    cr.speciesThermo().write(os);
    dictionary thermoDict(IStringStream(os.str())());

    wordList speciesList(thermoDict.toc());

    // Add elements
    forAll(speciesList, si)
    {
        dictionary elementsDict("elements");
        forAll(cr.specieComposition()[speciesList[si]], ei)
        {
            elementsDict.add
            (
                cr.specieComposition()[speciesList[si]][ei].name(),
                cr.specieComposition()[speciesList[si]][ei].nAtoms()
            );
        }

        thermoDict.subDict(speciesList[si]).add("elements", elementsDict);
    }

    thermoDict.write(OFstream(args[5])(), false);

    reactionsFile << nl;

    reactionsFile.writeKeyword("Tlow")
        << Reaction<gasHThermoPhysics>::TlowDefault
        << token::END_STATEMENT << nl;
    reactionsFile.writeKeyword("Thigh")
        << Reaction<gasHThermoPhysics>::ThighDefault
        << token::END_STATEMENT << nl << nl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
