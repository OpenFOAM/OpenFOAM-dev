/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

    argList::addOption
    (
        "precision",
        "label",
        "set the write precision"
    );

    argList args(argc, argv);

    label precision = 10;
    args.optionReadIfPresent("precision", precision);

    // Increase the precision of the output for JANAF coefficients
    Ostream::defaultPrecision(precision);

    bool newFormat = args.optionFound("newFormat");

    speciesTable species;
    chemkinReader cr(species, args[1], args[3], args[2], newFormat);
    const HashPtrTable<chemkinReader::thermoPhysics>& speciesThermo =
        cr.speciesThermo();

    dictionary thermoDict;
    thermoDict.add("species", cr.species());

    // Add the species thermo formatted entries
    {
        OStringStream os;
        speciesThermo.write(os);
        dictionary speciesThermoDict(IStringStream(os.str())());
        thermoDict.merge(speciesThermoDict);
    }

    // Temporary hack to splice the specie composition data into the thermo file
    // pending complete integration into the thermodynamics structure

    // Add elements
    forAllConstIter
    (
        HashPtrTable<chemkinReader::thermoPhysics>,
        speciesThermo,
        iter
    )
    {
        const word specieName(iter.key());

        dictionary elementsDict("elements");
        forAll(cr.specieComposition()[specieName], ei)
        {
            elementsDict.add
            (
                cr.specieComposition()[specieName][ei].name(),
                cr.specieComposition()[specieName][ei].nAtoms()
            );
        }

        thermoDict.subDict(specieName).add("elements", elementsDict);
    }

    thermoDict.write(OFstream(args[5])(), false);
    cr.reactions().write(OFstream(args[4])());

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
