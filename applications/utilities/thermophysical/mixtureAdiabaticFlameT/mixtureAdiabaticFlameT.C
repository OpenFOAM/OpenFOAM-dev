/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    mixtureAdiabaticFlameT

Description
    Calculates the adiabatic flame temperature for a given mixture
    at a given temperature.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "etcFiles.H"

#include "specie.H"
#include "perfectGas.H"
#include "thermo.H"
#include "janafThermo.H"
#include "absoluteEnthalpy.H"
#include "mixture.H"

using namespace Foam;

typedef species::thermo<janafThermo<perfectGas<specie>>, absoluteEnthalpy>
    thermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("properties dictionary");
    argList args(argc, argv);

    const fileName propertiesFileName(args[1]);

    // Construct properties dictionary
    IFstream propertiesFile(propertiesFileName);

    // Check propertiesFile stream is OK
    if (!propertiesFile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << propertiesFileName
            << abort(FatalError);
    }

    dictionary properties(propertiesFile);


    scalar P(readScalar(properties.lookup("P")));
    scalar T0(readScalar(properties.lookup("T0")));
    mixture rMix(properties.lookup("reactants"));
    mixture pMix(properties.lookup("products"));


    Info<< nl << "Reading thermodynamic data dictionary" << endl;

    fileName thermoDataFileName(findEtcFile("thermoData/thermoData"));

    // Construct properties dictionary
    IFstream thermoDataFile(thermoDataFileName);

    // Check thermoData stream is OK
    if (!thermoDataFile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << thermoDataFileName
            << abort(FatalError);
    }

    dictionary thermoData(thermoDataFile);


    thermo reactants
    (
        rMix[0].volFrac()*thermo(thermoData.subDict(rMix[0].name()))
    );

    for (label i = 1; i < rMix.size(); i++)
    {
        reactants = reactants
            + rMix[i].volFrac()*thermo(thermoData.subDict(rMix[i].name()));
    }


    thermo products
    (
        2*pMix[0].volFrac()*thermo(thermoData.subDict(pMix[0].name()))
    );

    for (label i = 1; i < pMix.size(); i++)
    {
        products = products
            + 2*pMix[i].volFrac()*thermo(thermoData.subDict(pMix[i].name()));
    }

    Info<< "Adiabatic flame temperature of mixture " << rMix.name() << " = "
         << products.THa(reactants.Ha(P, T0), P, 1000.0) << " K" << endl;

    return 0;
}


// ************************************************************************* //
