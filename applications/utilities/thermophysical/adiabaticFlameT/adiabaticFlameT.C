/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    adiabaticFlameT

Description
    Calculates the adiabatic flame temperature for a given fuel over a
    range of unburnt temperatures and equivalence ratios.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "etcFiles.H"
#include "dimensionedTypes.H"

#include "specie.H"
#include "perfectGas.H"
#include "thermo.H"
#include "janafThermo.H"
#include "absoluteEnthalpy.H"

using namespace Foam;

typedef species::thermo<janafThermo<perfectGas<specie>>, absoluteEnthalpy>
    thermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::validArgs.append("properties dictionary");
    argList args(argc, argv);

    const fileName propertiesDictName = args[1];

    // Construct control dictionary
    IFstream propertiesDict(propertiesDictName);

    // Check propertiesDict stream is OK
    if (!propertiesDict.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << propertiesDictName
            << exit(FatalError);
    }

    dictionary control(propertiesDict);


    scalar P(control.lookup<scalar>("P"));
    scalar T0(control.lookup<scalar>("T0"));
    const word fuelName(control.lookup("fuel"));
    scalar n(control.lookup<scalar>("n"));
    scalar m(control.lookup<scalar>("m"));


    Info<< nl << "Reading thermodynamic data dictionary" << endl;

    fileName thermoDataFileName(findEtcFile("thermoData/thermoData"));

    // Construct control dictionary
    IFstream thermoDataFile(thermoDataFileName);

    // Check thermoData stream is OK
    if (!thermoDataFile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << thermoDataFileName
            << exit(FatalError);
    }

    dictionary thermoData(thermoDataFile);


    scalar stoicO2 = n + m/4.0;
    scalar stoicN2 = (0.79/0.21)*stoicO2;
    scalar stoicCO2 = n;
    scalar stoicH2O = m/2.0;

    thermo FUEL
    (
        "fuel",
        thermo(thermoData.subDict(fuelName))
    );
    Info<< "fuel " << FUEL << ';' << endl;
    FUEL *= FUEL.W();

    thermo O2
    (
        "O2",
        thermo(thermoData.subDict("O2"))
    );
    O2 *= O2.W();

    thermo N2
    (
        "N2",
        thermo(thermoData.subDict("N2"))
    );
    N2 *= N2.W();

    thermo CO2
    (
        "CO2",
        thermo(thermoData.subDict("CO2"))
    );
    CO2 *= CO2.W();

    thermo H2O
    (
        "H2O",
        thermo(thermoData.subDict("H2O"))
    );
    H2O *= H2O.W();

    thermo oxidant
    (
        "oxidant",
        stoicO2*O2
      + stoicN2*N2
    );
    Info<< "oxidant " << (1/oxidant.Y())*oxidant << ';' << endl;

    dimensionedScalar stoichiometricAirFuelMassRatio
    (
        "stoichiometricAirFuelMassRatio",
        dimless,
        oxidant.Y()/FUEL.W()
    );

    Info<< "stoichiometricAirFuelMassRatio "
        << stoichiometricAirFuelMassRatio << ';' << endl;

    for (int i=0; i<300; i++)
    {
        scalar equiv = (i + 1)*0.01;
        scalar ft = 1/(1 + stoichiometricAirFuelMassRatio.value()/equiv);

        Info<< "phi = " << equiv << nl
            << "ft = " << ft << endl;

        scalar o2 = (1.0/equiv)*stoicO2;
        scalar n2 = (0.79/0.21)*o2;
        scalar fres = max(1.0 - 1.0/equiv, 0.0);
        scalar ores = max(1.0/equiv - 1.0, 0.0);
        scalar fburnt = 1.0 - fres;

        thermo reactants
        (
            "reactants",
            FUEL + (1.0/equiv)*oxidant
        );
        Info<< "reactants " << (1/reactants.Y())*reactants << ';' << endl;

        thermo burntProducts
        (
            "burntProducts",
          + (n2 - (0.79/0.21)*ores*stoicO2)*N2
          + fburnt*stoicCO2*CO2
          + fburnt*stoicH2O*H2O
        );
        Info<< "burntProducts "
            << (1/burntProducts.Y())*burntProducts << ';' << endl;

        thermo products
        (
            "products",
            fres*FUEL
          + n2*N2
          + fburnt*stoicCO2*CO2
          + fburnt*stoicH2O*H2O
          + ores*stoicO2*O2
        );

        Info<< "products " << (1/products.Y())*products << ';' << endl;

        scalar Tad = products.THa(reactants.Ha(P, T0), P, 1000.0);
        Info<< "Tad = " << Tad << nl << endl;
    }

    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //
