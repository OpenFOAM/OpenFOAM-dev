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
    equilibriumFlameT

Description
    Calculates the equilibrium flame temperature for a given fuel and
    pressure for a range of unburnt gas temperatures and equivalence
    ratios; the effects of dissociation on O2, H2O and CO2 are included.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "IFstream.H"
#include "etcFiles.H"
#include "IOmanip.H"
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
            << abort(FatalError);
    }

    dictionary control(propertiesDict);


    scalar P(control.lookup<scalar>("P"));
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
            << abort(FatalError);
    }

    dictionary thermoData(thermoDataFile);


    Info<< nl << "Reading thermodynamic data for relevant species"
        << nl << endl;

    // Reactants (mole-based)
    thermo FUEL(thermoData.subDict(fuelName)); FUEL *= FUEL.W();

    // Oxidant (mole-based)
    thermo O2(thermoData.subDict("O2")); O2 *= O2.W();
    thermo N2(thermoData.subDict("N2")); N2 *= N2.W();

    // Intermediates (mole-based)
    thermo H2(thermoData.subDict("H2")); H2 *= H2.W();

    // Products (mole-based)
    thermo CO2(thermoData.subDict("CO2")); CO2 *= CO2.W();
    thermo H2O(thermoData.subDict("H2O")); H2O *= H2O.W();
    thermo CO(thermoData.subDict("CO")); CO *= CO.W();


    // Product dissociation reactions

    thermo CO2BreakUp
    (
        CO2 == CO + 0.5*O2
    );

    thermo H2OBreakUp
    (
        H2O == H2 + 0.5*O2
    );


    // Stoichiometric number of moles of species for one mole of fuel
    scalar stoicO2 = n + m/4.0;
    scalar stoicN2 = (0.79/0.21)*(n + m/4.0);
    scalar stoicCO2 = n;
    scalar stoicH2O = m/2.0;

    // Oxidant
    thermo oxidant
    (
        "oxidant",
        stoicO2*O2
      + stoicN2*N2
    );

    dimensionedScalar stoichiometricAirFuelMassRatio
    (
        "stoichiometricAirFuelMassRatio",
        dimless,
        oxidant.Y()/FUEL.W()
    );

    Info<< "stoichiometricAirFuelMassRatio "
        << stoichiometricAirFuelMassRatio << ';' << endl;

    Info<< "Equilibrium flame temperature data ("
        << P/1e5 << " bar)" << nl << nl
        << setw(3) << "Phi"
        << setw(12) << "ft"
        << setw(7) << "T0"
        << setw(12) << "Tad"
        << setw(12) << "Teq"
        << setw(12) << "Terror"
        << setw(20) << "O2res (mole frac)" << nl
        << endl;


    // Loop over equivalence ratios
    for (int i=0; i<16; i++)
    {
        scalar equiv = 0.6 + i*0.05;
        scalar ft = 1/(1 + stoichiometricAirFuelMassRatio.value()/equiv);

    // Loop over initial temperatures
    for (int j=0; j<28; j++)
    {
        scalar T0 = 300.0 + j*100.0;

        // Number of moles of species for one mole of fuel
        scalar o2 = (1.0/equiv)*stoicO2;
        scalar n2 = (0.79/0.21)*o2;
        scalar fres = max(1.0 - 1.0/equiv, 0.0);
        scalar fburnt = 1.0 - fres;

        // Initial guess for number of moles of product species
        // ignoring product dissociation
        scalar oresInit = max(1.0/equiv - 1.0, 0.0)*stoicO2;
        scalar co2Init = fburnt*stoicCO2;
        scalar h2oInit = fburnt*stoicH2O;

        scalar ores = oresInit;
        scalar co2 = co2Init;
        scalar h2o = h2oInit;

        scalar co = 0.0;
        scalar h2 = 0.0;

        // Total number of moles in system
        scalar N = fres + n2 + co2 + h2o + ores;


        // Initial guess for adiabatic flame temperature
        scalar adiabaticFlameTemperature =
            T0
          + (fburnt/(1.0 + o2 + n2))/(1.0/(1.0 + (1.0 + 0.79/0.21)*stoicO2))
           *2000.0;

        scalar equilibriumFlameTemperature = adiabaticFlameTemperature;


        // Iteration loop for adiabatic flame temperature
        for (int j=0; j<20; j++)
        {
            if (j > 0)
            {
                co = co2*
                    min
                    (
                        CO2BreakUp.Kn(P, equilibriumFlameTemperature, N)
                       /::sqrt(max(ores, 0.001)),
                        1.0
                    );

                h2 = h2o*
                    min
                    (
                        H2OBreakUp.Kn(P, equilibriumFlameTemperature, N)
                       /::sqrt(max(ores, 0.001)),
                        1.0
                    );

                co2 = co2Init - co;
                h2o = h2oInit - h2;
                ores = oresInit + 0.5*co + 0.5*h2;
            }

            thermo reactants
            (
                FUEL + o2*O2 + n2*N2
            );

            thermo products
            (
                fres*FUEL + ores*O2 + n2*N2
              + co2*CO2 + h2o*H2O + co*CO + h2*H2
            );


            scalar equilibriumFlameTemperatureNew =
                products.THa(reactants.Ha(P, T0), P, adiabaticFlameTemperature);

            if (j==0)
            {
                adiabaticFlameTemperature = equilibriumFlameTemperatureNew;
            }
            else
            {
                equilibriumFlameTemperature = 0.5*
                (
                    equilibriumFlameTemperature
                  + equilibriumFlameTemperatureNew
                );
            }
        }

        Info<< setw(3) << equiv
            << setw(12) << ft
            << setw(7) << T0
            << setw(12) << adiabaticFlameTemperature
            << setw(12) << equilibriumFlameTemperature
            << setw(12)
            << adiabaticFlameTemperature - equilibriumFlameTemperature
            << setw(12) << ores/N
            << endl;
    }
    }

    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //
