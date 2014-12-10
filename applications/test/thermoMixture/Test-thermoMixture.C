/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    ThermoMixture

Description

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "IFstream.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "constTransport.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    typedef constTransport
    <
        species::thermo
        <
            hConstThermo<perfectGas<specie> >,
            sensibleEnthalpy
        >
    > ThermoType;

    dictionary dict(IFstream("thermoDict")());

    ThermoType t1(dict.subDict("specie1"));
    ThermoType t2(dict.subDict("specie2"));

    Info<< "Checking Cp of mixture of hConstThermo" << endl;

    Info<< "W 1, 2, (1 + 2) = " << t1.W() << " " << t2.W() << " "
        << (t1 + t2).W() << endl;

    Info<< "Cp 1, 2, 1 + 2 = " << t1.cp(1, 1) << " " << t2.cp(1, 1) << " "
        << (t1 + t2).cp(1, 1) << endl;

    ThermoType t3(t1);
    t3 += t2;
    Info<< "Cp (1 += 2) = " << t3.cp(1, 1) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
