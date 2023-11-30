/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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
    Test-thermoMixture

Description

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "IFstream.H"
#include "specie.H"
#include "perfectGas.H"
#include "hConstThermo.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    typedef hConstThermo<perfectGas<specie>> ThermoType;

    dictionary dict(IFstream("thermoDict")());

    ThermoType t1("specie1", dict.subDict("specie1"));
    ThermoType t2("specie2", dict.subDict("specie2"));

    Info<< "Checking Cp of mixture of hConstThermo:" << nl << endl;

    ThermoType t1Plus2(0.5*t1 + 0.5*t2);

    ThermoType t1PlusEq2(0.5*t1);
    t1PlusEq2 += 0.5*t2;

    Info<< "    W_1/W_2/W_{1+2}/W_{1+=2} = "
        << t1.W() << "/"
        << t2.W() << "/"
        << t1Plus2.W() << "/"
        << t1PlusEq2.W()
        << " [kg/kmol] " << nl << endl;

    Info<< "Cp_1/Cp_2/Cp_{1+2}/Cp_{1+=2} = "
        << t1.Cp(1, 1) << "/"
        << t2.Cp(1, 1) << "/"
        << t1Plus2.Cp(1, 1) << "/"
        << t1PlusEq2.Cp(1, 1)
        << " [J/kg/K]" << endl
        << "                             = "
        << t1.Cp(1, 1)*t1.W() << "/"
        << t2.Cp(1, 1)*t2.W() << "/"
        << t1Plus2.Cp(1, 1)*t1Plus2.W() << "/"
        << t1PlusEq2.Cp(1, 1)*t1PlusEq2.W()
        << " [J/kmol/K]" << nl << endl;

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
