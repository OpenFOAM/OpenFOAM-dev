/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Test-nonUniformTable1D

Description
    Tests the lookup of values of a linear function from and non-uniform
    table and reports any error.

\*---------------------------------------------------------------------------*/

#include "nonUniformTable1DThermophysicalFunction.H"
#include "IFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    dictionary dict(IFstream("thermoDict")());

    thermophysicalFunctions::nonUniformTable1D table(dict);

    const label n = 1000;
    const scalar T0 = table.values().first().first();
    const scalar Tn = table.values().last().first();
    const scalar deltaT = (Tn - T0)/n;

    for (int i = 0; i<n + 1; i++)
    {
        const scalar T = T0 + i*deltaT;
        const scalar f = table.f(0, T);

        if (mag(T - f) > small)
        {
            FatalError<< "failed" << exit(FatalError) << endl;
        }
    }

    Info<< "\nLookup of a linear function from a non-uniform table is correct\n"
        << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
