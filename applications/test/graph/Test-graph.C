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
    graphTest

Description
    Test program for making graphs

\*---------------------------------------------------------------------------*/

#include "graph.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    scalarField alpha(200);
    scalarField phie(alpha.size());
    scalarField phic(alpha.size());

    forAll(alpha, i)
    {
        alpha[i] = scalar(i)/50.0;
    }

    scalar R = 5.0;

    phie = (R - 1)/(sqrt(R/(alpha + 1.0e-6)) + 1.0) + 1.0;
    phic = (R - 1)/(sqrt(R*alpha) + 1.0) + 1.0;


    graph phi("@f! (R = 5)", "@a!", "@f!", alpha);
    phi.insert
    (
        "@f!&e!", new curve("@f!&e!", curve::curveStyle::CONTINUOUS, phie)
    );
    phi.insert
    (
        "@f!&c!", new curve("@f!&c!", curve::curveStyle::CONTINUOUS, phic)
    );
    phi.write("phi", "xmgr");

    Info<< "end" << endl;
}


// ************************************************************************* //
