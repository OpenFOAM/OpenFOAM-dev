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
#include "OFstream.H"
#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    scalarField x(100);
    scalarField r(x.size());

    forAll(x, i)
    {
        x[i] = -3 + 0.06*i;
    }

    scalarField b(0.5*(1.0 + erf(x)));
    scalarField c(1.0 - b);
    scalarField gradb((1/::sqrt(constant::mathematical::pi))*exp(-sqr(x)));
    scalarField lapb(-2*x*gradb);

    r = lapb*b*c/(gradb*gradb);


    graph("r", "x", "r", x, r).write("r", "xmgr");

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
