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

Description

\*---------------------------------------------------------------------------*/

#include "string.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    string a("a"), b("b"), c("c"), d("d"), e("e"),
           f("f"), g("g"), h("h"), i("i"), j("j");

    Info<< "1) a = b + c + d + e + f + g;\n";
    a = b + c + d + e + f + g;
    Info<< a << endl;

    {
        Info<< "2) a = e + f + g;\n";
        a = e + f + g;
        Info<< a << endl;
    }

    Info<< "3) a = h + i + j;\n";
    a = h + i + j;
    Info<< a << endl;

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
