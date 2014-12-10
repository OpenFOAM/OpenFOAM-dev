/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    test

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOmanip.H"
#include "scalar.H"
#include "List.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(void)
{
    string st("sfdsf  sdfs23df sdf32f .  sdfsdff23/2sf32");
    Info<< word(string::validate<word>(st)) << "END" << endl;

    string st1("1234567");

    Info<< label(st1.size()) << tab << string(word(st1)) << endl;

    Info<< setw(20) << setprecision(3) << 1.234234 << endl;

    Info<< hex << 255 << endl;

    Info.operator Foam::OSstream&() << "stop" << endl;
}

// ************************************************************************* //
