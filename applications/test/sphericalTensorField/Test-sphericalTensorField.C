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
    sphericalTensorFieldTest

\*---------------------------------------------------------------------------*/

#include "tensorField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    scalarField f1(1, 1);
    sphericalTensorField sf1(1, 1);
    sphericalTensorField sf2(1, 2);
    tensorField tf1(1, tensor::one);

    Info<< (tf1 & sf2) << endl;

    Info<< f1*sf1 << " " << sf1*3 << endl;

    Info<< ((sf1 + sf2) & (sf1 + sf2)) << endl;

    return 0;
}


// ************************************************************************* //
