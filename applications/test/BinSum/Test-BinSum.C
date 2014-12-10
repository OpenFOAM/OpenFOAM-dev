/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Test-BinSum

Description
    Test BinSum container

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "BinSum.H"
#include "IOstreams.H"
#include "Random.H"
#include "scalarField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Random rndGen(0);

    scalarField samples(10000000);
    forAll(samples, i)
    {
        samples[i] = rndGen.scalar01();
    }

    const scalar min = 0;
    const scalar max = 1;
    const scalar delta = 0.1;

    BinSum<scalar, scalarField> count(min, max, delta);
    BinSum<scalar, scalarField> sum(min, max, delta);

    forAll(samples, i)
    {
        count.add(samples[i], 1);
        sum.add(samples[i], samples[i]);
    }

    Info<< "sum    : " << sum << endl;
    Info<< "count  : " << count << endl;
    Info<< "average: " << sum/count << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
