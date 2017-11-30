/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    gen6DoF

Description
    Generate simple sinusoidal 6-DoF motion control-file.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "vector.H"
#include "Vector2D.H"
#include "Tuple2.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // End-time of the table
    const scalar endTime = 40;

    // Number of entries in the table
    const label nTimes = 100;

    // Amplitude of the translation [m]
    const vector transAmp(2, 3, 2);

    // Frequency of the translation [rad/s]
    const vector transOmega(0.5, 0.8, 0.4);

    // Amplitude of the rotation [deg]
    const vector rotAmp(30, 10, 10);

    // Frequency of the rotation [rad/s]
    const vector rotOmega(0.4, 0.7, 0.5);

    List<Tuple2<scalar,  Vector2D<vector>>> timeValues(nTimes);

    forAll(timeValues, i)
    {
        scalar t = (endTime*i)/(nTimes - 1);
        timeValues[i].first() = t;

        timeValues[i].second()[0] = vector
        (
            transAmp.x()*Foam::sin(transOmega.x()*t),
            transAmp.y()*Foam::sin(transOmega.y()*t),
            transAmp.z()*Foam::sin(transOmega.z()*t)
        );

        timeValues[i].second()[1] = vector
        (
            rotAmp.x()*Foam::sin(rotOmega.x()*t),
            rotAmp.y()*Foam::sin(rotOmega.y()*t),
            rotAmp.z()*Foam::sin(rotOmega.z()*t)
        );
    }

    {
        OFstream dataFile("6DoF.dat");
        dataFile << timeValues << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
