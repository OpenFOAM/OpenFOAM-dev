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

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "boolList.H"
#include "HashSet.H"
#include "StaticHashTable.H"
#include "cpuTime.H"
#include <vector>
#include "PackedBoolList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
    const label n = 100000000;
    const label nReport = 1000000;

    cpuTime timer;

    // test inserts
    // PackedBoolList
    PackedBoolList packed;
    for (label i = 0; i < n; i++)
    {
        if ((i % nReport) == 0 && i)
        {
            Info<< "i:" << i << " in " << timer.cpuTimeIncrement() << " s"
                <<endl;
        }
        packed[i] = 1;
    }
    Info<< "insert test: " << n << " elements in "
        << timer.cpuTimeIncrement() << " s\n\n";

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
