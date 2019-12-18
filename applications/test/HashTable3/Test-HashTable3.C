/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Description
    Test speeds for some HashTable operations

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "HashTable.H"
#include "HashPtrTable.H"
#include "Map.H"
#include "cpuTime.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    const label nLoops = 30;
    const label nBase  = 100000;
    const label nSize  = nLoops * nBase;

    cpuTime timer;

    // ie, a
    // Map<label> map(2 * nSize);
    // HashTable<label, label, Hash<label>> map(2 * nSize);
    HashTable<label, label, Hash<label>> map(2 * nSize);

    Info<< "Constructed map of size: " << nSize
        << " (size " << map.size() << " capacity " << map.capacity() << ") "
        << "  " << timer.cpuTimeIncrement() << " s\n\n";

    for (label i = 0; i < nSize; i++)
    {
        map.insert(i, i);
    }
    Info<< "Inserted " << nSize << " elements"
        << " (size " << map.size() << " capacity " << map.capacity() << ") "
        << timer.cpuTimeIncrement() << " s\n";

    label elemI = 0;
    for (label iLoop = 0; iLoop < nLoops; iLoop++)
    {
        for (label i = 0; i < nBase; i++)
        {
            map.erase(elemI++);
        }

        map.shrink();
        Info<< "loop " << iLoop << " - Erased " << nBase << " elements"
            << " (size " << map.size() << " capacity " << map.capacity() << ") "
            << timer.cpuTimeIncrement() << " s\n";
    }

    return 0;
}

// ************************************************************************* //
