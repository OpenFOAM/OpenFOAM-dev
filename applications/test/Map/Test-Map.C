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
    testMapIterators

Description
    For each time calculate the magnitude of velocity.

\*---------------------------------------------------------------------------*/

#include "Map.H"
#include <map>
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Map<bool> banana;

    banana.insert(5, true);

    // Taking a const iterator from find does not work!
    // Also, fails later on op==
    Map<bool>::const_iterator bananaIter = banana.find(5);

    // This works but now I can change the value.
    // Map<bool>::iterator bananaIter = banana.find(5);

    if (bananaIter == banana.end())
    {
        Info<< "not found" << endl;
    }
    else
    {
        Info<< "5 is " << bananaIter() << endl;
    }

    // Same with STL
    Info<< "Same with STL" << endl;

    std::map<label, bool> STLbanana;
    STLbanana[5] = true;
    std::map<label, bool>::const_iterator STLbananaIter = STLbanana.find(5);

    if (STLbananaIter == STLbanana.end())
    {
        Info<< "not found" << endl;
    }
    else
    {
        Info<< "5 is " << STLbananaIter->second << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
