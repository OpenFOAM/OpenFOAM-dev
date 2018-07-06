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

Description
    Miscellaneous tests for HashTable

\*---------------------------------------------------------------------------*/

#include "HashTable.H"
#include "HashPtrTable.H"
#include "Map.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    HashTable<label, Foam::string> table1
    {
        {"kjhk", 10},
        {"kjhk2", 12}
    };

    Info<< "table1: " << table1 << nl
        << "toc: " << table1.toc() << endl;

    HashTable<label, label, Hash<label>> table2
    {
        {3, 10},
        {5, 12},
        {7, 16}
    };

    Info<< "table2: " << table2 << nl
        << "toc: " << table2.toc() << endl;

    Map<label> table3(1);
    table3.transfer(table2);

    Info<< "table2: " << table2 << nl
        << "toc: " << table2.toc() << endl;

    Info<< "table3: " << table3 << nl
        << "toc: " << table3.toc() << endl;

    Map<label> table4(table3.xfer());

    Info<< "table3: " << table3 << nl
        << "toc: " << table3.toc() << endl;

    Info<< "table4: " << table4 << nl
        << "toc: " << table4.toc() << endl;

    HashPtrTable<label, Foam::string> ptable1(0);
    ptable1.insert("kjhkjh", new label(10));

    Info<< "PtrTable toc: " << ptable1.toc() << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
