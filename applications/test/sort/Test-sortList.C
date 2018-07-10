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

\*---------------------------------------------------------------------------*/

#include "SortableList.H"
#include "ListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    labelList orig(8);
    orig[0] = 7;
    orig[1] = 9;
    orig[2] = 1;
    orig[3] = 2;
    orig[4] = 4;
    orig[5] = 7;
    orig[6] = 4;
    orig[7] = 0;

    labelList order;

    labelList a(orig);
    sortedOrder(a, order);

    SortableList<label> aReverse(a.size());
    aReverse = a;

    Info<< "unsorted: " << a << endl;
    sort(a);
    Info<< "sorted:   " << a << endl;
    Info<< "indices:  " << order << endl;

    aReverse.reverseSort();
    Info<< "reverse sorted:   " << aReverse << endl;
    Info<< "reverse indices:  " << aReverse.indices() << endl;

    SortableList<label> b(orig);
    Info<< "unsorted: " << orig << endl;
    Info<< "sorted:   " << b << endl;
    Info<< "indices:  " << b.indices() << endl;

    Info<< "shrunk:   " << b.shrink() << endl;
    Info<< "indices:  " << b.indices() << endl;

    // repeat by assignment
    b = orig;
    Info<< "unsorted: " << b << endl;
    b.sort();

    Info<< "sorted:   " << b << endl;
    Info<< "indices:  " << b.indices() << endl;

    // find unique/duplicate values
    b = orig;

    Info<< "unsorted: " << b << endl;
    uniqueOrder(b, order);
    Info<< "unique:  " << order << endl;
    duplicateOrder(b, order);
    Info<< "duplicate:" << order << endl;

    // sort reverse
    Info<< "unsorted: " << b << endl;
    b.reverseSort();
    Info<< "rsort:    " << b << endl;
    Info<< "indices:  " << b.indices() << endl;

    // transfer assignment
    a = orig;
    b.transfer(a);
    Info<< "unsorted: " << b << endl;
    b.sort();

    Info<< "sorted:   " << b << endl;
    Info<< "indices:  " << b.indices() << endl;

    a.transfer(b);

    Info<< "plain:    " << a << endl;
    Info<< "sorted:   " << b << endl;
    Info<< "indices:  " << b.indices() << endl;

    // sort/duplicate/unique with identical values
    b.setSize(8);
    b = 5;

    Info<< "unsorted: " << b << endl;

    uniqueOrder(b, order);
    Info<< "unique:  " << order << endl;
    duplicateOrder(b, order);
    Info<< "duplicate:" << order << endl;
    b.sort();

    Info<< "sorted:   " << b << endl;
    Info<< "indices:  " << b.indices() << endl;

    // with a single value
    b.setSize(1);

    Info<< "unsorted: " << b << endl;
    uniqueOrder(b, order);
    Info<< "unique:  " << order << endl;
    duplicateOrder(b, order);
    Info<< "duplicate:" << order << endl;
    b.sort();

    Info<< "sorted:   " << b << endl;
    Info<< "indices:  " << b.indices() << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
