/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "OSspecific.H"

#include "IOstreams.H"
#include "SLList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    SLList<scalar> myList{2.1, 3.4};
    myList = {2.1, 3.4, 4.3};

    for (int i = 0; i<10; i++)
    {
        myList.append(1.3*i);
    }

    myList.append(100.3);
    myList.append(500.3);

    Info<< nl << "And again using STL iterator: " << nl << endl;

    forAllIter(SLList<scalar>, myList, iter)
    {
        Info<< "element:" << *iter << endl;
    }

    Info<< nl << "And again using STL const_iterator: " << nl << endl;

    const SLList<scalar>& const_myList = myList;

    forAllConstIter(SLList<scalar>, const_myList, iter)
    {
        Info<< "element:" << *iter << endl;
    }

    forAllIter(SLList<scalar>, myList, iter)
    {
        Info<< "Removing element:" << *iter << endl;
        myList.remove(iter);
    }

    forAllConstIter(SLList<scalar>, const_myList, iter)
    {
        Info<< "element:" << *iter << endl;
    }


    for (int i = 0; i<10; i++)
    {
        myList.append(1.3*i);
    }

    myList.append(100.3);
    myList.append(500.3);

    Info<< nl << "Testing transfer: " << nl << endl;
    Info<< "original: " << myList << endl;

    SLList<scalar> newList;
    newList.transfer(myList);

    Info<< nl << "source: " << myList << nl
        << nl << "target: " << newList << endl;


    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
