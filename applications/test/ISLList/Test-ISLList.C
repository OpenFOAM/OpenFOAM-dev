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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "IOstreams.H"
#include "ISLList.H"

using namespace Foam;

class Scalar
:
    public ISLList<Scalar>::link
{
public:

    scalar data_;

    Scalar()
    :
        data_(0)
    {}

    Scalar(scalar s)
    :
        data_(s)
    {}

    friend Ostream& operator<<(Ostream& os, const Scalar& s)
    {
        os  << s.data_;
        return os;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    ISLList<Scalar> myList;

    for (int i = 0; i<10; i++)
    {
        myList.append(new Scalar(1.3*i));
    }

    myList.append(new Scalar(100.3));
    myList.append(new Scalar(500.3));

    Info<< nl << "And again using STL iterator: " << nl << endl;

    forAllIter(SLList<scalar>, myList, iter)
    {
        Info<< "element:" << *iter << endl;
    }

    Info<< nl << "And again using STL const_iterator: " << nl << endl;

    const ISLList<Scalar>& const_myList = myList;

    forAllConstIter(SLList<scalar>, const_myList, iter)
    {
        Info<< "element:" << *iter << endl;
    }


    Info<< nl << "Testing transfer: " << nl << endl;
    Info<< "original: " << myList << endl;

    ISLList<Scalar> newList;
    newList.transfer(myList);

    Info<< nl << "source: " << myList << nl
        << nl << "target: " << newList << endl;

    Info<< nl << "Bye." << endl;
    return 0;
}


// ************************************************************************* //
