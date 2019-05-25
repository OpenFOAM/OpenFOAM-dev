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

\*---------------------------------------------------------------------------*/

#include "IndirectList.H"
#include "IOstreams.H"

using namespace Foam;

template<class ListType>
void printInfo(const ListType& lst)
{
    Info<< "addr: " << lst.addressing() << nl
        << "list: " << lst << nl
        << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<double> completeList(10);

    forAll(completeList, i)
    {
        completeList[i] = 0.1*i;
    }

    Info<< "raw : " << completeList << nl << endl;


    List<label> addresses(5);
    addresses[0] = 1;
    addresses[1] = 0;
    addresses[2] = 7;
    addresses[3] = 8;
    addresses[4] = 5;

    IndirectList<double> idl1(completeList, addresses);

    printInfo(idl1);

    addresses[4] = 1;
    addresses[3] = 0;
    addresses[2] = 7;
    addresses[1] = 8;
    addresses[0] = 5;

    idl1.resetAddressing(move(addresses));

    printInfo(idl1);

    // test copying
    UIndirectList<double> uidl1(idl1);
    IndirectList<double> idl2(uidl1);
    IndirectList<double> idl3(idl2);

    printInfo(uidl1);

    idl1.resetAddressing(List<label>());
//    idl2.resetAddressing(List<label>());

    Info<<"after resetAddressing:" << nl << endl;

    printInfo(uidl1);
    printInfo(idl1);
    printInfo(idl2);
    printInfo(idl3);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
