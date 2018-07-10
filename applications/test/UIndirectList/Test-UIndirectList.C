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

#include "UIndirectList.H"
#include "DynamicList.H"
#include "IOstreams.H"
#include "ListOps.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<double> completeList(10);

    forAll(completeList, i)
    {
        completeList[i] = 0.1*i;
    }

    List<label> addresses(5);
    addresses[0] = 1;
    addresses[1] = 0;
    addresses[2] = 7;
    addresses[3] = 8;
    addresses[4] = 5;

    UIndirectList<double> idl(completeList, addresses);

    Info<< idl << "\n";

    idl[1] = -666;

    Info<< "idl[1] changed: " << idl << endl;

    idl = -999;

    Info<< "idl changed: " << idl << endl;

    UIndirectList<double> idl2(idl);

    Info<< "idl2: " << idl2 << endl;


    {
        List<double> ident(idl.size());

        forAll(ident, i)
        {
            ident[i] = ident.size() - i;
        }
        idl = ident;
    }

    Info<< "idl assigned from UList: " << idl << endl;

    // test List operations

    List<double> flatList(UIndirectList<double>(completeList, addresses));
    Info<< "List constructed from UIndirectList: " << flatList << endl;

    flatList = UIndirectList<double>(completeList, addresses);
    Info<< "List assigned from UIndirectList: " << flatList << endl;

    flatList.append(UIndirectList<double>(completeList, addresses));
    Info<< "List::append(UIndirectList): " << flatList << endl;


    DynamicList<double> dynList(UIndirectList<double>(completeList, addresses));
    Info<< "DynamicList constructed from UIndirectList: " << dynList << endl;

    dynList.append(UIndirectList<double>(completeList, addresses));
    Info<< "DynamicList::append(UIndirectList): " << dynList << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
