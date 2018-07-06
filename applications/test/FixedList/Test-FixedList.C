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
    Test-FixedList

Description
    Simple tests and examples of use of FixedList

See also
    Foam::FixedList

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "FixedList.H"
#include "IFstream.H"
#include "OFstream.H"
#include "IPstream.H"
#include "OPstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList args(argc, argv);

    FixedList<label, 4> list;
    list[0] = 1;
    list[1] = 2;
    list[2] = 3;
    list[3] = 4;

    Info<< "list:" << list
        << " hash:" << FixedList<label, 4>::Hash<>()(list) << endl;

    label a[4] = {0, 1, 2, 3};
    FixedList<label, 4> list2(a);

    Info<< "list2:" << list2
        << " hash:" << FixedList<label, 4>::Hash<>()(list2) << endl;

    Info<< "list: " << list << nl
        << "list2: " << list2 << endl;
    list.swap(list2);
    Info<< "Swapped via the swap() method" << endl;
    Info<< "list: " << list << nl
        << "list2: " << list2 << endl;

    List<label> list3{0, 1, 2, 3};
    FixedList<label, 4> list4(list3.begin(), list3.end());
    Info<< "list3: " << list3 << nl
        << "list4: " << list4 << endl;

    list4 = {1, 2, 3, 5};
    Info<< "list4: " << list4 << nl;

    FixedList<label, 5> list5{0, 1, 2, 3, 4};
    Info<< "list5: " << list5 << endl;

    List<FixedList<label, 2>> list6{{0, 1}, {2, 3}};
    Info<< "list6: " << list6 << endl;

    if (Pstream::parRun())
    {
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            Serr<< "slave sending to master "
                << Pstream::masterNo() << endl;

            OPstream toMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            FixedList<label, 2> list3;
            list3[0] = 0;
            list3[1] = 1;
            toMaster << list3;
        }
        else
        {
            for
            (
                int slave = Pstream::firstSlave();
                slave <= Pstream::lastSlave();
                slave++
            )
            {
                Serr << "master receiving from slave " << slave << endl;
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);
                FixedList<label, 2> list3(fromSlave);

                Serr<< list3 << endl;
            }
        }
    }

    return 0;
}


// ************************************************************************* //
