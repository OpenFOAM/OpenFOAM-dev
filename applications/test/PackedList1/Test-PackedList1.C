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

#include "uLabel.H"
#include "IOstreams.H"
#include "PackedBoolList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Info<< "PackedList max_bits() = " << PackedList<>::max_bits() << nl;

    Info<< "\ntest allocation with value\n";
    PackedList<3> list1(5,1);
    list1.printInfo(Info, true);

    Info<< "\ntest assign uniform value\n";
    list1 = 3;
    list1.printInfo(Info, true);

    Info<< "\ntest assign uniform value (with overflow)\n";
    list1 = -1;
    list1.printInfo(Info, true);

    Info<< "\ntest zero\n";
    list1 = 0;
    list1.printInfo(Info, true);

    Info<< "\ntest set() with default argument (max_value)\n";
    list1.set(1);
    list1.set(3);
    list1.printInfo(Info, true);

    Info<< "\ntest unset() with in-range and out-of-range\n";
    list1.unset(3);
    list1.unset(100000);
    list1.printInfo(Info, true);

    Info<< "\ntest assign between references\n";
    list1[2] = 3;
    list1[4] = list1[2];
    list1.printInfo(Info, true);

    Info<< "\ntest assign between references, with chaining\n";
    list1[0] = 1;
    list1[4] = 1;
    list1.printInfo(Info, true);

    Info<< "\ntest assign between references, with chaining and auto-vivify\n";
    list1[1] = 2;
    list1[8] = 2;
    list1[10] = 2;
    list1[14] = 2;
    list1.printInfo(Info, true);


    Info<< "\ntest operator== between references\n";
    if (list1[1] == list1[8])
    {
        Info<< "[1] == [8] (expected)\n";
    }
    else
    {
        Info<< "[1] != [8] (unexpected)\n";
    }

    if (list1[0] != list1[1])
    {
        Info<< "[0] != [1] (expected)\n";
    }
    else
    {
        Info<< "[0] == [1] (unexpected)\n";
    }

    Info<< "\ntest operator== with iterator\n";
    {
        PackedList<3>::iterator iter = list1[1];

        if (iter != list1[8])
        {
            Info<< "iter != [8] (expected)\n";
        }
        else
        {
            Info<< "iter == [8] (unexpected)\n";
        }

        if (*iter != list1[8])
        {
            Info<< "*iter != [8] (unexpected)\n";
        }
        else
        {
            Info<< "*iter == [8] (expected)\n";
        }
    }


    {
        const PackedList<3>& constLst = list1;
        Info<< "\ntest operator[] const with out-of-range index\n";
        constLst.printInfo(Info, true);
        if (constLst[20])
        {
            Info<< "[20] is true (unexpected)\n";
        }
        else
        {
            Info<< "[20] is false (expected) list size should be unchanged "
                << "(const)\n";
        }
        constLst.printInfo(Info, true);

        Info<< "\ntest operator[] non-const with out-of-range index\n";
        if (list1[20])
        {
            Info<< "[20] is true (unexpected)\n";
        }
        else
        {
            Info<< "[20] is false (expected) but list was resized?? "
                << "(non-const)\n";
        }
        list1.printInfo(Info, true);
    }


    Info<< "\ntest operator[] with out-of-range index\n";
    if (!list1[20])
    {
        Info<< "[20] is false, as expected\n";
    }
    list1.printInfo(Info, true);

    Info<< "\ntest resize with value (without reallocation)\n";
    list1.resize(8, list1.max_value());
    list1.printInfo(Info, true);

    Info<< "\ntest flip() function\n";
    list1.flip();
    list1.printInfo(Info, true);

    Info<< "\nre-flip()\n";
    list1.flip();
    list1.printInfo(Info, true);

    Info<< "\ntest set() function\n";
    list1.set(1, 5);
    list1.printInfo(Info, true);

    Info<< "\ntest assign bool\n";
    list1 = false;
    list1.printInfo(Info, true);

    Info<< "\ntest assign bool\n";
    list1 = true;
    list1.printInfo(Info, true);

    Info<< "\ntest resize without value (with reallocation)\n";
    list1.resize(12);
    list1.printInfo(Info, true);

    Info<< "\ntest resize with value (with reallocation)\n";
    list1.resize(25, list1.max_value());
    list1.printInfo(Info, true);

    Info<< "\ntest resize smaller (should not touch allocation)\n";
    list1.resize(8);
    list1.printInfo(Info, true);

    Info<< "\ntest append() operation\n";
    list1.append(2);
    list1.append(3);
    list1.append(4);
    list1.printInfo(Info, true);

    Info<< "\ntest reserve() operation\n";
    list1.reserve(32);
    list1.printInfo(Info, true);

    Info<< "\ntest shrink() operation\n";
    list1.shrink();
    list1.printInfo(Info, true);

    Info<< "\ntest setCapacity() operation\n";
    list1.setCapacity(15);
    list1.printInfo(Info, true);

    Info<< "\ntest setCapacity() operation\n";
    list1.setCapacity(100);
    list1.printInfo(Info, true);

    Info<< "\ntest operator[] assignment\n";
    list1[16] = 5;
    list1.printInfo(Info, true);

    Info<< "\ntest operator[] assignment with auto-vivify\n";
    list1[36] = list1.max_value();
    list1.printInfo(Info, true);

    Info<< "\ntest setCapacity smaller\n";
    list1.setCapacity(24);
    list1.printInfo(Info, true);

    Info<< "\ntest resize much smaller\n";
    list1.resize(150);
    list1.printInfo(Info, true);

    Info<< "\ntest trim\n";
    list1.trim();
    list1.printInfo(Info, true);

    // add in some misc values
    list1[31] = 1;
    list1[32] = 2;
    list1[33] = 3;

    Info<< "\ntest iterator\n";
    PackedList<3>::iterator iter = list1.begin();
    Info<< "begin():";
    iter.printInfo(Info) << "\n";

    Info<< "iterator:" << iter() << "\n";
    iter() = 5;
    iter.printInfo(Info);
    list1.printInfo(Info, true);

    iter = list1[31];
    Info<< "iterator:" << iter() << "\n";
    iter.printInfo(Info);


    Info<< "\ntest get() method\n";
    Info<< "get(10):" << list1.get(10) << " and list[10]:" << list1[10] << "\n";
    list1.printInfo(Info, true);

    Info<< "\ntest iterator indexing\n";
    Info<< "cend() ";
    list1.cend().printInfo(Info) << "\n";

    {
        Info<< "\ntest assignment of iterator\n";
        list1.printInfo(Info, true);
        Info<< "cend()\n";
        list1.end().printInfo(Info);
        PackedList<3>::iterator cit = list1[100];
        Info<< "out-of-range: ";
        cit.printInfo(Info);
        cit = list1[15];
        Info<< "in-range: ";
        cit.printInfo(Info);
        Info<< "out-of-range: ";
        cit = list1[1000];
        cit.printInfo(Info);
    }


    for
    (
        PackedList<3>::iterator cit = list1[30];
        cit != list1.end();
        ++cit
    )
    {
        cit.printInfo(Info);
    }

    Info<< "\ntest operator[] auto-vivify\n";
    Info<< "size:" << list1.size() << "\n";

    const unsigned int val = list1[45];

    Info<< "list[45]:" << val << "\n";
    Info<< "size after read:" << list1.size() << "\n";

    list1[45] = list1.max_value();
    Info<< "size after write:" << list1.size() << "\n";
    Info<< "list[45]:" << list1[45] << "\n";
    list1[49] = list1[100];
    list1.printInfo(Info, true);


    Info<< "\ntest copy constructor + append\n";
    PackedList<3> list2(list1);
    list2.append(4);
    Info<< "source list:\n";
    list1.printInfo(Info, true);
    Info<< "destination list:\n";
    list2.printInfo(Info, true);

    Info<< "\ntest pattern that fills all bits\n";
    PackedList<4> list3(8, 8);

    label pos = list3.size() - 1;

    list3[pos--] = list3.max_value();
    list3[pos--] = 0;
    list3[pos--] = list3.max_value();
    list3.printInfo(Info, true);

    Info<< "removed final value: " << list3.remove() << endl;
    list3.printInfo(Info, true);

    Info<<"list: " << list3 << endl;


    List<bool> list4(16, false);
    {
        // fill with some values
        forAll(list4, i)
        {
            list4[i] = i % 3;
        }

        const UList<bool>& constLst = list4;
        Info<< "\ntest operator[] const with out-of-range index\n";
        Info<< constLst << endl;
        if (constLst[100])
        {
            Info<< "[100] is true (unexpected)\n";
        }
        else
        {
            Info<< "[100] is false (expected) "
                << "list size should be unchanged (const)\n";
        }
        Info<< constLst << endl;
    }


    PackedBoolList listb(list4);

    Info<< "copied from bool list " << endl;
    listb.printInfo(Info, true);

    {
        labelList indices = listb.used();

        Info<< "indices: " << indices << endl;
    }


    Info<< "\n\nDone.\n";

    return 0;
}


// ************************************************************************* //
