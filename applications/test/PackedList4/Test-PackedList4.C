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
#include "IStringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    PackedBoolList list1(20);
    // set every third one on
    forAll(list1, i)
    {
        list1[i] = !(i % 3);
    }

    Info<< "\nalternating bit pattern\n";
    list1.printInfo(Info, true);

    PackedBoolList list2 = ~list1;

    Info<< "\ncomplementary bit pattern\n";
    list2.printBits(Info);

    // set every other on
    forAll(list2, i)
    {
        list2[i] = !(i % 2);
    }

    Info<< "\nalternating bit pattern\n";
    list2.printBits(Info);

    list2.resize(28, false);
    list2.resize(34, true);
    list2.resize(40, false);
    for (label i=0; i < 4; ++i)
    {
        list2[i] = true;
    }

    Info<< "\nresized with false, 6 true + 6 false, bottom 4 bits true\n";
    list2.printInfo(Info, true);

    labelList list2Labels = list2.used();

    Info<< "\noperator|\n";

    list1.printBits(Info);
    list2.printBits(Info);
    Info<< "==\n";
    (list1 | list2).printBits(Info);

    Info<< "\noperator& : does trim\n";
    (list1 & list2).printBits(Info);

    Info<< "\noperator^\n";
    (list1 ^ list2).printBits(Info);


    Info<< "\noperator|=\n";
    {
        PackedBoolList list3 = list1;
        (list3 |= list2).printBits(Info);
    }

    Info<< "\noperator|= with labelUList\n";
    {
        PackedBoolList list3 = list1;
        (list3 |= list2Labels).printBits(Info);
    }

    Info<< "\noperator&=\n";
    {
        PackedBoolList list3 = list1;
        (list3 &= list2).printBits(Info);
    }

    Info<< "\noperator+=\n";
    {
        PackedBoolList list3 = list1;
        (list3 += list2).printBits(Info);
    }

    Info<< "\noperator+= with labelUList\n";
    {
        PackedBoolList list3 = list1;
        (list3 += list2Labels).printBits(Info);
    }

    Info<< "\noperator-=\n";
    {
        PackedBoolList list3 = list1;
        (list3 -= list2).printBits(Info);
    }

    Info<< "\noperator-= with labelUList\n";
    {
        PackedBoolList list3 = list1;
        (list3 -= list2Labels).printBits(Info);
    }

    PackedBoolList list4
    (
        IStringStream
        (
            "(1 n 1 n 1 n 1 1 off 0 0 f f 0 y yes y true y false on t)"
        )()
    );

    Info<< "\ntest Istream constructor\n";

    list4.printInfo(Info, true);
    Info<< list4 << " indices: " << list4.used()() <<endl;

    Info<< "\nassign from labelList\n";
    list4 = labelList
    (
        IStringStream
        (
            "(0 1 2 3 12 13 14 19 20 21)"
        )()
    );

    list4.printInfo(Info, true);
    Info<< list4 << " indices: " << list4.used()() <<endl;

    Info<< "\nassign from indices\n";
    list4.read
    (
        IStringStream
        (
            "{0 1 2 3 12 13 14 19 20 21}"
        )()
    );


    list4.printInfo(Info, true);
    Info<< list4 << " indices: " << list4.used()() <<endl;

    List<bool> boolLst(list4.size());
    forAll(list4, i)
    {
        boolLst[i] = list4[i];
    }

    Info<< "List<bool>: " << boolLst <<endl;


    // check roundabout assignments
    PackedList<2> pl2
    (
        IStringStream
        (
            "{(0 3)(1 3)(2 3)(3 3)(12 3)(13 3)(14 3)(19 3)(20 3)(21 3)}"
        )()
    );

    Info<< "roundabout assignment: " << pl2 << endl;

    list4.clear();
    forAll(pl2, i)
    {
        list4[i] = pl2[i];
    }

    list4.write(Info, true) << endl;

    list4.writeEntry("PackedBoolList", Info);

    return 0;
}


// ************************************************************************* //
