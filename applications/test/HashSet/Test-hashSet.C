/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "HashSet.H"
#include "Map.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    wordHashSet setA(0);
    HashTable<label, word> tableA;

    HashTable<nil> tableB;
    Map<label> mapA;

    setA.insert("kjhk");
    setA.insert("kjhk2");

    tableA.insert("value1", 1);
    tableA.insert("value2", 2);
    tableA.insert("value3", 3);

    tableB.insert("value4", nil());
    tableB.insert("value5", nil());
    tableB.insert("value6", nil());

    mapA.set(1, 1);
    mapA.set(2, 2);
    mapA.set(3, 3);
    mapA.set(4, 4);

    Info<< setA << endl;
    Info<< tableA << endl;
    Info<< mapA << endl;

    Info<< "create from HashSet: ";
    Info<< wordHashSet(setA) << endl;
    Info<< "create from HashTable<T>: ";
    Info<< wordHashSet(tableA) << endl;
    Info<< "create from HashTable<nil>: ";
    Info<< wordHashSet(tableB) << endl;

    Info<< "create from Map<label>: ";
    Info<< labelHashSet(mapA) << endl;

    Info<<"combined toc: "
        << (wordHashSet(setA) | wordHashSet(tableA) | wordHashSet(tableB))
        << nl;


    labelHashSet setB(1);
    setB.insert(11);
    setB.insert(42);

    Info<< "setB : " << setB << endl;

    labelHashSet setC(1);
    setC.insert(2008);
    setC.insert(1984);

    Info<< "setC : " << setC << endl;

    labelHashSet setD(1);
    setD.insert(11);
    setD.insert(100);
    setD.insert(49);
    setD.insert(36);
    setD.insert(2008);

    Info<< "setD : " << setD << endl;

    Info<< "setB == setC: " << (setB == setC) << endl;
    Info<< "setC != setD: " << (setC != setD) << endl;

    // test operations
    setB += setC;
    Info<< "setB += setC : " << setB << endl;

    setB &= setD;
    Info<< "setB &= setD : " << setB << endl;

    Info<< "setB : " << setB << endl;
    Info<< "setC : " << setC << endl;
    Info<< "setD : " << setD << endl;
    Info<< "setB ^ setC ^ setD : " << (setB ^ setC ^ setD) << endl;

    // test operator[]

    Info<< "setD : " << setD << endl;
    if (setD[0])
    {
        Info<< "setD has 0" << endl;
    }
    else
    {
        Info<< "setD has no 0" << endl;
    }


    if (setD[11])
    {
        Info<< "setD has 11" << endl;
    }
    else
    {
        Info<< "setD has no 0" << endl;
    }

    Info<< "setD : " << setD << endl;

    // this doesn't work (yet?)
    // setD[12] = true;

    List<label> someLst(10);
    forAll(someLst, elemI)
    {
        someLst[elemI] = elemI*elemI;
    }

    label added = setD.set(someLst);
    Info<< "added " << added << " from " << someLst.size() << endl;
    Info<< "setD : " << setD << endl;


    return 0;
}


// ************************************************************************* //
