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

Application

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "boolList.H"
#include "PackedBoolList.H"
#include "HashSet.H"
#include "StaticHashTable.H"
#include "cpuTime.H"
#include <vector>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
    const label n = 1000000;
    const label nIters = 1000;

    unsigned int sum = 0;

    PackedBoolList packed(n, 1);
    boolList unpacked(n, true);
    std::vector<bool> stlVector(n, true);

    labelHashSet emptyHash;
    labelHashSet fullHash(1000);
    for (label i = 0; i < n; i++)
    {
        fullHash.insert(i);
    }

    // fullStaticHash is really slow
    // give it lots of slots to help
    StaticHashTable<nil, label, Hash<label> > emptyStaticHash;
    StaticHashTable<nil, label, Hash<label> > fullStaticHash(100000);
    for (label i = 0; i < n; i++)
    {
        fullStaticHash.insert(i, nil());
    }

    emptyHash.printInfo(Info);
    fullHash.printInfo(Info);
    emptyStaticHash.printInfo(Info);
    fullStaticHash.printInfo(Info);


    cpuTime timer;

    for (label iter = 0; iter < nIters; ++iter)
    {
        packed.resize(40);
        packed.shrink();
        packed.resize(n, 1);
    }
    Info<< "resize/shrink/resize:" << timer.cpuTimeIncrement() << " s\n\n";

    // set every other bit on:
    Info<< "set every other bit on and count\n";
    packed.storage() = 0xAAAAAAAAu;

    // Count packed
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            sum += packed[i];
        }
    }
    Info<< "Counting brute-force:" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Count packed
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        sum += packed.count();
    }
    Info<< "Counting via count():" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Dummy addition
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += i + 1;
        }
    }
    Info<< "Dummy loop:" << timer.cpuTimeIncrement() << " s" << endl;
    Info<< "  sum " << sum << endl;

    //
    // Read
    //

    // Read stl
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        for (unsigned int i = 0; i < stlVector.size(); i++)
        {
            sum += stlVector[i];
        }
    }
    Info<< "Reading stl:" << timer.cpuTimeIncrement() << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Read unpacked
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += unpacked[i];
        }
    }
    Info<< "Reading unpacked:" << timer.cpuTimeIncrement() << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Read packed
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            sum += packed.get(i);
        }
    }
    Info<< "Reading packed using get:" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Read packed
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            sum += packed[i];
        }
    }
    Info<< "Reading packed using reference:" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Read via iterator
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAllIter(PackedBoolList, packed, it)
        {
            sum += it;
        }
    }
    Info<< "Reading packed using iterator:" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Read via iterator
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAllConstIter(PackedBoolList, packed, cit)
        {
            sum += cit();
        }
    }
    Info<< "Reading packed using const_iterator():" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Read empty hash
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += emptyHash.found(i);
        }
    }
    Info<< "Reading empty labelHashSet:" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Read full hash
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += fullHash.found(i);
        }
    }
    Info<< "Reading full labelHashSet:" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;


    // Read empty static hash
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += emptyStaticHash.found(i);
        }
    }
    Info<< "Reading empty StaticHash:" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;

#if 0
    // we can skip this test - it is usually quite slow
    // Read full static hash
    sum = 0;
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            sum += fullStaticHash.found(i);
        }
    }
    Info<< "Reading full StaticHash:" << timer.cpuTimeIncrement()
        << " s" << endl;
    Info<< "  sum " << sum << endl;
#endif

    Info<< "Starting write tests" << endl;

    //
    // Write
    //

    // Write stl
    for (label iter = 0; iter < nIters; ++iter)
    {
        for (unsigned int i = 0; i < stlVector.size(); i++)
        {
            stlVector[i] = true;
        }
    }
    Info<< "Writing stl:" << timer.cpuTimeIncrement() << " s" << endl;

    // Write unpacked
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(unpacked, i)
        {
            unpacked[i] = true;
        }
    }
    Info<< "Writing unpacked:" << timer.cpuTimeIncrement() << " s" << endl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            packed[i] = 1;
        }
    }
    Info<< "Writing packed using reference:" << timer.cpuTimeIncrement()
        << " s" << endl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAll(packed, i)
        {
            packed.set(i, 1);
        }
    }
    Info<< "Writing packed using set:" << timer.cpuTimeIncrement()
        << " s" << endl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        forAllIter(PackedBoolList, packed, it)
        {
            it() = 1;
        }
    }
    Info<< "Writing packed using iterator:" << timer.cpuTimeIncrement()
        << " s" << endl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        packed = 0;
    }
    Info<< "Writing packed uniform 0:" << timer.cpuTimeIncrement()
        << " s" << endl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        packed = 1;
    }
    Info<< "Writing packed uniform 1:" << timer.cpuTimeIncrement()
        << " s" << endl;


    PackedList<3> oddPacked(n, 3);

    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        packed = 0;
    }
    Info<< "Writing packed<3> uniform 0:" << timer.cpuTimeIncrement()
        << " s" << endl;


    // Write packed
    for (label iter = 0; iter < nIters; ++iter)
    {
        packed = 1;
    }
    Info<< "Writing packed<3> uniform 1:" << timer.cpuTimeIncrement()
        << " s" << endl;


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
