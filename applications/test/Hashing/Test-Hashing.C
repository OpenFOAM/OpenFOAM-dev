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
    testHashing

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"

#include "stringList.H"
#include "labelList.H"
#include "labelPair.H"
#include "edgeList.H"
#include "triFaceList.H"

#include "Hash.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    IFstream is("hashingTests");


    while (is.good())
    {
        const word listType(is);

        Info<< endl;
        IOobject::writeDivider(Info) << listType << endl;

        if (listType == "stringList")
        {
            Info<< "contiguous = " << contiguous<string>() << endl << endl;

            stringList lst(is);

            forAll(lst, i)
            {
                unsigned hash1 = string::hash()(lst[i]);

                Info<< hex << hash1 << ": " << lst[i] << endl;
            }

        }
        else if (listType == "labelList")
        {
            Info<<"contiguous = " << contiguous<label>() << endl << endl;

            labelList lst(is);

            forAll(lst, i)
            {
                // direct value
                unsigned hash1 = Hash<label>()(lst[i]);

                // hashed byte-wise
                unsigned hash2 = Hash<label>()(lst[i], 0);

                Info<< hex << hash1
                    << " (seeded: " << hash2 << ")"
                    << ": " << dec << lst[i] << endl;
            }

            if (contiguous<label>())
            {
                unsigned hash3 = Hasher
                (
                    lst.cdata(),
                    lst.size() * sizeof(label)
                );

                Info<<"contiguous hashed value " << hex << hash3 << endl;
            }
        }
        else if (listType == "labelListList")
        {
            List<List<label>> lst(is);

            forAll(lst, i)
            {
                unsigned hash1 = Hasher
                (
                    lst[i].cdata(),
                    lst[i].size() * sizeof(label)
                );

                Info<< hex << hash1
                    << ": " << dec << lst[i] << endl;
            }

        }
        else if (listType == "edgeList")
        {
            Info<<"contiguous = " << contiguous<edge>() << endl << endl;

            edgeList lst(is);

            forAll(lst, i)
            {
                unsigned hash1 = Hash<edge>()(lst[i]);

                // as FixedList
                unsigned hash2 = labelPair::Hash<>()(lst[i]);

                Info<< hex << hash1 << " (as FixedList: " << hash2
                    << "): " << dec << lst[i] << endl;
            }
        }
        else if (listType == "triFaceList")
        {
            Info<<"contiguous = " << contiguous<triFace>() << endl << endl;

            triFaceList lst(is);

            forAll(lst, i)
            {
                // direct value
                unsigned hash1 = Hash<triFace>()(lst[i]);
                unsigned hash2 = FixedList<label, 3>::Hash<>()(lst[i]);

                Info<< hex << hash1 << " (as FixedList: " << hash2
                    << "): " << dec << lst[i] << endl;
            }
        }
        else
        {
            Info<< "unknown type: " << listType << endl;
        }

    }

    return 0;
}


// ************************************************************************* //
