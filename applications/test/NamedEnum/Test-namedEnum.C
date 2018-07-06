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

#include "NamedEnum.H"
#include "IOstreams.H"

using namespace Foam;

class namedEnumTest
{
public:

    enum options
    {
        a,
        b,
        c
    };

    static const Foam::NamedEnum<options, 3> namedEnum;
};


template<>
const char* Foam::NamedEnum<namedEnumTest::options, 3>::names[] =
{
    "a",
    "b",
    "c"
};

const Foam::NamedEnum<namedEnumTest::options, 3> namedEnumTest::namedEnum;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< namedEnumTest::namedEnum["a"] << endl;
    Info<< namedEnumTest::namedEnum[namedEnumTest::a] << endl;

    namedEnumTest::options hmm(namedEnumTest::namedEnum.read(Sin));
    Info<< namedEnumTest::namedEnum[hmm] << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
