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
    testSHA1

Description


\*---------------------------------------------------------------------------*/

#include "OSHA1stream.H"
#include "IStringStream.H"
#include "dictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char * argv[])
{
    SHA1 sha;
    SHA1Digest shaDig;

    std::string str("The quick brown fox jumps over the lazy dog");
    Info<< shaDig << nl;
    Info<< SHA1("The quick brown fox jumps over the lazy dog") << nl;

    sha.append("The quick brown fox jumps over the lazy dog");
    Info<< sha << nl;

    sha.clear();
    sha.append("The quick brown fox jumps over the lazy dog");
    shaDig = sha;

    sha.append("\n");
    Info<< sha << nl;
    Info<< shaDig << nl;

    if (sha == shaDig)
    {
        Info<<"SHA1 digests are identical\n";
    }
    else
    {
        Info<<"SHA1 digests are different\n";
    }
    Info<<"lhs:" << sha << " rhs:" << shaDig << endl;

    // start over:
    sha.clear();
    sha.append(str);

    SHA1 sha_A = sha;

    sha.append("\n");

    Info<< "digest1: " << sha_A << nl;
    Info<< "digest2: " << sha << nl;

    // start over:
    sha.clear();
    sha.append("\"");
    sha.append(str);
    sha.append("\"");

    Info<< "digest3: " << sha << nl;

    // try the output buffer interface
    {
        OSHA1stream os;

        os  << str;
        Info<< os.digest() << endl;

        os  << str;
        Info<< os.digest() << endl;

        os.rewind();
        os  << "The quick brown fox jumps over the lazy dog";
        Info<< os.digest() << endl;

    }

    {
        dictionary dict
        (
            IStringStream
            (
                "parent { Default_Boundary_Region { type zeroGradient; } }"
                "inlet_1 { value inlet_1; }"
                "inlet_2 { value inlet_2; }"
                "inlet_3 { value inlet_3; }"
                "\"inlet_.*\" { value XXX; }"
            ) ()
        );

        Info<< "dict:" << endl;
        dict.write(Info, false);

        dictionary dict2(dict);

        OSHA1stream os;
        dict.write(os, false);
        Info<< os.digest() << endl;

        Info<< dict2.digest() << endl;
    }


    return 0;
}
