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
    prefixOSstreamTest

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "IOstreams.H"
#include "IStringStream.H"
#include "scalar.H"
#include "vector.H"
#include "ListOps.H"
#include "Pstream.H"
#include "argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList args(argc, argv);

    // Pout.prefix() = '[' + name(Pstream::myProcNo()) + "] ";

    List<vector> list(IStringStream("1 ((0 1 2))")());
    Info<< list << endl;

    List<vector> list2
    (
        IStringStream
        (
            "(\
                 (0 1 2)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
             )"
        )()
    );
    Pout<< list2 << endl;

    Info<< findIndex(list2, vector(3, 4, 5)) << endl;

    return 0;
}


// ************************************************************************* //
