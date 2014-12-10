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

\*---------------------------------------------------------------------------*/

#include "simpleMatrix.H"
#include "vector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    simpleMatrix<vector> hmm(3);

    hmm[0][0] = -3.0;
    hmm[0][1] = 10.0;
    hmm[0][2] = -4.0;
    hmm[1][0] = 2.0;
    hmm[1][1] = 3.0;
    hmm[1][2] = 10.0;
    hmm[2][0] = 2.0;
    hmm[2][1] = 6.0;
    hmm[2][2] = 1.0;

    hmm.source()[0] = vector(2.0, 1.0, 3.0);
    hmm.source()[1] = vector(1.0, 4.0, 3.0);
    hmm.source()[2] = vector(0.0, 5.0, 2.0);

    Info<< hmm << endl;
    Info<< hmm.solve() << endl;
    Info<< hmm << endl;
    Info<< hmm.LUsolve() << endl;
    Info<< hmm << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
