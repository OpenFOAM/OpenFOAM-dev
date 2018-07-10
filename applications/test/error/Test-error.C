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

#include "error.H"
#include "IOstreams.H"
#include "dictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    FatalError.throwExceptions();

    try
    {
        WarningInFunction << "warning 1" << endl;
        IOWarningInFunction(Serr) << "warning 2" << endl;

        dictionary dict;

        IOWarningInFunction(dict) << "warning 3" << endl;

        FatalErrorInFunction << "error 1" << endl;
        FatalErrorInFunction << "error 2" << exit(FatalError);
    }
    catch (Foam::error& fErr)
    {
        Serr<< "Caught Foam error " << fErr << nl << endl;
    }

    return 0;
}


// ************************************************************************* //
