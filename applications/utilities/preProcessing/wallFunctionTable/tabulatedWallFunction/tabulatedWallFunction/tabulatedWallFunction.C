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

#include "tabulatedWallFunction.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace tabulatedWallFunctions
    {
        defineTypeNameAndDebug(tabulatedWallFunction, 0);
        defineRunTimeSelectionTable(tabulatedWallFunction, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tabulatedWallFunctions::tabulatedWallFunction::tabulatedWallFunction
(
    const dictionary& dict,
    const polyMesh& mesh,
    const word& name
)
:
    dict_(dict),
    mesh_(mesh),
    coeffDict_(dict.subDict(name + "Coeffs")),
    invertedTableName_(dict.lookup("invertedTableName")),
    invertedTable_(invertedTableName_, mesh_, dict, true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tabulatedWallFunctions::tabulatedWallFunction::~tabulatedWallFunction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::tabulatedWallFunctions::tabulatedWallFunction::write()
{
    if (invertedTable_.log10())
    {
        invertedTable_.note() =
            "U+ as a function of log10(Re) computed using " + type();
    }
    else
    {
        invertedTable_.note() =
            "U+ as a function of Re computed using " + type();
    }

    Info<< "Writing inverted table to\n    " << invertedTable_.objectPath()
        << endl;

    invertedTable_.write();
}


// ************************************************************************* //
