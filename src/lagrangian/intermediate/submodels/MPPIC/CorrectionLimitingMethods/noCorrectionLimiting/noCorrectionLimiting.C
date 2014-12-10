/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "noCorrectionLimiting.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace CorrectionLimitingMethods
{
    defineTypeNameAndDebug(noCorrectionLimiting, 0);

    addToRunTimeSelectionTable
    (
        CorrectionLimitingMethod,
        noCorrectionLimiting,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CorrectionLimitingMethods::noCorrectionLimiting::noCorrectionLimiting
(
    const dictionary& dict
)
:
    CorrectionLimitingMethod(dict)
{}


Foam::CorrectionLimitingMethods::noCorrectionLimiting::noCorrectionLimiting
(
    const noCorrectionLimiting& cl
)
:
    CorrectionLimitingMethod(cl)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CorrectionLimitingMethods::noCorrectionLimiting::~noCorrectionLimiting()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector
Foam::CorrectionLimitingMethods::noCorrectionLimiting::limitedVelocity
(
    const vector uP,
    const vector dU,
    const vector uMean
) const
{
    return dU;
}


// ************************************************************************* //
