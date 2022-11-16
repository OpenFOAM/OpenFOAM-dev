/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable(breakupModel, exponential, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::breakupModels::exponential::exponential
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    breakupModel(popBal, dict),
    exponent_(dict.lookup<scalar>("exponent")),
    C_(dict.lookup<scalar>("C"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::breakupModels::exponential::setBreakupRate
(
    volScalarField& breakupRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];

    breakupRate.primitiveFieldRef() =
        C_*exp(exponent_*fi.x().value());
}


// ************************************************************************* //
