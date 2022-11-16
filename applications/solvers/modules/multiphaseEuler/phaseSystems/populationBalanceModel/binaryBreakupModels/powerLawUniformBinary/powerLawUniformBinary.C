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

#include "powerLawUniformBinary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{
    defineTypeNameAndDebug(powerLawUniformBinary, 0);
    addToRunTimeSelectionTable
    (
        binaryBreakupModel,
        powerLawUniformBinary,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::binaryBreakupModels::powerLawUniformBinary::
powerLawUniformBinary
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binaryBreakupModel(popBal, dict),
    power_(dict.lookup<scalar>("power"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::binaryBreakupModels::powerLawUniformBinary::
addToBinaryBreakupRate
(
    volScalarField& binaryBreakupRate,
    const label i,
    const label j
)
{
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    binaryBreakupRate.primitiveFieldRef() +=
        pow(fj.x().value(), power_)*2/fj.x().value();
}


// ************************************************************************* //
