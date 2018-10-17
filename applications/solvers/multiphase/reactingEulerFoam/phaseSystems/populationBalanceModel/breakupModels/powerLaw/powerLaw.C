/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "powerLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace breakupModels
{
    defineTypeNameAndDebug(powerLaw, 0);
    addToRunTimeSelectionTable(breakupModel, powerLaw, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::breakupModels::powerLaw::powerLaw
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    breakupModel(popBal, dict),
    power_(readScalar(dict.lookup("power")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::breakupModels::powerLaw::setBreakupRate
(
    volScalarField& breakupRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];

    breakupRate.primitiveFieldRef() = pow(fi.x().value(), power_);
}


// ************************************************************************* //
