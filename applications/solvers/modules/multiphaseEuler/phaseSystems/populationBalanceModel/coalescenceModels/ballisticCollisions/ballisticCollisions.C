/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2022 OpenFOAM Foundation
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

#include "ballisticCollisions.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(ballisticCollisions, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        ballisticCollisions,
        dictionary
    );
}
}
}

using Foam::constant::physicoChemical::k;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::ballisticCollisions::
ballisticCollisions
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::coalescenceModels::ballisticCollisions::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const label i,
    const label j
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    const volScalarField& T = popBal_.continuousPhase().thermo().T();
    const volScalarField& rhoP = fi.phase().rho();

    coalescenceRate +=
        sqrt(3*k*T/rhoP)*sqr(fi.d() + fj.d())
       *sqrt(1/pow3(fi.d()) + 1/pow3(fj.d()));
}


// ************************************************************************* //
