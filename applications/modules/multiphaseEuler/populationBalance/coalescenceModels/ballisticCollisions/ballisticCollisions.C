/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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
namespace populationBalance
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::coalescenceModels::ballisticCollisions::
ballisticCollisions
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalance::coalescenceModels::ballisticCollisions::rate
(
    const label i,
    const label j
) const
{
    using Foam::constant::physicoChemical::k;

    tmp<volScalarField> tdi = popBal_.d(i);
    const volScalarField::Internal& di = tdi();
    tmp<volScalarField> tdj = popBal_.d(j);
    const volScalarField::Internal& dj = tdj();

    const volScalarField::Internal& Tc = popBal_.continuousPhase().thermo().T();

    return
        sqrt(3*k*Tc/popBal_.phases()[i].rho()())
       *sqr(di + dj)
       *sqrt(1/pow3(di) + 1/pow3(dj));
}


// ************************************************************************* //
