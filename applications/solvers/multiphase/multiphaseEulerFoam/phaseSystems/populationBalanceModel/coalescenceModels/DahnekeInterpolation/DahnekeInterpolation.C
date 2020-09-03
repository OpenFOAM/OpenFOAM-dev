/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "DahnekeInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "BrownianCollisions.H"
#include "ballisticCollisions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(DahnekeInterpolation, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        DahnekeInterpolation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::DahnekeInterpolation::
DahnekeInterpolation
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    Brownian_(new BrownianCollisions(popBal, dict)),
    BrownianRate_
    (
        IOobject
        (
            "BrownianCollisionRate",
            popBal_.mesh().time().timeName(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar("BrownianCollisionRate", dimVolume/dimTime, Zero)
    ),
    ballistic_(new ballisticCollisions(popBal, dict)),
    ballisticRate_
    (
        IOobject
        (
            "ballisticCollisionRate",
            popBal_.mesh().time().timeName(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar("ballisticCollisionRate", dimVolume/dimTime, Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceModels::DahnekeInterpolation::precompute()
{
    Brownian_().precompute();
}


void
Foam::diameterModels::coalescenceModels::DahnekeInterpolation::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const label i,
    const label j
)
{
    BrownianRate_ = Zero;
    ballisticRate_ = Zero;

    Brownian_().addToCoalescenceRate(BrownianRate_, i, j);
    ballistic_().addToCoalescenceRate(ballisticRate_, i, j);

    const volScalarField KnD(BrownianRate_/(2*ballisticRate_));

    coalescenceRate += BrownianRate_*(1 + KnD)/(1 + 2*KnD + 2*sqr(KnD));
}


// ************************************************************************* //
