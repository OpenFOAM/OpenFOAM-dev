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
#include "brownianCollisions.H"
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
    brownian_(new brownianCollisions(popBal, dict)),
    brownianCollisionRate_
    (
        IOobject
        (
            "brownianCollisionRate",
            popBal_.mesh().time().timeName(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar("brownianCollisionRate", dimVolume/dimTime, 0.0)
    ),
    ballistic_(new ballisticCollisions(popBal, dict)),
    ballisticCollisionRate_
    (
        IOobject
        (
            "ballisticCollisionRate",
            popBal_.mesh().time().timeName(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar("ballisticCollisionRate", dimVolume/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::coalescenceModels::DahnekeInterpolation::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const label i,
    const label j
)
{
    brownianCollisionRate_ = Zero;
    ballisticCollisionRate_ = Zero;

    brownian_().addToCoalescenceRate(brownianCollisionRate_, i, j);
    ballistic_().addToCoalescenceRate(ballisticCollisionRate_, i, j);

    const volScalarField KnD
    (
        brownianCollisionRate_/(2.0*ballisticCollisionRate_)
    );

    coalescenceRate +=
        brownianCollisionRate_*(1.0 + KnD)/(1.0 + 2.0*KnD + 2.0*sqr(KnD));
}


// ************************************************************************* //
