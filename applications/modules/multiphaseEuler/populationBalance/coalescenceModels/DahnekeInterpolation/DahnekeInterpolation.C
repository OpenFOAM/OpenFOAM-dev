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

#include "DahnekeInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "BrownianCollisions.H"
#include "ballisticCollisions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
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

Foam::populationBalance::coalescenceModels::DahnekeInterpolation::
DahnekeInterpolation
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    Brownian_(new BrownianCollisions(popBal, dict)),
    ballistic_(new ballisticCollisions(popBal, dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::populationBalance::coalescenceModels::DahnekeInterpolation::precompute()
{
    Brownian_().precompute();
    ballistic_().precompute();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalance::coalescenceModels::DahnekeInterpolation::rate
(
    const label i,
    const label j
) const
{
    tmp<volScalarField::Internal> tBrownianRate = Brownian_().rate(i, j);
    tmp<volScalarField::Internal> tballisticRate = ballistic_().rate(i, j);

    const volScalarField::Internal KnD(tBrownianRate()/(2*tballisticRate));

    return tBrownianRate*(1 + KnD)/(1 + 2*KnD + 2*sqr(KnD));
}


// ************************************************************************* //
