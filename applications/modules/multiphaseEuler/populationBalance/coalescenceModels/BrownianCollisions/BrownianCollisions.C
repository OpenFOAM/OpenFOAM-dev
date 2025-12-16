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

#include "BrownianCollisions.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(BrownianCollisions, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        BrownianCollisions,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::coalescenceModels::BrownianCollisions::
BrownianCollisions
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    A1_(dict.lookupOrDefault<scalar>("A1", 2.514)),
    A2_(dict.lookupOrDefault<scalar>("A2", 0.8)),
    A3_(dict.lookupOrDefault<scalar>("A3", 0.55)),
    sigma_("sigma", dimLength, dict),
    lambda_
    (
        IOobject
        (
            "lambda",
            popBal_.time().name(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar
        (
            "lambda",
            dimLength,
            Zero
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::coalescenceModels::BrownianCollisions::
precompute()
{
    using Foam::constant::physicoChemical::k;
    using Foam::constant::mathematical::pi;

    const volScalarField::Internal& p =
        popBal_.continuousPhase().fluidThermo().p();

    const volScalarField::Internal& Tc = popBal_.continuousPhase().thermo().T();

    lambda_ = k*Tc/(sqrt(2.0)*pi*p*sqr(sigma_));
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalance::coalescenceModels::BrownianCollisions::rate
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

    tmp<volScalarField> tmuc(popBal_.continuousPhase().fluidThermo().mu());
    const volScalarField::Internal& muc = tmuc();

    const volScalarField::Internal Cci
    (
        1 + lambda_/di*(A1_ + A2_*exp(-A3_*di/lambda_))
    );

    const volScalarField::Internal Ccj
    (
        1 + lambda_/dj*(A1_ + A2_*exp(-A3_*dj/lambda_))
    );

    return 8*k*Tc/(3*muc)*(di + dj)*(Cci/di + Ccj/dj);
}


// ************************************************************************* //
