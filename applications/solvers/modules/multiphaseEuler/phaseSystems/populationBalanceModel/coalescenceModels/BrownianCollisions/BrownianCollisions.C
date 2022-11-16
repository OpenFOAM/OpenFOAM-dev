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

#include "BrownianCollisions.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
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

using Foam::constant::physicoChemical::k;
using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::BrownianCollisions::
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
            popBal_.time().timeName(),
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

void Foam::diameterModels::coalescenceModels::BrownianCollisions::precompute()
{
    const volScalarField& T = popBal_.continuousPhase().thermo().T();
    const volScalarField& p = popBal_.continuousPhase().thermo().p();

    lambda_ = k*T/(sqrt(2.0)*pi*p*sqr(sigma_));
}


void
Foam::diameterModels::coalescenceModels::BrownianCollisions::
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
    const volScalarField& mu = popBal_.continuousPhase().thermo().mu();

    const volScalarField Cci
    (
        1 + lambda_/fi.d()*(A1_ + A2_*exp(-A3_*fi.d()/lambda_))
    );

    const volScalarField Ccj
    (
        1 + lambda_/fj.d()*(A1_ + A2_*exp(-A3_*fj.d()/lambda_))
    );

    coalescenceRate += 8*k*T/(3*mu)*(fi.d() + fj.d())*(Cci/fi.d() + Ccj/fj.d());
}


// ************************************************************************* //
