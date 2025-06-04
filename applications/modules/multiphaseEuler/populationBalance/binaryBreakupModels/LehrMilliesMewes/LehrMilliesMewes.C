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

#include "LehrMilliesMewes.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace binaryBreakupModels
{
    defineTypeNameAndDebug(LehrMilliesMewes, 0);
    addToRunTimeSelectionTable
    (
        binaryBreakupModel,
        LehrMilliesMewes,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::binaryBreakupModels::LehrMilliesMewes::
LehrMilliesMewes
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binaryBreakupModel(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::binaryBreakupModels::LehrMilliesMewes::
addToBinaryBreakupRate
(
    volScalarField::Internal& binaryBreakupRate,
    const label i,
    const label j
)
{
    using Foam::constant::mathematical::pi;

    const dimensionedScalar& dSphi = popBal_.dSph(i);
    const dimensionedScalar& dSphj = popBal_.dSph(j);

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(j));
    const volScalarField::Internal& sigma = tsigma();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();

    volScalarField::Internal L
    (
        pow(sigma/rhoc, 3.0/5.0)/pow(epsilonc, 2.0/5.0)
    );

    // Reset of dimension to pure length to avoid problems in transcendental
    // functions due to small exponents
    L.dimensions().reset(dimLength);

    const volScalarField::Internal T
    (
        pow(sigma/rhoc, 2.0/5.0)/pow(epsilonc, 3.0/5.0)
    );

    binaryBreakupRate +=
        0.5*pow(dSphj/L, 5.0/3.0)
       *exp(-sqrt(2.0)/pow3(dSphj/L))
       *6/pow(pi, 1.5)/pow3(dSphi/L)
       *exp(-9.0/4.0*sqr(log(pow(2.0, 0.4)*dSphi/L)))
       /max(1 + erf(1.5*log(pow(2.0, 1.0/15.0)*dSphj/L)), small)
       /(T*pow3(L));
}


// ************************************************************************* //
