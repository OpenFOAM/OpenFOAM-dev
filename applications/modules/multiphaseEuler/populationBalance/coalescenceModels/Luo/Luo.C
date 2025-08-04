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

#include "Luo.H"
#include "phaseSystem.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "dispersedVirtualMassModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(Luo, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        Luo,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::coalescenceModels::Luo::
Luo
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    beta_("beta", dimless, dict, 2.05),
    C1_("C1", dimless, dict, 1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::coalescenceModels::Luo::addToCoalescenceRate
(
    volScalarField::Internal& coalescenceRate,
    const label i,
    const label j
)
{
    using Foam::constant::mathematical::pi;

    const dimensionedScalar& dSphi = popBal_.dSph(i);
    const dimensionedScalar& dSphj = popBal_.dSph(j);

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(i));
    const volScalarField::Internal& sigma = tsigma();

    const dispersedPhaseInterface interface
    (
        popBal_.phases()[i],
        popBal_.continuousPhase()
    );

    if
    (
       !popBal_.fluid().foundInterfacialModel
        <
            virtualMassModels::dispersedVirtualMassModel
        >(interface)
    )
    {
        FatalErrorInFunction
            << "A virtual mass model for " << popBal_.phases()[i].name()
            << " in " << popBal_.continuousPhase().name() << " is not "
            << "specified. This is required by the Luo coalescence model."
            << exit(FatalError);
    }

    tmp<volScalarField> tCvm =
        popBal_.fluid().lookupInterfacialModel
        <
            virtualMassModels::dispersedVirtualMassModel
        >(interface).Cvm();
    const volScalarField::Internal& Cvm = tCvm();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();

    const dimensionedScalar xi = dSphi/dSphj;

    const volScalarField::Internal uij
    (
        sqrt(beta_)*cbrt(epsilonc*dSphi)*sqrt(1 + pow(xi, -2.0/3.0))
    );

    coalescenceRate +=
        pi/4*sqr(dSphi + dSphj)*uij
       *exp
        (
          - C1_
           *sqrt(0.75*(1 + sqr(xi))*(1 + pow3(xi)))
           /(sqrt(popBal_.phases()[i].rho()()/rhoc + Cvm)*pow3(1 + xi))
           *sqrt(rhoc*dSphi*sqr(uij)/sigma)
        );
}


// ************************************************************************* //
