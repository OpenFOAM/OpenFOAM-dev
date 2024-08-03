/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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
namespace diameterModels
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

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::Luo::
Luo
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    beta_(dict.lookupOrDefault("beta", dimensionedScalar(dimless, 2.05))),
    C1_(dict.lookupOrDefault("C1", dimensionedScalar(dimless, 1)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceModels::Luo::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const label i,
    const label j
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];
    const phaseModel& continuousPhase = popBal_.continuousPhase();

    if
    (
        popBal_.fluid().foundInterfacialModel
        <virtualMassModels::dispersedVirtualMassModel>
        (
            dispersedPhaseInterface(fi.phase(), popBal_.continuousPhase())
        )
    )
    {
        const virtualMassModels::dispersedVirtualMassModel& vm =
            popBal_.fluid().lookupInterfacialModel
            <virtualMassModels::dispersedVirtualMassModel>
            (
                dispersedPhaseInterface(fi.phase(), popBal_.continuousPhase())
            );

        const dimensionedScalar xi = fi.dSph()/fj.dSph();

        const volScalarField uij
        (
            sqrt(beta_)
           *cbrt(popBal_.continuousTurbulence().epsilon()*fi.dSph())
           *sqrt(1 + pow(xi, -2.0/3.0))
        );

        coalescenceRate +=
            pi/4*sqr(fi.dSph() + fj.dSph())*uij
           *exp
            (
              - C1_
               *sqrt(0.75*(1 + sqr(xi))*(1 + pow3(xi)))
               /(
                    sqrt(fi.phase().rho()/continuousPhase.rho()
                  + vm.Cvm())*pow3(1 + xi)
                )
               *sqrt
                (
                    continuousPhase.rho()*fi.dSph()*sqr(uij)
                   /popBal_.sigmaWithContinuousPhase(fi.phase())
                )
            );
    }
    else
    {
        FatalErrorInFunction
            << "A virtual mass model for " << fi.phase().name() << " in "
            << popBal_.continuousPhase().name() << " is not specified. This is "
            << "required by the Luo coalescence model." << exit(FatalError);
    }
}


// ************************************************************************* //
