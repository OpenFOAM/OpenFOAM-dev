/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "Liao.H"
#include "fvcGrad.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{
    defineTypeNameAndDebug(Liao, 0);
    addToRunTimeSelectionTable
    (
        binaryBreakupModel,
        Liao,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::binaryBreakupModels::Liao::Liao
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binaryBreakupModel(popBal, dict),
    LiaoBase(popBal, dict),
    BTurb_(dict.lookupOrDefault("BTurb", dimensionedScalar(dimless, 1))),
    BShear_(dict.lookupOrDefault("BShear", dimensionedScalar(dimless, 1))),
    BEddy_(dict.lookupOrDefault("BEddy", dimensionedScalar(dimless, 1))),
    BFric_(dict.lookupOrDefault("BFric", dimensionedScalar(dimless, 0.25))),
    turbulence_(dict.lookup("turbulence")),
    laminarShear_(dict.lookup("laminarShear")),
    turbulentShear_(dict.lookup("turbulentShear")),
    interfacialFriction_(dict.lookup("interfacialFriction"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::binaryBreakupModels::Liao::precompute()
{
    LiaoBase::precompute();
}


void Foam::diameterModels::binaryBreakupModels::Liao::
addToBinaryBreakupRate
(
    volScalarField& binaryBreakupRate,
    const label i,
    const label j
)
{
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    dimensionedScalar dk(cbrt(pow3(fj.dSph()) - pow3(fi.dSph())));

    const volScalarField tauCrit1
    (
        6*popBal_.sigmaWithContinuousPhase(fj.phase())/fj.dSph()
       *(sqr(fi.dSph()/fj.dSph()) + sqr(dk/fj.dSph()) - 1)
    );

    const volScalarField tauCrit2
    (
        popBal_.sigmaWithContinuousPhase(fj.phase())/min(dk, fi.dSph())
    );

    const volScalarField tauCrit(max(tauCrit1, tauCrit2));

    if (turbulence_)
    {
        const volScalarField tauTurb
        (
            pos(fj.dSph() - kolmogorovLengthScale_)*BTurb_*continuousPhase.rho()
           *sqr(cbrt(popBal_.continuousTurbulence().epsilon()*fj.dSph()))
        );

        binaryBreakupRate +=
            pos(tauTurb - tauCrit)*1/fj.dSph()
           *sqrt(mag(tauTurb - tauCrit)/continuousPhase.rho())/fj.x();
    }

    if (laminarShear_)
    {
        const volScalarField tauShear
        (
            BShear_*continuousPhase.fluidThermo().mu()*shearStrainRate_
        );

        binaryBreakupRate +=
            pos(tauShear - tauCrit)*1/fj.dSph()
           *sqrt(mag(tauShear - tauCrit)/continuousPhase.rho())/fj.x();
    }

    if (turbulentShear_)
    {
        const volScalarField tauEddy
        (
            pos0(kolmogorovLengthScale_ - fj.dSph())
           *BEddy_*continuousPhase.fluidThermo().mu()*eddyStrainRate_
        );

        binaryBreakupRate +=
            pos(tauEddy - tauCrit)*1/fj.dSph()
           *sqrt(mag(tauEddy - tauCrit)/continuousPhase.rho())/fj.x();
    }

    if (interfacialFriction_)
    {
        const volScalarField tauFric
        (
            BFric_*0.5*continuousPhase.rho()*sqr(uTerminal_[j])*Cd_[j]
        );

        binaryBreakupRate +=
            pos(tauFric - tauCrit)*1/fj.dSph()
           *sqrt(mag(tauFric - tauCrit)/continuousPhase.rho())/fj.x();
    }
}


// ************************************************************************* //
