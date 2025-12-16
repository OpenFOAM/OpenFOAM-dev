/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace breakupModels
{
    defineTypeNameAndDebug(Liao, 0);
    addToRunTimeSelectionTable(breakupModel, Liao, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::breakupModels::Liao::Liao
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binary(popBal, dict),
    LiaoBase(popBal, dict),
    BTurb_("BTurb", dimless, dict, 1),
    BShear_("BShear", dimless, dict, 1),
    BEddy_("BEddy", dimless, dict, 1),
    BFric_("BFric", dimless, dict, 0.25),
    turbulence_(dict.lookup("turbulence")),
    laminarShear_(dict.lookup("laminarShear")),
    turbulentShear_(dict.lookup("turbulentShear")),
    interfacialFriction_(dict.lookup("interfacialFriction"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::breakupModels::Liao::precompute()
{
    LiaoBase::precompute();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalance::breakupModels::Liao::rate
(
    const label i,
    const label j
) const
{
    const dimensionedScalar& dSphi = popBal_.dSph(i);
    const dimensionedScalar& dSphj = popBal_.dSph(j);
    const dimensionedScalar& vj = popBal_.v(j);

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(i));
    const volScalarField::Internal& sigma = tsigma();

    tmp<volScalarField> tmuc(popBal_.continuousPhase().fluidThermo().mu());
    const volScalarField::Internal& muc = tmuc();

    const dimensionedScalar dk(cbrt(pow3(dSphj) - pow3(dSphi)));

    const volScalarField::Internal tauCrit1
    (
        6*sigma/dSphj*(sqr(dSphi/dSphj) + sqr(dk/dSphj) - 1)
    );

    const volScalarField::Internal tauCrit2
    (
        sigma/min(dk, dSphi)
    );

    const volScalarField::Internal tauCrit(max(tauCrit1, tauCrit2));

    tmp<volScalarField::Internal> tbinaryBreakupRate =
        volScalarField::Internal::New
        (
            "binaryBreakupRate",
            popBal_.mesh(),
            dimensionedScalar(inv(dimVolume*dimTime), scalar(0))
        );
    volScalarField::Internal& binaryBreakupRate = tbinaryBreakupRate.ref();

    if (turbulence_)
    {
        tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
        const volScalarField::Internal& epsilonc = tepsilonc();

        const volScalarField::Internal tauTurb
        (
            pos(dSphj - kolmogorovLengthScale_)*BTurb_*rhoc
           *sqr(cbrt(epsilonc*dSphj))
        );

        binaryBreakupRate +=
            pos(tauTurb - tauCrit)
           /dSphj
           *sqrt(mag(tauTurb - tauCrit)/rhoc)
           /vj;
    }

    if (laminarShear_)
    {
        const volScalarField::Internal tauShear
        (
            BShear_*muc*shearStrainRate_
        );

        binaryBreakupRate +=
            pos(tauShear - tauCrit)
           /dSphj
           *sqrt(mag(tauShear - tauCrit)/rhoc)
           /vj;
    }

    if (turbulentShear_)
    {
        const volScalarField::Internal tauEddy
        (
            pos0(kolmogorovLengthScale_ - dSphj)
           *BEddy_
           *muc
           *eddyStrainRate_
        );

        binaryBreakupRate +=
            pos(tauEddy - tauCrit)
           /dSphj
           *sqrt(mag(tauEddy - tauCrit)/rhoc)/vj;
    }

    if (interfacialFriction_)
    {
        const volScalarField::Internal tauFric
        (
            BFric_*0.5*rhoc*sqr(uTerminal_[j])*Cd_[j]
        );

        binaryBreakupRate +=
            pos(tauFric - tauCrit)
           /dSphj
           *sqrt(mag(tauFric - tauCrit)/rhoc)/vj;
    }

    return tbinaryBreakupRate;
}


// ************************************************************************* //
