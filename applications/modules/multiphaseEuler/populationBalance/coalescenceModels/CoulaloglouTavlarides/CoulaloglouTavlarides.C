/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

#include "CoulaloglouTavlarides.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(CoulaloglouTavlarides, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        CoulaloglouTavlarides,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::CoulaloglouTavlarides::
CoulaloglouTavlarides
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    C1_("C1", dimless, dict, 2.8),
    C2_("C2", inv(dimArea), dict, 1.83e9)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceModels::CoulaloglouTavlarides::
addToCoalescenceRate
(
    volScalarField::Internal& coalescenceRate,
    const label i,
    const label j
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(fi.phase()));
    const volScalarField::Internal& sigma = tsigma();

    tmp<volScalarField> tmuc(popBal_.continuousPhase().fluidThermo().mu());
    const volScalarField::Internal& muc = tmuc();
    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();

    coalescenceRate +=
        C1_*(pow(fi.x(), 2.0/3.0) + pow(fj.x(), 2.0/3.0))
       *sqrt(pow(fi.x(), 2.0/9.0) + pow(fj.x(), 2.0/9.0))
       *cbrt(epsilonc)
       /(1 + popBal_.alphas()())
       *exp
        (
          - C2_*muc*rhoc*epsilonc/sqr(sigma)
           /pow3(1 + popBal_.alphas()())
           *pow4(cbrt(fi.x())*cbrt(fj.x())/(cbrt(fi.x()) + cbrt(fj.x())))
        );
}


// ************************************************************************* //
