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

#include "LehrMilliesMewesCoalescence.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(LehrMilliesMewesCoalescence, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        LehrMilliesMewesCoalescence,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::coalescenceModels::LehrMilliesMewesCoalescence::
LehrMilliesMewesCoalescence
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    uCrit_("uCrit", dimVelocity, dict, 0.08),
    alphaMax_("alphaMax", dimless, dict, 0.6)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalance::coalescenceModels::LehrMilliesMewesCoalescence::rate
(
    const label i,
    const label j
) const
{
    using Foam::constant::mathematical::pi;

    const dimensionedScalar& dSphi = popBal_.dSph(i);
    const dimensionedScalar& dSphj = popBal_.dSph(j);
    const phaseModel& phasei = popBal_.phases()[i];
    const phaseModel& phasej = popBal_.phases()[j];

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();

    const volScalarField::Internal uChar
    (
        max
        (
            sqrt(2.0)
           *cbrt(epsilonc)
           *sqrt(cbrt(sqr(dSphi)) + cbrt(sqr(dSphj))),
            mag(phasei.U()()() - phasej.U()()())
        )
    );

    return
        pi/4
       *sqr(dSphi + dSphj)
       *min(uChar, uCrit_)
       *exp
        (
          - sqr(cbrt(alphaMax_)
           /cbrt(max(popBal_.alphas()(), phasei.residualAlpha())) - 1)
        );
}


// ************************************************************************* //
