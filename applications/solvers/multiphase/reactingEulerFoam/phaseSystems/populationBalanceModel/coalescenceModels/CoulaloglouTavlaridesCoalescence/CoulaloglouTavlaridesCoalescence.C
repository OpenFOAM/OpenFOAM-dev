/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "CoulaloglouTavlaridesCoalescence.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(CoulaloglouTavlaridesCoalescence, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        CoulaloglouTavlaridesCoalescence,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::CoulaloglouTavlaridesCoalescence::
CoulaloglouTavlaridesCoalescence
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    C1_(dimensionedScalar::lookupOrDefault("C1", dict, dimless, 2.8)),
    C2_(dimensionedScalar::lookupOrDefault("C2", dict, inv(dimArea), 1.83e9))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::diameterModels::coalescenceModels::CoulaloglouTavlaridesCoalescence::
addToCoalescenceRate
(
    volScalarField& coalescenceRate,
    const label i,
    const label j
)
{
    const phaseModel& continuousPhase = popBal_.continuousPhase();
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    coalescenceRate +=
        C1_*(pow(fi.x(), 2.0/3.0) + pow(fj.x(), 2.0/3.0))
       *sqrt(pow(fi.x(), 2.0/9.0) + pow(fj.x(), 2.0/9.0))
       *cbrt(popBal_.continuousTurbulence().epsilon())/(1 + popBal_.alphas())
       *exp
        (
          - C2_*continuousPhase.mu()*continuousPhase.rho()
           *popBal_.continuousTurbulence().epsilon()
           /sqr(popBal_.sigmaWithContinuousPhase(fi.phase()))
           /pow3(1 + popBal_.alphas())
           *pow4(cbrt(fi.x())*cbrt(fj.x())/(cbrt(fi.x()) + cbrt(fj.x())))
        );
}


// ************************************************************************* //
