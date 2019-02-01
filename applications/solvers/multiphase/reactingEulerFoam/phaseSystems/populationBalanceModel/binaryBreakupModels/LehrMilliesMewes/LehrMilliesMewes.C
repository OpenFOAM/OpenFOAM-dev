/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
#include "addToRunTimeSelectionTable.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
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

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::binaryBreakupModels::LehrMilliesMewes::LehrMilliesMewes
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binaryBreakupModel(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::binaryBreakupModels::LehrMilliesMewes::
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

    volScalarField L
    (
        pow
        (
            popBal_.sigmaWithContinuousPhase(fj.phase())/continuousPhase.rho(),
            3.0/5.0
        )
       /pow(popBal_.continuousTurbulence().epsilon(), 2.0/5.0)
    );

    // Reset of dimension to pure length to avoid problems in transcendental
    // functions due to small exponents
    L.dimensions().reset(dimLength);

    const volScalarField T
    (
        pow
        (
            popBal_.sigmaWithContinuousPhase(fj.phase())/continuousPhase.rho(),
            2.0/5.0
        )
       /pow(popBal_.continuousTurbulence().epsilon(), 3.0/5.0)
    );

    binaryBreakupRate +=
        0.5*pow(fj.d()/L, 5.0/3.0)
       *exp(-sqrt(2.0)/pow3(fj.d()/L))
       *6.0/pow(pi, 1.5)/pow3(fi.d()/L)
       *exp(-9.0/4.0*sqr(log(pow(2.0, 0.4)*fi.d()/L)))
       /max(1.0 + erf(1.5*log(pow(2.0, 1.0/15.0)*fj.d()/L)), SMALL)
       /(T*pow3(L));
}


// ************************************************************************* //
