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

#include "turbulentShear.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(turbulentShear, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        turbulentShear,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::turbulentShear::
turbulentShear
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict),
    C_("C", dimless, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceModels::turbulentShear::
addToCoalescenceRate
(
    volScalarField::Internal& coalescenceRate,
    const label i,
    const label j
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];
    const sizeGroup& fj = popBal_.sizeGroups()[j];

    tmp<volScalarField> tdi = fi.d();
    const volScalarField::Internal& di = tdi();
    tmp<volScalarField> tdj = fj.d();
    const volScalarField::Internal& dj = tdj();

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();
    tmp<volScalarField> tmu(popBal_.continuousPhase().fluidThermo().mu());
    const volScalarField::Internal muc = tmu();

    coalescenceRate += C_*sqrt(epsilonc*rhoc/muc)*pow3(di + dj);
}


// ************************************************************************* //
