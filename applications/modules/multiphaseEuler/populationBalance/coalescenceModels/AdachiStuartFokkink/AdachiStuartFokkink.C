/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

#include "AdachiStuartFokkink.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace coalescenceModels
{
    defineTypeNameAndDebug(AdachiStuartFokkink, 0);
    addToRunTimeSelectionTable
    (
        coalescenceModel,
        AdachiStuartFokkink,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::coalescenceModels::AdachiStuartFokkink::
AdachiStuartFokkink
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    coalescenceModel(popBal, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volInternalScalarField>
Foam::populationBalance::coalescenceModels::AdachiStuartFokkink::rate
(
    const label i,
    const label j
) const
{
    using Foam::constant::mathematical::pi;

    tmp<volScalarField> tdi = popBal_.d(i);
    const volInternalScalarField& di = tdi();
    tmp<volScalarField> tdj = popBal_.d(j);
    const volInternalScalarField& dj = tdj();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volInternalScalarField& epsilonc = tepsilonc();
    tmp<volScalarField> tnuc(popBal_.continuousPhase().fluidThermo().nu());
    const volInternalScalarField nuc = tnuc();

    return (4.0/3.0)*sqrt(0.3*pi*epsilonc/nuc)*pow3(di + dj);
}


// ************************************************************************* //
