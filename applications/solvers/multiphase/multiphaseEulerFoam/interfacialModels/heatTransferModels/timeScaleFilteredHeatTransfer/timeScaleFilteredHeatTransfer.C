/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "timeScaleFilteredHeatTransfer.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(timeScaleFilteredHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        heatTransferModel,
        timeScaleFilteredHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::timeScaleFilteredHeatTransfer::
timeScaleFilteredHeatTransfer
(
    const dictionary& dict,
    const phasePair& pair
)
:
    heatTransferModel(dict.subDict("heatTransferModel"), pair),
    heatTransferModel_
    (
        heatTransferModel::New(dict.subDict("heatTransferModel"), pair)
    ),
    minRelaxTime_("minRelaxTime", dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::timeScaleFilteredHeatTransfer::
~timeScaleFilteredHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::timeScaleFilteredHeatTransfer::K
(
    const scalar residualAlpha
) const
{
    const volScalarField limit
    (
        max(pair_.dispersed(), residualAlpha)
       *pair_.dispersed().thermo().Cp()
       *pair_.dispersed().rho()
       /minRelaxTime_
    );

    return min(heatTransferModel_->K(residualAlpha), limit);
}


// ************************************************************************* //
