/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "sphericalDiffusiveMassTransfer.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diffusiveMassTransferModels
{
    defineTypeNameAndDebug(sphericalDiffusiveMassTransfer, 0);
    addToRunTimeSelectionTable
    (
        diffusiveMassTransferModel,
        sphericalDiffusiveMassTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModels::sphericalDiffusiveMassTransfer::
sphericalDiffusiveMassTransfer
(
    const dictionary& dict,
    const phasePair& pair
)
:
    diffusiveMassTransferModel(dict, pair),
    Le_("Le", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModels::sphericalDiffusiveMassTransfer::
~sphericalDiffusiveMassTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diffusiveMassTransferModels::sphericalDiffusiveMassTransfer::K() const
{
    return 60*pair_.dispersed()/sqr(pair_.dispersed().d());
}


// ************************************************************************* //
