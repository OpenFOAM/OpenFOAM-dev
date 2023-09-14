/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "diffusiveMassTransferModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusiveMassTransferModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug
    (
        diffusiveMassTransferModel,
        0
    );
    defineSidedInterfacialModelTypeNameAndDebug
    (
        blendedDiffusiveMassTransferModel,
        0
    );
    defineRunTimeSelectionTable(diffusiveMassTransferModel, dictionary);
}

const Foam::dimensionSet Foam::diffusiveMassTransferModel::dimK(0, -2, 0, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModel::diffusiveMassTransferModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusiveMassTransferModel::~diffusiveMassTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::blendedDiffusiveMassTransferModel::K() const
{
    tmp<volScalarField> (diffusiveMassTransferModel::*k)() const =
        &diffusiveMassTransferModel::K;
    return evaluate(k, "K", diffusiveMassTransferModel::dimK, false);
}


// ************************************************************************* //
