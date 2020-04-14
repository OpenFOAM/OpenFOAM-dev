/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "MomentumTransportModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class Alpha,
    class Rho,
    class BasicMomentumTransportModel,
    class TransportModel
>
Foam::MomentumTransportModel
<
    Alpha,
    Rho,
    BasicMomentumTransportModel,
    TransportModel
>::MomentumTransportModel
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
:
    BasicMomentumTransportModel
    (
        rho,
        U,
        alphaRhoPhi,
        phi
    ),
    alpha_(alpha),
    transport_(transport)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template
<
    class Alpha,
    class Rho,
    class BasicMomentumTransportModel,
    class TransportModel
>
Foam::autoPtr
<
    Foam::MomentumTransportModel
    <
        Alpha,
        Rho,
        BasicMomentumTransportModel,
        TransportModel
    >
>
Foam::MomentumTransportModel
<
    Alpha,
    Rho,
    BasicMomentumTransportModel,
    TransportModel
>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
{
    const word modelType
    (
        momentumTransportModel::readModelDict
        (
            U.db(),
            alphaRhoPhi.group()
        ).lookup("simulationType")
    );

    Info<< "Selecting turbulence model type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown MomentumTransportModel type "
            << modelType << nl << nl
            << "Valid MomentumTransportModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<MomentumTransportModel>
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, transport)
    );
}


// ************************************************************************* //
