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

#include "PhaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::PhaseCompressibleMomentumTransportModel<TransportModel>::
PhaseCompressibleMomentumTransportModel
(
    const word& type,
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
:
    MomentumTransportModel
    <
        volScalarField,
        volScalarField,
        compressibleMomentumTransportModel,
        transportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::PhaseCompressibleMomentumTransportModel<TransportModel>>
Foam::PhaseCompressibleMomentumTransportModel<TransportModel>::New
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
{
    return autoPtr<PhaseCompressibleMomentumTransportModel>
    (
        static_cast<PhaseCompressibleMomentumTransportModel*>(
        MomentumTransportModel
        <
            volScalarField,
            volScalarField,
            compressibleMomentumTransportModel,
            transportModel
        >::New
        (
            alpha,
            rho,
            U,
            alphaRhoPhi,
            phi,
            transport
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volScalarField>
Foam::PhaseCompressibleMomentumTransportModel<TransportModel>::pPrime() const
{
    return volScalarField::New
    (
        IOobject::groupName("pPrime", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimPressure, 0)
    );
}


template<class TransportModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PhaseCompressibleMomentumTransportModel<TransportModel>::pPrimef() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("pPrimef", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimPressure, 0)
    );
}


// ************************************************************************* //
