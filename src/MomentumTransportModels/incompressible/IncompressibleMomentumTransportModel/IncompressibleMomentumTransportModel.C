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

#include "IncompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::IncompressibleMomentumTransportModel<TransportModel>::
IncompressibleMomentumTransportModel
(
    const word& type,
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport
)
:
    MomentumTransportModel
    <
        geometricOneField,
        geometricOneField,
        incompressibleMomentumTransportModel,
        TransportModel
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
Foam::autoPtr<Foam::IncompressibleMomentumTransportModel<TransportModel>>
Foam::IncompressibleMomentumTransportModel<TransportModel>::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const TransportModel& transport
)
{
    return autoPtr<IncompressibleMomentumTransportModel>
    (
        static_cast<IncompressibleMomentumTransportModel*>(
        MomentumTransportModel
        <
            geometricOneField,
            geometricOneField,
            incompressibleMomentumTransportModel,
            TransportModel
        >::New
        (
            geometricOneField(),
            geometricOneField(),
            U,
            phi,
            phi,
            transport
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleMomentumTransportModel<TransportModel>::devSigma() const
{
    return devTau();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleMomentumTransportModel<TransportModel>::divDevSigma
(
    volVectorField& U
) const
{
    return divDevTau(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleMomentumTransportModel<TransportModel>::
devTau() const
{
    NotImplemented;

    return devSigma();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleMomentumTransportModel<TransportModel>::
divDevTau
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevSigma(U);
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleMomentumTransportModel<TransportModel>::
divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;

    return divDevSigma(U);
}


// ************************************************************************* //
