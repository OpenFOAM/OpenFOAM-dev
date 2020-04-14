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

#include "PhaseIncompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>::
PhaseIncompressibleMomentumTransportModel
(
    const word& type,
    const volScalarField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel
)
:
    MomentumTransportModel
    <
        volScalarField,
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
        transportModel
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>>
Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>::New
(
    const volScalarField& alpha,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel
)
{
    return autoPtr<PhaseIncompressibleMomentumTransportModel>
    (
        static_cast<PhaseIncompressibleMomentumTransportModel*>(
        MomentumTransportModel
        <
            volScalarField,
            geometricOneField,
            incompressibleMomentumTransportModel,
            TransportModel
        >::New
        (
            alpha,
            geometricOneField(),
            U,
            alphaRhoPhi,
            phi,
            transportModel
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volScalarField>
Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>::pPrime() const
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
Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>::pPrimef() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("pPrimef", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimPressure, 0)
    );
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>::
devSigma() const
{
    return devTau();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>::divDevSigma
(
    volVectorField& U
) const
{
    return divDevTau(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>::
devTau() const
{
    NotImplemented;

    return devSigma();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleMomentumTransportModel<TransportModel>::
divDevTau
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevSigma(U);
}


// ************************************************************************* //
