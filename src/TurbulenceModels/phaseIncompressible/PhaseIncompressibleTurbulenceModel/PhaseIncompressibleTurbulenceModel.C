/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "PhaseIncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::
PhaseIncompressibleTurbulenceModel
(
    const volScalarField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel,
    const word& propertiesName
)
:
    TurbulenceModel
    <
        volScalarField,
        geometricOneField,
        incompressibleTurbulenceModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transportModel,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::PhaseIncompressibleTurbulenceModel<TransportModel> >
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::New
(
    const volScalarField& alpha,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel,
    const word& propertiesName
)
{
    return autoPtr<PhaseIncompressibleTurbulenceModel>
    (
        static_cast<PhaseIncompressibleTurbulenceModel*>(
        TurbulenceModel
        <
            volScalarField,
            geometricOneField,
            incompressibleTurbulenceModel,
            TransportModel
        >::New
        (
            alpha,
            geometricOneField(),
            U,
            alphaRhoPhi,
            phi,
            transportModel,
            propertiesName
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volScalarField>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::pPrime() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("pPrime", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("pPrimef", dimPressure, 0.0)
        )
    );
}


template<class TransportModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::pPrimef() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("pPrimef", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("pPrimef", dimPressure, 0.0)
        )
    );
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::devReff() const
{
    return devRhoReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U
) const
{
    return divDevRhoReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::
devRhoReff() const
{
    notImplemented
    (
        "PhaseIncompressibleTurbulenceModel<TransportModel>::"
        "devRhoReff()"
    );

    return devReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleTurbulenceModel<TransportModel>::
divDevRhoReff
(
    volVectorField& U
) const
{
    notImplemented
    (
        "PhaseIncompressibleTurbulenceModel<TransportModel>::"
        "divDevRhoReff(volVectorField& U)"
    );

    return divDevReff(U);
}


// ************************************************************************* //
