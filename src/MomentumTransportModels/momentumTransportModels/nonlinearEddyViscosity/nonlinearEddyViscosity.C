/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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

#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::nonlinearEddyViscosity<BasicMomentumTransportModel>::
nonlinearEddyViscosity
(
    const word& modelName,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    eddyViscosity<BasicMomentumTransportModel>
    (
        modelName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    nonlinearStress_
    (
        IOobject
        (
            this->groupName("nonlinearStress"),
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedSymmTensor
        (
            "nonlinearStress",
            sqr(dimVelocity),
            Zero
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::nonlinearEddyViscosity<BasicMomentumTransportModel>::sigma() const
{
    tmp<volSymmTensorField> tR
    (
        eddyViscosity<BasicMomentumTransportModel>::sigma()
    );
    tR.ref() += nonlinearStress_;
    return tR;
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::surfaceVectorField>
Foam::nonlinearEddyViscosity<BasicMomentumTransportModel>::devTau() const
{
    tmp<surfaceVectorField> tdevTau
    (
        eddyViscosity<BasicMomentumTransportModel>::devTau()
    );

    tdevTau.ref() += fvc::dotInterpolate
    (
        this->mesh().nf(),
        this->rho_*nonlinearStress_
    );

    return tdevTau;
}


template<class BasicMomentumTransportModel>
template<class RhoFieldType>
Foam::tmp<Foam::fvVectorMatrix>
Foam::nonlinearEddyViscosity<BasicMomentumTransportModel>::DivDevTau
(
    const RhoFieldType& rho,
    volVectorField& U
) const
{
    return
    (
        fvc::div(rho*nonlinearStress_)
      + eddyViscosity<BasicMomentumTransportModel>::DivDevTau(rho, U)
    );
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::nonlinearEddyViscosity<BasicMomentumTransportModel>::divDevTau
(
    volVectorField& U
) const
{
    return DivDevTau(this->rho_, U);
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::nonlinearEddyViscosity<BasicMomentumTransportModel>::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return DivDevTau(rho, U);
}


// ************************************************************************* //
