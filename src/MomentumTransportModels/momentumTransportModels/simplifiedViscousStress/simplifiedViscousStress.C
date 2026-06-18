/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "simplifiedViscousStress.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::simplifiedViscousStress<BasicMomentumTransportModel>::
simplifiedViscousStress
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
    BasicMomentumTransportModel
    (
        modelName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool Foam::simplifiedViscousStress<BasicMomentumTransportModel>::read()
{
    return BasicMomentumTransportModel::read();
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::surfaceVectorField>
Foam::simplifiedViscousStress<BasicMomentumTransportModel>::devTau() const
{
    const surfaceScalarField alphaRhoNuEff
    (
        fvc::interpolate(this->alpha_*this->rho_*this->nuEff())
    );

    return surfaceVectorField::New
    (
        this->groupName("devTau"),
       -alphaRhoNuEff*fvc::snGrad(this->U_)
    );
}


template<class BasicMomentumTransportModel>
template<class RhoFieldType>
Foam::tmp<Foam::fvVectorMatrix>
Foam::simplifiedViscousStress<BasicMomentumTransportModel>::DivDevTau
(
    const RhoFieldType& rho,
    volVectorField& U
) const
{
    const surfaceScalarField alphaRhoNuEff
    (
        fvc::interpolate(this->alpha_*rho*this->nuEff())
    );

    return -fvm::laplacian(alphaRhoNuEff, U);
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::simplifiedViscousStress<BasicMomentumTransportModel>::divDevTau
(
    volVectorField& U
) const
{
    return DivDevTau(this->rho_, U);
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::simplifiedViscousStress<BasicMomentumTransportModel>::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return DivDevTau(rho, U);
}


template<class BasicMomentumTransportModel>
void Foam::simplifiedViscousStress<BasicMomentumTransportModel>::correct()
{
    BasicMomentumTransportModel::correct();
}


// ************************************************************************* //
