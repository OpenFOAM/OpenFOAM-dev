/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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

#include "linearViscousStress.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
Foam::linearViscousStress<BasicMomentumTransportModel>::linearViscousStress
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
bool Foam::linearViscousStress<BasicMomentumTransportModel>::read()
{
    return BasicMomentumTransportModel::read();
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::surfaceVectorField>
Foam::linearViscousStress<BasicMomentumTransportModel>::devTau() const
{
    const surfaceScalarField alphaRhoNuEff
    (
        fvc::interpolate(this->alpha_*this->rho_*this->nuEff())
    );

    return surfaceVectorField::New
    (
        this->groupName("devTau"),
       -alphaRhoNuEff
       *(
           fvc::dotInterpolate(this->mesh().nf(), dev2(T(fvc::grad(this->U_))))
         + fvc::snGrad(this->U_)
        )
    );
}


template<class BasicMomentumTransportModel>
template<class RhoFieldType>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicMomentumTransportModel>::DivDevTau
(
    const RhoFieldType& rho,
    volVectorField& U
) const
{
    const surfaceScalarField alphaRhoNuEff
    (
        fvc::interpolate(this->alpha_*rho*this->nuEff())
    );

    const fvVectorMatrix divDevTauCorr
    (
        this->divDevTauCorr
        (
           -alphaRhoNuEff
           *fvc::dotInterpolate(this->mesh().Sf(), dev2(T(fvc::grad(U)))),
            U
        )
    );

    return divDevTauCorr - fvm::laplacian(alphaRhoNuEff, U);
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicMomentumTransportModel>::divDevTau
(
    volVectorField& U
) const
{
    return DivDevTau(this->rho_, U);
}


template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicMomentumTransportModel>::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return DivDevTau(rho, U);
}


template<class BasicMomentumTransportModel>
void Foam::linearViscousStress<BasicMomentumTransportModel>::correct()
{
    BasicMomentumTransportModel::correct();
}


// ************************************************************************* //
