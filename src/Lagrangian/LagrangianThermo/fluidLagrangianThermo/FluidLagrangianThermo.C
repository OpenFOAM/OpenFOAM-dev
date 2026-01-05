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

#include "FluidLagrangianThermo.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::FluidLagrangianThermo<BaseThermo>::~FluidLagrangianThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
void Foam::FluidLagrangianThermo<BaseThermo>::correct
(
    const LagrangianSubMesh& subMesh
)
{
    if (BaseThermo::debug) InfoInFunction << endl;

    const SubField<scalar> e = subMesh.sub(this->e_.primitiveField());

    SubField<scalar> p = subMesh.sub(this->p_.primitiveFieldRef());
    SubField<scalar> T = subMesh.sub(this->T_.primitiveFieldRef());
    SubField<scalar> psi = subMesh.sub(this->psi_.primitiveFieldRef());
    SubField<scalar> rho = subMesh.sub(this->rho_.primitiveFieldRef());
    SubField<scalar> Cv = subMesh.sub(this->Cv_.primitiveFieldRef());
    SubField<scalar> kappa = subMesh.sub(this->kappa_.primitiveFieldRef());
    SubField<scalar> mu = subMesh.sub(this->mu_.primitiveFieldRef());

    auto Yslicer = this->Yslicer();

    forAll(T, subi)
    {
        const label i = subMesh.start() + subi;

        auto composition = this->elementComposition(Yslicer, i);

        const typename BaseThermo::mixtureType::thermoMixtureType&
            thermoMixture = this->thermoMixture(composition);

        const typename BaseThermo::mixtureType::transportMixtureType&
            transportMixture =
            this->transportMixture(composition, thermoMixture);

        T[subi] = thermoMixture.Tes(e[subi], p[subi], T[subi]);

        psi[subi] = thermoMixture.psi(p[subi], T[subi]);
        rho[subi] = thermoMixture.rho(p[subi], T[subi]);
        Cv[subi] = thermoMixture.Cv(p[subi], T[subi]);

        kappa[subi] = transportMixture.kappa(p[subi], T[subi]);
        mu[subi] = transportMixture.mu(p[subi], T[subi]);
    }

    if (BaseThermo::debug) Info<< "    Finished" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::FluidLagrangianThermo<BaseThermo>::psi
(
    const LagrangianSubScalarField& T,
    const LagrangianInjection& injection
) const
{
    return
        this->LagrangianInjectionProperty
        (
            injection,
            T.mesh(),
            "psi",
            dimDensity/dimPressure,
            &BaseThermo::mixtureType::thermoMixture,
            &BaseThermo::mixtureType::thermoMixtureType::psi,
            this->p(injection, T.mesh())(),
            T
        );
}


template<class BaseThermo>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::FluidLagrangianThermo<BaseThermo>::mu
(
    const LagrangianSubScalarField& T,
    const LagrangianInjection& injection
) const
{
    typedef decltype(this->Yslicer(injection, T.mesh())) compositionType;

    const typename BaseThermo::mixtureType::transportMixtureType&
        (BaseThermo::mixtureType::*mixture)(const compositionType&) const =
        &BaseThermo::mixtureType::transportMixture;

    return
        this->LagrangianInjectionProperty
        (
            injection,
            T.mesh(),
            "mu",
            dimDynamicViscosity,
            mixture,
            &BaseThermo::mixtureType::transportMixtureType::mu,
            this->p(injection, T.mesh())(),
            T
        );
}


// ************************************************************************* //
