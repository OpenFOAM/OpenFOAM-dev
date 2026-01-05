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

#include "BasicLagrangianThermo.H"
#include "energyLagrangianScalarFieldSource.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::LagrangianInternalScalarField>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::
LagrangianInternalScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    Mixture mixture,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<LagrangianInternalScalarField> tPsi
    (
        LagrangianInternalScalarField::New
        (
            IOobject::groupName(psiName, this->group()),
            this->mesh(),
            psiDim
        )
    );
    LagrangianInternalScalarField& psi = tPsi.ref();

    auto Yslicer = this->Yslicer();

    forAll(psi, i)
    {
        auto composition = this->elementComposition(Yslicer, i);

        psi[i] = ((this->*mixture)(composition).*psiMethod)(args[i] ...);
    }

    return tPsi;
}


template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::
LagrangianSubScalarFieldProperty
(
    const LagrangianSubMesh& subMesh,
    const word& psiName,
    const dimensionSet& psiDim,
    Mixture mixture,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<LagrangianSubScalarField> tPsi
    (
        LagrangianSubScalarField::New
        (
            IOobject::groupName(psiName, this->group()),
            subMesh,
            psiDim
        )
    );
    LagrangianSubScalarField& psi = tPsi.ref();

    auto Yslicer = this->Yslicer();

    forAll(psi, subi)
    {
        const label i = subMesh.start() + subi;

        auto composition = this->elementComposition(Yslicer, i);

        psi[subi] = ((this->*mixture)(composition).*psiMethod)(args[i] ...);
    }

    return tPsi;
}


template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::
LagrangianInjectionProperty
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh,
    const word& psiName,
    const dimensionSet& psiDim,
    Mixture mixture,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<LagrangianSubScalarField> tPsi
    (
        LagrangianSubScalarField::New
        (
            IOobject::groupName(psiName, this->group()),
            subMesh,
            psiDim
        )
    );
    LagrangianSubScalarField& psi = tPsi.ref();

    auto Yslicer = this->Yslicer(injection, subMesh);

    forAll(psi, subi)
    {
        auto composition = this->injectionElementComposition(Yslicer, subi);

        psi[subi] = ((this->*mixture)(composition).*psiMethod)(args[subi] ...);
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::BasicLagrangianThermo
(
    const LagrangianMesh& mesh,
    const word& phaseName
)
:
    physicalProperties(mesh, phaseName),
    MixtureType(properties()),
    BasicThermoType
    (
        properties(),
        static_cast<const MixtureType&>(*this),
        mesh,
        phaseName
    ),
    e_
    (
        IOobject
        (
            IOobject::groupName
            (
                MixtureType::thermoType::heName(),
                phaseName
            ),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        LagrangianInternalScalarFieldProperty
        (
            "e",
            dimEnergy/dimMass,
            &MixtureType::thermoMixture,
            &MixtureType::thermoMixtureType::es,
            this->p_,
            this->T_
        )(),
        this->eBoundaryTypes(),
        this->eBoundaryBaseTypes(),
        this->template sourcesTypes<energyLagrangianScalarFieldSource>
        (
            this->T_
        ),
        this->T_.sources().errorLocation()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::
~BasicLagrangianThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
const Foam::IOdictionary&
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::properties() const
{
    return *this;
}


template<class MixtureType, class BasicThermoType>
Foam::IOdictionary&
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::properties()
{
    return *this;
}


template<class MixtureType, class BasicThermoType>
Foam::word
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::thermoName() const
{
    return MixtureType::thermoType::typeName();
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::W
(
    const LagrangianSubMesh& subMesh
) const
{
    return
        LagrangianSubScalarFieldProperty
        (
            subMesh,
            "W",
            dimMass/dimMoles,
            &MixtureType::thermoMixture,
            &MixtureType::thermoMixtureType::W
        );
}


template<class MixtureType, class BasicThermoType>
const Foam::LagrangianScalarDynamicField&
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::e() const
{
    return e_;
}


template<class MixtureType, class BasicThermoType>
Foam::LagrangianScalarDynamicField&
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::e()
{
    return e_;
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::rho
(
    const LagrangianSubScalarField& T,
    const LagrangianInjection& injection
) const
{
    return
        LagrangianInjectionProperty
        (
            injection,
            T.mesh(),
            "rho",
            dimDensity,
            &MixtureType::thermoMixture,
            &MixtureType::thermoMixtureType::rho,
            this->p(injection, T.mesh())(),
            T
        );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::e
(
    const LagrangianSubScalarField& T,
    const LagrangianInjection& injection
) const
{
    return
        LagrangianInjectionProperty
        (
            injection,
            T.mesh(),
            "e",
            dimEnergy/dimMass,
            &MixtureType::thermoMixture,
            &MixtureType::thermoMixtureType::es,
            this->p(injection, T.mesh())(),
            T
        );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::Cv
(
    const LagrangianSubScalarField& T,
    const LagrangianInjection& injection
) const
{
    return
        LagrangianInjectionProperty
        (
            injection,
            T.mesh(),
            "Cv",
            dimSpecificHeatCapacity,
            &MixtureType::thermoMixture,
            &MixtureType::thermoMixtureType::Cv,
            this->p(injection, T.mesh())(),
            T
        );
}


template<class MixtureType, class BasicThermoType>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::kappa
(
    const LagrangianSubScalarField& T,
    const LagrangianInjection& injection
) const
{
    typedef decltype(this->Yslicer(injection, T.mesh())) compositionType;

    const typename MixtureType::transportMixtureType&
        (MixtureType::*mixture)(const compositionType&) const =
        &MixtureType::transportMixture;

    return
        LagrangianInjectionProperty
        (
            injection,
            T.mesh(),
            "kappa",
            dimThermalConductivity,
            mixture,
            &MixtureType::transportMixtureType::kappa,
            this->p(injection, T.mesh())(),
            T
        );
}


template<class MixtureType, class BasicThermoType>
bool Foam::BasicLagrangianThermo<MixtureType, BasicThermoType>::read()
{
    if (physicalProperties::read())
    {
        MixtureType::read(*this);
        BasicThermoType::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
