/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "MulticomponentThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BaseThermo>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::MulticomponentThermo<BaseThermo>::volScalarFieldPropertyi
(
    const word& psiName,
    const dimensionSet& psiDim,
    Method psiMethod,
    const label speciei,
    const Args& ... args
) const
{
    const typename BaseThermo::mixtureType::thermoType& thermo =
        this->specieThermo(speciei);

    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName(psiName, this->T_.group()),
            this->T_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(psi, celli)
    {
        psi[celli] = (thermo.*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        forAll(psiBf[patchi], patchFacei)
        {
            psiBf[patchi][patchFacei] =
                (thermo.*psiMethod)
                (
                    args.boundaryField()[patchi][patchFacei] ...
                );
        }
    }

    return tPsi;
}


template<class BaseThermo>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField::Internal>
Foam::MulticomponentThermo<BaseThermo>::volInternalScalarFieldPropertyi
(
    const word& psiName,
    const dimensionSet& psiDim,
    Method psiMethod,
    const label speciei,
    const Args& ... args
) const
{
    const typename BaseThermo::mixtureType::thermoType& thermo =
        this->specieThermo(speciei);

    tmp<volScalarField::Internal> tPsi
    (
        volScalarField::Internal::New
        (
            IOobject::groupName(psiName, this->T_.group()),
            this->T_.mesh(),
            psiDim
        )
    );

    volScalarField::Internal& psi = tPsi.ref();

    forAll(psi, celli)
    {
        psi[celli] = (thermo.*psiMethod)(args[celli] ...);
    }

    return tPsi;
}


template<class BaseThermo>
template<class Method, class Arg, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::MulticomponentThermo<BaseThermo>::scalarFieldPropertyi
(
    Method psiMethod,
    const label speciei,
    const Arg& arg,
    const Args& ... args
) const
{
    const typename BaseThermo::mixtureType::thermoType& thermo =
        this->specieThermo(speciei);

    tmp<scalarField> tPsi(new scalarField(arg.size()));

    scalarField& psi = tPsi.ref();

    forAll(psi, i)
    {
        psi[i] = (thermo.*psiMethod)(arg[i], args[i] ...);
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::MulticomponentThermo<BaseThermo>::MulticomponentThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    BaseThermo(mesh, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::MulticomponentThermo<BaseThermo>::~MulticomponentThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::scalar Foam::MulticomponentThermo<BaseThermo>::WiValue
(
    const label speciei
) const
{
    return this->specieThermo(speciei).W();
}


template<class BaseThermo>
Foam::dimensionedScalar Foam::MulticomponentThermo<BaseThermo>::Wi
(
    const label speciei
) const
{
    return
        dimensionedScalar
        (
            "W",
            dimMass/dimMoles,
            this->specieThermo(speciei).W()
        );
}


template<class BaseThermo>
Foam::scalar Foam::MulticomponentThermo<BaseThermo>::rhoi
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).rho(p, T);
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::MulticomponentThermo<BaseThermo>::rhoi
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldPropertyi
    (
        "rho",
        dimDensity,
        &BaseThermo::mixtureType::thermoType::rho,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::scalar Foam::MulticomponentThermo<BaseThermo>::Cpi
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).Cp(p, T);
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::MulticomponentThermo<BaseThermo>::Cpi
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldPropertyi
    (
        "Cp",
        dimEnergy/dimMass/dimTemperature,
        &BaseThermo::mixtureType::thermoType::Cp,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::scalar Foam::MulticomponentThermo<BaseThermo>::hei
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).he(p, T);
}


template<class BaseThermo>
Foam::tmp<Foam::scalarField> Foam::MulticomponentThermo<BaseThermo>::hei
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return scalarFieldPropertyi
    (
        &BaseThermo::mixtureType::thermoType::he,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField::Internal>
Foam::MulticomponentThermo<BaseThermo>::hei
(
    const label speciei,
    const volScalarField::Internal& p,
    const volScalarField::Internal& T
) const
{
    return volInternalScalarFieldPropertyi
    (
        "he",
        dimEnergy/dimMass,
        &BaseThermo::mixtureType::thermoType::he,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::MulticomponentThermo<BaseThermo>::hei
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldPropertyi
    (
        "he",
        dimEnergy/dimMass,
        &BaseThermo::mixtureType::thermoType::he,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::scalar Foam::MulticomponentThermo<BaseThermo>::hsi
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).hs(p, T);
}


template<class BaseThermo>
Foam::tmp<Foam::scalarField> Foam::MulticomponentThermo<BaseThermo>::hsi
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return scalarFieldPropertyi
    (
        &BaseThermo::mixtureType::thermoType::hs,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField::Internal>
Foam::MulticomponentThermo<BaseThermo>::hsi
(
    const label speciei,
    const volScalarField::Internal& p,
    const volScalarField::Internal& T
) const
{
    return volInternalScalarFieldPropertyi
    (
        "hs",
        dimEnergy/dimMass,
        &BaseThermo::mixtureType::thermoType::hs,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::MulticomponentThermo<BaseThermo>::hsi
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldPropertyi
    (
        "hs",
        dimEnergy/dimMass,
        &BaseThermo::mixtureType::thermoType::hs,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::scalar Foam::MulticomponentThermo<BaseThermo>::hai
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).ha(p, T);
}


template<class BaseThermo>
Foam::tmp<Foam::scalarField> Foam::MulticomponentThermo<BaseThermo>::hai
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return scalarFieldPropertyi
    (
        &BaseThermo::mixtureType::thermoType::ha,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField::Internal>
Foam::MulticomponentThermo<BaseThermo>::hai
(
    const label speciei,
    const volScalarField::Internal& p,
    const volScalarField::Internal& T
) const
{
    return volInternalScalarFieldPropertyi
    (
        "ha",
        dimEnergy/dimMass,
        &BaseThermo::mixtureType::thermoType::ha,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::MulticomponentThermo<BaseThermo>::hai
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldPropertyi
    (
        "ha",
        dimEnergy/dimMass,
        &BaseThermo::mixtureType::thermoType::ha,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::dimensionedScalar Foam::MulticomponentThermo<BaseThermo>::hfi
(
    const label speciei
) const
{
    return
        dimensionedScalar
        (
            "hf",
            dimEnergy/dimMass,
            this->specieThermo(speciei).hf()
        );
}


template<class BaseThermo>
Foam::scalar Foam::MulticomponentThermo<BaseThermo>::hfiValue
(
    const label speciei
) const
{
    return this->specieThermo(speciei).hf();
}


template<class BaseThermo>
Foam::scalar Foam::MulticomponentThermo<BaseThermo>::kappai
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).kappa(p, T);
}


template<class BaseThermo>
Foam::tmp<Foam::volScalarField>
Foam::MulticomponentThermo<BaseThermo>::kappai
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldPropertyi
    (
        "kappa",
        dimThermalConductivity,
        &BaseThermo::mixtureType::thermoType::kappa,
        speciei,
        p,
        T
    );
}


// ************************************************************************* //
