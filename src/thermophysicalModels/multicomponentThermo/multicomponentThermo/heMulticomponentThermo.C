/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "heMulticomponentThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class HeThermo>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::heMulticomponentThermo<HeThermo>::volScalarFieldPropertyi
(
    const word& psiName,
    const dimensionSet& psiDim,
    Method psiMethod,
    const label speciei,
    const Args& ... args
) const
{
    const typename HeThermo::mixtureType::thermoType& thermo =
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


template<class HeThermo>
template<class Method, class Arg, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::heMulticomponentThermo<HeThermo>::scalarFieldPropertyi
(
    Method psiMethod,
    const label speciei,
    const Arg& arg,
    const Args& ... args
) const
{
    const typename HeThermo::mixtureType::thermoType& thermo =
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

template<class HeThermo>
Foam::heMulticomponentThermo<HeThermo>::heMulticomponentThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    HeThermo(mesh, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class HeThermo>
Foam::heMulticomponentThermo<HeThermo>::~heMulticomponentThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class HeThermo>
Foam::scalar Foam::heMulticomponentThermo<HeThermo>::Wi
(
    const label speciei
) const
{
    return this->specieThermo(speciei).W();
}


template<class HeThermo>
Foam::scalar Foam::heMulticomponentThermo<HeThermo>::hfi
(
    const label speciei
) const
{
    return this->specieThermo(speciei).Hf();
}


template<class HeThermo>
Foam::scalar Foam::heMulticomponentThermo<HeThermo>::rhoi
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).rho(p, T);
}


template<class HeThermo>
Foam::tmp<Foam::volScalarField>
Foam::heMulticomponentThermo<HeThermo>::rhoi
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
        &HeThermo::mixtureType::thermoType::rho,
        speciei,
        p,
        T
    );
}


template<class HeThermo>
Foam::scalar Foam::heMulticomponentThermo<HeThermo>::Cpi
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).Cp(p, T);
}


template<class HeThermo>
Foam::tmp<Foam::volScalarField>
Foam::heMulticomponentThermo<HeThermo>::Cpi
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
        &HeThermo::mixtureType::thermoType::Cp,
        speciei,
        p,
        T
    );
}


template<class HeThermo>
Foam::scalar Foam::heMulticomponentThermo<HeThermo>::hei
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).HE(p, T);
}


template<class HeThermo>
Foam::tmp<Foam::scalarField> Foam::heMulticomponentThermo<HeThermo>::hei
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return scalarFieldPropertyi
    (
        &HeThermo::mixtureType::thermoType::HE,
        speciei,
        p,
        T
    );
}


template<class HeThermo>
Foam::tmp<Foam::volScalarField>
Foam::heMulticomponentThermo<HeThermo>::hei
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
        &HeThermo::mixtureType::thermoType::HE,
        speciei,
        p,
        T
    );
}


template<class HeThermo>
Foam::scalar Foam::heMulticomponentThermo<HeThermo>::hsi
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).Hs(p, T);
}


template<class HeThermo>
Foam::tmp<Foam::scalarField> Foam::heMulticomponentThermo<HeThermo>::hsi
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return scalarFieldPropertyi
    (
        &HeThermo::mixtureType::thermoType::Hs,
        speciei,
        p,
        T
    );
}


template<class HeThermo>
Foam::tmp<Foam::volScalarField>
Foam::heMulticomponentThermo<HeThermo>::hsi
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
        &HeThermo::mixtureType::thermoType::Hs,
        speciei,
        p,
        T
    );
}


template<class HeThermo>
Foam::scalar Foam::heMulticomponentThermo<HeThermo>::hai
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).Ha(p, T);
}


template<class HeThermo>
Foam::tmp<Foam::scalarField> Foam::heMulticomponentThermo<HeThermo>::hai
(
    const label speciei,
    const scalarField& p,
    const scalarField& T
) const
{
    return scalarFieldPropertyi
    (
        &HeThermo::mixtureType::thermoType::Ha,
        speciei,
        p,
        T
    );
}


template<class HeThermo>
Foam::tmp<Foam::volScalarField>
Foam::heMulticomponentThermo<HeThermo>::hai
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
        &HeThermo::mixtureType::thermoType::Ha,
        speciei,
        p,
        T
    );
}


template<class HeThermo>
Foam::scalar Foam::heMulticomponentThermo<HeThermo>::kappai
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    return this->specieThermo(speciei).kappa(p, T);
}


template<class HeThermo>
Foam::tmp<Foam::volScalarField>
Foam::heMulticomponentThermo<HeThermo>::kappai
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    return volScalarFieldPropertyi
    (
        "kappa",
        dimPower/dimLength/dimTemperature,
        &HeThermo::mixtureType::thermoType::kappa,
        speciei,
        p,
        T
    );
}


// ************************************************************************* //
