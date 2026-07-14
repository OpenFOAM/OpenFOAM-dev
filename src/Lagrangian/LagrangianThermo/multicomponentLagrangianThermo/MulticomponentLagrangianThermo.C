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

#include "MulticomponentLagrangianThermo.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BaseThermo>
template<class Method, class ... Args>
Foam::tmp<Foam::LagrangianInternalScalarField>
Foam::MulticomponentLagrangianThermo<BaseThermo>::
LagrangianInternalScalarFieldPropertyi
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

    forAll(psi, i)
    {
        psi[i] = (thermo.*psiMethod)(args[i] ...);
    }

    return tPsi;
}


template<class BaseThermo>
template<class Method, class ... Args>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::MulticomponentLagrangianThermo<BaseThermo>::
LagrangianSubScalarFieldPropertyi
(
    const LagrangianSubMesh& subMesh,
    const word& psiName,
    const dimensionSet& psiDim,
    Method psiMethod,
    const label speciei,
    const Args& ... args
) const
{
    const typename BaseThermo::mixtureType::thermoType& thermo =
        this->specieThermo(speciei);

    tmp<LagrangianSubScalarField> tPsi
    (
        LagrangianSubScalarField::New
        (
            IOobject::groupName(subMesh.sub(psiName), this->group()),
            subMesh,
            psiDim
        )
    );
    LagrangianSubScalarField& psi = tPsi.ref();

    forAll(psi, subi)
    {
        psi[subi] = (thermo.*psiMethod)(args[subi] ...);
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::MulticomponentLagrangianThermo<BaseThermo>::
~MulticomponentLagrangianThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::dimensionedScalar Foam::MulticomponentLagrangianThermo<BaseThermo>::Wi
(
    const label speciei
) const
{
    return
        dimensionedScalar
        (
            "W",
            dimensions::mass/dimensions::moles,
            this->specieThermo(speciei).W()
        );
}


template<class BaseThermo>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::MulticomponentLagrangianThermo<BaseThermo>::rhoi
(
    const label speciei,
    const LagrangianSubScalarField& p,
    const LagrangianSubScalarField& T
) const
{
    return LagrangianSubScalarFieldPropertyi
    (
        p.mesh(),
        "rho",
        dimensions::density,
        &BaseThermo::mixtureType::thermoType::rho,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::MulticomponentLagrangianThermo<BaseThermo>::hsi
(
    const label speciei,
    const LagrangianSubScalarField& p,
    const LagrangianSubScalarField& T
) const
{
    return LagrangianSubScalarFieldPropertyi
    (
        p.mesh(),
        "hs",
        dimensions::specificEnergy,
        &BaseThermo::mixtureType::thermoType::hs,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::dimensionedScalar Foam::MulticomponentLagrangianThermo<BaseThermo>::hfi
(
    const label speciei
) const
{
    return
        dimensionedScalar
        (
            "hf",
            dimensions::specificEnergy,
            this->specieThermo(speciei).hf()
        );
}


template<class BaseThermo>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::MulticomponentLagrangianThermo<BaseThermo>::Cvi
(
    const label speciei,
    const LagrangianSubScalarField& p,
    const LagrangianSubScalarField& T
) const
{
    return LagrangianSubScalarFieldPropertyi
    (
        p.mesh(),
        "Cv",
        dimensions::specificHeatCapacity,
        &BaseThermo::mixtureType::thermoType::Cv,
        speciei,
        p,
        T
    );
}


template<class BaseThermo>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::MulticomponentLagrangianThermo<BaseThermo>::Cpi
(
    const label speciei,
    const LagrangianSubScalarField& p,
    const LagrangianSubScalarField& T
) const
{
    return LagrangianSubScalarFieldPropertyi
    (
        p.mesh(),
        "Cp",
        dimensions::specificHeatCapacity,
        &BaseThermo::mixtureType::thermoType::Cp,
        speciei,
        p,
        T
    );
}


// ************************************************************************* //
