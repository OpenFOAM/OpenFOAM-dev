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

#include "MulticomponentCarrierThermo.H"
#include "dimensions.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BaseThermo>
template<class Method, class ... Args>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::MulticomponentCarrierThermo<BaseThermo>::LagrangianSubScalarFieldPropertyi
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
            IOobject::groupName(subMesh.sub(psiName), this->phaseName()),
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::MulticomponentCarrierThermo<BaseThermo>::MulticomponentCarrierThermo
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const basicThermo& thermo
)
:
    BaseThermo(c, carriedCloud, thermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::MulticomponentCarrierThermo<BaseThermo>::~MulticomponentCarrierThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BaseThermo>
Foam::dimensionedScalar Foam::MulticomponentCarrierThermo<BaseThermo>::Wi
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
Foam::MulticomponentCarrierThermo<BaseThermo>::rhoi
(
    const label speciei,
    const LagrangianSubMesh& subMesh
) const
{
    return LagrangianSubScalarFieldPropertyi
    (
        subMesh,
        "rho",
        dimensions::density,
        &BaseThermo::mixtureType::thermoType::rho,
        speciei,
        this->p(subMesh)(),
        this->T(subMesh)()
    );
}


template<class BaseThermo>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::MulticomponentCarrierThermo<BaseThermo>::hsi
(
    const label speciei,
    const LagrangianSubMesh& subMesh
) const
{
    return LagrangianSubScalarFieldPropertyi
    (
        subMesh,
        "hsi",
        dimensions::specificEnergy,
        &BaseThermo::mixtureType::thermoType::hs,
        speciei,
        this->p(subMesh)(),
        this->T(subMesh)()
    );
}


template<class BaseThermo>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::MulticomponentCarrierThermo<BaseThermo>::hai
(
    const label speciei,
    const LagrangianSubMesh& subMesh
) const
{
    return LagrangianSubScalarFieldPropertyi
    (
        subMesh,
        "ha",
        dimensions::specificEnergy,
        &BaseThermo::mixtureType::thermoType::ha,
        speciei,
        this->p(subMesh)(),
        this->T(subMesh)()
    );
}


template<class BaseThermo>
Foam::dimensionedScalar Foam::MulticomponentCarrierThermo<BaseThermo>::hfi
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


// ************************************************************************* //
