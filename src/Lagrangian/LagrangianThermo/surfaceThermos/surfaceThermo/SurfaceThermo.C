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

#include "SurfaceThermo.H"
#include "LagrangianSubFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
template<class Mixture, class Method, class ... Args>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::SurfaceThermo<MixtureType, BasicThermoType>::
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
            IOobject::groupName(subMesh.sub(psiName), this->phaseName()),
            subMesh,
            psiDim
        )
    );
    LagrangianSubScalarField& psi = tPsi.ref();

    auto Yslicer = this->Yslicer(subMesh);

    forAll(psi, subi)
    {
        const label i = subMesh.start() + subi;

        auto composition = this->elementComposition(Yslicer, i);

        psi[subi] = ((this->*mixture)(composition).*psiMethod)(args[i] ...);
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
Foam::SurfaceThermo<MixtureType, BasicThermoType>::SurfaceThermo
(
    const basicLagrangianSubThermo& farThermo
)
:
    MixtureType(farThermo.properties()),
    BasicThermoType(farThermo)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
Foam::SurfaceThermo<MixtureType, BasicThermoType>::~SurfaceThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class MixtureType, class BasicThermoType>
Foam::word
Foam::SurfaceThermo<MixtureType, BasicThermoType>::mixtureName() const
{
    return MixtureType::typeName();
}


// ************************************************************************* //
