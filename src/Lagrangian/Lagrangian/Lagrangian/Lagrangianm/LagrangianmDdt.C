/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianmDdt.H"
#include "LagrangianDdtScheme.H"
#include "LagrangianMesh.H"
#include "LagrangianSubFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::Lagrangianm::initDdt
(
    const dimensionSet& mDims,
    const LagrangianSubSubField<Type>& psi,
    const bool instantaneousDdt
)
{
    const LagrangianMesh& mesh = psi.mesh().mesh();

    return
        Lagrangian::ddtScheme<Type>::New
        (
            mesh,
            mesh.schemes().ddt("ddt(" + psi.name() + ')')
        ).ref().LagrangianmInitDdt(mDims, psi, instantaneousDdt);
}


template<class Type, template<class> class PrimitiveField>
bool Foam::Lagrangianm::initDdt
(
    const dimensionSet& mDims,
    const LagrangianSubField<Type, PrimitiveField>& psi,
    const bool instantaneousDdt
)
{
    return Lagrangianm::initDdt(mDims, toSubField(psi)(), instantaneousDdt);
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::noDdt
(
    const LagrangianSubScalarField& deltaT,
    const dimensionSet& mDims,
    const LagrangianSubSubField<Type>& psi
)
{
    const LagrangianMesh& mesh = psi.mesh().mesh();

    return
        Lagrangian::ddtScheme<Type>::New
        (
            mesh,
            mesh.schemes().ddt("ddt(" + psi.name() + ')')
        ).ref().LagrangianmNoDdt(deltaT, mDims, psi);
}


template<class Type, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::noDdt
(
    const LagrangianSubScalarField& deltaT,
    const dimensionSet& mDims,
    const LagrangianSubField<Type, PrimitiveField>& psi
)
{
    return Lagrangianm::noDdt(deltaT, mDims, toSubField(psi)());
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::Ddt
(
    const LagrangianSubScalarField& deltaT,
    LagrangianSubSubField<Type>& psi
)
{
    const LagrangianMesh& mesh = psi.mesh().mesh();

    return
        Lagrangian::ddtScheme<Type>::New
        (
            mesh,
            mesh.schemes().ddt("ddt(" + psi.name() + ')')
        ).ref().LagrangianmDdt(deltaT, psi);
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::Ddt
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& m,
    LagrangianSubSubField<Type>& psi
)
{
    const LagrangianMesh& mesh = psi.mesh().mesh();

    return
        Lagrangian::ddtScheme<Type>::New
        (
            mesh,
            mesh.schemes().ddt("ddt(" + m.name() + ',' + psi.name() + ')')
        ).ref().LagrangianmDdt(deltaT, m, psi);
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::ddt
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubSubField<Type>& psi
)
{
    return Lagrangian::ddtScheme<Type>::Lagrangianmddt(deltaT, psi);
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::ddt
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& psi
)
{
    return Lagrangian::ddtScheme<Type>::Lagrangianmddt(deltaT, m, psi);
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::ddt0
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubSubField<Type>& psi
)
{
    return Lagrangian::ddtScheme<Type>::Lagrangianmddt0(deltaT, psi);
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::ddt0
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& psi
)
{
    return Lagrangian::ddtScheme<Type>::Lagrangianmddt0(deltaT, m, psi);
}


// ************************************************************************* //
