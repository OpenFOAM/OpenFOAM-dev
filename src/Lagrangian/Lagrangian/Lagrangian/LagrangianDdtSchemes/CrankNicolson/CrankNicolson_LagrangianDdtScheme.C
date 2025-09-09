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

#include "CrankNicolson_LagrangianDdtScheme.H"
#include "LagrangianFields.H"
#include "LagrangianEqn.H"
#include "LagrangianModels.H"
#include "NaNLagrangianFieldSource.H"
#include "internalLagrangianFieldSource.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::word Foam::Lagrangian::ddtSchemes::CrankNicolson<Type>::S0Name
(
    const LagrangianSubSubField<Type>& psi
) const
{
    return typedName("S0(" + psi.mesh().complete(psi.name()) + ")");
}


template<class Type>
Foam::word Foam::Lagrangian::ddtSchemes::CrankNicolson<Type>::deltaTSp0Name
(
    const LagrangianSubSubField<Type>& psi
) const
{
    return typedName("deltaTSp0(" + psi.mesh().complete(psi.name()) + ")");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::Lagrangian::ddtSchemes::CrankNicolson<Type>::LagrangianmInitDdt
(
    const dimensionSet& mDims,
    const LagrangianSubSubField<Type>& psi,
    const bool instantaneousDdt
)
{
    const LagrangianSubMesh& subMesh = psi.mesh();
    const LagrangianMesh& mesh = subMesh.mesh();

    const word S0Name = this->S0Name(psi);

    if (!mesh.foundObject<LagrangianDynamicField<Type>>(S0Name))
    {
        regIOobject::store
        (
            new LagrangianDynamicField<Type>
            (
                IOobject
                (
                    S0Name,
                    mesh.time().name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimensioned<Type>(mDims*psi.dimensions()/dimTime, Zero),
                wordList
                (
                    mesh.boundary().size(),
                    calculatedLagrangianPatchField<Type>::typeName
                ),
                wordList::null(),
                LagrangianModels::New(mesh).modelTypeFieldSourceTypes
                <
                    LagrangianInjection,
                    NaNLagrangianFieldSource<Type>,
                    LagrangianSource,
                    internalLagrangianFieldSource<Type>
                >()
            )
        );
    }

    if (!instantaneousDdt) return true;

    const word deltaTSp0Name = this->deltaTSp0Name(psi);

    if (!mesh.foundObject<LagrangianDynamicField<scalar>>(deltaTSp0Name))
    {
        regIOobject::store
        (
            new LagrangianDynamicField<scalar>
            (
                IOobject
                (
                    deltaTSp0Name,
                    mesh.time().name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimensioned<scalar>(mDims, Zero),
                wordList
                (
                    mesh.boundary().size(),
                    calculatedLagrangianPatchField<scalar>::typeName
                ),
                wordList::null(),
                LagrangianModels::New(mesh).modelTypeFieldSourceTypes
                <
                    LagrangianInjection,
                    NaNLagrangianFieldSource<scalar>,
                    LagrangianSource,
                    internalLagrangianFieldSource<scalar>
                >()
            )
        );
    }

    return true;
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtSchemes::CrankNicolson<Type>::LagrangianmNoDdt
(
    const LagrangianSubScalarField& deltaT,
    const dimensionSet& mDims,
    const LagrangianSubSubField<Type>& psi
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();
    const LagrangianMesh& mesh = subMesh.mesh();

    const word S0Name = this->S0Name(psi);
    const word deltaTSp0Name = this->deltaTSp0Name(psi);

    LagrangianDynamicField<Type>& S0 =
        mesh.lookupObjectRef<LagrangianDynamicField<Type>>(S0Name);
    LagrangianDynamicField<scalar>& deltaTSp0 =
        mesh.foundObject<LagrangianDynamicField<scalar>>(deltaTSp0Name)
      ? mesh.lookupObjectRef<LagrangianDynamicField<scalar>>(deltaTSp0Name)
      : NullObjectNonConstRef<LagrangianDynamicField<scalar>>();

    const bool isNone = subMesh.group() == LagrangianGroup::none;

    tmp<LagrangianEqn<Type>> tEqn
    (
        new LagrangianEqn<Type>
        (
            !isNone ? deltaT/2 : tmp<LagrangianSubScalarField>(deltaT),
            psi,
            deltaTSp0,
            S0
        )
    );

    if (!isNone) tEqn.ref().Su += subMesh.sub(S0);

    return tEqn;
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtSchemes::CrankNicolson<Type>::LagrangianmDdt
(
    const LagrangianSubScalarField& deltaT,
    LagrangianSubSubField<Type>& psi
)
{
    // Return the standard time-derivative equation, but add the
    // no-time-derivative matrix to it to get the time-step and trigger the
    // necessary side-effects, and also add an empty equation so that the name
    // and the non-const field reference are set
    return
        Lagrangian::ddtScheme<Type>::Lagrangianmddt(psi)
      + CrankNicolson<Type>::LagrangianmNoDdt(deltaT, dimless, psi)
      + tmp<LagrangianEqn<Type>>(new LagrangianEqn<Type>(psi));
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtSchemes::CrankNicolson<Type>::LagrangianmDdt
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& m,
    LagrangianSubSubField<Type>& psi
)
{
    // Return the standard time-derivative equation, but add the
    // no-time-derivative matrix to it to get the time-step and trigger the
    // necessary side-effects, and also add an empty equation so that the name
    // and the non-const field reference are set
    return
        Lagrangian::ddtScheme<Type>::Lagrangianmddt(m, psi)
      + CrankNicolson<Type>::LagrangianmNoDdt(deltaT, m.dimensions(), psi)
      + tmp<LagrangianEqn<Type>>(new LagrangianEqn<Type>(psi));
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::Lagrangian::ddtSchemes::CrankNicolson<Type>::LagrangiancDdt
(
    const LagrangianSubSubField<Type>& psi
)
{
    const LagrangianSubMesh& subMesh = psi.mesh();
    const LagrangianMesh& mesh = subMesh.mesh();

    return
      - subMesh.sub
        (
            mesh.lookupObject<LagrangianDynamicField<Type>>
            (
                S0Name(psi)
            )
        )
       /subMesh.sub
        (
            mesh.lookupObject<LagrangianDynamicField<scalar>>
            (
                deltaTSp0Name(psi)
            )
        );
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::Lagrangian::ddtSchemes::CrankNicolson<Type>::LagrangiancDdt
(
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& psi
)
{
    // This is an approximation if the time-coefficient (deltaTSp0) is not
    // "just" m. Some models may add to the time-coefficient; a virtual mass
    // model, for example, adds the added mass. The time-derivative of these
    // additions is not known and cannot therefore be factored out. This
    // approximation assumes that those additional terms are proportional to m
    // and neglects the time-derivative of the coefficient. For a
    // constant-coefficient virtual mass and constant relative density this
    // approach is exact.

    return m*LagrangiancDdt(psi);
}


// ************************************************************************* //
