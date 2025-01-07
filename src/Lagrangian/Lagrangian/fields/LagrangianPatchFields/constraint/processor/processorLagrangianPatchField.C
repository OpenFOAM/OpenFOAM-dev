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

#include "processorLagrangianPatch.H"
#include "processorLagrangianPatchField.H"
#include "transformField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<template<class> class GeoField>
void Foam::processorLagrangianPatchField<Type>::evaluate
(
    PstreamBuffers& pBufs,
    const GeoField<Type>& internalField
)
{
    UIPstream uips(processorPatch_.processorPatch().neighbProcNo(), pBufs);
    Field<Type> field(uips);

    const LagrangianMesh& mesh = processorPatch_.boundaryMesh().mesh();

    mesh.appendSpecifiedField<Type, GeoField>
    (
        processorPatch_.mesh(),
        const_cast<GeoField<Type>&>(internalField),
        field
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorLagrangianPatchField<Type>::processorLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(p, iIo),
    processorPatch_(refCast<const processorLagrangianPatch>(p))
{}


template<class Type>
Foam::processorLagrangianPatchField<Type>::processorLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianPatchField<Type>(p, iIo, dict),
    processorPatch_(refCast<const processorLagrangianPatch>(p))
{}


template<class Type>
Foam::processorLagrangianPatchField<Type>::processorLagrangianPatchField
(
    const processorLagrangianPatchField<Type>& ptf
)
:
    LagrangianPatchField<Type>(ptf),
    processorPatch_(ptf.processorPatch_)
{}


template<class Type>
Foam::processorLagrangianPatchField<Type>::processorLagrangianPatchField
(
    const processorLagrangianPatchField<Type>& ptf,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(ptf, iIo),
    processorPatch_(ptf.processorPatch_)
{}


template<class Type>
Foam::processorLagrangianPatchField<Type>::processorLagrangianPatchField
(
    const processorLagrangianPatchField<Type>& ptf,
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(ptf, p, iIo),
    processorPatch_(refCast<const processorLagrangianPatch>(p))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::processorLagrangianPatchField<Type>::initEvaluate
(
    PstreamBuffers& pBufs,
    const LagrangianScalarInternalDynamicField& fraction
)
{
    UOPstream(processorPatch_.processorPatch().neighbProcNo(), pBufs)()
        << this->primitiveSubField();
}


template<class Type>
void Foam::processorLagrangianPatchField<Type>::evaluate
(
    PstreamBuffers& pBufs,
    const LagrangianScalarInternalDynamicField& fraction
)
{
    if (notNull(this->internalField_))
    {
        return
            evaluate<LagrangianInternalDynamicField>
            (
                pBufs,
                this->internalField_
            );
    }

    if (notNull(this->internalNonDynamicField_))
    {
        return
            evaluate<LagrangianInternalField>
            (
                pBufs,
                this->internalNonDynamicField_
            );
    }

    FatalErrorInFunction
        << "Internal field " << this->internalIo_.name() << " is not of type "
        << LagrangianInternalDynamicField<Type>::typeName
        << exit(FatalError);
}


// ************************************************************************* //
