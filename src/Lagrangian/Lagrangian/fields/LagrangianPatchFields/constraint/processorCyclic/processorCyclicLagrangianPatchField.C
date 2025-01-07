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

#include "processorCyclicLagrangianPatchField.H"
#include "transformField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorCyclicLagrangianPatchField<Type>::
processorCyclicLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    processorLagrangianPatchField<Type>(p, iIo),
    processorCyclicPatch_(refCast<const processorCyclicLagrangianPatch>(p))
{}


template<class Type>
Foam::processorCyclicLagrangianPatchField<Type>::
processorCyclicLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
:
    processorLagrangianPatchField<Type>(p, iIo, dict),
    processorCyclicPatch_(refCast<const processorCyclicLagrangianPatch>(p))
{}


template<class Type>
Foam::processorCyclicLagrangianPatchField<Type>::
processorCyclicLagrangianPatchField
(
    const processorCyclicLagrangianPatchField<Type>& ptf
)
:
    processorLagrangianPatchField<Type>(ptf),
    processorCyclicPatch_(ptf.processorCyclicPatch_)
{}


template<class Type>
Foam::processorCyclicLagrangianPatchField<Type>::
processorCyclicLagrangianPatchField
(
    const processorCyclicLagrangianPatchField<Type>& ptf,
    const regIOobject& iIo
)
:
    processorLagrangianPatchField<Type>(ptf, iIo),
    processorCyclicPatch_(ptf.processorCyclicPatch_)
{}


template<class Type>
Foam::processorCyclicLagrangianPatchField<Type>::
processorCyclicLagrangianPatchField
(
    const processorCyclicLagrangianPatchField<Type>& ptf,
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    processorLagrangianPatchField<Type>(ptf, p, iIo),
    processorCyclicPatch_(refCast<const processorCyclicLagrangianPatch>(p))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::processorCyclicLagrangianPatchField<Type>::evaluate
(
    PstreamBuffers& pBufs,
    const LagrangianScalarInternalDynamicField& fraction
)
{
    processorLagrangianPatchField<Type>::evaluate(pBufs, fraction);

    LagrangianPatchField<Type>::operator=
    (
        processorCyclicPatch_.transform().transform
        (
            this->primitiveSubField().operator const Field<Type>&()
        )
    );
}


// ************************************************************************* //
