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

#include "nonConformalCyclicLagrangianPatchField.H"
#include "transformField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalCyclicLagrangianPatchField<Type>::
nonConformalCyclicLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(p, iIo),
    nonConformalCyclicPatch_
    (
        refCast<const nonConformalCyclicLagrangianPatch>(p)
    )
{}


template<class Type>
Foam::nonConformalCyclicLagrangianPatchField<Type>::
nonConformalCyclicLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianPatchField<Type>(p, iIo, dict),
    nonConformalCyclicPatch_
    (
        refCast<const nonConformalCyclicLagrangianPatch>(p)
    )
{}


template<class Type>
Foam::nonConformalCyclicLagrangianPatchField<Type>::
nonConformalCyclicLagrangianPatchField
(
    const nonConformalCyclicLagrangianPatchField<Type>& ptf
)
:
    LagrangianPatchField<Type>(ptf),
    nonConformalCyclicPatch_(ptf.nonConformalCyclicPatch_)
{}


template<class Type>
Foam::nonConformalCyclicLagrangianPatchField<Type>::
nonConformalCyclicLagrangianPatchField
(
    const nonConformalCyclicLagrangianPatchField<Type>& ptf,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(ptf, iIo),
    nonConformalCyclicPatch_(ptf.nonConformalCyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::nonConformalCyclicLagrangianPatchField<Type>::evaluate
(
    PstreamBuffers&,
    const LagrangianScalarInternalDynamicField& fraction
)
{
    LagrangianPatchField<Type>::operator=
    (
        nonConformalCyclicPatch_.transform().transform
        (
            this->primitiveSubField().operator const Field<Type>&()
        )
    );
}


// ************************************************************************* //
