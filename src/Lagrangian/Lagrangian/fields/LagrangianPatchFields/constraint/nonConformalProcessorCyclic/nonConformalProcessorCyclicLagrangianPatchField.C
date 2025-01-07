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

#include "nonConformalProcessorCyclicLagrangianPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalProcessorCyclicLagrangianPatchField<Type>::
nonConformalProcessorCyclicLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    processorCyclicLagrangianPatchField<Type>(p, iIo),
    nonConformalProcessorCyclicPatch_
    (
        refCast<const nonConformalProcessorCyclicLagrangianPatch>(p)
    )
{}


template<class Type>
Foam::nonConformalProcessorCyclicLagrangianPatchField<Type>::
nonConformalProcessorCyclicLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
:
    processorCyclicLagrangianPatchField<Type>(p, iIo, dict),
    nonConformalProcessorCyclicPatch_
    (
        refCast<const nonConformalProcessorCyclicLagrangianPatch>(p)
    )
{}


template<class Type>
Foam::nonConformalProcessorCyclicLagrangianPatchField<Type>::
nonConformalProcessorCyclicLagrangianPatchField
(
    const nonConformalProcessorCyclicLagrangianPatchField<Type>& ptf
)
:
    processorCyclicLagrangianPatchField<Type>(ptf),
    nonConformalProcessorCyclicPatch_(ptf.nonConformalProcessorCyclicPatch_)
{}


template<class Type>
Foam::nonConformalProcessorCyclicLagrangianPatchField<Type>::
nonConformalProcessorCyclicLagrangianPatchField
(
    const nonConformalProcessorCyclicLagrangianPatchField<Type>& ptf,
    const regIOobject& iIo
)
:
    processorCyclicLagrangianPatchField<Type>(ptf, iIo),
    nonConformalProcessorCyclicPatch_(ptf.nonConformalProcessorCyclicPatch_)
{}


template<class Type>
Foam::nonConformalProcessorCyclicLagrangianPatchField<Type>::
nonConformalProcessorCyclicLagrangianPatchField
(
    const nonConformalProcessorCyclicLagrangianPatchField<Type>& ptf,
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    processorCyclicLagrangianPatchField<Type>(ptf, p, iIo),
    nonConformalProcessorCyclicPatch_
    (
        refCast<const nonConformalProcessorCyclicLagrangianPatch>(p)
    )
{}


// ************************************************************************* //
