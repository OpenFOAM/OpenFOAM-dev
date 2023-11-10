/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "nonConformalProcessorCyclicPointPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalProcessorCyclicPointPatchField<Type>::
nonConformalProcessorCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    procPatch_(refCast<const nonConformalProcessorCyclicPointPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicPointPatchField<Type>::
nonConformalProcessorCyclicPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const nonConformalProcessorCyclicPointPatch>(p))
{
}


template<class Type>
Foam::nonConformalProcessorCyclicPointPatchField<Type>::
nonConformalProcessorCyclicPointPatchField
(
    const nonConformalProcessorCyclicPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const fieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    procPatch_
    (
        refCast<const nonConformalProcessorCyclicPointPatch>(ptf.patch())
    )
{}


template<class Type>
Foam::nonConformalProcessorCyclicPointPatchField<Type>::
nonConformalProcessorCyclicPointPatchField
(
    const nonConformalProcessorCyclicPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    procPatch_
    (
        refCast<const nonConformalProcessorCyclicPointPatch>(ptf.patch())
    )
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalProcessorCyclicPointPatchField<Type>::
~nonConformalProcessorCyclicPointPatchField()
{}


// ************************************************************************* //
