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

#include "internalLagrangianPatchField.H"
#include "transformField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::internalLagrangianPatchField<Type>::internalLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(p, iIo),
    internalPatch_(refCast<const internalLagrangianPatch>(p))
{}


template<class Type>
Foam::internalLagrangianPatchField<Type>::internalLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianPatchField<Type>(p, iIo, dict),
    internalPatch_(refCast<const internalLagrangianPatch>(p))
{}


template<class Type>
Foam::internalLagrangianPatchField<Type>::internalLagrangianPatchField
(
    const internalLagrangianPatchField<Type>& ptf
)
:
    LagrangianPatchField<Type>(ptf),
    internalPatch_(ptf.internalPatch_)
{}


template<class Type>
Foam::internalLagrangianPatchField<Type>::internalLagrangianPatchField
(
    const internalLagrangianPatchField<Type>& ptf,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(ptf, iIo),
    internalPatch_(ptf.internalPatch_)
{}


// ************************************************************************* //
