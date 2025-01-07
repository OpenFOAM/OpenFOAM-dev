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

#include "calculatedLagrangianPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::LagrangianPatchField<Type>::calculatedType()
{
    return calculatedLagrangianPatchField<Type>::typeName;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedLagrangianPatchField<Type>::calculatedLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(p, iIo)
{}


template<class Type>
Foam::calculatedLagrangianPatchField<Type>::calculatedLagrangianPatchField
(
    const LagrangianPatch& p,
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianPatchField<Type>(p, iIo, dict)
{}


template<class Type>
Foam::calculatedLagrangianPatchField<Type>::calculatedLagrangianPatchField
(
    const calculatedLagrangianPatchField<Type>& ptf
)
:
    LagrangianPatchField<Type>(ptf)
{}


template<class Type>
Foam::calculatedLagrangianPatchField<Type>::calculatedLagrangianPatchField
(
    const calculatedLagrangianPatchField<Type>& ptf,
    const regIOobject& iIo
)
:
    LagrangianPatchField<Type>(ptf, iIo)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
template<class Type2>
Foam::autoPtr<Foam::LagrangianPatchField<Type>>
Foam::LagrangianPatchField<Type>::NewCalculatedType
(
    const LagrangianPatchField<Type2>& pf
)
{
    typename LagrangianPatchConstructorTable::iterator patchTypeCstrIter =
        LagrangianPatchConstructorTablePtr_->find(pf.patch().type());

    if (patchTypeCstrIter != LagrangianPatchConstructorTablePtr_->end())
    {
        return autoPtr<LagrangianPatchField<Type>>
        (
            patchTypeCstrIter()
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
    else
    {
        return autoPtr<LagrangianPatchField<Type>>
        (
            new calculatedLagrangianPatchField<Type>
            (
                pf.patch(),
                Field<Type>::null()
            )
        );
    }
}


// ************************************************************************* //
