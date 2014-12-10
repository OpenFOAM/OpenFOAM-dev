/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "mappedFixedPushedInternalValueFvPatchField.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFixedPushedInternalValueFvPatchField<Type>::
mappedFixedPushedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mappedFixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::mappedFixedPushedInternalValueFvPatchField<Type>::
mappedFixedPushedInternalValueFvPatchField
(
    const mappedFixedPushedInternalValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedFixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::mappedFixedPushedInternalValueFvPatchField<Type>::
mappedFixedPushedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mappedFixedValueFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::mappedFixedPushedInternalValueFvPatchField<Type>::
mappedFixedPushedInternalValueFvPatchField
(
    const mappedFixedPushedInternalValueFvPatchField<Type>& ptf
)
:
    mappedFixedValueFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::mappedFixedPushedInternalValueFvPatchField<Type>::
mappedFixedPushedInternalValueFvPatchField
(
    const mappedFixedPushedInternalValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mappedFixedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFixedPushedInternalValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Retrieve the neighbour values and assign to this patch boundary field
    mappedFixedValueFvPatchField<Type>::updateCoeffs();

    // Assign the patch internal field to its boundary value
    Field<Type>& intFld = const_cast<Field<Type>&>(this->internalField());
    UIndirectList<Type>(intFld, this->patch().faceCells()) = *this;
}


template<class Type>
void Foam::mappedFixedPushedInternalValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    mappedFixedValueFvPatchField<Type>::write(os);
}


// ************************************************************************* //
