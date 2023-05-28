/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "mappedValueAndPatchInternalValueFvPatchField.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedValueAndPatchInternalValueFvPatchField<Type>::
mappedValueAndPatchInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mappedValueFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::mappedValueAndPatchInternalValueFvPatchField<Type>::
mappedValueAndPatchInternalValueFvPatchField
(
    const mappedValueAndPatchInternalValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::mappedValueAndPatchInternalValueFvPatchField<Type>::
mappedValueAndPatchInternalValueFvPatchField
(
    const mappedValueAndPatchInternalValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mappedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedValueAndPatchInternalValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(this->mappedValues(this->nbrPatchField()));

    UIndirectList<Type>
    (
        const_cast<Field<Type>&>(this->primitiveField()),
        this->patch().faceCells()
    ) = this->mappedValues(this->nbrPatchField().patchInternalField());

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedValueAndPatchInternalValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    mappedValueFvPatchField<Type>::write(os);
}


// ************************************************************************* //
