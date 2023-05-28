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

#include "mappedValueFvPatchField.H"
#include "mappedPolyPatch.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedValueFvPatchField<Type>::mappedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::mappedValueFvPatchField<Type>::mappedValueFvPatchField
(
    const mappedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::mappedValueFvPatchField<Type>::mappedValueFvPatchField
(
    const mappedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedValueFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField<Type>::map(ptf, mapper);
    mappedFvPatchField<Type>::clearOut();
}


template<class Type>
void Foam::mappedValueFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    fixedValueFvPatchField<Type>::reset(ptf);
    mappedFvPatchField<Type>::clearOut();
}


template<class Type>
void Foam::mappedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(this->mappedValues(this->nbrPatchField()));

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    mappedFvPatchField<Type>::write(os);

    writeEntry(os, "value", *this);
}


// ************************************************************************* //
