/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "mappedFieldFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFieldFvPatchField<Type>::mappedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase(p.patch()),
    mappedPatchFieldBase<Type>(*this, *this)
{}


template<class Type>
Foam::mappedFieldFvPatchField<Type>::mappedFieldFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchBase(p.patch(), dict),
    mappedPatchFieldBase<Type>(*this, *this, dict)
{}


template<class Type>
Foam::mappedFieldFvPatchField<Type>::mappedFieldFvPatchField
(
    const mappedFieldFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf)
{}


template<class Type>
Foam::mappedFieldFvPatchField<Type>::mappedFieldFvPatchField
(
    const mappedFieldFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFieldFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedFieldFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedFieldFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    fixedValueFvPatchField<Type>::reset(ptf);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedFieldFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(this->mappedField());

    if (debug)
    {
        Info<< "operating on field:" << this->internalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedFieldFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchBase::write(os);
    mappedPatchFieldBase<Type>::write(os);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
