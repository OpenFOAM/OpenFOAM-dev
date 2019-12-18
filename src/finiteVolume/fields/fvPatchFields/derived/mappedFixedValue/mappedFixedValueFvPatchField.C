/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "mappedFixedValueFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFixedValueFvPatchField<Type>::mappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this)
{}


template<class Type>
Foam::mappedFixedValueFvPatchField<Type>::mappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this, dict)
{}


template<class Type>
Foam::mappedFixedValueFvPatchField<Type>::mappedFixedValueFvPatchField
(
    const mappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchFieldBase<Type>(this->mapper(p, iF), *this, ptf)
{}


template<class Type>
Foam::mappedFixedValueFvPatchField<Type>::mappedFixedValueFvPatchField
(
    const mappedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    mappedPatchFieldBase<Type>(ptf)
{}


template<class Type>
Foam::mappedFixedValueFvPatchField<Type>::mappedFixedValueFvPatchField
(
    const mappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchFieldBase<Type>(this->mapper(this->patch(), iF), *this, ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::mappedPatchBase& Foam::mappedFixedValueFvPatchField<Type>::mapper
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.patch().name()
            << " of field " << iF.name()
            << " in file " << iF.objectPath()
            << exit(FatalError);
    }
    return refCast<const mappedPatchBase>(p.patch());
}


template<class Type>
void Foam::mappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(this->mappedField());

    if (debug)
    {
        Info<< "mapped on field:"
            << this->internalField().name()
            << " patch:" << this->patch().name()
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchFieldBase<Type>::write(os);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
