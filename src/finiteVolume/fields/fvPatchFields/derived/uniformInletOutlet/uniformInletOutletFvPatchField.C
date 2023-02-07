/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "uniformInletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phi")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict, false),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    uniformInletValue_(Function1<Type>::New("uniformInletValue", dict))
{
    this->refValue() =
        uniformInletValue_->value(this->db().time().userTimeValue());

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->refValue());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const uniformInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper, false), // Don't map
    phiName_(ptf.phiName_),
    uniformInletValue_(ptf.uniformInletValue_, false)
{
    // Evaluate refValue since not mapped
    this->refValue() =
        uniformInletValue_->value(this->db().time().userTimeValue());

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

    // Initialise the patch value to the refValue
    fvPatchField<Type>::operator=(this->refValue());

    mapper(*this, ptf);
}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const uniformInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_),
    uniformInletValue_(ptf.uniformInletValue_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->refValue() =
        uniformInletValue_->value(this->db().time().userTimeValue());

    const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            phiName_
        );

    this->valueFraction() = neg(phip);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (phiName_ != "phi")
    {
        writeEntry(os, "phi", phiName_);
    }
    writeEntry(os, this->uniformInletValue_());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchField<Type>::map(ptf, mapper);

    // Override
    this->refValue() =
        uniformInletValue_->value(this->db().time().userTimeValue());
}


template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    mixedFvPatchField<Type>::reset(ptf);

    // Override
    this->refValue() =
        uniformInletValue_->value(this->db().time().userTimeValue());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// ************************************************************************* //
