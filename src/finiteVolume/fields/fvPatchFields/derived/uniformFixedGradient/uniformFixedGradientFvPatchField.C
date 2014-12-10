/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "uniformFixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    uniformGradient_()
{}


template<class Type>
uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedGradientFvPatchField<Type>(p, iF, fld),
    uniformGradient_()
{}


template<class Type>
uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    uniformGradient_(ptf.uniformGradient_().clone().ptr())
{
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();
    this->gradient() = uniformGradient_->value(t);
}


template<class Type>
uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    uniformGradient_(DataEntry<Type>::New("uniformGradient", dict))
{
    if (dict.found("gradient"))
    {
        this->gradient() = Field<Type>("gradient", dict, p.size());
    }
    else
    {
        const scalar t = this->db().time().timeOutputValue();
        this->gradient() = uniformGradient_->value(t);
    }
}


template<class Type>
uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    uniformGradient_
    (
        ptf.uniformGradient_.valid()
      ? ptf.uniformGradient_().clone().ptr()
      : NULL
    )
{}


template<class Type>
uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    uniformGradient_
    (
        ptf.uniformGradient_.valid()
      ? ptf.uniformGradient_().clone().ptr()
      : NULL
    )
{
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();

    if (ptf.uniformGradient_.valid())
    {
        this->gradient() = uniformGradient_->value(t);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void uniformFixedGradientFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    this->gradient() = uniformGradient_->value(t);

    fixedGradientFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void uniformFixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    uniformGradient_->writeData(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
