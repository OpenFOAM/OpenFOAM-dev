/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    uniformGradient_()
{}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
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
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    uniformGradient_(Function1<Type>::New("uniformGradient", dict))
{
    this->evaluate();
}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    uniformGradient_(ptf.uniformGradient_, false)
{}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    uniformGradient_(ptf.uniformGradient_, false)
{
    // Evaluate the profile if defined
    if (ptf.uniformGradient_.valid())
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedGradientFvPatchField<Type>::updateCoeffs()
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
void Foam::uniformFixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    writeEntry(os, uniformGradient_());
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
