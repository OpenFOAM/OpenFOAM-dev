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
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
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
    mixedFvPatchField<Type>(p, iF),
    phiName_(ptf.phiName_),
    uniformInletValue_(ptf.uniformInletValue_, false)
{
    this->patchType() = ptf.patchType();

    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();
    this->refValue() = uniformInletValue_->value(t);
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;

    // Map value (unmapped get refValue)
    if (&iF && iF.size())
    {
        fvPatchField<Type>::operator=(this->refValue());
    }
    this->map(ptf, mapper);

}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    uniformInletValue_(DataEntry<Type>::New("uniformInletValue", dict))
{
    const scalar t = this->db().time().timeOutputValue();
    this->refValue() = uniformInletValue_->value(t);

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

    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::uniformInletOutletFvPatchField<Type>::uniformInletOutletFvPatchField
(
    const uniformInletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_),
    uniformInletValue_(ptf.uniformInletValue_, false)
{}


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

    const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            phiName_
        );

    this->valueFraction() = 1.0 - pos(phip);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    this->uniformInletValue_->writeData(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<Type>::autoMap(m);

    // Override
    const scalar t = this->db().time().timeOutputValue();
    this->refValue() = uniformInletValue_->value(t);
}


template<class Type>
void Foam::uniformInletOutletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);

    // Override
    const scalar t = this->db().time().timeOutputValue();
    this->refValue() = uniformInletValue_->value(t);
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
