/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "fixedNormalSlipFvPatchField.H"
#include "symmTransformField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF),
    fixedValue_(p.size(), Zero)
{}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF),
    fixedValue_("fixedValue", dict, p.size())
{
    evaluate();
}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper),
    fixedValue_(mapper(ptf.fixedValue_))
{}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF),
    fixedValue_(ptf.fixedValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    transformFvPatchField<Type>::autoMap(m);
    m(fixedValue_, fixedValue_);
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    transformFvPatchField<Type>::rmap(ptf, addr);

    const fixedNormalSlipFvPatchField<Type>& dmptf =
        refCast<const fixedNormalSlipFvPatchField<Type>>(ptf);

    fixedValue_.rmap(dmptf.fixedValue_, addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedNormalSlipFvPatchField<Type>::snGrad() const
{
    const vectorField nHat(this->patch().nf());
    const Field<Type> pif(this->patchInternalField());

    return
    (
        (nHat*(nHat & fixedValue_) + transform(I - sqr(nHat), pif)) - pif
    )*this->patch().deltaCoeffs();
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const vectorField nHat(this->patch().nf());

    Field<Type>::operator=
    (
        nHat*(nHat & fixedValue_)
      + transform(I - sqr(nHat), this->patchInternalField())
    );

    transformFvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedNormalSlipFvPatchField<Type>::snGradTransformDiag() const
{
    const vectorField nHat(this->patch().nf());
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::write(Ostream& os) const
{
    transformFvPatchField<Type>::write(os);
    writeEntry(os, "fixedValue", fixedValue_);
}


// ************************************************************************* //
