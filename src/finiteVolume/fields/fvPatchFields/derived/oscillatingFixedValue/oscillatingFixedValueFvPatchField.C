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

#include "oscillatingFixedValueFvPatchField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
scalar oscillatingFixedValueFvPatchField<Type>::currentScale() const
{
    const scalar t = this->db().time().timeOutputValue();
    const scalar a = amplitude_->value(t);
    const scalar f = frequency_->value(t);

    return 1.0 + a*sin(constant::mathematical::twoPi*f*t);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
oscillatingFixedValueFvPatchField<Type>::oscillatingFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    refValue_(p.size()),
    offset_(pTraits<Type>::zero),
    amplitude_(),
    frequency_(),
    curTimeIndex_(-1)
{}


template<class Type>
oscillatingFixedValueFvPatchField<Type>::oscillatingFixedValueFvPatchField
(
    const oscillatingFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    offset_(ptf.offset_),
    amplitude_(ptf.amplitude_().clone().ptr()),
    frequency_(ptf.frequency_().clone().ptr()),
    curTimeIndex_(-1)
{}


template<class Type>
oscillatingFixedValueFvPatchField<Type>::oscillatingFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    refValue_("refValue", dict, p.size()),
    offset_(dict.lookupOrDefault<Type>("offset", pTraits<Type>::zero)),
    amplitude_(DataEntry<scalar>::New("amplitude", dict)),
    frequency_(DataEntry<scalar>::New("frequency", dict)),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==
        (
            refValue_*currentScale()
          + offset_
        );
    }
}


template<class Type>
oscillatingFixedValueFvPatchField<Type>::oscillatingFixedValueFvPatchField
(
    const oscillatingFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    refValue_(ptf.refValue_),
    offset_(ptf.offset_),
    amplitude_(ptf.amplitude_().clone().ptr()),
    frequency_(ptf.frequency_().clone().ptr()),
    curTimeIndex_(-1)
{}


template<class Type>
oscillatingFixedValueFvPatchField<Type>::oscillatingFixedValueFvPatchField
(
    const oscillatingFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    offset_(ptf.offset_),
    amplitude_(ptf.amplitude_().clone().ptr()),
    frequency_(ptf.frequency_().clone().ptr()),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void oscillatingFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    refValue_.autoMap(m);
}


template<class Type>
void oscillatingFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const oscillatingFixedValueFvPatchField<Type>& tiptf =
        refCast<const oscillatingFixedValueFvPatchField<Type> >(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


template<class Type>
void oscillatingFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        fixedValueFvPatchField<Type>::operator==
        (
            refValue_*currentScale()
          + offset_
        );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void oscillatingFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
    amplitude_->writeData(os);
    frequency_->writeData(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
