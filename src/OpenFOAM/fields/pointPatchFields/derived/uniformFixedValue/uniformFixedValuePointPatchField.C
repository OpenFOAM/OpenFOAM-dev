/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "uniformFixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    uniformValue_()
{}


template<class Type>
uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, false),
    uniformValue_(DataEntry<Type>::New("uniformValue", dict))
{
    if (dict.found("value"))
    {
        fixedValuePointPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        const scalar t = this->db().time().timeOutputValue();
        fixedValuePointPatchField<Type>::operator=(uniformValue_->value(t));
    }
}


template<class Type>
uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper),
    uniformValue_(ptf.uniformValue_().clone().ptr())
{
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();
    fixedValuePointPatchField<Type>::operator=(uniformValue_->value(t));
}


template<class Type>
uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_().clone().ptr())
{}


template<class Type>
uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_().clone().ptr())
{
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();
    fixedValuePointPatchField<Type>::operator==(uniformValue_->value(t));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void uniformFixedValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    fixedValuePointPatchField<Type>::operator==(uniformValue_->value(t));

    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void uniformFixedValuePointPatchField<Type>::
write(Ostream& os) const
{
    // Note: write value
    fixedValuePointPatchField<Type>::write(os);
    uniformValue_->writeData(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
