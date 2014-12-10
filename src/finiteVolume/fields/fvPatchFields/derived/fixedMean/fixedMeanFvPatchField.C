/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "fixedMeanFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fixedMeanFvPatchField<Type>::fixedMeanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    meanValue_(pTraits<Type>::zero)
{}


template<class Type>
fixedMeanFvPatchField<Type>::fixedMeanFvPatchField
(
    const fixedMeanFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    meanValue_(ptf.meanValue_)
{}


template<class Type>
fixedMeanFvPatchField<Type>::fixedMeanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    meanValue_(pTraits<Type>(dict.lookup("meanValue")))
{}


template<class Type>
fixedMeanFvPatchField<Type>::fixedMeanFvPatchField
(
    const fixedMeanFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    meanValue_(ptf.meanValue_)
{}


template<class Type>
fixedMeanFvPatchField<Type>::fixedMeanFvPatchField
(
    const fixedMeanFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    meanValue_(ptf.meanValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void fixedMeanFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    Field<Type> newValues(this->patchInternalField());

    Type meanValuePsi =
        gSum(this->patch().magSf()*newValues)
       /gSum(this->patch().magSf());

    if (mag(meanValue_) > SMALL && mag(meanValuePsi)/mag(meanValue_) > 0.5)
    {
        newValues *= mag(meanValue_)/mag(meanValuePsi);
    }
    else
    {
        newValues += (meanValue_ - meanValuePsi);
    }

    this->operator==(newValues);

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void fixedMeanFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("meanValue") << meanValue_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
