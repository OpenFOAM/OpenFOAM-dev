/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool valuesRequired
)
:
    fvPatchField<Type>(p, iF, dict, false),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{
    if (valuesRequired)
    {
        if (dict.found("refValue"))
        {
            refValue_ = Field<Type>("refValue", dict, p.size());
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Essential entry 'refValue' missing"
                << exit(FatalIOError);
        }

        if (dict.found("refGradient"))
        {
            refGrad_ = Field<Type>("refGradient", dict, p.size());
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Essential entry 'refGradient' missing"
                << exit(FatalIOError);
        }

        if (dict.found("valueFraction"))
        {
            valueFraction_ = Field<scalar>("valueFraction", dict, p.size());
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Essential entry 'valueFraction' missing"
                << exit(FatalIOError);
        }

        evaluate();
    }
}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const mixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<Type>(ptf, p, iF, mapper, mappingRequired),
    refValue_(mapper(ptf.refValue_)),
    refGrad_(mapper(ptf.refGrad_)),
    valueFraction_(mapper(ptf.valueFraction_))
{}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const mixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mixedFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fvPatchFieldMapper& mapper
)
{
    fvPatchField<Type>::map(ptf, mapper);

    const mixedFvPatchField<Type>& mptf =
        refCast<const mixedFvPatchField<Type>>(ptf);

    mapper(refValue_, mptf.refValue_);
    mapper(refGrad_, mptf.refGrad_);
    mapper(valueFraction_, mptf.valueFraction_);
}


template<class Type>
void Foam::mixedFvPatchField<Type>::reset(const fvPatchField<Type>& ptf)
{
    fvPatchField<Type>::reset(ptf);

    const mixedFvPatchField<Type>& mptf =
        refCast<const mixedFvPatchField<Type>>(ptf);

    refValue_.reset(mptf.refValue_);
    refGrad_.reset(mptf.refGrad_);
    valueFraction_.reset(mptf.valueFraction_);
}


template<class Type>
void Foam::mixedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        valueFraction_*refValue_
      +
        (1.0 - valueFraction_)*
        (
            this->patchInternalField()
          + refGrad_/this->patch().deltaCoeffs()
        )
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::snGrad() const
{
    return
        valueFraction_
       *(refValue_ - this->patchInternalField())
       *this->patch().deltaCoeffs()
      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one)*(1.0 - valueFraction_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         valueFraction_*refValue_
       + (1.0 - valueFraction_)*refGrad_/this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*valueFraction_*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        valueFraction_*this->patch().deltaCoeffs()*refValue_
      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
void Foam::mixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "refValue", refValue_);
    writeEntry(os, "refGradient", refGrad_);
    writeEntry(os, "valueFraction", valueFraction_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
