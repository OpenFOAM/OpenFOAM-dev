/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2026 OpenFOAM Foundation
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

#include "functionalFixedValueFvPatchField.H"
#include "DimensionedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::functionalFixedValueFvPatchField<Type>::functionalFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    dimensionedValue_
    (
        IOobject
        (
            iF.name(),
            iF.time().name(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        p,
        iF.dimensions(),
        *this
    ),
    funcPtr_
    (
        DimensionedFieldFunction<DimensionedFvPatchField<Type>>::New
        (
            dict.subDict("function"),
            dimensionedValue_
        )
    )
{}


template<class Type>
Foam::functionalFixedValueFvPatchField<Type>::functionalFixedValueFvPatchField
(
    const functionalFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper, true),
    dimensionedValue_
    (
        ptf.dimensionedValue_,
        p,
        ptf.dimensionedValue_.dimensions(),
        *this
    ),
    funcPtr_(ptf.funcPtr_->clone(dimensionedValue_))
{}


template<class Type>
Foam::functionalFixedValueFvPatchField<Type>::functionalFixedValueFvPatchField
(
    const functionalFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, fvMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    dimensionedValue_
    (
        ptf.dimensionedValue_,
        ptf.patch(),
        ptf.dimensionedValue_.dimensions(),
        *this
    ),
    funcPtr_(ptf.funcPtr_->clone(dimensionedValue_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::functionalFixedValueFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    if (mapper.direct())
    {
        fixedValueFvPatchField<Type>::map(ptf, mapper);
    }
    else
    {
        this->setSize(this->patch().size());
        dimensionedValue_.reset(*this);
        funcPtr_->reset();
    }
}


template<class Type>
void Foam::functionalFixedValueFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    this->setSize(this->patch().size());
    dimensionedValue_.reset(*this);
    funcPtr_->reset();
}


template<class Type>
void Foam::functionalFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    funcPtr_->update();

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::functionalFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "function", *funcPtr_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
