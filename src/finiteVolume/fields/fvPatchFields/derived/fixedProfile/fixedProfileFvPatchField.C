/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "fixedProfileFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedProfileFvPatchField<Type>::fixedProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    profile_(),
    origin_(0),
    direction_(Zero)
{}


template<class Type>
Foam::fixedProfileFvPatchField<Type>::fixedProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFvPatchField<Type>(p, iF, fld),
    profile_(),
    origin_(0),
    direction_(Zero)
{}


template<class Type>
Foam::fixedProfileFvPatchField<Type>::fixedProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    profile_(Function1<Type>::New("profile", dict)),
    origin_(dict.lookup<scalar>("origin")),
    direction_(dict.lookup("direction"))
{
    if (mag(direction_) == 0)
    {
        FatalErrorInFunction
            << "The direction must be non-zero"
            << abort(FatalError);
    }

    // Normalise the direction
    direction_ /= mag(direction_);

    // Evaluate profile
    this->evaluate();
}


template<class Type>
Foam::fixedProfileFvPatchField<Type>::fixedProfileFvPatchField
(
    const fixedProfileFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(p, iF), // Don't map
    profile_(ptf.profile_, false),
    origin_(ptf.origin_),
    direction_(ptf.direction_)
{
    // Evaluate profile since value not mapped
    this->evaluate();
}


template<class Type>
Foam::fixedProfileFvPatchField<Type>::fixedProfileFvPatchField
(
    const fixedProfileFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    profile_(ptf.profile_, false),
    origin_(ptf.origin_),
    direction_(ptf.direction_)
{}


template<class Type>
Foam::fixedProfileFvPatchField<Type>::fixedProfileFvPatchField
(
    const fixedProfileFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    profile_(ptf.profile_, false),
    origin_(ptf.origin_),
    direction_(ptf.direction_)
{
    // Evaluate the profile if defined
    if (ptf.profile_.valid())
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedProfileFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalarField dirCmpt((direction_ & this->patch().Cf()) - origin_);

    fvPatchField<Type>::operator==(profile_->value(dirCmpt));

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::fixedProfileFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, profile_());
    writeEntry(os, "direction", direction_);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
