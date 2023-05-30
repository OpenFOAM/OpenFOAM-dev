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

#include "fixedGradientFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    gradient_(p.size(), Zero)
{}


template<class Type>
Foam::fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool gradientRequired
)
:
    fvPatchField<Type>(p, iF, dict, false),
    gradient_(p.size())
{
    if (gradientRequired)
    {
        if (dict.found("gradient"))
        {
            gradient_ = Field<Type>("gradient", dict, p.size());
            evaluate();
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Essential entry 'gradient' missing"
                << exit(FatalIOError);
        }
    }
}


template<class Type>
Foam::fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
    const fixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<Type>(ptf, p, iF, mapper, mappingRequired),
    gradient_(p.size())
{
    if (mappingRequired)
    {
        mapper(gradient_, ptf.gradient_);
    }
}


template<class Type>
Foam::fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
    const fixedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedGradientFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fvPatchFieldMapper& mapper
)
{
    fvPatchField<Type>::map(ptf, mapper);

    const fixedGradientFvPatchField<Type>& fgptf =
        refCast<const fixedGradientFvPatchField<Type>>(ptf);

    mapper(gradient_, fgptf.gradient_);
}


template<class Type>
void Foam::fixedGradientFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::reset(ptf);

    const fixedGradientFvPatchField<Type>& fgptf =
        refCast<const fixedGradientFvPatchField<Type>>(ptf);

    gradient_.reset(fgptf.gradient_);
}


template<class Type>
void Foam::fixedGradientFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedGradientFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>(new Field<Type>(this->size(), pTraits<Type>::one));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedGradientFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return gradient()/this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedGradientFvPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedGradientFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return gradient();
}


template<class Type>
void Foam::fixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "gradient", gradient_);
}


// ************************************************************************* //
