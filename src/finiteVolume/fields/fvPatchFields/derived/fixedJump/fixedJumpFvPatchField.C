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

#include "fixedJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicFvPatchField<Type>(p, iF),
    jump_(this->size(), pTraits<Type>::zero)
{}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fixedJumpFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpCyclicFvPatchField<Type>(ptf, p, iF, mapper),
    jump_(ptf.jump_, mapper)
{}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicFvPatchField<Type>(p, iF),
    jump_(p.size(), pTraits<Type>::zero)
{
    if (this->cyclicPatch().owner())
    {
        jump_ = Field<Type>("jump", dict, p.size());
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::blocking);
    }
}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fixedJumpFvPatchField<Type>& ptf
)
:
    jumpCyclicFvPatchField<Type>(ptf),
    jump_(ptf.jump_)
{}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fixedJumpFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicFvPatchField<Type>(ptf, iF),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fixedJumpFvPatchField<Type>::jump() const
{
    if (this->cyclicPatch().owner())
    {
        return jump_;
    }
    else
    {
        return refCast<const fixedJumpFvPatchField<Type> >
        (
            this->neighbourPatchField()
        ).jump();
    }
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpCyclicFvPatchField<Type>::autoMap(m);
    jump_.autoMap(m);
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    jumpCyclicFvPatchField<Type>::rmap(ptf, addr);

    const fixedJumpFvPatchField<Type>& tiptf =
        refCast<const fixedJumpFvPatchField<Type> >(ptf);
    jump_.rmap(tiptf.jump_, addr);
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("patchType") << this->interfaceFieldType()
        << token::END_STATEMENT << nl;

    if (this->cyclicPatch().owner())
    {
        jump_.writeEntry("jump", os);
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
