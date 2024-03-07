/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::Field<Type>& Foam::fixedJumpFvPatchField<Type>::jumpRef()
{
    if (!this->cyclicPatch().owner())
    {
        FatalErrorInFunction
            << "The jump field is not available on the neighbour patch"
            << exit(FatalError);
    }

    return jump_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool jumpRequired
)
:
    jumpCyclicFvPatchField<Type>(p, iF, dict),
    jump_(p.size())
{
    if (jumpRequired)
    {
        if (this->cyclicPatch().owner())
        {
            jump_ = Field<Type>("jump", dict, p.size());
        }

        this->evaluateNoUpdateCoeffs();
    }
}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fixedJumpFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    jumpCyclicFvPatchField<Type>(ptf, p, iF, mapper),
    jump_(mapper(ptf.jump_))
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
Foam::tmp<Foam::Field<Type>> Foam::fixedJumpFvPatchField<Type>::jump() const
{
    if (this->cyclicPatch().owner())
    {
        return -jump_;
    }
    else
    {
        const fixedJumpFvPatchField<Type>& nbrField =
            refCast<const fixedJumpFvPatchField<Type>>(this->nbrPatchField());

        return nbrField.jump_;
    }
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    jumpCyclicFvPatchField<Type>::map(ptf, mapper);

    const fixedJumpFvPatchField<Type>& tiptf =
        refCast<const fixedJumpFvPatchField<Type>>(ptf);
    mapper(jump_, tiptf.jump_);
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    jumpCyclicFvPatchField<Type>::reset(ptf);

    const fixedJumpFvPatchField<Type>& tiptf =
        refCast<const fixedJumpFvPatchField<Type>>(ptf);
    jump_.reset(tiptf.jump_);
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    if (this->cyclicPatch().owner())
    {
        writeEntry(os, "jump", jump_);
    }
}


// ************************************************************************* //
