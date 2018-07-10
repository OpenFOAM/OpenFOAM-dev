/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "fixedJumpAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedJumpAMIFvPatchField<Type>::fixedJumpAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicAMIFvPatchField<Type>(p, iF),
    jump_(this->size(), Zero)
{}


template<class Type>
Foam::fixedJumpAMIFvPatchField<Type>::fixedJumpAMIFvPatchField
(
    const fixedJumpAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpCyclicAMIFvPatchField<Type>(ptf, p, iF, mapper),
    jump_(ptf.jump_, mapper)
{}


template<class Type>
Foam::fixedJumpAMIFvPatchField<Type>::fixedJumpAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicAMIFvPatchField<Type>(p, iF),
    jump_(p.size(), Zero)
{
    if (this->cyclicAMIPatch().owner())
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
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::fixedJumpAMIFvPatchField<Type>::fixedJumpAMIFvPatchField
(
    const fixedJumpAMIFvPatchField<Type>& ptf
)
:
    jumpCyclicAMIFvPatchField<Type>(ptf),
    jump_(ptf.jump_)
{}


template<class Type>
Foam::fixedJumpAMIFvPatchField<Type>::fixedJumpAMIFvPatchField
(
    const fixedJumpAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicAMIFvPatchField<Type>(ptf, iF),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fixedJumpAMIFvPatchField<Type>::jump() const
{
    if (this->cyclicAMIPatch().owner())
    {
        return jump_;
    }
    else
    {
        const fixedJumpAMIFvPatchField& nbrPatch =
            refCast<const fixedJumpAMIFvPatchField<Type>>
            (
                this->neighbourPatchField()
            );

        if (this->cyclicAMIPatch().applyLowWeightCorrection())
        {
            return this->cyclicAMIPatch().interpolate
            (
                nbrPatch.jump(),
                Field<Type>(this->size(), Zero)
            );
        }
        else
        {
            return this->cyclicAMIPatch().interpolate(nbrPatch.jump());
        }
    }
}


template<class Type>
void Foam::fixedJumpAMIFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpCyclicAMIFvPatchField<Type>::autoMap(m);
    jump_.autoMap(m);
}


template<class Type>
void Foam::fixedJumpAMIFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    jumpCyclicAMIFvPatchField<Type>::rmap(ptf, addr);

    const fixedJumpAMIFvPatchField<Type>& tiptf =
        refCast<const fixedJumpAMIFvPatchField<Type>>(ptf);
    jump_.rmap(tiptf.jump_, addr);
}


template<class Type>
void Foam::fixedJumpAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("patchType") << this->interfaceFieldType()
        << token::END_STATEMENT << nl;

    if (this->cyclicAMIPatch().owner())
    {
        jump_.writeEntry("jump", os);
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
