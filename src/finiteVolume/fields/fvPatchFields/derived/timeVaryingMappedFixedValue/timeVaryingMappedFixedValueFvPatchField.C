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

#include "timeVaryingMappedFixedValueFvPatchField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    fieldMapper_(p, dict, iF.name())
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        // Note: we use evaluate() here to trigger updateCoeffs followed
        //       by re-setting of fvatchfield::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const timeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    fieldMapper_(ptf.fieldMapper_)
{}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const timeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    fieldMapper_(ptf.fieldMapper_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::map
(
    const fvPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchField<Type>::map(ptf, mapper);
    fieldMapper_.map
    (
        refCast<const timeVaryingMappedFixedValueFvPatchField<Type>>(ptf)
       .fieldMapper_,
        mapper
    );
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::reset
(
    const fvPatchField<Type>& ptf
)
{
    fixedValueFvPatchField<Type>::reset(ptf);
    fieldMapper_.reset
    (
        refCast<const timeVaryingMappedFixedValueFvPatchField<Type>>(ptf)
       .fieldMapper_
    );
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(fieldMapper_.map());

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    fieldMapper_.write(os);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
