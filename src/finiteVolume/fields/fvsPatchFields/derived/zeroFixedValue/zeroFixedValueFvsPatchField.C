/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "zeroFixedValueFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroFixedValueFvsPatchField<Type>::zeroFixedValueFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fixedValueFvsPatchField<Type>(p, iF)
{
    this->operator==(Zero);
}


template<class Type>
Foam::zeroFixedValueFvsPatchField<Type>::zeroFixedValueFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvsPatchField<Type>(p, iF)
{
    this->operator==(Zero);
}


template<class Type>
Foam::zeroFixedValueFvsPatchField<Type>::zeroFixedValueFvsPatchField
(
    const zeroFixedValueFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvsPatchField<Type>(ptf, p, iF, mapper)
{
    this->operator==(Zero);
}


template<class Type>
Foam::zeroFixedValueFvsPatchField<Type>::zeroFixedValueFvsPatchField
(
    const zeroFixedValueFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fixedValueFvsPatchField<Type>(ptf, iF)
{
    this->operator==(Zero);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::zeroFixedValueFvsPatchField<Type>::map
(
    const fvsPatchField<Type>& ptf,
    const fieldMapper& mapper
)
{
    this->resize(this->patch().size(), Zero);
}


template<class Type>
void Foam::zeroFixedValueFvsPatchField<Type>::write(Ostream& os) const
{
    fvsPatchField<Type>::write(os);
}


// ************************************************************************* //
