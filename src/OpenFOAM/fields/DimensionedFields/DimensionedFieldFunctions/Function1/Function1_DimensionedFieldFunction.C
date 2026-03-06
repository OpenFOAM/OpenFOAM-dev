/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "Function1_DimensionedFieldFunction.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::Function1<DimensionedFieldType>::
Function1
(
    const dictionary& dict,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dict, field),
    funcPtr_
    (
        Foam::Function1<typename DimensionedFieldType::Type_>::New
        (
            "function",
            dimLength,
            field.dimensions(),
            dict
        )
    ),
    direction_(normalised(dict.lookup<vector>("direction")))
{
    evaluate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Function1<DimensionedFieldType>::
evaluate()
{
    this->field_.primitiveFieldRef() =
        funcPtr_->value(direction_ & this->field_.mesh().C());
}


template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Function1<DimensionedFieldType>::write
(
    Ostream& os
) const
{
    writeEntry(os, dimLength, this->field_.dimensions(), funcPtr_());
    writeEntry(os, "direction", direction_);
}


// ************************************************************************* //
