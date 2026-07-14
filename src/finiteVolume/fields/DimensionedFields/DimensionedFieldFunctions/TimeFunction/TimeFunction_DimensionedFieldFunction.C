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

#include "TimeFunction_DimensionedFieldFunction.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::TimeFunction<DimensionedFieldType>::
TimeFunction
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
            field.time().userUnits(),
            field.dimensions(),
            dict
        )
    )
{}


template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::TimeFunction<DimensionedFieldType>::
TimeFunction
(
    const TimeFunction& dff,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dff, field),
    funcPtr_(dff.funcPtr_, false)
{}


template<class DimensionedFieldType>
Foam::autoPtr<Foam::DimensionedFieldFunction<DimensionedFieldType>>
Foam::DimensionedFieldFunctions::TimeFunction<DimensionedFieldType>::clone
(
    DimensionedFieldType& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
    (
        new TimeFunction<DimensionedFieldType>(*this, field)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::TimeFunction<DimensionedFieldType>::
evaluate()
{
    this->field_.primitiveFieldRef() =
        funcPtr_->value(this->field_.time().value());
}


template<class DimensionedFieldType>
bool Foam::DimensionedFieldFunctions::TimeFunction<DimensionedFieldType>::
update()
{
    evaluate();
    return true;
}


template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::TimeFunction<DimensionedFieldType>::write
(
    Ostream& os
) const
{
    writeEntry(os, dimensions::length, this->field_.dimensions(), funcPtr_());
}


// ************************************************************************* //
