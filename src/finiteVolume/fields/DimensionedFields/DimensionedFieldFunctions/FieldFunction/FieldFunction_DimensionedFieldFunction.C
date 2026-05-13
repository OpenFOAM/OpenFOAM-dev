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

#include "FieldFunction_DimensionedFieldFunction.H"
#include "GeometricField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::FieldFunction<DimensionedFieldType>::
FieldFunction
(
    const dictionary& dict,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dict, field),
    fieldName_(dict.lookup<word>("field")),
    funcPtr_(nullptr),
    dictPtr_(new dictionary(dict))
{}


template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::FieldFunction<DimensionedFieldType>::
FieldFunction
(
    const FieldFunction& dff,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dff, field),
    fieldName_(dff.fieldName_),
    funcPtr_(dff.funcPtr_, false),
    dictPtr_(nullptr)
{}


template<class DimensionedFieldType>
Foam::autoPtr<Foam::DimensionedFieldFunction<DimensionedFieldType>>
Foam::DimensionedFieldFunctions::FieldFunction<DimensionedFieldType>::clone
(
    DimensionedFieldType& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
    (
        new FieldFunction<DimensionedFieldType>(*this, field)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::FieldFunction<DimensionedFieldType>::
evaluate()
{
    typedef
        DimensionedField<scalar, typename DimensionedFieldType::GeoMesh_>
        scalarFieldType;

    tmp<scalarFieldType> tx =
        this->field_.mesh().template lookupField<scalar>(fieldName_);

    if (dictPtr_.valid())
    {
        funcPtr_.set
        (
            Foam::Function1<typename DimensionedFieldType::Type_>::New
            (
                "function",
                tx().dimensions(),
                this->field_.dimensions(),
                dictPtr_()
            ).ptr()
        );

        dictPtr_.clear();
    }

    this->field_.primitiveFieldRef() = funcPtr_->value(tx().primitiveField());
}


template<class DimensionedFieldType>
bool Foam::DimensionedFieldFunctions::FieldFunction<DimensionedFieldType>::
update()
{
    evaluate();
    return true;
}


template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::FieldFunction<DimensionedFieldType>::write
(
    Ostream& os
) const
{
    typedef
        DimensionedField<scalar, typename DimensionedFieldType::GeoMesh_>
        scalarFieldType;

    writeEntry(os, "field", fieldName_);

    if (dictPtr_.valid())
    {
        writeEntry(os, "function", dictPtr_());
    }
    else
    {
        tmp<scalarFieldType> tx =
            this->field_.mesh().template lookupField<scalar>(fieldName_);

        writeEntry
        (
            os,
            tx().dimensions(),
            this->field_.dimensions(),
            funcPtr_()
        );
    }
}


// ************************************************************************* //
