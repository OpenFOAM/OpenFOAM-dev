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

#include "UFunctionalDimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::UFunctionalDimensionedField<Type, GeoMesh>::UFunctionalDimensionedField
(
    const word& name,
    const word& funcName,
    const GeoMesh& mesh,
    const dimensionSet& dimensions,
    Field<Type>& field,
    const dictionary& dict
)
:
    SlicedDimensionedField<Type, GeoMesh>
    (
        IOobject
        (
            word(name, '_', funcName),
            mesh.time().name(),
            mesh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensions,
        field
    ),
    field_(field),
    funcName_(funcName),
    funcPtr_
    (
        dict.isDict(funcName)
      ? DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>::New
        (
            dict.subDict(funcName),
            *this
        )
      : autoPtr<DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>>
        (
            nullptr
        )
    )
{
    if (funcPtr_.valid())
    {
        if (mesh.time().completeCase())
        {
            funcPtr_->evaluate();
        }
    }
    else
    {
        field_ = Field<Type>(funcName_, dimensions, dict, mesh.size());
    }
}


template<class Type, class GeoMesh>
Foam::UFunctionalDimensionedField<Type, GeoMesh>::UFunctionalDimensionedField
(
    const word& name,
    const word& funcName,
    const GeoMesh& mesh,
    const dimensionSet& dimensions,
    Field<Type>& field,
    const dictionary& dict,
    const Type& defaultValue
)
:
    SlicedDimensionedField<Type, GeoMesh>
    (
        IOobject
        (
            word(name, '_', funcName),
            mesh.time().name(),
            mesh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensions,
        field
    ),
    field_(field),
    funcName_(funcName),
    funcPtr_
    (
        dict.isDict(funcName)
      ? DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>::New
        (
            dict.subDict(funcName),
            *this
        )
      : autoPtr<DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>>
        (
            nullptr
        )
    )
{
    if (funcPtr_.valid())
    {
        if (mesh.time().completeCase())
        {
            funcPtr_->evaluate();
        }
    }
    else if (dict.found(funcName))
    {
        field_ = Field<Type>(funcName_, dimensions, dict, mesh.size());
        SlicedDimensionedField<Type, GeoMesh>::reset(field_);
    }
    else
    {
        field_ = Field<Type>(mesh.size(), defaultValue);
        SlicedDimensionedField<Type, GeoMesh>::reset(field_);
    }
}


template<class Type, class GeoMesh>
Foam::UFunctionalDimensionedField<Type, GeoMesh>::UFunctionalDimensionedField
(
    const UFunctionalDimensionedField<Type, GeoMesh>& udff,
    const GeoMesh& mesh,
    Field<Type>& field
)
:
    SlicedDimensionedField<Type, GeoMesh>
    (
        udff,
        mesh,
        udff.dimensions(),
        field
    ),
    field_(field),
    funcName_(udff.funcName_),
    funcPtr_
    (
        udff.funcPtr_.valid()
      ? udff.funcPtr_->clone(*this)
      : autoPtr<DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>>
        (
            nullptr
        )
    )
{}


template<class Type, class GeoMesh>
Foam::UFunctionalDimensionedField<Type, GeoMesh>::UFunctionalDimensionedField
(
    const UFunctionalDimensionedField<Type, GeoMesh>& udff,
    Field<Type>& field
)
:
    SlicedDimensionedField<Type, GeoMesh>
    (
        udff,
        udff.mesh(),
        udff.dimensions(),
        field
    ),
    field_(field),
    funcName_(udff.funcName_),
    funcPtr_
    (
        udff.funcPtr_.valid()
      ? udff.funcPtr_->clone(*this)
      : autoPtr<DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>>
        (
            nullptr
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
void Foam::UFunctionalDimensionedField<Type, GeoMesh>::map(const bool evaluate)
{
    SlicedDimensionedField<Type, GeoMesh>::reset(field_);
    if (evaluate && funcPtr_.valid())
    {
        funcPtr_->reset();
    }
}


template<class Type, class GeoMesh>
void Foam::UFunctionalDimensionedField<Type, GeoMesh>::reset()
{
    SlicedDimensionedField<Type, GeoMesh>::reset(field_);

    if (funcPtr_.valid())
    {
        funcPtr_->reset();
    }
}


template<class Type, class GeoMesh>
bool Foam::UFunctionalDimensionedField<Type, GeoMesh>::update()
{
    if (funcPtr_.valid())
    {
        return funcPtr_->update();
    }
    else
    {
        return false;
    }
}


template<class Type, class GeoMesh>
void Foam::UFunctionalDimensionedField<Type, GeoMesh>::write(Ostream& os) const
{
    if (funcPtr_.valid())
    {
        funcPtr_->write(os);
    }
}


template<class Type, class GeoMesh>
void Foam::writeEntry
(
    Ostream& os,
    const UFunctionalDimensionedField<Type, GeoMesh>& udff
)
{
    if (udff.funcPtr_.valid())
    {
        writeEntry(os, udff.funcName_, *udff.funcPtr_);
    }
    else
    {
        writeEntry(os, udff.funcName_, udff.field_);
    }
}


// ************************************************************************* //
