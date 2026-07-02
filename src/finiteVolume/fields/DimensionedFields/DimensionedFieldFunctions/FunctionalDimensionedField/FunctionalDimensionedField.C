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

#include "FunctionalDimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::FunctionalDimensionedField<Type, GeoMesh>::FunctionalDimensionedField
(
    const word& name,
    const word& funcName,
    const GeoMesh& mesh,
    const dimensionSet& dimensions,
    const dictionary& dict
)
:
    DimensionedField<Type, GeoMesh>
    (
        IOobject
        (
            name + '_' + funcName,
            mesh.time().name(),
            mesh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensions,
        false
    ),
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
        this->primitiveFieldRef() =
            Field<Type>(funcName_, dimensions, dict, mesh.size());
    }
}


template<class Type, class GeoMesh>
Foam::FunctionalDimensionedField<Type, GeoMesh>::FunctionalDimensionedField
(
    const word& name,
    const word& funcName,
    const GeoMesh& mesh,
    const dimensionSet& dimensions,
    const dictionary& dict,
    const Type& defaultValue
)
:
    DimensionedField<Type, GeoMesh>
    (
        IOobject
        (
            name + '_' + funcName,
            mesh.time().name(),
            mesh.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensions,
        defaultValue
    ),
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
        this->primitiveFieldRef() =
            Field<Type>(funcName_, dimensions, dict, mesh.size());
    }
    else
    {
        this->primitiveFieldRef() =
            Field<Type>(mesh.size(), defaultValue);
    }
}


template<class Type, class GeoMesh>
Foam::FunctionalDimensionedField<Type, GeoMesh>::FunctionalDimensionedField
(
    const FunctionalDimensionedField<Type, GeoMesh>& udff,
    const GeoMesh& mesh
)
:
    DimensionedField<Type, GeoMesh>
    (
        udff,
        mesh,
        udff.dimensions(),
        false
    ),
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
Foam::FunctionalDimensionedField<Type, GeoMesh>::FunctionalDimensionedField
(
    const FunctionalDimensionedField<Type, GeoMesh>& udff
)
:
    DimensionedField<Type, GeoMesh>(udff),
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
void Foam::FunctionalDimensionedField<Type, GeoMesh>::map(const bool evaluate)
{
    InfoInFunction << endl;
    this->setSize(this->mesh().size());

    if (evaluate && funcPtr_.valid())
    {
        funcPtr_->reset();
    }
}


template<class Type, class GeoMesh>
void Foam::FunctionalDimensionedField<Type, GeoMesh>::reset()
{
    this->setSize(this->mesh().size());

    if (funcPtr_.valid())
    {
        funcPtr_->reset();
    }
}


template<class Type, class GeoMesh>
bool Foam::FunctionalDimensionedField<Type, GeoMesh>::update()
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
void Foam::FunctionalDimensionedField<Type, GeoMesh>::write(Ostream& os) const
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
    const FunctionalDimensionedField<Type, GeoMesh>& udff
)
{
    if (udff.funcPtr_.valid())
    {
        writeEntry(os, udff.funcName_, *udff.funcPtr_);
    }
    else
    {
        writeEntry(os, udff.funcName_, udff.primitiveField());
    }
}


// ************************************************************************* //
