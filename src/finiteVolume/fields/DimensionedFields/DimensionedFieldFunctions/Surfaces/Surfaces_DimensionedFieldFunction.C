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

#include "Surfaces_DimensionedFieldFunction.H"
#include "searchableSurfacesInsideFraction.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Surfaces<DimensionedFieldType>::
readSurfacesAndValues()
{
    const objectRegistry& db = this->field_.mesh().db();

    surfacesPtr_.set
    (
        new searchableSurfaceList
        (
            IOobject
            (
                searchableSurfaceList::typeName,
                db.time().constant(),
                searchableSurface::geometryDir(db.time()),
                db,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            surfacesDict_,
            false
        )
    );

    values_.setSize(surfacesPtr_().size());

    forAll(surfacesPtr_(), i)
    {
        values_.set
        (
            i,
            new dimensioned<Type>
            (
                "value",
                this->field_.dimensions(),
                surfacesDict_.subDict(surfacesPtr_()[i].name())
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::Surfaces<DimensionedFieldType>::Surfaces
(
    const dictionary& dict,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dict, field),
    defaultValue_("defaultValue", this->field_.dimensions(), dict),
    surfacesDict_(dict.subDict("surfaces"))
{
    readSurfacesAndValues();
}


template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::Surfaces<DimensionedFieldType>::Surfaces
(
    const Surfaces& dff,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dff, field),
    defaultValue_(dff.defaultValue_),
    surfacesDict_(dff.surfacesDict_)
{
    readSurfacesAndValues();
}


template<class DimensionedFieldType>
Foam::autoPtr
<
    Foam::DimensionedFieldFunction<DimensionedFieldType>
> Foam::DimensionedFieldFunctions::Surfaces<DimensionedFieldType>::clone
(
    DimensionedFieldType& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
    (
        new Surfaces<DimensionedFieldType>(*this, field)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Surfaces<DimensionedFieldType>::evaluate()
{
    this->field_ = defaultValue_;

    forAll(surfacesPtr_(), i)
    {
        this->field_ +=
            DimensionedField
            <
                scalar,
                typename DimensionedFieldType::GeoMesh_
            >::New
            (
                IOobject::groupName("alpha", surfacesPtr_()[i].name()),
                this->field_.mesh(),
                dimless,
                searchableSurfaces::insideFraction
                (
                    surfacesPtr_()[i],
                    this->field_.mesh().poly()
                )
            )
           *(values_[i] - this->field_);
    }
}


template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Surfaces<DimensionedFieldType>::write
(
    Ostream& os
) const
{
    writeEntry(os, defaultValue_);
    writeEntry(os, surfacesDict_.dictName(), surfacesDict_);
}


// ************************************************************************* //
