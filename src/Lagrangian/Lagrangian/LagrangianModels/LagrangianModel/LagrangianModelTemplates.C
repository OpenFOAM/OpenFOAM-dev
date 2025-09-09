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

#include "LagrangianModel.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class AlphaRhoFieldType, class ... AlphaRhoFieldTypes>
Foam::word Foam::LagrangianModel::fieldName
(
    const AlphaRhoFieldType& alphaRhoField,
    const AlphaRhoFieldTypes& ... alphaRhoFields
)
{
    return fieldName(alphaRhoFields ...);
}


template
<
    class Type,
    class GeoMesh,
    template<class> class PrimitiveField
>
const Foam::word& Foam::LagrangianModel::fieldName
(
    const GeometricField<Type, GeoMesh, PrimitiveField>& field
)
{
    return field.name();
}


template
<
    class Type,
    class GeoMesh,
    template<class> class PrimitiveField
>
const Foam::word& Foam::LagrangianModel::fieldName
(
    const DimensionedField<Type, GeoMesh, PrimitiveField>& field
)
{
    return field.name();
}


template
<
    class Type,
    template<class> class PrimitiveField
>
Foam::word Foam::LagrangianModel::fieldName
(
    const DimensionedField<Type, LagrangianSubMesh, PrimitiveField>& field
)
{
    return field.mesh().complete(field.name());
}


template<class AlphaRhoFieldType, class ... AlphaRhoFieldTypes>
Foam::word Foam::LagrangianModel::fieldsName
(
    const AlphaRhoFieldType& alphaRhoField,
    const AlphaRhoFieldTypes& ... alphaRhoFields
)
{
    return fieldName(alphaRhoField) + '*' + fieldsName(alphaRhoFields ...);
}


template<class AlphaRhoFieldType>
Foam::word Foam::LagrangianModel::fieldsName
(
    const AlphaRhoFieldType& alphaRhoField
)
{
    return fieldName(alphaRhoField);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PrimitiveField,
    template<class> class PrimitiveEqnField
>
bool Foam::LagrangianModel::addsSupToField
(
    const LagrangianSubField<Type, PrimitiveField>& field,
    const LagrangianSubField<Type, PrimitiveEqnField>& eqnField
) const
{
    return addsSupToField(fieldName(field.name()), fieldName(eqnField.name()));
}


// ************************************************************************* //

