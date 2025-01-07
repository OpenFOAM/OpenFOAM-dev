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


template<class AlphaRhoFieldType>
Foam::word Foam::LagrangianModel::fieldName
(
    const AlphaRhoFieldType& alphaRhoField
)
{
    const word group = alphaRhoField.group();
    const word member = alphaRhoField.member();
    const string::size_type i = member.find(':');
    return IOobject::groupName(member(i), group);
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

template<class Type, template<class> class PrimitiveField>
bool Foam::LagrangianModel::addsSupToField
(
    const LagrangianSubField<Type, PrimitiveField>& field
) const
{
    return addsSupToField(fieldName(field.name()));
}


// ************************************************************************* //

