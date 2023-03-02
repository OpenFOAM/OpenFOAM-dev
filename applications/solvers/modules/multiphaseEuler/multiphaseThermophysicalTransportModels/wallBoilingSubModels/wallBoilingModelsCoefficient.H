/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

Namespace
    Foam::wallBoilingModels

\*---------------------------------------------------------------------------*/

#ifndef wallBoilingModelsCoefficient_H
#define wallBoilingModelsCoefficient_H

#include "dimensionedTypes.H"
#include "volFieldsFwd.H"
#include "UniformDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

template<class FieldType>
struct coefficient;

template<class Type>
struct coefficient<Field<Type>>
{
    template<class Type2>
    static const Type2& value(const dimensioned<Type2>& v)
    {
        return v.value();
    }

    template<class Type2>
    static const Type2& value(const UniformDimensionedField<Type2>& v)
    {
        return v.value();
    }
};

template<class Type>
struct coefficient<VolField<Type>>
{
    template<class Type2>
    static const dimensioned<Type2>& value(const dimensioned<Type2>& v)
    {
        return v;
    }

    template<class Type2>
    static const dimensioned<Type2>& value
    (
        const UniformDimensionedField<Type2>& v
    )
    {
        return v;
    }
};

} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
