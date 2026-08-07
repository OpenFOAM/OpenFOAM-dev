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

Namespace
    Foam::expression

Description
    Specialise expression::Operate for dimensionSet so that non-dimensioned
    types access the dimensionless set and dimensioned types access their
    dimensions.

    !!! This could/should ultimately be put directly into dimensionSet.H

\*---------------------------------------------------------------------------*/

#ifndef expressionDimensionSet_H
#define expressionDimensionSet_H

#include "dimensionSet.H"

#include "expressionOperate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
struct Operate
<
    Type,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<Type>
>
{
    using type = const dimensionSet&;

    static type access(const Type& t, const dimensionSet&)
    {
        return dimless;
    }
};

template<class Type>
struct Operate
<
    dimensioned<Type>,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<dimensioned<Type>>
>
{
    using type = const dimensionSet&;

    static type access(const dimensioned<Type>& dt, const dimensionSet&)
    {
        return dt.dimensions();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
