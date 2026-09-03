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
    Specialise expression::Operate for dimensioned<Type> so that it returns its
    value when an expression::Value argument is given. Also specialise
    expression::Name for dimensioned<Type> so that just the name is included
    into an expression string, rather than the full string representation
    (i.e., name, dimensions and value).

    !!! This could/should ultimately be put directly into dimensionedType.H

\*---------------------------------------------------------------------------*/

#ifndef expressionDimensionedType_H
#define expressionDimensionedType_H

#include "dimensionedType.H"

#include "expressionOperate.H"
#include "expressionName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                            Class Value Declaration
\*---------------------------------------------------------------------------*/

struct Value
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class ... Args>
struct Operate
<
    dimensioned<Type>,
    Arguments<Value, Args ...>,
    EnableIfNotExpressionable<dimensioned<Type>>
>
{
    using type = const Type&;

    static type access
    (
        const dimensioned<Type>& dt,
        const Value&,
        const Args& ... args
    )
    {
        return dt.value();
    }
};

template<class Type>
struct Name
<
    dimensioned<Type>,
    EnableIfNotExpressionable<dimensioned<Type>>
>
{
    static word value(const dimensioned<Type>& dt)
    {
        return dt.name();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
