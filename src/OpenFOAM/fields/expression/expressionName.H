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

Classes
    Foam::expression::Operation
    Foam::expression::Name

Description
    Classes for constructing string representations of expressions.

    The Foam::expression::Operation trait provides the base operator type and
    its precedence for a given expression.

    The Foam::expression::Name class constructs the actual name. These are
    specialised for unary and binary operations in order to construct compound
    names. Foam::expression::Name implements a streamed output for any
    non-expression and a type and a pointer address for any expressionable
    container. This should be specialised for any container type for which a
    user-readable name is available.

\*---------------------------------------------------------------------------*/

#ifndef expressionName_H
#define expressionName_H

#include "expressionExpressionable.H"
#include "OStringStream.H"
#include "typeName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                         Class Operation Declaration
\*---------------------------------------------------------------------------*/

template<class Expression>
struct Operation
{
    using type = void;

    static const signed char precedence =
        std::numeric_limits<signed char>::min();
};

/*---------------------------------------------------------------------------*\
                           Class Name Declaration
\*---------------------------------------------------------------------------*/

template<class Expression, class = void>
struct Name;

template<class Expression>
struct Name<Expression, EnableIfNotExpressionable<Expression>>
{
    static word value(const Expression& e)
    {
        OStringStream oss;
        oss << e;
        return oss.str();
    }
};

template<class Expression>
struct Name<Expression, EnableIfExpressionable<Expression>>
{
    static word value(const Expression& f)
    {
        std::ostringstream oss;
        oss << std::hex << reinterpret_cast<uintptr_t>(&(*f.begin()));
        return "[" + typeName<Expression, word>() + ":" + oss.str() + "]";
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Function to get the name of an expression. More convenient than using the
//  classes above directly, as the expression type is deduced.
template<class Expression>
word name(const Expression& e)
{
    return Name<Expression>::value(e);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
