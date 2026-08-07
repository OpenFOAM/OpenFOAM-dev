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

Class
    Foam::expression::Operate

Description
    Class which provides the ability to access and evaluate an expression
    element-by-element. This must be specialised for expressionable containers
    in order for them to support field algebra.

\*---------------------------------------------------------------------------*/

#ifndef expressionOperate_H
#define expressionOperate_H

#include "expressionExpressionable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                          Class Operate Declaration
\*---------------------------------------------------------------------------*/

template<class Expression, class = Arguments<>, class = void>
struct Operate;

template<class Expression, class ... Args>
struct Operate
<
    Expression,
    Arguments<Args ...>,
    EnableIfNotExpressionable<Expression>
>
{
    using type = const Expression&;

    static type access(const Expression& e, const Args& ... args)
    {
        return e;
    }

    static bool compare(const Expression& e, const Args& ... args)
    {
        return false;
    }

    static void advance(const Expression& e, Args& ... args)
    {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Access part of a given expression. E.g., given an index as an argument,
//  return an element. Or given a type representing the start or end, return
//  an iterator.
template<class Expression, class ... Args>
typename Operate<Expression, Arguments<Args ...>>::type
access(const Expression& e, const Args& ... args)
{
    return Operate<Expression, Arguments<Args ...>>::access(e, args ...);
}

//- Compare parts of a given expression
template<class Expression, class ... Args>
bool compare(const Expression& e, const Args& ... args)
{
    return Operate<Expression, Arguments<Args ...>>::compare(e, args ...);
}

//- Advance the parts of a given expression to the next element
template<class Expression, class ... Args>
void advance(const Expression& e, Args& ... args)
{
    return Operate<Expression, Arguments<Args ...>>::advance(e, args ...);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
