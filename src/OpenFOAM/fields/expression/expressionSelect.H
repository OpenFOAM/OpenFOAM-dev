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
    Foam::expression::Select

Description
    Simple wrapper that modifies expression access so that it constructs and
    returns a (tmp) container of a given type, rather than producing a
    sub-expression. Used to construct from expressions, rather than just to
    assign to them. This will need overloading if the type in question doesn't
    have a basic construct-from-expression constructor, and/or if the type is
    abstract and some sort of default selection is needed.

\*---------------------------------------------------------------------------*/

#ifndef expressionSelect_H
#define expressionSelect_H

#include "expressionExpressionable.H"
#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                            Class Select Declaration
\*---------------------------------------------------------------------------*/

template<class Type, const bool Ptr, class = void>
struct Select;

template<class Type, bool Ptr>
struct Select<Type, Ptr, EnableIfExpressionable<Type>>
{
    using type = tmp<Type>;

    template<class Expression>
    static tmp<Type> New(const Expression& e)
    {
        return tmp<Type>(new Type(e));
    }
};

template<class Type>
struct Select<Type, true, EnableIfNotExpressionable<Type>>
{
    using type = autoPtr<Type>;

    template<class Expression>
    static type New(const Expression& e)
    {
        return autoPtr<Type>(new Type(e));
    }
};

template<class Type>
struct Select<Type, false, EnableIfNotExpressionable<Type>>
{
    using type = Type;

    template<class Expression>
    static type New(const Expression& e)
    {
        return Type(e);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a new object constructed from the result of an expression access
template<class Type, const bool Ptr, class Expression, class ... Args>
typename Select<Type, Ptr>::type New(const Expression& e, const Args& ... args)
{
    return Select<Type, Ptr>::template New<Expression>(e, args ...);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
