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
#include "expressionBound.H"
#include "expressionBase.H"

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


/*---------------------------------------------------------------------------*\
                          Class OperateIterable Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class = Arguments<>>
struct OperateIterable;

//- Define access by index
template<class Container>
struct OperateIterable<Container, Arguments<label>>
{
    using type = const typename Container::value_type&;

    static type access(const Container& f, const label i)
    {
        return f[i];
    }
};

//- Define access to the start and end iterators
template<class Container, bool Index>
struct OperateIterable<Container, Arguments<Bound<Index>>>
{
    using type = typename Container::const_iterator;

    static type access(const Container& f, const Bound<false>&)
    {
        return f.begin();
    }

    static type access(const Container& f, const Bound<true>&)
    {
        return f.end();
    }
};

//- Define access by dereferencing an iterator. Also define comparison and
//  incrementation of iterators.
template<class Container, class ... Args>
struct OperateIterable
<
    Container,
    Arguments<typename Container::const_iterator, Args ...>
>
{
    using type = const typename Container::value_type&;

    static type access
    (
        const Container& f,
        const typename Container::const_iterator& fi
    )
    {
        return *fi;
    }

    static bool compare
    (
        const Container& f,
        const typename Container::const_iterator& fi,
        const typename Container::const_iterator& fj
    )
    {
        return fi == fj;
    }

    static void advance
    (
        const Container& f,
        typename Container::const_iterator& fi
    )
    {
        ++ fi;
    }
};

//- Define access to the base as the container itself. An iterable container
//  *is* typically the expressionable base class. It does not *have* an
//  expressionable base class.
template<class Container>
struct OperateIterable<Container, Arguments<Base>>
{
    using type = const Container&;

    static type access(const Container& f, const Base&)
    {
        return f;
    }
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
