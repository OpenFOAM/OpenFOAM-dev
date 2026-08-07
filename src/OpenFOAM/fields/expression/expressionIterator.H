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
    Foam::expression::Bounds
    Foam::expression::IteratorBase
    Foam::expression::Iterator
    Foam::expression::IteratorNew

Description
    Iterator class which facilitates iteration through an expression

\*---------------------------------------------------------------------------*/

#ifndef expressionIterator_H
#define expressionIterator_H

#include "expressionOperate.H"
#include "expressionSelect.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                            Class Bound Declaration
\*---------------------------------------------------------------------------*/

template<bool>
struct Bound
{};

/*---------------------------------------------------------------------------*\
                           Class IteratorBase Declaration
\*---------------------------------------------------------------------------*/

template<template<class> class Iterator, class Expression>
struct IteratorBase
{
protected:

    const Expression& e_;

    using eiType = typename Operate<Expression, Arguments<Bound<false>>>::type;

    eiType ei_;

public:

    template<bool Index>
    inline IteratorBase
    (
        const Expression& e,
        const Bound<Index>&
    )
    :
        e_(e),
        ei_
        (
            Operate<Expression, Arguments<Bound<Index>>>::access
            (
                e,
                Bound<Index>()
            )
        )
    {}

    inline Iterator<Expression>& operator++()
    {
        expression::advance(e_, ei_);
        return static_cast<Iterator<Expression>&>(*this);
    }

    inline Iterator<Expression> operator++(int)
    {
        Iterator<Expression> ei0(e_, ei_);
        expression::advance(e_, ei_);
        return ei0;
    }

    inline bool operator==(const IteratorBase<Iterator, Expression>& iter) const
    {
        return expression::compare(e_, ei_, iter.ei_);
    }

    inline bool operator!=(const IteratorBase<Iterator, Expression>& iter) const
    {
        return !expression::compare(e_, ei_, iter.ei_);
    }
};

/*---------------------------------------------------------------------------*\
                           Class Iterator Declaration
\*---------------------------------------------------------------------------*/

template<class Expression>
struct Iterator
:
    public IteratorBase<Iterator, Expression>
{
    using IteratorBase<Iterator, Expression>::IteratorBase;

    typedef typename Operate
    <
        Expression,
        Arguments<typename IteratorBase<Iterator, Expression>::eiType>
    >::type value_type;

    value_type operator*() const
    {
        return expression::access(this->e_, this->ei_);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return an iterator pointing to the start of the expression
template<class Expression>
Iterator<Expression> begin(const Expression& e)
{
    return Iterator<Expression>(e, Bound<false>());
}

//- Return an iterator pointing to the end of the expression
template<class Expression>
Iterator<Expression> end(const Expression& e)
{
    return Iterator<Expression>(e, Bound<true>());
}

/*---------------------------------------------------------------------------*\
                          Class IteratorNew Declaration
\*---------------------------------------------------------------------------*/

template<class Type, bool Ptr, class Expression>
struct IteratorNew
:
    public IteratorBase<Iterator, Expression>
{
    using IteratorBase<Iterator, Expression>::IteratorBase;

    typedef typename Select<Type, Ptr>::type value_type;

    typename Select<Type, Ptr>::type operator*() const
    {
        return
            expression::New<Type, Ptr>
            (
                expression::access(this->e_, this->ei_)
            );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return an iterator pointing to the start of the expression
template<class Type, bool Ptr, class Expression>
IteratorNew<Type, Ptr, Expression> beginNew(const Expression& e)
{
    return IteratorNew<Type, Ptr, Expression>(e, Bound<false>());
}

//- Return an iterator pointing to the end of the expression
template<class Type, bool Ptr, class Expression>
IteratorNew<Type, Ptr, Expression> endNew(const Expression& e)
{
    return IteratorNew<Type, Ptr, Expression>(e, Bound<true>());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Expression>
struct std::iterator_traits
<
    Foam::expression::Iterator<Expression>
>
{
    typedef Foam::label difference_type;
    typedef typename
        Foam::expression::Iterator<Expression>::value_type
        value_type;
    typedef void pointer_type;
    typedef void reference;
    typedef std::forward_iterator_tag iterator_category;
};

template<class Type, bool Ptr, class Expression>
struct std::iterator_traits
<
    Foam::expression::IteratorNew<Type, Ptr, Expression>
>
{
    typedef Foam::label difference_type;
    typedef typename
        Foam::expression::IteratorNew<Type, Ptr, Expression>::value_type
        value_type;
    typedef void pointer_type;
    typedef void reference;
    typedef std::forward_iterator_tag iterator_category;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
