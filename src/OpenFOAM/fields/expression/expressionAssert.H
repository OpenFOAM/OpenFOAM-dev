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
    Functions which check the consistency of the contents of expressions

\*---------------------------------------------------------------------------*/

#ifndef expressionAssert_H
#define expressionAssert_H

#include "expressionContainers.H"
#include "expressionName.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Assert that all containers in an expression have the same property value
template<class Property, class Expression>
void assertSameAllContainerProperty(const Expression& e)
{
    typedef AllContainers<Expression> allContainers;

    FixedList<typename Property::type, allContainers::n> ps =
        getAll<Property>(e);

    for (unsigned i = 1; i < allContainers::n; ++ i)
    {
        if (ps[0] != ps[i])
        {
            FatalErrorInFunction
                << "expression '" << name(e)
                << "' of type '" << typeName<Expression>()
                << "' contains containers with different "
                << Property::names() << exit(FatalError);
        }
    }
}

//- Assert that all containers in an expression have the same property value as
//  that of the given container
template<class Property, class Container, class Expression>
void assertSameAllContainerProperty
(
    const Container& c,
    const Expression& e
)
{
    assertSameAllContainerProperty<Property>(e);

    typedef AllContainers<Expression> allContainers;

    FixedList<typename Property::type, allContainers::n> ps =
        getAll<Property>(e);

    if (allContainers::n)
    {
        if (Property::get(c) != ps[0])
        {
            FatalErrorInFunction
                << "expression '" << name(e)
                << "' of type '" << typeName<Expression>()
                << "' has a different " << Property::name()
                << " to container '" << name(c)
                << "' of type '" << typeName<Container>()
                << "'" << exit(FatalError);
        }
    }
}

//- Assert that all containers of a type in an expression have the same
//  property value
template<class Container, class Property, class Expression>
void assertSameTypeContainerProperty(const Expression& e)
{
    typedef TypeContainers<Expression, Container> typeContainers;

    FixedList<typename Property::type, typeContainers::n> ps =
        getTypes<Container, Property>(e);

    for (unsigned i = 1; i < typeContainers::n; ++ i)
    {
        if (ps[0] != ps[i])
        {
            FatalErrorInFunction
                << "expression '" << name(e)
                << "' of type '" << typeName<Expression>()
                << "' contains containers of type '"
                << typeName<Container>() << "' with different "
                << Property::names() << exit(FatalError);
        }
    }
}

//- Assert that all containers of a type in an expression have the same
//  property value as that of the given container
template<class Property, class Container, class Expression>
void assertSameTypeContainerProperty
(
    const Container& c,
    const Expression& e
)
{
    assertSameTypeContainerProperty<Container, Property>(e);

    typedef TypeContainers<Expression, Container> typeContainers;

    FixedList<typename Property::type, typeContainers::n> ps =
        getTypes<Container, Property>(e);

    if (typeContainers::n)
    {
        if (Property::get(c) != ps[0])
        {
            FatalErrorInFunction
                << "expression '" << name(e)
                << "' of type '" << typeName<Expression>()
                << "' has a different " << Property::name()
                << " to container '" << name(c)
                << "' of type '" << typeName<Container>()
                << "'" << exit(FatalError);
        }
    }
}

//- Assert that two containers have the same property value
template<class Property, class Container>
void assertSameTypeContainerProperty
(
    const Container& a,
    const Container& b
)
{
    if (Property::get(a) != Property::get(b))
    {
        FatalErrorInFunction
            << "containers '" << name(a)
            << "' and '" << name(b)
            << "' of type '" << typeName<Container>()
            << "' have different " << Property::names()
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
