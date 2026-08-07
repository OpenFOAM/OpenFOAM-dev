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
    Foam::expression::Reference
    Foam::expression::Pointer
    Foam::expression::Size
    Foam::expression::Mesh
    Foam::expression::MeshPointer
    Foam::expression::AllContainers
    Foam::expression::TypeContainers

Description
    Classes which provide access to the containers referenced within an
    expression, and to properties of those containers. Properties are
    identified by utility classes which also provide generic access to the
    property and (sometimes) a string representation thereof. The AllContainers
    and TypeContainers classes then apply these utility classes to expressions
    containing multiple such containers to access (potentially) multiple
    property values. AllContainers access all the container types in the
    expression, whilst TypeContainers only access those of a specific type.
    These classes should not need specialising further.

\*---------------------------------------------------------------------------*/

#ifndef expressionContainers_H
#define expressionContainers_H

#include "FixedList.H"
#include "IOobject.H"
#include "nil.H"

#include "expressionBase.H"
#include "expressionExpressionable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                          Class Reference Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct Reference
{
    using type = const Type&;

    static type get(const Type& c)
    {
        return c;
    }
};

/*---------------------------------------------------------------------------*\
                           Class Pointer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct Pointer
{
    using type = const Type*;

    static type get(const Type& c)
    {
        return &c;
    }
};

/*---------------------------------------------------------------------------*\
                            Class Size Declaration
\*---------------------------------------------------------------------------*/

struct Size
{
    using type = label;

    static const char* name()
    {
        return "size";
    }

    static const char* names()
    {
        return "sizes";
    }

    template<class Container>
    static type get(const Container& c)
    {
        return c.size();
    }
};

/*---------------------------------------------------------------------------*\
                            Class IO Declaration
\*---------------------------------------------------------------------------*/

struct IO
{
    using type = const IOobject&;

    template<class Container>
    static type get(const Container& c)
    {
        return c.io();
    }
};

/*---------------------------------------------------------------------------*\
                            Class Mesh Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct Mesh
{
    using type = const Type&;

    template<class Container>
    static type get(const Container& c)
    {
        return c.mesh();
    }
};

/*---------------------------------------------------------------------------*\
                         Class MeshPointer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct MeshPointer
{
    using type = const Type*;

    static const char* name()
    {
        return "mesh";
    }

    static const char* names()
    {
        return "meshes";
    }

    template<class Container>
    static type get(const Container& c)
    {
        return &c.mesh();
    }
};

/*---------------------------------------------------------------------------*\
                            Class Patch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct Patch
{
    using type = const Type&;

    template<class Container>
    static type get(const Container& c)
    {
        return c.patch();
    }
};

/*---------------------------------------------------------------------------*\
                        Class AllContainers Declaration
\*---------------------------------------------------------------------------*/

template<class Expression, class = void>
struct AllContainers;

template<class Expression>
struct AllContainers<Expression, EnableIfNotExpressionable<Expression>>
{
    static const unsigned n = 0;

    using type = nil;

    template<class Property, unsigned, unsigned N>
    static void get
    (
        const Expression& e,
        FixedList<typename Property::type, N>& result
    )
    {}
};

template<class Expression>
struct AllContainers<Expression, EnableIfExpressionable<Expression>>
{
    static const unsigned n = 1;

    using type = Expression;

    template<class Property, unsigned>
    static typename Property::type get(const Expression& e)
    {
        return Property::get(e);
    }

    template<class Property, unsigned I, unsigned N>
    static void get
    (
        const Expression& e,
        FixedList<typename Property::type, N>& result
    )
    {
        result.template elmt<I>() = Property::get(e);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Get a property from the first container in an expression
template<class Property, class Expression>
typename Property::type getAll(const Expression& e)
{
    typedef AllContainers<Expression> allContainers;

    return allContainers::template get<Property, 0>(e);
}

//- Get a list of the properties from all the containers in an expression
template<class Property, class Expression>
FixedList<typename Property::type, AllContainers<Expression>::n>
getAlls(const Expression& e)
{
    typedef AllContainers<Expression> allContainers;

    FixedList<typename Property::type, allContainers::n> result;

    allContainers::template get<Property, 0>(e, result);

    return result;
}

/*---------------------------------------------------------------------------*\
                        Class TypeContainers Declaration
\*---------------------------------------------------------------------------*/

template<class Expression, class Type, class = void>
struct TypeContainers;

template<class Expression, class Type>
struct TypeContainers<Expression, Type, EnableIfNotBaseOf<Type, Expression>>
{
    static const unsigned n = 0;

    template<class Property, unsigned, unsigned N>
    static void get
    (
        const Expression& e,
        FixedList<typename Property::type, N>& result
    )
    {}
};

template<class Expression, class Type>
struct TypeContainers<Expression, Type, EnableIfBaseOf<Type, Expression>>
{
    static const unsigned n = 1;

    template<class Property, unsigned>
    static typename Property::type get(const Expression& e)
    {
        return Property::get(e);
    }

    template<class Property, unsigned I, unsigned N>
    static void get
    (
        const Expression& e,
        FixedList<typename Property::type, N>& result
    )
    {
        result.template elmt<I>() = Property::get(e);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Get a property from the first container of a specific type in an expression
template
<
    class Container,
    class Property = Reference<Container>,
    class Expression
>
typename Property::type getType(const Expression& e)
{
    typedef TypeContainers<Expression, Container> typeContainers;

    return typeContainers::template get<Property, 0>(e);
}

//- Get a list of the properties from all the containers of a specific type in
//  an expression
template
<
    class Container,
    class Property,
    class Expression
>
FixedList<typename Property::type, TypeContainers<Expression, Container>::n>
getTypes(const Expression& e)
{
    typedef TypeContainers<Expression, Container> typeContainers;

    FixedList<typename Property::type, typeContainers::n> result;

    typeContainers::template get<Property, 0>(e, result);

    return result;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
