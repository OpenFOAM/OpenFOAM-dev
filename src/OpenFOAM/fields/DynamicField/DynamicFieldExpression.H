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

\*---------------------------------------------------------------------------*/

#include "expressions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A DynamicField is expressionable
template<class Type>
struct Expressionable<DynamicField<Type>>
:
    public std::integral_constant
    <
        signed char,
        std::max(static_cast<signed char>(0), Expressionable<Type>::value)
      + static_cast<signed char>(1)
    >
{};

//- Access to a DynamicField by index is expressionable if the type is
template<class Type>
struct Expressionable<DynamicField<Type>, Arguments<label>>
:
    public Expressionable<Type>
{};

//- Access to a DynamicField by iterator is expressionable if the type is
template<class Type>
struct Expressionable
<
    DynamicField<Type>,
    Arguments<typename DynamicField<Type>::const_iterator>
>
:
    public Expressionable<Type>
{};

//- Define access to a value in a DynamicField by index
template<class Type>
struct Operate<DynamicField<Type>, Arguments<label>>
{
    using type = const Type&;

    static type access(const DynamicField<Type>& f, const label i)
    {
        return f[i];
    }
};

//- Define access to the start and end iterators of a DynamicField
template<class Type, bool Index>
struct Operate<DynamicField<Type>, Arguments<Bound<Index>>>
{
    using type = typename DynamicField<Type>::const_iterator;

    static type access(const DynamicField<Type>& f, const Bound<false>&)
    {
        return f.begin();
    }

    static type access(const DynamicField<Type>& f, const Bound<true>&)
    {
        return f.end();
    }
};

//- Define access to a value in a DynamicField by dereferencing an iterator.
//  Also define comparison and incrementation of iterators.
template<class Type, class ... Args>
struct Operate
<
    DynamicField<Type>,
    Arguments<typename DynamicField<Type>::const_iterator, Args ...>
>
{
    using type = const Type&;

    static type access
    (
        const DynamicField<Type>& f,
        const typename DynamicField<Type>::const_iterator& fi
    )
    {
        return *fi;
    }

    static bool compare
    (
        const DynamicField<Type>& f,
        const typename DynamicField<Type>::const_iterator& fi,
        const typename DynamicField<Type>::const_iterator& fj
    )
    {
        return fi == fj;
    }

    static void advance
    (
        const DynamicField<Type>& f,
        typename DynamicField<Type>::const_iterator& fi
    )
    {
        ++ fi;
    }
};

//- Define access to the base expression of a DynamicField as itself
template<class Type>
struct Operate<DynamicField<Type>, Arguments<Base>>
{
    using type = const DynamicField<Type>&;

    static type access(const DynamicField<Type>& f, const Base&)
    {
        return f;
    }
};

//- Define the types returned by DynamicField expressions as the same as for
//  Field expressions
template<class Type>
struct Result<DynamicField<Type>>
:
    public Result<Field<Type>>
{};

//- Define reduction of DynamicField expressions as the same as for Field
//  expressions
template<class Type>
struct Reduce<DynamicField<Type>>
:
    public Reduce<Field<Type>>
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //
