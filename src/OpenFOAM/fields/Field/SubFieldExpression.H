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

//- A SubField is expressionable
template<class Type>
struct Expressionable<SubField<Type>>
:
    public std::integral_constant
    <
        char,
        std::max(char(0), Expressionable<Type>::value) + char(1)
    >
{};

//- Access to a SubField by index is expressionable if the type is
template<class Type>
struct Expressionable<SubField<Type>, Arguments<label>>
:
    public Expressionable<Type>
{};

//- Access to a SubField by iterator is expressionable if the type is
template<class Type>
struct Expressionable
<
    SubField<Type>,
    Arguments<typename SubField<Type>::const_iterator>
>
:
    public Expressionable<Type>
{};

//- Define access to a value in a SubField by index
template<class Type>
struct Operate<SubField<Type>, Arguments<label>>
{
    using type = const Type&;

    static type access(const SubField<Type>& f, const label i)
    {
        return f[i];
    }
};

//- Define access to the start and end iterators of a SubField
template<class Type, bool Index>
struct Operate<SubField<Type>, Arguments<Bound<Index>>>
{
    using type = typename SubField<Type>::const_iterator;

    static type access(const SubField<Type>& f, const Bound<false>&)
    {
        return f.begin();
    }

    static type access(const SubField<Type>& f, const Bound<true>&)
    {
        return f.end();
    }
};

//- Define access to a value in a SubField by dereferencing an iterator. Also
//  define comparison and incrementation of iterators.
template<class Type, class ... Args>
struct Operate
<
    SubField<Type>,
    Arguments<typename SubField<Type>::const_iterator, Args ...>
>
{
    using type = const Type&;

    static type access
    (
        const SubField<Type>& f,
        const typename SubField<Type>::const_iterator& fi
    )
    {
        return *fi;
    }

    static bool compare
    (
        const SubField<Type>& f,
        const typename SubField<Type>::const_iterator& fi,
        const typename SubField<Type>::const_iterator& fj
    )
    {
        return fi == fj;
    }

    static void advance
    (
        const SubField<Type>& f,
        typename SubField<Type>::const_iterator& fi
    )
    {
        ++ fi;
    }
};

//- Define access to the base expression of a SubField as itself
template<class Type>
struct Operate<SubField<Type>, Arguments<Base>>
{
    using type = const SubField<Type>&;

    static type access(const SubField<Type>& f, const Base&)
    {
        return f;
    }
};

//- Define the types returned by SubField expressions as the same as for
//  Field expressions
template<class Type>
struct Result<SubField<Type>>
:
    public Result<Field<Type>>
{};

//- Define reduction of SubField expressions as the same as for Field
//  expressions
template<class Type>
struct Reduce<SubField<Type>>
:
    public Reduce<Field<Type>>
{};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //
