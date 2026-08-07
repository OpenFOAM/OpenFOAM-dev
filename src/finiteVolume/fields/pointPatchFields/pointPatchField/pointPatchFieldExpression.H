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

//- A pointPatchField is expressionable (yes, really ...)
template<class Type>
struct Expressionable<pointPatchField<Type>>
:
    public std::integral_constant<char, 1>
{};

//- Define access to a pointPatchField as returning a Field, but do not
//  actually provide an implementation. This just facilitates generation of the
//  identity value for reduction.
template<class Type>
struct Operate<pointPatchField<Type>, Arguments<>>
{
    using type = const Field<Type>&;
};

//- Define the types returned by pointPatchField expressions
template<class Type>
struct Result<pointPatchField<Type>>
{
    template<class Expression>
    using operateType =
        typename Operate<std::decay_t<Expression>, Arguments<>>::type;

    template<class Expression>
    using elementType =
        typename Result
        <
            Field<Type>
        >::template elementType<operateType<Expression>>;

    template<class Expression>
    using type = nil;
};

//- Define reduction of a pointPatchField as returning the identity for the
//  given operation
template<class Type>
struct Reduce<pointPatchField<Type>>
{
    template<class Op, class Expression>
    static typename Result
    <
        pointPatchField<Type>
    >::template elementType<Expression> value(const Expression& e)
    {
        return
            Op::template identity
            <
                typename Result
                <
                    pointPatchField<Type>
                >::template elementType<Expression>
            >();
    }
};

//- Implement selection of a pointPatchField within an expression. Construct a
//  calculated field type on non-constrained patch types, and pass through a
//  reference to the internal field if this was provided.
template<class Type>
struct Select
<
    pointPatchField<Type>,
    true,
    EnableIfExpressionable<pointPatchField<Type>>
>
{
    using type = autoPtr<pointPatchField<Type>>;

    template<class Expression>
    static type New
    (
        const Expression& e,
        const DimensionedField<Type, pointMesh, Field>& iF
    )
    {
        autoPtr<pointPatchField<Type>> tResult =
            pointPatchField<Type>::New
            (
                pointPatchField<Type>::calculatedType(),
                getAll<Patch<pointPatch>>(e),
                iF
            );

        tResult() = e;

        return tResult;
    }

    template<class Expression>
    static type New(const Expression& e)
    {
        return
            New
            (
                e,
                NullObjectRef<DimensionedField<Type, pointMesh, Field>>()
            );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //
