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
#include "expressionDimensionSet.H"
#include "expressionDimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A DimensionedField is expressionable
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Expressionable<DimensionedField<Type, GeoMesh, PrimitiveField>>
:
    public std::integral_constant<signed char, 1>
{};

//- Access to the dimensions of a DimensionedField is not expressionable. It
//  returns a non-expressionable object (i.e., a dimensionSet).
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Expressionable
<
    DimensionedField<Type, GeoMesh, PrimitiveField>,
    Arguments<dimensionSet>
>
:
    public std::integral_constant<signed char, 0>
{};

//- Define access to the base primitive field of a DimensionedField. Note
//  that expression::Value is specified as well as expression::Base. This is
//  not important here, but it is important when accessing dimensioned
//  constants in the same layer. The presence of expression::Value will ensure
//  that their dimensions are removed and only their value is retained in the
//  accessed expression.
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Operate
<
    DimensionedField<Type, GeoMesh, PrimitiveField>,
    Arguments<Value, Base>
>
{
    using type = const PrimitiveField<Type>&;

    static type access
    (
        const DimensionedField<Type, GeoMesh, PrimitiveField>& f,
        const Value&,
        const Base&
    )
    {
        return f.primitiveField();
    }
};

//- Access the dimensions of a DimensionedField
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Operate
<
    DimensionedField<Type, GeoMesh, PrimitiveField>,
    Arguments<dimensionSet>
>
{
    using type = const dimensionSet&;

    static type access
    (
        const DimensionedField<Type, GeoMesh, PrimitiveField>& f,
        const dimensionSet&
    )
    {
        return f.dimensions();
    }
};

//- Define the result types for a DimensionedField as the dimensioned
//  equivalents of the primitive types
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Result<DimensionedField<Type, GeoMesh, PrimitiveField>>
{
    template<class Expression>
    using primitiveType =
        typename Result
        <
            PrimitiveField<Type>
        >::template elementType
        <
            typename Operate
            <
                std::decay_t<Expression>,
                Arguments<Value, Base>
            >::type
        >;

    template<class Expression>
    using elementType = dimensioned<primitiveType<Expression>>;

    template<class Expression>
    using type =
        DimensionedField<primitiveType<Expression>, GeoMesh, Field>;
};

//- Implement a reduction of DimensionedField which reduces the primitive
//  field (locally and then globally) and then combines the result with the
//  name and dimensions of the expression to return a dimensioned type
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Reduce<DimensionedField<Type, GeoMesh, PrimitiveField>>
{
    // Global (parallel) reductions should be applied to the result by default
    static const bool globalReduce = true;

    template<class Op, class Expression>
    static typename Result
    <
        DimensionedField<Type, GeoMesh, PrimitiveField>
    >::template elementType<Expression> value(const Expression& e)
    {
        return
            typename Result
            <
                DimensionedField<Type, GeoMesh, PrimitiveField>
            >::template elementType<Expression>
            (
                Op::name(name(e)),
                access(e, NullObjectRef<dimensionSet>()),
                Reduce<PrimitiveField<Type>>::template value<Op>
                (
                    access(e, Value(), Base())
                )
            );
    }
};

//- The name of a DimensionedField is stored in its IOobject, so use this
//  rather than the default string representation of the type and pointer
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Name
<
    DimensionedField<Type, GeoMesh, PrimitiveField>,
    EnableIfExpressionable
    <
        DimensionedField<Type, GeoMesh, PrimitiveField>
    >
>
{
    static word value
    (
        const DimensionedField<Type, GeoMesh, PrimitiveField>& f
    )
    {
        return f.name();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //
