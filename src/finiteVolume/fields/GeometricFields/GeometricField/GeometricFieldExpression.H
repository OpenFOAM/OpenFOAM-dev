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

//- Type used to identify access to the underlying internal field
struct GeometricField_InternalField {};

//- Type used to identify access to the underlying boundary field
struct GeometricField_BoundaryField {};

namespace expression
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A GeometricField is expressionable
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Expressionable<GeometricField<Type, GeoMesh, PrimitiveField>>
:
    public std::integral_constant<signed char, 1>
{};

//- Define access to the internal field of a GeometricField
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Operate
<
    GeometricField<Type, GeoMesh, PrimitiveField>,
    Arguments<GeometricField_InternalField>
>
{
    using type = const DimensionedField<Type, GeoMesh, PrimitiveField>&;

    static type access
    (
        const GeometricField<Type, GeoMesh, PrimitiveField>& f,
        const GeometricField_InternalField&
    )
    {
        return f.internalField();
    }
};

//- Define access to the boundary field of a GeometricField
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Operate
<
    GeometricField<Type, GeoMesh, PrimitiveField>,
    Arguments<Value, GeometricField_BoundaryField>
>
{
    using type =
        const
        typename GeometricField<Type, GeoMesh, PrimitiveField>::BoundaryField&;

    static type access
    (
        const GeometricField<Type, GeoMesh, PrimitiveField>& f,
        const Value&,
        const GeometricField_BoundaryField&
    )
    {
        return f.boundaryField();
    }
};

//- Define the result types for a GeometricField as the dimensioned
//  equivalents of the primitive types
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Result<GeometricField<Type, GeoMesh, PrimitiveField>>
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
                std::decay_t
                <
                    typename Operate
                    <
                        std::decay_t<Expression>,
                        Arguments<GeometricField_InternalField>
                    >::type
                >,
                Arguments<Value, Base>
            >::type
        >;

    template<class Expression>
    using elementType = dimensioned<primitiveType<Expression>>;

    template<class Expression>
    using type = GeometricField<primitiveType<Expression>, GeoMesh, Field>;
};

//- Implement a reduction of GeometricField which reduces both the primitive
//  and boundary fields (locally and then globally) and then combines the
//  result with the name and dimensions of the expression to return a
//  dimensioned type
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Reduce<GeometricField<Type, GeoMesh, PrimitiveField>>
{
    // Global (parallel) reductions should be applied to the result by default
    static const bool globalReduce = true;

    template<class Op, class Expression>
    static typename Result
    <
        GeometricField<Type, GeoMesh, PrimitiveField>
    >::template elementType<Expression> value
    (
        const Expression& e,
        const bool globalReduce,
        const label comm
    )
    {
        using BoundaryField =
            typename GeometricField<Type, GeoMesh, PrimitiveField>::
            BoundaryField;

        auto x =
            Op::value
            (
                Reduce<PrimitiveField<Type>>::template value<Op>
                (
                    access
                    (
                        access(e, GeometricField_InternalField()),
                        Value(),
                        Base()
                    ),
                    globalReduce,
                    comm
                ),
                Reduce<BoundaryField>::template value<Op>
                (
                    access
                    (
                        e,
                        Value(),
                        GeometricField_BoundaryField()
                    ),
                    globalReduce,
                    comm
                )
            );

        return
            typename Result
            <
                GeometricField<Type, GeoMesh, PrimitiveField>
            >::template elementType<Expression>
            (
                Op::name(name(e)),
                access
                (
                    access(e, GeometricField_InternalField()),
                    NullObjectRef<dimensionSet>()
                ),
                x
            );
    }

    template<class Expression>
    static typename Result
    <
        GeometricField<Type, GeoMesh, PrimitiveField>
    >::template elementType<Expression> average
    (
        const Expression& e,
        const bool globalReduce,
        const label comm
    )
    {
        return Reduce<DimensionedField<Type, GeoMesh, PrimitiveField>>::
            average
            (
                access(e, GeometricField_InternalField()),
                globalReduce,
                comm
            );
    }

    template<class Expression, class WeightExpression>
    static typename Result
    <
        GeometricField<Type, GeoMesh, PrimitiveField>
    >::template elementType<Expression> weightedAverage
    (
        const Expression& e,
        const WeightExpression& we,
        const bool globalReduce,
        const label comm
    )
    {
        return Reduce<DimensionedField<Type, GeoMesh, PrimitiveField>>::
            weightedAverage
            (
                access(e, GeometricField_InternalField()),
                access(we, GeometricField_InternalField()),
                globalReduce,
                comm
            );
    }
};

//- The name of a GeometricField is stored in its IOobject, so use this
//  rather than the default string representation of the type and pointer
template<class Type, class GeoMesh, template<class> class PrimitiveField>
struct Name
<
    GeometricField<Type, GeoMesh, PrimitiveField>,
    EnableIfExpressionable
    <
        GeometricField<Type, GeoMesh, PrimitiveField>
    >
>
{
    static word value
    (
        const GeometricField<Type, GeoMesh, PrimitiveField>& f
    )
    {
        return f.name();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expression
} // End namespace Foam

// ************************************************************************* //
