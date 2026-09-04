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
    Specialise expression::Operate for dimensionSet so that non-dimensioned
    types access the dimensionless set and dimensioned types access their
    dimensions.

    !!! This could/should ultimately be put directly into dimensionSet.H

\*---------------------------------------------------------------------------*/

#ifndef expressionDimensionSet_H
#define expressionDimensionSet_H

#include "dimensionedType.H"

#include "expressionOperate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expression
{

/*---------------------------------------------------------------------------*\
                     Class DimensionsOrUniformExponent Declaration
\*---------------------------------------------------------------------------*/

template<class Exponent>
class DimensionsOrUniformExponent
{
    const dimensionSet& dimensions_;

    const Exponent exponent_;

public:

    DimensionsOrUniformExponent(const Exponent& exponent)
    :
        dimensions_(dimless),
        exponent_(exponent)
    {}

    DimensionsOrUniformExponent(const dimensioned<Exponent>& dt)
    :
        dimensions_(dt.dimensions()),
        exponent_(dt.value())
    {}

    operator const dimensionSet&() const
    {
        return dimensions_;
    }

    friend dimensionSet pow
    (
        const dimensionSet& dimensions,
        const DimensionsOrUniformExponent<Exponent>& doue
    )
    {
        doue.dimensions_ = dimless;
        return pow(dimensions, doue.exponent_);
    }
};

/*---------------------------------------------------------------------------*\
                  Class dimensionsOrVariableExponent Declaration
\*---------------------------------------------------------------------------*/

class dimensionsOrVariableExponent
{
    const dimensionSet dimensions_;

public:

    dimensionsOrVariableExponent(const dimensionSet& dimensions)
    :
        dimensions_(dimensions)
    {}

    operator const dimensionSet&() const
    {
        return dimensions_;
    }

    #define UNARY_DOVE_FUNCTION(Func)                                          \
                                                                               \
        friend dimensionsOrVariableExponent Func                               \
        (                                                                      \
            const dimensionsOrVariableExponent& dove                           \
        )                                                                      \
        {                                                                      \
            return Func(dove.dimensions_);                                     \
        }

    #define UNARY_DOVE_OPERATOR(OpSym)                                         \
        UNARY_DOVE_FUNCTION(operator OpSym)

    #define BINARY_DOVE_FUNCTION(Func)                                         \
                                                                               \
        friend dimensionsOrVariableExponent Func                               \
        (                                                                      \
            const dimensionsOrVariableExponent& a,                             \
            const dimensionsOrVariableExponent& b                              \
        )                                                                      \
        {                                                                      \
            return Func(a.dimensions_, b.dimensions_);                         \
        }                                                                      \
                                                                               \
        friend dimensionsOrVariableExponent Func                               \
        (                                                                      \
            const dimensionsOrVariableExponent& dove,                          \
            const dimensionSet& dimensions                                     \
        )                                                                      \
        {                                                                      \
            return Func(dove.dimensions_, dimensions);                         \
        }                                                                      \
                                                                               \
        friend dimensionsOrVariableExponent Func                               \
        (                                                                      \
            const dimensionSet& dimensions,                                    \
            const dimensionsOrVariableExponent& dove                           \
        )                                                                      \
        {                                                                      \
            return Func(dimensions, dove.dimensions_);                         \
        }

    #define BINARY_DOVE_OPERATOR(OpSym)                                        \
        BINARY_DOVE_FUNCTION(operator OpSym)

    UNARY_DOVE_FUNCTION(sqr)
    UNARY_DOVE_FUNCTION(mag)
    UNARY_DOVE_FUNCTION(magSqr)

    UNARY_DOVE_OPERATOR(-)

    BINARY_DOVE_FUNCTION(max)
    BINARY_DOVE_FUNCTION(min)

    BINARY_DOVE_OPERATOR(+)
    BINARY_DOVE_OPERATOR(-)
    BINARY_DOVE_OPERATOR(*)
    BINARY_DOVE_OPERATOR(/)
    BINARY_DOVE_OPERATOR(&)
    BINARY_DOVE_OPERATOR(^)
    BINARY_DOVE_OPERATOR(&&)

    #undef UNARY_DOVE_FUNCTION
    #undef UNARY_DOVE_OPERATOR
    #undef BINARY_DOVE_FUNCTION
    #undef BINARY_DOVE_OPERATOR

    friend const dimensionSet& pow
    (
        const dimensionSet& dimensions,
        const dimensionsOrVariableExponent& dove
    )
    {
        dimensions = dimless;
        dove.dimensions_ = dimless;
        return dimless;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
struct Operate
<
    Type,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<Type>
>
{
    using type = const dimensionSet&;

    static type access(const Type& t, const dimensionSet&)
    {
        return dimless;
    }
};

template<>
struct Operate
<
    int32_t,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<int32_t>
>
{
    using type = DimensionsOrUniformExponent<int32_t>;

    static type access(const int32_t l, const dimensionSet&)
    {
        return type(l);
    }
};

template<>
struct Operate
<
    int64_t,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<int64_t>
>
{
    using type = DimensionsOrUniformExponent<int64_t>;

    static type access(const int64_t l, const dimensionSet&)
    {
        return type(l);
    }
};

template<>
struct Operate
<
    float,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<float>
>
{
    using type = DimensionsOrUniformExponent<float>;

    static type access(const float s, const dimensionSet&)
    {
        return type(s);
    }
};

template<>
struct Operate
<
    double,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<double>
>
{
    using type = DimensionsOrUniformExponent<double>;

    static type access(const double s, const dimensionSet&)
    {
        return type(s);
    }
};

template<>
struct Operate
<
    long double,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<double>
>
{
    using type = DimensionsOrUniformExponent<long double>;

    static type access(const long double s, const dimensionSet&)
    {
        return type(s);
    }
};

template<class Type>
struct Operate
<
    dimensioned<Type>,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<dimensioned<Type>>
>
{
    using type = const dimensionSet&;

    static type access(const dimensioned<Type>& dt, const dimensionSet&)
    {
        return dt.dimensions();
    }
};

template<>
struct Operate
<
    dimensionedScalar,
    Arguments<dimensionSet>,
    EnableIfNotExpressionable<dimensionedScalar>
>
{
    using type = DimensionsOrUniformExponent<scalar>;

    static type access(const dimensionedScalar& ds, const dimensionSet&)
    {
        return type(ds);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace expression

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
