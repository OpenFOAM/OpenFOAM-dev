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

Description
    Unary and binary vector functions and operators used by the expression
    templating system

\*---------------------------------------------------------------------------*/

#ifndef vectorExpressions_H
#define vectorExpressions_H

#include "expressions.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Unary vector functions

UnaryFunction(cmptMax)
UnaryFunction(cmptMin)
UnaryFunction(cmptAv)
UnaryFunction(cmptMag)

#define UnaryVectorOpFunction(func)                                            \
    UnaryOpFunction(vector, vector, func)

UnaryVectorOpFunction(normalised)
UnaryVectorOpFunction(perpendicular)

#undef UnaryVectorOpFunction

/*---------------------------------------------------------------------------*\
                         Class powProductOp Declaration
\*---------------------------------------------------------------------------*/

// Powers of a Form
// Equivalent to outer-products between the Form and itself
template<direction r>
struct powProductOp;

// Form^0 = 1.0
template<>
struct powProductOp<0>
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& f)
    {
        return "pow<0>(" + f + ')';
    }

    template<class Form, class Cmpt, direction Ncmpts>
    static typename powProduct<Form, 0>::type value
    (
        const VectorSpace<Form, Cmpt, Ncmpts>& v
    )
    {
        return 1.0;
    }

    static const dimensionSet& value(const dimensionSet&)
    {
        return dimless;
    }
};

// Form^1 = Form
template<>
struct powProductOp<1>
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& f)
    {
        return "pow<1>(" + f + ')';
    }

    template<class Form, class Cmpt, direction Ncmpts>
    static typename powProduct<Form, 1>::type value
    (
        const VectorSpace<Form, Cmpt, Ncmpts>& v
    )
    {
        return static_cast<const Form&>(v);
    }

    static dimensionSet value(const dimensionSet& a)
    {
        return a;
    }
};

// Form^2 = sqr(Form)
template<>
struct powProductOp<2>
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& f)
    {
        return "pow<2>(" + f + ')';
    }

    template<class Form, class Cmpt, direction Ncmpts>
    static typename powProduct<Form, 2>::type value
    (
        const VectorSpace<Form, Cmpt, Ncmpts>& v
    )
    {
        return sqr(static_cast<const Form&>(v));
    }

    static dimensionSet value(const dimensionSet& a)
    {
        return sqr(a);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Pow product unary function

template<direction r, class Form, class Cmpt, direction Ncmpts>
typename powProduct<Form, r>::type pow
(
    const VectorSpace<Form, Cmpt, Ncmpts>& v
)
{
    return powProductOp<r>::value(static_cast<const Form&>(v));
}

template<class Type, direction r>
dimensioned<typename powProduct<Type, r>::type> pow
(
    const dimensioned<Type>& dt
)
{
    return dimensioned<typename powProduct<Type, r>::type>
    (
        "pow(" + dt.name() + ',' + name(r) + ')',
        pow(dt.dimensions(), r),
        powProductOp<r>::value(dt.value())
    );
}


// template
// <
//     direction r,
//     class X,
//     class = EnableIfExpressionable<UnaryExpression<X, powProductOp<r>>>
// >
// UnaryExpression<const X&, powProductOp<r>>
// pow(const X& x)
// {
//     return UnaryExpression<const X&, powProductOp<r>>(x);
// }

#define UnaryPowProduct(XPar, XArg, XVar)                                      \
    template                                                                   \
    <                                                                          \
        direction r,                                                           \
        class X,                                                               \
        class = EnableIfExpressionable<UnaryExpression<X, powProductOp<r>>>    \
    >                                                                          \
    UnaryExpression<XPar, powProductOp<r>> pow(XArg)                           \
    {                                                                          \
        return UnaryExpression<XPar, powProductOp<r>>(XVar);                   \
    }                                                                          \

UnaryPowProduct(const X&, const X& x, x)
UnaryPowProduct(const X&, X& x, x)
UnaryPowProduct(X, X&& x, std::move(x))
UnaryPowProduct(const X&, const tmp<X>& x, x)
UnaryPowProduct(const X&, tmp<X>& x, x)
UnaryPowProduct(const X&, tmp<X>&& x, x)

#undef UnaryPowProduct

// void hmm()
// {
//     const Field<vector> f(1, vector::zero);
//     Info << pow<2>(f);
// }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Binary vector functions

BinaryFunction(cmptMultiply);
BinaryFunction(cmptDivide);

/*---------------------------------------------------------------------------*\
                         Class componentOp Declaration
\*---------------------------------------------------------------------------*/

struct componentOp
{
    static const signed char precedence =
        expression::Operation<void>::precedence;

    static const bool associative = false;

    static word name(const word& vs, const word& d)
    {
        return "component(" + vs + ',' + d + ')';
    }

    template<class VSType, class Dtype>
    static auto value(const VSType& vs, const Dtype& d)
    {
        return component(vs, d);
    }
};

BinaryFunction(component);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
