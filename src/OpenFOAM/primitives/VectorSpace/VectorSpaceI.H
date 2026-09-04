/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "error.H"
#include "products.H"
#include "VectorSpaceOp.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Cmpt, direction Ncmpts>
inline VectorSpace<Form, Cmpt, Ncmpts>::VectorSpace()
{}


template<class Form, class Cmpt, direction Ncmpts>
inline VectorSpace<Form, Cmpt, Ncmpts>::VectorSpace(const Foam::zero)
{
    VectorSpaceOp<Ncmpts, 0>::op([](Cmpt& x) { x = Zero; }, *this);
}


template<class Form, class Cmpt, direction Ncmpts>
template<class Form2, class Cmpt2>
inline VectorSpace<Form, Cmpt, Ncmpts>::VectorSpace
(
    const VectorSpace<Form2, Cmpt2, Ncmpts>& vs
)
{
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt2& b) { a = b; },
        *this, vs
    );
}


template<class Form, class Cmpt, direction Ncmpts>
template<class SubVector, direction BStart>
inline
VectorSpace<Form, Cmpt, Ncmpts>::ConstBlock<SubVector, BStart>::ConstBlock
(
    const vsType& vs
)
:
    vs_(vs)
{
    static_assert
    (
        vsType::nComponents >= BStart + nComponents,
        "Requested block size > VectorSpace size"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt, direction Ncmpts>
inline direction VectorSpace<Form, Cmpt, Ncmpts>::size()
{
    return Ncmpts;
}


template<class Form, class Cmpt, direction Ncmpts>
inline const Cmpt& VectorSpace<Form, Cmpt, Ncmpts>::component
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}


template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt& VectorSpace<Form, Cmpt, Ncmpts>::component
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}


template<class Form, class Cmpt, direction Ncmpts>
inline void VectorSpace<Form, Cmpt, Ncmpts>::component
(
    Cmpt& c,
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    c = v_[d];
}


template<class Form, class Cmpt, direction Ncmpts>
inline void VectorSpace<Form, Cmpt, Ncmpts>::replace
(
    const direction d,
    const Cmpt& c
)
{
    #ifdef FULLDEBUG
    if (d >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    v_[d] = c;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form VectorSpace<Form, Cmpt, Ncmpts>::uniform(const Cmpt& s)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op([&s](Cmpt& x) { x = s; }, v);
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
template<class SubVector, direction BStart>
inline const typename VectorSpace<Form, Cmpt, Ncmpts>::template
    ConstBlock<SubVector, BStart>
VectorSpace<Form, Cmpt, Ncmpts>::block() const
{
    return *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt, direction Ncmpts>
inline const Cmpt& VectorSpace<Form, Cmpt, Ncmpts>::operator[]
(
    const direction d
) const
{
    #ifdef FULLDEBUG
    if (d >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}


template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt& VectorSpace<Form, Cmpt, Ncmpts>::operator[]
(
    const direction d
)
{
    #ifdef FULLDEBUG
    if (d >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return v_[d];
}


template<class Form, class Cmpt, direction Ncmpts>
template<class SubVector, direction BStart>
inline const Cmpt&
VectorSpace<Form, Cmpt, Ncmpts>::
ConstBlock<SubVector, BStart>::operator[]
(
    const direction i
) const
{
    #ifdef FULLDEBUG
    if (i >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }
    #endif

    return vs_[BStart + i];
}


template<class Form, class Cmpt, direction Ncmpts>
template<class SubVector, direction BStart>
inline const Cmpt&
VectorSpace<Form, Cmpt, Ncmpts>::
ConstBlock<SubVector, BStart>::operator()
(
    const direction i,
    const direction j
) const
{
    #ifdef FULLDEBUG
    if (i >= Ncmpts)
    {
        FatalErrorInFunction
            << "index out of range"
            << abort(FatalError);
    }

    if (j != 0)
    {
        FatalErrorInFunction
            << "j != 0"
            << abort(FatalError);
    }
    #endif

    return vs_[BStart + i];
}


template<class Form, class Cmpt, direction Ncmpts>
inline void VectorSpace<Form, Cmpt, Ncmpts>::operator+=
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& x, const Cmpt& y) { x += y; },
        *this, vs
    );
}


template<class Form, class Cmpt, direction Ncmpts>
inline void VectorSpace<Form, Cmpt, Ncmpts>::operator-=
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& x, const Cmpt& y) { x -= y; },
        *this, vs
    );
}


template<class Form, class Cmpt, direction Ncmpts>
inline void VectorSpace<Form, Cmpt, Ncmpts>::operator=(const Foam::zero)
{
    VectorSpaceOp<Ncmpts, 0>::op([](Cmpt& x) { x = Zero; }, *this);
}


template<class Form, class Cmpt, direction Ncmpts>
template<class Type, class Enable>
inline void VectorSpace<Form, Cmpt, Ncmpts>::operator*=
(
    const Type s
)
{
    VectorSpaceOp<Ncmpts, 0>::op([s](Cmpt& x) { x *= s; }, *this);
}


template<class Form, class Cmpt, direction Ncmpts>
template<class Type, class Enable>
inline void VectorSpace<Form, Cmpt, Ncmpts>::operator/=
(
    const Type s
)
{
    VectorSpaceOp<Ncmpts, 0>::op([&s](Cmpt& x) { x /= s; }, *this);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt& setComponent
(
    VectorSpace<Form, Cmpt, Ncmpts>& vs,
    const direction d
)
{
    return vs.component(d);
}


template<class Form, class Cmpt, direction Ncmpts>
inline const Cmpt& component
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs,
    const direction d
)
{
    return vs.component(d);
}


// Powers of a Form
// Equivalent to outer-products between the Form and itself
// Form^0 = 1.0
// template<class Form, class Cmpt, direction Ncmpts>
// inline typename powProduct<Form, 0>::type pow
// (
//     const VectorSpace<Form, Cmpt, Ncmpts>&,
//     typename powProduct<Form, 0>::type =
//         pTraits<typename powProduct<Form, 0>::type>::zero
// )
// {
//     return 1.0;
// }

// // Form^1 = Form
// template<class Form, class Cmpt, direction Ncmpts>
// inline typename powProduct<Form, 1>::type pow
// (
//     const VectorSpace<Form, Cmpt, Ncmpts>& v,
//     typename powProduct<Form, 1>::type =
//         pTraits<typename powProduct<Form, 1>::type>::zero
// )
// {
//     return static_cast<const Form&>(v);
// }


// // Form^2 = sqr(Form)
// template<class Form, class Cmpt, direction Ncmpts>
// inline typename powProduct<Form, 2>::type pow
// (
//     const VectorSpace<Form, Cmpt, Ncmpts>& v,
//     typename powProduct<Form, 2>::type =
//         pTraits<typename powProduct<Form, 2>::type>::zero
// )
// {
//     return sqr(static_cast<const Form&>(v));
// }


template<class Form, class Cmpt, direction Ncmpts>
inline auto magSqr
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    auto ms = magSqr(vs.v_[0]);
    VectorSpaceOp<Ncmpts, 1>::op
    (
        [&ms](const Cmpt& x) { ms += magSqr(x); },
        vs
    );
    return ms;
}


template<class Form, class Cmpt, direction Ncmpts>
inline auto mag
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    return sqrt(magSqr(static_cast<const Form&>(vs)));
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form normalised
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    const auto magVs = mag(vs);
    return magVs > 0 ? vs/magVs : Zero;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptMultiply
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b, const Cmpt& c) { a = cmptMultiply(b, c); },
        v,
        vs1,
        vs2
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptPow
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b, const Cmpt& c) { a = cmptPow(b, c); },
        v,
        vs1,
        vs2
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptDivide
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b, const Cmpt& c) { a = cmptDivide(b, c); },
        v,
        vs1,
        vs2
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form stabilise
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs,
    const Cmpt& small
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [&small](Cmpt& a, const Cmpt& b) { a = stabilise(b, small); },
        v,
        vs
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt cmptMax
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    Cmpt cMax = vs.v_[0];
    VectorSpaceOp<Ncmpts, 1>::op
    (
        [&cMax](const Cmpt& c) { cMax = max(c, cMax); },
        vs
    );
    return cMax;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt cmptMin
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    Cmpt cMin = vs.v_[0];
    VectorSpaceOp<Ncmpts, 1>::op
    (
        [&cMin](const Cmpt& c) { cMin = min(c, cMin); },
        vs
    );
    return cMin;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt cmptSum
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    Cmpt sum = vs.v_[0];
    VectorSpaceOp<Ncmpts, 1>::op
    (
        [&sum](const Cmpt& c) { sum += c; },
        vs
    );
    return sum;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt cmptAv
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    return cmptSum(vs)/Ncmpts;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt cmptProduct
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    Cmpt product = vs.v_[0];
    VectorSpaceOp<Ncmpts, 1>::op
    (
        [&product](const Cmpt& c) { product *= c; },
        vs
    );
    return product;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptSqr
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b) { a = sqr(b); },
        v,
        vs
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form cmptMag
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b) { a = mag(b); },
        v,
        vs
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form max
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b, const Cmpt& c) { a = max(b, c); },
        v, vs1, vs2
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form min
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b, const Cmpt& c) { a = min(b, c); },
        v, vs1, vs2
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form minMod
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b, const Cmpt& c) { a = minMod(b, c); },
        v, vs1, vs2
    );
    return v;
}


template
<
    class Form1, class Cmpt1, direction Ncmpts1,
    class Form2, class Cmpt2, direction Ncmpts2
>
inline typename innerProduct<Form1, Form2>::type dot
(
    const VectorSpace<Form1, Cmpt1, Ncmpts1>& t1,
    const VectorSpace<Form2, Cmpt2, Ncmpts2>& t2
)
{
    return static_cast<const Form1&>(t1) & static_cast<const Form2&>(t2);
}


template<class Form, class Cmpt, direction Ncmpts>
inline direction findMax(const VectorSpace<Form, Cmpt, Ncmpts>& vs)
{
    direction index = 0;
    for (direction i=1; i<Ncmpts; ++i)
    {
        index = vs[index] > vs[i] ? index : i;
    }
    return index;
}


template<class Form, class Cmpt, direction Ncmpts>
inline direction findMin(const VectorSpace<Form, Cmpt, Ncmpts>& vs)
{
    direction index = 0;
    for (direction i=1; i<Ncmpts; ++i)
    {
        index = vs[index] < vs[i] ? index : i;
    }
    return index;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator-
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b) { a = -b; },
        v, vs
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Form operator+
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b, const Cmpt& c) { a = b + c; },
        v, vs1, vs2
    );
    return v;
}

template<class Form, class Cmpt, direction Ncmpts>
inline Form operator-
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [](Cmpt& a, const Cmpt& b, const Cmpt& c) { a = b - c; },
        v, vs1, vs2
    );
    return v;
}


template
<
    class Form,
    class Cmpt,
    direction Ncmpts,
    class Type,
    class Enable = EnableRank0<Type>
>
inline Form operator*
(
    Type s,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [&s](Cmpt& a, const Cmpt& b) { a = s*b; },
        v, vs
    );
    return v;
}


template
<
    class Form,
    class Cmpt,
    direction Ncmpts,
    class Type,
    class Enable = EnableRank0<Type>
>
inline Form operator*
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs,
    Type s
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [&s](Cmpt& a, const Cmpt& b) { a = b*s; },
        v, vs
    );
    return v;
}


template
<
    class Form,
    class Cmpt,
    direction Ncmpts,
    class Type,
    class Enable = EnableRank0<Type>
>
inline Form operator/
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs,
    Type s
)
{
    Form v;
    VectorSpaceOp<Ncmpts, 0>::op
    (
        [&s](Cmpt& a, const Cmpt& b) { a = b/s; },
        v, vs
    );
    return v;
}


template<class Form, class Cmpt, direction Ncmpts>
inline Cmpt operator&&
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    Cmpt ddProd = vs1.v_[0]*vs2.v_[0];
    VectorSpaceOp<Ncmpts, 1>::op
    (
        [&ddProd](const Cmpt& a, const Cmpt& b) { ddProd += a*b; },
        vs1, vs2
    );
    return ddProd;
}


template<class Form, class Cmpt, direction Ncmpts>
inline bool operator==
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    bool eq = true;
    for (direction i=0; i<Ncmpts; ++i)
    {
        if (!(eq &= (equal(vs1.v_[i], vs2.v_[i])))) break;
    }
    return eq;
}


template<class Form, class Cmpt, direction Ncmpts>
inline bool operator!=
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    return !(vs1 == vs2);
}


template<class Form, class Cmpt, direction Ncmpts>
inline bool operator>
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    bool gt = true;
    for (direction i=0; i<Ncmpts; ++i)
    {
        if (!(gt &= vs1.v_[i] > vs2.v_[i])) break;
    }
    return gt;
}


template<class Form, class Cmpt, direction Ncmpts>
inline bool operator<
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    bool lt = true;
    for (direction i=0; i<Ncmpts; ++i)
    {
        if (!(lt &= vs1.v_[i] < vs2.v_[i])) break;
    }
    return lt;
}


template<class Form, class Cmpt, direction Ncmpts>
inline bool operator>=
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    return !(vs1 < vs2);
}


template<class Form, class Cmpt, direction Ncmpts>
inline bool operator<=
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs1,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs2
)
{
    return !(vs1 > vs2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
