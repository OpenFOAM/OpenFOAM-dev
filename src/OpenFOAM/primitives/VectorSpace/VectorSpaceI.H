/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "VectorSpaceM.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline VectorSpace<Form, Cmpt, nCmpt>::VectorSpace()
{}


template<class Form, class Cmpt, int nCmpt>
inline VectorSpace<Form, Cmpt, nCmpt>::VectorSpace
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    VectorSpaceOps<nCmpt,0>::eqOp(*this, vs, eqOp<Cmpt>());
}


template<class Form, class Cmpt, int nCmpt>
template<class Form2, class Cmpt2>
inline VectorSpace<Form, Cmpt, nCmpt>::VectorSpace
(
    const VectorSpace<Form2, Cmpt2, nCmpt>& vs
)
{
    VectorSpaceOps<nCmpt,0>::eqOp(*this, vs, eqOp<Cmpt>());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline label VectorSpace<Form, Cmpt, nCmpt>::size() const
{
    return nCmpt;
}


template<class Form, class Cmpt, int nCmpt>
inline const Cmpt& VectorSpace<Form, Cmpt, nCmpt>::component
(
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn
        (
            "VectorSpace<Form, Cmpt, nCmpt>::component(direction) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt& VectorSpace<Form, Cmpt, nCmpt>::component
(
    const direction d
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn("VectorSpace<Form, Cmpt, nCmpt>::component(direction)")
            << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::component
(
    Cmpt& c,
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn
        (
            "VectorSpace<Form, Cmpt, nCmpt>::component(Cmpt&, direction) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    c = v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::replace
(
    const direction d,
    const Cmpt& c
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn
        (
            "VectorSpace<Form, Cmpt, nCmpt>::"
            "replace(direction, const Cmpt&) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    v_[d] = c;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline const Cmpt& VectorSpace<Form, Cmpt, nCmpt>::operator[]
(
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn
        (
            "VectorSpace<Form, Cmpt, nCmpt>::operator[](direction d) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt& VectorSpace<Form, Cmpt, nCmpt>::operator[]
(
    const direction d
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt)
    {
        FatalErrorIn("VectorSpace<Form, Cmpt, nCmpt>::operator[](direction d)")
            << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    VectorSpaceOps<nCmpt,0>::eqOp(*this, vs, eqOp<Cmpt>());
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator+=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    VectorSpaceOps<nCmpt,0>::eqOp(*this, vs, plusEqOp<Cmpt>());
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator-=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    VectorSpaceOps<nCmpt,0>::eqOp(*this, vs, minusEqOp<Cmpt>());
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator*=
(
    const scalar s
)
{
    VectorSpaceOps<nCmpt,0>::eqOpS(*this, s, multiplyEqOp2<Cmpt, scalar>());
}


template<class Form, class Cmpt, int nCmpt>
inline void VectorSpace<Form, Cmpt, nCmpt>::operator/=
(
    const scalar s
)
{
    VectorSpaceOps<nCmpt,0>::eqOpS(*this, s, divideEqOp2<Cmpt, scalar>());
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline Cmpt& setComponent
(
    VectorSpace<Form, Cmpt, nCmpt>& vs,
    const direction d
)
{
    return vs.component(d);
}


template<class Form, class Cmpt, int nCmpt>
inline const Cmpt& component
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs,
    const direction d
)
{
    return vs.component(d);
}


// Powers of a Form
// Equivalent to outer-products between the Form and itself
// Form^0 = 1.0
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 0>::type pow
(
    const VectorSpace<Form, Cmpt, nCmpt>&,
    typename powProduct<Form, 0>::type
    = pTraits<typename powProduct<Form, 0>::type>::zero
)
{
    return 1.0;
}

// Form^1 = Form
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 1>::type pow
(
    const VectorSpace<Form, Cmpt, nCmpt>& v,
    typename powProduct<Form, 1>::type
  = pTraits<typename powProduct<Form, 1>::type>::zero
)
{
    return static_cast<const Form&>(v);
}


// Form^2 = sqr(Form)
template<class Form, class Cmpt, int nCmpt>
inline typename powProduct<Form, 2>::type pow
(
    const VectorSpace<Form, Cmpt, nCmpt>& v,
    typename powProduct<Form, 2>::type
  = pTraits<typename powProduct<Form, 2>::type>::zero
)
{
    return sqr(static_cast<const Form&>(v));
}


template<class Form, class Cmpt, int nCmpt>
inline scalar magSqr
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    scalar ms = magSqr(vs.v_[0]);
    VectorSpaceOps<nCmpt,1>::SeqOp(ms, vs, plusEqMagSqrOp2<scalar, Cmpt>());
    return ms;
}


template<class Form, class Cmpt, int nCmpt>
inline scalar mag
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    return ::sqrt(magSqr(static_cast<const Form&>(vs)));
}


template<class Form, class Cmpt, int nCmpt>
inline VectorSpace<Form, Cmpt, nCmpt> cmptMultiply
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::op(v, vs1, vs2, cmptMultiplyOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline VectorSpace<Form, Cmpt, nCmpt> cmptDivide
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::op(v, vs1, vs2, cmptDivideOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline VectorSpace<Form, Cmpt, nCmpt> stabilise
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs,
    const Cmpt& small
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::opVS(v, vs, small, stabiliseOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptMax
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Cmpt cMax = vs.v_[0];
    VectorSpaceOps<nCmpt,1>::SeqOp(cMax, vs, maxEqOp<Cmpt>());
    return cMax;
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptMin
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Cmpt cMin = vs.v_[0];
    VectorSpaceOps<nCmpt,1>::SeqOp(cMin, vs, minEqOp<Cmpt>());
    return cMin;
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptSum
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Cmpt sum = vs.v_[0];
    VectorSpaceOps<nCmpt,1>::SeqOp(sum, vs, plusEqOp<Cmpt>());
    return sum;
}


template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptAv
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    return cmptSum(vs)/nCmpt;
}

template<class Form, class Cmpt, int nCmpt>
inline Cmpt cmptProduct
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Cmpt product = vs.v_[0];
    VectorSpaceOps<nCmpt,1>::SeqOp(product, vs, multiplyEqOp<Cmpt>());
    return product;
}


template<class Form, class Cmpt, int nCmpt>
inline Form cmptMag
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::eqOp(v, vs, eqMagOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form max
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::op(v, vs1, vs2, maxOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form min
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::op(v, vs1, vs2, minOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form minMod
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::op(v, vs1, vs2, minModOp<Cmpt>());
    return v;
}


template<class Type>
inline Type dot(const scalar s, const Type& t)
{
    return s * t;
}


template<class Type>
inline Type dot(const Type& t, const scalar s)
{
    return t * s;
}


template
<
    class Form1, class Cmpt1, int nCmpt1,
    class Form2, class Cmpt2, int nCmpt2
>
inline typename innerProduct<Form1, Form2>::type dot
(
    const VectorSpace<Form1, Cmpt1, nCmpt1>& t1,
    const VectorSpace<Form2, Cmpt2, nCmpt2>& t2
)
{
    return static_cast<const Form1&>(t1) & static_cast<const Form2&>(t2);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt, int nCmpt>
inline Form operator-
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::eqOp(v, vs, eqMinusOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator+
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::op(v, vs1, vs2, plusOp<Cmpt>());
    return v;
}

template<class Form, class Cmpt, int nCmpt>
inline Form operator-
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::op(v, vs1, vs2, minusOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator*
(
    scalar s,
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::opSV(v, s, vs, multiplyOp3<Cmpt, scalar, Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator*
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs,
    scalar s
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::opVS(v, vs, s, multiplyOp3<Cmpt, Cmpt, scalar>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator/
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs,
    scalar s
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::opVS(v, vs, s, divideOp3<Cmpt, Cmpt, scalar>());
    return v;
}

/*
template<class Form, class Cmpt, int nCmpt>
inline Form operator/
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::op(v, vs1, vs2, divideOp<Cmpt>());
    return v;
}


template<class Form, class Cmpt, int nCmpt>
inline Form operator/
(
    scalar s,
    const VectorSpace<Form, Cmpt, nCmpt>& vs
)
{
    Form v;
    VectorSpaceOps<nCmpt,0>::opSV(v, s, vs, divideOp2<scalar, Cmpt>());
    return v;
}
*/


template<class Form, class Cmpt, int nCmpt>
inline Cmpt operator&&
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    Cmpt ddProd = vs1.v_[0]*vs2.v_[0];
    for (int i=1; i<nCmpt; ++i)
    {
        ddProd += vs1.v_[i]*vs2.v_[i];
    }
    return ddProd;
}


template<class Form, class Cmpt, int nCmpt>
inline bool operator==
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    bool eq = true;
    for (int i=0; i<nCmpt; ++i)
    {
        if (!(eq &= (equal(vs1.v_[i], vs2.v_[i])))) break;
    }
    return eq;
}


template<class Form, class Cmpt, int nCmpt>
inline bool operator!=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    return !(vs1 == vs2);
}


template<class Form, class Cmpt, int nCmpt>
inline bool operator>
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    bool gt = true;
    for (int i=0; i<nCmpt; ++i)
    {
        if (!(gt &= vs1.v_[i] > vs2.v_[i])) break;
    }
    return gt;
}


template<class Form, class Cmpt, int nCmpt>
inline bool operator<
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    bool lt = true;
    for (int i=0; i<nCmpt; ++i)
    {
        if (!(lt &= vs1.v_[i] < vs2.v_[i])) break;
    }
    return lt;
}


template<class Form, class Cmpt, int nCmpt>
inline bool operator>=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    return !(vs1 < vs2);
}


template<class Form, class Cmpt, int nCmpt>
inline bool operator<=
(
    const VectorSpace<Form, Cmpt, nCmpt>& vs1,
    const VectorSpace<Form, Cmpt, nCmpt>& vs2
)
{
    return !(vs1 > vs2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
