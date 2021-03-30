/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    3D tensor transformation operations.

\*---------------------------------------------------------------------------*/

#ifndef transform_H
#define transform_H

#include "tensor.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Rotational transformation tensor from unit vector n1 to n2
inline tensor rotationTensor
(
    const vector& n1,
    const vector& n2
)
{
    const scalar c = n1 & n2;
    const vector n3 = n1 ^ n2;
    const scalar magSqrN3 = magSqr(n3);

    // n1 and n2 define a plane n3
    if (magSqrN3 > small)
    {
        // Return rotational transformation tensor in the n3-plane using
        // Rodrigues' rotation formula
        return
            c*I
          + (n2*n1 - n1*n2)
          + (1 - c)*sqr(n3)/magSqrN3;
    }
    // n1 and n2 are contradirectional
    else if (c < 0)
    {
        // Return rotational transformation tensor in a plane with arbitrary
        // normal perpendicular to n1 and n2 using a subset of Rodrigues'
        // rotation formula
        const vector n4 = perpendicular(n1);
        const scalar magSqrN4 = magSqr(n4);
        return - I + 2*sqr(n4)/magSqrN4;
    }
    // n1 and n2 are codirectional
    else
    {
        // Return null transformation tensor
        return I;
    }
}


//- Rotational transformation tensor about the x-axis by omega radians
//  The rotation is defined in a right-handed coordinate system
//  i.e. clockwise with respect to the axis from -ve to +ve
//  (looking along the axis).
inline tensor Rx(const scalar& omega)
{
    const scalar s = sin(omega);
    const scalar c = cos(omega);
    return tensor
    (
        1,  0,  0,
        0,  c, -s,
        0, s,  c
    );
}


//- Rotational transformation tensor about the y-axis by omega radians
//  The rotation is defined in a right-handed coordinate system
//  i.e. clockwise with respect to the axis from -ve to +ve
//  (looking along the axis).
inline tensor Ry(const scalar& omega)
{
    const scalar s = sin(omega);
    const scalar c = cos(omega);
    return tensor
    (
        c,  0,  s,
        0,  1,  0,
       -s,  0,  c
    );
}


//- Rotational transformation tensor about the z-axis by omega radians
//  The rotation is defined in a right-handed coordinate system
//  i.e. clockwise with respect to the axis from -ve to +ve
//  (looking along the axis).
inline tensor Rz(const scalar& omega)
{
    const scalar s = sin(omega);
    const scalar c = cos(omega);
    return tensor
    (
        c, -s,  0,
        s,  c,  0,
        0,  0,  1
    );
}


//- Rotational transformation tensor about axis a by omega radians
//  The rotation is defined in a right-handed coordinate system
//  i.e. clockwise with respect to the axis from -ve to +ve
//  (looking along the axis).
inline tensor Ra(const vector& a, const scalar omega)
{
    const scalar s = sin(omega);
    const scalar c = cos(omega);

    return tensor
    (
        sqr(a.x())*(1 - c)  + c,
        a.y()*a.x()*(1 - c) - a.z()*s,
        a.x()*a.z()*(1 - c) + a.y()*s,

        a.x()*a.y()*(1 - c) + a.z()*s,
        sqr(a.y())*(1 - c)  + c,
        a.y()*a.z()*(1 - c) - a.x()*s,

        a.x()*a.z()*(1 - c) - a.y()*s,
        a.y()*a.z()*(1 - c) + a.x()*s,
        sqr(a.z())*(1 - c)  + c
    );
}


inline label transform(const tensor&, const bool i)
{
    return i;
}


inline label transform(const tensor&, const label i)
{
    return i;
}


inline scalar transform(const tensor&, const scalar s)
{
    return s;
}


template<class Cmpt>
inline Vector<Cmpt> transform(const tensor& tt, const Vector<Cmpt>& v)
{
    return tt & v;
}


template<class Cmpt>
inline Tensor<Cmpt> transform(const tensor& tt, const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.xx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.xy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.xz(),

        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.yx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.yy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.yz(),

        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.zx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.zy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.zz(),

        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.xx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.xy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.xz(),

        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.yx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.yy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.yz(),

        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.zx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.zy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.zz(),

        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.xx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.xy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.xz(),

        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.yx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.yy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.yz(),

        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.zx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.zy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.zz()
    );
}


template<class Cmpt>
inline SphericalTensor<Cmpt> transform
(
    const tensor& tt,
    const SphericalTensor<Cmpt>& st
)
{
    return st;
}


template<class Cmpt>
inline SymmTensor<Cmpt> transform(const tensor& tt, const SymmTensor<Cmpt>& st)
{
    return SymmTensor<Cmpt>
    (
        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.xx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.xy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.xz(),

        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.yx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.yy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.yz(),

        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.zx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.zy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.zz(),

        (tt.yx()*st.xx() + tt.yy()*st.xy() + tt.yz()*st.xz())*tt.yx()
      + (tt.yx()*st.xy() + tt.yy()*st.yy() + tt.yz()*st.yz())*tt.yy()
      + (tt.yx()*st.xz() + tt.yy()*st.yz() + tt.yz()*st.zz())*tt.yz(),

        (tt.yx()*st.xx() + tt.yy()*st.xy() + tt.yz()*st.xz())*tt.zx()
      + (tt.yx()*st.xy() + tt.yy()*st.yy() + tt.yz()*st.yz())*tt.zy()
      + (tt.yx()*st.xz() + tt.yy()*st.yz() + tt.yz()*st.zz())*tt.zz(),

        (tt.zx()*st.xx() + tt.zy()*st.xy() + tt.zz()*st.xz())*tt.zx()
      + (tt.zx()*st.xy() + tt.zy()*st.yy() + tt.zz()*st.yz())*tt.zy()
      + (tt.zx()*st.xz() + tt.zy()*st.yz() + tt.zz()*st.zz())*tt.zz()
    );
}


template<class Type1, class Type2>
inline Type1 transformMask(const Type2& t)
{
    return t;
}


template<>
inline sphericalTensor transformMask<sphericalTensor>(const tensor& t)
{
    return sph(t);
}


template<>
inline symmTensor transformMask<symmTensor>(const tensor& t)
{
    return symm(t);
}


//- Estimate angle of vec in coordinate system (e0, e1, e0^e1).
//  Is guaranteed to return increasing number but is not correct
//  angle. Used for sorting angles.  All input vectors need to be normalized.
//
//  Calculates scalar which increases with angle going from e0 to vec in
//  the coordinate system e0, e1, e0^e1
//
//  Jumps from 2*pi -> 0 at -small so hopefully parallel vectors with small
//  rounding errors should still get the same quadrant.
//
inline scalar pseudoAngle
(
    const vector& e0,
    const vector& e1,
    const vector& vec
)
{
    scalar cos = vec & e0;
    scalar sin = vec & e1;

    if (sin < -small)
    {
        return (3.0 + cos)*constant::mathematical::piByTwo;
    }
    else
    {
        return (1.0 - cos)*constant::mathematical::piByTwo;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
