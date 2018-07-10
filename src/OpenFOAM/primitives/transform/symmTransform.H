/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    3D symmetric tensor transformation operations.

\*---------------------------------------------------------------------------*/

#ifndef symmTransform_H
#define symmTransform_H

#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline scalar transform(const symmTensor&, const scalar s)
{
    return s;
}


template<class Cmpt>
inline Vector<Cmpt> transform(const symmTensor& stt, const Vector<Cmpt>& v)
{
    return stt & v;
}


template<class Cmpt>
inline Tensor<Cmpt> transform(const symmTensor& stt, const Tensor<Cmpt>& t)
{
    // return stt & t & stt.T();
    return Tensor<Cmpt>
    (
        (stt.xx()*t.xx() + stt.xy()*t.yx() + stt.xz()*t.zx())*stt.xx()
      + (stt.xx()*t.xy() + stt.xy()*t.yy() + stt.xz()*t.zy())*stt.xy()
      + (stt.xx()*t.xz() + stt.xy()*t.yz() + stt.xz()*t.zz())*stt.xz(),

        (stt.xx()*t.xx() + stt.xy()*t.yx() + stt.xz()*t.zx())*stt.xy()
      + (stt.xx()*t.xy() + stt.xy()*t.yy() + stt.xz()*t.zy())*stt.yy()
      + (stt.xx()*t.xz() + stt.xy()*t.yz() + stt.xz()*t.zz())*stt.yz(),

        (stt.xx()*t.xx() + stt.xy()*t.yx() + stt.xz()*t.zx())*stt.xz()
      + (stt.xx()*t.xy() + stt.xy()*t.yy() + stt.xz()*t.zy())*stt.yz()
      + (stt.xx()*t.xz() + stt.xy()*t.yz() + stt.xz()*t.zz())*stt.zz(),

        (stt.xy()*t.xx() + stt.yy()*t.yx() + stt.yz()*t.zx())*stt.xx()
      + (stt.xy()*t.xy() + stt.yy()*t.yy() + stt.yz()*t.zy())*stt.xy()
      + (stt.xy()*t.xz() + stt.yy()*t.yz() + stt.yz()*t.zz())*stt.xz(),

        (stt.xy()*t.xx() + stt.yy()*t.yx() + stt.yz()*t.zx())*stt.xy()
      + (stt.xy()*t.xy() + stt.yy()*t.yy() + stt.yz()*t.zy())*stt.yy()
      + (stt.xy()*t.xz() + stt.yy()*t.yz() + stt.yz()*t.zz())*stt.yz(),

        (stt.xy()*t.xx() + stt.yy()*t.yx() + stt.yz()*t.zx())*stt.xz()
      + (stt.xy()*t.xy() + stt.yy()*t.yy() + stt.yz()*t.zy())*stt.yz()
      + (stt.xy()*t.xz() + stt.yy()*t.yz() + stt.yz()*t.zz())*stt.zz(),

        (stt.xz()*t.xx() + stt.yz()*t.yx() + stt.zz()*t.zx())*stt.xx()
      + (stt.xz()*t.xy() + stt.yz()*t.yy() + stt.zz()*t.zy())*stt.xy()
      + (stt.xz()*t.xz() + stt.yz()*t.yz() + stt.zz()*t.zz())*stt.xz(),

        (stt.xz()*t.xx() + stt.yz()*t.yx() + stt.zz()*t.zx())*stt.xy()
      + (stt.xz()*t.xy() + stt.yz()*t.yy() + stt.zz()*t.zy())*stt.yy()
      + (stt.xz()*t.xz() + stt.yz()*t.yz() + stt.zz()*t.zz())*stt.yz(),

        (stt.xz()*t.xx() + stt.yz()*t.yx() + stt.zz()*t.zx())*stt.xz()
      + (stt.xz()*t.xy() + stt.yz()*t.yy() + stt.zz()*t.zy())*stt.yz()
      + (stt.xz()*t.xz() + stt.yz()*t.yz() + stt.zz()*t.zz())*stt.zz()
    );
}


template<class Cmpt>
inline SphericalTensor<Cmpt> transform
(
    const symmTensor& stt,
    const SphericalTensor<Cmpt>& st
)
{
    return st;
}


template<class Cmpt>
inline SymmTensor<Cmpt> transform
(
    const symmTensor& stt,
    const SymmTensor<Cmpt>& st
)
{
    return SymmTensor<Cmpt>
    (
        (stt.xx()*st.xx() + stt.xy()*st.xy() + stt.xz()*st.xz())*stt.xx()
      + (stt.xx()*st.xy() + stt.xy()*st.yy() + stt.xz()*st.yz())*stt.xy()
      + (stt.xx()*st.xz() + stt.xy()*st.yz() + stt.xz()*st.zz())*stt.xz(),

        (stt.xx()*st.xx() + stt.xy()*st.xy() + stt.xz()*st.xz())*stt.xy()
      + (stt.xx()*st.xy() + stt.xy()*st.yy() + stt.xz()*st.yz())*stt.yy()
      + (stt.xx()*st.xz() + stt.xy()*st.yz() + stt.xz()*st.zz())*stt.yz(),

        (stt.xx()*st.xx() + stt.xy()*st.xy() + stt.xz()*st.xz())*stt.xz()
      + (stt.xx()*st.xy() + stt.xy()*st.yy() + stt.xz()*st.yz())*stt.yz()
      + (stt.xx()*st.xz() + stt.xy()*st.yz() + stt.xz()*st.zz())*stt.zz(),

        (stt.xy()*st.xx() + stt.yy()*st.xy() + stt.yz()*st.xz())*stt.xy()
      + (stt.xy()*st.xy() + stt.yy()*st.yy() + stt.yz()*st.yz())*stt.yy()
      + (stt.xy()*st.xz() + stt.yy()*st.yz() + stt.yz()*st.zz())*stt.yz(),

        (stt.xy()*st.xx() + stt.yy()*st.xy() + stt.yz()*st.xz())*stt.xz()
      + (stt.xy()*st.xy() + stt.yy()*st.yy() + stt.yz()*st.yz())*stt.yz()
      + (stt.xy()*st.xz() + stt.yy()*st.yz() + stt.yz()*st.zz())*stt.zz(),

        (stt.xz()*st.xx() + stt.yz()*st.xy() + stt.zz()*st.xz())*stt.xz()
      + (stt.xz()*st.xy() + stt.yz()*st.yy() + stt.zz()*st.yz())*stt.yz()
      + (stt.xz()*st.xz() + stt.yz()*st.yz() + stt.zz()*st.zz())*stt.zz()
    );
}


template<>
inline sphericalTensor transformMask<sphericalTensor>(const symmTensor& st)
{
    return sph(st);
}


template<>
inline symmTensor transformMask<symmTensor>(const symmTensor& st)
{
    return st;
}


template<>
inline tensor transformMask<tensor>(const symmTensor& st)
{
    return st;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
