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

#include "SymmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Tensor<Cmpt>::Tensor()
{}


template<class Cmpt>
template<class Cmpt2>
inline Tensor<Cmpt>::Tensor(const VectorSpace<Tensor<Cmpt2>, Cmpt2, 9>& vs)
:
    VectorSpace<Tensor<Cmpt>, Cmpt, 9>(vs)
{}


template<class Cmpt>
inline Tensor<Cmpt>::Tensor(const SphericalTensor<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = 0; this->v_[XZ] = 0;
    this->v_[YX] = 0; this->v_[YY] = st.ii(); this->v_[YZ] = 0;
    this->v_[ZX] = 0; this->v_[ZY] = 0; this->v_[ZZ] = st.ii();
}


template<class Cmpt>
inline Tensor<Cmpt>::Tensor(const SymmTensor<Cmpt>& st)
{
    this->v_[XX] = st.xx(); this->v_[XY] = st.xy(); this->v_[XZ] = st.xz();
    this->v_[YX] = st.xy(); this->v_[YY] = st.yy(); this->v_[YZ] = st.yz();
    this->v_[ZX] = st.xz(); this->v_[ZY] = st.yz(); this->v_[ZZ] = st.zz();
}


template<class Cmpt>
inline Tensor<Cmpt>::Tensor(const Vector<Vector<Cmpt> >& tr)
{
    this->v_[XX] = tr.x().x();
    this->v_[XY] = tr.x().y();
    this->v_[XZ] = tr.x().z();

    this->v_[YX] = tr.y().x();
    this->v_[YY] = tr.y().y();
    this->v_[YZ] = tr.y().z();

    this->v_[ZX] = tr.z().x();
    this->v_[ZY] = tr.z().y();
    this->v_[ZZ] = tr.z().z();
}


template<class Cmpt>
inline Tensor<Cmpt>::Tensor
(
    const Vector<Cmpt>& x,
    const Vector<Cmpt>& y,
    const Vector<Cmpt>& z
)
{
    this->v_[XX] = x.x(); this->v_[XY] = x.y(); this->v_[XZ] = x.z();
    this->v_[YX] = y.x(); this->v_[YY] = y.y(); this->v_[YZ] = y.z();
    this->v_[ZX] = z.x(); this->v_[ZY] = z.y(); this->v_[ZZ] = z.z();
}


template<class Cmpt>
inline Tensor<Cmpt>::Tensor
(
    const Cmpt txx, const Cmpt txy, const Cmpt txz,
    const Cmpt tyx, const Cmpt tyy, const Cmpt tyz,
    const Cmpt tzx, const Cmpt tzy, const Cmpt tzz
)
{
    this->v_[XX] = txx; this->v_[XY] = txy; this->v_[XZ] = txz;
    this->v_[YX] = tyx; this->v_[YY] = tyy; this->v_[YZ] = tyz;
    this->v_[ZX] = tzx; this->v_[ZY] = tzy; this->v_[ZZ] = tzz;
}


template<class Cmpt>
inline Tensor<Cmpt>::Tensor(Istream& is)
:
    VectorSpace<Tensor<Cmpt>, Cmpt, 9>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline Vector<Cmpt> Tensor<Cmpt>::x() const
{
    return Vector<Cmpt>(this->v_[XX], this->v_[XY], this->v_[XZ]);
}


template<class Cmpt>
inline Vector<Cmpt> Tensor<Cmpt>::y() const
{
    return Vector<Cmpt>(this->v_[YX], this->v_[YY], this->v_[YZ]);
}


template<class Cmpt>
inline Vector<Cmpt>  Tensor<Cmpt>::z() const
{
    return Vector<Cmpt>(this->v_[ZX], this->v_[ZY], this->v_[ZZ]);
}


template<class Cmpt>
inline Vector<Cmpt>  Tensor<Cmpt>::vectorComponent(const direction cmpt) const
{
    switch (cmpt)
    {
        case 0:
            return x();
        break;
        case 1:
            return y();
        break;
        case 2:
            return z();
        break;
    }
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::xx() const
{
    return this->v_[XX];
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::xy() const
{
    return this->v_[XY];
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::xz() const
{
    return this->v_[XZ];
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::yx() const
{
    return this->v_[YX];
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::yy() const
{
    return this->v_[YY];
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::yz() const
{
    return this->v_[YZ];
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::zx() const
{
    return this->v_[ZX];
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::zy() const
{
    return this->v_[ZY];
}


template<class Cmpt>
inline const Cmpt&  Tensor<Cmpt>::zz() const
{
    return this->v_[ZZ];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::xx()
{
    return this->v_[XX];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::xy()
{
    return this->v_[XY];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::xz()
{
    return this->v_[XZ];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::yx()
{
    return this->v_[YX];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::yy()
{
    return this->v_[YY];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::yz()
{
    return this->v_[YZ];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::zx()
{
    return this->v_[ZX];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::zy()
{
    return this->v_[ZY];
}


template<class Cmpt>
inline Cmpt& Tensor<Cmpt>::zz()
{
    return this->v_[ZZ];
}


template<class Cmpt>
inline Tensor<Cmpt> Tensor<Cmpt>::T() const
{
    return Tensor<Cmpt>
    (
        xx(), yx(), zx(),
        xy(), yy(), zy(),
        xz(), yz(), zz()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline void Tensor<Cmpt>::operator=(const SphericalTensor<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = 0; this->v_[XZ] = 0;
    this->v_[YX] = 0; this->v_[YY] = st.ii(); this->v_[YZ] = 0;
    this->v_[ZX] = 0; this->v_[ZY] = 0; this->v_[ZZ] = st.ii();
}


template<class Cmpt>
inline void Tensor<Cmpt>::operator=(const SymmTensor<Cmpt>& st)
{
    this->v_[XX] = st.xx(); this->v_[XY] = st.xy(); this->v_[XZ] = st.xz();
    this->v_[YX] = st.xy(); this->v_[YY] = st.yy(); this->v_[YZ] = st.yz();
    this->v_[ZX] = st.xz(); this->v_[ZY] = st.yz(); this->v_[ZZ] = st.zz();
}


template<class Cmpt>
inline void Tensor<Cmpt>::operator=(const Vector<Vector<Cmpt> >& tr)
{
    this->v_[XX] = tr.x().x();
    this->v_[XY] = tr.x().y();
    this->v_[XZ] = tr.x().z();

    this->v_[YX] = tr.y().x();
    this->v_[YY] = tr.y().y();
    this->v_[YZ] = tr.y().z();

    this->v_[ZX] = tr.z().x();
    this->v_[ZY] = tr.z().y();
    this->v_[ZZ] = tr.z().z();
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline Vector<Cmpt> operator*(const Tensor<Cmpt>& t)
{
    return Vector<Cmpt>(t.yz(), -t.xz(), t.xy());
}


template<class Cmpt>
inline Tensor<Cmpt> operator*(const Vector<Cmpt>& v)
{
    return Tensor<Cmpt>
    (
             0, -v.z(),   v.y(),
         v.z(),      0,  -v.x(),
        -v.y(),  v.x(),       0
    );
}


template<class Cmpt>
inline typename innerProduct<Tensor<Cmpt>, Tensor<Cmpt> >::type
operator&(const Tensor<Cmpt>& t1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        t1.xx()*t2.xx() + t1.xy()*t2.yx() + t1.xz()*t2.zx(),
        t1.xx()*t2.xy() + t1.xy()*t2.yy() + t1.xz()*t2.zy(),
        t1.xx()*t2.xz() + t1.xy()*t2.yz() + t1.xz()*t2.zz(),

        t1.yx()*t2.xx() + t1.yy()*t2.yx() + t1.yz()*t2.zx(),
        t1.yx()*t2.xy() + t1.yy()*t2.yy() + t1.yz()*t2.zy(),
        t1.yx()*t2.xz() + t1.yy()*t2.yz() + t1.yz()*t2.zz(),

        t1.zx()*t2.xx() + t1.zy()*t2.yx() + t1.zz()*t2.zx(),
        t1.zx()*t2.xy() + t1.zy()*t2.yy() + t1.zz()*t2.zy(),
        t1.zx()*t2.xz() + t1.zy()*t2.yz() + t1.zz()*t2.zz()
    );
}


template<class Cmpt>
inline typename innerProduct<Tensor<Cmpt>, Vector<Cmpt> >::type
operator&(const Tensor<Cmpt>& t, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>
    (
        t.xx()*v.x() + t.xy()*v.y() + t.xz()*v.z(),
        t.yx()*v.x() + t.yy()*v.y() + t.yz()*v.z(),
        t.zx()*v.x() + t.zy()*v.y() + t.zz()*v.z()
    );
}


template<class Cmpt>
inline typename innerProduct<Vector<Cmpt>, Tensor<Cmpt> >::type
operator&(const Vector<Cmpt>& v, const Tensor<Cmpt>& t)
{
    return Vector<Cmpt>
    (
        v.x()*t.xx() + v.y()*t.yx() + v.z()*t.zx(),
        v.x()*t.xy() + v.y()*t.yy() + v.z()*t.zy(),
        v.x()*t.xz() + v.y()*t.yz() + v.z()*t.zz()
    );
}


template<class Cmpt>
inline typename outerProduct<Vector<Cmpt>, Vector<Cmpt> >::type
operator*(const Vector<Cmpt>& v1, const Vector<Cmpt>& v2)
{
    return Tensor<Cmpt>
    (
        v1.x()*v2.x(), v1.x()*v2.y(), v1.x()*v2.z(),
        v1.y()*v2.x(), v1.y()*v2.y(), v1.y()*v2.z(),
        v1.z()*v2.x(), v1.z()*v2.y(), v1.z()*v2.z()
    );
}


template<class Cmpt>
inline typename innerProduct<Vector<Cmpt>, Tensor<Cmpt> >::type
operator/(const Vector<Cmpt>& v, const Tensor<Cmpt>& t)
{
    return inv(t) & v;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the trace of a tensor
template<class Cmpt>
inline Cmpt tr(const Tensor<Cmpt>& t)
{
    return t.xx() + t.yy() + t.zz();
}


//- Return the spherical part of a tensor
template<class Cmpt>
inline SphericalTensor<Cmpt> sph(const Tensor<Cmpt>& t)
{
    return (1.0/3.0)*tr(t);
}


//- Return the symmetric part of a tensor
template<class Cmpt>
inline SymmTensor<Cmpt> symm(const Tensor<Cmpt>& t)
{
    return SymmTensor<Cmpt>
    (
        t.xx(), 0.5*(t.xy() + t.yx()), 0.5*(t.xz() + t.zx()),
                t.yy(),                0.5*(t.yz() + t.zy()),
                                       t.zz()
    );
}


//- Return twice the symmetric part of a tensor
template<class Cmpt>
inline SymmTensor<Cmpt> twoSymm(const Tensor<Cmpt>& t)
{
    return SymmTensor<Cmpt>
    (
        2*t.xx(), (t.xy() + t.yx()), (t.xz() + t.zx()),
                  2*t.yy(),          (t.yz() + t.zy()),
                                     2*t.zz()
    );
}


//- Return the skew-symmetric part of a tensor
template<class Cmpt>
inline Tensor<Cmpt> skew(const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        0.0, 0.5*(t.xy() - t.yx()), 0.5*(t.xz() - t.zx()),
        0.5*(t.yx() - t.xy()), 0.0, 0.5*(t.yz() - t.zy()),
        0.5*(t.zx() - t.xz()), 0.5*(t.zy() - t.yz()), 0.0
    );
}


//- Return the skew-symmetric part of a symmetric tensor
template<class Cmpt>
inline const Tensor<Cmpt>& skew(const SymmTensor<Cmpt>& st)
{
    return Tensor<Cmpt>::zero;
}


//- Return the deviatoric part of a tensor
template<class Cmpt>
inline Tensor<Cmpt> dev(const Tensor<Cmpt>& t)
{
    return t - SphericalTensor<Cmpt>::oneThirdI*tr(t);
}


//- Return the deviatoric part of a tensor
template<class Cmpt>
inline Tensor<Cmpt> dev2(const Tensor<Cmpt>& t)
{
    return t - SphericalTensor<Cmpt>::twoThirdsI*tr(t);
}


//- Return the determinant of a tensor
template<class Cmpt>
inline Cmpt det(const Tensor<Cmpt>& t)
{
    return
    (
        t.xx()*t.yy()*t.zz() + t.xy()*t.yz()*t.zx()
      + t.xz()*t.yx()*t.zy() - t.xx()*t.yz()*t.zy()
      - t.xy()*t.yx()*t.zz() - t.xz()*t.yy()*t.zx()
    );
}


//- Return the cofactor tensor of a tensor
template<class Cmpt>
inline Tensor<Cmpt> cof(const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        t.yy()*t.zz() - t.zy()*t.yz(),
        t.zx()*t.yz() - t.yx()*t.zz(),
        t.yx()*t.zy() - t.yy()*t.zx(),

        t.xz()*t.zy() - t.xy()*t.zz(),
        t.xx()*t.zz() - t.xz()*t.zx(),
        t.xy()*t.zx() - t.xx()*t.zy(),

        t.xy()*t.yz() - t.xz()*t.yy(),
        t.yx()*t.xz() - t.xx()*t.yz(),
        t.xx()*t.yy() - t.yx()*t.xy()
    );
}


//- Return the inverse of a tensor given the determinant
template<class Cmpt>
inline Tensor<Cmpt> inv(const Tensor<Cmpt>& t, const Cmpt dett)
{
    return Tensor<Cmpt>
    (
        t.yy()*t.zz() - t.zy()*t.yz(),
        t.xz()*t.zy() - t.xy()*t.zz(),
        t.xy()*t.yz() - t.xz()*t.yy(),

        t.zx()*t.yz() - t.yx()*t.zz(),
        t.xx()*t.zz() - t.xz()*t.zx(),
        t.yx()*t.xz() - t.xx()*t.yz(),

        t.yx()*t.zy() - t.yy()*t.zx(),
        t.xy()*t.zx() - t.xx()*t.zy(),
        t.xx()*t.yy() - t.yx()*t.xy()
    )/dett;
}


//- Return the inverse of a tensor
template<class Cmpt>
inline Tensor<Cmpt> inv(const Tensor<Cmpt>& t)
{
    return inv(t, det(t));
}


//- Return the 1st invariant of a tensor
template<class Cmpt>
inline Cmpt invariantI(const Tensor<Cmpt>& t)
{
    return tr(t);
}


//- Return the 2nd invariant of a tensor
template<class Cmpt>
inline Cmpt invariantII(const Tensor<Cmpt>& t)
{
    return
    (
        0.5*sqr(tr(t))
      - 0.5*
        (
           t.xx()*t.xx() + t.xy()*t.xy() + t.xz()*t.xz()
         + t.yx()*t.yx() + t.yy()*t.yy() + t.yz()*t.yz()
         + t.zx()*t.zx() + t.zy()*t.zy() + t.zz()*t.zz()
        )
    );
}


//- Return the 3rd invariant of a tensor
template<class Cmpt>
inline Cmpt invariantIII(const Tensor<Cmpt>& t)
{
    return det(t);
}


// * * * * * * * * * Mixed Tensor SphericalTensor Operators  * * * * * * * * //

template<class Cmpt>
inline Tensor<Cmpt>
operator+(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii() + t2.xx(), t2.xy(),            t2.xz(),
        t2.yx(),            st1.ii() + t2.yy(), t2.yz(),
        t2.zx(),            t2.zy(),            st1.ii() + t2.zz()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator+(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() + st2.ii(), t1.xy(),            t1.xz(),
        t1.yx(),            t1.yy() + st2.ii(), t1.yz(),
        t1.zx(),            t1.zy(),            t1.zz() + st2.ii()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator-(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii() - t2.xx(), -t2.xy(),            -t2.xz(),
       -t2.yx(),             st1.ii() - t2.yy(), -t2.yz(),
       -t2.zx(),            -t2.zy(),             st1.ii() - t2.zz()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator-(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() - st2.ii(), t1.xy(),            t1.xz(),
        t1.yx(),            t1.yy() - st2.ii(), t1.yz(),
        t1.zx(),            t1.zy(),            t1.zz() - st2.ii()
    );
}


//- Inner-product between a spherical tensor and a tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii()*t2.xx(), st1.ii()*t2.xy(), st1.ii()*t2.xz(),
        st1.ii()*t2.yx(), st1.ii()*t2.yy(), st1.ii()*t2.yz(),
        st1.ii()*t2.zx(), st1.ii()*t2.zy(), st1.ii()*t2.zz()
    );
}


//- Inner-product between a tensor and a spherical tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx()*st2.ii(), t1.xy()*st2.ii(), t1.xz()*st2.ii(),
        t1.yx()*st2.ii(), t1.yy()*st2.ii(), t1.yz()*st2.ii(),
        t1.zx()*st2.ii(), t1.zy()*st2.ii(), t1.zz()*st2.ii()
    );
}


//- Double-dot-product between a spherical tensor and a tensor
template<class Cmpt>
inline Cmpt
operator&&(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return(st1.ii()*t2.xx() + st1.ii()*t2.yy() + st1.ii()*t2.zz());
}


//- Double-dot-product between a tensor and a spherical tensor
template<class Cmpt>
inline Cmpt
operator&&(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return(t1.xx()*st2.ii() + t1.yy()*st2.ii() + t1.zz()*st2.ii());
}


template<class Cmpt>
class typeOfSum<SphericalTensor<Cmpt>, Tensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class typeOfSum<Tensor<Cmpt>, SphericalTensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, Tensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<Tensor<Cmpt>, SphericalTensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


// * * * * * * * * * * Mixed Tensor SymmTensor Operators * * * * * * * * * * //

template<class Cmpt>
inline Tensor<Cmpt>
operator+(const SymmTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.xx() + t2.xx(), st1.xy() + t2.xy(), st1.xz() + t2.xz(),
        st1.xy() + t2.yx(), st1.yy() + t2.yy(), st1.yz() + t2.yz(),
        st1.xz() + t2.zx(), st1.yz() + t2.zy(), st1.zz() + t2.zz()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator+(const Tensor<Cmpt>& t1, const SymmTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() + st2.xx(), t1.xy() + st2.xy(), t1.xz() + st2.xz(),
        t1.yx() + st2.xy(), t1.yy() + st2.yy(), t1.yz() + st2.yz(),
        t1.zx() + st2.xz(), t1.zy() + st2.yz(), t1.zz() + st2.zz()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator-(const SymmTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.xx() - t2.xx(), st1.xy() - t2.xy(), st1.xz() - t2.xz(),
        st1.xy() - t2.yx(), st1.yy() - t2.yy(), st1.yz() - t2.yz(),
        st1.xz() - t2.zx(), st1.yz() - t2.zy(), st1.zz() - t2.zz()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator-(const Tensor<Cmpt>& t1, const SymmTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() - st2.xx(), t1.xy() - st2.xy(), t1.xz() - st2.xz(),
        t1.yx() - st2.xy(), t1.yy() - st2.yy(), t1.yz() - st2.yz(),
        t1.zx() - st2.xz(), t1.zy() - st2.yz(), t1.zz() - st2.zz()
    );
}


//- Inner-product between a symmetric tensor and a tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const SymmTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.xx()*t2.xx() + st1.xy()*t2.yx() + st1.xz()*t2.zx(),
        st1.xx()*t2.xy() + st1.xy()*t2.yy() + st1.xz()*t2.zy(),
        st1.xx()*t2.xz() + st1.xy()*t2.yz() + st1.xz()*t2.zz(),

        st1.xy()*t2.xx() + st1.yy()*t2.yx() + st1.yz()*t2.zx(),
        st1.xy()*t2.xy() + st1.yy()*t2.yy() + st1.yz()*t2.zy(),
        st1.xy()*t2.xz() + st1.yy()*t2.yz() + st1.yz()*t2.zz(),

        st1.xz()*t2.xx() + st1.yz()*t2.yx() + st1.zz()*t2.zx(),
        st1.xz()*t2.xy() + st1.yz()*t2.yy() + st1.zz()*t2.zy(),
        st1.xz()*t2.xz() + st1.yz()*t2.yz() + st1.zz()*t2.zz()
    );
}


//- Inner-product between a tensor and a symmetric tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const Tensor<Cmpt>& t1, const SymmTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx()*st2.xx() + t1.xy()*st2.xy() + t1.xz()*st2.xz(),
        t1.xx()*st2.xy() + t1.xy()*st2.yy() + t1.xz()*st2.yz(),
        t1.xx()*st2.xz() + t1.xy()*st2.yz() + t1.xz()*st2.zz(),

        t1.yx()*st2.xx() + t1.yy()*st2.xy() + t1.yz()*st2.xz(),
        t1.yx()*st2.xy() + t1.yy()*st2.yy() + t1.yz()*st2.yz(),
        t1.yx()*st2.xz() + t1.yy()*st2.yz() + t1.yz()*st2.zz(),

        t1.zx()*st2.xx() + t1.zy()*st2.xy() + t1.zz()*st2.xz(),
        t1.zx()*st2.xy() + t1.zy()*st2.yy() + t1.zz()*st2.yz(),
        t1.zx()*st2.xz() + t1.zy()*st2.yz() + t1.zz()*st2.zz()
    );
}


//- Double-dot-product between a symmetric tensor and a tensor
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return
    (
        st1.xx()*t2.xx() + st1.xy()*t2.xy() + st1.xz()*t2.xz() +
        st1.xy()*t2.yx() + st1.yy()*t2.yy() + st1.yz()*t2.yz() +
        st1.xz()*t2.zx() + st1.yz()*t2.zy() + st1.zz()*t2.zz()
    );
}


//- Double-dot-product between a tensor and a symmetric tensor
template<class Cmpt>
inline Cmpt
operator&&(const Tensor<Cmpt>& t1, const SymmTensor<Cmpt>& st2)
{
    return
    (
        t1.xx()*st2.xx() + t1.xy()*st2.xy() + t1.xz()*st2.xz() +
        t1.yx()*st2.xy() + t1.yy()*st2.yy() + t1.yz()*st2.yz() +
        t1.zx()*st2.xz() + t1.zy()*st2.yz() + t1.zz()*st2.zz()
    );
}


template<class Cmpt>
class typeOfSum<SymmTensor<Cmpt>, Tensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class typeOfSum<Tensor<Cmpt>, SymmTensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SymmTensor<Cmpt>, Tensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<Tensor<Cmpt>, SymmTensor<Cmpt> >
{
public:

    typedef Tensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
