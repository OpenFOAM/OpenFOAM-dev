/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "SphericalTensor.H"
#include "SymmTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline DiagTensor<Cmpt>::DiagTensor()
{}


template<class Cmpt>
template<class Cmpt2>
inline DiagTensor<Cmpt>::DiagTensor
(
    const VectorSpace<DiagTensor<Cmpt2>, Cmpt2, 3>& vs
)
:
    VectorSpace<DiagTensor<Cmpt>, Cmpt, 3>(vs)
{}


template<class Cmpt>
inline DiagTensor<Cmpt>::DiagTensor
(
    const Cmpt& vxx,
    const Cmpt& vyy,
    const Cmpt& vzz
)
{
    this->v_[XX] = vxx;
    this->v_[YY] = vyy;
    this->v_[ZZ] = vzz;
}


template<class Cmpt>
inline DiagTensor<Cmpt>::DiagTensor(Istream& is)
:
    VectorSpace<DiagTensor<Cmpt>, Cmpt, 3>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt&  DiagTensor<Cmpt>::xx() const
{
    return this->v_[XX];
}

template<class Cmpt>
inline const Cmpt&  DiagTensor<Cmpt>::yy() const
{
    return this->v_[YY];
}

template<class Cmpt>
inline const Cmpt&  DiagTensor<Cmpt>::zz() const
{
    return this->v_[ZZ];
}


template<class Cmpt>
inline Cmpt& DiagTensor<Cmpt>::xx()
{
    return this->v_[XX];
}

template<class Cmpt>
inline Cmpt& DiagTensor<Cmpt>::yy()
{
    return this->v_[YY];
}

template<class Cmpt>
inline Cmpt& DiagTensor<Cmpt>::zz()
{
    return this->v_[ZZ];
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline Tensor<Cmpt>
operator+(const DiagTensor<Cmpt>& dt1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        dt1.xx() + t2.xx(), t2.xy(),            t2.xz(),
        t2.yx(),            dt1.yy() + t2.yy(), t2.yz(),
        t2.zx(),            t2.zy(),            dt1.zz() + t2.zz()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator+(const Tensor<Cmpt>& t1, const DiagTensor<Cmpt>& dt2)
{
    return Tensor<Cmpt>
    (
        t1.xx() + dt2.xx(), t1.xy(),            t1.xz(),
        t1.yx(),            t1.yy() + dt2.yy(), t1.yz(),
        t1.zx(),            t1.zy(),            t1.zz() + dt2.zz()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator-(const DiagTensor<Cmpt>& dt1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        dt1.xx() - t2.xx(), -t2.xy(),            -t2.xz(),
       -t2.yx(),             dt1.yy() - t2.yy(), -t2.yz(),
       -t2.zx(),            -t2.zy(),             dt1.zz() - t2.zz()
    );
}


template<class Cmpt>
inline Tensor<Cmpt>
operator-(const Tensor<Cmpt>& t1, const DiagTensor<Cmpt>& dt2)
{
    return Tensor<Cmpt>
    (
        t1.xx() - dt2.xx(), t1.xy(),            t1.xz(),
        t1.yx(),            t1.yy() - dt2.yy(), t1.yz(),
        t1.zx(),            t1.zy(),            t1.zz() - dt2.zz()
    );
}


//- Inner-product between two diagonal tensors
template<class Cmpt>
inline DiagTensor<Cmpt>
operator&(const DiagTensor<Cmpt>& dt1, const DiagTensor<Cmpt>& dt2)
{
    return DiagTensor<Cmpt>
    (
        dt1.xx()*dt2.xx(),
        dt1.yy()*dt2.yy(),
        dt1.zz()*dt2.zz()
    );
}


//- Inner-product between a diagonal tensor and a tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const DiagTensor<Cmpt>& dt1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        dt1.xx()*t2.xx(),
        dt1.xx()*t2.xy(),
        dt1.xx()*t2.xz(),

                          dt1.yy()*t2.yx(),
                          dt1.yy()*t2.yy(),
                          dt1.yy()*t2.yz(),

                                            dt1.zz()*t2.zx(),
                                            dt1.zz()*t2.zy(),
                                            dt1.zz()*t2.zz()
    );
}


//- Inner-product between a tensor and a diagonal tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const Tensor<Cmpt>& t1, const DiagTensor<Cmpt>& dt2)
{
    return Tensor<Cmpt>
    (
        t1.xx()*dt2.xx(),
                          t1.xy()*dt2.yy(),
                                            t1.xz()*dt2.zz(),

        t1.yx()*dt2.xx(),
                          t1.yy()*dt2.yy(),
                                            t1.yz()*dt2.zz(),

        t1.zx()*dt2.xx(),
                          t1.zy()*dt2.yy(),
                                            t1.zz()*dt2.zz()
    );
}


//- Inner-product between a diagonal tensor and a vector
template<class Cmpt>
inline Vector<Cmpt>
operator&(const DiagTensor<Cmpt>& dt, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>
    (
        dt.xx()*v.x(),
                       dt.yy()*v.y(),
                                      dt.zz()*v.z()
    );
}


//- Inner-product between a vector and a diagonal tensor
template<class Cmpt>
inline Vector<Cmpt>
operator&(const Vector<Cmpt>& v, const DiagTensor<Cmpt>& dt)
{
    return Vector<Cmpt>
    (
        v.x()*dt.xx(),
                       v.y()*dt.yy(),
                                      v.z()*dt.zz()
    );
}


//- Division of a scalar by a diagonalTensor
template<class Cmpt>
inline DiagTensor<Cmpt>
operator/(const scalar s, const DiagTensor<Cmpt>& dt)
{
    return DiagTensor<Cmpt>(s/dt.xx(), s/dt.yy(), s/dt.zz());
}


//- Division of a vector by a diagonalTensor
template<class Cmpt>
inline Vector<Cmpt>
operator/(const Vector<Cmpt> v, const DiagTensor<Cmpt>& dt)
{
    return Vector<Cmpt>(v.x()/dt.xx(), v.y()/dt.yy(), v.z()/dt.zz());
}


//- Return the trace of a diagonal tensor
template<class Cmpt>
inline Cmpt tr(const DiagTensor<Cmpt>& dt)
{
    return dt.xx() + dt.yy() + dt.zz();
}


//- Return the spherical part of a diagonal tensor
template<class Cmpt>
inline SphericalTensor<Cmpt> sph(const DiagTensor<Cmpt>& dt)
{
    return 0.5*tr(dt);
}


//- Return the determinant of a diagonal tensor
template<class Cmpt>
inline Cmpt det(const DiagTensor<Cmpt>& t)
{
    return t.xx()*t.yy()*t.zz();
}


//- Return the inverse of a diagonal tensor
template<class Cmpt>
inline DiagTensor<Cmpt> inv(const DiagTensor<Cmpt>& dt)
{
    return DiagTensor<Cmpt>(1.0/dt.xx(), 1.0/dt.yy(), 1.0/dt.zz());
}


//- Return the diagonal of a symmetric tensor as a diagonal tensor
template<class Cmpt>
inline DiagTensor<Cmpt> diag(const SymmTensor<Cmpt>& t)
{
    return DiagTensor<Cmpt>(t.xx(), t.yy(), t.zz());
}


//- Return the diagonal of a tensor as a diagonal tensor
template<class Cmpt>
inline DiagTensor<Cmpt> diag(const Tensor<Cmpt>& t)
{
    return DiagTensor<Cmpt>(t.xx(), t.yy(), t.zz());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
