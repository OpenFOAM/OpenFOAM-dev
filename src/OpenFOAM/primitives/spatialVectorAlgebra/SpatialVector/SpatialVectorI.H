/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::SpatialVector()
{}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::SpatialVector(const Foam::zero)
:
    SpatialVector::vsType(Zero)
{}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::SpatialVector
(
    const typename SpatialVector::vsType& vs
)
:
    SpatialVector::vsType(vs)
{}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::SpatialVector
(
    const Vector<Cmpt>& w,
    const Vector<Cmpt>& l
)
{
    this->v_[0] = w.x();
    this->v_[1] = w.y();
    this->v_[2] = w.z();
    this->v_[3] = l.x();
    this->v_[4] = l.y();
    this->v_[5] = l.z();
}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::SpatialVector
(
    const Cmpt& v0,
    const Cmpt& v1,
    const Cmpt& v2,
    const Cmpt& v3,
    const Cmpt& v4,
    const Cmpt& v5
)
{
    this->v_[0] = v0;
    this->v_[1] = v1;
    this->v_[2] = v2;
    this->v_[3] = v3;
    this->v_[4] = v4;
    this->v_[5] = v5;
}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::SpatialVector(Istream& is)
:
    SpatialVector::vsType(is)
{}


template<class Cmpt>
inline Foam::SpatialVector<Cmpt>::dual::dual(const SpatialVector& v)
:
    v_(v)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::SpatialVector<Cmpt>::wx() const
{
    return this->v_[WX];
}


template<class Cmpt>
inline const Cmpt& Foam::SpatialVector<Cmpt>::wy() const
{
    return this->v_[WY];
}


template<class Cmpt>
inline const Cmpt& Foam::SpatialVector<Cmpt>::wz() const
{
    return this->v_[WZ];
}


template<class Cmpt>
inline const Cmpt& Foam::SpatialVector<Cmpt>::lx() const
{
    return this->v_[LX];
}


template<class Cmpt>
inline const Cmpt& Foam::SpatialVector<Cmpt>::ly() const
{
    return this->v_[LY];
}


template<class Cmpt>
inline const Cmpt& Foam::SpatialVector<Cmpt>::lz() const
{
    return this->v_[LZ];
}


template<class Cmpt>
inline Cmpt& Foam::SpatialVector<Cmpt>::wx()
{
    return this->v_[WX];
}


template<class Cmpt>
inline Cmpt& Foam::SpatialVector<Cmpt>::wy()
{
    return this->v_[WY];
}


template<class Cmpt>
inline Cmpt& Foam::SpatialVector<Cmpt>::wz()
{
    return this->v_[WZ];
}


template<class Cmpt>
inline Cmpt& Foam::SpatialVector<Cmpt>::lx()
{
    return this->v_[LX];
}


template<class Cmpt>
inline Cmpt& Foam::SpatialVector<Cmpt>::ly()
{
    return this->v_[LY];
}


template<class Cmpt>
inline Cmpt& Foam::SpatialVector<Cmpt>::lz()
{
    return this->v_[LZ];
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SpatialVector<Cmpt>::w() const
{
    return Vector<Cmpt>(this->v_[0], this->v_[1], this->v_[2]);
}

template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SpatialVector<Cmpt>::l() const
{
    return Vector<Cmpt>(this->v_[3], this->v_[4], this->v_[5]);
}


template<class Cmpt>
const Foam::SpatialVector<Cmpt>& Foam::SpatialVector<Cmpt>::dual::v() const
{
    return v_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline typename Foam::SpatialVector<Cmpt>::dual
Foam::SpatialVector<Cmpt>::operator*() const
{
    return dual(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Return the cross-product between two spatial vectors
template<class Cmpt>
inline SpatialVector<Cmpt> operator^
(
    const SpatialVector<Cmpt>& u,
    const SpatialVector<Cmpt>& v
)
{
    return SpatialVector<Cmpt>
    (
       -u.wz()*v.wy() + u.wy()*v.wz(),
        u.wz()*v.wx() - u.wx()*v.wz(),
       -u.wy()*v.wx() + u.wx()*v.wy(),
       -u.lz()*v.wy() + u.ly()*v.wz() - u.wz()*v.ly() + u.wy()*v.lz(),
        u.lz()*v.wx() - u.lx()*v.wz() + u.wz()*v.lx() - u.wx()*v.lz(),
       -u.ly()*v.wx() + u.lx()*v.wy() - u.wy()*v.lx() + u.wx()*v.ly()
    );
}


//- Return the dual cross-product between two spatial vectors
template<class Cmpt>
inline SpatialVector<Cmpt> operator^
(
    const SpatialVector<Cmpt>& v,
    const typename SpatialVector<Cmpt>::dual& df
)
{
    const SpatialVector<Cmpt>& f = df.v();

    return SpatialVector<Cmpt>
    (
       -v.wz()*f.wy() + v.wy()*f.wz() - v.lz()*f.ly() + v.ly()*f.lz(),
        v.wz()*f.wx() - v.wx()*f.wz() + v.lz()*f.lx() - v.lx()*f.lz(),
       -v.wy()*f.wx() + v.wx()*f.wy() - v.ly()*f.lx() + v.lx()*f.ly(),
       -v.wz()*f.ly() + v.wy()*f.lz(),
        v.wz()*f.lx() - v.wx()*f.lz(),
       -v.wy()*f.lx() + v.wx()*f.ly()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
