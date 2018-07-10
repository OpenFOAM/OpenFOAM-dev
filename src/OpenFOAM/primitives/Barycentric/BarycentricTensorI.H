/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
inline Foam::BarycentricTensor<Cmpt>::BarycentricTensor()
{}


template<class Cmpt>
inline Foam::BarycentricTensor<Cmpt>::BarycentricTensor(const Foam::zero)
:
    BarycentricTensor::msType(Zero)
{}


template<class Cmpt>
inline Foam::BarycentricTensor<Cmpt>::BarycentricTensor
(
    const Barycentric<Cmpt>& x,
    const Barycentric<Cmpt>& y,
    const Barycentric<Cmpt>& z
)
{
    this->v_[XA] = x.a();
    this->v_[XB] = x.b();
    this->v_[XC] = x.c();
    this->v_[XD] = x.d();

    this->v_[YA] = y.a();
    this->v_[YB] = y.b();
    this->v_[YC] = y.c();
    this->v_[YD] = y.d();

    this->v_[ZA] = z.a();
    this->v_[ZB] = z.b();
    this->v_[ZC] = z.c();
    this->v_[ZD] = z.d();
}


template<class Cmpt>
inline Foam::BarycentricTensor<Cmpt>::BarycentricTensor
(
    const Vector<Cmpt>& a,
    const Vector<Cmpt>& b,
    const Vector<Cmpt>& c,
    const Vector<Cmpt>& d
)
{
    this->v_[XA] = a.x();
    this->v_[XB] = b.x();
    this->v_[XC] = c.x();
    this->v_[XD] = d.x();

    this->v_[YA] = a.y();
    this->v_[YB] = b.y();
    this->v_[YC] = c.y();
    this->v_[YD] = d.y();

    this->v_[ZA] = a.z();
    this->v_[ZB] = b.z();
    this->v_[ZC] = c.z();
    this->v_[ZD] = d.z();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Barycentric<Cmpt> Foam::BarycentricTensor<Cmpt>::x() const
{
    return
        Barycentric<Cmpt>
        (
            this->v_[XA],
            this->v_[XB],
            this->v_[XC],
            this->v_[XD]
        );
}


template<class Cmpt>
inline Foam::Barycentric<Cmpt> Foam::BarycentricTensor<Cmpt>::y() const
{
    return
        Barycentric<Cmpt>
        (
            this->v_[YA],
            this->v_[YB],
            this->v_[YC],
            this->v_[YD]
        );
}


template<class Cmpt>
inline Foam::Barycentric<Cmpt> Foam::BarycentricTensor<Cmpt>::z() const
{
    return
        Barycentric<Cmpt>
        (
            this->v_[ZA],
            this->v_[ZB],
            this->v_[ZC],
            this->v_[ZD]
        );
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::BarycentricTensor<Cmpt>::a() const
{
    return Vector<Cmpt>(this->v_[XA], this->v_[YA], this->v_[ZA]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::BarycentricTensor<Cmpt>::b() const
{
    return Vector<Cmpt>(this->v_[XB], this->v_[YB], this->v_[ZB]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::BarycentricTensor<Cmpt>::c() const
{
    return Vector<Cmpt>(this->v_[XC], this->v_[YC], this->v_[ZC]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::BarycentricTensor<Cmpt>::d() const
{
    return Vector<Cmpt>(this->v_[XD], this->v_[YD], this->v_[ZD]);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline Vector<Cmpt> operator&
(
    const BarycentricTensor<Cmpt>& T,
    const Barycentric<Cmpt>& b
)
{
    return Vector<Cmpt>(T.x() & b, T.y() & b, T.z() & b);
}


template<class Cmpt>
inline Barycentric<Cmpt> operator&
(
    const Vector<Cmpt>& v,
    const BarycentricTensor<Cmpt>& T
)
{
    return Barycentric<Cmpt>(v & T.a(), v & T.b(), v & T.c(), v & T.d());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
