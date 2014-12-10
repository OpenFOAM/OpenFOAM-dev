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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Vector<Cmpt>::Vector()
{}


template<class Cmpt>
template<class Cmpt2>
inline Vector<Cmpt>::Vector(const VectorSpace<Vector<Cmpt2>, Cmpt2, 3>& vs)
:
    VectorSpace<Vector<Cmpt>, Cmpt, 3>(vs)
{}


template<class Cmpt>
inline Vector<Cmpt>::Vector(const Cmpt& vx, const Cmpt& vy, const Cmpt& vz)
{
    this->v_[X] = vx;
    this->v_[Y] = vy;
    this->v_[Z] = vz;
}


template<class Cmpt>
inline Vector<Cmpt>::Vector(Istream& is)
:
    VectorSpace<Vector<Cmpt>, Cmpt, 3>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt&  Vector<Cmpt>::x() const
{
    return this->v_[X];
}

template<class Cmpt>
inline const Cmpt&  Vector<Cmpt>::y() const
{
    return this->v_[Y];
}

template<class Cmpt>
inline const Cmpt&  Vector<Cmpt>::z() const
{
    return this->v_[Z];
}


template<class Cmpt>
inline Cmpt& Vector<Cmpt>::x()
{
    return this->v_[X];
}

template<class Cmpt>
inline Cmpt& Vector<Cmpt>::y()
{
    return this->v_[Y];
}

template<class Cmpt>
inline Cmpt& Vector<Cmpt>::z()
{
    return this->v_[Z];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Vector<Cmpt>& Vector<Cmpt>::centre
(
    const Foam::List<Vector<Cmpt> >&
)const
{
    return *this;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline typename innerProduct<Vector<Cmpt>, Vector<Cmpt> >::type
operator&(const Vector<Cmpt>& v1, const Vector<Cmpt>& v2)
{
    return Cmpt(v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z());
}


template<class Cmpt>
inline Vector<Cmpt> operator^(const Vector<Cmpt>& v1, const Vector<Cmpt>& v2)
{
    return Vector<Cmpt>
    (
        (v1.y()*v2.z() - v1.z()*v2.y()),
        (v1.z()*v2.x() - v1.x()*v2.z()),
        (v1.x()*v2.y() - v1.y()*v2.x())
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
