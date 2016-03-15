/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
inline Foam::SpatialVector<Cmpt>::SpatialVector(const Foam::zero z)
:
    SpatialVector::vsType(z)
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
    const Vector<Cmpt>& angular,
    const Vector<Cmpt>& linear
)
{
    this->v_[0] = angular.x();
    this->v_[1] = angular.y();
    this->v_[2] = angular.z();
    this->v_[3] = linear.x();
    this->v_[4] = linear.y();
    this->v_[5] = linear.z();
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SpatialVector<Cmpt>::angular() const
{
    return Vector<Cmpt>(this->v_[0], this->v_[1], this->v_[2]);
}

template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SpatialVector<Cmpt>::linear() const
{
    return Vector<Cmpt>(this->v_[3], this->v_[4], this->v_[5]);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline void Foam::SpatialVector<Cmpt>::operator=(const Foam::zero z)
{
    SpatialVector::vsType::operator=(z);
}


// ************************************************************************* //
