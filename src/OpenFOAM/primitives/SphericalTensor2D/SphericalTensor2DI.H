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

\*---------------------------------------------------------------------------*/

#include "Vector2D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D()
{}


template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D(const Foam::zero)
:
    SphericalTensor2D::vsType(Zero)
{}


template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D
(
    const VectorSpace<SphericalTensor2D<Cmpt>, Cmpt, 1>& vs
)
:
    SphericalTensor2D::vsType(vs)
{}


template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D(const Cmpt& stii)
{
    this->v_[II] = stii;
}


template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D(Istream& is)
:
    SphericalTensor2D::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::SphericalTensor2D<Cmpt>::ii() const
{
    return this->v_[II];
}


template<class Cmpt>
inline Cmpt& Foam::SphericalTensor2D<Cmpt>::ii()
{
    return this->v_[II];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two spherical tensors
template<class Cmpt>
inline SphericalTensor2D<Cmpt>
operator&
(
    const SphericalTensor2D<Cmpt>& st1,
    const SphericalTensor2D<Cmpt>& st2
)
{
    return SphericalTensor2D<Cmpt>(st1.ii()*st2.ii());
}


//- Inner-product between a spherical tensor and a vector
template<class Cmpt>
inline Vector2D<Cmpt>
operator&(const SphericalTensor2D<Cmpt>& st, const Vector2D<Cmpt>& v)
{
    return Vector2D<Cmpt>
    (
        st.ii()*v.x(),
                       st.ii()*v.y()
    );
}


//- Inner-product between a vector and a spherical tensor
template<class Cmpt>
inline Vector2D<Cmpt>
operator&(const Vector2D<Cmpt>& v, const SphericalTensor2D<Cmpt>& st)
{
    return Vector2D<Cmpt>
    (
        v.x()*st.ii(),
                       v.y()*st.ii()
    );
}


//- Division of a scalar by a sphericalTensor2D
template<class Cmpt>
inline SphericalTensor2D<Cmpt>
operator/(const scalar s, const SphericalTensor2D<Cmpt>& st)
{
    return SphericalTensor2D<Cmpt>(s/st.ii());
}


//- Return the trace of a spherical tensor
template<class Cmpt>
inline Cmpt tr(const SphericalTensor2D<Cmpt>& st)
{
    return 2*st.ii();
}


//- Return the spherical part of a spherical tensor, i.e. itself
template<class Cmpt>
inline SphericalTensor2D<Cmpt> sph(const SphericalTensor2D<Cmpt>& st)
{
    return st;
}


//- Return the determinant of a spherical tensor
template<class Cmpt>
inline Cmpt det(const SphericalTensor2D<Cmpt>& st)
{
    return st.ii()*st.ii();
}


//- Return the inverse of a symmetric tensor
template<class Cmpt>
inline SphericalTensor2D<Cmpt> inv(const SphericalTensor2D<Cmpt>& st)
{
    return SphericalTensor2D<Cmpt>(1.0/st.ii());
}


template<class Cmpt>
class outerProduct<SphericalTensor2D<Cmpt>, Cmpt>
{
public:

    typedef SphericalTensor2D<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SphericalTensor2D<Cmpt>>
{
public:

    typedef SphericalTensor2D<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor2D<Cmpt>, SphericalTensor2D<Cmpt>>
{
public:

    typedef SphericalTensor2D<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor2D<Cmpt>, Vector2D<Cmpt>>
{
public:

    typedef Vector2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector2D<Cmpt>, SphericalTensor2D<Cmpt>>
{
public:

    typedef Vector2D<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
