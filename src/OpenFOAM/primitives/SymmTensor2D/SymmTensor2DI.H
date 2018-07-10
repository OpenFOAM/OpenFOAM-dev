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
#include "Tensor2D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::SymmTensor2D<Cmpt>::SymmTensor2D()
{}


template<class Cmpt>
inline Foam::SymmTensor2D<Cmpt>::SymmTensor2D(const Foam::zero)
:
    SymmTensor2D::vsType(Zero)
{}


template<class Cmpt>
inline Foam::SymmTensor2D<Cmpt>::SymmTensor2D
(
    const VectorSpace<SymmTensor2D<Cmpt>, Cmpt, 3>& vs
)
:
    SymmTensor2D::vsType(vs)
{}


template<class Cmpt>
inline Foam::SymmTensor2D<Cmpt>::SymmTensor2D(const SphericalTensor2D<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = 0;
                            this->v_[YY] = st.ii();
}


template<class Cmpt>
inline Foam::SymmTensor2D<Cmpt>::SymmTensor2D
(
    const Cmpt txx, const Cmpt txy,
                    const Cmpt tyy
)
{
    this->v_[XX] = txx; this->v_[XY] = txy;
                        this->v_[YY] = tyy;
}


template<class Cmpt>
inline Foam::SymmTensor2D<Cmpt>::SymmTensor2D(Istream& is)
:
    SymmTensor2D::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor2D<Cmpt>::xx() const
{
    return this->v_[XX];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor2D<Cmpt>::xy() const
{
    return this->v_[XY];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor2D<Cmpt>::yy() const
{
    return this->v_[YY];
}


template<class Cmpt>
inline Cmpt& Foam::SymmTensor2D<Cmpt>::xx()
{
    return this->v_[XX];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor2D<Cmpt>::xy()
{
    return this->v_[XY];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor2D<Cmpt>::yy()
{
    return this->v_[YY];
}


template<class Cmpt>
inline const Foam::SymmTensor2D<Cmpt>& Foam::SymmTensor2D<Cmpt>::T() const
{
    return *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline void Foam::SymmTensor2D<Cmpt>::operator=
(
    const SphericalTensor2D<Cmpt>& st
)
{
    this->v_[XX] = st.ii(); this->v_[XY] = 0;
                            this->v_[YY] = st.ii();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two symmetric tensors
template<class Cmpt>
inline Tensor2D<Cmpt>
operator&(const SymmTensor2D<Cmpt>& st1, const SymmTensor2D<Cmpt>& st2)
{
    return Tensor2D<Cmpt>
    (
        st1.xx()*st2.xx() + st1.xy()*st2.xy(),
        st1.xx()*st2.xy() + st1.xy()*st2.yy(),

        st1.xy()*st2.xx() + st1.yy()*st2.xy(),
        st1.xy()*st2.xy() + st1.yy()*st2.yy()
    );
}


//- Double-dot-product between a symmetric tensor and a symmetric tensor
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor2D<Cmpt>& st1, const SymmTensor2D<Cmpt>& st2)
{
    return
    (
        st1.xx()*st2.xx() + 2*st1.xy()*st2.xy()
                          +   st1.yy()*st2.yy()
    );
}


//- Inner-product between a symmetric tensor and a vector
template<class Cmpt>
inline Vector2D<Cmpt>
operator&(const SymmTensor2D<Cmpt>& st, const Vector2D<Cmpt>& v)
{
    return Vector2D<Cmpt>
    (
        st.xx()*v.x() + st.xy()*v.y(),
        st.xy()*v.x() + st.yy()*v.y()
    );
}


//- Inner-product between a vector and a symmetric tensor
template<class Cmpt>
inline Vector2D<Cmpt>
operator&(const Vector2D<Cmpt>& v, const SymmTensor2D<Cmpt>& st)
{
    return Vector2D<Cmpt>
    (
        v.x()*st.xx() + v.y()*st.xy(),
        v.x()*st.xy() + v.y()*st.yy()
    );
}


//- Inner-sqr of a symmetric tensor
template<class Cmpt>
inline SymmTensor2D<Cmpt>
innerSqr(const SymmTensor2D<Cmpt>& st)
{
    return SymmTensor2D<Cmpt>
    (
        st.xx()*st.xx() + st.xy()*st.xy(),
        st.xx()*st.xy() + st.xy()*st.yy(),
        st.xy()*st.xy() + st.yy()*st.yy()
    );
}


template<class Cmpt>
inline Cmpt magSqr(const SymmTensor2D<Cmpt>& st)
{
    return
    (
        magSqr(st.xx()) + 2*magSqr(st.xy())
                        +   magSqr(st.yy())
    );
}


//- Return the trace of a symmetric tensor
template<class Cmpt>
inline Cmpt tr(const SymmTensor2D<Cmpt>& st)
{
    return st.xx() + st.yy();
}


//- Return the spherical part of a symmetric tensor
template<class Cmpt>
inline SphericalTensor2D<Cmpt> sph(const SymmTensor2D<Cmpt>& st)
{
    return (1.0/2.0)*tr(st);
}


//- Return the symmetric part of a symmetric tensor, i.e. itself
template<class Cmpt>
inline const SymmTensor2D<Cmpt>& symm(const SymmTensor2D<Cmpt>& st)
{
    return st;
}


//- Return twice the symmetric part of a symmetric tensor
template<class Cmpt>
inline SymmTensor2D<Cmpt> twoSymm(const SymmTensor2D<Cmpt>& st)
{
    return 2*st;
}


//- Return the deviatoric part of a symmetric tensor
template<class Cmpt>
inline SymmTensor2D<Cmpt> dev(const SymmTensor2D<Cmpt>& st)
{
    return st - SphericalTensor2D<Cmpt>::oneThirdI*tr(st);
}


//- Return the deviatoric part of a symmetric tensor
template<class Cmpt>
inline SymmTensor2D<Cmpt> dev2(const SymmTensor2D<Cmpt>& st)
{
    return st - SphericalTensor2D<Cmpt>::twoThirdsI*tr(st);
}


//- Return the determinant of a symmetric tensor
template<class Cmpt>
inline Cmpt det(const SymmTensor2D<Cmpt>& st)
{
    return
    (
        st.xx()*st.yy() - st.xy()*st.xy()
    );
}


//- Return the cofactor symmetric tensor of a symmetric tensor
template<class Cmpt>
inline SymmTensor2D<Cmpt> cof(const SymmTensor2D<Cmpt>& st)
{
    return SymmTensor2D<Cmpt>
    (
        st.yy(), -st.xy(),
                  st.xx()
    );
}


//- Return the inverse of a symmetric tensor give the determinant
template<class Cmpt>
inline SymmTensor2D<Cmpt> inv(const SymmTensor2D<Cmpt>& st, const Cmpt detst)
{
    return cof(st)/detst;
}


//- Return the inverse of a symmetric tensor
template<class Cmpt>
inline SymmTensor2D<Cmpt> inv(const SymmTensor2D<Cmpt>& st)
{
    return inv(st, det(st));
}


//- Return the 1st invariant of a symmetric tensor
template<class Cmpt>
inline Cmpt invariantI(const SymmTensor2D<Cmpt>& st)
{
    return tr(st);
}


//- Return the 2nd invariant of a symmetric tensor
template<class Cmpt>
inline Cmpt invariantII(const SymmTensor2D<Cmpt>& st)
{
    return
    (
        0.5*sqr(tr(st))
      - 0.5*
        (
           st.xx()*st.xx() + st.xy()*st.xy()
         + st.xy()*st.xy() + st.yy()*st.yy()
        )
    );
}


//- Return the 3rd invariant of a symmetric tensor
template<class Cmpt>
inline Cmpt invariantIII(const SymmTensor2D<Cmpt>& st)
{
    return det(st);
}


template<class Cmpt>
inline SymmTensor2D<Cmpt>
operator+(const SphericalTensor2D<Cmpt>& spt1, const SymmTensor2D<Cmpt>& st2)
{
    return SymmTensor2D<Cmpt>
    (
        spt1.ii() + st2.xx(), st2.xy(),
                              spt1.ii() + st2.yy()
    );
}


template<class Cmpt>
inline SymmTensor2D<Cmpt>
operator+(const SymmTensor2D<Cmpt>& st1, const SphericalTensor2D<Cmpt>& spt2)
{
    return SymmTensor2D<Cmpt>
    (
        st1.xx() + spt2.ii(), st1.xy(),
                              st1.yy() + spt2.ii()
    );
}


template<class Cmpt>
inline SymmTensor2D<Cmpt>
operator-(const SphericalTensor2D<Cmpt>& spt1, const SymmTensor2D<Cmpt>& st2)
{
    return SymmTensor2D<Cmpt>
    (
        spt1.ii() - st2.xx(), -st2.xy(),
                               spt1.ii() - st2.yy()
    );
}


template<class Cmpt>
inline SymmTensor2D<Cmpt>
operator-(const SymmTensor2D<Cmpt>& st1, const SphericalTensor2D<Cmpt>& spt2)
{
    return SymmTensor2D<Cmpt>
    (
        st1.xx() - spt2.ii(), st1.xy(),
                              st1.yy() - spt2.ii()
    );
}


//- Inner-product between a spherical symmetric tensor and a symmetric tensor
template<class Cmpt>
inline SymmTensor2D<Cmpt>
operator&(const SphericalTensor2D<Cmpt>& spt1, const SymmTensor2D<Cmpt>& st2)
{
    return SymmTensor2D<Cmpt>
    (
        spt1.ii()*st2.xx(), spt1.ii()*st2.xy(),
                            spt1.ii()*st2.yy()
    );
}


//- Inner-product between a tensor and a spherical tensor
template<class Cmpt>
inline SymmTensor2D<Cmpt>
operator&(const SymmTensor2D<Cmpt>& st1, const SphericalTensor2D<Cmpt>& spt2)
{
    return SymmTensor2D<Cmpt>
    (
        st1.xx()*spt2.ii(), st1.xy()*spt2.ii(),
                            st1.yy()*spt2.ii()
    );
}


//- Double-dot-product between a spherical tensor and a symmetric tensor
template<class Cmpt>
inline Cmpt
operator&&(const SphericalTensor2D<Cmpt>& spt1, const SymmTensor2D<Cmpt>& st2)
{
    return(spt1.ii()*st2.xx() + spt1.ii()*st2.yy());
}


//- Double-dot-product between a tensor and a spherical tensor
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor2D<Cmpt>& st1, const SphericalTensor2D<Cmpt>& spt2)
{
    return(st1.xx()*spt2.ii() + st1.yy()*spt2.ii());
}


template<class Cmpt>
inline SymmTensor2D<Cmpt> sqr(const Vector2D<Cmpt>& v)
{
    return SymmTensor2D<Cmpt>
    (
        v.x()*v.x(), v.x()*v.y(),
                     v.y()*v.y()
    );
}


template<class Cmpt>
class outerProduct<SymmTensor2D<Cmpt>, Cmpt>
{
public:

    typedef SymmTensor2D<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SymmTensor2D<Cmpt>>
{
public:

    typedef SymmTensor2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor2D<Cmpt>, SymmTensor2D<Cmpt>>
{
public:

    typedef Tensor2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor2D<Cmpt>, Vector2D<Cmpt>>
{
public:

    typedef Vector2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector2D<Cmpt>, SymmTensor2D<Cmpt>>
{
public:

    typedef Vector2D<Cmpt> type;
};


template<class Cmpt>
class typeOfSum<SphericalTensor2D<Cmpt>, SymmTensor2D<Cmpt>>
{
public:

    typedef SymmTensor2D<Cmpt> type;
};

template<class Cmpt>
class typeOfSum<SymmTensor2D<Cmpt>, SphericalTensor2D<Cmpt>>
{
public:

    typedef SymmTensor2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SphericalTensor2D<Cmpt>, SymmTensor2D<Cmpt>>
{
public:

    typedef SymmTensor2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor2D<Cmpt>, SphericalTensor2D<Cmpt>>
{
public:

    typedef SymmTensor2D<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
