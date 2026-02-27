/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "boundSphere.H"

// * * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * //

inline Foam::boundSphere::boundSphere(const point& c, const scalar rSqr)
:
    c_(c),
    rSqr_(rSqr)
{}


// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

template<class PointList>
inline Foam::boundSphere::PointUList<PointList>
Foam::boundSphere::pointUList(const PointList& ps)
{
    return ps;
}


inline void Foam::boundSphere::sortBoundaryPis
(
    FixedList<label, 4>& boundaryPis
)
{
    for (label iter = 0; iter < 3; ++ iter)
    {
        for (label i = 0; i < 3; ++ i)
        {
            if
            (
                boundaryPis[i] == -1 ? true
              : boundaryPis[i + 1] == -1 ? false
              : boundaryPis[i] > boundaryPis[i + 1]
            )
            {
                Swap(boundaryPis[i], boundaryPis[i + 1]);
            }
        }
    }
}


template<class Remote>
inline void Foam::boundSphere::sortBoundaryPis
(
    FixedList<Remote, 4>& boundaryPis
)
{
    for (label iter = 0; iter < 3; ++ iter)
    {
        for (label i = 0; i < 3; ++ i)
        {
            if
            (
                boundaryPis[i].proci == -1 ? true
              : boundaryPis[i + 1].proci == -1 ? false
              : boundaryPis[i].proci > boundaryPis[i + 1].proci ? true
              : boundaryPis[i].proci < boundaryPis[i + 1].proci ? false
              : boundaryPis[i].elementi > boundaryPis[i + 1].elementi
            )
            {
                Swap(boundaryPis[i], boundaryPis[i + 1]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::boundSphere::boundSphere()
:
    c_(point::uniform(NaN)),
    rSqr_(-vGreat)
{}


inline Foam::boundSphere Foam::boundSphere::cr(const point& c, const scalar r)
{
    return r < 0 ? boundSphere() : boundSphere(c, sqr(r));
}


inline Foam::boundSphere Foam::boundSphere::cr(const Tuple2<point, scalar>& cr)
{
    return
        cr.second() < 0
      ? boundSphere()
      : boundSphere(cr.first(), sqr(cr.second()));
}


inline Foam::boundSphere Foam::boundSphere::crSqr
(
    const point& c,
    const scalar rSqr
)
{
    return rSqr < 0 ? boundSphere() : boundSphere(c, rSqr);
}


inline Foam::boundSphere Foam::boundSphere::crSqr
(
    const Tuple2<point, scalar>& crSqr
)
{
    return
        crSqr.second() < 0
      ? boundSphere()
      : boundSphere(crSqr.first(), crSqr.second());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::point& Foam::boundSphere::c() const
{
    return c_;
}


inline Foam::scalar Foam::boundSphere::rSqr() const
{
    return rSqr_;
}


inline Foam::scalar Foam::boundSphere::r() const
{
    return sqrt(rSqr_);
}


inline bool Foam::boundSphere::valid() const
{
    return rSqr_ >= 0;
}


inline bool Foam::boundSphere::contains(const point& p) const
{
    return valid() ? magSqr(p - c_) <= rSqr_ : false;
}


inline void Foam::boundSphere::inflate(const scalar s)
{
    rSqr_ *= sqr(1 + s);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline bool Foam::boundSphere::overlap
(
    const boundSphere& a,
    const boundSphere& b
)
{
    if (!a.valid() || !b.valid()) return false;

    return magSqr(a.c_ - b.c_) < sqr(a.r() + b.r());
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::intersect
(
    const PointList& ps,
    const FixedList<label, 4>& pis,
    const label nPs
)
{
    return implementation::intersect(pointUList(ps), pis, nPs);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::trivial
(
    const PointList& ps,
    const FixedList<label, 4>& pis,
    const label nPs
)
{
    FixedList<label, 4> boundaryPis({-1, -1, -1, -1});

    return implementation::trivial(pointUList(ps), pis, nPs, boundaryPis);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::trivial
(
    const PointList& ps,
    const FixedList<label, 4>& pis,
    const label nPs,
    FixedList<label, 4>& boundaryPis
)
{
    return implementation::trivial(pointUList(ps), pis, nPs, boundaryPis);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::bruteForce
(
    const PointList& ps
)
{
    FixedList<label, 4> boundaryPis({-1, -1, -1, -1});

    return implementation::bruteForce(pointUList(ps), boundaryPis);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::bruteForce
(
    const PointList& ps,
    FixedList<label, 4>& boundaryPis
)
{
    return implementation::bruteForce(pointUList(ps), boundaryPis);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::local
(
    const PointList& ps
)
{
    randomGenerator rndGen(0);
    FixedList<label, 4> boundaryPis({-1, -1, -1, -1});

    return implementation::local(pointUList(ps), rndGen, boundaryPis);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::local
(
    const PointList& ps,
    randomGenerator& rndGen
)
{
    FixedList<label, 4> boundaryPis({-1, -1, -1, -1});

    return implementation::local(pointUList(ps), rndGen, boundaryPis);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::local
(
    const PointList& ps,
    FixedList<label, 4>& boundaryPis
)
{
    randomGenerator rndGen(0);

    return implementation::local(pointUList(ps), rndGen, boundaryPis);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::local
(
    const PointList& ps,
    randomGenerator& rndGen,
    FixedList<label, 4>& boundaryPis
)
{
    return implementation::local(pointUList(ps), rndGen, boundaryPis);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::global
(
    const PointList& ps,
    const bool strict
)
{
    randomGenerator rndGen(0);
    FixedList<RemoteData<point>, 4> boundaryPis;

    return implementation::global(pointUList(ps), rndGen, boundaryPis, strict);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::global
(
    const PointList& ps,
    randomGenerator& rndGen,
    const bool strict
)
{
    FixedList<RemoteData<point>, 4> boundaryPis;

    return implementation::global(pointUList(ps), rndGen, boundaryPis, strict);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::global
(
    const PointList& ps,
    FixedList<RemoteData<point>, 4>& boundaryPis,
    const bool strict
)
{
    randomGenerator rndGen(0);

    return implementation::global(pointUList(ps), rndGen, boundaryPis, strict);
}


template<class PointList>
inline Foam::boundSphere Foam::boundSphere::global
(
    const PointList& ps,
    randomGenerator& rndGen,
    FixedList<RemoteData<point>, 4>& boundaryPis,
    const bool strict
)
{
    return implementation::global(pointUList(ps), rndGen, boundaryPis, strict);
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

inline bool Foam::operator==(const boundSphere& a, const boundSphere& b)
{
    return a.c_ == b.c_ && a.rSqr_ == b.rSqr_;
}


inline bool Foam::operator!=(const boundSphere& a, const boundSphere& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * IOstream operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, boundSphere& s)
{
    is.readBegin("boundSphere");
    is >> s.c_ >> s.rSqr_;
    is.readEnd("boundSphere");

    // Check state of Istream
    is.check("operator>>(Istream&, boundSphere&)");

    if (is.format() == IOstream::ASCII && s.valid())
    {
        s.rSqr_ = sqr(s.rSqr_);
    }

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const boundSphere& s)
{
    os << s.c_ << token::SPACE;

    if (os.format() == IOstream::ASCII && s.valid())
    {
        os << sqrt(s.rSqr_);
    }
    else
    {
        os << s.rSqr_;
    }

    return os;
}


// ************************************************************************* //
