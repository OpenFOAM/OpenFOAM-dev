/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "triangle.H"
#include "IOstreams.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::tetrahedron
(
    const Point& a,
    const Point& b,
    const Point& c,
    const Point& d
)
:
    a_(a),
    b_(b),
    c_(c),
    d_(d)
{}


template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::tetrahedron
(
    const UList<Point>& points,
    const FixedList<label, 4>& indices
)
:
    a_(points[indices[0]]),
    b_(points[indices[1]]),
    c_(points[indices[2]]),
    d_(points[indices[3]])
{}


template<class Point, class PointRef>
inline Foam::tetrahedron<Point, PointRef>::tetrahedron(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline const Point& Foam::tetrahedron<Point, PointRef>::a() const
{
    return a_;
}


template<class Point, class PointRef>
inline const Point& Foam::tetrahedron<Point, PointRef>::b() const
{
    return b_;
}


template<class Point, class PointRef>
inline const Point& Foam::tetrahedron<Point, PointRef>::c() const
{
    return c_;
}


template<class Point, class PointRef>
inline const Point& Foam::tetrahedron<Point, PointRef>::d() const
{
    return d_;
}


template<class Point, class PointRef>
inline Foam::triPointRef Foam::tetrahedron<Point, PointRef>::tri
(
    const label facei
) const
{
    // Warning. Ordering of faces needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell

    if (facei == 0)
    {
        return triPointRef(b_, c_, d_);
    }
    else if (facei == 1)
    {
        return triPointRef(a_, d_, c_);
    }
    else if (facei == 2)
    {
        return triPointRef(a_, b_, d_);
    }
    else if (facei == 3)
    {
        return triPointRef(a_, c_, b_);
    }
    else
    {
        FatalErrorInFunction
            << "index out of range 0 -> 3. facei = " << facei
            << abort(FatalError);
        return triPointRef(b_, c_, d_);
    }
}


template<class Point, class PointRef>
inline Foam::vector Foam::tetrahedron<Point, PointRef>::Sa() const
{
    return triangle<Point, PointRef>(b_, c_, d_).area();
}


template<class Point, class PointRef>
inline Foam::vector Foam::tetrahedron<Point, PointRef>::Sb() const
{
    return triangle<Point, PointRef>(a_, d_, c_).area();
}


template<class Point, class PointRef>
inline Foam::vector Foam::tetrahedron<Point, PointRef>::Sc() const
{
    return triangle<Point, PointRef>(a_, b_, d_).area();
}


template<class Point, class PointRef>
inline Foam::vector Foam::tetrahedron<Point, PointRef>::Sd() const
{
    return triangle<Point, PointRef>(a_, c_, b_).area();
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::centre() const
{
    return 0.25*(a_ + b_ + c_ + d_);
}


template<class Point, class PointRef>
inline Foam::scalar Foam::tetrahedron<Point, PointRef>::mag() const
{
    return (1.0/6.0)*(((b_ - a_) ^ (c_ - a_)) & (d_ - a_));
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::circumCentre() const
{
    vector a = b_ - a_;
    vector b = c_ - a_;
    vector c = d_ - a_;

    scalar lambda = magSqr(c) - (a & c);
    scalar mu = magSqr(b) - (a & b);

    vector ba = b ^ a;
    vector ca = c ^ a;

    vector num = lambda*ba - mu*ca;
    scalar denom = (c & ba);

    if (Foam::mag(denom) < rootVSmall)
    {
        // Degenerate tetrahedron, returning centre instead of circumCentre.

        return centre();
    }

    return a_ + 0.5*(a + num/denom);
}


template<class Point, class PointRef>
inline Foam::scalar Foam::tetrahedron<Point, PointRef>::circumRadius() const
{
    vector a = b_ - a_;
    vector b = c_ - a_;
    vector c = d_ - a_;

    scalar lambda = magSqr(c) - (a & c);
    scalar mu = magSqr(b) - (a & b);

    vector ba = b ^ a;
    vector ca = c ^ a;

    vector num = lambda*ba - mu*ca;
    scalar denom = (c & ba);

    if (Foam::mag(denom) < rootVSmall)
    {
        // Degenerate tetrahedron, returning great for circumRadius.
        return great;
    }

    return Foam::mag(0.5*(a + num/denom));
}


template<class Point, class PointRef>
inline Foam::scalar Foam::tetrahedron<Point, PointRef>::quality() const
{
    return
        mag()
       /(
            8.0/(9.0*sqrt(3.0))
           *pow3(min(circumRadius(), great))
          + rootVSmall
        );
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::randomPoint
(
    Random& rndGen
) const
{
    return barycentricToPoint(barycentric01(rndGen));
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::randomPoint
(
    cachedRandom& rndGen
) const
{
    return barycentricToPoint(barycentric01(rndGen));
}


template<class Point, class PointRef>
inline Point Foam::tetrahedron<Point, PointRef>::barycentricToPoint
(
    const barycentric& bary
) const
{
    return bary[0]*a_ + bary[1]*b_ + bary[2]*c_ + bary[3]*d_;
}


template<class Point, class PointRef>
inline Foam::barycentric Foam::tetrahedron<Point, PointRef>::pointToBarycentric
(
    const point& pt
) const
{
    barycentric bary;
    pointToBarycentric(pt, bary);
    return bary;
}


template<class Point, class PointRef>
inline Foam::scalar Foam::tetrahedron<Point, PointRef>::pointToBarycentric
(
    const point& pt,
    barycentric& bary
) const
{
    // Reference:
    // http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)

    vector e0(a_ - d_);
    vector e1(b_ - d_);
    vector e2(c_ - d_);

    tensor t
    (
        e0.x(), e1.x(), e2.x(),
        e0.y(), e1.y(), e2.y(),
        e0.z(), e1.z(), e2.z()
    );

    scalar detT = det(t);

    if (Foam::mag(detT) < small)
    {
        // Degenerate tetrahedron, returning 1/4 barycentric coordinates

        bary = barycentric(0.25, 0.25, 0.25, 0.25);

        return detT;
    }

    vector res = inv(t, detT) & (pt - d_);

    bary[0] = res.x();
    bary[1] = res.y();
    bary[2] = res.z();
    bary[3] = 1 - cmptSum(res);

    return detT;
}


template<class Point, class PointRef>
inline Foam::pointHit Foam::tetrahedron<Point, PointRef>::nearestPoint
(
    const point& p
) const
{
    // Adapted from:
    // Real-time collision detection, Christer Ericson, 2005, p142-144

    // Assuming initially that the point is inside all of the
    // halfspaces, so closest to itself.

    point closestPt = p;

    scalar minOutsideDistance = vGreat;

    bool inside = true;

    if (((p - b_) & Sa()) >= 0)
    {
        // p is outside halfspace plane of tri
        pointHit info = triangle<Point, PointRef>(b_, c_, d_).nearestPoint(p);

        inside = false;

        if (info.distance() < minOutsideDistance)
        {
            closestPt = info.rawPoint();

            minOutsideDistance = info.distance();
        }
    }

    if (((p - a_) & Sb()) >= 0)
    {
        // p is outside halfspace plane of tri
        pointHit info = triangle<Point, PointRef>(a_, d_, c_).nearestPoint(p);

        inside = false;

        if (info.distance() < minOutsideDistance)
        {
            closestPt = info.rawPoint();

            minOutsideDistance = info.distance();
        }
    }

    if (((p - a_) & Sc()) >= 0)
    {
        // p is outside halfspace plane of tri
        pointHit info = triangle<Point, PointRef>(a_, b_, d_).nearestPoint(p);

        inside = false;

        if (info.distance() < minOutsideDistance)
        {
            closestPt = info.rawPoint();

            minOutsideDistance = info.distance();
        }
    }

    if (((p - a_) & Sd()) >= 0)
    {
        // p is outside halfspace plane of tri
        pointHit info = triangle<Point, PointRef>(a_, c_, b_).nearestPoint(p);

        inside = false;

        if (info.distance() < minOutsideDistance)
        {
            closestPt = info.rawPoint();

            minOutsideDistance = info.distance();
        }
    }

    // If the point is inside, then the distance to the closest point
    // is zero
    if (inside)
    {
        minOutsideDistance = 0;
    }

    return pointHit
    (
        inside,
        closestPt,
        minOutsideDistance,
        !inside
    );
}


template<class Point, class PointRef>
bool Foam::tetrahedron<Point, PointRef>::inside(const point& pt) const
{
    // For robustness, assuming that the point is in the tet unless
    // "definitively" shown otherwise by obtaining a positive dot
    // product greater than a tolerance of small.

    // The tet is defined: tet(Cc, tetBasePt, pA, pB) where the area
    // vectors and base points for the half-space planes are:
    // area[0] = Sa();
    // area[1] = Sb();
    // area[2] = Sc();
    // area[3] = Sd();
    // planeBase[0] = tetBasePt = b_
    // planeBase[1] = ptA       = c_
    // planeBase[2] = tetBasePt = b_
    // planeBase[3] = tetBasePt = b_

    vector n = Zero;

    {
        // 0, a
        const point& basePt = b_;

        n = Sa();
        n /= (Foam::mag(n) + vSmall);

        if (((pt - basePt) & n) > small)
        {
            return false;
        }
    }

    {
        // 1, b
        const point& basePt = c_;

        n = Sb();
        n /= (Foam::mag(n) + vSmall);

        if (((pt - basePt) & n) > small)
        {
            return false;
        }
    }

    {
        // 2, c
        const point& basePt = b_;

        n = Sc();
        n /= (Foam::mag(n) + vSmall);

        if (((pt - basePt) & n) > small)
        {
            return false;
        }
    }

    {
        // 3, d
        const point& basePt = b_;

        n = Sd();
        n /= (Foam::mag(n) + vSmall);

        if (((pt - basePt) & n) > small)
        {
            return false;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    tetrahedron<Point, PointRef>& t
)
{
    is.readBegin("tetrahedron");
    is  >> t.a_ >> t.b_ >> t.c_ >> t.d_;
    is.readEnd("tetrahedron");

    is.check("Istream& operator>>(Istream&, tetrahedron&)");

    return is;
}


template<class Point, class PointRef>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tetrahedron<Point, PointRef>& t
)
{
    os  << nl
        << token::BEGIN_LIST
        << t.a_ << token::SPACE
        << t.b_ << token::SPACE
        << t.c_ << token::SPACE
        << t.d_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
