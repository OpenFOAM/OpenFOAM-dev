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

#include "IOstreams.H"
#include "pointHit.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::triangle
(
    const Point& a,
    const Point& b,
    const Point& c
)
:
    a_(a),
    b_(b),
    c_(c)
{}


template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::triangle
(
    const UList<Point>& points,
    const FixedList<label, 3>& indices
)
:
    a_(points[indices[0]]),
    b_(points[indices[1]]),
    c_(points[indices[2]])
{}



template<class Point, class PointRef>
inline Foam::triangle<Point, PointRef>::triangle(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline const Point& Foam::triangle<Point, PointRef>::a() const
{
    return a_;
}

template<class Point, class PointRef>
inline const Point& Foam::triangle<Point, PointRef>::b() const
{
    return b_;
}

template<class Point, class PointRef>
inline const Point& Foam::triangle<Point, PointRef>::c() const
{
    return c_;
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::centre() const
{
    return (1.0/3.0)*(a_ + b_ + c_);
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::mag() const
{
    return Foam::mag(normal());
}


template<class Point, class PointRef>
inline Foam::vector Foam::triangle<Point, PointRef>::normal() const
{
    return 0.5*((b_ - a_)^(c_ - a_));
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::circumCentre() const
{
    scalar d1 =  (c_ - a_) & (b_ - a_);
    scalar d2 = -(c_ - b_) & (b_ - a_);
    scalar d3 =  (c_ - a_) & (c_ - b_);

    scalar c1 = d2*d3;
    scalar c2 = d3*d1;
    scalar c3 = d1*d2;

    scalar c = c1 + c2 + c3;

    if (Foam::mag(c) < ROOTVSMALL)
    {
        // Degenerate triangle, returning centre instead of circumCentre.

        return centre();
    }

    return
    (
        ((c2 + c3)*a_ + (c3 + c1)*b_ + (c1 + c2)*c_)/(2*c)
    );
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::circumRadius() const
{
    scalar d1 =  (c_ - a_) & (b_ - a_);
    scalar d2 = -(c_ - b_) & (b_ - a_);
    scalar d3 =  (c_ - a_) & (c_ - b_);

    scalar denom = d2*d3 + d3*d1 + d1*d2;

    if (Foam::mag(denom) < VSMALL)
    {
        // Degenerate triangle, returning GREAT for circumRadius.

        return GREAT;
    }
    else
    {
        scalar a = (d1 + d2)*(d2 + d3)*(d3 + d1) / denom;

        return 0.5*Foam::sqrt(min(GREAT, max(0, a)));
    }
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::quality() const
{
    scalar c = circumRadius();

    if (c < ROOTVSMALL)
    {
        // zero circumRadius, something has gone wrong.
        return SMALL;
    }

    return mag()/(Foam::sqr(c)*3.0*sqrt(3.0)/4.0);
}


template<class Point, class PointRef>
inline Foam::scalar Foam::triangle<Point, PointRef>::sweptVol
(
    const triangle& t
) const
{
    return (1.0/12.0)*
    (
        ((t.a_ - a_) & ((b_ - a_)^(c_ - a_)))
      + ((t.b_ - b_) & ((c_ - b_)^(t.a_ - b_)))
      + ((c_ - t.c_) & ((t.b_ - t.c_)^(t.a_ - t.c_)))

      + ((t.a_ - a_) & ((b_ - a_)^(c_ - a_)))
      + ((b_ - t.b_) & ((t.a_ - t.b_)^(t.c_ - t.b_)))
      + ((c_ - t.c_) & ((b_ - t.c_)^(t.a_ - t.c_)))
    );
}


template<class Point, class PointRef>
inline Foam::tensor Foam::triangle<Point, PointRef>::inertia
(
    PointRef refPt,
    scalar density
) const
{
    Point aRel = a_ - refPt;
    Point bRel = b_ - refPt;
    Point cRel = c_ - refPt;

    tensor V
    (
        aRel.x(), aRel.y(), aRel.z(),
        bRel.x(), bRel.y(), bRel.z(),
        cRel.x(), cRel.y(), cRel.z()
    );

    scalar a = Foam::mag((b_ - a_)^(c_ - a_));

    tensor S = 1/24.0*(tensor::one + I);

    return
    (
        a*I/24.0*
        (
            (aRel & aRel)
          + (bRel & bRel)
          + (cRel & cRel)
          + ((aRel + bRel + cRel) & (aRel + bRel + cRel))
        )
      - a*(V.T() & S & V)
    )
   *density;
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::randomPoint(Random& rndGen) const
{
    // Generating Random Points in Triangles
    // by Greg Turk
    // from "Graphics Gems", Academic Press, 1990
    // http://tog.acm.org/GraphicsGems/gems/TriPoints.c

    scalar s = rndGen.scalar01();

    scalar t = sqrt(rndGen.scalar01());

    return (1 - t)*a_ + (1 - s)*t*b_ + s*t*c_;
}


template<class Point, class PointRef>
inline Point Foam::triangle<Point, PointRef>::randomPoint
(
    cachedRandom& rndGen
) const
{
    // Generating Random Points in Triangles
    // by Greg Turk
    // from "Graphics Gems", Academic Press, 1990
    // http://tog.acm.org/GraphicsGems/gems/TriPoints.c

    scalar s = rndGen.sample01<scalar>();

    scalar t = sqrt(rndGen.sample01<scalar>());

    return (1 - t)*a_ + (1 - s)*t*b_ + s*t*c_;
}


template<class Point, class PointRef>
Foam::scalar Foam::triangle<Point, PointRef>::barycentric
(
    const point& pt,
    List<scalar>& bary
) const
{
    // From:
    // Real-time collision detection, Christer Ericson, 2005, p47-48

    vector v0 = b_ - a_;
    vector v1 = c_ - a_;
    vector v2 = pt - a_;

    scalar d00 = v0 & v0;
    scalar d01 = v0 & v1;
    scalar d11 = v1 & v1;
    scalar d20 = v2 & v0;
    scalar d21 = v2 & v1;

    scalar denom = d00*d11 - d01*d01;

    if (Foam::mag(denom) < SMALL)
    {
        // Degenerate triangle, returning 1/3 barycentric coordinates.

        bary = List<scalar>(3, 1.0/3.0);

        return denom;
    }

    bary.setSize(3);

    bary[1] = (d11*d20 - d01*d21)/denom;
    bary[2] = (d00*d21 - d01*d20)/denom;
    bary[0] = 1.0 - bary[1] - bary[2];

    return denom;
}


template<class Point, class PointRef>
inline Foam::pointHit Foam::triangle<Point, PointRef>::ray
(
    const point& p,
    const vector& q,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    // Express triangle in terms of baseVertex (point a_) and
    // two edge vectors
    vector E0 = b_ - a_;
    vector E1 = c_ - a_;

    // Initialize intersection to miss.
    pointHit inter(p);

    vector n(0.5*(E0 ^ E1));
    scalar area = Foam::mag(n);

    if (area < VSMALL)
    {
        // Ineligible miss.
        inter.setMiss(false);

        // The miss point is the nearest point on the triangle. Take any one.
        inter.setPoint(a_);

        // The distance to the miss is the distance between the
        // original point and plane of intersection. No normal so use
        // distance from p to a_
        inter.setDistance(Foam::mag(a_ - p));

        return inter;
    }

    vector q1 = q/Foam::mag(q);

    if (dir == intersection::CONTACT_SPHERE)
    {
        n /= area;

        return ray(p, q1 - n, alg, intersection::VECTOR);
    }

    // Intersection point with triangle plane
    point pInter;
    // Is intersection point inside triangle
    bool hit;
    {
        // Reuse the fast ray intersection routine below in FULL_RAY
        // mode since the original intersection routine has rounding problems.
        pointHit fastInter = intersection(p, q1, intersection::FULL_RAY);
        hit = fastInter.hit();

        if (hit)
        {
            pInter = fastInter.rawPoint();
        }
        else
        {
            // Calculate intersection of ray with triangle plane
            vector v = a_ - p;
            pInter = p + (q1&v)*q1;
        }
    }

    // Distance to intersection point
    scalar dist = q1 & (pInter - p);

    const scalar planarPointTol =
        Foam::min
        (
            Foam::min
            (
                Foam::mag(E0),
                Foam::mag(E1)
            ),
            Foam::mag(c_ - b_)
        )*intersection::planarTol();

    bool eligible =
        alg == intersection::FULL_RAY
     || (alg == intersection::HALF_RAY && dist > -planarPointTol)
     || (
            alg == intersection::VISIBLE
         && ((q1 & normal()) < -VSMALL)
        );

    if (hit && eligible)
    {
        // Hit. Set distance to intersection.
        inter.setHit();
        inter.setPoint(pInter);
        inter.setDistance(dist);
    }
    else
    {
        // Miss or ineligible hit.
        inter.setMiss(eligible);

        // The miss point is the nearest point on the triangle
        inter.setPoint(nearestPoint(p).rawPoint());

        // The distance to the miss is the distance between the
        // original point and plane of intersection
        inter.setDistance(Foam::mag(pInter - p));
    }


    return inter;
}


// From "Fast, Minimum Storage Ray/Triangle Intersection"
// Moeller/Trumbore.
template<class Point, class PointRef>
inline Foam::pointHit Foam::triangle<Point, PointRef>::intersection
(
    const point& orig,
    const vector& dir,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    const vector edge1 = b_ - a_;
    const vector edge2 = c_ - a_;

    // begin calculating determinant - also used to calculate U parameter
    const vector pVec = dir ^ edge2;

    // if determinant is near zero, ray lies in plane of triangle
    const scalar det = edge1 & pVec;

    // Initialise to miss
    pointHit intersection(false, vector::zero, GREAT, false);

    if (alg == intersection::VISIBLE)
    {
        // Culling branch
        if (det < ROOTVSMALL)
        {
            // Ray on wrong side of triangle. Return miss
            return intersection;
        }
    }
    else if (alg == intersection::HALF_RAY || alg == intersection::FULL_RAY)
    {
        // Non-culling branch
        if (det > -ROOTVSMALL && det < ROOTVSMALL)
        {
            // Ray parallel to triangle. Return miss
            return intersection;
        }
    }

    const scalar inv_det = 1.0 / det;

    /* calculate distance from a_ to ray origin */
    const vector tVec = orig-a_;

    /* calculate U parameter and test bounds */
    const scalar u = (tVec & pVec)*inv_det;

    if (u < -tol || u > 1.0+tol)
    {
        // return miss
        return intersection;
    }

    /* prepare to test V parameter */
    const vector qVec = tVec ^ edge1;

    /* calculate V parameter and test bounds */
    const scalar v = (dir & qVec) * inv_det;

    if (v < -tol || u + v > 1.0+tol)
    {
        // return miss
        return intersection;
    }

    /* calculate t, scale parameters, ray intersects triangle */
    const scalar t = (edge2 & qVec) * inv_det;

    if (alg == intersection::HALF_RAY && t < -tol)
    {
        // Wrong side of orig. Return miss
        return intersection;
    }

    intersection.setHit();
    intersection.setPoint(a_ + u*edge1 + v*edge2);
    intersection.setDistance(t);

    return intersection;
}


template<class Point, class PointRef>
Foam::pointHit Foam::triangle<Point, PointRef>::nearestPointClassify
(
    const point& p,
    label& nearType,
    label& nearLabel
) const
{
    // Adapted from:
    // Real-time collision detection, Christer Ericson, 2005, p136-142

    // Check if P in vertex region outside A
    vector ab = b_ - a_;
    vector ac = c_ - a_;
    vector ap = p - a_;

    scalar d1 = ab & ap;
    scalar d2 = ac & ap;

    if (d1 <= 0.0 && d2 <= 0.0)
    {
        // barycentric coordinates (1, 0, 0)

        nearType = POINT;
        nearLabel = 0;
        return pointHit(false, a_, Foam::mag(a_ - p), true);
    }

    // Check if P in vertex region outside B
    vector bp = p - b_;
    scalar d3 = ab & bp;
    scalar d4 = ac & bp;

    if (d3 >= 0.0 && d4 <= d3)
    {
        // barycentric coordinates (0, 1, 0)

        nearType = POINT;
        nearLabel = 1;
        return pointHit(false, b_, Foam::mag(b_ - p), true);
    }

    // Check if P in edge region of AB, if so return projection of P onto AB
    scalar vc = d1*d4 - d3*d2;

    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        if ((d1 - d3) < ROOTVSMALL)
        {
            // Degenerate triangle, for d1 = d3, a_ and b_ are likely coincident
            nearType = POINT;
            nearLabel = 0;
            return pointHit(false, a_, Foam::mag(a_ - p), true);
        }

        // barycentric coordinates (1-v, v, 0)
        scalar v = d1/(d1 - d3);

        point nearPt =  a_ + v*ab;
        nearType = EDGE;
        nearLabel = 0;
        return pointHit(false, nearPt, Foam::mag(nearPt - p), true);
    }

    // Check if P in vertex region outside C
    vector cp = p - c_;
    scalar d5 = ab & cp;
    scalar d6 = ac & cp;

    if (d6 >= 0.0 && d5 <= d6)
    {
        // barycentric coordinates (0, 0, 1)

        nearType = POINT;
        nearLabel = 2;
        return pointHit(false, c_, Foam::mag(c_ - p), true);
    }

    // Check if P in edge region of AC, if so return projection of P onto AC
    scalar vb = d5*d2 - d1*d6;

    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        if ((d2 - d6) < ROOTVSMALL)
        {
            // Degenerate triangle, for d2 = d6, a_ and c_ are likely coincident
            nearType = POINT;
            nearLabel = 0;
            return pointHit(false, a_, Foam::mag(a_ - p), true);
        }

        // barycentric coordinates (1-w, 0, w)
        scalar w = d2/(d2 - d6);

        point nearPt = a_ + w*ac;
        nearType = EDGE;
        nearLabel = 2;
        return pointHit(false, nearPt, Foam::mag(nearPt - p), true);
    }

    // Check if P in edge region of BC, if so return projection of P onto BC
    scalar va = d3*d6 - d5*d4;

    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        if (((d4 - d3) + (d5 - d6)) < ROOTVSMALL)
        {
            // Degenerate triangle, for (d4 - d3) = (d6 - d5), b_ and c_ are
            // likely coincident
            nearType = POINT;
            nearLabel = 1;
            return pointHit(false, b_, Foam::mag(b_ - p), true);
        }

        // barycentric coordinates (0, 1-w, w)
        scalar w = (d4 - d3)/((d4 - d3) + (d5 - d6));

        point nearPt = b_ + w*(c_ - b_);
        nearType = EDGE;
        nearLabel = 1;
        return pointHit(false, nearPt, Foam::mag(nearPt - p), true);
    }

    // P inside face region. Compute Q through its barycentric
    // coordinates (u, v, w)

    if ((va + vb + vc) < ROOTVSMALL)
    {
        // Degenerate triangle, return the centre because no edge or points are
        // closest
        point nearPt = centre();
        nearType = NONE,
        nearLabel = -1;
        return pointHit(true, nearPt, Foam::mag(nearPt - p), false);
    }

    scalar denom = 1.0/(va + vb + vc);
    scalar v = vb * denom;
    scalar w = vc * denom;

    // = u*a + v*b + w*c, u = va*denom = 1.0 - v - w

    point nearPt = a_ + ab*v + ac*w;
    nearType = NONE,
    nearLabel = -1;
    return pointHit(true, nearPt, Foam::mag(nearPt - p), false);
}


template<class Point, class PointRef>
inline Foam::pointHit Foam::triangle<Point, PointRef>::nearestPoint
(
    const point& p
) const
{
    // Dummy labels
    label nearType = -1;
    label nearLabel = -1;

    return nearestPointClassify(p, nearType, nearLabel);
}


template<class Point, class PointRef>
inline bool Foam::triangle<Point, PointRef>::classify
(
    const point& p,
    label& nearType,
    label& nearLabel
) const
{
    return nearestPointClassify(p, nearType, nearLabel).hit();
}


template<class Point, class PointRef>
inline Foam::pointHit Foam::triangle<Point, PointRef>::nearestPoint
(
    const linePointRef& ln,
    pointHit& lnInfo
) const
{
    vector q = ln.vec();
    pointHit triInfo
    (
        triangle<Point, PointRef>::intersection
        (
            ln.start(),
            q,
            intersection::FULL_RAY
        )
    );

    if (triInfo.hit())
    {
        // Line hits triangle. Find point on line.
        if (triInfo.distance() > 1)
        {
            // Hit beyond endpoint
            lnInfo.setMiss(true);
            lnInfo.setPoint(ln.end());
            scalar dist = Foam::mag(triInfo.hitPoint()-lnInfo.missPoint());
            lnInfo.setDistance(dist);
            triInfo.setMiss(true);
            triInfo.setDistance(dist);
        }
        else if (triInfo.distance() < 0)
        {
            // Hit beyond startpoint
            lnInfo.setMiss(true);
            lnInfo.setPoint(ln.start());
            scalar dist = Foam::mag(triInfo.hitPoint()-lnInfo.missPoint());
            lnInfo.setDistance(dist);
            triInfo.setMiss(true);
            triInfo.setDistance(dist);
        }
        else
        {
            // Hit on line
            lnInfo.setHit();
            lnInfo.setPoint(triInfo.hitPoint());
            lnInfo.setDistance(0.0);
            triInfo.setDistance(0.0);
        }
    }
    else
    {
        // Line skips triangle. See which triangle edge it gets closest to

        point nearestEdgePoint;
        point nearestLinePoint;
        //label minEdgeIndex = 0;
        scalar minDist = ln.nearestDist
        (
            linePointRef(a_, b_),
            nearestLinePoint,
            nearestEdgePoint
        );

        {
            point linePoint;
            point triEdgePoint;
            scalar dist = ln.nearestDist
            (
                linePointRef(b_, c_),
                linePoint,
                triEdgePoint
            );
            if (dist < minDist)
            {
                minDist = dist;
                nearestEdgePoint = triEdgePoint;
                nearestLinePoint = linePoint;
                //minEdgeIndex = 1;
            }
        }

        {
            point linePoint;
            point triEdgePoint;
            scalar dist = ln.nearestDist
            (
                linePointRef(c_, a_),
                linePoint,
                triEdgePoint
            );
            if (dist < minDist)
            {
                minDist = dist;
                nearestEdgePoint = triEdgePoint;
                nearestLinePoint = linePoint;
                //minEdgeIndex = 2;
            }
        }

        lnInfo.setDistance(minDist);
        triInfo.setDistance(minDist);
        triInfo.setMiss(false);
        triInfo.setPoint(nearestEdgePoint);

        // Convert point on line to pointHit
        if (Foam::mag(nearestLinePoint-ln.start()) < SMALL)
        {
            lnInfo.setMiss(true);
            lnInfo.setPoint(ln.start());
        }
        else if (Foam::mag(nearestLinePoint-ln.end()) < SMALL)
        {
            lnInfo.setMiss(true);
            lnInfo.setPoint(ln.end());
        }
        else
        {
            lnInfo.setHit();
            lnInfo.setPoint(nearestLinePoint);
        }
    }
    return triInfo;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    triangle<Point, PointRef>& t
)
{
    is.readBegin("triangle");
    is  >> t.a_ >> t.b_ >> t.c_;
    is.readEnd("triangle");

    is.check("Istream& operator>>(Istream&, triangle&)");
    return is;
}


template<class Point, class PointRef>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const triangle<Point, PointRef>& t
)
{
    os  << nl
        << token::BEGIN_LIST
        << t.a_ << token::SPACE
        << t.b_ << token::SPACE
        << t.c_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
