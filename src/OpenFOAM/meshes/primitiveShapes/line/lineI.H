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

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::line<Point, PointRef>::line(const Point& start, const Point& end)
:
    a_(start),
    b_(end)
{}


template<class Point, class PointRef>
inline Foam::line<Point, PointRef>::line
(
    const UList<Point>& points,
    const FixedList<label, 2>& indices
)
:
    a_(points[indices[0]]),
    b_(points[indices[1]])
{}


template<class Point, class PointRef>
inline Foam::line<Point, PointRef>::line(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline PointRef Foam::line<Point, PointRef>::start() const
{
    return a_;
}

template<class Point, class PointRef>
inline PointRef Foam::line<Point, PointRef>::end() const
{
    return b_;
}


template<class Point, class PointRef>
inline Point Foam::line<Point, PointRef>::centre() const
{
    return 0.5*(a_ + b_);
}


template<class Point, class PointRef>
inline Foam::scalar Foam::line<Point, PointRef>::mag() const
{
    return ::Foam::mag(vec());
}


template<class Point, class PointRef>
inline Point Foam::line<Point, PointRef>::vec() const
{
    return b_ - a_;
}


template<class Point, class PointRef>
Foam::PointHit<Point> Foam::line<Point, PointRef>::nearestDist
(
    const Point& p
) const
{
    Point v = vec();

    Point w(p - a_);

    scalar c1 = v & w;

    if (c1 <= 0)
    {
        return PointHit<Point>(false, a_, Foam::mag(p - a_), true);
    }

    scalar c2 = v & v;

    if (c2 <= c1)
    {
        return PointHit<Point>(false, b_, Foam::mag(p - b_), true);
    }

    scalar b = c1/c2;

    Point pb(a_ + b*v);

    return PointHit<Point>(true, pb, Foam::mag(p - pb), false);
}


template<class Point, class PointRef>
Foam::scalar Foam::line<Point, PointRef>::nearestDist
(
    const line<Point, const Point&>& edge,
    Point& thisPt,
    Point& edgePt
) const
{
    // From Mathworld Line-Line distance/(Gellert et al. 1989, p. 538).
    Point a(end() - start());
    Point b(edge.end() - edge.start());
    Point c(edge.start() - start());

    Point crossab = a ^ b;
    scalar magCrossSqr = magSqr(crossab);

    if (magCrossSqr > vSmall)
    {
        scalar s = ((c ^ b) & crossab)/magCrossSqr;
        scalar t = ((c ^ a) & crossab)/magCrossSqr;

        // Check for end points outside of range 0..1
        if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
        {
            // Both inside range 0..1
            thisPt = start() + a*s;
            edgePt = edge.start() + b*t;
        }
        else
        {
            // Do brute force. Distance of everything to everything.
            // Can quite possibly be improved!

            // From edge endpoints to *this
            PointHit<Point> this0(nearestDist(edge.start()));
            PointHit<Point> this1(nearestDist(edge.end()));
            scalar thisDist = min(this0.distance(), this1.distance());

            // From *this to edge
            PointHit<Point> edge0(edge.nearestDist(start()));
            PointHit<Point> edge1(edge.nearestDist(end()));
            scalar edgeDist = min(edge0.distance(), edge1.distance());

            if (thisDist < edgeDist)
            {
                if (this0.distance() < this1.distance())
                {
                    thisPt = this0.rawPoint();
                    edgePt = edge.start();
                }
                else
                {
                    thisPt = this1.rawPoint();
                    edgePt = edge.end();
                }
            }
            else
            {
                if (edge0.distance() < edge1.distance())
                {
                    thisPt = start();
                    edgePt = edge0.rawPoint();
                }
                else
                {
                    thisPt = end();
                    edgePt = edge1.rawPoint();
                }
            }
        }
    }
    else
    {
        // Parallel lines. Find overlap of both lines by projecting onto
        // direction vector (now equal for both lines).

        scalar edge0 = edge.start() & a;
        scalar edge1 = edge.end() & a;
        bool edgeOrder = edge0 < edge1;

        scalar minEdge = (edgeOrder ? edge0 : edge1);
        scalar maxEdge = (edgeOrder ? edge1 : edge0);
        const Point& minEdgePt = (edgeOrder ? edge.start() : edge.end());
        const Point& maxEdgePt = (edgeOrder ? edge.end() : edge.start());

        scalar this0 = start() & a;
        scalar this1 = end() & a;
        bool thisOrder = this0 < this1;

        scalar minThis = min(this0, this1);
        scalar maxThis = max(this1, this0);
        const Point& minThisPt = (thisOrder ? start() : end());
        const Point& maxThisPt = (thisOrder ? end() : start());

        if (maxEdge < minThis)
        {
            // edge completely below *this
            edgePt = maxEdgePt;
            thisPt = minThisPt;
        }
        else if (maxEdge < maxThis)
        {
            // maxEdge inside interval of *this
            edgePt = maxEdgePt;
            thisPt = nearestDist(edgePt).rawPoint();
        }
        else
        {
            // maxEdge outside. Check if minEdge inside.
            if (minEdge < minThis)
            {
                // Edge completely envelops this. Take any this point and
                // determine nearest on edge.
                thisPt = minThisPt;
                edgePt = edge.nearestDist(thisPt).rawPoint();
            }
            else if (minEdge < maxThis)
            {
                // minEdge inside this interval.
                edgePt = minEdgePt;
                thisPt = nearestDist(edgePt).rawPoint();
            }
            else
            {
                // minEdge outside this interval
                edgePt = minEdgePt;
                thisPt = maxThisPt;
            }
        }
    }

    return Foam::mag(thisPt - edgePt);
}


template<class Point, class PointRef>
bool Foam::line<Point, PointRef>::insideBoundBox(const Point& p) const
{
    if
    (
        ( p.x() < min(a_.x(), b_.x()) || p.x() > max(a_.x(), b_.x()) )
     || ( p.y() < min(a_.y(), b_.y()) || p.y() > max(a_.y(), b_.y()) )
     || ( p.z() < min(a_.z(), b_.z()) || p.z() > max(a_.z(), b_.z()) )
    )
    {
        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    line<Point, PointRef>& l
)
{
    is.readBegin("line");
    is  >> l.a_ >> l.b_;
    is.readEnd("line");

    is.check("Istream& operator>>(Istream&, line&)");
    return is;
}


template<class Point, class PointRef>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const line<Point, PointRef>& l
)
{
    os  << token::BEGIN_LIST
        << l.a_ << token::SPACE
        << l.b_
        << token::END_LIST;
    return os;
}


// ************************************************************************* //
