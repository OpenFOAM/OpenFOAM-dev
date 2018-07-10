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

#include "tetrahedron.H"
#include "triPointRef.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
Foam::pointHit Foam::tetrahedron<Point, PointRef>::containmentSphere
(
    const scalar tol
) const
{
    // (Probably very inefficient) minimum containment sphere calculation.
    // From http://www.imr.sandia.gov/papers/imr11/shewchuk2.pdf:
    // Sphere ctr is smallest one of
    // - tet circumcentre
    // - triangle circumcentre
    // - edge mids

    const scalar fac = 1 + tol;

    // Halve order of tolerance for comparisons of sqr.
    const scalar facSqr = Foam::sqrt(fac);


    // 1. Circumcentre itself.

    pointHit pHit(circumCentre());
    pHit.setHit();
    scalar minRadiusSqr = magSqr(pHit.rawPoint() - a_);


    // 2. Try circumcentre of tet triangles. Create circumcircle for triFace and
    // check if 4th point is inside.

    {
        point ctr = triPointRef(a_, b_, c_).circumCentre();
        scalar radiusSqr = magSqr(ctr - a_);

        if
        (
            radiusSqr < minRadiusSqr
         && Foam::magSqr(d_-ctr) <= facSqr*radiusSqr
        )
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }
    {
        point ctr = triPointRef(a_, b_, d_).circumCentre();
        scalar radiusSqr = magSqr(ctr - a_);

        if
        (
            radiusSqr < minRadiusSqr
         && Foam::magSqr(c_-ctr) <= facSqr*radiusSqr
        )
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }
    {
        point ctr = triPointRef(a_, c_, d_).circumCentre();
        scalar radiusSqr = magSqr(ctr - a_);

        if
        (
            radiusSqr < minRadiusSqr
         && Foam::magSqr(b_-ctr) <= facSqr*radiusSqr
        )
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }
    {
        point ctr = triPointRef(b_, c_, d_).circumCentre();
        scalar radiusSqr = magSqr(ctr - b_);

        if
        (
            radiusSqr < minRadiusSqr
         && Foam::magSqr(a_-ctr) <= facSqr*radiusSqr
        )
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }


    // 3. Try midpoints of edges

    // mid of edge A-B
    {
        point ctr = 0.5*(a_ + b_);
        scalar radiusSqr = magSqr(a_ - ctr);
        scalar testRadSrq = facSqr*radiusSqr;

        if
        (
            radiusSqr < minRadiusSqr
         && magSqr(c_-ctr) <= testRadSrq
         && magSqr(d_-ctr) <= testRadSrq)
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }

    // mid of edge A-C
    {
        point ctr = 0.5*(a_ + c_);
        scalar radiusSqr = magSqr(a_ - ctr);
        scalar testRadSrq = facSqr*radiusSqr;

        if
        (
            radiusSqr < minRadiusSqr
         && magSqr(b_-ctr) <= testRadSrq
         && magSqr(d_-ctr) <= testRadSrq
        )
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }

    // mid of edge A-D
    {
        point ctr = 0.5*(a_ + d_);
        scalar radiusSqr = magSqr(a_ - ctr);
        scalar testRadSrq = facSqr*radiusSqr;

        if
        (
            radiusSqr < minRadiusSqr
         && magSqr(b_-ctr) <= testRadSrq
         && magSqr(c_-ctr) <= testRadSrq
        )
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }

    // mid of edge B-C
    {
        point ctr = 0.5*(b_ + c_);
        scalar radiusSqr = magSqr(b_ - ctr);
        scalar testRadSrq = facSqr*radiusSqr;

        if
        (
            radiusSqr < minRadiusSqr
         && magSqr(a_-ctr) <= testRadSrq
         && magSqr(d_-ctr) <= testRadSrq
        )
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }

    // mid of edge B-D
    {
        point ctr = 0.5*(b_ + d_);
        scalar radiusSqr = magSqr(b_ - ctr);
        scalar testRadSrq = facSqr*radiusSqr;

        if
        (
            radiusSqr < minRadiusSqr
         && magSqr(a_-ctr) <= testRadSrq
         && magSqr(c_-ctr) <= testRadSrq)
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }

    // mid of edge C-D
    {
        point ctr = 0.5*(c_ + d_);
        scalar radiusSqr = magSqr(c_ - ctr);
        scalar testRadSrq = facSqr*radiusSqr;

        if
        (
            radiusSqr < minRadiusSqr
         && magSqr(a_-ctr) <= testRadSrq
         && magSqr(b_-ctr) <= testRadSrq
        )
        {
            pHit.setMiss(false);
            pHit.setPoint(ctr);
            minRadiusSqr = radiusSqr;
        }
    }


    pHit.setDistance(sqrt(minRadiusSqr));

    return pHit;
}


template<class Point, class PointRef>
void Foam::tetrahedron<Point, PointRef>::gradNiSquared
(
    scalarField& buffer
) const
{
    // Change of sign between face area vector and gradient
    // does not matter because of square

    // Warning: Added a mag to produce positive coefficients even if
    // the tetrahedron is twisted inside out.  This is pretty
    // dangerous, but essential for mesh motion.
    scalar magVol = Foam::mag(mag());

    buffer[0] = (1.0/9.0)*magSqr(Sa())/magVol;
    buffer[1] = (1.0/9.0)*magSqr(Sb())/magVol;
    buffer[2] = (1.0/9.0)*magSqr(Sc())/magVol;
    buffer[3] = (1.0/9.0)*magSqr(Sd())/magVol;
}


template<class Point, class PointRef>
void Foam::tetrahedron<Point, PointRef>::gradNiDotGradNj
(
    scalarField& buffer
) const
{
    // Warning. Ordering of edges needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell

    // Warning: Added a mag to produce positive coefficients even if
    // the tetrahedron is twisted inside out.  This is pretty
    // dangerous, but essential for mesh motion.

    // Double change of sign between face area vector and gradient

    scalar magVol = Foam::mag(mag());
    vector sa = Sa();
    vector sb = Sb();
    vector sc = Sc();
    vector sd = Sd();

    buffer[0] = (1.0/9.0)*(sa & sb)/magVol;
    buffer[1] = (1.0/9.0)*(sa & sc)/magVol;
    buffer[2] = (1.0/9.0)*(sa & sd)/magVol;
    buffer[3] = (1.0/9.0)*(sd & sb)/magVol;
    buffer[4] = (1.0/9.0)*(sb & sc)/magVol;
    buffer[5] = (1.0/9.0)*(sd & sc)/magVol;
}


template<class Point, class PointRef>
void Foam::tetrahedron<Point, PointRef>::gradNiGradNi
(
    tensorField& buffer
) const
{
    // Change of sign between face area vector and gradient
    // does not matter because of square

    scalar magVol = Foam::mag(mag());

    buffer[0] = (1.0/9.0)*sqr(Sa())/magVol;
    buffer[1] = (1.0/9.0)*sqr(Sb())/magVol;
    buffer[2] = (1.0/9.0)*sqr(Sc())/magVol;
    buffer[3] = (1.0/9.0)*sqr(Sd())/magVol;
}


template<class Point, class PointRef>
void Foam::tetrahedron<Point, PointRef>::gradNiGradNj
(
    tensorField& buffer
) const
{
    // Warning. Ordering of edges needs to be the same for a tetrahedron
    // class, a tetrahedron cell shape model and a tetCell

    // Double change of sign between face area vector and gradient

    scalar magVol = Foam::mag(mag());
    vector sa = Sa();
    vector sb = Sb();
    vector sc = Sc();
    vector sd = Sd();

    buffer[0] = (1.0/9.0)*(sa * sb)/magVol;
    buffer[1] = (1.0/9.0)*(sa * sc)/magVol;
    buffer[2] = (1.0/9.0)*(sa * sd)/magVol;
    buffer[3] = (1.0/9.0)*(sd * sb)/magVol;
    buffer[4] = (1.0/9.0)*(sb * sc)/magVol;
    buffer[5] = (1.0/9.0)*(sd * sc)/magVol;
}


template<class Point, class PointRef>
Foam::boundBox Foam::tetrahedron<Point, PointRef>::bounds() const
{
    return
        boundBox
        (
            min(a(), min(b(), min(c(), d()))),
            max(a(), max(b(), max(c(), d())))
        );
}


// ************************************************************************* //
