/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
#include "labelPair.H"
#include "triPointRef.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool Foam::isValidBoundSphere(const Tuple2<point, scalar>& sphere)
{
    return sphere.second() >= 0;
}


template<class PointField>
Foam::Tuple2<Foam::point, Foam::scalar>
Foam::intersectBoundSphere
(
    const PointField& ps,
    const FixedList<label, 4>& pis,
    const label nPs
)
{
    switch (nPs)
    {
        case 0:
            return Tuple2<point, scalar>(point::uniform(NaN), -vGreat);

        case 1:
            return Tuple2<point, scalar>(ps[pis[0]], 0);

        case 2:
            // Return a mid point-centred sphere
            return Tuple2<point, scalar>
            (
                (ps[pis[0]] + ps[pis[1]])/2,
                mag(ps[pis[0]] - ps[pis[1]])/2
            );

        case 3:
        {
            // Return the sphere that intersects the triangle
            return
                triPointRef
                (
                    ps[pis[0]],
                    ps[pis[1]],
                    ps[pis[2]]
                ).circumCircle();
        }

        case 4:
        {
            // Return the sphere that intersects the tetrahedron
            return
                tetPointRef
                (
                    ps[pis[0]],
                    ps[pis[1]],
                    ps[pis[2]],
                    ps[pis[3]]
                ).circumSphere();
        }

        default:
            FatalErrorInFunction
                << "Cannot compute the intersect bounding sphere of more than "
                << "four points" << exit(FatalError);
            return Tuple2<point, scalar>(point::uniform(NaN), NaN);
    }
}


template<class PointField>
Foam::Tuple2<Foam::point, Foam::scalar>
Foam::trivialBoundSphere
(
    const PointField& ps,
    const FixedList<label, 4>& pis,
    const label nPs
)
{
    static const scalar tol = 1 + 2*sqrt(small*rootSmall);

    // Search for spheres that intersect sub-sets of the points that bound all
    // the other points
    switch (nPs)
    {
        case 0:
        case 1:
        case 2:
            break;

        case 3:
        {
            // Test the spheres that intersect the edges
            for (label i = 0; i < 3; ++ i)
            {
                const point& p0 = ps[pis[i]];
                const point& p1 = ps[pis[(i + 1) % 3]];
                const point& pOpp = ps[pis[(i + 2) % 3]];

                const point c = (p0 + p1)/2;
                const scalar rSqr = magSqr(p0 - p1)/4;

                if (magSqr(pOpp - c) <= tol*rSqr)
                {
                    return Tuple2<point, scalar>(c, Foam::sqrt(rSqr));
                }
            }

            break;
        }

        case 4:
        {
            // Test the spheres that intersect the edges
            static const FixedList<labelPair, 6> p01s =
                {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
            static const FixedList<labelPair, 6> pOpp01s =
                {{2, 3}, {1, 3}, {1, 2}, {0, 3}, {0, 2}, {0, 1}};
            for (label i = 0; i < 6; ++ i)
            {
                const point& p0 = ps[pis[p01s[i].first()]];
                const point& p1 = ps[pis[p01s[i].second()]];
                const point& pOpp0 = ps[pis[pOpp01s[i].first()]];
                const point& pOpp1 = ps[pis[pOpp01s[i].second()]];

                const point c = (p0 + p1)/2;
                const scalar rSqr = magSqr(p0 - p1)/4;

                if
                (
                    magSqr(pOpp0 - c) <= tol*rSqr
                 && magSqr(pOpp1 - c) <= tol*rSqr
                )
                {
                    return Tuple2<point, scalar>(c, Foam::sqrt(rSqr));
                }
            }

            // Test the spheres that intersect the triangles
            point minC = point::uniform(NaN);
            scalar minR = vGreat;
            for (label i = 0; i < 4; ++ i)
            {
                const triPointRef tri
                (
                    ps[pis[i]],
                    ps[pis[(i + 1) % 4]],
                    ps[pis[(i + 2) % 4]]
                );
                const point& pOpp = ps[pis[(i + 3) % 4]];

                const Tuple2<point, scalar> circ = tri.circumCircle();
                const point& c = circ.first();
                const scalar r = circ.second();

                if (magSqr(pOpp - c) <= tol*sqr(r) && r < minR)
                {
                    minC = c;
                    minR = r;
                }
            }
            if (minR != vGreat)
            {
                return Tuple2<point, scalar>(minC, minR);
            }

            break;
        }

        default:
            FatalErrorInFunction
                << "Cannot compute the trivial bounding sphere of more than "
                << "four points" << exit(FatalError);
            return Tuple2<point, scalar>(point::uniform(NaN), NaN);
    }

    // Return the sphere that intersects the points
    return intersectBoundSphere(ps, pis, nPs);
}


template<class PointField>
Foam::Tuple2<Foam::point, Foam::scalar> Foam::weizlBoundSphere
(
    const PointField& ps,
    List<label>& pis,
    const label nPs,
    FixedList<label, 4>& boundaryPis,
    const label nBoundaryPs
)
{
    static const scalar tol = 1 + 2*sqrt(small*rootSmall);

    Tuple2<point, scalar> sphere;

    if (nBoundaryPs != 0)
    {
        sphere = intersectBoundSphere(ps, boundaryPis, nBoundaryPs);
    }

    if (nBoundaryPs == 4)
    {
        return sphere;
    }

    for (label i = 0; i < nPs; ++ i)
    {
        if
        (
            (nBoundaryPs == 0 && i == 0)
         || (magSqr(ps[pis[i]] - sphere.first()) > tol*sqr(sphere.second()))
        )
        {
            boundaryPis[nBoundaryPs] = pis[i];

            sphere = weizlBoundSphere(ps, pis, i, boundaryPis, nBoundaryPs + 1);

            // Move the limiting point to the start of the list so that the
            // sphere grows as quickly as possible in the recursive calls
            Swap(pis[0], pis[i]);
        }
    }

    return sphere;
}


template<class PointField>
Foam::Tuple2<Foam::point, Foam::scalar>
Foam::boundSphere(const PointField& ps, Random& rndGen)
{
    if (ps.size() <= 4)
    {
        static const FixedList<label, 4> pis({0, 1, 2, 3});
        return trivialBoundSphere(ps, pis, ps.size());
    }
    else
    {
        labelList pis(identityMap(ps.size()));
        rndGen.permute(pis);
        FixedList<label, 4> boundaryPis({-1, -1, -1, -1});
        return weizlBoundSphere(ps, pis, ps.size(), boundaryPis, 0);
    }
}


template<class PointField>
Foam::Tuple2<Foam::point, Foam::scalar>
Foam::boundSphere(const PointField& ps)
{
    if (ps.size() <= 4)
    {
        static const FixedList<label, 4> pis({0, 1, 2, 3});
        return trivialBoundSphere(ps, pis, ps.size());
    }
    else
    {
        labelList pis(identityMap(ps.size()));
        Random(0).permute(pis);
        FixedList<label, 4> boundaryPis({-1, -1, -1, -1});
        return weizlBoundSphere(ps, pis, ps.size(), boundaryPis, 0);
    }
}


// ************************************************************************* //

