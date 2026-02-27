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
#include "labelPair.H"
#include "triPointRef.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class PointList>
Foam::boundSphere Foam::boundSphere::implementation::intersect
(
    const PointList& ps,
    const FixedList<label, 4>& pis,
    const label nPs
)
{
    switch (nPs)
    {
        case 0:
            return boundSphere(point::uniform(NaN), -vGreat);

        case 1:
            return boundSphere(ps[pis[0]], 0);

        case 2:
        {
            // Return a sphere centred on the mid-point between the two points
            // and with radius half the distance between the points
            const point c = (ps[pis[0]] + ps[pis[1]])/2;
            return boundSphere
            (
                c,
                max
                (
                    magSqr(ps[pis[0]] - c),
                    magSqr(ps[pis[1]] - c)
                )
            );
        }

        case 3:
        {
            // If the triangle is not degenerate then return the sphere for
            // which a great circle is the triangles circum-circle
            const triPointRef tri
            (
                ps[pis[0]],
                ps[pis[1]],
                ps[pis[2]]
            );

            const boundSphere s = crSqr(tri.circumCircleSqr());

            // If the quality is less than small then return an invalid sphere
            static const scalar ratio = small*3.0/4.0*sqrt(scalar(3));
            if (!s.valid() || magSqr(tri.area()) < sqr(ratio)*sqr(s.rSqr_))
            {
                return boundSphere();
            }

            // Re-calculate the radius to make sure that the sphere contains
            // all the points. This makes Welzl's algorithm robust to the
            // presence of coincident points.
            return
                boundSphere
                (
                    s.c(),
                    max
                    (
                        magSqr(tri.a() - s.c()),
                        max
                        (
                            magSqr(tri.b() - s.c()),
                            magSqr(tri.c() - s.c())
                        )
                    )
                );
        }

        case 4:
        {
            // If the tetrahedron is not degenerate then return the
            // cicrum-sphere
            const tetPointRef tet
            (
                ps[pis[0]],
                ps[pis[1]],
                ps[pis[2]],
                ps[pis[3]]
            );

            const boundSphere s = crSqr(tet.circumSphereSqr());

            // If the quality is less than small then return an invalid sphere
            static const scalar ratio = small*8.0/27.0*sqrt(scalar(3));
            if (!s.valid() || sqr(tet.mag()) < sqr(ratio)*pow3(s.rSqr_))
            {
                return boundSphere();
            }

            // Re-calculate the radius to make sure that the sphere contains
            // all the points. This makes Welzl's algorithm robust to the
            // presence of coincident points.
            return
                boundSphere
                (
                    s.c(),
                    max
                    (
                        magSqr(tet.a() - s.c()),
                        max
                        (
                            magSqr(tet.b() - s.c()),
                            max
                            (
                                magSqr(tet.c() - s.c()),
                                magSqr(tet.d() - s.c())
                            )
                        )
                    )
                );
        }

        default:
            FatalErrorInFunction
                << "Cannot compute the intersect bounding sphere of more than "
                << "four points" << exit(FatalError);
            return boundSphere(point::uniform(NaN), NaN);
    }
}


template<class PointList>
Foam::boundSphere Foam::boundSphere::implementation::trivial
(
    const PointList& ps,
    const FixedList<label, 4>& pis,
    const label nPs,
    FixedList<label, 4>& boundaryPis
)
{
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
                const label pi0 = pis[i];
                const label pi1 = pis[(i + 1) % 3];
                const label piOpp = pis[(i + 2) % 3];

                const boundSphere s = intersect(ps, {pi0, pi1, -1, -1}, 2);

                if (s.contains(ps[piOpp]))
                {
                    boundaryPis[0] = pis[i];
                    boundaryPis[1] = pis[(i + 1) % 3];
                    return s;
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
                const label pi0 = pis[p01s[i].first()];
                const label pi1 = pis[p01s[i].second()];
                const label piOpp0 = pis[pOpp01s[i].first()];
                const label piOpp1 = pis[pOpp01s[i].second()];

                const boundSphere s = intersect(ps, {pi0, pi1, -1, -1}, 2);

                if (s.contains(ps[piOpp0]) && s.contains(ps[piOpp1]))
                {
                    boundaryPis[0] = pis[p01s[i].first()];
                    boundaryPis[1] = pis[p01s[i].second()];
                    return s;
                }
            }

            // Test the spheres that intersect the triangles
            boundSphere minS(point::uniform(NaN), vGreat);
            label minSi = -1;
            for (label i = 0; i < 4; ++ i)
            {
                const label pi0 = pis[i];
                const label pi1 = pis[(i + 1) % 4];
                const label pi2 = pis[(i + 2) % 4];
                const label piOpp = pis[(i + 3) % 4];

                const boundSphere s = intersect(ps, {pi0, pi1, pi2, -1}, 3);

                if (s.contains(ps[piOpp]) && s.rSqr_ < minS.rSqr_)
                {
                    minS = s;
                    minSi = i;
                }
            }

            if (minSi != -1)
            {
                boundaryPis[0] = pis[minSi];
                boundaryPis[1] = pis[(minSi + 1) % 4];
                boundaryPis[2] = pis[(minSi + 2) % 4];
                return minS;
            }

            break;
        }

        default:
            FatalErrorInFunction
                << "Cannot compute the trivial bounding sphere of more than "
                << "four points" << exit(FatalError);
            return boundSphere();
    }

    // Set all the points as on the boundary
    for (label i = 0; i < nPs; ++ i)
    {
        boundaryPis[i] = pis[i];
    }
    sortBoundaryPis(boundaryPis);

    // Return the sphere that intersects all the points
    return intersect(ps, pis, nPs);
}


template<class PointList>
Foam::boundSphere Foam::boundSphere::implementation::bruteForce
(
    const PointList& ps,
    FixedList<label, 4>& boundaryPis
)
{
    // Update the currently stored bounding sphere if the supplied set of
    // points defines a smaller bounding sphere
    auto update = []
    (
        const PointList& ps,
        const FixedList<label, 4>& pis,
        const label nPs,
        boundSphere& s,
        FixedList<label, 4>& boundaryPis
    )
    {
        // Construct the intersect sphere for these points
        const boundSphere sStar = intersect(ps, pis, nPs);

        // Determine if the sphere bounds all the points, excluding any
        // actually on the boundary so as to prevent round-off error issues
        bool bounds = true;
        for (label pi = 0; bounds && pi < ps.size(); ++ pi)
        {
            label i = 0;
            for (; i < nPs; ++ i) if (pi == pis[i]) break;
            if (i != nPs) continue;

            if (!sStar.contains(ps[pi])) bounds = false;
        }

        // If this sphere bounds everything and is smaller than the current
        // sphere then update
        if (bounds && sStar.rSqr() < s.rSqr())
        {
            s = sStar;
            boundaryPis = pis;
        }
    };

    // Return an invalid sphere if there are no points
    if (ps.size() == 0) return boundSphere();

    // Initialise a very large sphere
    boundSphere s = boundSphere::crSqr(point::uniform(NaN), vGreat);
    boundaryPis = -1;

    // Search through all combinations of points to find the minimum sphere
    FixedList<label, 4> pis({-1, -1, -1, -1});
    for (pis[0] = 0; pis[0] < ps.size(); ++ pis[0])
    {
        update(ps, pis, 1, s, boundaryPis);
        for (pis[1] = pis[0] + 1; pis[1] < ps.size(); ++ pis[1])
        {
            update(ps, pis, 2, s, boundaryPis);
            for (pis[2] = pis[1] + 1; pis[2] < ps.size(); ++ pis[2])
            {
                update(ps, pis, 3, s, boundaryPis);
                for (pis[3] = pis[2] + 1; pis[3] < ps.size(); ++ pis[3])
                {
                    update(ps, pis, 4, s, boundaryPis);
                }
                pis[3] = -1;
            }
            pis[2] = -1;
        }
        pis[1] = -1;
    }

    return s;
}


template<class PointList, class DLPermutation>
Foam::boundSphere Foam::boundSphere::implementation::Welzl
(
    const PointList& ps,
    DLPermutation& pis,
    const typename DLPermutation::const_iterator& pisEnd,
    FixedList<label, 4>& boundaryPis,
    const label nBoundaryPs
)
{
    boundSphere s = intersect(ps, boundaryPis, nBoundaryPs);

    if (pisEnd == pis.rend() || nBoundaryPs == 4) return s;

    for (label boundaryi = nBoundaryPs; boundaryi < 4; ++ boundaryi)
    {
        boundaryPis[boundaryi] = -1;
    }

    for
    (
        typename DLPermutation::const_iterator piIter = pis.begin();
        piIter != pisEnd;
        ++ piIter
    )
    {
        if (s.contains(ps[*piIter])) continue;

        boundaryPis[nBoundaryPs] = *piIter;

        s = Welzl(ps, pis, piIter, boundaryPis, nBoundaryPs + 1);

        // Move this point up to the start
        if (piIter != pis.begin())
        {
            const label i = *piIter;
            piIter = piIter.previous();
            pis.raise(i);
        }
    }

    return s;
}


template<class PointList>
Foam::boundSphere Foam::boundSphere::implementation::local
(
    const PointList& ps,
    randomGenerator& rndGen,
    FixedList<label, 4>& boundaryPis
)
{
    if (ps.size() <= 4)
    {
        static const FixedList<label, 4> pis({0, 1, 2, 3});
        boundaryPis = -1;
        return trivial(ps, pis, ps.size(), boundaryPis);
    }
    else
    {
        dLPermutation pis(ps.size(), rndGen);
        boundaryPis = -1;
        boundSphere s = Welzl(ps, pis, pis.end(), boundaryPis, 0);
        sortBoundaryPis(boundaryPis);
        return s;
    }
}


template<class PointList>
Foam::boundSphere Foam::boundSphere::implementation::global
(
    const PointList& ps,
    randomGenerator& rndGen,
    FixedList<RemoteData<point>, 4>& boundaryPs,
    const bool strict
)
{
    // Run the local algorithm if serial
    if (!Pstream::parRun())
    {
        FixedList<label, 4> boundaryPis({-1, -1, -1, -1});

        const boundSphere s = local(ps, rndGen, boundaryPis);

        forAll(boundaryPis, boundaryi)
        {
            if (boundaryPis[boundaryi] != -1)
            {
                boundaryPs[boundaryi] =
                    RemoteData<point>
                    (
                        Pstream::myProcNo(),
                        boundaryPis[boundaryi],
                        ps[boundaryPis[boundaryi]]
                    );
            }
        }

        return s;
    }

    LocalAndRemotePoints<PointList> larps(ps);

    dynamicDLPermutation pis(ps.size(), rndGen);
    FixedList<label, 4> boundaryPis({-1, -1, -1, -1});

    while (true)
    {
        // Calculate the sphere for the current set of points
        const boundSphere s = Welzl(larps, pis, pis.end(), boundaryPis, 0);

        // Communicate the boundary points
        List<FixedList<RemoteData<point>, 4>> procBoundaryPs(Pstream::nProcs());
        {
            label boundaryi = 0;
            for (; boundaryi < 4 && boundaryPis[boundaryi] != -1; ++ boundaryi)
            {
                procBoundaryPs[Pstream::myProcNo()][boundaryi] =
                    larps(boundaryPis[boundaryi]);
            }

            sortBoundaryPis(procBoundaryPs[Pstream::myProcNo()]);

            Pstream::gatherList(procBoundaryPs);
            Pstream::scatterList(procBoundaryPs);
        }

        // Determine if all boundary points are the same on all processors. If
        // they are then the sphere is complete.
        bool complete = true;
        forAll(procBoundaryPs, proci)
        {
            if (proci == Pstream::myProcNo()) continue;

            for (label boundaryi = 0; boundaryi < 4; ++ boundaryi)
            {
                const remote& r = procBoundaryPs[proci][boundaryi];

                if (r != procBoundaryPs[Pstream::myProcNo()][boundaryi])
                {
                    complete = false;
                }

                if (!complete) break;
            }

            if (!complete) break;
        }

        if (complete)
        {
            boundaryPs = procBoundaryPs[Pstream::myProcNo()];
            return s;
        }

        // Add remote boundary points to the local lists and iterate
        complete = true;
        forAll(procBoundaryPs, proci)
        {
            if (proci == Pstream::myProcNo()) continue;

            for (label boundaryi = 0; boundaryi < 4; ++ boundaryi)
            {
                const RemoteData<point>& rp = procBoundaryPs[proci][boundaryi];

                if (rp.proci == -1) break;
                if (larps.found(rp)) continue;

                if (s.contains(rp.data)) continue;

                larps.insertNoCheck(rp);
                pis.prepend();

                complete = false;
            }
        }

        // If nothing was added on any processor then the state has stopped
        // changing. If we are being strict, then this should never happen, so
        // raise an error. However, this can also happen due to ambiguous cases
        // for which multiple solutions are valid to within round-off error
        // (e.g., a cube), in which case the different ordering on different
        // processes can result in a different (but still valid) solution. So,
        // if we are not being strict then just return the sphere and the
        // boundary points from the master process.
        if (returnReduce(complete, andOp<bool>()))
        {
            if (strict)
            {
                FatalErrorInFunction
                    << "Global bound sphere iteration did not converge"
                    << exit(FatalError);
            }

            boundaryPs = procBoundaryPs[Pstream::master()];
            return s;
        }
    }

    return boundSphere();
}


#define InstantiateBoundSphere(PointList)                                      \
                                                                               \
    template Foam::boundSphere Foam::boundSphere::implementation::intersect    \
    (                                                                          \
        const PointList& ps,                                                   \
        const FixedList<label, 4>& pis,                                        \
        const label nPs                                                        \
    );                                                                         \
                                                                               \
    template Foam::boundSphere Foam::boundSphere::implementation::trivial      \
    (                                                                          \
        const PointList& ps,                                                   \
        const FixedList<label, 4>& pis,                                        \
        const label nPs,                                                       \
        FixedList<label, 4>& boundaryPis                                       \
    );                                                                         \
                                                                               \
    template Foam::boundSphere Foam::boundSphere::implementation::bruteForce   \
    (                                                                          \
        const PointList& ps,                                                   \
        FixedList<label, 4>& boundaryPis                                       \
    );                                                                         \
                                                                               \
    template Foam::boundSphere Foam::boundSphere::implementation::local        \
    (                                                                          \
        const PointList& ps,                                                   \
        randomGenerator& rndGen,                                               \
        FixedList<label, 4>& boundaryPis                                       \
    );                                                                         \
                                                                               \
    template Foam::boundSphere Foam::boundSphere::implementation::global       \
    (                                                                          \
        const PointList& ps,                                                   \
        randomGenerator& rndGen,                                               \
        FixedList<RemoteData<point>, 4>& boundaryPis,                          \
        const bool strict                                                      \
    );


InstantiateBoundSphere(UList<point>);
InstantiateBoundSphere(UIndirectList<point>);


#undef InstantiateBoundSphere


// ************************************************************************* //

