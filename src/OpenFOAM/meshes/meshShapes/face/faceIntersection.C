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

#include "face.H"
#include "pointHit.H"
#include "triPointRef.H"
#include "line.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::pointHit Foam::face::ray
(
    const point& p,
    const vector& n,
    const pointField& meshPoints,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    // Return potential intersection with face with a ray starting
    // at p, direction n (does not need to be normalized)
    // Does face-center decomposition and returns triangle intersection
    // point closest to p.

    // In case of miss the point is the nearest point intersection of the
    // face plane and the ray and the distance is the distance between the
    // intersection point and the nearest point on the face

    // If the face is a triangle, do a direct calculation
    if (size() == 3)
    {
        return triPointRef
        (
            meshPoints[operator[](0)],
            meshPoints[operator[](1)],
            meshPoints[operator[](2)]
        ).ray(p, n, alg, dir);
    }

    point ctr = Foam::average(points(meshPoints));

    scalar nearestHitDist = great;

    scalar nearestMissDist = great;
    bool eligible = false;

    // Initialize to miss, distance = great
    pointHit nearest(p);

    const labelList& f = *this;

    label nPoints = size();

    point nextPoint = ctr;

    for (label pI = 0; pI < nPoints; pI++)
    {
        nextPoint = meshPoints[f[fcIndex(pI)]];

        // Note: for best accuracy, centre point always comes last
        //
        pointHit curHit = triPointRef
        (
            meshPoints[f[pI]],
            nextPoint,
            ctr
        ).ray(p, n, alg, dir);

        if (curHit.hit())
        {
            if (Foam::mag(curHit.distance()) < Foam::mag(nearestHitDist))
            {
                nearestHitDist = curHit.distance();
                nearest.setHit();
                nearest.setPoint(curHit.hitPoint());
            }
        }
        else if (!nearest.hit())
        {
            // Miss and no hit yet. Update miss statistics.
            if (curHit.eligibleMiss())
            {
                eligible = true;

                // Miss distance is the distance between the plane intersection
                // point and the nearest point of the triangle
                scalar missDist =
                    Foam::mag
                    (
                        p + curHit.distance()*n
                      - curHit.missPoint()
                    );

                if (missDist < nearestMissDist)
                {
                    nearestMissDist = missDist;
                    nearest.setDistance(curHit.distance());
                    nearest.setPoint(curHit.missPoint());
                }
            }
        }
    }

    if (nearest.hit())
    {
        nearest.setDistance(nearestHitDist);
    }
    else
    {
        // Haven't hit a single face triangle
        nearest.setMiss(eligible);
    }

    return nearest;
}


Foam::pointHit Foam::face::intersection
(
    const point& p,
    const vector& q,
    const point& ctr,
    const pointField& meshPoints,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    // If the face is a triangle, do a direct calculation
    if (size() == 3)
    {
        return triPointRef
        (
            meshPoints[operator[](0)],
            meshPoints[operator[](1)],
            meshPoints[operator[](2)]
        ).intersection(p, q, alg, tol);
    }

    scalar nearestHitDist = vGreat;

    // Initialize to miss, distance = great
    pointHit nearest(p);

    const labelList& f = *this;

    forAll(f, pI)
    {
        // Note: for best accuracy, centre point always comes last
        pointHit curHit = triPointRef
        (
            meshPoints[f[pI]],
            meshPoints[f[fcIndex(pI)]],
            ctr
        ).intersection(p, q, alg, tol);

        if (curHit.hit())
        {
            if (Foam::mag(curHit.distance()) < Foam::mag(nearestHitDist))
            {
                nearestHitDist = curHit.distance();
                nearest.setHit();
                nearest.setPoint(curHit.hitPoint());
            }
        }
    }

    if (nearest.hit())
    {
        nearest.setDistance(nearestHitDist);
    }

    return nearest;
}


Foam::pointHit Foam::face::nearestPoint
(
    const point& p,
    const pointField& meshPoints
) const
{
    // Dummy labels
    label nearType = -1;
    label nearLabel = -1;

    return nearestPointClassify(p, meshPoints, nearType, nearLabel);
}


Foam::pointHit Foam::face::nearestPointClassify
(
    const point& p,
    const pointField& meshPoints,
    label& nearType,
    label& nearLabel
) const
{
    // If the face is a triangle, do a direct calculation
    if (size() == 3)
    {
        return triPointRef
        (
            meshPoints[operator[](0)],
            meshPoints[operator[](1)],
            meshPoints[operator[](2)]
        ).nearestPointClassify(p, nearType, nearLabel);
    }

    const face& f = *this;
    point ctr = centre(meshPoints);

    // Initialize to miss, distance=great
    pointHit nearest(p);

    nearType = -1;
    nearLabel = -1;

    label nPoints = f.size();

    point nextPoint = ctr;

    for (label pI = 0; pI < nPoints; pI++)
    {
        nextPoint = meshPoints[f[fcIndex(pI)]];

        label tmpNearType = -1;
        label tmpNearLabel = -1;

        // Note: for best accuracy, centre point always comes last
        triPointRef tri
        (
            meshPoints[f[pI]],
            nextPoint,
            ctr
        );

        pointHit curHit = tri.nearestPointClassify
        (
            p,
            tmpNearType,
            tmpNearLabel
        );

        if (Foam::mag(curHit.distance()) < Foam::mag(nearest.distance()))
        {
            nearest.setDistance(curHit.distance());

            // Assume at first that the near type is NONE on the
            // triangle (i.e. on the face of the triangle) then it is
            // therefore also for the face.

            nearType = NONE;

            if (tmpNearType == triPointRef::EDGE && tmpNearLabel == 0)
            {
                // If the triangle edge label is 0, then this is also
                // an edge of the face, if not, it is on the face

                nearType = EDGE;

                nearLabel = pI;
            }
            else if (tmpNearType == triPointRef::POINT && tmpNearLabel < 2)
            {
                // If the triangle point label is 0 or 1, then this is
                // also a point of the face, if not, it is on the face

                nearType = POINT;

                nearLabel = pI + tmpNearLabel;
            }

            if (curHit.hit())
            {
                nearest.setHit();
                nearest.setPoint(curHit.hitPoint());
            }
            else
            {
                // In nearest point, miss is always eligible
                nearest.setMiss(true);
                nearest.setPoint(curHit.missPoint());
            }
        }
    }

    return nearest;
}


// ************************************************************************* //
