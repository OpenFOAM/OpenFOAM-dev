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

#include "faceAreaIntersect.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<faceAreaIntersect::triangulationMode, 2>::names[] =
    {
        "fan",
        "mesh"
    };
}

const Foam::NamedEnum<Foam::faceAreaIntersect::triangulationMode, 2>
    Foam::faceAreaIntersect::triangulationModeNames_;

Foam::scalar Foam::faceAreaIntersect::tol = 1e-6;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::faceAreaIntersect::triSliceWithPlane
(
    const triPoints& tri,
    const plane& p,
    FixedList<triPoints, 10>& tris,
    label& nTris,
    const scalar len
)
{
    // distance to cutting plane
    FixedList<scalar, 3> d;

    // determine how many of the points are above the cutting plane
    label nCoPlanar = 0;
    label nPos = 0;
    label posI = -1;
    label negI = -1;
    label copI = -1;
    forAll(tri, i)
    {
        d[i] = ((tri[i] - p.refPoint()) & p.normal());

        if (mag(d[i]) < tol*len)
        {
            nCoPlanar++;
            copI = i;
            d[i] = 0.0;
        }
        else
        {
            if (d[i] > 0)
            {
                nPos++;
                posI = i;
            }
            else
            {
                negI = i;
            }
        }
    }


    // Determine triangle area contribution

    if
    (
        (nPos == 3)
     || ((nPos == 2) && (nCoPlanar == 1))
     || ((nPos == 1) && (nCoPlanar == 2))
    )
    {
        /*
                /\          _____
               /  \         \   /          /\
              /____\         \ /          /  \
            __________    ____v____    __/____\__

            all points above cutting plane
            - add complete triangle to list
        */

        tris[nTris++] = tri;
    }
    else if ((nPos == 2) && (nCoPlanar == 0))
    {
        /*
            i1________i2
              \      /
             --\----/--
                \  /
                 \/
                 i0

            2 points above plane, 1 below
            - resulting quad above plane split into 2 triangles
            - forget triangle below plane
        */

        // point under the plane
        label i0 = negI;

        // indices of remaining points
        label i1 = d.fcIndex(i0);
        label i2 = d.fcIndex(i1);

        // determine the two intersection points
        point p01 = planeIntersection(d, tri, i0, i1);
        point p02 = planeIntersection(d, tri, i0, i2);

        // forget triangle below plane
        // - decompose quad above plane into 2 triangles and add to list
        setTriPoints(tri[i1], tri[i2], p02, nTris, tris);
        setTriPoints(tri[i1], p02, p01, nTris, tris);
    }
    else if (nPos == 1)
    {
        // point above the plane
        label i0 = posI;

        if (nCoPlanar == 0)
        {
            /*
                     i0
                     /\
                    /  \
                 --/----\--
                  /______\
                i2        i1

                1 point above plane, 2 below
                - keep triangle above intersection plane
                - forget quad below plane
            */

            // indices of remaining points
            label i1 = d.fcIndex(i0);
            label i2 = d.fcIndex(i1);

            // determine the two intersection points
            point p01 = planeIntersection(d, tri, i1, i0);
            point p02 = planeIntersection(d, tri, i2, i0);

            // add triangle above plane to list
            setTriPoints(tri[i0], p01, p02, nTris, tris);
        }
        else
        {
            /*
                  i0
                  |\
                  | \
                __|__\_i2_
                  |  /
                  | /
                  |/
                  i1

                1 point above plane, 1 on plane, 1 below
                - keep triangle above intersection plane
            */

            // point indices
            label i1 = negI;
            label i2 = copI;

            // determine the intersection point
            point p01 = planeIntersection(d, tri, i1, i0);

            // add triangle above plane to list - clockwise points
            if (d.fcIndex(i0) == i1)
            {
                setTriPoints(tri[i0], p01, tri[i2], nTris, tris);
            }
            else
            {
                setTriPoints(tri[i0], tri[i2], p01, nTris, tris);
            }
        }
    }
    else
    {
        /*
            _________    __________    ___________
                             /\          \    /
               /\           /  \          \  /
              /  \         /____\          \/
             /____\

            all points below cutting plane - forget
        */
    }
}


Foam::scalar Foam::faceAreaIntersect::triangleIntersect
(
    const triPoints& src,
    const triPoints& tgt,
    const vector& n
)
{
    // Work storage
    FixedList<triPoints, 10> workTris1;
    label nWorkTris1 = 0;

    FixedList<triPoints, 10> workTris2;
    label nWorkTris2 = 0;

    // cut source triangle with all inwards pointing faces of target triangle
    // - triangles in workTris1 are inside target triangle

    scalar t = sqrt(triArea(src));

    // edge 0
    {
        // cut triangle src with plane and put resulting sub-triangles in
        // workTris1 list

        scalar s = mag(tgt[1] - tgt[0]);
        plane pl0(tgt[0], tgt[1], tgt[1] + s*n);
        triSliceWithPlane(src, pl0, workTris1, nWorkTris1, t);
    }

    if (nWorkTris1 == 0)
    {
        return 0.0;
    }

    // edge1
    {
        // cut workTris1 with plane and put resulting sub-triangles in
        // workTris2 list (re-use tris storage)

        scalar s = mag(tgt[2] - tgt[1]);
        plane pl1(tgt[1], tgt[2], tgt[2] + s*n);

        nWorkTris2 = 0;

        for (label i = 0; i < nWorkTris1; i++)
        {
            triSliceWithPlane(workTris1[i], pl1, workTris2, nWorkTris2, t);
        }

        if (nWorkTris2 == 0)
        {
            return 0.0;
        }
    }

    // edge2
    {
        // cut workTris2 with plane and put resulting sub-triangles in
        // workTris1 list (re-use workTris1 storage)

        scalar s = mag(tgt[2] - tgt[0]);
        plane pl2(tgt[2], tgt[0], tgt[0] + s*n);

        nWorkTris1 = 0;

        for (label i = 0; i < nWorkTris2; i++)
        {
            triSliceWithPlane(workTris2[i], pl2, workTris1, nWorkTris1, t);
        }

        if (nWorkTris1 == 0)
        {
            return 0.0;
        }
        else
        {
            // calculate area of sub-triangles
            scalar area = 0.0;
            for (label i = 0; i < nWorkTris1; i++)
            {
                area += triArea(workTris1[i]);
            }

            return area;
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::faceAreaIntersect::faceAreaIntersect
(
    const pointField& pointsA,
    const pointField& pointsB,
    const bool reverseB
)
:
    pointsA_(pointsA),
    pointsB_(pointsB),
    reverseB_(reverseB)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::faceAreaIntersect::calc
(
    const face& faceA,
    const face& faceB,
    const vector& n,
    const triangulationMode& triMode
)
{
    // split faces into triangles
    DynamicList<face> trisA;
    DynamicList<face> trisB;

    switch (triMode)
    {
        case tmFan:
        {
            triangleFan(faceA, trisA);
            triangleFan(faceB, trisB);

            break;
        }
        case tmMesh:
        {
            faceA.triangles(pointsA_, trisA);
            faceB.triangles(pointsB_, trisB);

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::faceAreaIntersect::calc"
                "("
                    "const face&, "
                    "const face&, "
                    "const vector&, "
                    "const triangulationMode&"
                ")"
            )   << "Unknown triangulation mode enumeration"
                << abort(FatalError);
        }
    }

    // intersect triangles
    scalar totalArea = 0.0;
    forAll(trisA, tA)
    {
        triPoints tpA = getTriPoints(pointsA_, trisA[tA], false);

//        if (triArea(tpA) > ROOTVSMALL)
        {
            forAll(trisB, tB)
            {
                triPoints tpB = getTriPoints(pointsB_, trisB[tB], !reverseB_);

//                if (triArea(tpB) > ROOTVSMALL)
                {
                    totalArea += triangleIntersect(tpA, tpB, n);
                }
            }
        }
    }

    return totalArea;
}


// ************************************************************************* //
