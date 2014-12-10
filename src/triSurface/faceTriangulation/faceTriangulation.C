/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "faceTriangulation.H"
#include "plane.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::faceTriangulation::edgeRelTol = 1e-6;


// Edge to the right of face vertex i
Foam::label Foam::faceTriangulation::right(const label, label i)
{
    return i;
}


// Edge to the left of face vertex i
Foam::label Foam::faceTriangulation::left(const label size, label i)
{
    return i ? i-1 : size-1;
}


// Calculate (normalized) edge vectors.
// edges[i] gives edge between point i+1 and i.
Foam::tmp<Foam::vectorField> Foam::faceTriangulation::calcEdges
(
    const face& f,
    const pointField& points
)
{
    tmp<vectorField> tedges(new vectorField(f.size()));
    vectorField& edges = tedges();

    forAll(f, i)
    {
        point thisPt = points[f[i]];
        point nextPt = points[f[f.fcIndex(i)]];

        vector vec(nextPt - thisPt);
        vec /= mag(vec) + VSMALL;

        edges[i] = vec;
    }

    return tedges;
}


// Calculates half angle components of angle from e0 to e1
void Foam::faceTriangulation::calcHalfAngle
(
    const vector& normal,
    const vector& e0,
    const vector& e1,
    scalar& cosHalfAngle,
    scalar& sinHalfAngle
)
{
    // truncate cos to +-1 to prevent negative numbers
    scalar cos = max(-1, min(1, e0 & e1));

    scalar sin = (e0 ^ e1) & normal;

    if (sin < -ROOTVSMALL)
    {
        // 3rd or 4th quadrant
        cosHalfAngle = - Foam::sqrt(0.5*(1 + cos));
        sinHalfAngle = Foam::sqrt(0.5*(1 - cos));
    }
    else
    {
        // 1st or 2nd quadrant
        cosHalfAngle = Foam::sqrt(0.5*(1 + cos));
        sinHalfAngle = Foam::sqrt(0.5*(1 - cos));
    }
}


// Calculate intersection point between edge p1-p2 and ray (in 2D).
// Return true and intersection point if intersection between p1 and p2.
Foam::pointHit Foam::faceTriangulation::rayEdgeIntersect
(
    const vector& normal,
    const point& rayOrigin,
    const vector& rayDir,
    const point& p1,
    const point& p2,
    scalar& posOnEdge
)
{
    // Start off from miss
    pointHit result(p1);

    // Construct plane normal to rayDir and intersect
    const vector y = normal ^ rayDir;

    posOnEdge = plane(rayOrigin, y).normalIntersect(p1, (p2-p1));

    // Check intersection to left of p1 or right of p2
    if ((posOnEdge < 0) || (posOnEdge > 1))
    {
        // Miss
    }
    else
    {
        // Check intersection behind rayOrigin
        point intersectPt = p1 + posOnEdge * (p2 - p1);

        if (((intersectPt - rayOrigin) & rayDir) < 0)
        {
            // Miss
        }
        else
        {
            // Hit
            result.setHit();
            result.setPoint(intersectPt);
            result.setDistance(mag(intersectPt - rayOrigin));
        }
    }
    return result;
}


// Return true if triangle given its three points (anticlockwise ordered)
// contains point
bool Foam::faceTriangulation::triangleContainsPoint
(
    const vector& n,
    const point& p0,
    const point& p1,
    const point& p2,
    const point& pt
)
{
    scalar area01Pt = triPointRef(p0, p1, pt).normal() & n;
    scalar area12Pt = triPointRef(p1, p2, pt).normal() & n;
    scalar area20Pt = triPointRef(p2, p0, pt).normal() & n;

    if ((area01Pt > 0) && (area12Pt > 0) && (area20Pt > 0))
    {
        return true;
    }
    else if ((area01Pt < 0) && (area12Pt < 0) && (area20Pt < 0))
    {
        FatalErrorIn("triangleContainsPoint") << abort(FatalError);
        return false;
    }
    else
    {
        return false;
    }
}


// Starting from startIndex find diagonal. Return in index1, index2.
// Index1 always startIndex except when convex polygon
void Foam::faceTriangulation::findDiagonal
(
    const pointField& points,
    const face& f,
    const vectorField& edges,
    const vector& normal,
    const label startIndex,
    label& index1,
    label& index2
)
{
    const point& startPt = points[f[startIndex]];

    // Calculate angle at startIndex
    const vector& rightE = edges[right(f.size(), startIndex)];
    const vector leftE = -edges[left(f.size(), startIndex)];

    // Construct ray which bisects angle
    scalar cosHalfAngle = GREAT;
    scalar sinHalfAngle = GREAT;
    calcHalfAngle(normal, rightE, leftE, cosHalfAngle, sinHalfAngle);
    vector rayDir
    (
        cosHalfAngle*rightE
      + sinHalfAngle*(normal ^ rightE)
    );
    // rayDir should be normalized already but is not due to rounding errors
    // so normalize.
    rayDir /= mag(rayDir) + VSMALL;


    //
    // Check all edges (apart from rightE and leftE) for nearest intersection
    //

    label faceVertI = f.fcIndex(startIndex);

    pointHit minInter(false, vector::zero, GREAT, true);
    label minIndex = -1;
    scalar minPosOnEdge = GREAT;

    for (label i = 0; i < f.size() - 2; i++)
    {
        scalar posOnEdge;
        pointHit inter =
            rayEdgeIntersect
            (
                normal,
                startPt,
                rayDir,
                points[f[faceVertI]],
                points[f[f.fcIndex(faceVertI)]],
                posOnEdge
            );

        if (inter.hit() && inter.distance() < minInter.distance())
        {
            minInter = inter;
            minIndex = faceVertI;
            minPosOnEdge = posOnEdge;
        }

        faceVertI = f.fcIndex(faceVertI);
    }


    if (minIndex == -1)
    {
        //WarningIn("faceTriangulation::findDiagonal")
        //    << "Could not find intersection starting from " << f[startIndex]
        //    << " for face " << f << endl;

        index1 = -1;
        index2 = -1;
        return;
    }

    const label leftIndex = minIndex;
    const label rightIndex = f.fcIndex(minIndex);

    // Now ray intersects edge from leftIndex to rightIndex.
    // Check for intersection being one of the edge points. Make sure never
    // to return two consecutive points.

    if (mag(minPosOnEdge) < edgeRelTol && f.fcIndex(startIndex) != leftIndex)
    {
        index1 = startIndex;
        index2 = leftIndex;

        return;
    }
    if
    (
        mag(minPosOnEdge - 1) < edgeRelTol
     && f.fcIndex(rightIndex) != startIndex
    )
    {
        index1 = startIndex;
        index2 = rightIndex;

        return;
    }

    // Select visible vertex that minimizes
    // angle to bisection. Visibility checking by checking if inside triangle
    // formed by startIndex, leftIndex, rightIndex

    const point& leftPt = points[f[leftIndex]];
    const point& rightPt = points[f[rightIndex]];

    minIndex = -1;
    scalar maxCos = -GREAT;

    // all vertices except for startIndex and ones to left and right of it.
    faceVertI = f.fcIndex(f.fcIndex(startIndex));
    for (label i = 0; i < f.size() - 3; i++)
    {
        const point& pt = points[f[faceVertI]];

        if
        (
            (faceVertI == leftIndex)
         || (faceVertI == rightIndex)
         || (triangleContainsPoint(normal, startPt, leftPt, rightPt, pt))
        )
        {
            // pt inside triangle (so perhaps visible)
            // Select based on minimal angle (so guaranteed visible).
            vector edgePt0 = pt - startPt;
            edgePt0 /= mag(edgePt0);

            scalar cos = rayDir & edgePt0;
            if (cos > maxCos)
            {
                maxCos = cos;
                minIndex = faceVertI;
            }
        }
        faceVertI = f.fcIndex(faceVertI);
    }

    if (minIndex == -1)
    {
        // no vertex found. Return startIndex and one of the intersected edge
        // endpoints.
        index1 = startIndex;

        if (f.fcIndex(startIndex) != leftIndex)
        {
            index2 = leftIndex;
        }
        else
        {
            index2 = rightIndex;
        }

        return;
    }

    index1 = startIndex;
    index2 = minIndex;
}


// Find label of vertex to start splitting from. Is:
//     1] flattest concave angle
//     2] flattest convex angle if no concave angles.
Foam::label Foam::faceTriangulation::findStart
(
    const face& f,
    const vectorField& edges,
    const vector& normal
)
{
    const label size = f.size();

    scalar minCos = GREAT;
    label minIndex = -1;

    forAll(f, fp)
    {
        const vector& rightEdge = edges[right(size, fp)];
        const vector leftEdge = -edges[left(size, fp)];

        if (((rightEdge ^ leftEdge) & normal) < ROOTVSMALL)
        {
            scalar cos = rightEdge & leftEdge;
            if (cos < minCos)
            {
                minCos = cos;
                minIndex = fp;
            }
        }
    }

    if (minIndex == -1)
    {
        // No concave angle found. Get flattest convex angle
        minCos = GREAT;

        forAll(f, fp)
        {
            const vector& rightEdge = edges[right(size, fp)];
            const vector leftEdge = -edges[left(size, fp)];

            scalar cos = rightEdge & leftEdge;
            if (cos < minCos)
            {
                minCos = cos;
                minIndex = fp;
            }
        }
    }

    return minIndex;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Split face f into triangles. Handles all simple (convex & concave)
// polygons.
bool Foam::faceTriangulation::split
(
    const bool fallBack,
    const pointField& points,
    const face& f,
    const vector& normal,
    label& triI
)
{
    const label size = f.size();

    if (size <= 2)
    {
        WarningIn
        (
            "split(const bool, const pointField&, const face&"
            ", const vector&, label&)"
        )   << "Illegal face:" << f
            << " with points " << UIndirectList<point>(points, f)()
            << endl;

        return false;
    }
    else if (size == 3)
    {
        // Triangle. Just copy.
        triFace& tri = operator[](triI++);
        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];

        return true;
    }
    else
    {
        // General case. Start splitting for -flattest concave angle
        // -or flattest convex angle if no concave angles.

        tmp<vectorField> tedges(calcEdges(f, points));
        const vectorField& edges = tedges();

        label startIndex = findStart(f, edges, normal);

        // Find diagonal to split face across
        label index1 = -1;
        label index2 = -1;

        forAll(f, iter)
        {
            findDiagonal
            (
                points,
                f,
                edges,
                normal,
                startIndex,
                index1,
                index2
            );

            if (index1 != -1 && index2 != -1)
            {
                // Found correct diagonal
                break;
            }

            // Try splitting from next startingIndex.
            startIndex = f.fcIndex(startIndex);
        }

        if (index1 == -1 || index2 == -1)
        {
            if (fallBack)
            {
                // Do naive triangulation. Find smallest angle to start
                // triangulating from.
                label maxIndex = -1;
                scalar maxCos = -GREAT;

                forAll(f, fp)
                {
                    const vector& rightEdge = edges[right(size, fp)];
                    const vector leftEdge = -edges[left(size, fp)];

                    scalar cos = rightEdge & leftEdge;
                    if (cos > maxCos)
                    {
                        maxCos = cos;
                        maxIndex = fp;
                    }
                }

                WarningIn
                (
                    "split(const bool, const pointField&, const face&"
                    ", const vector&, label&)"
                )   << "Cannot find valid diagonal on face " << f
                    << " with points " << UIndirectList<point>(points, f)()
                    << nl
                    << "Returning naive triangulation starting from "
                    << f[maxIndex] << " which might not be correct for a"
                    << " concave or warped face" << endl;


                label fp = f.fcIndex(maxIndex);

                for (label i = 0; i < size-2; i++)
                {
                    label nextFp = f.fcIndex(fp);

                    triFace& tri = operator[](triI++);
                    tri[0] = f[maxIndex];
                    tri[1] = f[fp];
                    tri[2] = f[nextFp];

                    fp = nextFp;
                }

                return true;
            }
            else
            {
                WarningIn
                (
                    "split(const bool, const pointField&, const face&"
                    ", const vector&, label&)"
                )   << "Cannot find valid diagonal on face " << f
                    << " with points " << UIndirectList<point>(points, f)()
                    << nl
                    << "Returning empty triFaceList" << endl;

                return false;
            }
        }


        // Split into two subshapes.
        //     face1: index1 to index2
        //     face2: index2 to index1

        // Get sizes of the two subshapes
        label diff = 0;
        if (index2 > index1)
        {
            diff = index2 - index1;
        }
        else
        {
            // folded round
            diff = index2 + size - index1;
        }

        label nPoints1 = diff + 1;
        label nPoints2 = size - diff + 1;

        if (nPoints1 == size || nPoints2 == size)
        {
            FatalErrorIn
            (
                "split(const bool, const pointField&, const face&"
                ", const vector&, label&)"
            )   << "Illegal split of face:" << f
                << " with points " << UIndirectList<point>(points, f)()
                << " at indices " << index1 << " and " << index2
                << abort(FatalError);
        }


        // Collect face1 points
        face face1(nPoints1);

        label faceVertI = index1;
        for (int i = 0; i < nPoints1; i++)
        {
            face1[i] = f[faceVertI];
            faceVertI = f.fcIndex(faceVertI);
        }

        // Collect face2 points
        face face2(nPoints2);

        faceVertI = index2;
        for (int i = 0; i < nPoints2; i++)
        {
            face2[i] = f[faceVertI];
            faceVertI = f.fcIndex(faceVertI);
        }

        // Decompose the split faces
        //Pout<< "Split face:" << f << " into " << face1 << " and " << face2
        //    << endl;
        //string oldPrefix(Pout.prefix());
        //Pout.prefix() = "  " + oldPrefix;

        bool splitOk =
            split(fallBack, points, face1, normal, triI)
         && split(fallBack, points, face2, normal, triI);

        //Pout.prefix() = oldPrefix;

        return splitOk;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::faceTriangulation::faceTriangulation()
:
    triFaceList()
{}


// Construct from components
Foam::faceTriangulation::faceTriangulation
(
    const pointField& points,
    const face& f,
    const bool fallBack
)
:
    triFaceList(f.size()-2)
{
    vector avgNormal = f.normal(points);
    avgNormal /= mag(avgNormal) + VSMALL;

    label triI = 0;

    bool valid = split(fallBack, points, f, avgNormal, triI);

    if (!valid)
    {
        setSize(0);
    }
}


// Construct from components
Foam::faceTriangulation::faceTriangulation
(
    const pointField& points,
    const face& f,
    const vector& n,
    const bool fallBack
)
:
    triFaceList(f.size()-2)
{
    label triI = 0;

    bool valid = split(fallBack, points, f, n, triI);

    if (!valid)
    {
        setSize(0);
    }
}


// Construct from Istream
Foam::faceTriangulation::faceTriangulation(Istream& is)
:
    triFaceList(is)
{}


// ************************************************************************* //
