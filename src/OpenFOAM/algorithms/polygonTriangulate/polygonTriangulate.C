/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "polygonTriangulate.H"
#include "tensor2D.H"

// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

template<class PointField>
Foam::scalar Foam::polygonTriangulate::area
(
    const triFace& triPoints,
    const PointField& points,
    const vector& normal
)
{
    const point& o = points[triPoints[0]];
    const vector ao = points[triPoints[1]] - o;
    const vector ab = points[triPoints[2]] - o;

    return (ao ^ ab) & normal;
}


template<class PointField>
Foam::scalar Foam::polygonTriangulate::quality
(
    const triFace& triPoints,
    const PointField& points,
    const vector& normal
)
{
    const scalar An = area(triPoints, points, normal);

    scalar PSqr = 0;
    forAll(triPoints, triEdgei)
    {
        const point& p0 = points[triPoints[triEdgei]];
        const point& p1 = points[triPoints[(triEdgei + 1) % 3]];
        PSqr += magSqr(p0 - p1);
    }

    static const scalar equilateralAnByPSqr = sqrt(scalar(3))/12;

    return PSqr != 0 ? An/PSqr/equilateralAnByPSqr : 0;
}


template<class PointField>
bool Foam::polygonTriangulate::intersection
(
    const edge& edgePointsA,
    const edge& edgePointsB,
    const PointField& points,
    const vector& normal
)
{
    const vector tau0 = perpendicular(normal);
    const vector tau1 = normal ^ tau0;

    const point& pointA0 = points[edgePointsA[0]];
    const point& pointA1 = points[edgePointsA[1]];
    const point& pointB0 = points[edgePointsB[0]];
    const point& pointB1 = points[edgePointsB[1]];

    const point2D point2DA0(tau0 & pointA0, tau1 & pointA0);
    const point2D point2DA1(tau0 & pointA1, tau1 & pointA1);
    const point2D point2DB0(tau0 & pointB0, tau1 & pointB0);
    const point2D point2DB1(tau0 & pointB1, tau1 & pointB1);

    const tensor2D M =
        tensor2D(point2DA0 - point2DA1, point2DB1 - point2DB0).T();

    const scalar detM = det(M);
    const tensor2D detMInvM = cof(M);

    const vector2D magDetMTu = sign(detM)*detMInvM & (point2DA0 - point2DB0);

    return 0 <= cmptMin(magDetMTu) && cmptMax(magDetMTu) <= mag(detM);
}


template<class PointField>
Foam::label Foam::polygonTriangulate::nIntersections
(
    const edge& edgePoints,
    const PointField& points,
    const vector& normal
)
{
    label n = 0;

    forAll(points, pi)
    {
        const edge otherEdgePoints(pi, points.fcIndex(pi));

        if (edgePoints.commonVertex(otherEdgePoints) == -1)
        {
            n += intersection(edgePoints, otherEdgePoints, points, normal);
        }
    }

    return n;
}


template<class PointField>
Foam::scalar Foam::polygonTriangulate::angle
(
    const label pointi,
    const PointField& points,
    const vector& normal
)
{
    const point& o = points[pointi];
    const vector oa = points[points.rcIndex(pointi)] - o;
    const vector ob = points[points.fcIndex(pointi)] - o;

    const vector oaNegNNOa = oa - normal*(normal & oa);
    const vector obNegNNOb = ob - normal*(normal & ob);

    return
        atan2(normal & (oaNegNNOa ^ obNegNNOb), - oaNegNNOa & obNegNNOb)
      + constant::mathematical::pi;
}


template<class PointField>
bool Foam::polygonTriangulate::ear
(
    const label pointi,
    const PointField& points,
    const vector& normal
)
{
    const point& o = points[pointi];
    const vector oa = points[points.rcIndex(pointi)] - o;
    const vector ob = points[points.fcIndex(pointi)] - o;

    const tensor A = tensor(ob, oa, normal).T();
    const tensor T(A.y() ^ A.z(), A.z() ^ A.x(), A.x() ^ A.y());
    const scalar detA = det(A);

    for
    (
        label pointj = points.fcIndex(points.fcIndex(pointi));
        pointj != points.rcIndex(pointi);
        pointj = points.fcIndex(pointj)
    )
    {
        const vector detAY = (points[pointj] - o) & T;

        if (detAY.x() > 0 && detAY.y() > 0 && detAY.x() + detAY.y() < detA)
        {
            return false;
        }
    }

    return true;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::List<Foam::point> Foam::polygonTriangulate::randomPolygon
(
    Random& rndGen,
    const label n,
    const scalar error
)
{
    // Get random points on a unit disk
    List<point> points(n);
    forAll(points, pointi)
    {
        const scalar theta =
            2*constant::mathematical::pi*rndGen.sample01<scalar>();
        const scalar r = sqrt(rndGen.sample01<scalar>());
        points[pointi] = point(r*cos(theta), r*sin(theta), 0);
    }

    // Initialise intersected polygon
    labelList pointis(identityMap(n));

    // Reorder until non-intersected
    bool incomplete = true;
    while (incomplete)
    {
        incomplete = false;

        for (label piA0 = 0; piA0 < n; ++ piA0)
        {
            const label piA1 = points.fcIndex(piA0);

            for (label piB0 = piA0 + 2; piB0 < (piA0 == 0 ? n - 1 : n); ++ piB0)
            {
                const label piB1 = points.fcIndex(piB0);

                if
                (
                    intersection
                    (
                        edge(pointis[piA0], pointis[piA1]),
                        edge(pointis[piB0], pointis[piB1]),
                        points,
                        vector(0, 0, 1)
                    )
                )
                {
                    SubList<label> subOrder(pointis, piB0 + 1 - piA1, piA1);
                    inplaceReverseList(subOrder);
                    incomplete = true;
                }

                if (incomplete) break;
            }

            if (incomplete) break;
        }
    }

    // Add noise and potential self-intersection
    forAll(points, pointi)
    {
        const scalar theta =
            2*constant::mathematical::pi*rndGen.sample01<scalar>();
        const scalar r = sqrt(rndGen.sample01<scalar>());
        points[pointi] =
            (1 - error)*points[pointi]
          + error*point(r*cos(theta), r*sin(theta), rndGen.sample01<scalar>());
    }

    // Reorder and return
    return List<point>(points, pointis);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class PointField>
void Foam::polygonTriangulate::optimiseTriangulation
(
    const label trii,
    const PointField& points,
    const vector& normal,
    UList<triFace>& triPoints,
    UList<FixedList<label, 3>>& triEdges,
    UList<labelPair>& edgeTris
)
{
    const triFace& t = triPoints[trii];
    const scalar q = quality(t, points, normal);

    forAll(triEdges[trii], triEdgei)
    {
        const label edgei = triEdges[trii][triEdgei];

        const label otherTrii = edgeTris[edgei][edgeTris[edgei][0] == trii];

        if (otherTrii == -1) continue;

        const label otherTriEdgei = findIndex(triEdges[otherTrii], edgei);

        const triFace& otherT = triPoints[otherTrii];
        const scalar otherQ = quality(otherT, points, normal);

        const label piA = triPoints[trii][(triEdgei + 1) % 3];
        const label piB = triPoints[trii][(triEdgei + 2) % 3];
        const label piC = triPoints[otherTrii][(otherTriEdgei + 1) % 3];
        const label piD = triPoints[otherTrii][(otherTriEdgei + 2) % 3];

        const triFace flipT0(piA, piB, piD);
        const triFace flipT1(piC, piD, piB);

        if
        (
            triPointsSet_.found(flipT0)
         || triPointsSet_.found(flipT1)
        ) continue;

        const scalar flipQ0 = quality(flipT0, points, normal);
        const scalar flipQ1 = quality(flipT1, points, normal);

        if (min(q, otherQ) < min(flipQ0, flipQ1))
        {
            triPointsSet_.set(t);
            triPointsSet_.set(otherT);

            const label eiA = triEdges[trii][(triEdgei + 1) % 3];
            const label eiB = triEdges[trii][(triEdgei + 2) % 3];
            const label eiC = triEdges[otherTrii][(otherTriEdgei + 1) % 3];
            const label eiD = triEdges[otherTrii][(otherTriEdgei + 2) % 3];

            triPoints[trii] = {piA, piB, piD};
            triEdges[trii] = {eiA, edgei, eiD};
            edgeTris[eiD][edgeTris[eiD][0] != otherTrii] = trii;

            triPoints[otherTrii] = {piC, piD, piB};
            triEdges[otherTrii] = {eiC, edgei, eiB};
            edgeTris[eiB][edgeTris[eiB][0] != trii] = otherTrii;

            optimiseTriangulation
            (
                trii,
                points,
                normal,
                triPoints,
                triEdges,
                edgeTris
            );
            optimiseTriangulation
            (
                otherTrii,
                points,
                normal,
                triPoints,
                triEdges,
                edgeTris
            );

            break;
        }
    }
}


void Foam::polygonTriangulate::simpleTriangulate
(
    const UList<point>& allPoints,
    const vector& normal,
    UList<triFace>& triPoints,
    UList<FixedList<label, 3>>& triEdges,
    UList<labelPair>& edgeTris,
    const bool optimal
)
{
    // Get the polygon size
    const label n = allPoints.size();

    // Clear the triangulation
    triPoints = triFace(-1, -1, -1);
    triEdges = FixedList<label, 3>({-1, -1, -1});
    edgeTris = labelPair(-1, -1);

    // Initialise workspace
    pointis_.setSize(n);
    edges_.setSize(n);
    angle_.setSize(n);
    ear_.setSize(n);
    forAll(pointis_, i)
    {
        pointis_[i] = i;
        edges_[i] = i;
        angle_[i] = angle(i, allPoints, normal);
        ear_[i] = ear(i, allPoints, normal);
    }

    // Create the indirect point list for the current un-triangulated polygon
    UIndirectList<point> points(allPoints, pointis_);

    // Generate triangles
    forAll(triPoints, trii)
    {
        // Find the ear with the smallest angle
        scalar minEarAngle = vGreat;
        label minEarAnglei = -1;
        for (label i = 0; i < pointis_.size(); ++ i)
        {
            if (ear_[i] && minEarAngle > angle_[i])
            {
                minEarAngle = angle_[i];
                minEarAnglei = i;
            }
        }

        // If nothing is an ear, then this face is degenerate. Default to
        // the smallest angle, ignoring the ear test.
        if (minEarAnglei == -1)
        {
            minEarAnglei = findMin(angle_);
        }

        // Get the adjacent points
        label minEarAnglei0 = pointis_.rcIndex(minEarAnglei);
        label minEarAnglei1 = pointis_.fcIndex(minEarAnglei);

        // Add a triangle
        triPoints[trii] = triFace
        (
            pointis_[minEarAnglei0],
            pointis_[minEarAnglei],
            pointis_[minEarAnglei1]
        );
        triEdges[trii] = FixedList<label, 3>
        ({
            edges_[minEarAnglei0],
            edges_[minEarAnglei],
            trii == triPoints.size() - 1 ? edges_[minEarAnglei1] : n + trii
        });
        forAll(triEdges[trii], triEdgei)
        {
            const label edgei = triEdges[trii][triEdgei];
            edgeTris[edgei][edgeTris[edgei][0] != -1] = trii;
        }

        // Remove the cut off point from the polygon
        edges_[minEarAnglei0] = triEdges[trii][2];
        for (label i = minEarAnglei; i < pointis_.size() - 1; ++ i)
        {
            pointis_[i] = pointis_[i + 1];
            edges_[i] = edges_[i + 1];
            angle_[i] = angle_[i + 1];
            ear_[i] = ear_[i + 1];
        }
        pointis_.resize(pointis_.size() - 1);
        edges_.resize(edges_.size() - 1);
        angle_.resize(angle_.size() - 1);
        ear_.resize(ear_.size() - 1);

        // Update connected points
        if (minEarAnglei0 > minEarAnglei) -- minEarAnglei0;
        angle_[minEarAnglei0] = angle(minEarAnglei0, points, normal);
        ear_[minEarAnglei0] = ear(minEarAnglei0, points, normal);
        if (minEarAnglei1 > minEarAnglei) -- minEarAnglei1;
        angle_[minEarAnglei1] = angle(minEarAnglei1, points, normal);
        ear_[minEarAnglei1] = ear(minEarAnglei1, points, normal);

        // Optimise the quality
        if (optimal)
        {
            triPointsSet_.clear();
            optimiseTriangulation
            (
                trii,
                allPoints,
                normal,
                triPoints,
                triEdges,
                edgeTris
            );
        }
    }
}


void Foam::polygonTriangulate::partitionTriangulate
(
    const UList<point>& points,
    const vector& normal,
    const label spanPointi,
    const label spanEdgei,
    UList<triFace>& triPoints,
    UList<FixedList<label, 3>>& triEdges,
    UList<labelPair>& edgeTris,
    const bool simple,
    const bool optimal
)
{
    // Get the polygon size
    const label n = points.size();

    // Get the points of the spanning triangle
    const label piA = (spanEdgei + 1) % n;
    const label piP = spanPointi;
    const label piB = spanEdgei;

    // Determine the sizes of the polygons connected to either side of the
    // spanning triangle
    const label nA = (piP - piA + n) % n + 1;
    const label nB = (piB - piP + n) % n + 1;

    // Get the edges of the spanning triangle
    const label eiA = nA >= 3 ? n : piA;
    const label eiB = nB >= 3 ? n + (nA >= 3) : piP;
    const label eiE = piB;

    // Add the spanning triangle
    triPoints[0] = triFace(piA, piP, piB);
    triEdges[0] = FixedList<label, 3>({eiA, eiB, eiE});

    // Rotate the point list so that it can be passed to sub-triangulations
    // without indirect addressing
    inplaceRotateList(const_cast<UList<point>&>(points), - piA);

    // Triangulate the sub-polygon connected to edge A
    if (nA >= 3)
    {
        // Subset the polygon and triangulate
        SubList<triFace> triPointsA(triPoints, nA - 2, 1);
        SubList<FixedList<label, 3>> triEdgesA(triEdges, nA - 2, 1);
        SubList<labelPair> edgeTrisA(edgeTris, 2*nA - 3);
        triangulate
        (
            SubList<point>(points, nA),
            normal,
            triPointsA,
            triEdgesA,
            edgeTrisA,
            simple,
            optimal
        );

        // Map the point labels back to the full polygon
        forAll(triPointsA, triAi)
        {
            forAll(triPointsA[triAi], triPointAi)
            {
                label& pointi = triPointsA[triAi][triPointAi];
                pointi = (pointi + piA) % n;
            }
        }

        // Map the edge labels back to the full polygon
        forAll(triEdgesA, triAi)
        {
            forAll(triEdgesA[triAi], triEdgeAi)
            {
                label& edgei = triEdgesA[triAi][triEdgeAi];
                edgei =
                    edgei >= nA ? max(eiA, eiB) + edgei - nA + 1
                  : edgei == nA - 1 ? eiA
                  : (piA + edgei) % n;
            }
        }

        // Map the tri labels back to the full polygon
        forAll(edgeTrisA, edgeAi)
        {
            forAll(edgeTrisA[edgeAi], edgeTriAi)
            {
                label& trii = edgeTrisA[edgeAi][edgeTriAi];
                trii = trii == -1 ? -1 : trii + 1;
            }
        }
    }

    // Triangulate the sub-polygon connected to edge B
    if (nB >= 3)
    {
        // Subset the polygon and triangulate
        SubList<triFace> triPointsB(triPoints, nB - 2, nA - 1);
        SubList<FixedList<label, 3>> triEdgesB(triEdges, nB - 2, nA - 1);
        SubList<labelPair> edgeTrisB(edgeTris, 2*nB - 3, 2*nA - 3);
        triangulate
        (
            SubList<point>(points, nB, nA - 1),
            normal,
            triPointsB,
            triEdgesB,
            edgeTrisB,
            simple,
            optimal
        );

        // Map the point labels back to the full polygon
        forAll(triPointsB, triBi)
        {
            forAll(triPointsB[triBi], triPointBi)
            {
                label& pointi = triPointsB[triBi][triPointBi];
                pointi = (pointi + piA + nA - 1) % n;
            }
        }

        // Map the edge labels back to the full polygon
        forAll(triEdgesB, triBi)
        {
            forAll(triEdgesB[triBi], triEdgeBi)
            {
                label& edgei = triEdgesB[triBi][triEdgeBi];
                edgei =
                    edgei >= nB ? max(eiA, eiB) + nA - 3 + edgei - nB + 1
                  : edgei == nB - 1 ? eiB
                  : (piA + nA - 1 + edgei) % n;
            }
        }

        // Map the tri labels back to the full polygon
        forAll(edgeTrisB, edgeBi)
        {
            forAll(edgeTrisB[edgeBi], edgeTriBi)
            {
                label& trii = edgeTrisB[edgeBi][edgeTriBi];
                trii = trii == -1 ? -1 : trii + nA - 1;
            }
        }
    }

    // Rotate the point list back
    inplaceRotateList(const_cast<UList<point>&>(points), piA);

    // Reorder the edge-tris
    {
        // Swap the A-internal-edges and B-outer-edges to get all outer
        // edges and all internal edges in contiguous blocks
        SubList<labelPair> l(edgeTris, nA + nB - 3, nA - 1);
        inplaceRotateList(l, - nA + 2);
    }
    {
        // Rotate the outer edges right by one so that the intersection
        // edge (currently at the end) moves adjacent to the outer edges
        SubList<labelPair> l(edgeTris, nA + nB - 3, nA + nB - 2);
        inplaceRotateList(l, 1);
    }
    {
        // Rotate the outer edges so that they are in sequence with the
        // original point ordering
        SubList<labelPair> l(edgeTris, n);
        inplaceRotateList(l, - n + piA);
    }
    if (nA >= 3 && nB >= 3)
    {
        // Move edge-B leftwards adjacent to edge-A as this is created
        // before any of the other internal edges
        const label fromi = n + nA - 2, toi = n + 1;
        const labelPair temp = edgeTris[fromi];
        for (label i = fromi; i > toi; -- i)
        {
            edgeTris[i] = edgeTris[i-1];
        }
        edgeTris[toi] = temp;
    }

    // Add associations to the intersection triangle
    edgeTris[eiA][edgeTris[eiA][0] != -1] = 0;
    edgeTris[eiB][edgeTris[eiB][0] != -1] = 0;
    edgeTris[eiE] = labelPair(0, -1);
}


void Foam::polygonTriangulate::complexTriangulate
(
    const UList<point>& points,
    const vector& normal,
    UList<triFace>& triPoints,
    UList<FixedList<label, 3>>& triEdges,
    UList<labelPair>& edgeTris,
    const bool optimal
)
{
    // Get the polygon size
    const label n = points.size();

    // Detect self intersections. When one is found, remove one of the
    // intersected edges by adding a spanning triangle and recurse. Pick the
    // spanning triangle that results in the smallest negative area.
    for (label piA0 = 0; piA0 < n; ++ piA0)
    {
        const label piA1 = points.fcIndex(piA0);

        for (label piB0 = piA0 + 2; piB0 < (piA0 == 0 ? n - 1 : n); ++ piB0)
        {
            const label piB1 = points.fcIndex(piB0);

            if
            (
                intersection
                (
                    edge(piA0, piA1),
                    edge(piB0, piB1),
                    points,
                    normal
                )
            )
            {
                // Get a bunch of unique spanning triangles that remove one of
                // the intersected edges
                FixedList<labelPair, 8> spanPointAndEdgeis
                ({
                    labelPair(piA0, piB0),
                    labelPair(piA1, piB0),
                    labelPair(piB0, piA0),
                    labelPair(piB1, piA0),
                    labelPair(points.rcIndex(piA0), piA0),
                    labelPair(points.fcIndex(piA1), piA0),
                    labelPair(points.rcIndex(piB0), piB0),
                    labelPair(points.fcIndex(piB1), piB0)
                });
                forAll(spanPointAndEdgeis, spani)
                {
                    for (label spanj = spani + 1; spanj < 8; ++ spanj)
                    {
                        if
                        (
                            spanPointAndEdgeis[spani]
                         == spanPointAndEdgeis[spanj]
                        )
                        {
                            spanPointAndEdgeis[spanj] = labelPair(-1, -1);
                        }
                    }
                }

                // Find the spanning triangle that results in the smallest
                // negative triangulation area
                scalar spanSumNegA = - vGreat;
                label spanPointi = -1, spanEdgei = -1;
                forAll(spanPointAndEdgeis, spani)
                {
                    if (spanPointAndEdgeis[spani] == labelPair(-1, -1))
                        continue;

                    const label piP = spanPointAndEdgeis[spani].first();
                    const label ei = spanPointAndEdgeis[spani].second();
                    const label piA0 = points.rcIndex(ei);
                    const label piA = ei;
                    const label piB = points.fcIndex(ei);
                    const label piB1 = points.fcIndex(piB);

                    // Only consider this spanning triangle if it decreases
                    // the number of intersections
                    const label netNIntersections =
                      - nIntersections(edge(piA, piB), points, normal)
                      + (piP == piA0 ? -1 : +1)
                       *nIntersections(edge(piP, piA), points, normal)
                      + (piP == piB1 ? -1 : +1)
                       *nIntersections(edge(piP, piB), points, normal);
                    if (netNIntersections >= 0) continue;

                    // Triangulate
                    partitionTriangulate
                    (
                        points,
                        normal,
                        piP,
                        ei,
                        triPoints,
                        triEdges,
                        edgeTris,
                        false,
                        false
                    );

                    // Compute the negative area. If it is smaller than the
                    // previous value then update the spanning triangle. If it
                    // is zero then break and use this triangle.
                    scalar sumNegA = 0;
                    forAll(triPoints, trii)
                    {
                        const scalar a = area(triPoints[trii], points, normal);
                        sumNegA += negPart(a);
                    }
                    if (sumNegA > spanSumNegA)
                    {
                        spanSumNegA = sumNegA;
                        spanPointi = piP;
                        spanEdgei = ei;
                    }
                    if (spanSumNegA == 0)
                    {
                        break;
                    }
                }

                // If a suitable spanning triangle was found then partition the
                // triangulation
                if (spanPointi != -1)
                {
                    partitionTriangulate
                    (
                        points,
                        normal,
                        spanPointi,
                        spanEdgei,
                        triPoints,
                        triEdges,
                        edgeTris,
                        false,
                        optimal
                    );

                    return;
                }
            }
        }
    }

    // If there are no intersections then do simple polygon triangulation
    simpleTriangulate(points, normal, triPoints, triEdges, edgeTris, optimal);
}


void Foam::polygonTriangulate::triangulate
(
    const UList<point>& points,
    const vector& normal,
    UList<triFace>& triPoints,
    UList<FixedList<label, 3>>& triEdges,
    UList<labelPair>& edgeTris,
    const bool simple,
    const bool optimal
)
{
    if (simple)
    {
        simpleTriangulate
        (
            points,
            normal,
            triPoints,
            triEdges,
            edgeTris,
            optimal
        );
    }
    else
    {
        complexTriangulate
        (
            points,
            normal,
            triPoints,
            triEdges,
            edgeTris,
            optimal
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polygonTriangulate::polygonTriangulate()
:
    pointis_(),
    edges_(),
    angle_(),
    ear_(),

    triPointsSet_(),

    points_(),
    triPoints_(),
    triEdges_(),
    edgeTris_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polygonTriangulate::~polygonTriangulate()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::UList<Foam::triFace>& Foam::polygonTriangulate::triangulate
(
    const UList<point>& points,
    const vector& normal,
    const bool simple,
    const bool optimal
)
{
    // Get the polygon size
    const label n = points.size();

    // Resize workspace
    triPoints_.resize(n - 2);
    triEdges_.resize(n - 2);
    edgeTris_.resize(2*n - 3);

    // Run
    triangulate
    (
        points,
        normal,
        triPoints_,
        triEdges_,
        edgeTris_,
        simple,
        optimal
    );

    return triPoints_;
}


const Foam::UList<Foam::triFace>& Foam::polygonTriangulate::triangulate
(
    const UList<point>& points,
    const bool simple,
    const bool optimal
)
{
    return
        triangulate
        (
            points,
            normalised(face::area(points)),
            simple,
            optimal
        );
}


// ************************************************************************* //
