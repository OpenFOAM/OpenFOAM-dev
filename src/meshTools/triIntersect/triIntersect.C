/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "triIntersect.H"
#include "boundBox.H"
#include "cubicEqn.H"
#include "unitConversion.H"
#include "OFstream.H"
#include "tensor2D.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace triIntersect
{


//- The maximum dot product between a source point normal and a target plane
//  considered to be a valid, forward projection
const scalar maxDot = - cos(degToRad(80));


//- Divide two numbers and protect the result from overflowing
scalar protectedDivide(const scalar a, const scalar b)
{
    return mag(a)/vGreat < mag(b) ? a/b : sign(a)*sign(b)*vGreat;
}


//- Divide two numbers, protect the result from overflowing, and clip the
//  result between 0 and 1
scalar protectedDivideAndClip01(const scalar a, const scalar b)
{
    const scalar bStar = max(mag(a), mag(b));

    return bStar == 0 ? 0 : max(a*sign(b), 0)/bStar;
}


//- Print 3x3 FixedListLists on one line
template <class Type>
Ostream& operator<<(Ostream& os, const FixedList<FixedList<Type, 3>, 3>& l)
{
    os << token::BEGIN_LIST;
    forAll(l, i)
    {
        if (i) os << token::SPACE;
        os << l[i];
    }
    os << token::END_LIST;
    return os;
}


//- Clip the given vector between values of 0 and 1, and also clip one minus
//  its component sum. Clipping is applied to groups of components. It is done
//  by moving the value linearly towards the value where all components in the
//  group, and one minus their sum, share the same value.
vector clipped01(const vector x, const FixedList<label, 3> groups)
{
    vector y(x);

    for (label group = 0; group < groups[findMax(groups)] + 1; ++ group)
    {
        label n = 1;
        forAll(x, i)
        {
            if (groups[i] == group)
            {
                n ++;
            }
        }

        if (n == 1)
        {
            continue;
        }
        else if (n == 2)
        {
            forAll(x, i)
            {
                if (groups[i] == group)
                {
                    y[i] = min(max(x[i], 0), 1);
                }
            }
        }
        else
        {
            scalar xn = 1;
            forAll(x, i)
            {
                if (groups[i] == group)
                {
                    xn -= x[i];
                }
            }

            scalar phi = 0;
            forAll(x, i)
            {
                if (groups[i] == group)
                {
                    if (x[i] < 0)
                    {
                        phi = max(phi, n*x[i]/(n*x[i] - 1));
                    }
                }
            }
            if (xn < 0)
            {
                phi = max(phi, n*xn/(n*xn - 1));
            }

            forAll(x, i)
            {
                if (groups[i] == group)
                {
                    y[i] = min(max((1 - phi)*x[i] + phi/n, 0), 1);
                }
            }
        }
    }

    return y;
}


//- Solve a projection equation given a value of the t variable
vector solveProjectionGivenT
(
    const vector& C,
    const vector& Ct,
    const vector& Cu,
    const vector& Cv,
    const vector& Ctu,
    const vector& Ctv,
    const FixedList<label, 3> groups,
    const scalar t
)
{
    // Solve a least squares problem for u and v
    const vector CCtt = C + Ct*t;
    const vector CuCtut = Cu + Ctu*t;
    const vector CvCtvt = Cv + Ctv*t;

    const tensor2D A
    (
        CuCtut & CuCtut, CuCtut & CvCtvt,
        CvCtvt & CuCtut, CvCtvt & CvCtvt
    );
    const vector2D B
    (
        - CuCtut & CCtt,
        - CvCtvt & CCtt
    );

    const scalar detA = det(A);
    const vector2D detAuv = cof(A) & B;

    const vector tuv
    (
        t,
        protectedDivide(detAuv.x(), detA),
        protectedDivide(detAuv.y(), detA)
    );

    // Apply group clipping
    return clipped01(tuv, groups);
}


//- Solve a projection equation
Tuple2<bool, vector> solveProjection
(
    const vector& C,
    const vector& Ct,
    const vector& Cu,
    const vector& Cv,
    const vector& Ctu,
    const vector& Ctv,
    const FixedList<label, 3> groups
)
{
    // Solve the cubic projection equation for t
    const Roots<3> tRoots =
        cubicEqn
        (
            (Ct ^ Ctu) & Ctv,
            ((C ^ Ctu) & Ctv) + ((Ct ^ Cu) & Ctv) + ((Ct ^ Ctu) & Cv),
            ((C ^ Cu) & Ctv) + ((C ^ Ctu) & Cv) + ((Ct ^ Cu) & Cv),
            (C ^ Cu) & Cv
        ).roots();

    // Solve the remaining problem for u and v
    label nTuvs = 0;
    FixedList<vector, 3> tuvs;
    forAll(tRoots, tRooti)
    {
        if (tRoots.type(tRooti) != rootType::real) continue;

        if (mag(tRoots[tRooti]) > great) continue;

        const vector tuv =
            solveProjectionGivenT
            (
                C,
                Ct,
                Cu,
                Cv,
                Ctu,
                Ctv,
                {-1, -1, -1},
                tRoots[tRooti]
            );

        if (cmptMax(cmptMag(tuv)) > rootVGreat) continue;

        tuvs[nTuvs ++] = tuv;
    }

    // Apply clipping
    FixedList<scalar, 3> tuvClippage(NaN);
    for (label i = 0; i < nTuvs; ++ i)
    {
        const vector tuvOld = tuvs[i];
        tuvs[i] = clipped01(tuvs[i], groups);
        tuvClippage[i] = cmptSum(cmptMag(tuvs[i] - tuvOld));
    }

    // Sort so the least clipped roots come first
    for (label i = 0; i < nTuvs - 1; ++ i)
    {
        for (label j = 0; j < nTuvs - 1; ++ j)
        {
            if (tuvClippage[j] > tuvClippage[j + 1])
            {
                Swap(tuvs[j], tuvs[j + 1]);
                Swap(tuvClippage[j], tuvClippage[j + 1]);
            }
        }
    }

    // Analyse each clipped solution value versus estimated error. If the
    // value is small relative to the error, then return success.
    for (label i = 0; i < nTuvs; ++ i)
    {
        const scalar t = tuvs[i].x(), u = tuvs[i].y(), v = tuvs[i].z();

        const scalar magSqrF =
            magSqr(C + Ct*t + Cu*u + Cv*v + Ctu*t*u + Ctv*t*v);
        const scalar magSqrErrF =
            magSqr
            (
                1*cmptMag(C)
              + 2*cmptMag(Ct)*mag(t)
              + 2*cmptMag(Cu)*mag(u)
              + 2*cmptMag(Cv)*mag(v)
              + 3*cmptMag(Ctu)*mag(t)*mag(u)
              + 3*cmptMag(Ctv)*mag(t)*mag(v)
            )*small;

        if (magSqrF < magSqrErrF)
        {
            return Tuple2<bool, vector>(true, tuvs[i]);
        }
    }

    // No suitable roots were found. Return failure.
    return Tuple2<bool, vector>(false, vector::uniform(NaN));
}


//- Calculate the non-dimensional offsets of the source points from the target
//  edges. These values are considered indicative only. The calculation is not
//  as reliable as that for the target point offsets from the source edges (see
//  below). The target point offsets should take precedence where possible.
FixedList<FixedList<scalar, 3>, 3> srcTgtEdgeOffset
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const FixedList<bool, 3>& srcOwns,
    const FixedList<point, 3>& tgtPs,
    const FixedList<bool, 3>& tgtOwns
)
{
    FixedList<FixedList<scalar, 3>, 3> result;

    const triPointRef tgtTri(tgtPs[0], tgtPs[1], tgtPs[2]);
    const vector tgtN = tgtTri.normal();

    // Check whether the intersections occur in a forward direction
    FixedList<bool, 3> srcFwd;
    forAll(srcPs, srcPi)
    {
        srcFwd[srcPi] = (srcNs[srcPi] & tgtN) < maxDot;
    }

    // For all forward projecting source points, determine the offset from the
    // target edges
    forAll(srcPs, srcPi)
    {
        if (srcFwd[srcPi])
        {
            forAll(tgtPs, tgtEi)
            {
                const label tgtPi0 = tgtEi, tgtPi1 = (tgtEi + 1) % 3;
                result[srcPi][tgtEi] =
                    srcPointTgtEdgeOffset
                    (
                        srcPs[srcPi],
                        srcNs[srcPi],
                        {tgtPs[tgtPi0], tgtPs[tgtPi1]},
                        tgtOwns[tgtEi]
                    );
            }
        }
    }

    // For all backward projecting source points, initialise to outside
    // everything
    forAll(srcPs, srcPi)
    {
        if (!srcFwd[srcPi])
        {
            result[srcPi] = {-vGreat, -vGreat, -vGreat};
        }
    }

    // For source edges with one forward projecting and one backward
    // projecting point, compute the point and normal where the projection
    // direction changes and use this to determine the offset of the backward
    // projecting point
    forAll(srcPs, srcEi)
    {
        const label srcPi0 = srcEi, srcPi1 = (srcEi + 1) % 3;

        if (srcFwd[srcPi0] != srcFwd[srcPi1])
        {
            // Get the source edge parameter and normal at the asymptote where
            // the normal switches sign relative to the target
            const scalar srcT =
                  protectedDivide
                  (
                      maxDot - (srcNs[srcPi0] & tgtN),
                      (srcNs[srcPi1] - srcNs[srcPi0]) & tgtN
                  );
            const point srcP = (1 - srcT)*srcPs[srcPi0] + srcT*srcPs[srcPi1];
            const vector srcN = (1 - srcT)*srcNs[srcPi0] + srcT*srcNs[srcPi1];

            forAll(tgtPs, tgtEi)
            {
                const label tgtPi0 = tgtEi, tgtPi1 = (tgtEi + 1) % 3;
                result[srcFwd[srcPi0] ? srcPi1 : srcPi0][tgtEi] =
                    srcPointTgtEdgeOffset
                    (
                        srcP,
                        srcN,
                        {tgtPs[tgtPi0], tgtPs[tgtPi1]},
                        tgtOwns[tgtEi]
                    );
            }
        }
    }

    return result;
}


//- Calculate the non-dimensional offsets the target points from the source
//  edges. These values are considered definitive, and should take precedence
//  over the source point target edge offsets.
FixedList<FixedList<scalar, 3>, 3> tgtSrcEdgeOffset
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const FixedList<bool, 3>& srcOwns,
    const FixedList<point, 3>& tgtPs,
    const FixedList<bool, 3>& tgtOwns
)
{
    FixedList<FixedList<scalar, 3>, 3> result;

    // For all target points, determine the offset from each source edge
    forAll(tgtPs, tgtPi)
    {
        forAll(srcPs, srcEi)
        {
            const label srcPi0 = srcEi, srcPi1 = (srcEi + 1) % 3;
            result[tgtPi][srcEi] =
                srcEdgeTgtPointOffset
                (
                    {srcPs[srcPi0], srcPs[srcPi1]},
                    {srcNs[srcPi0], srcNs[srcPi1]},
                    tgtPs[tgtPi],
                    srcOwns[srcEi]
                );
        }
    }

    return result;
}


//- Construct point-inside/outside-edge topology from a set of point-edge
//  offsets. Uses the sign of the offsets.
FixedList<FixedList<label, 3>, 3> thisInOtherEdge
(
    const FixedList<FixedList<scalar, 3>, 3>& thisOtherEdgeOffset
)
{
    FixedList<FixedList<label, 3>, 3> result;

    // Determine the edge association from the sign of the offset
    forAll(thisOtherEdgeOffset, thisPi)
    {
        forAll(thisOtherEdgeOffset[thisPi], otherEi)
        {
            result[thisPi][otherEi] =
                thisOtherEdgeOffset[thisPi][otherEi] > 0 ? +1 : -1;
        }
    }

    return result;
}


//- Construct point-inside/outside-triangle topology from a set of
//  point-inside/outside-edge topology
FixedList<label, 3> thisInOtherTri
(
    const FixedList<FixedList<label, 3>, 3>& thisInOtherEdge
)
{
    FixedList<label, 3> result;

    // Combine edge associations to get triangle associations
    forAll(thisInOtherEdge, thisPi)
    {
        result[thisPi] = count(thisInOtherEdge[thisPi], 1) == 3 ? +1 : -1;
    }

    return result;
}


//- Construct target-point-inside/outside-source-triangle topology from a set
//  of target-point-inside/outside-source-edge topology, and some additional
//  geometric information to handle cases where the source normal direction
//  is reversed relative to the target triangle
FixedList<label, 3> tgtInSrcTri
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const FixedList<point, 3>& tgtPs,
    const FixedList<FixedList<label, 3>, 3>& tgtInSrcEdge
)
{
    const triPointRef tgtTri(tgtPs[0], tgtPs[1], tgtPs[2]);
    const vector tgtN = tgtTri.normal();

    // Combine edge associations to get triangle associations
    FixedList<label, 3> result = thisInOtherTri(tgtInSrcEdge);

    // Filter to only include forward intersections
    forAll(tgtInSrcEdge, tgtPi)
    {
        if (result[tgtPi] == 1)
        {
            const barycentric2D srcTs =
                srcTriTgtPointIntersection(srcPs, srcNs, tgtPs[tgtPi]);
            const vector srcN = srcTriInterpolate(srcTs, srcNs);
            const bool tgtFwd = (srcN & tgtN) < maxDot;
            result[tgtPi] = tgtFwd ? +1 : -1;
        }
    }

    return result;
}


//- Override results of the srcInTgt/tgtInSrc calculations with explicit
//  connections between points on either side
void thisIsOther
(
    const FixedList<label, 3>& thisOtherPis,
    FixedList<FixedList<label, 3>, 3>& thisInOtherEdge,
    FixedList<label, 3>& thisInOtherTri
)
{
    forAll(thisOtherPis, thisPi)
    {
        const label otherPi = thisOtherPis[thisPi];

        if (otherPi != -1)
        {
            const label otherEi0 = (otherPi + 2) % 3, otherEi1 = otherPi;

            thisInOtherTri[thisPi] = 0;

            thisInOtherEdge[thisPi][otherEi0] = 0;
            thisInOtherEdge[thisPi][otherEi1] = 0;
        }
    }
}


//- Calculate whether the points of the given source triangle project inside or
//  outside the opposing target triangle and its edges
void srcInTgt
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const FixedList<bool, 3>& srcOwns,
    const FixedList<label, 3>& srcTgtPis,
    const FixedList<point, 3>& tgtPs,
    const FixedList<bool, 3>& tgtOwns,
    FixedList<FixedList<label, 3>, 3>& srcInTgtEdge,
    FixedList<label, 3>& srcInTgtTri
)
{
    const FixedList<FixedList<scalar, 3>, 3>& srcTgtEdgeOffset =
        triIntersect::srcTgtEdgeOffset(srcPs, srcNs, srcOwns, tgtPs, tgtOwns);

    srcInTgtEdge = thisInOtherEdge(srcTgtEdgeOffset);

    srcInTgtTri = thisInOtherTri(srcInTgtEdge);

    thisIsOther(srcTgtPis, srcInTgtEdge, srcInTgtTri);
}


//- Calculate whether the points of the given target triangle project inside or
//  outside the opposing source triangle and its edges
void tgtInSrc
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const FixedList<bool, 3>& srcOwns,
    const FixedList<point, 3>& tgtPs,
    const FixedList<bool, 3>& tgtOwns,
    const FixedList<label, 3>& tgtSrcPis,
    FixedList<FixedList<label, 3>, 3>& tgtInSrcEdge,
    FixedList<label, 3>& tgtInSrcTri
)
{
    const FixedList<FixedList<scalar, 3>, 3>& tgtSrcEdgeOffset =
        triIntersect::tgtSrcEdgeOffset(srcPs, srcNs, srcOwns, tgtPs, tgtOwns);

    tgtInSrcEdge = thisInOtherEdge(tgtSrcEdgeOffset);

    tgtInSrcTri = triIntersect::tgtInSrcTri(srcPs, srcNs, tgtPs, tgtInSrcEdge);

    thisIsOther(tgtSrcPis, tgtInSrcEdge, tgtInSrcTri);
}


//- Order intersection locations into a polygon
bool orderLocations
(
    const UList<location>& locations,
    bool isSrcEdge,
    const label i0,
    label& nVisited,
    boolList& visited,
    labelList& order
)
{
    // Mark this location as visited
    order[nVisited ++] = i0;
    visited[i0] = true;

    // Get the index of the edge attached to this point
    const location& l0 = locations[i0];
    const label ei0 =
        isSrcEdge
      ? (l0.isSrcPoint() ? l0.srcPointi() : l0.srcEdgei())
      : (l0.isTgtPoint() ? (l0.tgtPointi() + 2) % 3 : l0.tgtEdgei());

    // Terminate if connected back to the first location
    {
        const label i1 = order.first();

        const location& l1 = locations[i1];

        if
        (
            i0 != order.first()
         && !(l1.isSrcNotTgtPoint() && !isSrcEdge)
         && !(l1.isTgtNotSrcPoint() && isSrcEdge)
        )
        {
            const location& l1 = locations[i1];
            const label ei1 =
                isSrcEdge
              ? (l1.isSrcPoint() ? (l1.srcPointi() + 2) % 3 : l1.srcEdgei())
              : (l1.isTgtPoint() ? l1.tgtPointi() : l1.tgtEdgei());

            if (ei0 == ei1)
            {
                return true;
            }
        }
    }

    // Search for the next connected location and recurse if found
    forAll(locations, i1)
    {
        if (!visited[i1])
        {
            const location& l1 = locations[i1];

            if
            (
                !(l1.isSrcNotTgtPoint() && !isSrcEdge)
             && !(l1.isTgtNotSrcPoint() && isSrcEdge)
            )
            {
                const label ei1 =
                    isSrcEdge
                  ? (l1.isSrcPoint() ? (l1.srcPointi() + 2) % 3 : l1.srcEdgei())
                  : (l1.isTgtPoint() ? l1.tgtPointi() : l1.tgtEdgei());

                if (ei0 == ei1)
                {
                    auto branch = [&](const bool isSrcEdge)
                    {
                        return orderLocations
                        (
                            locations,
                            isSrcEdge,
                            i1,
                            nVisited,
                            visited,
                            order
                        );
                    };

                    if
                    (
                        (
                            !l1.isSrcAndTgtPoint()
                          && branch(l1.isIntersection() != isSrcEdge)
                        )
                     || (
                            l1.isSrcAndTgtPoint()
                         && (branch(true) || branch(false))
                        )
                    )
                    {
                        return true;
                    }
                }
            }
        }
    }

    // This branch failed to find a connected location. Un-visit this location.
    order[-- nVisited] = -1;
    visited[i0] = false;

    return false;
}


//- Construct the intersection topology
bool generateLocations
(
    const FixedList<label, 3>& tgtSrcPis,
    const FixedList<FixedList<label, 3>, 3>& srcInTgtEdge,
    const FixedList<FixedList<label, 3>, 3>& tgtInSrcEdge,
    const FixedList<label, 3>& srcInTgtTri,
    const FixedList<label, 3>& tgtInSrcTri,
    DynamicList<location>& pointLocations
)
{
    // Step 1: Process trivial rejection cases

    // If the entire target triangle is outside or on the same source edge
    // then there can be no intersection.
    forAll(srcInTgtEdge, srcEi)
    {
        bool outside = true;
        forAll(tgtInSrcEdge, tgtPi)
        {
            if (tgtInSrcEdge[tgtPi][srcEi] == 1)
            {
                outside = false;
                break;
            }
        }
        if (outside)
        {
            return true;
        }
    }

    // If all source points are outside all target edges this indicates
    // that the triangles are oppositely oriented, in which case there can
    // also be no intersection.
    if (count(srcInTgtEdge, {-1, -1, -1}) == 3)
    {
        return true;
    }

    // Step 2: Define point addition/checking functions

    // Add crossing point locations, inserting source points as necessary
    auto addPointLocations = [&pointLocations]
    (
        const location l1,
        const location l2 = location(),
        const bool add = true
    )
    {
        if (!pointLocations.empty())
        {
            const location l0 = pointLocations.last();

            if (l0.isIntersection() || l0.isSrcAndTgtPoint())
            {
                const label srcEi0 =
                    l0.isIntersection()
                  ? l0.srcEdgei()
                  : (l0.srcPointi() + 2) % 3;
                const label tgtEi0 =
                    l0.isIntersection()
                  ? l0.tgtEdgei()
                  : l0.tgtPointi();

                const label srcEi1 =
                    l1.isIntersection()
                  ? l1.srcEdgei()
                  : l1.srcPointi();
                const label tgtEi1 =
                    l1.isIntersection()
                  ? l1.tgtEdgei()
                  : (l1.tgtPointi() + 2) % 3;

                if
                (
                    (l0.isIntersection() && l1.isIntersection())
                 || tgtEi0 != tgtEi1
                )
                {
                    for
                    (
                        label srcEj = srcEi0;
                        srcEj != srcEi1;
                        srcEj = (srcEj + 2) % 3
                    )
                    {
                        pointLocations.append(location::srcPoint(srcEj));
                    }
                }
            }
        }

        if (add)
        {
            pointLocations.append(l1);

            if (!l2.isNull())
            {
                pointLocations.append(l2);
            }
        }
    };

    // One target point is within the source triangle and one is not
    auto inTriToOut = [&addPointLocations,&srcInTgtEdge]
    (
        const label tgtEi,
        const label tgtOutSrcEi1,
        const label tgtOutSrcPi1,
        const bool reverse
    )
    {
        const label srcEi =
            tgtOutSrcEi1 != -1
          ? tgtOutSrcEi1
          : srcInTgtEdge[tgtOutSrcPi1][tgtEi] == 1
          ? (tgtOutSrcPi1 + 2*reverse) % 3
          : (tgtOutSrcPi1 + 2*!reverse) % 3;

        addPointLocations(location::intersection(srcEi, tgtEi));
    };

    // One target point is a source point and the other is outside a source edge
    auto isPointToOutEdge = [&addPointLocations,&srcInTgtEdge]
    (
        const label tgtEi,
        const label tgtIsSrcPi0,
        const label tgtOutSrcEi1,
        const bool reverse
    )
    {
        const label srcEi0Next = (tgtIsSrcPi0 + 2*reverse) % 3;
        const label srcEi0Opp = (tgtIsSrcPi0 + 1) % 3;

        if
        (
            srcInTgtEdge[(tgtIsSrcPi0 + 1) % 3][tgtEi] == -1
         && srcInTgtEdge[(tgtIsSrcPi0 + 2) % 3][tgtEi] == -1
        )
        {
            return false;
        }

        if (tgtOutSrcEi1 == srcEi0Next)
        {
            return false;
        }

        if (tgtOutSrcEi1 == srcEi0Opp)
        {
            addPointLocations(location::intersection(srcEi0Opp, tgtEi));
        }

        return true;
    };

    // One target point is a source point and the other is outside a source
    // corner
    auto isPointToOutCorner = []
    (
        const label tgtEi,
        const label tgtIsSrcPi0,
        const label tgtOutSrcPi1,
        const bool reverse
    )
    {
        return tgtOutSrcPi1 != (tgtIsSrcPi0 + 1 + reverse) % 3;
    };

    // Both target points are outside source edges
    auto outEdgeToOutEdge = [&addPointLocations,&srcInTgtEdge]
    (
        const label tgtEi,
        const label tgtOutSrcEi0,
        const label tgtOutSrcEi1
    )
    {
        const label srcPi = (5 - tgtOutSrcEi0 - tgtOutSrcEi1) % 3;

        if
        (
            (tgtOutSrcEi0 != (tgtOutSrcEi1 + 1) % 3)
         && (srcInTgtEdge[srcPi][tgtEi] != 1)
        )
        {
            return false;
        }

        if
        (
            (tgtOutSrcEi0 == (tgtOutSrcEi1 + 1) % 3)
         != (srcInTgtEdge[srcPi][tgtEi] == 1)
        )
        {
            addPointLocations
            (
                location::intersection(tgtOutSrcEi0, tgtEi),
                location::intersection(tgtOutSrcEi1, tgtEi)
            );
        }

        return true;
    };

    // One target point is outside a source edge and the other is outside a
    // source corner
    auto outEdgeToOutCorner = [&addPointLocations,&srcInTgtEdge]
    (
        const label tgtEi,
        const label tgtOutSrcEi0,
        const label tgtOutSrcPi1,
        const bool reverse
    )
    {
        if (tgtOutSrcEi0 == tgtOutSrcPi1)
        {
            return !reverse;
        }

        if ((tgtOutSrcEi0 + 1) % 3 == tgtOutSrcPi1)
        {
            return reverse;
        }

        const label srcPi1Prev = (tgtOutSrcPi1 + 1 + !reverse) % 3;

        if (srcInTgtEdge[srcPi1Prev][tgtEi] == -1)
        {
            return false;
        }

        const label srcPi1Next = (tgtOutSrcPi1 + 1 + reverse) % 3;

        if (srcInTgtEdge[srcPi1Next][tgtEi] == 1)
        {
            return true;
        }

        location l1 = location::intersection(tgtOutSrcEi0, tgtEi);

        const label srcEi =
            srcInTgtEdge[tgtOutSrcPi1][tgtEi] == 1
          ? (tgtOutSrcPi1 + 2*reverse) % 3
          : (tgtOutSrcPi1 + 2*!reverse) % 3;

        location l2 = location::intersection(srcEi, tgtEi);

        if (reverse)
        {
            Swap(l1, l2);
        }

        addPointLocations(l1, l2);

        return true;
    };

    // Both target points are outside source corners
    auto outCornerToOutCorner = []
    (
        const label tgtEi,
        const label tgtOutSrcPi0,
        const label tgtOutSrcPi1
    )
    {
        return tgtOutSrcPi0 != (tgtOutSrcPi1 + 2) % 3;
    };

    // Step 3: Walk around the target edges to form the intersection polygon
    for (label tgtEi = 0; tgtEi < 3; tgtEi ++)
    {
        const label tgtPi0 = tgtEi, tgtPi1 = (tgtEi + 1) % 3;

        const bool tgtInSrcTri0 = tgtInSrcTri[tgtPi0] == 1;
        const bool tgtInSrcTri1 = tgtInSrcTri[tgtPi1] == 1;

        const label tgtIsSrcPi0 =
            tgtInSrcTri[tgtPi0] == 0 ? tgtSrcPis[tgtPi0] : -1;
        const label tgtIsSrcPi1 =
            tgtInSrcTri[tgtPi1] == 0 ? tgtSrcPis[tgtPi1] : -1;

        const label tgtOutSrcEi0 =
            count(tgtInSrcEdge[tgtPi0], -1) == 1
          ? findIndex(tgtInSrcEdge[tgtPi0], -1)
          : -1;
        const label tgtOutSrcEi1 =
            count(tgtInSrcEdge[tgtPi1], -1) == 1
          ? findIndex(tgtInSrcEdge[tgtPi1], -1)
          : -1;

        const label tgtOutSrcPi0 =
            count(tgtInSrcEdge[tgtPi0], -1) == 2
          ? (findIndex(tgtInSrcEdge[tgtPi0], 1) + 2) % 3
          : -1;
        const label tgtOutSrcPi1 =
            count(tgtInSrcEdge[tgtPi1], -1) == 2
          ? (findIndex(tgtInSrcEdge[tgtPi1], 1) + 2) % 3
          : -1;

        // Add the first point if it within or part of the source triangle
        if (tgtInSrcTri0)
        {
            pointLocations.append(location::tgtPoint(tgtPi0));
        }
        if (tgtIsSrcPi0 != -1)
        {
            addPointLocations(location::srcTgtPoint(tgtIsSrcPi0, tgtPi0));
        }

        // Add crossings
        if
        (
            (tgtInSrcTri0 && tgtInSrcTri1)
         || (tgtOutSrcEi0 != -1 && tgtOutSrcEi0 == tgtOutSrcEi1)
         || (tgtOutSrcPi0 != -1 && tgtOutSrcPi0 == tgtOutSrcPi1)
        )
        {
            // Both target points are in the same source quadrant. There is
            // nothing to check or to add.
        }
        else if
        (
            (tgtInSrcTri0 && tgtIsSrcPi1 != -1)
         || (tgtIsSrcPi0 != -1 && tgtInSrcTri1)
        )
        {
            // One target point is within the source triangle and one is a
            // source point. There is nothing to check or to add.
        }
        else if (tgtInSrcTri0 && (tgtOutSrcEi1 != -1 || tgtOutSrcPi1 != -1))
        {
            // The first target point is within the source triangle and the
            // second is outside
            inTriToOut(tgtEi, tgtOutSrcEi1, tgtOutSrcPi1, 0);
        }
        else if ((tgtOutSrcEi0 != -1 || tgtOutSrcPi0 != -1) && tgtInSrcTri1)
        {
            // (reverse of previous clause)
            inTriToOut(tgtEi, tgtOutSrcEi0, tgtOutSrcPi0, 1);
        }
        else if (tgtIsSrcPi0 != -1 && tgtIsSrcPi1 != -1)
        {
            // Both target points are source points. Check the ordering is
            // compatible with an intersection.
            if (tgtIsSrcPi0 != (tgtIsSrcPi1 + 1) % 3)
            {
                pointLocations.clear();
                return false;
            }
        }
        else if (tgtIsSrcPi0 != -1 && tgtOutSrcEi1 != -1)
        {
            // The first target point is a source point and the second is
            // outside a source edge
            if (!isPointToOutEdge(tgtEi, tgtIsSrcPi0, tgtOutSrcEi1, 0))
            {
                pointLocations.clear();
                return false;
            }
        }
        else if (tgtOutSrcEi0 != -1 && tgtIsSrcPi1 != -1)
        {
            // (reverse of previous clause)
            if (!isPointToOutEdge(tgtEi, tgtIsSrcPi1, tgtOutSrcEi0, 1))
            {
                pointLocations.clear();
                return false;
            }
        }
        else if (tgtIsSrcPi0 != -1 && tgtOutSrcPi1 != -1)
        {
            // The first target point is a source point and the second is
            // outside a source corner
            if (!isPointToOutCorner(tgtEi, tgtIsSrcPi0, tgtOutSrcPi1, 0))
            {
                pointLocations.clear();
                return false;
            }
        }
        else if (tgtOutSrcPi0 != -1 && tgtIsSrcPi1 != -1)
        {
            // (reverse of previous clause)
            if (!isPointToOutCorner(tgtEi, tgtIsSrcPi1, tgtOutSrcPi0, 1))
            {
                pointLocations.clear();
                return false;
            }
        }
        else if (tgtOutSrcEi0 != -1 && tgtOutSrcEi1 != -1)
        {
            // Both target points are outside source edges
            if (!outEdgeToOutEdge(tgtEi, tgtOutSrcEi0, tgtOutSrcEi1))
            {
                pointLocations.clear();
                return false;
            }
        }
        else if (tgtOutSrcEi0 != -1 && tgtOutSrcPi1 != -1)
        {
            // The first target point is outside a source edge and the
            // second is outside a source corner
            if (!outEdgeToOutCorner(tgtEi, tgtOutSrcEi0, tgtOutSrcPi1, 0))
            {
                pointLocations.clear();
                return false;
            }
        }
        else if (tgtOutSrcPi0 != -1 && tgtOutSrcEi1 != -1)
        {
            // (reverse of previous clause)
            if (!outEdgeToOutCorner(tgtEi, tgtOutSrcEi1, tgtOutSrcPi0, 1))
            {
                pointLocations.clear();
                return false;
            }
        }
        else if (tgtOutSrcPi0 != -1 && tgtOutSrcPi1 != -1)
        {
            // Both target points are outside source corners
            if (!outCornerToOutCorner(tgtEi, tgtOutSrcPi0, tgtOutSrcPi1))
            {
                pointLocations.clear();
                return false;
            }
        }
        else
        {
            // A target point is outside all source edges. The projection
            // has collapsed.
            pointLocations.clear();
            return false;
        }
    }

    // Step 4: Complete the polygon by adding any remaining source points that
    // were not traversed during the walk of the target edges
    if (!pointLocations.empty())
    {
        const location& l = pointLocations.first();

        if (l.isIntersection() || l.isSrcAndTgtPoint())
        {
            addPointLocations(l, location(), false);
        }
    }
    else
    {
        forAllReverse(srcInTgtEdge, srcPi)
        {
            pointLocations.append(location::srcPoint(srcPi));
        }
    }

    // Step 5: The above walk was done around the target triangle, but the
    // result should be ordered in the direction of the source triangle, so the
    // list of locations must be reversed
    inplaceReverseList(pointLocations);

    return true;
}


//- Construct the intersection geometry
void generateGeometry
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const FixedList<point, 3>& tgtPs,
    DynamicList<point>& srcPoints,
    DynamicList<vector>& srcPointNormals,
    DynamicList<point>& tgtPoints,
    const DynamicList<location>& pointLocations
)
{
    srcPoints.resize(pointLocations.size());
    srcPointNormals.resize(pointLocations.size());
    tgtPoints.resize(pointLocations.size());

    forAll(pointLocations, pointi)
    {
        const location& l = pointLocations[pointi];

        if (l.isSrcAndTgtPoint())
        {
            const point& srcP = srcPs[l.srcPointi()];
            const vector& srcN = srcNs[l.srcPointi()];
            const point& tgtP = tgtPs[l.tgtPointi()];

            srcPoints[pointi] = srcP;
            srcPointNormals[pointi] = srcN;
            tgtPoints[pointi] = tgtP;
        }
        else if (l.isSrcPoint())
        {
            const point& srcP = srcPs[l.srcPointi()];
            const vector& srcN = srcNs[l.srcPointi()];
            barycentric2D tgtTs =
                srcPointTgtTriIntersection(srcP, srcN, tgtPs);

            // Force inside the target triangle
            if (cmptMin(tgtTs) < 0)
            {
                const direction iMin = findMin(tgtTs);
                const direction iMax = findMax(tgtTs);
                const direction iMid = 3 - iMin - iMax;

                if (tgtTs[iMid] < 0)
                {
                    tgtTs[iMin] = 0;
                    tgtTs[iMax] = 1;
                    tgtTs[iMid] = 0;
                }
                else
                {
                    const scalar t = tgtTs[iMax] + tgtTs[iMid];
                    tgtTs[iMin] = 0;
                    tgtTs[iMax] /= t;
                    tgtTs[iMid] /= t;
                }
            }

            srcPoints[pointi] = srcP;
            srcPointNormals[pointi] = srcN;
            tgtPoints[pointi] = tgtTriInterpolate(tgtTs, tgtPs);
        }
        else if (l.isTgtPoint())
        {
            const point& tgtP = tgtPs[l.tgtPointi()];
            const barycentric2D srcTs =
                srcTriTgtPointIntersection(srcPs, srcNs, tgtP);

            srcPoints[pointi] = srcTriInterpolate(srcTs, srcPs);
            srcPointNormals[pointi] = srcTriInterpolate(srcTs, srcNs);
            tgtPoints[pointi] = tgtP;
        }
        else // if (l.isIntersection())
        {
            const label srcPi0 = l.srcEdgei(), srcPi1 = (srcPi0 + 1) % 3;
            const label tgtPi0 = l.tgtEdgei(), tgtPi1 = (tgtPi0 + 1) % 3;

            const Pair<scalar> ts =
                srcEdgeTgtEdgeIntersection
                (
                    {srcPs[srcPi0], srcPs[srcPi1]},
                    {srcNs[srcPi0], srcNs[srcPi1]},
                    {tgtPs[tgtPi0], tgtPs[tgtPi1]}
                );
            const scalar srcT = ts.first(), tgtT = ts.second();

            srcPoints[pointi] =
                (1 - srcT)*srcPs[srcPi0] + srcT*srcPs[srcPi1];
            srcPointNormals[pointi] =
                (1 - srcT)*srcNs[srcPi0] + srcT*srcNs[srcPi1];
            tgtPoints[pointi] =
                (1 - tgtT)*tgtPs[tgtPi0] + tgtT*tgtPs[tgtPi1];
        }
    }
}


} // End namespace triIntersect
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::triIntersect::writeTriProjection
(
    const word& name,
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const label nEdge,
    const label nNormal,
    const scalar lNormal
)
{
    scalar lengthScale = 0;
    for (label i = 0; i < 3; ++ i)
    {
        lengthScale = max(lengthScale, mag(srcPs[i] - srcPs[(i + 1) % 3]));
    }

    const label nu = nEdge, nv = nNormal;
    const scalar u0 = 0, u1 = 1;
    const scalar v0 = -lNormal/2*lengthScale, v1 = lNormal/2*lengthScale;

    pointField ps(3*(nu + 1)*(nv + 1));
    for (label i = 0; i < 3; ++ i)
    {
        const point& p0 = srcPs[i], p1 = srcPs[(i + 1) % 3];
        const vector& n0 = srcNs[i], n1 = srcNs[(i + 1) % 3];

        for (label iu = 0; iu <= nu; ++ iu)
        {
            const scalar u = u0 + (u1 - u0)*scalar(iu)/nu;
            for (label iv = 0; iv <= nv; ++ iv)
            {
                const scalar v = v0 + (v1 - v0)*scalar(iv)/nv;
                const vector x = p0 + (p1 - p0)*u + (n0 + (n1 - n0)*u)*v;
                ps[i*(nu + 1)*(nv + 1) + iu*(nv + 1) + iv] = x;
            }
        }
    }

    faceList fs(3*nu*nv);
    for (label i = 0; i < 3; ++ i)
    {
        for (label iu = 0; iu < nu; ++ iu)
        {
            for (label iv = 0; iv < nv; ++ iv)
            {
                fs[i*nu*nv + iu*nv + iv] =
                    face
                    ({
                        i*(nu + 1)*(nv + 1) + (nv + 1)*iu + iv,
                        i*(nu + 1)*(nv + 1) + (nv + 1)*iu + iv + 1,
                        i*(nu + 1)*(nv + 1) + (nv + 1)*(iu + 1) + iv + 1,
                        i*(nu + 1)*(nv + 1) + (nv + 1)*(iu + 1) + iv
                    });
            }
        }
    }

    Info<< indent << "Writing face to " << name + ".vtk" << endl;
    vtkWritePolyData::write
    (
        name + ".vtk",
        name,
        false,
        ps,
        labelList(),
        labelListList(),
        fs
    );
}


Foam::scalar Foam::triIntersect::srcEdgeTgtPointOffset
(
    const Pair<point>& srcPs,
    const Pair<vector>& srcNs,
    const point& tgtP
)
{
    const tensor A(srcPs[1] - srcPs[0], srcNs[0], srcNs[1]);
    const scalar detA = det(A);

    const tensor T(A.y()^A.z(), A.z()^A.x(), A.x()^A.y());

    const vector detAY = T & (tgtP - srcPs[0]);
    const scalar Yx = protectedDivideAndClip01(detAY.x(), detA);

    const scalar offset = Yx*detAY.y() - (1 - Yx)*detAY.z();

    return offset == 0 ? - vSmall : offset;
}


Foam::scalar Foam::triIntersect::srcEdgeTgtPointOffset
(
    const Pair<point>& srcPs,
    const Pair<vector>& srcNs,
    const point& tgtP,
    const bool srcDirection
)
{
    if (srcDirection)
    {
        return srcEdgeTgtPointOffset(srcPs, srcNs, tgtP);
    }
    else
    {
        return - srcEdgeTgtPointOffset(reverse(srcPs), reverse(srcNs), tgtP);
    }
}


Foam::scalar Foam::triIntersect::srcPointTgtEdgeOffset
(
    const point& srcP,
    const vector& srcN,
    const Pair<point>& tgtPs
)
{
    const tensor A(srcN, tgtPs[1] - tgtPs[0], srcN^(tgtPs[1] - tgtPs[0]));
    const scalar detA = det(A);

    const tensor T(A.y()^A.z(), A.z()^A.x(), A.x()^A.y());

    const vector detAY = T & (tgtPs[0] - srcP);

    const scalar offset = protectedDivide(detAY.z(), detA);

    return offset == 0 ? - vSmall : offset;
}


Foam::scalar Foam::triIntersect::srcPointTgtEdgeOffset
(
    const point& srcP,
    const vector& srcN,
    const Pair<point>& tgtPs,
    const bool tgtDirection
)
{
    if (tgtDirection)
    {
        return srcPointTgtEdgeOffset(srcP, srcN, tgtPs);
    }
    else
    {
        return - srcPointTgtEdgeOffset(srcP, srcN, reverse(tgtPs));
    }
}


Foam::Pair<Foam::scalar> Foam::triIntersect::srcEdgeTgtEdgeIntersection
(
    const Pair<point>& srcPs,
    const Pair<vector>& srcNs,
    const Pair<point>& tgtPs
)
{
    scalar srcT, tgtT;

    #ifndef BISECT
    const Tuple2<bool, vector> solution =
        solveProjection
        (
            srcPs[0] - tgtPs[0],
            srcPs[1] - srcPs[0],
            tgtPs[0] - tgtPs[1],
            srcNs[0],
            Zero,
            srcNs[1] - srcNs[0],
            {0, 1, -1}
        );
    if (solution.first())
    {
        // If the analytical solution succeeds, then use the result
        srcT = solution.second().x();
        tgtT = solution.second().y();
    }
    else
    #endif
    {
        // If the analytical solution fails, then solve by bisection

        // !!! This method, whilst elegant, isn't sufficiently robust. The
        // srcPointTgtEdgeOffset calculation is not as reliable as
        // srcEdgeTgtPointOffset. So, we can't bisect the source edge to get
        // srcT and then reuse the solveProjectionGivenT stuff. We need to
        // bisect the target edge to get tgtT, and then have a specific process
        // for calculating srcT from tgtT.
        /*
        scalar srcT0 = 0, srcT1 = 1;

        const scalar o0 = srcPointTgtEdgeOffset(srcPs[0], srcNs[0], tgtPs);
        const scalar o1 = srcPointTgtEdgeOffset(srcPs[1], srcNs[1], tgtPs);

        const scalar s = o0 > o1 ? +1 : -1;

        for (label i = 0; i < ceil(std::log2(1/small)); ++ i)
        {
            const scalar srcT = (srcT0 + srcT1)/2;
            const vector srcP = (1 - srcT)*srcPs[0] + srcT*srcPs[1];
            const vector srcN = (1 - srcT)*srcNs[0] + srcT*srcNs[1];

            const scalar o = s*srcPointTgtEdgeOffset(srcP, srcN, tgtPs);

            if (o > 0)
            {
                srcT0 = srcT;
            }
            else
            {
                srcT1 = srcT;
            }
        }

        srcT = (srcT0 + srcT1)/2;
        tgtT =
            solveProjectionGivenT
            (
                srcPs[0] - tgtPs[0],
                srcPs[1] - srcPs[0],
                tgtPs[0] - tgtPs[1],
                srcNs[0],
                Zero,
                srcNs[1] - srcNs[0],
                {0, 1, -1},
                srcT
            ).y();
        */

        // !!! This method appears robust
        scalar tgtT0 = 0, tgtT1 = 1;

        const scalar o0 = srcEdgeTgtPointOffset(srcPs, srcNs, tgtPs[0]);
        const scalar o1 = srcEdgeTgtPointOffset(srcPs, srcNs, tgtPs[1]);

        const scalar s = o0 > o1 ? +1 : -1;

        for (label i = 0; i < ceil(std::log2(1/small)); ++ i)
        {
            const scalar tgtT = (tgtT0 + tgtT1)/2;
            const vector tgtP = (1 - tgtT)*tgtPs[0] + tgtT*tgtPs[1];

            const scalar o = s*srcEdgeTgtPointOffset(srcPs, srcNs, tgtP);

            if (o > 0)
            {
                tgtT0 = tgtT;
            }
            else
            {
                tgtT1 = tgtT;
            }
        }

        tgtT = (tgtT0 + tgtT1)/2;

        // Solve for the corresponding source edge coordinate
        const vector srcDP = srcPs[1] - srcPs[0];
        const vector tgtP = (1 - tgtT)*tgtPs[0] + tgtT*tgtPs[1];

        const tensor A(srcDP, srcNs[0], srcNs[1]);
        const scalar detA = det(A);

        const vector Tx = A.y()^A.z();

        const scalar magDetAYx = sign(detA)*(Tx & (tgtP - srcPs[0]));

        const scalar srcTStar = protectedDivideAndClip01(magDetAYx, mag(detA));

        const vector srcN = (1 - srcTStar)*srcNs[0] + srcTStar*srcNs[1];

        const vector srcDPPerpN = srcDP - (srcDP & srcN)*srcN;

        srcT =
            protectedDivideAndClip01
            (
                (tgtP - srcPs[0]) & srcDPPerpN,
                srcDP & srcDPPerpN
            );
    }

    /*
    // Check that the points match
    {
        const point srcP = (1 - srcT)*srcPs[0] + srcT*srcPs[1];
        const point srcN = (1 - srcT)*srcNs[0] + srcT*srcNs[1];
        const point tgtP = (1 - tgtT)*tgtPs[0] + tgtT*tgtPs[1];
        const scalar srcU = ((tgtP - srcP) & srcN)/magSqr(srcN);
        const point srcQ = srcP + srcU*srcN;
        Info<< "srcT=" << srcT << ", tgtT=" << tgtT
            << ", err=" << magSqr(srcQ - tgtP) << endl;
    }
    */

    return Pair<scalar>(srcT, tgtT);
}


Foam::barycentric2D Foam::triIntersect::srcTriTgtPointIntersection
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const point& tgtP
)
{
    auto srcPN = []
    (
        const FixedList<vector, 3>& srcPNs,
        const vector2D& srcTs
    )
    {
        const scalar srcT0 = 1 - srcTs.x() - srcTs.y();
        return srcT0*srcPNs[0] + srcTs.x()*srcPNs[1] + srcTs.y()*srcPNs[2];
    };

    auto srcEdgePNs = [srcPN]
    (
        const FixedList<vector, 3>& srcPNs,
        const vector2D& srcTs0,
        const vector2D& srcTs1
    )
    {
        return Pair<vector>(srcPN(srcPNs, srcTs0), srcPN(srcPNs, srcTs1));
    };

    auto offset = [&](const vector2D& srcTs0, const vector2D& srcTs1)
    {
        return srcEdgeTgtPointOffset
        (
            srcEdgePNs(srcPs, srcTs0, srcTs1),
            srcEdgePNs(srcNs, srcTs0, srcTs1),
            tgtP
        );
    };

    const scalar oA = offset(vector2D(0, 0), vector2D(1, 0));
    const scalar oB = offset(vector2D(1, 0), vector2D(0, 1));
    const scalar oC = offset(vector2D(0, 1), vector2D(0, 0));
    const FixedList<scalar, 3> offsets({oA, oB, oC});

    // If inside the triangle (or outside if the triangle is inverted) ...
    if (offsets[findMin(offsets)] >= 0 || offsets[findMax(offsets)] <= 0)
    {
        scalar srcT, srcU;

        #ifndef BISECT
        const Tuple2<bool, vector> solution =
            solveProjection
            (
                srcPs[0] - tgtP,
                srcNs[0],
                srcPs[1] - srcPs[0],
                srcPs[2] - srcPs[0],
                srcNs[1] - srcNs[0],
                srcNs[2] - srcNs[0],
                {-1, 0, 0}
            );
        if (solution.first())
        {
            // If the analytical solution succeeds, then use the result
            srcT = solution.second().y();
            srcU = solution.second().z();
        }
        else
        #endif
        {
            // If the analytical solution fails, then solve by bisection
            const scalar sign = offsets[findMin(offsets)] > 0 ? +1 : -1;

            vector2D srcTsA(0, 0), srcTsB(1, 0), srcTsC(0, 1);

            for (label i = 0; i < ceil(std::log2(1/small)); ++ i)
            {
                const vector2D srcTsAB = (srcTsA + srcTsB)/2;
                const vector2D srcTsBC = (srcTsB + srcTsC)/2;
                const vector2D srcTsCA = (srcTsC + srcTsA)/2;

                const scalar oA = sign*offset(srcTsCA, srcTsAB);
                const scalar oB = sign*offset(srcTsAB, srcTsBC);
                const scalar oC = sign*offset(srcTsBC, srcTsCA);
                const FixedList<scalar, 3> offsets({oA, oB, oC});

                const label offsetMini = findMin(offsets);
                if (offsets[offsetMini] > 0)
                {
                    srcTsA = srcTsAB;
                    srcTsB = srcTsBC;
                    srcTsC = srcTsCA;
                }
                else if (offsetMini == 0)
                {
                    srcTsC = srcTsCA;
                    srcTsB = srcTsAB;
                }
                else if (offsetMini == 1)
                {
                    srcTsA = srcTsAB;
                    srcTsC = srcTsBC;
                }
                else if (offsetMini == 2)
                {
                    srcTsB = srcTsBC;
                    srcTsA = srcTsCA;
                }
            }

            srcT = (srcTsA[0] + srcTsB[0] + srcTsC[0])/3;
            srcU = (srcTsA[1] + srcTsB[1] + srcTsC[1])/3;
        }

        return barycentric2D(1 - srcT - srcU, srcT, srcU);
    }

    // If outside an edge ...
    forAll(srcPs, srcEi)
    {
        const label srcEi0 = (srcEi + 2) % 3, srcEi1 = (srcEi + 1) % 3;

        if
        (
            offsets[srcEi] <= 0
         && offsets[srcEi0] >= 0
         && offsets[srcEi1] >= 0
        )
        {
            const label srcPi0 = srcEi, srcPi1 = (srcEi + 1) % 3;
            const label srcPiOpp = (srcEi + 2) % 3;

            scalar srcT, srcU;

            #ifndef BISECT
            const Tuple2<bool, vector> solution =
                solveProjection
                (
                    srcPs[srcPi0] - tgtP,
                    srcPs[srcPi1] - srcPs[srcPi0],
                    srcPs[srcPi0] - srcPs[srcPiOpp],
                    srcNs[srcPi0],
                    srcPs[srcPi1] - srcPs[srcPi0],
                    srcNs[srcPi1] - srcNs[srcPi0],
                    {0, -1, -1}
                );
            if (solution.first())
            {
                // If the analytical solution succeeds, then use the result
                srcT = solution.second().x();
                srcU = solution.second().y();
            }
            else
            #endif
            {
                // If the analytical solution fails, then solve by bisection
                const vector2D srcTsOpp(srcPiOpp == 1, srcPiOpp == 2);
                const vector2D srcTs0(srcPi0 == 1, srcPi0 == 2);
                const vector2D srcTs1(srcPi1 == 1, srcPi1 == 2);

                scalar srcT0 = 0, srcT1 = 1;

                for (label i = 0; i < ceil(std::log2(1/small)); ++ i)
                {
                    const scalar srcT = (srcT0 + srcT1)/2;
                    const vector2D srcTs01(srcTs0*(1 - srcT) + srcTs1*srcT);

                    const scalar o = offset(srcTsOpp, srcTs01);

                    if (o > 0)
                    {
                        srcT0 = srcT;
                    }
                    else
                    {
                        srcT1 = srcT;
                    }
                }

                srcT = (srcT0 + srcT1)/2;
                srcU =
                    solveProjectionGivenT
                    (
                        srcPs[srcPi0] - tgtP,
                        srcPs[srcPi1] - srcPs[srcPi0],
                        srcPs[srcPi0] - srcPs[srcPiOpp],
                        srcNs[srcPi0],
                        srcPs[srcPi1] - srcPs[srcPi0],
                        srcNs[srcPi1] - srcNs[srcPi0],
                        {0, -1, -1},
                        srcT
                    ).y();
            }

            // Convert to the triangle's coordinate system
            barycentric2D y;
            y[srcPiOpp] = - srcU;
            y[srcPi0] = (1 + srcU)*(1 - srcT);
            y[srcPi1] = (1 + srcU)*srcT;
            return y;
        }
    }

    // If outside a corner ...
    forAll(srcPs, srcPi)
    {
        const label srcEiOpp = (srcPi + 1) % 3;
        const label srcEi0 = (srcPi + 2) % 3, srcEi1 = srcPi;

        if
        (
            offsets[srcEiOpp] >= 0
         && offsets[srcEi0] <= 0
         && offsets[srcEi1] <= 0
        )
        {
            // Solve for the intersection coordinates directly
            const label srcPi0 = (srcPi + 2) % 3, srcPi1 = (srcPi + 1) % 3;

            const tensor A
            (
                srcPs[srcPi] - srcPs[srcPi0],
                srcPs[srcPi] - srcPs[srcPi1],
                srcNs[srcPi]
            );
            const vector T0(A.y()^A.z()), T1(A.z()^A.x());
            const scalar detA = A.x() & T0;

            const scalar srcT =
                protectedDivide(T0 & (tgtP - srcPs[srcPi]), detA);
            const scalar srcU =
                protectedDivide(T1 & (tgtP - srcPs[srcPi]), detA);

            // Convert to the triangle's coordinate system
            barycentric2D y;
            y[srcPi0] = - srcT;
            y[srcPi1] = - srcU;
            y[srcPi] = 1 + srcT + srcU;
            return y;
        }
    }

    // Above logic means we should never reach here
    FatalErrorInFunction
        << "Point " << tgtP << " could not be classified within triangle "
        << srcPs << " with projection normals " << srcNs << exit(FatalError);
    return barycentric2D::uniform(NaN);
}


Foam::barycentric2D Foam::triIntersect::srcPointTgtTriIntersection
(
    const point& srcP,
    const vector& srcN,
    const FixedList<point, 3>& tgtPs
)
{
    const tensor A(tgtPs[1] - tgtPs[0], tgtPs[2] - tgtPs[0], - srcN);
    const scalar detA = det(A);

    const vector T0(A.y()^A.z()), T1(A.z()^A.x());
    const tensor T(- T0 - T1, T0, T1);

    const vector detAY = (T & (srcP - tgtPs[0])) + vector(detA, 0, 0);

    const scalar maxMagDetAY = mag(detAY[findMax(cmptMag(detAY))]);

    vector y =
        maxMagDetAY/vGreat < mag(detA)
      ? detAY/detA
      : detAY/maxMagDetAY*vGreat;

    return barycentric2D(y.x(), y.y(), y.z());
}


void Foam::triIntersect::intersectTris
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const FixedList<bool, 3>& srcOwns,
    const FixedList<label, 3>& srcTgtPis,
    const FixedList<point, 3>& tgtPs,
    const FixedList<bool, 3>& tgtOwns,
    const FixedList<label, 3>& tgtSrcPis,
    DynamicList<point>& srcPoints,
    DynamicList<vector>& srcPointNormals,
    DynamicList<point>& tgtPoints,
    DynamicList<location>& pointLocations,
    const bool debug,
    const word& writePrefix
)
{
    const bool write = writePrefix != word::null;

    if (debug || write)
    {
        Info<< indent << "Intersecting triangles" << incrIndent << endl;
    }

    if (write)
    {
        writePolygon(writePrefix + "_srcTri", srcPs);
        writePolygon(writePrefix + "_tgtTri", tgtPs);
        writeTriProjection(writePrefix + "_srcPrj", srcPs, srcNs);
    }

    // Determine what source points lie within target edges and vice-versa
    FixedList<FixedList<label, 3>, 3> srcInTgtEdge, tgtInSrcEdge;
    FixedList<label, 3> srcInTgtTri, tgtInSrcTri;
    srcInTgt
    (
        srcPs, srcNs, srcOwns, srcTgtPis,
        tgtPs, tgtOwns,
        srcInTgtEdge,
        srcInTgtTri
    );
    tgtInSrc
    (
        srcPs, srcNs, srcOwns,
        tgtPs, tgtOwns, tgtSrcPis,
        tgtInSrcEdge,
        tgtInSrcTri
    );

    if (debug)
    {
        if (count(srcTgtPis, -1) != 3)
        {
            Info<< indent << "srcTgtPis=" << srcTgtPis << endl;
        }
        Info<< indent << "srcInTgtTri=" << srcInTgtTri << endl
            << indent << "srcInTgtEdge=" << srcInTgtEdge << endl;
        if (count(tgtSrcPis, -1) != 3)
        {
            Info<< indent << "tgtSrcPis=" << tgtSrcPis << endl;
        }
        Info<< indent << "tgtInSrcTri=" << tgtInSrcTri << endl
            << indent << "tgtInSrcEdge=" << tgtInSrcEdge << endl;
    }

    // Generate the locations
    generateLocations
    (
        tgtSrcPis,
        srcInTgtEdge,
        tgtInSrcEdge,
        srcInTgtTri,
        tgtInSrcTri,
        pointLocations
    );

    // Generate the geometry
    generateGeometry
    (
        srcPs,
        srcNs,
        tgtPs,
        srcPoints,
        srcPointNormals,
        tgtPoints,
        pointLocations
    );

    if (write)
    {
        writePolygon(writePrefix + "_srcIctFace", srcPoints);
        writePolygon(writePrefix + "_tgtIctFace", tgtPoints);
    }

    if (debug || write)
    {
        Info<< decrIndent;
    }
}

// ************************************************************************* //
