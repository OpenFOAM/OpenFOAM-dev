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

#include "intersectionPatchToPatch.H"
#include "triIntersect.H"
#include "vtkWritePolyData.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace patchToPatches
{
    defineTypeNameAndDebug(intersection, 0);
    addToRunTimeSelectionTable(patchToPatch, intersection, bool);

    int intersection::debugSrcFacei =
        debug::debugSwitch((intersection::typeName + "SrcFace").c_str(), -1);
    int intersection::debugTgtFacei =
        debug::debugSwitch((intersection::typeName + "TgtFace").c_str(), -1);
}
}


// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

template<class Type>
Foam::FixedList<Type, 3>
Foam::patchToPatches::intersection::triPointValues
(
    const triFace& t,
    const UList<Type>& values
)
{
    FixedList<Type, 3> result;
    forAll(t, i)
    {
        result[i] = values[t[i]];
    }
    return result;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::treeBoundBox Foam::patchToPatches::intersection::srcBoxStatic
(
    const face& srcFace,
    const pointField& srcPoints,
    const vectorField& srcPointNormals
)
{
    static DynamicList<point> ps;

    ps.clear();

    const scalar l = sqrt(mag(srcFace.area(srcPoints)));

    forAll(srcFace, srcFacePointi)
    {
        const label srcPointi = srcFace[srcFacePointi];

        const point& p = srcPoints[srcPointi];
        const vector& n = srcPointNormals[srcPointi];

        ps.append(p - l/2*n);
        ps.append(p + l/2*n);
    }

    return treeBoundBox(ps);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::patchToPatches::intersection::srcBox
(
    const face& srcFace,
    const pointField& srcPoints,
    const vectorField& srcPointNormals
) const
{
    return srcBoxStatic(srcFace, srcPoints, srcPointNormals);
}


bool Foam::patchToPatches::intersection::intersectFaces
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const label srcFacei,
    const label tgtFacei
)
{
    // Quick rejection based on bound box
    const treeBoundBox srcFaceBox =
        srcBox
        (
            srcPatch.localFaces()[srcFacei],
            srcPatch.localPoints(),
            srcPointNormals
        );
    const treeBoundBox tgtFaceBox(tgtPatch.points(), tgtPatch[tgtFacei]);
    if (!srcFaceBox.overlaps(tgtFaceBox)) return false;

    // Construct face triangulations on demand
    if (srcTriPoints_[srcFacei].empty())
    {
        triEngine_.triangulate
        (
            UIndirectList<point>
            (
                srcPatch.localPoints(),
                srcPatch.localFaces()[srcFacei]
            )
        );

        srcTriPoints_[srcFacei] =
            triEngine_.triPoints(srcPatch.localFaces()[srcFacei]);
        srcTriFaceEdges_[srcFacei] = triEngine_.triEdges();
    }
    if (tgtTriPoints_[tgtFacei].empty())
    {
        triEngine_.triangulate
        (
            UIndirectList<point>
            (
                tgtPatch.localPoints(),
                tgtPatch.localFaces()[tgtFacei]
            )
        );

        tgtTriPoints_[tgtFacei] =
            triEngine_.triPoints(tgtPatch.localFaces()[tgtFacei]);
        tgtTriFaceEdges_[tgtFacei] = triEngine_.triEdges();
    }

    // Construct and initialise workspace
    bool srcCouples = false;
    couple srcCouple;
    srcFaceEdgePart_.resize(srcPatch[srcFacei].size());
    forAll(srcFaceEdgePart_, srcFaceEdgei)
    {
        const edge e =
            srcPatch.localFaces()[srcFacei].faceEdge(srcFaceEdgei);
        const vector eC = e.centre(srcPatch.localPoints());
        srcFaceEdgePart_[srcFaceEdgei] = part(Zero, eC);
    }

    bool tgtCouples = false;
    couple tgtCouple;
    tgtFaceEdgePart_.resize(tgtPatch[tgtFacei].size());
    forAll(tgtFaceEdgePart_, tgtFaceEdgei)
    {
        const edge e =
            tgtPatch.localFaces()[tgtFacei].faceEdge(tgtFaceEdgei);
        const vector eC = e.centre(tgtPatch.localPoints());
        tgtFaceEdgePart_[tgtFaceEdgei] = part(Zero, eC);
    }

    part errorPart(Zero, srcPatch.faceCentres()[srcFacei]);

    // Cache the face area magnitudes
    const scalar srcMagA = mag(srcPatch.faceAreas()[srcFacei]);
    const scalar tgtMagA = mag(tgtPatch.faceAreas()[tgtFacei]);

    // Determine whether or not to debug this tri intersection
    const bool debugTriIntersect =
        (debugSrcFacei != -1 || debugTgtFacei != -1)
     && (debugSrcFacei == -1 || debugSrcFacei == srcFacei)
     && (debugTgtFacei == -1 || debugTgtFacei == tgtFacei);

    // Loop the face triangles and compute the intersections
    bool anyCouples = false;
    forAll(srcTriPoints_[srcFacei], srcFaceTrii)
    {
        const triFace& srcT = srcTriPoints_[srcFacei][srcFaceTrii];

        const FixedList<point, 3> srcPs =
            triPointValues(srcT, srcPatch.localPoints());
        const FixedList<vector, 3> srcNs =
            triPointValues(srcT, srcPointNormals);

        forAll(tgtTriPoints_[tgtFacei], tgtFaceTrii)
        {
            const triFace tgtT =
                reverse()
              ? tgtTriPoints_[tgtFacei][tgtFaceTrii].reverseFace()
              : tgtTriPoints_[tgtFacei][tgtFaceTrii];

            const FixedList<point, 3> tgtPs =
                triPointValues(tgtT, tgtPatch.localPoints());

            // Do tri-intersection
            ictSrcPoints_.clear();
            ictSrcPointNormals_.clear();
            ictTgtPoints_.clear();
            ictPointLocations_.clear();
            triIntersect::intersectTris
            (
                srcPs,
                srcNs,
                {false, false, false},
                {-1, -1, -1},
                tgtPs,
                {false, false, false},
                {-1, -1, -1},
                ictSrcPoints_,
                ictSrcPointNormals_,
                ictTgtPoints_,
                ictPointLocations_,
                debugTriIntersect,
                debugTriIntersect
              ? word
                (
                    typeName
                  + "_srcFace=" + Foam::name(srcFacei)
                  + "_tgtFace=" + Foam::name(tgtFacei)
                  + "_intersection=" + Foam::name
                    (srcFaceTrii*tgtTriPoints_[tgtFacei].size() + tgtFaceTrii)
                )
              : word::null
            );

            // If there is no intersection then continue
            if (ictPointLocations_.empty())
            {
                continue;
            }

            // Mark that there has been an intersection
            anyCouples = true;

            // Compute the intersection geometry
            const part ictSrcPart(ictSrcPoints_);
            const part ictTgtPart(ictTgtPoints_);

            // If the intersection is below tolerance then continue
            if
            (
                mag(ictSrcPart.area) < small*srcMagA
             || mag(ictTgtPart.area) < small*tgtMagA
            )
            {
                continue;
            }

            // Mark that the source and target faces intersect
            srcCouples = tgtCouples = true;

            // Store the intersection geometry
            srcCouple += ictSrcPart;
            srcCouple.nbr += ictTgtPart;
            if (reverse())
            {
                tgtCouple += ictTgtPart;
                tgtCouple.nbr += ictSrcPart;
            }
            else
            {
                tgtCouple -= ictTgtPart;
                tgtCouple.nbr -= ictSrcPart;
            }

            // Store the intersection polygons for debugging
            const label debugSrcPoint0 = debugPoints_.size();
            const label debugTgtPoint0 =
                debugPoints_.size() + ictSrcPoints_.size();
            if (debug)
            {
                debugPoints_.append(ictSrcPoints_);
                debugPoints_.append(ictTgtPoints_);
                debugFaces_.append
                (
                    debugSrcPoint0 + identityMap(ictSrcPoints_.size())
                );
                debugFaceSrcFaces_.append(srcFacei);
                debugFaceTgtFaces_.append(tgtFacei);
                debugFaceSides_.append(1);
                debugFaces_.append
                (
                    debugTgtPoint0 + identityMap(ictTgtPoints_.size())
                );
                debugFaceSrcFaces_.append(srcFacei);
                debugFaceTgtFaces_.append(tgtFacei);
                debugFaceSides_.append(-1);
            }

            // Store edge and error areas
            forAll(ictPointLocations_, i0)
            {
                const label i1 = ictPointLocations_.fcIndex(i0);

                // Get the locations on each end of this edge of the
                // intersection polygon
                const triIntersect::location l0 = ictPointLocations_[i0];
                const triIntersect::location l1 = ictPointLocations_[i1];

                // Get the geometry for the projection of this edge
                const part ictEdgePart
                (
                    FixedList<point, 4>
                    ({
                        ictSrcPoints_[i0],
                        ictSrcPoints_[i1],
                        ictTgtPoints_[i1],
                        ictTgtPoints_[i0]
                    })
                );

                // Store the "side" of the intersection that this edge
                // corresponds to
                label ictEdgeSide = -labelMax;

                // If this edge corresponds to an edge of the source
                // triangle
                if
                (
                    l0.isSrcNotTgtPoint()
                 || l1.isSrcNotTgtPoint()
                 || (
                        l0.isIntersection()
                     && l1.isIntersection()
                     && l0.srcEdgei() == l1.srcEdgei()
                    )
                )
                {
                    const label srcEi =
                        l0.isSrcPoint() ? l0.srcPointi()
                      : l1.isSrcPoint() ? (l1.srcPointi() + 2) % 3
                      : l0.srcEdgei();

                    const label srcFaceEdgei =
                        srcTriFaceEdges_[srcFacei][srcFaceTrii][srcEi];

                    if (srcFaceEdgei < srcPatch[srcFacei].size())
                    {
                        srcFaceEdgePart_[srcFaceEdgei] += ictEdgePart;
                        ictEdgeSide = 1;
                    }
                    else
                    {
                        errorPart += ictEdgePart;
                        ictEdgeSide = 0;
                    }
                }

                // If this edge corresponds to an edge of the target
                // triangle
                else if
                (
                    l0.isTgtNotSrcPoint()
                 || l1.isTgtNotSrcPoint()
                 || (
                        l0.isIntersection()
                     && l1.isIntersection()
                     && l0.tgtEdgei() == l1.tgtEdgei()
                    )
                )
                {
                    const label tgtEi =
                        l0.isTgtPoint() ? (l0.tgtPointi() + 2) % 3
                      : l1.isTgtPoint() ? l1.tgtPointi()
                      : l0.tgtEdgei();

                    const label tgtFaceEdgei =
                        tgtTriFaceEdges_[tgtFacei][tgtFaceTrii]
                        [reverse() ? 2 - tgtEi : tgtEi];

                    if (tgtFaceEdgei < tgtPatch[tgtFacei].size())
                    {
                        tgtFaceEdgePart_[tgtFaceEdgei] +=
                            reverse() ? ictEdgePart : -ictEdgePart;
                        ictEdgeSide = -1;
                    }
                    else
                    {
                        errorPart += ictEdgePart;
                        ictEdgeSide = 0;
                    }
                }

                // No other location combinations should be possible for an
                // intersection without any shared points
                else
                {
                    FatalErrorInFunction
                        << "Tri-intersection topology not recognised. "
                        << "This is a bug." << exit(FatalError);
                }

                // Store the projected edge quadrilateral for debugging
                if (debug)
                {
                    debugFaces_.append
                    (
                        labelList
                        ({
                            debugSrcPoint0 + i0,
                            debugSrcPoint0 + i1,
                            debugTgtPoint0 + i1,
                            debugTgtPoint0 + i0
                        })
                    );
                    debugFaceSrcFaces_.append(srcFacei);
                    debugFaceTgtFaces_.append(tgtFacei);
                    debugFaceSides_.append(ictEdgeSide);
                }
            }
        }
    }

    // If the source face couples the target, then store the intersection
    if (srcCouples)
    {
        srcLocalTgtFaces_[srcFacei].append(tgtFacei);
        srcCouples_[srcFacei].append(srcCouple);
    }

    // If any intersection has occurred then store the edge and error parts
    if (anyCouples)
    {
        forAll(srcFaceEdgeParts_[srcFacei], srcFaceEdgei)
        {
            srcFaceEdgeParts_[srcFacei][srcFaceEdgei] +=
                srcFaceEdgePart_[srcFaceEdgei];
        }
        srcErrorParts_[srcFacei] +=
            reverse() ? sum(tgtFaceEdgePart_) : -sum(tgtFaceEdgePart_);
        srcErrorParts_[srcFacei] += errorPart;
    }

    // If the target face couples the source, then store in the intersection
    if (tgtCouples)
    {
        tgtLocalSrcFaces_[tgtFacei].append(srcFacei);
        tgtCouples_[tgtFacei].append(tgtCouple);
    }

    return anyCouples;
}


void Foam::patchToPatches::intersection::initialise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    patchToPatch::initialise
    (
        srcPatch,
        srcPointNormals,
        srcPointNormals0,
        tgtPatch
    );

    srcCouples_.resize(srcPatch.size());
    forAll(srcLocalTgtFaces_, i)
    {
        srcCouples_[i].clear();
    }

    srcEdgeParts_.resize(srcPatch.nEdges());
    forAll(srcEdgeParts_, srcEdgei)
    {
        const edge& e = srcPatch.edges()[srcEdgei];
        const point c = e.centre(srcPatch.localPoints());
        srcEdgeParts_[srcEdgei] = part(Zero, c);
    }

    srcErrorParts_.resize(srcPatch.size());
    forAll(srcErrorParts_, srcFacei)
    {
        srcErrorParts_[srcFacei] =
            part(Zero, srcPatch.faceCentres()[srcFacei]);
    }

    tgtCouples_.resize(tgtPatch.size());
    forAll(tgtLocalSrcFaces_, i)
    {
        tgtCouples_[i].clear();
    }

    srcTriPoints_ = List<triFaceList>(srcPatch.size());
    srcTriFaceEdges_ = List<List<FixedList<label, 3>>>(srcPatch.size());
    tgtTriPoints_ = List<triFaceList>(tgtPatch.size());
    tgtTriFaceEdges_ = List<List<FixedList<label, 3>>>(tgtPatch.size());

    srcFaceEdgeParts_.resize(srcPatch.size());
    forAll(srcFaceEdgeParts_, srcFacei)
    {
        srcFaceEdgeParts_[srcFacei].resize(srcPatch[srcFacei].size());
        forAll(srcFaceEdgeParts_[srcFacei], srcFaceEdgei)
        {
            const label srcEdgei =
                srcPatch.faceEdges()[srcFacei][srcFaceEdgei];
            srcFaceEdgeParts_[srcFacei][srcFaceEdgei] = srcEdgeParts_[srcEdgei];
        }
    }

    if (debug)
    {
        debugPoints_.clear();
        debugFaces_.clear();
        debugFaceSrcFaces_.clear();
        debugFaceTgtFaces_.clear();
        debugFaceSides_.clear();
    }
}


Foam::labelList Foam::patchToPatches::intersection::finaliseLocal
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch
)
{
    const labelList newToOldLocalTgtFace =
        patchToPatch::finaliseLocal
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch
        );

    tgtCouples_ = List<DynamicList<couple>>(tgtCouples_, newToOldLocalTgtFace);

    return newToOldLocalTgtFace;
}


void Foam::patchToPatches::intersection::rDistributeTgt
(
    const primitiveOldTimePatch& tgtPatch
)
{
    patchToPatch::rDistributeTgt(tgtPatch);

    patchToPatchTools::rDistributeListList
    (
        tgtPatch.size(),
        tgtMapPtr_(),
        tgtCouples_
    );
}


Foam::label Foam::patchToPatches::intersection::finalise
(
    const primitiveOldTimePatch& srcPatch,
    const vectorField& srcPointNormals,
    const vectorField& srcPointNormals0,
    const primitiveOldTimePatch& tgtPatch,
    const transformer& tgtToSrc
)
{
    const label nCouples =
        patchToPatch::finalise
        (
            srcPatch,
            srcPointNormals,
            srcPointNormals0,
            tgtPatch,
            tgtToSrc
        );

    // Convert face-edge-parts to edge-parts
    labelList srcEdgeNParts(srcEdgeParts_.size(), 0);
    forAll(srcEdgeParts_, srcEdgei)
    {
        const edge& e = srcPatch.edges()[srcEdgei];

        srcEdgeParts_[srcEdgei] = part();

        forAll(srcPatch.edgeFaces()[srcEdgei], i)
        {
            const label srcFacei = srcPatch.edgeFaces()[srcEdgei][i];
            const label srcFaceEdgei =
                findIndex(srcPatch.faceEdges()[srcFacei], srcEdgei);

            const edge fe =
                srcPatch.localFaces()[srcFacei].faceEdge(srcFaceEdgei);

            if (edge::compare(e, fe) > 0)
            {
                srcEdgeParts_[srcEdgei] +=
                    srcFaceEdgeParts_[srcFacei][srcFaceEdgei];
            }
            else
            {
                srcEdgeParts_[srcEdgei] -=
                    srcFaceEdgeParts_[srcFacei][srcFaceEdgei];
            }

            srcEdgeNParts[srcEdgei] ++;
        }
    }
    forAll(srcEdgeParts_, srcEdgei)
    {
        srcEdgeParts_[srcEdgei].area /= srcEdgeNParts[srcEdgei];
    }

    // Add the difference between the face-edge-part and the edge-part into the
    // face-error-parts
    forAll(srcEdgeParts_, srcEdgei)
    {
        const edge& e = srcPatch.edges()[srcEdgei];

        forAll(srcPatch.edgeFaces()[srcEdgei], i)
        {
            const label srcFacei = srcPatch.edgeFaces()[srcEdgei][i];
            const label srcFaceEdgei =
                findIndex(srcPatch.faceEdges()[srcFacei], srcEdgei);

            const edge fe =
                srcPatch.localFaces()[srcFacei].faceEdge(srcFaceEdgei);

            if (edge::compare(e, fe) > 0)
            {
                srcErrorParts_[srcFacei] -= srcEdgeParts_[srcEdgei];
            }
            else
            {
                srcErrorParts_[srcFacei] += srcEdgeParts_[srcEdgei];
            }

            srcErrorParts_[srcFacei] +=
                srcFaceEdgeParts_[srcFacei][srcFaceEdgei];
        }
    }

    // Transform the target couples back to the target side
    if (!isNull(tgtToSrc))
    {
        forAll(tgtCouples_, tgtFacei)
        {
            forAll(tgtCouples_[tgtFacei], i)
            {
                couple& c = tgtCouples_[tgtFacei][i];

                c.area = tgtToSrc.invTransform(c.area);
                c.centre = tgtToSrc.invTransformPosition(c.centre);
                c.nbr.area = tgtToSrc.invTransform(c.nbr.area);
                c.nbr.centre = tgtToSrc.invTransformPosition(c.nbr.centre);
            }
        }
    }

    // Calculate coverage and total areas on both sides
    auto coverage = []
    (
        const primitivePatch& patch,
        const List<DynamicList<couple>>& couples,
        scalar& area,
        scalar& coupleArea,
        List<scalar>& coverage
    )
    {
        area = 0;
        coupleArea = 0;
        coverage.resize(patch.size());

        forAll(patch, facei)
        {
            const scalar magA = mag(patch.faceAreas()[facei]);

            vector aCouple = Zero;
            forAll(couples[facei], i)
            {
                aCouple += couples[facei][i].area;
            }
            const scalar magACouple = mag(aCouple);

            area += magA;
            coupleArea += magACouple;
            coverage[facei] = magACouple/magA;
        }

        reduce(area, sumOp<scalar>());
        reduce(coupleArea, sumOp<scalar>());
    };
    scalar srcArea = 0, srcCoupleArea = 0;
    scalar tgtArea = 0, tgtCoupleArea = 0;
    coverage(srcPatch, srcCouples_, srcArea, srcCoupleArea, srcCoverage_);
    coverage(tgtPatch, tgtCouples_, tgtArea, tgtCoupleArea, tgtCoverage_);

    // Clear the triangulation workspace
    srcTriPoints_.clear();
    srcTriFaceEdges_.clear();
    tgtTriPoints_.clear();
    tgtTriFaceEdges_.clear();

    // Clear face-edge-parts
    srcFaceEdgePart_.clear();
    tgtFaceEdgePart_.clear();
    srcFaceEdgeParts_.clear();

    // Checking and reporting
    if (nCouples != 0)
    {
        scalarField srcOpenness(srcPatch.size());
        scalarField srcError(srcPatch.size());
        scalarField srcDepth(srcPatch.size());
        scalarField srcAngle(srcPatch.size());
        forAll(srcPatch, srcFacei)
        {
            const vector& a = srcPatch.faceAreas()[srcFacei];
            const scalar magA = mag(a);
            const point& c = srcPatch.faceCentres()[srcFacei];

            couple Cpl(part(Zero, c), part(Zero, c));
            forAll(srcCouples_[srcFacei], srcTgtFacei)
            {
                const couple& cpl = srcCouples_[srcFacei][srcTgtFacei];

                Cpl += cpl;
                Cpl.nbr += cpl.nbr;
            }

            vector projectionA = Zero;
            scalar projectionV = 0;
            forAll(srcCouples_[srcFacei], srcTgtFacei)
            {
                const couple& cpl = srcCouples_[srcFacei][srcTgtFacei];

                projectionA += cpl.nbr.area;
                projectionV +=
                    - (cpl.area/3 & (cpl.centre - Cpl.centre))
                    + (cpl.nbr.area/3 & (cpl.nbr.centre - Cpl.centre));
            }
            forAll(srcPatch.faceEdges()[srcFacei], srcFaceEdgei)
            {
                const label srcEdgei =
                    srcPatch.faceEdges()[srcFacei][srcFaceEdgei];

                const edge& e = srcPatch.edges()[srcEdgei];
                const edge fe =
                    srcPatch.localFaces()[srcFacei].faceEdge(srcFaceEdgei);

                const scalar sign = edge::compare(e, fe);

                projectionA += sign*srcEdgeParts_[srcEdgei].area;
                projectionV +=
                    sign*srcEdgeParts_[srcEdgei].area/3
                  & (srcEdgeParts_[srcEdgei].centre - Cpl.centre);
            }
            projectionA += srcErrorParts_[srcFacei].area;
            projectionV +=
                srcErrorParts_[srcFacei].area/3
              & (srcErrorParts_[srcFacei].centre - Cpl.centre);

            const vector aHat = normalised(a);
            const vector aOppHat = normalised(a - Cpl.area + Cpl.nbr.area);
            srcAngle[srcFacei] =
                radToDeg(acos(min(max(aHat & aOppHat, -1), +1)));
            srcOpenness[srcFacei] = mag(projectionA - Cpl.area)/magA;
            srcError[srcFacei] = mag(srcErrorParts_[srcFacei].area)/magA;
            srcDepth[srcFacei] = mag(projectionV)/pow3(sqrt(magA));
        }

        reduce(tgtArea, sumOp<scalar>());
        reduce(tgtCoupleArea, sumOp<scalar>());

        Info<< indent << "Source min/average/max coverage = "
            << gMin(srcCoverage_) << '/' << srcCoupleArea/srcArea << '/'
            << gMax(srcCoverage_) << endl
            << indent << "Target min/average/max coverage = "
            << gMin(tgtCoverage_) << '/' << tgtCoupleArea/tgtArea << '/'
            << gMax(tgtCoverage_) << endl
            << indent << "Source average openness/error/depth/angle = "
            << gAverage(srcOpenness) << '/' << gAverage(srcError) << '/'
            << gAverage(srcDepth) << '/' << gAverage(srcAngle) << endl
            << indent << "Source max openness/error/depth/angle = "
            << gMax(srcOpenness) << '/' << gMax(srcError) << '/'
            << gMax(srcDepth) << '/' << gMax(srcAngle) << endl;

        if (debug)
        {
            word name = patchToPatch::typeName + '_' + typeName;

            if (Pstream::parRun())
            {
                name += "_proc" + Foam::name(Pstream::myProcNo());
            }

            Info<< indent << "Writing intersected faces to "
                << name + ".vtk" << endl;
            vtkWritePolyData::write
            (
                name + ".vtk",
                name,
                false,
                debugPoints_,
                labelList(),
                labelListList(),
                debugFaces_,
                "srcFace", false, Field<label>(debugFaceSrcFaces_),
                "tgtFace", false, Field<label>(debugFaceTgtFaces_),
                "side", false, Field<label>(debugFaceSides_)
            );

            debugPoints_.clear();
            debugFaces_.clear();
            debugFaceSrcFaces_.clear();
            debugFaceTgtFaces_.clear();
            debugFaceSides_.clear();

            Info<< indent << "Writing source patch to "
                << name + "_srcPatch.vtk" << endl;
            vtkWritePolyData::write
            (
                name + "_srcPatch" + ".vtk",
                name + "_srcPatch",
                false,
                srcPatch.localPoints(),
                labelList(),
                labelListList(),
                srcPatch.localFaces(),
                "coverage", false, scalarField(srcCoverage_),
                "openness", false, srcOpenness,
                "error", false, srcError,
                "depth", false, srcDepth,
                "angle", false, srcAngle,
                "normals", true, srcPointNormals
            );

            Info<< indent << "Writing target patch to "
                << name + "_tgtPatch.vtk" << endl;
            vtkWritePolyData::write
            (
                name + "_tgtPatch" + ".vtk",
                name + "_tgtPatch",
                false,
                tgtPatch.localPoints(),
                labelList(),
                labelListList(),
                tgtPatch.localFaces(),
                "coverage", false, scalarField(tgtCoverage_)
            );
        }
    }

    return nCouples;
}


Foam::tmpNrc<Foam::List<Foam::DynamicList<Foam::scalar>>>
Foam::patchToPatches::intersection::srcWeights() const
{
    List<DynamicList<scalar>>* resultPtr
    (
        new List<DynamicList<scalar>>(srcCouples_.size())
    );
    List<DynamicList<scalar>>& result = *resultPtr;

    forAll(srcCouples_, srcFacei)
    {
        result[srcFacei].resize(srcCouples_[srcFacei].size());
        scalar aSum = 0;

        forAll(srcCouples_[srcFacei], i)
        {
            const scalar a = mag(srcCouples_[srcFacei][i].area);
            result[srcFacei][i] = a;
            aSum += a;
        }

        forAll(srcCouples_[srcFacei], i)
        {
            result[srcFacei][i] *=
                min(max(srcCoverage_[srcFacei], small), scalar(1))/aSum;
        }
    }

    return tmpNrc<List<DynamicList<scalar>>>(resultPtr);
}


Foam::tmpNrc<Foam::List<Foam::DynamicList<Foam::scalar>>>
Foam::patchToPatches::intersection::tgtWeights() const
{
    List<DynamicList<scalar>>* resultPtr
    (
        new List<DynamicList<scalar>>(tgtCouples_.size())
    );
    List<DynamicList<scalar>>& result = *resultPtr;

    forAll(tgtCouples_, tgtFacei)
    {
        result[tgtFacei].resize(tgtCouples_[tgtFacei].size());
        scalar aSum = 0;

        forAll(tgtCouples_[tgtFacei], i)
        {
            const scalar a = mag(tgtCouples_[tgtFacei][i].area);
            result[tgtFacei][i] = a;
            aSum += a;
        }

        forAll(tgtCouples_[tgtFacei], i)
        {
            result[tgtFacei][i] *=
                min(max(tgtCoverage_[tgtFacei], small), scalar(1))/aSum;
        }
    }

    return tmpNrc<List<DynamicList<scalar>>>(resultPtr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchToPatches::intersection::intersection(const bool reverse)
:
    patchToPatch(reverse),

    srcCouples_(),
    srcCoverage_(),
    srcEdgeParts_(),
    srcErrorParts_(),
    tgtCouples_(),
    tgtCoverage_(),

    triEngine_(),

    srcTriPoints_(),
    srcTriFaceEdges_(),
    tgtTriPoints_(),
    tgtTriFaceEdges_(),

    ictSrcPoints_(),
    ictSrcPointNormals_(),
    ictTgtPoints_(),
    ictPointLocations_(),

    srcFaceEdgePart_(),
    tgtFaceEdgePart_(),

    srcFaceEdgeParts_(),

    debugPoints_(),
    debugFaces_(),
    debugFaceSrcFaces_(),
    debugFaceTgtFaces_(),
    debugFaceSides_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchToPatches::intersection::~intersection()
{}


// ************************************************************************* //
