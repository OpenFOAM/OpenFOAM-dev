/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

Description
    Functions with which to perform an intersection of a pair of triangles; the
    source and target. The source triangle is specified with three point
    locations and three point normals. It is projected along it's point normals
    onto the target triangle in order to calculate the intersection. The target
    triangle is specified as three point locations only.

SourceFiles
    triIntersect.C

\*---------------------------------------------------------------------------*/

#ifndef triIntersect_H
#define triIntersect_H

#include "barycentric2D.H"
#include "DynamicField.H"
#include "primitivePatch.H"
#include "triFace.H"
#include "triIntersectLocation.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace triIntersect
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write a VTK file of a polygon
template<class PointField>
void writePolygon
(
    const word& name,
    const PointField& ps
);

//- Write a VTK file of a triangle projection
void writeTriProjection
(
    const word& name,
    const FixedList<point, 3>& srcTriPs,
    const FixedList<vector, 3>& srcTriNs,
    const label nEdge = 20,
    const label nNormal = 20,
    const scalar lNormal = 0.5
);

//- Calculate the signed offset of a target point in relation to a projected
//  source edge
scalar srcEdgeTgtPointOffset
(
    const Pair<point>& srcPs,
    const Pair<vector>& srcNs,
    const point& tgtP
);

//- Calculate the signed offset of a target point in relation to a projected
//  source edge. Additional ordering argument can be used to account for
//  edge orientation and thereby guarantee consistency of the result.
scalar srcEdgeTgtPointOffset
(
    const Pair<point>& srcPs,
    const Pair<vector>& srcNs,
    const point& tgtP,
    const bool srcDirection
);

//- Calculate the signed offset of a projected source point in relation to a
//  target edge
scalar srcPointTgtEdgeOffset
(
    const point& srcP,
    const vector& srcN,
    const Pair<point>& tgtPs
);

//- Calculate the signed offset of a projected source point in relation to a
//  target edge. Additional ordering argument can be used to account for
//  edge orientation and thereby guarantee consistency of the result.
scalar srcPointTgtEdgeOffset
(
    const point& srcPs,
    const vector& srcNs,
    const Pair<point>& tgtP,
    const bool tgtDirection
);

//- Calculate the intersection of a target edge with a source edge's
//  projected surface. Return the local coordinates of the intersection along
//  both the source and target edges.
Pair<scalar> srcEdgeTgtEdgeIntersection
(
    const Pair<point>& srcPs,
    const Pair<vector>& srcNs,
    const Pair<point>& tgtPs
);

//- Calculate the intersection of a target point with a source triangle's
//  projected volume. Return the local coordinates of the intersection within
//  the source triangle.
barycentric2D srcTriTgtPointIntersection
(
    const FixedList<point, 3>& srcPs,
    const FixedList<vector, 3>& srcNs,
    const point& tgtP
);

//- Use the coordinates obtained from srcTriTgtPointIntersection to interpolate
//  a property within a source triangle's projected volume
template<class Type>
Type srcTriInterpolate
(
    const barycentric2D& y,
    const FixedList<Type, 3>& tgtPsis
);

//- Calculate the intersection of a projected source point with a target
//  triangle. Return the local coordinates of the intersection within the
//  target triangle.
barycentric2D srcPointTgtTriIntersection
(
    const point& srcP,
    const vector& srcN,
    const FixedList<point, 3>& tgtPs
);

//- Use the coordinates obtained from srcPointTgtTriIntersection to interpolate
//  a property within a target triangle
template<class Type>
Type tgtTriInterpolate
(
    const barycentric2D& y,
    const FixedList<Type, 3>& tgtPsis
);

//- Construct the intersection of a source triangle's projected volume and a
//  target triangle
void intersectTris
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
    const word& writePrefix=word::null
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace triIntersect
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "triIntersectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
