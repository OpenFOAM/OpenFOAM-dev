/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "conformalVoronoiMesh.H"
#include "vectorTools.H"
#include "triangle.H"
#include "tetrahedron.H"
#include "ConstCirculator.H"
#include "DelaunayMeshTools.H"
#include "OBJstream.H"

using namespace Foam::vectorTools;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::createEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
) const
{
    if (foamyHexMeshControls().circulateEdges())
    {
        createEdgePointGroupByCirculating(feMesh, edHit, pts);
    }
    else
    {
        label edgeI = edHit.index();

        extendedFeatureEdgeMesh::edgeStatus edStatus =
            feMesh.getEdgeStatus(edgeI);

        switch (edStatus)
        {
            case extendedFeatureEdgeMesh::EXTERNAL:
            {
                createExternalEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::INTERNAL:
            {
                createInternalEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::FLAT:
            {
                createFlatEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::OPEN:
            {
                createOpenEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::MULTIPLE:
            {
                createMultipleEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::NONE:
            {
                break;
            }
        }
    }
}


bool Foam::conformalVoronoiMesh::meshableRegion
(
    const plane::side side,
    const extendedFeatureEdgeMesh::sideVolumeType volType
) const
{
    switch (volType)
    {
        case extendedFeatureEdgeMesh::INSIDE:
        {
            return (side == plane::FLIP) ? true : false;
        }
        case extendedFeatureEdgeMesh::OUTSIDE:
        {
            return (side == plane::NORMAL) ? true : false;
        }
        case extendedFeatureEdgeMesh::BOTH:
        {
            return true;
        }
        case extendedFeatureEdgeMesh::NEITHER:
        {
            return false;
        }
        default:
        {
            return false;
        }
    }
}


bool Foam::conformalVoronoiMesh::regionIsInside
(
    const extendedFeatureEdgeMesh::sideVolumeType volTypeA,
    const vector& normalA,
    const extendedFeatureEdgeMesh::sideVolumeType volTypeB,
    const vector& normalB,
    const vector& masterPtVec
) const
{
    plane::side sideA
    (
        ((masterPtVec & normalA) <= 0) ? plane::FLIP : plane::NORMAL
    );

    plane::side sideB
    (
        ((masterPtVec & normalB) <= 0) ? plane::FLIP : plane::NORMAL
    );

    const bool meshableRegionA = meshableRegion(sideA, volTypeA);
    const bool meshableRegionB = meshableRegion(sideB, volTypeB);

    if (meshableRegionA == meshableRegionB)
    {
        return meshableRegionA;
    }
    else
    {
        WarningInFunction
            << ""
            << endl;

        return false;
    }
}


void Foam::conformalVoronoiMesh::createEdgePointGroupByCirculating
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
)  const
{
    typedef Foam::indexedVertexEnum::vertexType         vertexType;
    typedef extendedFeatureEdgeMesh::sideVolumeType     sideVolumeType;

    const Foam::point& edgePt = edHit.hitPoint();
    const label edgeI = edHit.index();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const vector& edDir = feMesh.edgeDirections()[edgeI];
    const labelList& edNormalIs = feMesh.edgeNormals()[edgeI];
    const labelList& feNormalDirections = feMesh.normalDirections()[edgeI];

    const List<sideVolumeType>& normalVolumeTypes = feMesh.normalVolumeTypes();

    ConstCirculator<labelList> circ(edNormalIs);
    ConstCirculator<labelList> circNormalDirs(feNormalDirections);

    Map<Foam::point> masterPoints;
    Map<vertexType> masterPointsTypes;
    Map<plane> masterPointReflectionsPrev;
    Map<plane> masterPointReflectionsNext;

//    Info<< "Edge = " << edHit << ", edDir = " << edDir << endl;
//    Info<< "    edNorms = " << edNormalIs.size() << endl;

    bool addedMasterPreviously = false;
    label initialRegion = -1;

    if (circ.size()) do
    {
        const sideVolumeType volType = normalVolumeTypes[circ()];
        const sideVolumeType nextVolType = normalVolumeTypes[circ.next()];

        const vector& normal = feNormals[circ()];
        const vector& nextNormal = feNormals[circ.next()];

        vector normalDir = (normal ^ edDir);
        normalDir *= circNormalDirs()/mag(normalDir);

        vector nextNormalDir = (nextNormal ^ edDir);
        nextNormalDir *= circNormalDirs.next()/mag(nextNormalDir);

//        Info<< "    " << circ() << " " << circ.next() << nl
//            << "    " << circNormalDirs() << " " << circNormalDirs.next()
//            << nl
//            << "    normals = " << normal << " " << nextNormal << nl
//            << " normalDirs = " << normalDir << " " << nextNormalDir << nl
//            << "      cross = " << (normalDir ^ nextNormalDir) << nl
//            << "    "
//            << extendedFeatureEdgeMesh::sideVolumeTypeNames_[volType] << " "
//            << extendedFeatureEdgeMesh::sideVolumeTypeNames_[nextVolType]
//            << endl;

        // Calculate master point
        vector masterPtVec(normalDir + nextNormalDir);
        masterPtVec /= mag(masterPtVec) + small;

        if
        (
            ((normalDir ^ nextNormalDir) & edDir) < small
         || mag(masterPtVec) < small
        )
        {
//            Info<< "    IGNORE REGION" << endl;
            addedMasterPreviously = false;

            if
            (
                circ.size() == 2
             && mag((normal & nextNormal) - 1) < small
            )
            {
                const vector n = 0.5*(normal + nextNormal);

                const vector s = ppDist*(edDir ^ n);

                if (volType == extendedFeatureEdgeMesh::BOTH)
                {
                    createBafflePointPair(ppDist, edgePt + s, n, true, pts);
                    createBafflePointPair(ppDist, edgePt - s, n, true, pts);
                }
                else
                {
                    WarningInFunction
                        << "Failed to insert flat/open edge as volType is "
                        << extendedFeatureEdgeMesh::sideVolumeTypeNames_
                           [
                               volType
                           ]
                        << endl;
                }

                break;
            }

            continue;
        }

        const Foam::point masterPt = edgePt + ppDist*masterPtVec;

        // Check that region is inside or outside
        const bool inside =
            regionIsInside
            (
                volType,
                normal,
                nextVolType,
                nextNormal,
                masterPtVec
            );

        // Specialize for size = 1 && baffle
        if (mag((normalDir & nextNormalDir) - 1) < small)
        {
            if (inside)
            {
//                Info<< "Specialize for size 1 and baffle" << endl;

                vector s = ppDist*(edDir ^ normal);

                plane facePlane(edgePt, normal);

                Foam::point pt1 = edgePt + s + ppDist*normal;
                Foam::point pt2 = edgePt - s + ppDist*normal;

                Foam::point pt3 = facePlane.mirror(pt1);
                Foam::point pt4 = facePlane.mirror(pt2);

                pts.append(Vb(pt1, Vb::vtInternalFeatureEdge));
                pts.append(Vb(pt2, Vb::vtInternalFeatureEdge));
                pts.append(Vb(pt3, Vb::vtInternalFeatureEdge));
                pts.append(Vb(pt4, Vb::vtInternalFeatureEdge));

                break;
            }
            else
            {
                WarningInFunction
                    << "Faces are parallel but master point is not inside"
                    << endl;
            }
        }

        if (!addedMasterPreviously)
        {
            if (initialRegion == -1)
            {
                initialRegion = circ.nRotations();
            }

            addedMasterPreviously = true;

            masterPoints.insert(circ.nRotations(), masterPt);
            masterPointsTypes.insert
            (
                circ.nRotations(),
                inside
              ? Vb::vtInternalFeatureEdge
              : Vb::vtExternalFeatureEdge
            );

            masterPointReflectionsPrev.insert
            (
                circ.nRotations(),
                plane(edgePt, normal)
            );

            masterPointReflectionsNext.insert
            (
                circ.nRotations(),
                plane(edgePt, nextNormal)
            );
        }
        else if (addedMasterPreviously)
        {
            addedMasterPreviously = true;

            masterPointReflectionsNext.erase(circ.nRotations() - 1);

            // Shift the master point to be normal to the plane between it and
            // the previous master point
            // Should be the intersection of the normal and the plane with the
            // new master point in it.

            plane p(masterPoints[circ.nRotations() - 1], normalDir);
            plane::ray r(edgePt, masterPt - edgePt);

            scalar cutPoint = p.normalIntersect(r);

            masterPoints.insert
            (
                circ.nRotations(),
                edgePt + cutPoint*(masterPt - edgePt)
            );

            masterPointsTypes.insert
            (
                circ.nRotations(),
                inside
              ? Vb::vtInternalFeatureEdge
              : Vb::vtExternalFeatureEdge
            );

            masterPointReflectionsNext.insert
            (
                circ.nRotations(),
                plane(edgePt, nextNormal)
            );
        }

        if
        (
            masterPoints.size() > 1
         && inside
         && circ.nRotations() == circ.size() - 1
        )
        {
            if (initialRegion == 0)
            {
                plane p(masterPoints[initialRegion], nextNormalDir);
                plane::ray r(edgePt, masterPt - edgePt);

                scalar cutPoint = p.normalIntersect(r);

                masterPoints[circ.nRotations()] =
                    edgePt + cutPoint*(masterPt - edgePt);

                // Remove the first reflection plane if we are no longer
                // circulating

                masterPointReflectionsPrev.erase(initialRegion);
                masterPointReflectionsNext.erase(circ.nRotations());
            }
            else
            {

            }
        }
    }
    while
    (
        circ.circulate(CirculatorBase::CLOCKWISE),
        circNormalDirs.circulate(CirculatorBase::CLOCKWISE)
    );


    forAllConstIter(Map<Foam::point>, masterPoints, iter)
    {
        const Foam::point& pt = masterPoints[iter.key()];
        const vertexType ptType = masterPointsTypes[iter.key()];

//        Info<< "    Adding Master " << iter.key() << " " << pt << " "
//            << indexedVertexEnum::vertexTypeNames_[ptType] << endl;

        pts.append(Vb(pt, ptType));

        const vertexType reflectedPtType =
        (
            ptType == Vb::vtInternalFeatureEdge
          ? Vb::vtExternalFeatureEdge
          : Vb::vtInternalFeatureEdge
        );

        if (masterPointReflectionsPrev.found(iter.key()))
        {
            const Foam::point reflectedPt =
                masterPointReflectionsPrev[iter.key()].mirror(pt);

//            Info<< "        Adding Prev " << reflectedPt << " "
//                << indexedVertexEnum::vertexTypeNames_[reflectedPtType]
//                << endl;

            pts.append(Vb(reflectedPt, reflectedPtType));
        }

        if (masterPointReflectionsNext.found(iter.key()))
        {
            const Foam::point reflectedPt =
               masterPointReflectionsNext[iter.key()].mirror(pt);

//            Info<< "        Adding Next " << reflectedPt << " "
//                << indexedVertexEnum::vertexTypeNames_[reflectedPtType]
//                << endl;

            pts.append(Vb(reflectedPt, reflectedPtType));
        }
    }

//    pts.append(Vb(edgePt, Vb::vtExternalFeatureEdge));
}


void Foam::conformalVoronoiMesh::createExternalEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
) const
{
    const Foam::point& edgePt = edHit.hitPoint();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];
    const List<extendedFeatureEdgeMesh::sideVolumeType>& normalVolumeTypes =
        feMesh.normalVolumeTypes();

    // As this is an external edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    const extendedFeatureEdgeMesh::sideVolumeType& volTypeA =
        normalVolumeTypes[edNormalIs[0]];

    const extendedFeatureEdgeMesh::sideVolumeType& volTypeB =
        normalVolumeTypes[edNormalIs[1]];

    if (areParallel(nA, nB))
    {
        // The normals are nearly parallel, so this is too sharp a feature to
        // conform to.
        return;
    }

    // Normalised distance of reference point from edge point
    vector refVec((nA + nB)/(1 + (nA & nB)));

    if (magSqr(refVec) > sqr(5.0))
    {
        // Limit the size of the conformation
        ppDist *= 5.0/mag(refVec);

        // Pout<< nl << "createExternalEdgePointGroup limit "
        //     << "edgePt " << edgePt << nl
        //     << "refVec " << refVec << nl
        //     << "mag(refVec) " << mag(refVec) << nl
        //     << "ppDist " << ppDist << nl
        //     << "nA " << nA << nl
        //     << "nB " << nB << nl
        //     << "(nA & nB) " << (nA & nB) << nl
        //     << endl;
    }

    // Convex. So refPt will be inside domain and hence a master point
    Foam::point refPt = edgePt - ppDist*refVec;

    // Insert the master point pairing the the first slave

    if (!geometryToConformTo_.inside(refPt))
    {
        return;
    }

    pts.append
    (
        Vb
        (
            refPt,
            vertexCount() + pts.size(),
            Vb::vtInternalFeatureEdge,
            Pstream::myProcNo()
        )
    );

    // Insert the slave points by reflecting refPt in both faces.
    // with each slave referring to the master

    Foam::point reflectedA = refPt + 2*ppDist*nA;
    pts.append
    (
        Vb
        (
            reflectedA,
            vertexCount() + pts.size(),
            (
                volTypeA == extendedFeatureEdgeMesh::BOTH
              ? Vb::vtInternalFeatureEdge
              : Vb::vtExternalFeatureEdge
            ),
            Pstream::myProcNo()
        )
    );

    Foam::point reflectedB = refPt + 2*ppDist*nB;
    pts.append
    (
        Vb
        (
            reflectedB,
            vertexCount() + pts.size(),
            (
                volTypeB == extendedFeatureEdgeMesh::BOTH
              ? Vb::vtInternalFeatureEdge
              : Vb::vtExternalFeatureEdge
            ),
            Pstream::myProcNo()
        )
    );

    ptPairs_.addPointPair
    (
        pts[pts.size() - 3].index(),
        pts[pts.size() - 1].index()
    );

    ptPairs_.addPointPair
    (
        pts[pts.size() - 3].index(),
        pts[pts.size() - 2].index()
    );
}


void Foam::conformalVoronoiMesh::createInternalEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
) const
{
    const Foam::point& edgePt = edHit.hitPoint();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];
    const List<extendedFeatureEdgeMesh::sideVolumeType>& normalVolumeTypes =
        feMesh.normalVolumeTypes();

    // As this is an external edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    const extendedFeatureEdgeMesh::sideVolumeType& volTypeA =
        normalVolumeTypes[edNormalIs[0]];

//    const extendedFeatureEdgeMesh::sideVolumeType& volTypeB =
//        normalVolumeTypes[edNormalIs[1]];

    if (areParallel(nA, nB))
    {
        // The normals are nearly parallel, so this is too sharp a feature to
        // conform to.
        return;
    }

    // Normalised distance of reference point from edge point
    vector refVec((nA + nB)/(1 + (nA & nB)));

    if (magSqr(refVec) > sqr(5.0))
    {
        // Limit the size of the conformation
        ppDist *= 5.0/mag(refVec);

        // Pout<< nl << "createInternalEdgePointGroup limit "
        //     << "edgePt " << edgePt << nl
        //     << "refVec " << refVec << nl
        //     << "mag(refVec) " << mag(refVec) << nl
        //     << "ppDist " << ppDist << nl
        //     << "nA " << nA << nl
        //     << "nB " << nB << nl
        //     << "(nA & nB) " << (nA & nB) << nl
        //     << endl;
    }

    // Concave. master and reflected points inside the domain.
    Foam::point refPt = edgePt - ppDist*refVec;

    // Generate reflected master to be outside.
    Foam::point reflMasterPt = refPt + 2*(edgePt - refPt);

    // Reflect reflMasterPt in both faces.
    Foam::point reflectedA = reflMasterPt - 2*ppDist*nA;

    Foam::point reflectedB = reflMasterPt - 2*ppDist*nB;

    scalar totalAngle =
        radToDeg(constant::mathematical::pi + radAngleBetween(nA, nB));

    // Number of quadrants the angle should be split into
    int nQuads = int(totalAngle/foamyHexMeshControls().maxQuadAngle()) + 1;

    // The number of additional master points needed to obtain the
    // required number of quadrants.
    int nAddPoints = min(max(nQuads - 2, 0), 2);

    // Add number_of_vertices() at insertion of first vertex to all numbers:
    // Result for nAddPoints 1 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     2
    // reflectedB   1     2
    // reflMasterPt 2     0

    // Result for nAddPoints 1 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     3
    // reflectedB   1     3
    // refPt        2     3
    // reflMasterPt 3     0

    // Result for nAddPoints 2 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     4
    // reflectedB   1     4
    // reflectedAa  2     4
    // reflectedBb  3     4
    // reflMasterPt 4     0

    if
    (
        !geometryToConformTo_.inside(reflectedA)
     || !geometryToConformTo_.inside(reflectedB)
    )
    {
        return;
    }

    // Master A is inside.
    pts.append
    (
        Vb
        (
            reflectedA,
            vertexCount() + pts.size(),
            Vb::vtInternalFeatureEdge,
            Pstream::myProcNo()
        )
    );

    // Master B is inside.
    pts.append
    (
        Vb
        (
            reflectedB,
            vertexCount() + pts.size(),
            Vb::vtInternalFeatureEdge,
            Pstream::myProcNo()
        )
    );

    // Slave is outside.
    pts.append
    (
        Vb
        (
            reflMasterPt,
            vertexCount() + pts.size(),
            (
                volTypeA == extendedFeatureEdgeMesh::BOTH
              ? Vb::vtInternalFeatureEdge
              : Vb::vtExternalFeatureEdge
            ),
            Pstream::myProcNo()
        )
    );

    ptPairs_.addPointPair
    (
        pts[pts.size() - 2].index(),
        pts[pts.size() - 1].index()
    );

    ptPairs_.addPointPair
    (
        pts[pts.size() - 3].index(),
        pts[pts.size() - 1].index()
    );

    if (nAddPoints == 1)
    {
        // One additional point is the reflection of the slave point,
        // i.e. the original reference point
        pts.append
        (
            Vb
            (
                refPt,
                vertexCount() + pts.size(),
                Vb::vtInternalFeatureEdge,
                Pstream::myProcNo()
            )
        );
    }
    else if (nAddPoints == 2)
    {
        Foam::point reflectedAa = refPt + ppDist*nB;
        pts.append
        (
            Vb
            (
                reflectedAa,
                vertexCount() + pts.size(),
                Vb::vtInternalFeatureEdge,
                Pstream::myProcNo()
            )
        );

        Foam::point reflectedBb = refPt + ppDist*nA;
        pts.append
        (
            Vb
            (
                reflectedBb,
                vertexCount() + pts.size(),
                Vb::vtInternalFeatureEdge,
                Pstream::myProcNo()
            )
        );
    }
}


void Foam::conformalVoronoiMesh::createFlatEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
) const
{
    const Foam::point& edgePt = edHit.hitPoint();

    const scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];
    const List<extendedFeatureEdgeMesh::sideVolumeType>& normalVolumeTypes =
        feMesh.normalVolumeTypes();

    // As this is a flat edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    // Average normal to remove any bias to one side, although as this
    // is a flat edge, the normals should be essentially the same
    const vector n = 0.5*(nA + nB);

    // Direction along the surface to the control point, sense of edge
    // direction not important, as +s and -s can be used because this
    // is a flat edge
    vector s = ppDist*(feMesh.edgeDirections()[edHit.index()] ^ n);

    if (normalVolumeTypes[edNormalIs[0]] == extendedFeatureEdgeMesh::OUTSIDE)
    {
        createPointPair(ppDist, edgePt + s, -n, true, pts);
        createPointPair(ppDist, edgePt - s, -n, true, pts);
    }
    else if (normalVolumeTypes[edNormalIs[0]] == extendedFeatureEdgeMesh::BOTH)
    {
        createBafflePointPair(ppDist, edgePt + s, n, true, pts);
        createBafflePointPair(ppDist, edgePt - s, n, true, pts);
    }
    else
    {
        createPointPair(ppDist, edgePt + s, n, true, pts);
        createPointPair(ppDist, edgePt - s, n, true, pts);
    }
}


void Foam::conformalVoronoiMesh::createOpenEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
) const
{
    // Assume it is a baffle and insert flat edge point pairs
    const Foam::point& edgePt = edHit.hitPoint();

    const scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    if (edNormalIs.size() == 1)
    {
//        Info<< "Inserting open edge point group around " << edgePt << endl;
//        Info<< "    ppDist    = " << ppDist << nl
//            << "    edNormals = " << edNormalIs
//            << endl;

        const vector& n = feNormals[edNormalIs[0]];

        vector s = ppDist*(feMesh.edgeDirections()[edHit.index()] ^ n);

        plane facePlane(edgePt, n);

        const label initialPtsSize = pts.size();

        if
        (
            !geometryToConformTo_.inside(edgePt)
        )
        {
            return;
        }

        createBafflePointPair(ppDist, edgePt - s, n, true, pts);
        createBafflePointPair(ppDist, edgePt + s, n, false, pts);

        for (label ptI = initialPtsSize; ptI < pts.size(); ++ptI)
        {
            pts[ptI].type() = Vb::vtInternalFeatureEdge;
        }
    }
    else
    {
        Info<< "NOT INSERTING OPEN EDGE POINT GROUP WITH MORE THAN 1 "
            << "EDGE NORMAL, NOT IMPLEMENTED" << endl;
    }
}


void Foam::conformalVoronoiMesh::createMultipleEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
) const
{
//    Info<< "NOT INSERTING MULTIPLE EDGE POINT GROUP, NOT IMPLEMENTED" << endl;

    const Foam::point& edgePt = edHit.hitPoint();

    const scalar ppDist = pointPairDistance(edgePt);

    const vector edDir = feMesh.edgeDirections()[edHit.index()];

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];
    const labelList& normalDirs = feMesh.normalDirections()[edHit.index()];

    const List<extendedFeatureEdgeMesh::sideVolumeType>& normalVolumeTypes =
        feMesh.normalVolumeTypes();

    labelList nNormalTypes(4, label(0));

    forAll(edNormalIs, edgeNormalI)
    {
        const extendedFeatureEdgeMesh::sideVolumeType sType =
            normalVolumeTypes[edNormalIs[edgeNormalI]];

        nNormalTypes[sType]++;
    }

    if (nNormalTypes[extendedFeatureEdgeMesh::BOTH] == 4)
    {
        label masterEdgeNormalIndex = -1;

        forAll(edNormalIs, edgeNormalI)
        {
            const extendedFeatureEdgeMesh::sideVolumeType sType =
                normalVolumeTypes[edNormalIs[edgeNormalI]];

            if (sType == extendedFeatureEdgeMesh::BOTH)
            {
                masterEdgeNormalIndex = edgeNormalI;
                break;
            }
        }

        const vector& n = feNormals[edNormalIs[masterEdgeNormalIndex]];

        label nDir = normalDirs[masterEdgeNormalIndex];

        vector normalDir =
            (feNormals[edNormalIs[masterEdgeNormalIndex]] ^ edDir);
        normalDir *= nDir/mag(normalDir);

        Foam::point pt1 = edgePt + ppDist*normalDir + ppDist*n;
        Foam::point pt2 = edgePt + ppDist*normalDir - ppDist*n;

        plane plane3(edgePt, normalDir);

        Foam::point pt3 = plane3.mirror(pt1);
        Foam::point pt4 = plane3.mirror(pt2);

        pts.append
        (
            Vb
            (
                pt1,
                vertexCount() + pts.size(),
                Vb::vtInternalSurface,
                Pstream::myProcNo()
            )
        );
        pts.append
        (
            Vb
            (
                pt2,
                vertexCount() + pts.size(),
                Vb::vtInternalSurface,
                Pstream::myProcNo()
            )
        );

        ptPairs_.addPointPair
        (
            pts[pts.size() - 2].index(), // external 0 -> slave
            pts[pts.size() - 1].index()
        );

        pts.append
        (
            Vb
            (
                pt3,
                vertexCount() + pts.size(),
                Vb::vtInternalSurface,
                Pstream::myProcNo()
            )
        );

        ptPairs_.addPointPair
        (
            pts[pts.size() - 3].index(), // external 0 -> slave
            pts[pts.size() - 1].index()
        );

        pts.append
        (
            Vb
            (
                pt4,
                vertexCount() + pts.size(),
                Vb::vtInternalSurface,
                Pstream::myProcNo()
            )
        );

        ptPairs_.addPointPair
        (
            pts[pts.size() - 3].index(), // external 0 -> slave
            pts[pts.size() - 1].index()
        );

        ptPairs_.addPointPair
        (
            pts[pts.size() - 2].index(), // external 0 -> slave
            pts[pts.size() - 1].index()
        );
    }
    else if
    (
        nNormalTypes[extendedFeatureEdgeMesh::BOTH] == 1
     && nNormalTypes[extendedFeatureEdgeMesh::INSIDE] == 2
    )
    {
        label masterEdgeNormalIndex = -1;

        forAll(edNormalIs, edgeNormalI)
        {
            const extendedFeatureEdgeMesh::sideVolumeType sType =
                normalVolumeTypes[edNormalIs[edgeNormalI]];

            if (sType == extendedFeatureEdgeMesh::BOTH)
            {
                masterEdgeNormalIndex = edgeNormalI;
                break;
            }
        }

        const vector& n = feNormals[edNormalIs[masterEdgeNormalIndex]];

        label nDir = normalDirs[masterEdgeNormalIndex];

        vector normalDir =
            (feNormals[edNormalIs[masterEdgeNormalIndex]] ^ edDir);
        normalDir *= nDir/mag(normalDir);

        const label nextNormalI =
            (masterEdgeNormalIndex + 1) % edNormalIs.size();
        if ((normalDir & feNormals[edNormalIs[nextNormalI]]) > 0)
        {
            normalDir *= -1;
        }

        Foam::point pt1 = edgePt + ppDist*normalDir + ppDist*n;
        Foam::point pt2 = edgePt + ppDist*normalDir - ppDist*n;

        plane plane3(edgePt, normalDir);

        Foam::point pt3 = plane3.mirror(pt1);
        Foam::point pt4 = plane3.mirror(pt2);

        pts.append
        (
            Vb
            (
                pt1,
                vertexCount() + pts.size(),
                Vb::vtInternalSurface,
                Pstream::myProcNo()
            )
        );
        pts.append
        (
            Vb
            (
                pt2,
                vertexCount() + pts.size(),
                Vb::vtInternalSurface,
                Pstream::myProcNo()
            )
        );

        ptPairs_.addPointPair
        (
            pts[pts.size() - 2].index(), // external 0 -> slave
            pts[pts.size() - 1].index()
        );

        pts.append
        (
            Vb
            (
                pt3,
                vertexCount() + pts.size(),
                Vb::vtExternalSurface,
                Pstream::myProcNo()
            )
        );

        ptPairs_.addPointPair
        (
            pts[pts.size() - 3].index(), // external 0 -> slave
            pts[pts.size() - 1].index()
        );

        pts.append
        (
            Vb
            (
                pt4,
                vertexCount() + pts.size(),
                Vb::vtExternalSurface,
                Pstream::myProcNo()
            )
        );

        ptPairs_.addPointPair
        (
            pts[pts.size() - 3].index(), // external 0 -> slave
            pts[pts.size() - 1].index()
        );
    }


//    // As this is a flat edge, there are two normals by definition
//    const vector& nA = feNormals[edNormalIs[0]];
//    const vector& nB = feNormals[edNormalIs[1]];
//
//    // Average normal to remove any bias to one side, although as this
//    // is a flat edge, the normals should be essentially the same
//    const vector n = 0.5*(nA + nB);
//
//    // Direction along the surface to the control point, sense of edge
//    // direction not important, as +s and -s can be used because this
//    // is a flat edge
//    vector s = ppDist*(feMesh.edgeDirections()[edHit.index()] ^ n);
//
//    createBafflePointPair(ppDist, edgePt + s, n, true, pts);
//    createBafflePointPair(ppDist, edgePt - s, n, true, pts);
}


void Foam::conformalVoronoiMesh::insertFeaturePoints(bool distribute)
{
    Info<< nl
        << "Inserting feature points" << endl;

    const label preFeaturePointSize(number_of_vertices());

    if (Pstream::parRun() && distribute)
    {
        ftPtConformer_.distribute(decomposition());
    }

    const List<Vb>& ftPtVertices = ftPtConformer_.featurePointVertices();

    // Insert the created points directly as already distributed.
    Map<label> oldToNewIndices =
        this->DelaunayMesh<Delaunay>::insertPoints(ftPtVertices, true);

    ftPtConformer_.reIndexPointPairs(oldToNewIndices);

    label nFeatureVertices = number_of_vertices() - preFeaturePointSize;
    reduce(nFeatureVertices, sumOp<label>());

    Info<< "    Inserted " << nFeatureVertices << " feature vertices" << endl;
}


//Foam::scalar Foam::conformalVoronoiMesh::pyramidVolume
//(
//    const Foam::point& apex,
//    const Foam::point& a,
//    const Foam::point& b,
//    const Foam::point& c,
//    const bool printInfo
//) const
//{
//    triPointRef tri(a, b, c);
//
//    tetPointRef tet(tri.a(), tri.b(), tri.c(), apex);
//
//    scalar volume = tet.mag();
//
////    scalar volume = (1.0/3.0)*constant::mathematical::pi;
////
////    K::Circle_3 circle(toPoint(a), toPoint(b), toPoint(c));
////
////    scalar height = mag(topoint(circle.center()) - apex);
////
////    volume *= circle.squared_radius()*height;
//
//    if (printInfo)
//    {
//        Info<< "Calculating volume of pyramid..." << nl
//            << "    Apex      : " << apex << nl
//            << "    Point a   : " << a << nl
//            << "    Point b   : " << b << nl
//            << "    Point c   : " << c << nl
//            << "        Center  : " << tri.centre() << nl
//            << "        Volume  : " << volume << endl;
//    }
//
//    return volume;
//}


//void Foam::conformalVoronoiMesh::createPyramidMasterPoint
//(
//    const Foam::point& apex,
//    const vectorField& edgeDirections,
//    Foam::point& masterPoint,
//    vectorField& norms
//) const
//{
//    pointField basePoints(edgeDirections.size() + 1);
//
//    forAll(edgeDirections, eI)
//    {
//        basePoints[eI] = edgeDirections[eI] + apex;
//    }
//
//    basePoints[edgeDirections.size() + 1] = apex;
//
//    face f(identity(edgeDirections.size()));
//
//    pyramidPointFaceRef p(f, apex);
//
//    const scalar ppDist = pointPairDistance(apex);
//
//
//    vector unitDir = f.centre();
//    unitDir /= mag(unitDir);
//
//    masterPoint = apex + ppDist*unitDir;
//
//    norms.setSize(edgeDirections.size());
//
//    forAll(norms, nI)
//    {
//        norms[nI] =
//    }
//}


//void Foam::conformalVoronoiMesh::createConvexConcaveFeaturePoints
//(
//    DynamicList<Foam::point>& pts,
//    DynamicList<label>& indices,
//    DynamicList<label>& types
//)
//{
//    const PtrList<extendedFeatureEdgeMesh>& feMeshes
//    (
//        geometryToConformTo_.features()
//    );
//
//    forAll(feMeshes, i)
//    {
//        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);
//
//        for
//        (
//            label ptI = feMesh.convexStart();
//            ptI < feMesh.mixedStart();
//            ptI++
//        )
//        {
//            const Foam::point& apex = feMesh.points()[ptI];
//
//            if (!positionOnThisProc(apex))
//            {
//                continue;
//            }
//
//            const vectorField& featPtEdgeDirections
//                = feMesh.featurePointEdgeDirections(ptI);
//
//            Foam::point masterPoint;
//            vectorField tetNorms;
//
//            createPyramidMasterPoint
//            (
//                apex,
//                featPtEdgeDirections,
//                masterPoint,
//                tetNorms
//            );
//
//
//
//            // Result when the points are eventually inserted (example n = 4)
//            // Add number_of_vertices() at insertion of first vertex to all
//            // numbers:
//            // pt           index type
//            // internalPt   0     1
//            // externalPt0  1     0
//            // externalPt1  2     0
//            // externalPt2  3     0
//            // externalPt3  4     0
//
//            // Result when the points are eventually inserted (example n = 5)
//            // Add number_of_vertices() at insertion of first vertex to all
//            // numbers:
//            // pt           index type
//            // internalPt0  0     5
//            // internalPt1  1     5
//            // internalPt2  2     5
//            // internalPt3  3     5
//            // internalPt4  4     5
//            // externalPt   5     4
//
//            if (geometryToConformTo_.inside(masterPoint))
//            {
//
//            }
//            else
//            {
//
//            }
//
//            pts.append(masterPoint);
//            indices.append(0);
//            types.append(1);
//
//            label internalPtIndex = -1;
//
//            forAll(tetNorms, nI)
//            {
//                const vector& n = tetNorms[nI];
//
//                Foam::point reflectedPoint
//                    = reflectPoint(featPt, masterPoint, n);
//
//                pts.append(reflectedPoint);
//                indices.append(0);
//                types.append(internalPtIndex--);
//            }
//        }
//    }
//}
