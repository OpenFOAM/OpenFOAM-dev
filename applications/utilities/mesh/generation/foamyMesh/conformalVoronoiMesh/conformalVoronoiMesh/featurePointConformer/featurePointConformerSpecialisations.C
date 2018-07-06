/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "featurePointConformer.H"
#include "vectorTools.H"
#include "pointFeatureEdgesTypes.H"
#include "conformalVoronoiMesh.H"
#include "pointConversion.H"

using namespace Foam::vectorTools;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::featurePointConformer::createSpecialisedFeaturePoint
(
    const extendedFeatureEdgeMesh& feMesh,
    const labelList& pEds,
    const pointFeatureEdgesTypes& pFEdgesTypes,
    const List<extendedFeatureEdgeMesh::edgeStatus>& allEdStat,
    const label ptI,
    DynamicList<Vb>& pts
) const
{
    if
    (
        !pFEdgesTypes.found(extendedFeatureEdgeMesh::EXTERNAL)
     || !pFEdgesTypes.found(extendedFeatureEdgeMesh::INTERNAL)
    )
    {
        return false;
    }

    if
    (
        pFEdgesTypes[extendedFeatureEdgeMesh::EXTERNAL] == 2
     && pFEdgesTypes[extendedFeatureEdgeMesh::INTERNAL] == 1
     && pEds.size() == 3
    )
    {
        if (debug) Info<< "nExternal == 2 && nInternal == 1" << endl;

        const Foam::point& featPt = feMesh.points()[ptI];

        if
        (
            Pstream::parRun()
         && !foamyHexMesh_.decomposition().positionOnThisProcessor(featPt)
        )
        {
            return false;
        }

        label nVert = foamyHexMesh_.number_of_vertices();

        const label initialNumOfPoints = pts.size();

        const scalar ppDist = foamyHexMesh_.pointPairDistance(featPt);

        const vectorField& normals = feMesh.normals();

        const labelListList& edgeNormals = feMesh.edgeNormals();

        label concaveEdgeI = -1;
        labelList convexEdgesI(2, label(-1));
        label nConvex = 0;

        forAll(pEds, i)
        {
            const extendedFeatureEdgeMesh::edgeStatus& eS = allEdStat[i];

            if (eS == extendedFeatureEdgeMesh::INTERNAL)
            {
                concaveEdgeI = pEds[i];
            }
            else if (eS == extendedFeatureEdgeMesh::EXTERNAL)
            {
                convexEdgesI[nConvex++] = pEds[i];
            }
            else if (eS == extendedFeatureEdgeMesh::FLAT)
            {
                WarningInFunction
                    << "Edge " << eS << " is flat"
                    << endl;
            }
            else
            {
                FatalErrorInFunction
                    << "Edge " << eS << " not concave/convex"
                    << exit(FatalError);
            }
        }

        const vector& concaveEdgePlaneANormal =
            normals[edgeNormals[concaveEdgeI][0]];

        const vector& concaveEdgePlaneBNormal =
            normals[edgeNormals[concaveEdgeI][1]];

        // Intersect planes parallel to the concave edge planes offset
        // by ppDist and the plane defined by featPt and the edge vector.
        plane planeA
        (
            featPt + ppDist*concaveEdgePlaneANormal,
            concaveEdgePlaneANormal
        );

        plane planeB
        (
            featPt + ppDist*concaveEdgePlaneBNormal,
            concaveEdgePlaneBNormal
        );

        const vector& concaveEdgeDir = feMesh.edgeDirection
        (
            concaveEdgeI,
            ptI
        );

        // Todo,needed later but want to get rid of this.
        const Foam::point concaveEdgeLocalFeatPt =
            featPt + ppDist*concaveEdgeDir;

        // Finding the nearest point on the intersecting line to the edge
        // point. Floating point errors often occur using planePlaneIntersect

        plane planeF(concaveEdgeLocalFeatPt, concaveEdgeDir);

        const Foam::point concaveEdgeExternalPt = planeF.planePlaneIntersect
        (
            planeA,
            planeB
        );

        // Redefine planes to be on the feature surfaces to project through

        planeA = plane(featPt, concaveEdgePlaneANormal);

        planeB = plane(featPt, concaveEdgePlaneBNormal);

        const Foam::point internalPtA =
            concaveEdgeExternalPt
          - 2.0*planeA.distance(concaveEdgeExternalPt)
            *concaveEdgePlaneANormal;

        pts.append
        (
            Vb
            (
                internalPtA,
                foamyHexMesh_.vertexCount() + pts.size(),
                Vb::vtInternalFeaturePoint,
                Pstream::myProcNo()
            )
        );

        const label internalPtAIndex(pts.last().index());

        const Foam::point internalPtB =
            concaveEdgeExternalPt
          - 2.0*planeB.distance(concaveEdgeExternalPt)
            *concaveEdgePlaneBNormal;

        pts.append
        (
            Vb
            (
                internalPtB,
                foamyHexMesh_.vertexCount() + pts.size(),
                Vb::vtInternalFeaturePoint,
                Pstream::myProcNo()
            )
        );

        const label internalPtBIndex(pts.last().index());

        // Add the external points

        Foam::point externalPtD;
        Foam::point externalPtE;

        vector convexEdgePlaneCNormal(Zero);
        vector convexEdgePlaneDNormal(Zero);

        const labelList& concaveEdgeNormals = edgeNormals[concaveEdgeI];
        const labelList& convexEdgeANormals = edgeNormals[convexEdgesI[0]];
        const labelList& convexEdgeBNormals = edgeNormals[convexEdgesI[1]];

        forAll(concaveEdgeNormals, edgeNormalI)
        {
            bool convexEdgeA = false;
            bool convexEdgeB = false;

            forAll(convexEdgeANormals, edgeAnormalI)
            {
                const vector& concaveNormal
                    = normals[concaveEdgeNormals[edgeNormalI]];
                const vector& convexNormal
                    = normals[convexEdgeANormals[edgeAnormalI]];

                if (debug)
                {
                    Info<< "Angle between vectors = "
                        << degAngleBetween(concaveNormal, convexNormal) << endl;
                }

                // Need a looser tolerance, because sometimes adjacent triangles
                // on the same surface will be slightly out of alignment.
                if (areParallel(concaveNormal, convexNormal, tolParallel))
                {
                    convexEdgeA = true;
                }
            }

            forAll(convexEdgeBNormals, edgeBnormalI)
            {
                const vector& concaveNormal
                    = normals[concaveEdgeNormals[edgeNormalI]];
                const vector& convexNormal
                    = normals[convexEdgeBNormals[edgeBnormalI]];

                if (debug)
                {
                    Info<< "Angle between vectors = "
                        << degAngleBetween(concaveNormal, convexNormal) << endl;
                }

                // Need a looser tolerance, because sometimes adjacent triangles
                // on the same surface will be slightly out of alignment.
                if (areParallel(concaveNormal, convexNormal, tolParallel))
                {
                    convexEdgeB = true;
                }
            }

            if ((convexEdgeA && convexEdgeB) || (!convexEdgeA && !convexEdgeB))
            {
                WarningInFunction
                    << "Both or neither of the convex edges share the concave "
                    << "edge's normal."
                    << " convexEdgeA = " << convexEdgeA
                    << " convexEdgeB = " << convexEdgeB
                    << endl;

                // Remove points that have just been added before returning
                for (label i = 0; i < 2; ++i)
                {
                    pts.remove();
                    nVert--;
                }

                return false;
            }

            if (convexEdgeA)
            {
                forAll(convexEdgeANormals, edgeAnormalI)
                {
                    const vector& concaveNormal
                        = normals[concaveEdgeNormals[edgeNormalI]];
                    const vector& convexNormal
                        = normals[convexEdgeANormals[edgeAnormalI]];

                    if
                    (
                        !areParallel(concaveNormal, convexNormal, tolParallel)
                    )
                    {
                        convexEdgePlaneCNormal = convexNormal;

                        plane planeC(featPt, convexEdgePlaneCNormal);

                        externalPtD =
                            internalPtA
                          + 2.0*planeC.distance(internalPtA)
                           *convexEdgePlaneCNormal;

                        pts.append
                        (
                            Vb
                            (
                                externalPtD,
                                foamyHexMesh_.vertexCount() + pts.size(),
                                Vb::vtExternalFeaturePoint,
                                Pstream::myProcNo()
                            )
                        );

                        ftPtPairs_.addPointPair
                        (
                            internalPtAIndex,
                            pts.last().index()
                        );
                    }
                }
            }

            if (convexEdgeB)
            {
                forAll(convexEdgeBNormals, edgeBnormalI)
                {
                    const vector& concaveNormal
                        = normals[concaveEdgeNormals[edgeNormalI]];
                    const vector& convexNormal
                        = normals[convexEdgeBNormals[edgeBnormalI]];

                    if
                    (
                        !areParallel(concaveNormal, convexNormal, tolParallel)
                    )
                    {
                        convexEdgePlaneDNormal = convexNormal;

                        plane planned(featPt, convexEdgePlaneDNormal);

                        externalPtE =
                            internalPtB
                          + 2.0*planned.distance(internalPtB)
                           *convexEdgePlaneDNormal;

                        pts.append
                        (
                            Vb
                            (
                                externalPtE,
                                foamyHexMesh_.vertexCount() + pts.size(),
                                Vb::vtExternalFeaturePoint,
                                Pstream::myProcNo()
                            )
                        );

                        ftPtPairs_.addPointPair
                        (
                            internalPtBIndex,
                            pts.last().index()
                        );
                    }
                }
            }
        }

        pts.append
        (
            Vb
            (
                concaveEdgeExternalPt,
                foamyHexMesh_.vertexCount() + pts.size(),
                Vb::vtExternalFeaturePoint,
                Pstream::myProcNo()
            )
        );

        ftPtPairs_.addPointPair
        (
            internalPtBIndex,
            pts.last().index()
        );

        ftPtPairs_.addPointPair
        (
            internalPtAIndex,
            pts.last().index()
        );

        const label concaveEdgeExternalPtIndex(pts.last().index());

        const scalar totalAngle = radToDeg
        (
            constant::mathematical::pi
          + radAngleBetween(concaveEdgePlaneANormal, concaveEdgePlaneBNormal)
        );

        if (totalAngle > foamyHexMeshControls_.maxQuadAngle())
        {
            // Add additional mitreing points
            // scalar angleSign = 1.0;


            vector convexEdgesPlaneNormal =
                0.5*(convexEdgePlaneCNormal + convexEdgePlaneDNormal);

            plane planeM(featPt, convexEdgesPlaneNormal);

//            if
//            (
//                geometryToConformTo_.outside
//                (
//                    featPt - convexEdgesPlaneNormal*ppDist
//                )
//            )
//            {
//                angleSign = -1.0;
//            }

//            scalar phi =
//                angleSign*acos(concaveEdgeDir & -convexEdgesPlaneNormal);
//
//            scalar guard =
//            (
//                1.0 + sin(phi)*ppDist/mag
//                (
//                    concaveEdgeLocalFeatPt - concaveEdgeExternalPt
//                )
//            )/cos(phi) - 1.0;

            const Foam::point internalPtF =
                concaveEdgeExternalPt
            //+ (2.0 + guard)*(concaveEdgeLocalFeatPt - concaveEdgeExternalPt);
              + 2.0*(concaveEdgeLocalFeatPt - concaveEdgeExternalPt);

            pts.append
            (
                Vb
                (
                    internalPtF,
                    foamyHexMesh_.vertexCount() + pts.size(),
                    Vb::vtInternalFeaturePoint,
                    Pstream::myProcNo()
                )
            );

            const label internalPtFIndex(pts.last().index());

            ftPtPairs_.addPointPair
            (
                concaveEdgeExternalPtIndex,
                pts.last().index()
            );

            const Foam::point externalPtG =
                internalPtF
              + 2.0*planeM.distance(internalPtF)*convexEdgesPlaneNormal;

            pts.append
            (
                Vb
                (
                    externalPtG,
                    foamyHexMesh_.vertexCount() + pts.size(),
                    Vb::vtExternalFeaturePoint,
                    Pstream::myProcNo()
                )
            );

            ftPtPairs_.addPointPair
            (
                internalPtFIndex,
                pts.last().index()
            );
        }

        if (debug)
        {
            for (label ptI = initialNumOfPoints; ptI < pts.size(); ++ptI)
            {
                Info<< "Point " << ptI << " : ";
                meshTools::writeOBJ(Info, topoint(pts[ptI].point()));
            }
        }

        return true;
    }
    else if
    (
        pFEdgesTypes[extendedFeatureEdgeMesh::EXTERNAL] == 1
     && pFEdgesTypes[extendedFeatureEdgeMesh::INTERNAL] == 2
     && pEds.size() == 3
    )
    {
        if (debug)
        {
            Info<< "nExternal == 1 && nInternal == 2" << endl;
        }

        const Foam::point& featPt = feMesh.points()[ptI];

        if
        (
            Pstream::parRun()
         && !foamyHexMesh_.decomposition().positionOnThisProcessor(featPt)
        )
        {
            return false;
        }

        label nVert = foamyHexMesh_.number_of_vertices();

        const label initialNumOfPoints = pts.size();

        const scalar ppDist = foamyHexMesh_.pointPairDistance(featPt);

        const vectorField& normals = feMesh.normals();

        const labelListList& edgeNormals = feMesh.edgeNormals();

        label convexEdgeI = -1;
        labelList concaveEdgesI(2, label(-1));
        label nConcave = 0;

        forAll(pEds, i)
        {
            const extendedFeatureEdgeMesh::edgeStatus& eS = allEdStat[i];

            if (eS == extendedFeatureEdgeMesh::EXTERNAL)
            {
                convexEdgeI = pEds[i];
            }
            else if (eS == extendedFeatureEdgeMesh::INTERNAL)
            {
                concaveEdgesI[nConcave++] = pEds[i];
            }
            else if (eS == extendedFeatureEdgeMesh::FLAT)
            {
                WarningInFunction
                    << "Edge " << eS << " is flat"
                    << endl;
            }
            else
            {
                FatalErrorInFunction
                    << "Edge " << eS << " not concave/convex"
                    << exit(FatalError);
            }
        }

        const vector& convexEdgePlaneANormal =
            normals[edgeNormals[convexEdgeI][0]];

        const vector& convexEdgePlaneBNormal =
            normals[edgeNormals[convexEdgeI][1]];

        // Intersect planes parallel to the concave edge planes offset
        // by ppDist and the plane defined by featPt and the edge vector.
        plane planeA
        (
            featPt - ppDist*convexEdgePlaneANormal,
            convexEdgePlaneANormal
        );

        plane planeB
        (
            featPt - ppDist*convexEdgePlaneBNormal,
            convexEdgePlaneBNormal
        );

        const vector& convexEdgeDir = feMesh.edgeDirection
        (
            convexEdgeI,
            ptI
        );

        // Todo,needed later but want to get rid of this.
        const Foam::point convexEdgeLocalFeatPt =
            featPt + ppDist*convexEdgeDir;

        // Finding the nearest point on the intersecting line to the edge
        // point. Floating point errors often occur using planePlaneIntersect

        plane planeF(convexEdgeLocalFeatPt, convexEdgeDir);

        const Foam::point convexEdgeExternalPt = planeF.planePlaneIntersect
        (
            planeA,
            planeB
        );

        // Redefine planes to be on the feature surfaces to project through

        planeA = plane(featPt, convexEdgePlaneANormal);

        planeB = plane(featPt, convexEdgePlaneBNormal);

        const Foam::point internalPtA =
            convexEdgeExternalPt
          + 2.0*planeA.distance(convexEdgeExternalPt)
           *convexEdgePlaneANormal;

        pts.append
        (
            Vb
            (
                internalPtA,
                foamyHexMesh_.vertexCount() + pts.size(),
                Vb::vtExternalFeaturePoint,
                Pstream::myProcNo()
            )
        );

        const label internalPtAIndex(pts.last().index());

        const Foam::point internalPtB =
            convexEdgeExternalPt
          + 2.0*planeB.distance(convexEdgeExternalPt)
           *convexEdgePlaneBNormal;

        pts.append
        (
            Vb
            (
                internalPtB,
                foamyHexMesh_.vertexCount() + pts.size(),
                Vb::vtExternalFeaturePoint,
                Pstream::myProcNo()
            )
        );

        const label internalPtBIndex(pts.last().index());

        // Add the internal points

        Foam::point externalPtD;
        Foam::point externalPtE;

        vector concaveEdgePlaneCNormal(Zero);
        vector concaveEdgePlaneDNormal(Zero);

        const labelList& convexEdgeNormals = edgeNormals[convexEdgeI];
        const labelList& concaveEdgeANormals = edgeNormals[concaveEdgesI[0]];
        const labelList& concaveEdgeBNormals = edgeNormals[concaveEdgesI[1]];

        forAll(convexEdgeNormals, edgeNormalI)
        {
            bool concaveEdgeA = false;
            bool concaveEdgeB = false;

            forAll(concaveEdgeANormals, edgeAnormalI)
            {
                const vector& convexNormal
                    = normals[convexEdgeNormals[edgeNormalI]];
                const vector& concaveNormal
                    = normals[concaveEdgeANormals[edgeAnormalI]];

                if (debug)
                {
                    Info<< "Angle between vectors = "
                        << degAngleBetween(convexNormal, concaveNormal) << endl;
                }

                // Need a looser tolerance, because sometimes adjacent triangles
                // on the same surface will be slightly out of alignment.
                if (areParallel(convexNormal, concaveNormal, tolParallel))
                {
                    concaveEdgeA = true;
                }
            }

            forAll(concaveEdgeBNormals, edgeBnormalI)
            {
                const vector& convexNormal
                    = normals[convexEdgeNormals[edgeNormalI]];
                const vector& concaveNormal
                    = normals[concaveEdgeBNormals[edgeBnormalI]];

                if (debug)
                {
                    Info<< "Angle between vectors = "
                        << degAngleBetween(convexNormal, concaveNormal) << endl;
                }

                // Need a looser tolerance, because sometimes adjacent triangles
                // on the same surface will be slightly out of alignment.
                if (areParallel(convexNormal, concaveNormal, tolParallel))
                {
                    concaveEdgeB = true;
                }
            }

            if
            (
                (concaveEdgeA && concaveEdgeB)
             || (!concaveEdgeA && !concaveEdgeB)
            )
            {
                WarningInFunction
                    << "Both or neither of the concave edges share the convex "
                    << "edge's normal."
                    << " concaveEdgeA = " << concaveEdgeA
                    << " concaveEdgeB = " << concaveEdgeB
                    << endl;

                // Remove points that have just been added before returning
                for (label i = 0; i < 2; ++i)
                {
                    pts.remove();
                    nVert--;
                }

                return false;
            }

            if (concaveEdgeA)
            {
                forAll(concaveEdgeANormals, edgeAnormalI)
                {
                    const vector& convexNormal
                        = normals[convexEdgeNormals[edgeNormalI]];
                    const vector& concaveNormal
                        = normals[concaveEdgeANormals[edgeAnormalI]];

                    if
                    (
                        !areParallel(convexNormal, concaveNormal, tolParallel)
                    )
                    {
                        concaveEdgePlaneCNormal = concaveNormal;

                        plane planeC(featPt, concaveEdgePlaneCNormal);

                        externalPtD =
                            internalPtA
                          - 2.0*planeC.distance(internalPtA)
                           *concaveEdgePlaneCNormal;

                        pts.append
                        (
                            Vb
                            (
                                externalPtD,
                                foamyHexMesh_.vertexCount() + pts.size(),
                                Vb::vtInternalFeaturePoint,
                                Pstream::myProcNo()
                            )
                        );

                        ftPtPairs_.addPointPair
                        (
                            internalPtAIndex,
                            pts.last().index()
                        );
                    }
                }
            }

            if (concaveEdgeB)
            {
                forAll(concaveEdgeBNormals, edgeBnormalI)
                {
                    const vector& convexNormal
                        = normals[convexEdgeNormals[edgeNormalI]];
                    const vector& concaveNormal
                        = normals[concaveEdgeBNormals[edgeBnormalI]];

                    if
                    (
                        !areParallel(convexNormal, concaveNormal, tolParallel)
                    )
                    {
                        concaveEdgePlaneDNormal = concaveNormal;

                        plane planned(featPt, concaveEdgePlaneDNormal);

                        externalPtE =
                            internalPtB
                          - 2.0*planned.distance(internalPtB)
                           *concaveEdgePlaneDNormal;

                        pts.append
                        (
                            Vb
                            (
                                externalPtE,
                                foamyHexMesh_.vertexCount() + pts.size(),
                                Vb::vtInternalFeaturePoint,
                                Pstream::myProcNo()
                            )
                        );

                        ftPtPairs_.addPointPair
                        (
                            internalPtBIndex,
                            pts.last().index()
                        );
                    }
                }
            }
        }

        pts.append
        (
            Vb
            (
                convexEdgeExternalPt,
                foamyHexMesh_.vertexCount() + pts.size(),
                Vb::vtInternalFeaturePoint,
                Pstream::myProcNo()
            )
        );

        ftPtPairs_.addPointPair
        (
            internalPtBIndex,
            pts.last().index()
        );

        ftPtPairs_.addPointPair
        (
            internalPtAIndex,
            pts.last().index()
        );

        const scalar totalAngle = radToDeg
        (
            constant::mathematical::pi
          + radAngleBetween(convexEdgePlaneANormal, convexEdgePlaneBNormal)
        );

        if (totalAngle > foamyHexMeshControls_.maxQuadAngle())
        {
            // Add additional mitreing points
            // scalar angleSign = 1.0;


            vector convexEdgesPlaneNormal =
                0.5*(concaveEdgePlaneCNormal + concaveEdgePlaneDNormal);

            plane planeM(featPt, convexEdgesPlaneNormal);

//            if
//            (
//                geometryToConformTo_.outside
//                (
//                    featPt - convexEdgesPlaneNormal*ppDist
//                )
//            )
//            {
//                angleSign = -1.0;
//            }

//            scalar phi =
//                angleSign*acos(concaveEdgeDir & -convexEdgesPlaneNormal);
//
//            scalar guard =
//            (
//                1.0 + sin(phi)*ppDist/mag
//                (
//                    concaveEdgeLocalFeatPt - concaveEdgeExternalPt
//                )
//            )/cos(phi) - 1.0;

            const Foam::point internalPtF =
                convexEdgeExternalPt
            //+ (2.0 + guard)*(concaveEdgeLocalFeatPt - concaveEdgeExternalPt);
              + 2.0*(convexEdgeLocalFeatPt - convexEdgeExternalPt);

            pts.append
            (
                Vb
                (
                    internalPtF,
                    foamyHexMesh_.vertexCount() + pts.size(),
                    Vb::vtExternalFeaturePoint,
                    Pstream::myProcNo()
                )
            );

            ftPtPairs_.addPointPair
            (
                pts[pts.size() - 2].index(),
                pts.last().index()
            );

            const Foam::point externalPtG =
                internalPtF
              - 2.0*planeM.distance(internalPtF)*convexEdgesPlaneNormal;

            pts.append
            (
                Vb
                (
                    externalPtG,
                    foamyHexMesh_.vertexCount() + pts.size(),
                    Vb::vtInternalFeaturePoint,
                    Pstream::myProcNo()
                )
            );

            ftPtPairs_.addPointPair
            (
                pts[pts.size() - 2].index(),
                pts.last().index()
            );
        }

        if (debug)
        {
            for (label ptI = initialNumOfPoints; ptI < pts.size(); ++ptI)
            {
                Info<< "Point " << ptI << " "
                    << indexedVertexEnum::vertexTypeNames_[pts[ptI].type()]
                    << " : ";

                meshTools::writeOBJ(Info, topoint(pts[ptI].point()));
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
